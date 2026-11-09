//============================================================================
/* Copyright (c) The Technical University of Denmark
 * Developed at the Computational Visualization Center (CVC), The University of Texas at Austin
 * Author: Ojaswa Sharma
 * E-mail: os@imm.dtu.dk
 * File: multiphaseseg.cu
 * Main multi-phase levelset implementation and host functions to the CUDA device functions
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License version 2.1 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 */
//============================================================================


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <MultiphaseSegmentation/multiphaseseg.h>
#include <MultiphaseSegmentation/private/multiphaseseg_global.cu.h>
#include <MultiphaseSegmentation/private/pde_update_kernels.cu.h>
#include <MultiphaseSegmentation/private/init_interface_kernels.cu.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include<map>
#include<utility>
#include<vector>

using namespace MultiphaseSegmentation;
using namespace std;
 
MPSegmentation::MPSegmentation(bool doInitCUDA):HOSegmentation(doInitCUDA) 
{
  //MultiPhase stuff
  nImplicits = 3; // Let's do it for this. Nevertheless this parameter should be specified by the user. No. of resulting classes = 2^nImplicits.
}

bool MPSegmentation::allocateHost(float *vol, float **phi)
{
  debugPrint("Allocating host memory... ");
  h_vol = vol;
  h_PHI = phi;//Memory for phi should already be allocated by the caller.
  if(!h_vol) return false;
  if(!h_PHI) return false;
  h_phi = *(h_PHI);
  
/*  int x, y, z;
  float l_I;
  int nxy = volSize.width*volSize.height;

    for(int i=0; i< nelem; i++) {
    z = i/nxy;
    y = (i - z*nxy)/volSize.width; 
    x = (i - z*nxy - y*volSize.width);
	debugPrint("at beginning z,y,x = %d %d %d\n", z, y, x);
//    l_class = 0;
    l_I = VOXEL_AT(h_vol, x, y, z);
	debugPrint("l_I = %f \n", l_I);
  }
*/ 
  //Allocate extra buffer for PHI
  h_BUFFER[0] = h_PHI; currentPhi = 0;
  h_BUFFER[1] = new float*[nImplicits];
  for(int i=0; i<nImplicits; i++)
  {
    h_BUFFER[1][i] = new float[nelem];
    if(h_BUFFER[1][i] == NULL)
    {
      debugPrint("\nFailed to allocate %lu bytes of memory for phi(%d)! Not enough memory to hold %d implicit surfaces.\nReduce the number of classes.\n", nelem*sizeof(float), i, nImplicits);
      for(int j=i; j>=0; j--)
        delete [](h_BUFFER[1][j]);
      return false;
    }     
  }

  //Average values for all possible classes
  c_avg = new float[1<<nImplicits];//Average values. One for each class.

  //Allocate memory for the spare subvolume on the host --needed by computations by the parent class.
  h_subvolSpare =  new float[subvolSize.width*subvolSize.height*subvolSize.depth];
  if(h_subvolSpare == NULL) {
    debugPrint("\nFailed to allocate %lu bytes of memory\n", subvolSize.width*subvolSize.height*subvolSize.depth*sizeof(float));
    for(int i=0; i<nImplicits; i++)
      delete [](h_BUFFER[1][i]);
    delete [](h_BUFFER[1]);
    return false;
  } 
  //Print total memory allocated.
  unsigned long l_totalmem = (nImplicits*nelem + 1<<nImplicits + subvolSize.width*subvolSize.height*subvolSize.depth)*sizeof(float);
  float l_totalmemf = l_totalmem/1024.0;
  if(l_totalmemf < 1.0)
    debugPrint("%lu bytes.\n", l_totalmem);
  else {
    l_totalmemf /=1024.0; 
    if(l_totalmemf < 1.0)
      debugPrint("%.3f KB.\n", l_totalmemf*1024.0);
    else {
      l_totalmemf /=1024.0; 
      if(l_totalmemf < 1.0)
        debugPrint("%.3f MB.\n", l_totalmemf*1024.0);
      else
        debugPrint("%.3f GB.\n", l_totalmemf);
    }
  }
  return true;
}

bool MPSegmentation::freeHost()
{
  for(int i=0; i<nImplicits; i++)
    delete [](h_BUFFER[1][i]);
  delete [](h_BUFFER[1]);
  delete []c_avg;
  delete []h_subvolSpare;
  return true;
}

bool MPSegmentation::allocateDevice()
{
  //Create 3D array
  debugPrint("Allocating device memory... \n");
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
  printMemInfo();

  //Perform a quick check on required mem!
  unsigned int freemem, totalmem, peakrequiredmem;
  if(!getFreeGPUMem(freemem, totalmem)) return false;
  int subvolmem = subvolSize.width*subvolSize.height*subvolSize.depth*sizeof(float);
  peakrequiredmem = subvolmem*(5 + nImplicits);  

  if(peakrequiredmem > totalmem) {
    debugPrint("On GPU peak required memory (%u bytes) is more than the total GPU memory (%u bytes). Please reduce the number of classes.\n", peakrequiredmem, totalmem);
    return false;  
  }
  if(peakrequiredmem > freemem) {
    debugPrint("On GPU peak required memory (%u bytes) is more than the available (%u bytes). Please free some GPU memory (perhaps by a video mode switch.)\n", peakrequiredmem, freemem);
    return false;
  }

  //Split nImplicits into nImplicitsX, nImplicitsy, nImplicitsZ so that the max texture size requirement is satisfed.
  debugPrint("For %d phi: ", nImplicits);
  unsigned int maxSubvolUnits = MAX_3D_TEX_SIZE/subvolDim;//  =16. so a total of 16^3 phi possible => 2^(16^3) classes!!! <limited by texture memory though>
  computeOptimalFactors(nImplicits, maxSubvolUnits, nImplicitsX, nImplicitsY, nImplicitsZ);
  cudaExtent phi_subvolSize = make_cudaExtent(subvolDim*nImplicitsX, subvolDim*nImplicitsY, subvolDim*nImplicitsZ);
  //Allocate mem and bind phiArray to 3D texture
  CUDA_SAFE_CALL(cudaMalloc3DArray(&d_phiArray, &channelDesc, phi_subvolSize));
  CHK_CUDA_ERR("\n");
  texPhi.normalized = false;
  texPhi.filterMode = cudaFilterModePoint;
  texPhi.addressMode[0] = cudaAddressModeWrap; 
  texPhi.addressMode[1] = cudaAddressModeWrap; 
  texPhi.addressMode[2] = cudaAddressModeWrap; 
  CUDA_SAFE_CALL(cudaBindTextureToArray(texPhi, d_phiArray, channelDesc));
  cuda_err = cudaGetLastError();
  if(cuda_err != cudaSuccess) {
    debugPrint("initCuda: %s\n", cudaGetErrorString(cuda_err)) ;
    return(false); 
  }
  debugPrint("allocated %d (= %d x %d x %d optimal factors) implicits.\n", nImplicitsX*nImplicitsY*nImplicitsZ, nImplicitsX, nImplicitsY, nImplicitsZ);

  //Allocate mem and bind volArray to 3D texture
  CUDA_SAFE_CALL(cudaMalloc3DArray(&d_volArray, &channelDesc, subvolSize));
  CHK_CUDA_ERR("\n");
  texVol.normalized = false;
  texVol.filterMode = cudaFilterModePoint;
  texVol.addressMode[0] = cudaAddressModeWrap; 
  texVol.addressMode[1] = cudaAddressModeWrap; 
  texVol.addressMode[2] = cudaAddressModeWrap; 
  CUDA_SAFE_CALL(cudaBindTextureToArray(texVol, d_volArray, channelDesc));
  cuda_err = cudaGetLastError();
  if(cuda_err != cudaSuccess) {
    freeDevice(1);
    debugPrint("initCuda: %s\n", cudaGetErrorString(cuda_err)) ;
    return(false); 
  }
  debugPrint("Allocated vol.\n");
  
  //Allocate memory for 3D device memory for kernel output
  cudaExtent pitchedVolSize = make_cudaExtent(subvolDim*sizeof(float), subvolDim, subvolDim); 
  CUDA_SAFE_CALL(cudaMalloc3D(&d_volPPtr, pitchedVolSize));
  cuda_err = cudaGetLastError();
  if(cuda_err != cudaSuccess) { 
    freeDevice(2);
    debugPrint("initCuda: %s\n", cudaGetErrorString(cuda_err)) ;
    return(false); 
  }
  debugPrint("Allocated kernel output array.\n");
  
  CUDA_SAFE_CALL( cudaMalloc((void**)&d_c_avg, (1<<nImplicits)*sizeof(float)));
  texAverageValues.normalized = false;
  texAverageValues.filterMode = cudaFilterModePoint;
  texAverageValues.addressMode[0] = cudaAddressModeWrap; 
  CUDA_SAFE_CALL(cudaBindTexture(0, texAverageValues, d_c_avg, channelDesc));
  cuda_err = cudaGetLastError();
  if(cuda_err != cudaSuccess) { 
    freeDevice(3);
    debugPrint("initCuda: %s\n", cudaGetErrorString(cuda_err)) ;
    return(false); 
  }
  debugPrint("Allocated average values array.\n");
  
  printMemInfo();
  return true;
}

bool MPSegmentation::freeDevice(int level)
{
  switch(level)
  {
    case -1:
      CUDA_SAFE_CALL(cudaUnbindTexture(texAverageValues));
      CUDA_SAFE_CALL(cudaFree(d_c_avg));
    case 3:
      CUDA_SAFE_CALL(cudaFree(d_volPPtr.ptr));
    case 2:
      CUDA_SAFE_CALL(cudaUnbindTexture(texVol)); //unbind texture
      CUDA_SAFE_CALL(cudaFreeArray(d_volArray));
    case 1: 
      CUDA_SAFE_CALL(cudaUnbindTexture(texPhi)); //unbind texture
      CUDA_SAFE_CALL(cudaFreeArray(d_phiArray));
    }
  return true;
}

bool MPSegmentation::initCuda()
{
  if(DoInitializeCUDA)
  {
    cuInit(0); //Any driver API can also be called now on...
    cudaSetDevice(cutGetMaxGflopsDeviceId());
  }
  CUT_SAFE_CALL(cutCreateTimer(&timer));
  CUT_SAFE_CALL(cutResetTimer(timer));
  return true;
}

bool MPSegmentation::runSolver(float *vol, float **phi, int _width, int _height, int _depth, MPLevelSetParams *MPLSParams)
{
  setParameters (_width, _height, _depth, MPLSParams);
  //Set up various volume sizes
  volSize = make_cudaExtent(datainfo.n[0], datainfo.n[1], datainfo.n[2]);
  nelem = volSize.width*volSize.height*volSize.depth;
  // set the volume size and number of subvolumes
  subvolIndicesExtents = make_cudaExtent(iDivUp(volSize.width-2, subvolDim-2), iDivUp(volSize.height-2, subvolDim-2), iDivUp(volSize.depth-2, subvolDim-2));
  subvolSize = make_cudaExtent(subvolDim, subvolDim, subvolDim);

  if(!allocateHost(vol, phi)) return false;
  if(!initCuda()) return false;
  if(!allocateDevice()) return false; 

  debugPrint("Total volume: %zdx%zdx%zd\n", volSize.width, volSize.height, volSize.depth);
  debugPrint("Total input volume: %dx%dx%d\n", datainfo.n_input[0], datainfo.n_input[1], datainfo.n_input[2]);
  debugPrint("Total subvolumes: %zdx%zdx%zd\n", subvolIndicesExtents.width, subvolIndicesExtents.height, subvolIndicesExtents.depth);

  bool result = solverMain();
  debugPrint("Exiting solver...\n");

  freeDevice();
  CUT_SAFE_CALL(cutDeleteTimer(timer));
  freeHost();
  return result;
}

int MPSegmentation::solverMain()//entry point into the CUDA solver. 
{
  double t0, t1;
  t0 = getTime();
  normalizeVolume();
  initInterfaceMulti();

  if(nIter > 0) {
    debugPrint("Solving the PDE...\n");
    for (int i = 0; i<nIter; i++)
      if(!solve()) return false;
  }
  t1 = getTime();
  debugPrint("Total evolution time: %f sec. Average: %f sec/iteration in %d iterations.\n", t1 - t0, (t1-t0)/float(nIter), nIter);

  return true;
}

bool MPSegmentation::solve()
{
  //Solve the PDE
  if(iter_num < nIter) {
    iter_num++;
    debugPrint("Iteration: %d\n", iter_num);
    if(!computeSDTEikonal()) return false;
    if(!computeAverageIntensities()) return false;
    if(!PDEUpdate()) return false; //update phi and volume via level set PDE's
  }
  debugPrint("\tSolver done.\n");
  return true;
}

bool MPSegmentation::initInterfaceMulti()
{
  /// Initialize level set interface (to a cuboid)
  debugPrint("Initializing the level set interface... ");
  float maxDist = -float(multi_init_dr);
  int l_box = 2*(multi_init_r + multi_init_dr) + multi_init_s;
  int l_nx = volSize.width/l_box;
  int l_ny = volSize.height/l_box;
  int l_nz = volSize.depth/l_box;
  float l_offx = float(volSize.width - l_nx*l_box); 
  float l_offy = float(volSize.height - l_ny*l_box); 
  float l_offz = float(volSize.depth - l_nz*l_box); 
  float radius_reach = float(multi_init_s)/2. + float(multi_init_r + multi_init_dr);
  float rad0 = float(multi_init_r);
  float rad1 = float(multi_init_r + multi_init_dr); 
  float l_boxf = float(l_box);

  //Device memory to store randomized centers
  float4 *d_centers = 0;
  float4 *h_centers = new float4[l_nx*l_ny*l_nz];  
  CUDA_SAFE_CALL(cudaMalloc((void**)&d_centers, l_nx*l_ny*l_nz*sizeof(float4)));

  const dim3 SDTblockSize(SDTBlockDim, SDTBlockDim, SDTBlockDim);
  const dim3 SDTgridSize(iDivUp(subvolDimSDT, SDTBlockDim)*iDivUp(subvolDimSDT, SDTBlockDim)*iDivUp(subvolDimSDT, SDTBlockDim), 1);//linearized grid, since, we cannot have a 3D grid
  cudaExtent logicalGridSize = make_cudaExtent(iDivUp(subvolDimSDT, SDTBlockDim), iDivUp(subvolDimSDT, SDTBlockDim), iDivUp(subvolDimSDT, SDTBlockDim));
  char _msg[256];
  //The 1 voxel border is not initialized by our cuda kernel. Do it on the CPU :p. Outside is -ve
  for(int _x=0;_x<volSize.width; _x++) for(int _y=0; _y<volSize.height; _y++) { *(h_phi+ _x + _y*volSize.width) = maxDist; *(h_phi + _x + _y*volSize.width + (volSize.depth-1)*volSize.width*volSize.height) = maxDist; }
  for(int _x=0;_x<volSize.width; _x++) for(int _z=0; _z<volSize.depth; _z++) { *(h_phi+ _x + _z*volSize.width*volSize.height) = maxDist; *(h_phi + _x + (volSize.height-1)*volSize.width + _z*volSize.width*volSize.height) = maxDist; }
  for(int _z=0;_z<volSize.depth; _z++) for(int _y=0; _y<volSize.height; _y++) { *(h_phi + _y*volSize.width + _z*volSize.width*volSize.height) = maxDist; *(h_phi + volSize.width - 1 + _y*volSize.width + _z*volSize.width*volSize.height) = maxDist; }

  CUT_SAFE_CALL(cutStartTimer(timer));
  cudaExtent subvolIdx = make_cudaExtent(0, 0, 0);
  cudaExtent copyvol = make_cudaExtent(subvolDimSDT-2, subvolDimSDT-2, subvolDimSDT-2);
  cudaPos src_off = make_cudaPos(1, 1, 1);
  cudaPos dst_off = make_cudaPos(0, 0, 0);
  float3 offset = make_float3(l_offx, l_offy, l_offz);
  float2 radii = make_float2(rad0, rad1);
  int3 l_n = make_int3(l_nx, l_ny, l_nz);

  for (int j=0; j<nImplicits; j++) {
    h_phi = h_BUFFER[currentPhi][j];
    debugPrint(".");
    for(int k=0; k<l_nx*l_ny*l_nz; k++) { //Allocate random centers 

      h_centers[k].x = (drand48()*2. - 1.)*float(multi_init_s - 1);
      h_centers[k].y = (drand48()*2. - 1.)*float(multi_init_s - 1);
      h_centers[k].z = (drand48()*2. - 1.)*float(multi_init_s - 1);
    }
    CUDA_SAFE_CALL(cudaMemcpy((void *)d_centers, (void *)h_centers, l_nx*l_ny*l_nz*sizeof(float4), cudaMemcpyHostToDevice));
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float4>();
    CUDA_SAFE_CALL(cudaBindTexture(0, texCoord, d_centers, channelDesc)); //Make it available to the kernel
    for(int _z = 0; _z< subvolIndicesExtents.depth; _z++)
      for(int _y = 0; _y< subvolIndicesExtents.height; _y++)
        for(int _x = 0; _x< subvolIndicesExtents.width; _x++)
        {
          subvolIdx.width = _x; subvolIdx.height = _y; subvolIdx.depth = _z;
          d_init_interface_multi_smooth<<<SDTgridSize, SDTblockSize>>>(d_volPPtr, logicalGridSize, subvolIdx, l_n, maxDist, subvolDimSDT, SDTBlockDim, offset, l_boxf, radius_reach, radii);	//change from subvolDim to subvolDimSDT by zq
          CUT_CHECK_ERROR("Kernel failed");
          sprintf(_msg,  "---%d-%d-%d", _x, _y, _z);
          CHK_CUDA_ERR(_msg);
          cudaThreadSynchronize();

          //copy the results to main volume
          //Watchout for boundaries at the end. We don't want to write to volume outside allocated array!!
          adjustDnloadSubvolSize(_x, _y, _z, 2, copyvol);
          dst_off.x = 1 + _x*(subvolDimSDT-2); dst_off.y = 1 + _y*(subvolDimSDT-2); dst_off.z = 1 + _z*(subvolDimSDT-2);
          copy3DMemToHost(d_volPPtr, h_phi, copyvol, volSize, src_off, dst_off); 
          CHK_CUDA_ERR("->");
        }
  }
  CUT_SAFE_CALL(cutStopTimer(timer));
  milliseconds = cutGetAverageTimerValue(timer);
  CUT_SAFE_CALL(cutResetTimer(timer));
  delete []h_centers;
  CUDA_SAFE_CALL(cudaUnbindTexture(texCoord));    
  CUDA_SAFE_CALL(cudaFree(d_centers));  
  debugPrint(" %f ms\n", milliseconds);

  return true;
}

bool MPSegmentation::computeSDTEikonal()
{
  //Note: While using the parent method computeSDTEikonal, there is some inconsistency with binding texture to the respective cudaArray. A solution is to bind texture again in this routine in the parent class. This has been done in MSLevelSet class. -os/2010-10-10
  int l_subvolDim = subvolDim;
  subvolDim = subvolDimSDT;//Use custom subvoldim.
  debugPrint(" subvolDim is %d %d\n", subvolDim, SDTBlockDim);
  int old_currentPhi = currentPhi;





  for (int i=0; i< nImplicits; i++) {
    debugPrint("\t[%d] ", i);
    h_buffer[currentPhi] = h_BUFFER[currentPhi][i];
    h_buffer[(currentPhi+1)%2] = h_BUFFER[(currentPhi+1)%2][i];
    if(!MSLevelSet::computeSDTEikonal()) return false;  
    currentPhi = old_currentPhi; //MSLevelset::computeSDTEikonal() increments this counter. This is to reset the effect.
  }
  currentPhi = (currentPhi+1)%2; //Actually increment the counter
  subvolDim = l_subvolDim; //Restore original subvol dimension

  /*
  int x, y, z;
  float l_I;
  int nxy = volSize.width*volSize.height;

    for(int i=0; i< nelem; i++) {
    z = i/nxy;
    y = (i - z*nxy)/volSize.width; 
    x = (i - z*nxy - y*volSize.width);
	debugPrint("After compute SDTEikonal z,y,x = %d %d %d\n", z, y, x);
//    l_class = 0;
    l_I = VOXEL_AT(h_vol, x, y, z);
	debugPrint("l_I = %f \n", l_I);
  }
//*/ 


  return true;
}

bool MPSegmentation::computeAverageIntensities()
{
  //Note: As opposd to the single phi case, this routine does not calculate averages using Heaviside function.
  //Compute average intensity value inside and outside the volume
  debugPrint("\tComputing average intensities... ");
  unsigned int l_nClasses = 1<<nImplicits;
  unsigned int *l_c_count= new unsigned int[l_nClasses];
  h_PHI = h_BUFFER[currentPhi];


 


  double t0, t1;
  t0 = getTime();
  unsigned int l_class;
  float l_I = 0.0;
  for(unsigned int j = 0; j<l_nClasses; j++) { *(l_c_count+j)  = 0; c_avg[j] = 0.0; }

  unsigned long i = 0;
  unsigned long nxy = volSize.width*volSize.height;
  unsigned int x, y, z, k;
#pragma omp parallel for private(i, x, y, z, l_I, l_class, k) 
  for(i=0; i< nelem; i++) {
    z = i/nxy;
    y = (i - z*nxy)/volSize.width; 
    x = (i - z*nxy - y*volSize.width);
 //	debugPrint("z,y,x = %d %d %d\n", z, y, x);
    l_class = 0;
    l_I = VOXEL_AT(h_vol, x, y, z);
//	debugPrint("l_I = %f \n", l_I);

    for(k=0; k<nImplicits; k++)
      if(VOXEL_AT(*(h_PHI+k), x, y, z) >= 0.0f) l_class |= (0x1<<(nImplicits - 1 - k)); //Ordering of phi index in a byte (or int): MSB--> [0][1][2][3]...[nImplicits - 1] <--LSB
#pragma omp critical
    {
      c_avg[l_class] +=l_I;
      l_c_count[l_class]++;
    }
  }

//  debugPrint(" c_avg[l_class] = %f, l_c_count: = %d\n", c_avg[0],  l_c_count[0]); 
  //Divide to get actual average values...
  for(unsigned int j = 0; j<l_nClasses; j++) {
     if(l_c_count[j] == 0)
       c_avg[j] = 0.0f;
     else
       c_avg[j] /=(float)l_c_count[j];   
      //print value...
      debugPrint("%sClass[%d]<", (j%2 == 0)?"\n\t\t":"\t\t", j);
      for(unsigned int k=nImplicits; k>0; k--) debugPrint("%u", (j>>(k - 1))&0x1);
      debugPrint(", %8u> = %16f", l_c_count[j], c_avg[j]);
  }
  delete []l_c_count;
          
  t1 = getTime();
  debugPrint("\n\t\tTotal time %f sec.\n", t1 - t0);

  return true;
}

bool MPSegmentation::multiPhasePDEUpdate(bool onlyborder)
{
  float mu_factor = mu/(h*h);
  float l_vol;
  float l_extEnergy, l_Hprod, l_intEnergy;
  float l_p1, l_p0; //terms (u0 - c<>) for inside and outside c<> values for a particular region with once except all bits fixed. e.g. if the 3rd bit is varying: c_{111} and c_{011}, c_{101} and c_{001}, etc.
  unsigned int l_terms = 1<<(nImplicits-1);
  unsigned int l_bsz = sizeof(nImplicits)<<4;//number of bits in nImplicit's data type: unsigned int here. Multiply by 8
  unsigned int l_index0, l_index1;
  unsigned int l_bit;
  float l_delta_dirac;
  float *l_H_vals, *l_phi_vals, *l_phi_updates;
  
  //Update PDEs simultaneously...
  unsigned long idx = 0;
  unsigned long nxy = volSize.width*volSize.height;
  unsigned int x, y, z;
#pragma omp parallel private(idx, x, y, z, l_vol, l_H_vals, l_phi_vals, l_phi_updates, l_extEnergy, l_intEnergy, l_index0, l_index1, l_p0, l_p1, l_Hprod, l_bit, l_delta_dirac)
  {
    l_H_vals = new float[nImplicits];//cache float values for any x, y, z to avoid multiple memory reads.
    l_phi_vals = new float[nImplicits];//cache float values for any x, y, z to avoid multiple memory reads.  
    l_phi_updates = new float[nImplicits];//update values
#pragma omp for
    for(idx=0; idx< nelem; idx++) {
      z = idx/nxy;
      y = (idx - z*nxy)/volSize.width; 
      x = (idx - z*nxy - y*volSize.width);
      if(onlyborder && x > 1 && x < (volSize.width-2) && y > 1 && y < (volSize.height-2) && z > 1 && z < (volSize.depth-2)) continue; //Early reject
      l_vol = VOXEL_AT(h_vol, x, y, z);
      for(int j=0; j<nImplicits; j++)
      {
        l_phi_vals[j] = VOXEL_AT(h_BUFFER[currentPhi][j], x, y, z);
        l_H_vals[j] = 0.5*(1.0 + (2.0/PI)*atanf(l_phi_vals[j]/epsilon));
      }

      if(nImplicits == 2) {
        float t11 = l_vol - c_avg[3]; t11 *= t11;
        float t01 = l_vol - c_avg[1]; t01 *= t01;
        float t10 = l_vol - c_avg[2]; t10 *= t10;
        float t00 = l_vol - c_avg[0]; t00 *= t00;
        l_phi_updates[0] = (t11 - t01)*l_H_vals[1] + (t10 - t00)*(1. - l_H_vals[1]);
        l_phi_updates[1] = (t11 - t10)*l_H_vals[0] + (t01 - t00)*(1. - l_H_vals[0]);
        for(int i=0; i<nImplicits; i++) {
          h_phi = h_BUFFER[currentPhi][i];
          l_intEnergy = mu_factor*divergenceOfNormalizedGradient_FiniteDifference(x, y, z);
          l_delta_dirac = (epsilon/M_PI)/(epsilon*epsilon + l_phi_vals[i]*l_phi_vals[i]);	  
          l_phi_updates[i] = delta_t*l_delta_dirac*(l_intEnergy - l_phi_updates[i]); //l_phi_updates now contains the actual update to phi
          VOXEL_AT(h_BUFFER[(currentPhi+1)%2][i], x, y, z) = VOXEL_AT(h_BUFFER[currentPhi][i], x, y, z) + l_phi_updates[i];
        }
      } else if(nImplicits == 3) {
        float t111 = l_vol - c_avg[7]; t111 *= t111;
        float t011 = l_vol - c_avg[3]; t011 *= t011;
        float t110 = l_vol - c_avg[6]; t110 *= t110;
        float t010 = l_vol - c_avg[2]; t010 *= t010;
        float t101 = l_vol - c_avg[5]; t101 *= t101;
        float t001 = l_vol - c_avg[1]; t001 *= t001;
        float t100 = l_vol - c_avg[4]; t100 *= t100;
        float t000 = l_vol - c_avg[0]; t000 *= t000;
        //Store extrenal energy terms in l_phi_updates variable!
        l_phi_updates[0] = (t111 - t011)*l_H_vals[1]*l_H_vals[2] + (t110 - t010)*l_H_vals[1]*(1. - l_H_vals[2]) + (t101 - t001)*(1. - l_H_vals[1])*l_H_vals[2] + (t100 - t000)*(1. - l_H_vals[1])*(1. - l_H_vals[2]);
        l_phi_updates[1] = (t111 - t101)*l_H_vals[0]*l_H_vals[2] + (t110 - t100)*l_H_vals[0]*(1. - l_H_vals[2]) + (t011 - t001)*(1. - l_H_vals[0])*l_H_vals[2] + (t010 - t000)*(1. - l_H_vals[0])*(1. - l_H_vals[2]);
        l_phi_updates[2] = (t111 - t110)*l_H_vals[0]*l_H_vals[1] + (t101 - t100)*l_H_vals[0]*(1. - l_H_vals[1]) + (t011 - t010)*(1. - l_H_vals[0])*l_H_vals[1] + (t001 - t000)*(1. - l_H_vals[0])*(1. - l_H_vals[1]);
        for(int i=0; i<nImplicits; i++) {
          h_phi = h_BUFFER[currentPhi][i];
          l_intEnergy = mu_factor*divergenceOfNormalizedGradient_FiniteDifference(x, y, z);
          l_delta_dirac = (epsilon/M_PI)/(epsilon*epsilon + l_phi_vals[i]*l_phi_vals[i]);	  
          l_phi_updates[i] = delta_t*l_delta_dirac*(l_intEnergy - l_phi_updates[i]); //l_phi_updates now contains the actual update to phi
          VOXEL_AT(h_BUFFER[(currentPhi+1)%2][i], x, y, z) = VOXEL_AT(h_BUFFER[currentPhi][i], x, y, z) + l_phi_updates[i];
        }
      } else if(nImplicits == 4) {
        float t1111 = l_vol - c_avg[15]; t1111 *= t1111;
        float t0111 = l_vol - c_avg[7] ; t0111 *= t0111;
        float t1011 = l_vol - c_avg[11]; t1011 *= t1011;
        float t0011 = l_vol - c_avg[3] ; t0011 *= t0011;
        float t1101 = l_vol - c_avg[13]; t1101 *= t1101;
        float t0101 = l_vol - c_avg[5] ; t0101 *= t0101;
        float t1001 = l_vol - c_avg[9] ; t1001 *= t1001;
        float t0001 = l_vol - c_avg[1] ; t0001 *= t0001;

        float t1110 = l_vol - c_avg[14]; t1110 *= t1110;
        float t0110 = l_vol - c_avg[6] ; t0110 *= t0110;
        float t1010 = l_vol - c_avg[10]; t1010 *= t1010;
        float t0010 = l_vol - c_avg[2] ; t0010 *= t0010;
        float t1100 = l_vol - c_avg[12]; t1100 *= t1100;
        float t0100 = l_vol - c_avg[4] ; t0100 *= t0100;
        float t1000 = l_vol - c_avg[8] ; t1000 *= t1000;
        float t0000 = l_vol - c_avg[0] ; t0000 *= t0000;
        
        //Store extrenal energy terms in l_phi_updates variable!
	l_phi_updates[0] = (t1111 - t0111)*l_H_vals[1]*l_H_vals[2]*l_H_vals[3] + (t1110 - t0110)*l_H_vals[1]*l_H_vals[2]*(1. - l_H_vals[3]) + (t1101 - t0101)*l_H_vals[1]*(1. - l_H_vals[2])*l_H_vals[3] +
			   (t1100 - t0100)*l_H_vals[1]*(1. - l_H_vals[2])*(1. - l_H_vals[3]) + (t1011 - t0011)*(1. - l_H_vals[1])*l_H_vals[2]*l_H_vals[3] + 
			   (t1010 - t0010)*(1. - l_H_vals[1])*l_H_vals[2]*(1. - l_H_vals[3]) + (t1001 - t0001)*(1. - l_H_vals[1])*(1. - l_H_vals[2])*l_H_vals[3] +
			   (t1000 - t0000)*(1. - l_H_vals[1])*(1. - l_H_vals[2])*(1. - l_H_vals[3]);
        l_phi_updates[1] = (t1111 - t1011)*l_H_vals[0]*l_H_vals[2]*l_H_vals[3] + (t1110 - t1010)*l_H_vals[0]*l_H_vals[2]*(1. - l_H_vals[3]) + (t1101 - t1001)*l_H_vals[0]*(1. - l_H_vals[2])*l_H_vals[3] + 
			   (t1100 - t1000)*l_H_vals[0]*(1. - l_H_vals[2])*(1. - l_H_vals[3]) + (t0111 - t0011)*(1. - l_H_vals[0])*l_H_vals[2]*l_H_vals[3] +
			   (t0110 - t0010)*(1. - l_H_vals[0])*l_H_vals[2]*(1. - l_H_vals[3]) + (t0101 - t0001)*(1. - l_H_vals[0])*(1. - l_H_vals[2])*l_H_vals[3] + 
			   (t0100 - t0000)*(1. - l_H_vals[0])*(1. - l_H_vals[2])*(1. - l_H_vals[3]);  
        l_phi_updates[2] = (t1111 - t1101)*l_H_vals[0]*l_H_vals[1]*l_H_vals[3] + (t1110 - t1100)*l_H_vals[0]*l_H_vals[1]*(1. - l_H_vals[3]) + (t1011 - t1001)*l_H_vals[0]*(1. - l_H_vals[1])*l_H_vals[3] + 
			   (t1010 - t1000)*l_H_vals[0]*(1. - l_H_vals[1])*(1. - l_H_vals[3]) + (t0111 - t0101)*(1. - l_H_vals[0])*l_H_vals[1]*l_H_vals[3] + 
			   (t0110 - t0100)*(1. - l_H_vals[0])*l_H_vals[1]*(1. - l_H_vals[3]) + (t0011 - t0001)*(1. - l_H_vals[0])*(1. - l_H_vals[1])*l_H_vals[3] + 
			   (t0010 - t0000)*(1. - l_H_vals[0])*(1. - l_H_vals[1])*(1. - l_H_vals[3]);
	l_phi_updates[3] = (t1111 - t1110)*l_H_vals[0]*l_H_vals[1]*l_H_vals[2] + (t1101 - t1100)*l_H_vals[0]*l_H_vals[1]*(1. - l_H_vals[2]) + (t1011 - t1010)*l_H_vals[0]*(1. - l_H_vals[1])*l_H_vals[2] + 
			   (t1001 - t1000)*l_H_vals[0]*(1. - l_H_vals[1])*(1. - l_H_vals[2]) + (t0111 - t0110)*(1. - l_H_vals[0])*l_H_vals[1]*l_H_vals[2] + 
			   (t0101 - t0100)*(1. - l_H_vals[0])*l_H_vals[1]*(1. - l_H_vals[2]) + (t0011 - t0010)*(1. - l_H_vals[0])*(1. - l_H_vals[1])*l_H_vals[2] + 
			   (t0001 - t0000)*(1. - l_H_vals[0])*(1. - l_H_vals[1])*(1. - l_H_vals[2]);
        for(int i=0; i<nImplicits; i++) {
          h_phi = h_BUFFER[currentPhi][i];
          l_intEnergy = mu_factor*divergenceOfNormalizedGradient_FiniteDifference(x, y, z);
          l_delta_dirac = (epsilon/M_PI)/(epsilon*epsilon + l_phi_vals[i]*l_phi_vals[i]);	  
          l_phi_updates[i] = delta_t*l_delta_dirac*(l_intEnergy - l_phi_updates[i]); //l_phi_updates now contains the actual update to phi
          VOXEL_AT(h_BUFFER[(currentPhi+1)%2][i], x, y, z) = VOXEL_AT(h_BUFFER[currentPhi][i], x, y, z) + l_phi_updates[i];
        }
      } else {
        //Solve PDE
        for(int i=0; i<nImplicits; i++)
        {
          // Compute the external energy term. This term has 2^{nImplicits -1} terms
          l_extEnergy = 0.0;

          for(unsigned int k=0; k<l_terms; k++)//This loop enumerates possible combinations of (nImplicits-1) bits
          {
            l_index0 = ((k>>i)<<(i+1)) | ((i==0)?0:((k<<(l_bsz-i))>>(l_bsz-i)));//insert 0 at ith place
            l_index1 = ((k>>i)<<(i+1)) | (0x1<<i) | ((i==0)?0:((k<<(l_bsz-i))>>(l_bsz-i)));//insert 1 at ith place
            l_p0 = l_vol - c_avg[l_index0];
            l_p1 = l_vol - c_avg[l_index1];	    
            l_p0 *= l_p0;
            l_p1 *= l_p1;
            l_Hprod = 1.0;
            for (unsigned int j=0; j<nImplicits; j++)//bit position
            {
              if(j==i) continue; //There is no multiplier for H_i
              l_bit = (l_index0>>j)&0x1;//Get the jth bit
              l_Hprod *= (float)l_bit + (1.0 - 2.0*float(l_bit))*l_H_vals[j];///check -os 2010-10-11
            }
            l_extEnergy += (l_p1 - l_p0)*l_Hprod;
          }

          // Compute the external energy term.	  
          h_phi = h_BUFFER[currentPhi][i];
          l_intEnergy = mu_factor*divergenceOfNormalizedGradient_FiniteDifference(x, y, z);

          // Compute update
          l_delta_dirac = (epsilon/M_PI)/(epsilon*epsilon + l_phi_vals[i]*l_phi_vals[i]);	  
          l_phi_updates[i] = delta_t*l_delta_dirac*(l_intEnergy - l_extEnergy);
          VOXEL_AT(h_BUFFER[(currentPhi+1)%2][i], x, y, z) = VOXEL_AT(h_BUFFER[currentPhi][i], x, y, z) + l_phi_updates[i];
        }
      }
    }
    delete []l_phi_vals; delete []l_H_vals; delete []l_phi_updates;
  }
  return true;
}

// The bordering scheme changes here since we need 1 voxel border around subvol such that part of the actual
// volume boundary has a border with zero values. Thus an effectvive border of zeros around the full volume.
// This is needed only for the coeff texture. For the other two textures we don't care since computation is dependent 
// only on the current (x, y, z). Effective computational volume: (subvolDim-2)^3
bool MPSegmentation::PDEUpdate()
{
  /// Update phi by single time step.
  ///Compute the higher and lower order terms.
  debugPrint("\tUpdating the PDE (level set and input volume)");
#if UPDATE_MULTIPHASE_PDE_ON_CPU 
  debugPrint("[on CPU]... ");
  double t0, t1;
  t0 = getTime();
  multiPhasePDEUpdate(false);
  t1 = getTime();
  milliseconds = (t1 - t0)*1000.0;
#else
  int l_subvolDim = subvolDim;
  subvolDim = subvolDimPDE;//Use custom subvol size.
  debugPrint("[on GPU]... ");
  float mu_factor = mu/(h*h);
  const dim3 PDEblockSize(PDEBlockDim, PDEBlockDim, PDEBlockDim);
  const dim3 PDEgridSize(iDivUp(subvolDim, PDEBlockDim)*iDivUp(subvolDim, PDEBlockDim)*iDivUp(subvolDim, PDEBlockDim), 1);//linearized grid, since, we cannot have a 3D grid
  cudaExtent logicalGridSize = make_cudaExtent(iDivUp(subvolDim, PDEBlockDim), iDivUp(subvolDim, PDEBlockDim), iDivUp(subvolDim, PDEBlockDim));
  //Upload average values to the GPU global mem
  CUDA_SAFE_CALL( cudaMemcpy( (void*)d_c_avg, (void*)c_avg, (1<<nImplicits)*sizeof(float), cudaMemcpyHostToDevice));
  
  cudaExtent copyvol = make_cudaExtent(subvolDim-4, subvolDim-4, subvolDim-4);
  cudaExtent PDEsubvolIndicesExtents = make_cudaExtent(iDivUp(volSize.width, subvolDim-4), iDivUp(volSize.height, subvolDim-4), iDivUp(volSize.depth, subvolDim-4));  
  cudaPos src_off = make_cudaPos(2, 2, 2);
  
  cudaPos dst_off = make_cudaPos(0, 0, 0);
  char _msg[256];
  float maxDist = sqrtf(float(volSize.width*volSize.width + volSize.height*volSize.height + volSize.depth*volSize.depth));
  
  CUT_SAFE_CALL(cutStartTimer(timer));  
  cudaExtent copyvol_upload = make_cudaExtent(0,0,0);
  cudaPos offset_upload = make_cudaPos(0, 0, 0);

  //We copy the boundary to the texture, so the extra volume uploaded must have the same value as the outside of the interface initialization.
  for(int i=0; i< subvolDim*subvolDim*subvolDim; i++)      
    *(h_subvolSpare + i) = -maxDist;	
      
  for(int _z = 0; _z< PDEsubvolIndicesExtents.depth; _z++)
    for(int _y = 0; _y< PDEsubvolIndicesExtents.height; _y++)
      for(int _x = 0; _x< PDEsubvolIndicesExtents.width; _x++)
      {
        adjustUploadSubvolSize(_x, _y, _z, 4, copyvol_upload);
        //load texture values for coeff, phi and intensity values
        // Watchout for boundaries at the end. We don't want to read garbage values outside the allocated host volume!!
        if(zeroOutArray)
        {
          //zero out Arrays first so they do not have undesired values stored in them from previous load
          //N.B.: there seems to be no API to memset a cudaArray, so we do it indirectly
          cudaMemset3D(d_volPPtr, 0, make_cudaExtent(subvolDim*sizeof(float), subvolDim, subvolDim));// zero out the device mem
          copy3DMemToArray(d_volPPtr, d_volArray); 
          //Fill the phi texture with -maxDist so that it doesn't have undesirable values when bundary volume parts are copied.
          for(unsigned int i=0; i<nImplicits; i++)
          {
            int __x = i%nImplicitsX;
            int __y = (i/nImplicitsX)%nImplicitsY;
            int __z = i/(nImplicitsX*nImplicitsY);
            copy3DHostToArray(h_subvolSpare, d_phiArray, subvolSize, subvolSize, make_cudaPos(0, 0, 0), make_cudaPos(__x*subvolDim, __y*subvolDim, __z*subvolDim)); 
          }	              
        }

        for(unsigned int i=0; i<nImplicits; i++)// Load all phi functions to GPU in subvolume size
        {
          //Available blocks in GPU mem is: nImplicitsZ*nImplicitsY*nImplicitsX >=nImplicits. Just fill it linearly. X first , then Y and then Z.
          int __x = i%nImplicitsX;
          int __y = (i/nImplicitsX)%nImplicitsY;
          int __z = i/(nImplicitsX*nImplicitsY);
          copy3DHostToArray(h_BUFFER[currentPhi][i], d_phiArray, copyvol_upload, volSize, make_cudaPos(_x*(subvolDim-4), _y*(subvolDim-4), _z*(subvolDim-4)), make_cudaPos(__x*subvolDim, __y*subvolDim, __z*subvolDim));
          CHK_CUDA_ERR("->");
        }
        //Upload volume
        copy3DHostToArray(h_vol, d_volArray, copyvol_upload, volSize, make_cudaPos(_x*(subvolDim-4), _y*(subvolDim-4), _z*(subvolDim-4)));
        CHK_CUDA_ERR("->");
		
        //update all the PDE's 
        int3 texPhi_partition;
        for(unsigned int i=0; i<nImplicits; i++)
        {
          //PDE evolution
          texPhi_partition.x = nImplicitsX;
          texPhi_partition.y = nImplicitsY;
          texPhi_partition.z = nImplicitsZ;
          MultiPhase_PDE_update<<<PDEgridSize, PDEblockSize>>>(d_volPPtr, nImplicits, texPhi_partition, i, logicalGridSize, mu_factor, epsilon, delta_t, subvolDim, PDEBlockDim);
          CUT_CHECK_ERROR("Kernel failed");
          sprintf(_msg,  "---%d-%d-%d", _x, _y, _z);
          CHK_CUDA_ERR(_msg);     	
          cudaThreadSynchronize();	  	
	  
          // Copy results back to h_phi. 
          // Watchout for boundaries at the end. We don't want to write to volume outside allocated array!!
          adjustDnloadSubvolSize(_x, _y, _z, 4, copyvol);
          dst_off.x = 2+ _x*(subvolDim-4); dst_off.y = 2 + _y*(subvolDim-4); dst_off.z = 2 + _z*(subvolDim-4);
          copy3DMemToHost(d_volPPtr, h_BUFFER[(currentPhi+1)%2][i], copyvol, volSize, src_off, dst_off); 
          CHK_CUDA_ERR("<-%s");	  
        }
      }  

  CUT_SAFE_CALL(cutStopTimer(timer));
  milliseconds = cutGetAverageTimerValue(timer);
  CUT_SAFE_CALL(cutResetTimer(timer));

  //Compute border correctly on CPU
  double t0, t1;
  t0 = getTime();
  multiPhasePDEUpdate(true);
  t1 = getTime();
  milliseconds += (t1 - t0)*1000.0;
  subvolDim = l_subvolDim; //Restore the original subvol dimension.
#endif
  currentPhi = (currentPhi+1)%2;
  h_PHI = h_BUFFER[currentPhi];
  debugPrint("%f ms\n", milliseconds);  

  return true;
}

//Problem:
// Given an integer i > 0, find integers (optimal factors of i) a, b, c > 0, <= n, such that
// Q = min(x*y*z - i)
bool MPSegmentation::computeOptimalFactors(unsigned int i, unsigned int n, unsigned int &a, unsigned int &b, unsigned int &c) //Solves an integer programming problem.
{
  unsigned int Q;
  bool _2dSearchFailed = false;
  bool _3dSearchFailed = false;
  multimap <unsigned int, Factors> possibilities;

  if(i > n*n*n)
  {
    debugPrint("computeOptimalFactors: i cannot be more than n^3!\n");
    return false;
  }

  //In this optimization, we'll prefer filling up the X dimension first, then Y and then Z. This is just a heuristic.
  if(i<=n)
  {
    a = i;
    b = 1;
    c = 1;
    Q = 0;
  }
  else if(i<=n*n) //planar optimization
  {
    c = 1;
    //debugPrint("Performing optimal search in 2D...\n");
    int _a, _b;
    for(_a = n; _a>=1; _a--)
    {
      _b = (i%_a == 0)? (i/_a) : (i/_a + 1);
      if (_b > n)
      {
        if(_a == 1) _2dSearchFailed = true;
        continue;
      }
      Q = _a*_b - i;
      //debugPrint("a = %2u, b = %2u, c = %2u; Q = %3u\n", _a, _b, z, Q);
      if(Q == 0) break; 
      else {
        Factors dec = {_a, _b, c};
        possibilities.insert(pair<unsigned int, Factors>(Q,dec));
      }
    }
   
    if(Q > 0) _2dSearchFailed = true;
    else {
      a = _a;
      b = _b;
    }
  }

  if(_2dSearchFailed || (i>n*n))
  {
    //if(_2dSearchFailed) debugPrint("2D search failed, doing a 3D search...\n");
    //else debugPrint("Performing a 3D search...\n");
    unsigned _a, _b, _c, _ab;
    for(_b = n; _b>=1; _b--)
      for(_a = n; _a>=1; _a--)
      {
        _ab = _a*_b;
        if(_2dSearchFailed && _ab > n) continue;
        _c = (i%_ab == 0)? (i/_ab) : (i/_ab + 1);
        if (_c > n)
          continue;
        Q = _ab*_c - i;
        //debugPrint("a = %2u, b = %2u, c = %2u; Q = %3u\n", _a, _b, _c, Q);
        if(Q == 0) break; 
        else {
          Factors dec = {_a, _b, _c};
          possibilities.insert(pair<unsigned int, Factors>(Q,dec));
        }
      }
    if(Q > 0) _3dSearchFailed = true;
    else {
      a = _a;
      b = _b;
      c = _c;
    }
  }
  if(_2dSearchFailed || _3dSearchFailed) {
    //debugPrint("Q = 0 is not possible. Searching next minimum Q...\n");
    multimap<unsigned int, Factors>::iterator p;
    p = possibilities.begin();
    Q = p->first;
    a = p->second.x;
    b = p->second.y;
    c = p->second.z;
  }
  //debugPrint("computeOptimalFactors: Optimal solution: %2u = %2u x %2u x %2u; Q = %3u\n", i, a, b, c, Q);
  possibilities.clear();
  return true;
}

void MPSegmentation::copy3DHostToArray(float *_src, cudaArray *_dst, cudaExtent copy_extent, cudaExtent src_extent, cudaPos src_offset, cudaPos dst_offset)
{
  cudaMemcpy3DParms copyParams = {0};
  float *h_source = _src + src_offset.x + src_offset.y*volSize.width + src_offset.z*volSize.width*volSize.height;
  copyParams.srcPtr = make_cudaPitchedPtr((void*)h_source, src_extent.width*sizeof(float), src_extent.width, src_extent.height);
  copyParams.dstArray = _dst;
  copyParams.kind = cudaMemcpyHostToDevice;
  copyParams.extent = copy_extent;
  copyParams.dstPos = dst_offset;
  CUDA_SAFE_CALL(cudaMemcpy3D(&copyParams));
  CUT_CHECK_ERROR("Host -> Array Memcpy failed\n");
}
    
void MPSegmentation::setParameters(int _width, int _height, int _depth, MPLevelSetParams *MPLSParams)
{
  MSLevelSet::setParameters(_width, _height, _depth, dynamic_cast<MSLevelSetParams*>(MPLSParams));
  nImplicits = MPLSParams->nImplicits;
  debugPrint("---Multi-phase parameters---\n");
  debugPrint("Specified nImplicits: %d\n", nImplicits);
  subvolDim = MPLSParams->subvolDim;
  subvolDimSDT = MPLSParams->subvolDimSDT;
  subvolDimAvg = MPLSParams->subvolDimAvg;
  subvolDimPDE = MPLSParams->subvolDimPDE;
  debugPrint("subvolDim: SDT = %d, Avg = %d, PDE = %d\n", subvolDimSDT, subvolDimAvg, subvolDimPDE);
}

float MPSegmentation::divergenceOfNormalizedGradient_FiniteDifference( int x,  int y,  int z)
{
  float grad_left[3], grad_right[3], div;
  float l_phi = VOXEL_AT(h_phi, x, y, z);
  float grad_norm;
  if(x > 1 && x < (volSize.width-2)) {
	  grad_left[0] = (l_phi - VOXEL_AT(h_phi, x - 2, y, z))/2.0;
	  grad_right[0] = (VOXEL_AT(h_phi, x + 2, y, z) - l_phi)/2.0;
  } else if (x == 1) {
	  grad_left[0] = l_phi - VOXEL_AT(h_phi, x - 1, y, z);//One sided difference
	  grad_right[0] = (VOXEL_AT(h_phi, x + 2, y, z) - l_phi)/2.0;
  } else if(x == (volSize.width-2)) {
	  grad_left[0] = (l_phi - VOXEL_AT(h_phi, x - 2, y, z))/2.0;
	  grad_right[0] = VOXEL_AT(h_phi, x+1, y, z) - l_phi;
  } else if(x == 0) {
	  grad_left[0] = 0.;
	  grad_right[0] = (VOXEL_AT(h_phi, x + 2, y, z) - l_phi)/2.0;
  } else {
	  grad_left[0] = (l_phi - VOXEL_AT(h_phi, x - 2, y, z))/2.0;
	  grad_right[0] = 0.;
  }

  if(y > 1 && y < (volSize.height-2)) {
	  grad_left[1] = (l_phi - VOXEL_AT(h_phi, x, y - 2, z))/2.0;
	  grad_right[1] = (VOXEL_AT(h_phi, x, y + 2, z) - l_phi)/2.0;
  } else if (y == 1) {
	  grad_left[1] = l_phi - VOXEL_AT(h_phi, x, y - 1, z);//One sided difference
	  grad_right[1] = (VOXEL_AT(h_phi, x, y + 2, z) - l_phi)/2.0;
  } else if(y == (volSize.height-2)) {
	  grad_left[1] = (l_phi - VOXEL_AT(h_phi, x, y - 2, z))/2.0;
	  grad_right[1] = VOXEL_AT(h_phi, x, y + 1, z) - l_phi;
  } else if(y == 0) {
	  grad_left[1] = 0.;
	  grad_right[1] = (VOXEL_AT(h_phi, x, y + 2, z) - l_phi)/2.0;
  } else {
	  grad_left[1] = (l_phi - VOXEL_AT(h_phi, x, y - 2, z))/2.0;
	  grad_right[1] = 0.;
  }

  if(z > 1 && z < (volSize.depth-2)) {
	  grad_left[2] = (l_phi - VOXEL_AT(h_phi, x, y, z - 2))/2.0;
	  grad_right[2] = (VOXEL_AT(h_phi, x, y, z + 2) - l_phi)/2.0;
  } else if (z == 1) {
	  grad_left[2] = l_phi - VOXEL_AT(h_phi, x, y, z - 1);//One sided difference
	  grad_right[2] = (VOXEL_AT(h_phi, x, y, z + 2) - l_phi)/2.0;
  } else if(z == (volSize.depth-2)) {
	  grad_left[2] = (l_phi - VOXEL_AT(h_phi, x, y, z - 2))/2.0;
	  grad_right[2] = VOXEL_AT(h_phi, x, y, z + 1) - l_phi;
  } else if(z == 0) {
	  grad_left[2] = 0.;
	  grad_right[2] = (VOXEL_AT(h_phi, x, y, z + 2) - l_phi)/2.0;
  } else {
	  grad_left[2] = (l_phi - VOXEL_AT(h_phi, x, y, z - 2))/2.0;
	  grad_right[2] = 0.;
  }

  grad_norm = sqrt(grad_left[0]*grad_left[0] + grad_left[1]*grad_left[1] + grad_left[2]*grad_left[2] + TINY);
  grad_left[0] /= grad_norm; grad_left[1] /= grad_norm; grad_left[2] /= grad_norm;

  grad_norm = sqrt(grad_right[0]*grad_right[0] + grad_right[1]*grad_right[1] + grad_right[2]*grad_right[2] + TINY);
  grad_right[0] /= grad_norm; grad_right[1] /= grad_norm; grad_right[2] /= grad_norm;

  div= (grad_right[0] - grad_left[0] + grad_right[1] - grad_left[1] + grad_right[2] - grad_left[2])/2.0;
  return div;
}
