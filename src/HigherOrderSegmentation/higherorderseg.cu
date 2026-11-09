//============================================================================
/* Copyright (c) The Technical University of Denmark
 * Developed at the Computational Visualization Center (CVC), The University of Texas at Austin
 * Author: Ojaswa Sharma
 * E-mail: os@imm.dtu.dk
 * File: higherorderseg.cu
 * Main higher order levelset segmentation implementation and host functions to the CUDA device functions
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

#include <HigherOrderSegmentation/higherorderseg.h>
#include <HigherOrderSegmentation/private/higherorderseg_global.cu.h>
#include <HigherOrderSegmentation/private/cubic_coefficients_kernels.cu.h>
#include <HigherOrderSegmentation/private/pde_update_kernels.cu.h>

#include <cutil_inline.h>

using namespace HigherOrderSegmentation;

HOSegmentation::HOSegmentation(bool doInitCUDA): MSLevelSet(doInitCUDA)
{
  CCblockDim = 512;
  d_coeffArray = NULL;
  h_coeff = NULL;
}

bool HOSegmentation::runSolver(float *vol, int _width, int _height, int _depth, MSLevelSetParams *MSLSParams,
               		   void (*evolutionCallback)(const float *vol, int dimx, int dimy, int dimz))
{
  h_vol = vol;
  if(!h_vol) return false;
  
  setParameters (_width, _height, _depth, MSLSParams);

  bool result = solverMain();
  return result;
}

int HOSegmentation::solverMain()//entry point into the CUDA solver. 
{
  if(DoInitializeCUDA)
  {
    cuInit(0); //Any driver API can also be called now on...
    cudaSetDevice(cutGetMaxGflopsDeviceId());
  }
  CUT_SAFE_CALL(cutCreateTimer(&timer));
  CUT_SAFE_CALL(cutResetTimer(timer));
  
  // set the volume size and number of subvolumes
  volSize = make_cudaExtent(datainfo.n[0], datainfo.n[1], datainfo.n[2]);
  subvolIndicesExtents = make_cudaExtent(iDivUp(volSize.width-2, subvolDim-2), iDivUp(volSize.height-2, subvolDim-2), iDivUp(volSize.depth-2, subvolDim-2));
  subvolSize = make_cudaExtent(subvolDim, subvolDim, subvolDim);

  if(!initCuda()) return false;//Perhaps memory allcation failed!
  
  // Allocate two buffers for host memory for phi (Signed Distance function)
  h_buffer[0] = (float*)malloc(volSize.width*volSize.height*volSize.depth*sizeof(float));
  if(h_buffer[0] == NULL)
  {
    debugPrint("Failed to allocate %ld bytes of memory\n", volSize.width*volSize.height*volSize.depth*sizeof(float));
    freeHostMem();
    return false;
  } 
  h_buffer[1] = (float*)malloc(volSize.width*volSize.height*volSize.depth*sizeof(float));
  if(h_buffer[1] == NULL)
  {
    debugPrint("Failed to allocate %ld bytes of memory\n", volSize.width*volSize.height*volSize.depth*sizeof(float));
    freeHostMem();
    return false;
  } 
  h_phi = h_buffer[0];
  currentPhi = 0;

  //Allocate memory for the spare subvolume on the host --needed by computations by the parent class.
  h_subvolSpare = (float*)malloc((subvolDim)*(subvolDim)*(subvolDim)*sizeof(float));
  if(h_subvolSpare == NULL)
  {
    debugPrint("Failed to allocate %ld bytes of memory\n", (subvolDim)*(subvolDim)*(subvolDim)*sizeof(float));
    freeHostMem();
    return false;
  } 
  
  // Allocate host memory for cubic coefficients
  h_coeff = (float*)malloc(volSize.width*volSize.height*volSize.depth*sizeof(float));
  if(h_coeff == NULL)
  {
    debugPrint("Failed to allocate %ld bytes of memory\n", volSize.width*volSize.height*volSize.depth*sizeof(float));
    freeHostMem();
    return false;
  }
  
  // Allocate cuda array and global mem for subvol. Bind texture to cuda array
  debugPrint("Total volume: %zdx%zdx%zd\n", volSize.width, volSize.height, volSize.depth);
  debugPrint("Total input volume: %udx%udx%ud\n", datainfo.n_input[0], datainfo.n_input[1], datainfo.n_input[2]);
  debugPrint("Total subvolumes: %zdx%zdx%zd\n", subvolIndicesExtents.width, subvolIndicesExtents.height, subvolIndicesExtents.depth);

  //scale the input volume
  normalizeVolume();

  initInterfaceMulti();// Initialize level set interface for the phi function
  
  //Compute 9 derivatives  (fx, fy, fz, fxx, fyy, fzz, fxy, fyz, fzx)  at 27 points of the cubic spline function to compute higher order 
  //term in PDEUpdate. This should be done only once.
  float *ptr;
  float _x, _y, _z;
  for(int k=0; k<3; k++)
    for(int j=0; j<3; j++)
      for(int i=0; i<3; i++)
      {
       ptr = h_cubicDerivatives + (i+3*j+9*k)*9;
        _x = (float)i - 1.0; _y = (float)j - 1.0; _z = (float)k - 1.0;
        *(ptr)   = evalTriCubic_Dx(_x, _y, _z);
        *(ptr+1) = evalTriCubic_Dy(_x, _y, _z);
        *(ptr+2) = evalTriCubic_Dz(_x, _y, _z);
        *(ptr+3) = evalTriCubic_Dxx(_x, _y, _z);
        *(ptr+4) = evalTriCubic_Dyy(_x, _y, _z);
        *(ptr+5) = evalTriCubic_Dzz(_x, _y, _z);
        *(ptr+6) = evalTriCubic_Dxy(_x, _y, _z);
        *(ptr+7) = evalTriCubic_Dyz(_x, _y, _z);
        *(ptr+8) = evalTriCubic_Dzx(_x, _y, _z);
      }

  //transfer this to the GPU constant memory.
  CUDA_SAFE_CALL(cudaMemcpyToSymbol(d_cubicDerivatives, h_cubicDerivatives, 27*9*sizeof(float), 0, cudaMemcpyHostToDevice));
  CHK_CUDA_ERR("<-Constant mem.");   
    
  if(nIter > 0)
  {
    debugPrint("Solving the PDE...\n");
    for (int i = 0; i<nIter; i++)
      if(!solve()) return false;
  }
   
  //Copy h_phi to h_vol so that volRover can display it.
  memcpy(h_vol, h_phi, volSize.width*volSize.height*volSize.depth*sizeof(float));
  //memcpy(h_vol, h_coeff, volSize.width*volSize.height*volSize.depth*sizeof(float));
  cleanUp();
  return true;
}

bool HOSegmentation::solve()
{
  //Solve the PDE
  if(iter_num < nIter)
  {
    iter_num++;
    debugPrint("Iteration: %d\n\t", iter_num);
    if(!computeSDTEikonal()) return false; // re-initialize phi to a signed distance function
    debugPrint("\t");
    if(!computeAverageIntensities()) return false;
    debugPrint("\t");
    if(!computeCubicCoefficients()) return false; //Compute cubic coefficients for phi.
    debugPrint("\t");
    if(!PDEUpdate()) return false; //update phi and vol via level set PDE's
  }
  debugPrint("Solver done.\n");
  return true;

}

void HOSegmentation::cleanUp()
{
  freeHostMem();
  //Free GPU memory (if any allocated within this class)
  MSLevelSet::cleanUp();
}

void HOSegmentation::freeHostMem()
{
  //Calling parent class's cleanup trues to free h_phi again. Beware of this. Setting h_phi to NULL comes to rescue. (TODO: Improve the memory management on host sides between two classes.)
  if(h_buffer[0]) {free(h_buffer[0]); h_buffer[0] = NULL;}
  if(h_buffer[1]) {free(h_buffer[1]); h_buffer[1] = NULL;}
  if(h_subvolSpare) { free(h_subvolSpare); h_subvolSpare = NULL;}
}

bool HOSegmentation::initCuda()
{
  unsigned int freemem, totalmem, peakrequiredmem;
  if(!getFreeGPUMem(freemem, totalmem)) return false;
  unsigned int subvolmem = subvolSize.width*subvolSize.height*subvolSize.depth*sizeof(float);
  unsigned int max_vol_plane_mem_req_after_freeup = MAX(volSize.width*volSize.height, MAX(volSize.height*volSize.depth, volSize.depth*volSize.width))*sizeof(float) - 3*subvolmem;// Might be -ve, but okay!
  unsigned int max_req_coeff_mem = MAX(max_vol_plane_mem_req_after_freeup, subvolmem);
  peakrequiredmem = subvolmem*6 + max_req_coeff_mem;//d_phi, d_vol, d_coeff, spare1, 2, 3, d_volPPtr
  if(peakrequiredmem > freemem) {
    debugPrint("On GPU peak required memory (%u bytes) is more than the available (%u bytes). Please free some GPU memory (perhaps by a video mode switch.)\n", peakrequiredmem, freemem);
    return false;
  }
  bool retval =  MSLevelSet::initCuda();
  
  texCoeff.normalized = false; //Set correct texture mode.
  texCoeff.filterMode = cudaFilterModePoint;
  texCoeff.addressMode[0] = cudaAddressModeWrap; 
  texCoeff.addressMode[1] = cudaAddressModeWrap; 
  texCoeff.addressMode[2] = cudaAddressModeWrap; 

  return retval;
}

bool HOSegmentation::getFreeGPUMem(unsigned int &free, unsigned int &total)
{
  //Figure out how much GPU memory is remaining...
  CUresult res;
  CUdevice dev;
  CUcontext ctx;
  cuDeviceGet(&dev,cudaDevice);
  cuCtxCreate(&ctx, 0, dev);
  size_t _free, _total;
  res = cuMemGetInfo(&_free, &_total);
  free = _free;
  total = _total;
  if(res != CUDA_SUCCESS) {
    debugPrint("Failed to get results from cuMemGetInfo (status = %x)", res);
    return false;
  }
  cuCtxDetach(ctx);
  return true;
}

//The component in gpuSlice set to 0 is determined as per available memory
bool HOSegmentation::getAvailbleGPUSliceSize(cudaExtent &gpuSlice)
{
  if(gpuSlice.width != 0 && gpuSlice.height != 0 && gpuSlice.depth != 0) {
    debugPrint("At least one dimension must have 0 size in input argument while requesting GPU slice size.\n");
    return false;
  }
  unsigned int free, total;
  if(!getFreeGPUMem(free, total)) return false;
  
  //Figure out optimal volume slice size. say, Width x Height x l_slice, if gpuSlice.depth == 0
  int plane_size;
  int l_size;
  
  if(gpuSlice.width == 0) {
    plane_size = gpuSlice.height*gpuSlice.depth;
    l_size =  free/(2*plane_size*sizeof(float));
    if(l_size > volSize.width) l_size = volSize.width;
    gpuSlice.width = l_size;
  }
  else if(gpuSlice.height == 0) {
    plane_size = gpuSlice.width*gpuSlice.depth;
    l_size =  free/(2*plane_size*sizeof(float));
    if(l_size > volSize.height) l_size = volSize.height;
    gpuSlice.height = l_size;
  }  
  else if(gpuSlice.depth == 0) {
    plane_size = gpuSlice.height*gpuSlice.width;
    l_size =  free/(2*plane_size*sizeof(float));
    if(l_size > volSize.depth) l_size = volSize.depth;    
    gpuSlice.depth = l_size;
  }
    
  if(l_size < 1) //consider the case when exisiting cuda mem needs to be freed! We need at least once slice to process.
  {
    debugPrint("Not enought GPU memory for cubic coefficient compuation.\n");
    return false;
  }  
  return true;
}

bool HOSegmentation::allocateSliceOnGPU(cudaArray*& d_coeffArray, cudaPitchedPtr &d_coeffPPtr, cudaExtent coeffSubvolSize)
{
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
  //Allocate memory and bind coeffArray to 3D texture
  CUDA_SAFE_CALL(cudaMalloc3DArray(&d_coeffArray, &channelDesc, coeffSubvolSize));
  CHK_CUDA_ERR("\n");
  //Allocate memory for 3D device memory for cubic coefficient output
  cudaExtent pitchedVolSize = make_cudaExtent(coeffSubvolSize.width*sizeof(float), coeffSubvolSize.height, coeffSubvolSize.depth); 
  CUDA_SAFE_CALL(cudaMalloc3D(&d_coeffPPtr, pitchedVolSize));
  CHK_CUDA_ERR("\n");
  
  return true;
}

float HOSegmentation::evaluateCubicSplineAtGridPoint(int x, int y, int z)
{
  float yz_matrix[9], z_vector[3], result;

  if((x > 0 && x< (volSize.width-1)) && (y>0 && y<(volSize.height-1)) && (z>0 && z<(volSize.depth-1)))
    for(int k=-1; k<2; k++)
      for(int j=-1; j<2; j++)
        yz_matrix[j + 1 + (k + 1)*3] = (VOXEL_AT(h_coeff, x-1, y+j, z+k) + 4.0*VOXEL_AT(h_coeff, x, y+j, z+k) + VOXEL_AT(h_coeff, x+1, y+j, z+k))/6.0;
  else
    for(int k=-1; k<2; k++)
      for(int j=-1; j<2; j++)
        yz_matrix[j + 1 + (k + 1)*3] = (SAFE_VOXEL_AT(h_coeff, x-1, y+j, z+k) + 4.0*SAFE_VOXEL_AT(h_coeff, x, y+j, z+k) + SAFE_VOXEL_AT(h_coeff, x+1, y+j, z+k))/6.0;   

  for(int k=-1; k<2; k++)
    z_vector[k+1] = (yz_matrix[(k+1)*3] + 4.0*yz_matrix[1 + (k+1)*3] + yz_matrix[2 + (k+1)*3])/6.0;
    
  result = (z_vector[0] + 4.0*z_vector[1] + z_vector[2])/6.0;        
    
  return result;
}

float HOSegmentation::evaluateCubicSplineGenericDerivativeAtGridPoint(int x, int y, int z, float *xGridSplineValues, float *yGridSplineValues, float *zGridSplineValues)// GridSPlineValues are values (may be derivatives) at grid points. Is is vector of 3 floats.
{
  float yz_matrix[9], z_vector[3], result;

  if((x > 0 && x< (volSize.width-1)) && (y>0 && y<(volSize.height-1)) && (z>0 && z<(volSize.depth-1)))
    for(int k=-1; k<2; k++)
      for(int j=-1; j<2; j++)
        yz_matrix[j + 1 + (k + 1)*3] = xGridSplineValues[0]*VOXEL_AT(h_coeff, x-1, y+j, z+k) + xGridSplineValues[1]*VOXEL_AT(h_coeff, x, y+j, z+k) + xGridSplineValues[2]*VOXEL_AT(h_coeff, x+1, y+j, z+k);
  else
    for(int k=-1; k<2; k++)
      for(int j=-1; j<2; j++)
        yz_matrix[j + 1 + (k + 1)*3] = xGridSplineValues[0]*SAFE_VOXEL_AT(h_coeff, x-1, y+j, z+k) + xGridSplineValues[1]*SAFE_VOXEL_AT(h_coeff, x, y+j, z+k) + xGridSplineValues[2]*SAFE_VOXEL_AT(h_coeff, x+1, y+j, z+k);   

  for(int k=-1; k<2; k++)
    z_vector[k+1] = yGridSplineValues[0]*yz_matrix[(k+1)*3] + yGridSplineValues[1]*yz_matrix[1 + (k+1)*3] + yGridSplineValues[2]*yz_matrix[2 + (k+1)*3];
    
  result = zGridSplineValues[0]*z_vector[0] + zGridSplineValues[1]*z_vector[1] + zGridSplineValues[2]*z_vector[2];        
    
  return result;
}

float* HOSegmentation::accessX(float *arr, int y, int z, int k) 
{
  return (arr + k + (y+z*volSize.height)*volSize.width);
}

float* HOSegmentation::accessY(float *arr, int z, int x, int k)
{
  return (arr + x + (k+z*volSize.height)*volSize.width);
}

float* HOSegmentation::accessZ(float *arr, int x, int y, int k)
{
  return (arr + x + (y+k*volSize.height)*volSize.width);
}

void HOSegmentation::cubicCoeff1D(int n, int x1, int x2, float z1_2n, float* c_plus, float K, float z1, float *readBuffer, float *writeBuffer, float* (HOSegmentation::*accessDim)(float*, int, int, int))
{
  float z1_k, ck_minus;
  c_plus[0] = 0.0;
  if(n<K)
    for(int _k=1; _k<n; _k++)
    {
      z1_k = powf(z1, (float)_k);
      c_plus[0] += *((this->*accessDim)(readBuffer, x1, x2, _k - 1))*(z1_k - z1_2n/z1_k);
    }
    else //Use Horner's scheme.
      for(int _k=1; _k<n; _k++)
        c_plus[0] = ( *((this->*accessDim)(readBuffer, x1, x2, n - 1 - _k)) + c_plus[0])*z1;	  	
    c_plus[0] = -c_plus[0]/(1.0 - z1_2n);      
    
  for(int _k =1; _k<n; _k++)
    c_plus[_k] = *((this->*accessDim)(readBuffer, x1, x2, _k - 1)) + z1*c_plus[_k-1];
	
  ck_minus = 0.0; // We don't need to store these values. This would be c_minus[n] = 0, at start
  for(int _k=(n-1); _k>=1; _k--)
  {
    ck_minus = z1*(ck_minus - c_plus[_k]);  
    *((this->*accessDim)(writeBuffer, x1, x2, _k - 1)) = 6.0*ck_minus;		
  }    
}

bool HOSegmentation::computeCubicCoefficients()
{
  ///Compute cubc coefficients for the signed distance field
  debugPrint("Computing cubic coefficients for phi... ");
  float z1 = sqrtf(3.0) - 2.0;  
  float K = logf((float)EPSILON)/logf(fabs(z1));
  
#if COMPUTE_COEFFICIENTS_ON_CPU
  clock_t t0, t1;
  t0 = clock();
  float z1_2n;//, z1_k, ck_minus;
  float *c_plus = NULL;

  //---------------------------------------------------------------  
  debugPrint("Z ");
  z1_2n = powf(z1, 2.0*((float)volSize.depth + 1.0));
  c_plus = (float*)malloc((volSize.depth+1)*sizeof(float));  
   
  for(int _y = 0; _y< volSize.height; _y++)
    for(int _x = 0; _x< volSize.width; _x++)
      cubicCoeff1D(volSize.depth+1, _x, _y, z1_2n, c_plus, K, z1, h_phi, h_coeff, &HOSegmentation::accessZ);
  free(c_plus); c_plus = NULL;
  //---------------------------------------------------------------        
  debugPrint("Y ");
  z1_2n = powf(z1, 2.0*((float)volSize.height + 1.0));
  c_plus = (float*)malloc((volSize.height + 1)*sizeof(float));  
   
  for(int _z = 0; _z< volSize.depth; _z++)
    for(int _x = 0; _x< volSize.width; _x++)
      cubicCoeff1D(volSize.height+1, _z, _x, z1_2n, c_plus, K, z1, h_coeff, h_coeff, &HOSegmentation::accessY); 
  free(c_plus); c_plus = NULL;    
  //---------------------------------------------------------------  
  debugPrint("X ");
  z1_2n = powf(z1, 2.0*((float)volSize.width + 1.0));
  c_plus = (float*)malloc((volSize.width + 1)*sizeof(float));  
   
  for(int _z = 0; _z< volSize.depth; _z++)
    for(int _y = 0; _y< volSize.height; _y++)
       cubicCoeff1D(volSize.width+1, _y, _z, z1_2n, c_plus, K, z1, h_coeff, h_coeff, &HOSegmentation::accessX); 
  free(c_plus); c_plus = NULL; 
  
  t1 = clock();
  milliseconds = (float)(t1-t0)*1000.0/CLOCKS_PER_SEC;  
#else
  CUT_SAFE_CALL(cutStartTimer(timer));
  cudaPitchedPtr d_coeffPPtr;// subvolume in device memory - used for cubic coefficient computations
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
  CUDA_SAFE_CALL(cudaFree(d_spare1PPtr.ptr));//Free some GPU mem.
  CUDA_SAFE_CALL(cudaFree(d_spare2PPtr.ptr));
  CUDA_SAFE_CALL(cudaFree(d_spare3PPtr.ptr));
  char _msg[256];
  const dim3 CCblockSize(CCblockDim, 1, 1); //1D block
  dim3 CCgridSize(1, 1, 1);// 1D grid

  ///===================================================================================
  ///===================================================================================
  //1. Perform coeff. computation in Z-direction
  debugPrint("\tComputing along Z:\n");
  //Get the optimal volume slice size to process
  cudaExtent coeffSubvolSize = make_cudaExtent(volSize.width, 0, volSize.depth);//Let's slice along Y
  if(!getAvailbleGPUSliceSize(coeffSubvolSize))
  {
    debugPrint("Failed to compute cubic coefficients.\n");
    return false;
  }
  debugPrint("\tSlice dim: %zd x %zd x %zd\n", coeffSubvolSize.width, coeffSubvolSize.height, coeffSubvolSize.depth);
  
  //Allocate memory and bind coeffArray to 3D texture
  allocateSliceOnGPU(d_coeffArray, d_coeffPPtr, coeffSubvolSize);
  CUDA_SAFE_CALL(cudaBindTextureToArray(texCoeff, d_coeffArray, channelDesc)); //bind texture
  CHK_CUDA_ERR("\n");  
  int nslices = iDivUp(volSize.height, coeffSubvolSize.height);
  cudaExtent copyvol;
  
  //Compute coefficients in 1 direction  
  for(int _yslice = 0; _yslice < nslices; _yslice++)
  {
    //Copy subvol to texture array. This copying is safe in the sense that bondary is taken care of.
    copyvol = coeffSubvolSize;    
    zeroOutArray = false;
    if(_yslice == (nslices -1))
    {
      copyvol.height = (volSize.height)%(coeffSubvolSize.height);
      if(copyvol.height == 0)
        copyvol.height = coeffSubvolSize.height;
      else
        zeroOutArray = true;     
    }
    if(zeroOutArray)
    {
      cudaExtent pitchedVolSize = make_cudaExtent(coeffSubvolSize.width*sizeof(float), coeffSubvolSize.height, coeffSubvolSize.depth);     
      CUDA_SAFE_CALL(cudaMemset3D(d_coeffPPtr, 0, pitchedVolSize));// zero out the device mem
      copy3DMemToArray(d_coeffPPtr, d_coeffArray, coeffSubvolSize);
    }
    copy3DHostToArray(h_phi, d_coeffArray, copyvol , volSize , make_cudaPos(0, _yslice*coeffSubvolSize.height, 0));//Values taken from h_phi only the first time.
    CHK_CUDA_ERR("->");
    //The launch configuration is such that: threads/block = 512, blocks/grid = ceil(W*H/512)
    //CCgridSize.x = iDivUp(coeffSubvolSize.width*coeffSubvolSize.height, CCblockDim);
    CCgridSize.x = iDivUp(copyvol.width*copyvol.height, CCblockDim); 
    if(volSize.depth < K)
      d_cubic_coefficients_1DZ<<<CCgridSize, CCblockSize>>>(d_coeffPPtr, copyvol);    
    else
      d_cubic_coefficients_1DZ_Fast<<<CCgridSize, CCblockSize>>>(d_coeffPPtr, copyvol);    
    CUT_CHECK_ERROR("Kernel failed");
    sprintf(_msg,  "---%d---", _yslice);
    CHK_CUDA_ERR(_msg);     	
    cudaThreadSynchronize();
    //Copy results to h_coeff host volume
    copy3DMemToHost(d_coeffPPtr, h_coeff, copyvol, coeffSubvolSize, make_cudaPos(0, 0, 0), make_cudaPos(0, _yslice*coeffSubvolSize.height, 0));      
    CHK_CUDA_ERR("<-");
  }
  
  CUDA_SAFE_CALL(cudaUnbindTexture(texCoeff)); //unbind texture  
  CHK_CUDA_ERR("\n");    
  CUDA_SAFE_CALL(cudaFreeArray(d_coeffArray));//free already allocated mem
  CHK_CUDA_ERR("\n");    
  CUDA_SAFE_CALL(cudaFree(d_coeffPPtr.ptr));
  CHK_CUDA_ERR("\n");  
  
  ///===================================================================================
  ///===================================================================================
  //2. Perform coeff. computation in Y-direction
   debugPrint("\tComputing along Y:\n");
  //Get the optimal volume slice size to process
  coeffSubvolSize = make_cudaExtent(volSize.width, volSize.height, 0);//Let's slice along Z
  if(!getAvailbleGPUSliceSize(coeffSubvolSize))
  {
    debugPrint("Failed to compute cubic coefficients.\n");
    return false;
  }
  debugPrint("\tSlice dim: %zd x %zd x %zd\n", coeffSubvolSize.width, coeffSubvolSize.height, coeffSubvolSize.depth);
  
  //Allocate memory and bind coeffArray to 3D texture
  allocateSliceOnGPU(d_coeffArray, d_coeffPPtr, coeffSubvolSize);
  CUDA_SAFE_CALL(cudaBindTextureToArray(texCoeff, d_coeffArray, channelDesc)); //bind texture
  CHK_CUDA_ERR("\n");  
  nslices = iDivUp(volSize.depth, coeffSubvolSize.depth);
  
  //Compute coefficients in 1 direction  
  for(int _zslice = 0; _zslice < nslices; _zslice++)
  {
    //Copy subvol to texture array. This copying is safe in the sense that bondary is taken care of.
    copyvol = coeffSubvolSize;    
    zeroOutArray = false;
    if(_zslice == (nslices -1))
    {
      copyvol.depth = (volSize.depth)%(coeffSubvolSize.depth);
      if(copyvol.depth == 0)
        copyvol.depth = coeffSubvolSize.depth;
      else
        zeroOutArray = true;     
    }
    if(zeroOutArray)
    {
      cudaExtent pitchedVolSize = make_cudaExtent(coeffSubvolSize.width*sizeof(float), coeffSubvolSize.height, coeffSubvolSize.depth);         
      CUDA_SAFE_CALL(cudaMemset3D(d_coeffPPtr, 0, pitchedVolSize));// zero out the device mem
      copy3DMemToArray(d_coeffPPtr, d_coeffArray, coeffSubvolSize);
    }
    copy3DHostToArray(h_coeff, d_coeffArray, copyvol , volSize , make_cudaPos(0, 0, _zslice*coeffSubvolSize.depth));
    CHK_CUDA_ERR("->");
    
    //The launch configuration is such that: threads/block = 512, blocks/grid = ceil(W*D/512)
    //CCgridSize.x = iDivUp(coeffSubvolSize.width*coeffSubvolSize.depth, CCblockDim);
    CCgridSize.x = iDivUp(copyvol.width*copyvol.depth, CCblockDim);     
    if(volSize.height < K)       
      d_cubic_coefficients_1DY<<<CCgridSize, CCblockSize>>>(d_coeffPPtr, coeffSubvolSize);    
    else
      d_cubic_coefficients_1DY_Fast<<<CCgridSize, CCblockSize>>>(d_coeffPPtr, coeffSubvolSize);          
    CUT_CHECK_ERROR("Kernel failed");
    sprintf(_msg,  "---%d---", _zslice);
    CHK_CUDA_ERR(_msg);     	
    cudaThreadSynchronize();
    //Copy results to h_coeff host volume
    copy3DMemToHost(d_coeffPPtr, h_coeff, copyvol, volSize, make_cudaPos(0, 0, 0), make_cudaPos(0, 0, _zslice*coeffSubvolSize.depth));
    CHK_CUDA_ERR("<-");
  }
  ///===================================================================================
  ///===================================================================================
  //3. Perform coeff. computation in X-direction  
  debugPrint("\tComputing along X:\n");
  //Let's slice along Z again - this means that the slice size doesn't change.
  debugPrint("\tSlice dim: %zd x %zd x %zd\n", coeffSubvolSize.width, coeffSubvolSize.height, coeffSubvolSize.depth);
  
  //Compute coefficients in 1 direction  
  for(int _zslice = 0; _zslice < nslices; _zslice++)
  {
    //Copy subvol to texture array. This copying is safe in the sense that bondary is taken care of.
    copyvol = coeffSubvolSize;    
    zeroOutArray = false;
    if(_zslice == (nslices -1))
    {
      copyvol.depth = (volSize.depth)%(coeffSubvolSize.depth);
      if(copyvol.depth == 0)
        copyvol.depth = coeffSubvolSize.depth;
      else
        zeroOutArray = true;     
    }
    if(zeroOutArray)
    {
      cudaExtent pitchedVolSize = make_cudaExtent(coeffSubvolSize.width*sizeof(float), coeffSubvolSize.height, coeffSubvolSize.depth);         
      CUDA_SAFE_CALL(cudaMemset3D(d_coeffPPtr, 0, pitchedVolSize));// zero out the device mem
      copy3DMemToArray(d_coeffPPtr, d_coeffArray, coeffSubvolSize);
    }
    copy3DHostToArray(h_coeff, d_coeffArray, copyvol , volSize , make_cudaPos(0, 0, _zslice*coeffSubvolSize.depth));
    CHK_CUDA_ERR("->");
    
    //The launch configuration is such that: threads/block = 512, blocks/grid = ceil(H*D/512)
    //CCgridSize.x = iDivUp(coeffSubvolSize.height*coeffSubvolSize.depth, CCblockDim);
    CCgridSize.x = iDivUp(copyvol.height*copyvol.depth, CCblockDim);
    if(volSize.width < K)                
      d_cubic_coefficients_1DX<<<CCgridSize, CCblockSize>>>(d_coeffPPtr, coeffSubvolSize);    
    else
      d_cubic_coefficients_1DX_Fast<<<CCgridSize, CCblockSize>>>(d_coeffPPtr, coeffSubvolSize);          
    CUT_CHECK_ERROR("Kernel failed");
    sprintf(_msg,  "---%d---", _zslice);
    CHK_CUDA_ERR(_msg);     	
    cudaThreadSynchronize();
    //Copy results to h_coeff host volume
    copy3DMemToHost(d_coeffPPtr, h_coeff, copyvol, volSize, make_cudaPos(0, 0, 0), make_cudaPos(0, 0, _zslice*coeffSubvolSize.depth));
    CHK_CUDA_ERR("<-");
  }
  
  CUDA_SAFE_CALL(cudaUnbindTexture(texCoeff)); //unbind texture  
  CHK_CUDA_ERR("\n");    
  CUDA_SAFE_CALL(cudaFreeArray(d_coeffArray));//free already allocated mem
  CHK_CUDA_ERR("\n");    
  CUDA_SAFE_CALL(cudaFree(d_coeffPPtr.ptr));
  CHK_CUDA_ERR("\n");  
  ///===================================================================================
  ///===================================================================================  

  //Allocate freed spare mem.
  cudaExtent pitchedVolSize = make_cudaExtent(subvolDim*sizeof(float), subvolDim, subvolDim); 
  CUDA_SAFE_CALL(cudaMalloc3D(&d_spare1PPtr, pitchedVolSize));
  CHK_CUDA_ERR("\n");
  CUDA_SAFE_CALL(cudaMalloc3D(&d_spare2PPtr, pitchedVolSize));
  CHK_CUDA_ERR("\n");
  CUDA_SAFE_CALL(cudaMalloc3D(&d_spare3PPtr, pitchedVolSize));
  CHK_CUDA_ERR("\n");
  
  CUT_SAFE_CALL(cutStopTimer(timer));
  milliseconds = cutGetAverageTimerValue(timer);
  CUT_SAFE_CALL(cutResetTimer(timer));
#endif  
  debugPrint("%f ms\n", milliseconds);
  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@     
  /*FILE *fperr; // Redirect console to file.
  bool err_redirected = false;
  if ((fperr = freopen("stderr.txt","w",stderr)) != NULL)
   err_redirected = true;
   
  debugPrint("Testing computed coefficients...\n");
  float result; int counter = 0;
  for(int _x = 0; _x< volSize.width; _x++)  
    for(int _y = 0; _y< volSize.height; _y++){
      for(int _z = 0; _z< volSize.depth; _z++){
        result = evaluateCubicSplineAtGridPoint(_x, _y, _z);
        debugPrint("%8.5f ", result);	
      }
      debugPrint("\n");      
      for(int _z = 0; _z< volSize.depth; _z++)
        debugPrint("%8.5f ", VOXEL_AT(h_phi, _x, _y, _z));
      debugPrint("\n");	
      debugPrint("<%d %d> ", _x, _y);	
      for(int _z = 0; _z< volSize.depth; _z++)
        debugPrint("%8.5f ", VOXEL_AT(h_coeff, _x, _y, _z));
      debugPrint("\n");	
    }
  debugPrint("Testing complete.");
  if (err_redirected) freopen("CON","w",stderr);*/
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  return true;
}

float HOSegmentation::divergenceOfNormalizedGradient_FiniteDifference( int x,  int y,  int z)
{
  float grad_left[3], grad_right[3], div;
  float l_phi;
  float grad_norm;
  if((x>1 && x< (volSize.width-2)) && (y>1 && y<(volSize.height-2)) && (z>1 && z<(volSize.depth-2)))
  {
	l_phi = VOXEL_AT(h_phi, x, y, z);
	grad_left[0] = (l_phi - VOXEL_AT(h_phi, x - 2, y, z))/2.0;
	grad_left[1] = (l_phi - VOXEL_AT(h_phi, x, y - 2, z))/2.0;
	grad_left[2] = (l_phi - VOXEL_AT(h_phi, x, y, z - 2))/2.0;
	grad_right[0] = (VOXEL_AT(h_phi, x + 2, y, z) - l_phi)/2.0;
	grad_right[1] = (VOXEL_AT(h_phi, x, y + 2, z) - l_phi)/2.0;
	grad_right[2] = (VOXEL_AT(h_phi, x, y, z + 2) - l_phi)/2.0;	
  }
  else
  {
	l_phi = SAFE_VOXEL_AT(h_phi, x, y, z);
	grad_left[0] = (l_phi - SAFE_VOXEL_AT(h_phi, x - 2, y, z))/2.0;
	grad_left[1] = (l_phi - SAFE_VOXEL_AT(h_phi, x, y - 2, z))/2.0;
	grad_left[2] = (l_phi - SAFE_VOXEL_AT(h_phi, x, y, z - 2))/2.0;
	grad_right[0] = (SAFE_VOXEL_AT(h_phi, x + 2, y, z) - l_phi)/2.0;
	grad_right[1] = (SAFE_VOXEL_AT(h_phi, x, y + 2, z) - l_phi)/2.0;
	grad_right[2] = (SAFE_VOXEL_AT(h_phi, x, y, z + 2) - l_phi)/2.0;	
  }
	grad_norm = sqrt(grad_left[0]*grad_left[0] + grad_left[1]*grad_left[1] + grad_left[2]*grad_left[2] + TINY);
	grad_left[0] /= grad_norm; grad_left[1] /= grad_norm; grad_left[2] /= grad_norm;
	grad_norm = sqrt(grad_right[0]*grad_right[0] + grad_right[1]*grad_right[1] + grad_right[2]*grad_right[2] + TINY);
	grad_right[0] /= grad_norm; grad_right[1] /= grad_norm; grad_right[2] /= grad_norm;

	div = (grad_right[0] - grad_left[0] + grad_right[1] - grad_left[1] + grad_right[2] - grad_left[2])/2.0;
	return div;
}

float HOSegmentation::divergenceOfNormalizedGradient_CubicSpline( int x,  int y,  int z)
{
  float phix=0.0, phiy=0.0, phiz=0.0, phixx=0.0, phiyy=0.0, phizz=0.0, phixy=0.0, phiyz=0.0, phizx=0.0;  
  float coeff = 0.0;
  float *ptr;
  if((x>0 && x< (volSize.width-1)) && (y>0 && y<(volSize.height-1)) && (z>0 && z<(volSize.depth-1)))
    for(int k=0; k<3; k++)
      for(int j=0; j<3; j++)
        for(int i=0; i<3; i++)
        {
          ptr = h_cubicDerivatives + (i+3*j+9*k)*9;	
          coeff = VOXEL_AT(h_coeff,x+i-1, y+j-1, z+k-1);
          phix += coeff*ptr[0];
          phiy += coeff*ptr[1];
          phiz += coeff*ptr[2];
          phixx += coeff*ptr[3];
          phiyy += coeff*ptr[4];
          phizz += coeff*ptr[5];
          phixy += coeff*ptr[6];
          phiyz += coeff*ptr[7];
          phizx += coeff*ptr[8];
        }    
  else
    for(int k=0; k<3; k++)
      for(int j=0; j<3; j++)
        for(int i=0; i<3; i++)
        {
          ptr = h_cubicDerivatives + (i+3*j+9*k)*9;	
          coeff = SAFE_VOXEL_AT(h_coeff, (int)x+i-1, (int)y+j-1, (int)z+k-1);
          phix += coeff*ptr[0];
          phiy += coeff*ptr[1];
          phiz += coeff*ptr[2];
          phixx += coeff*ptr[3];
          phiyy += coeff*ptr[4];
          phizz += coeff*ptr[5];
          phixy += coeff*ptr[6];
          phiyz += coeff*ptr[7];
          phizx += coeff*ptr[8];
        }
  
  //Reuse variable coeff to store ||grad(phi)||^2
  coeff = phix*phix + phiy*phiy + phiz*phiz + TINY;
  coeff *=sqrtf(coeff);
  //real divNormGrad = phix*phix*phixx + phiy*phiy*phiyy + phiz*phiz*phizz;
  //divNormGrad += 2.0*(phix*phiy*phixy + phiy*phiz*phiyz + phiz*phix*phizx);
  //divNormGrad /= coeff;
  //divNormGrad = (phixx + phiyy + phizz) - divNormGrad;
  //divNormGrad /=sqrtf(coeff);
  float divNormGrad = phix*phix*(phiyy + phizz) + phiy*phiy*(phixx + phizz) + phiz*phiz*(phixx + phiyy) - 2.0*(phix*phiy*phixy + phiy*phiz*phiyz + phiz*phix*phizx);
  divNormGrad /=coeff;
  
//  float CS[] = {1.0/6.0, 4.0/6.0, 1.0/6.0};
//  float CSx[]= {-1.0/2.0, 0.0, 1.0/2.0};
//  float CSxx[] = {1.0, -2.0, 1.0};
//  phix = evaluateCubicSplineGenericDerivativeAtGridPoint(x, y, z, CSx, CS, CS);
//  phiy = evaluateCubicSplineGenericDerivativeAtGridPoint(x, y, z, CS, CSx, CS);	
//  phiz = evaluateCubicSplineGenericDerivativeAtGridPoint(x, y, z, CS, CS, CSx);
//  phixx = evaluateCubicSplineGenericDerivativeAtGridPoint(x, y, z, CSxx, CS, CS);
//  phiyy = evaluateCubicSplineGenericDerivativeAtGridPoint(x, y, z, CS, CSxx, CS);
//  phizz = evaluateCubicSplineGenericDerivativeAtGridPoint(x, y, z, CS, CS, CSxx);
//  phixy = evaluateCubicSplineGenericDerivativeAtGridPoint(x, y, z, CSx, CSx, CS);
//  phiyz = evaluateCubicSplineGenericDerivativeAtGridPoint(x, y, z, CS, CSx, CSx);
//  phizx = evaluateCubicSplineGenericDerivativeAtGridPoint(x, y, z, CSx, CS, CSx);
	      
  return divNormGrad;  	      
}

// The bordering scheme changes here since we need 1 voxel border around subvol such that part of the actual
// volume boundary has a border with zero values. Thus an effectvive border of zeros around the full volume.
// This is needed only for the coeff texture. For the other two textures we don't care since computation is dependent 
// only on the current (x, y, z). Effective computational volume: (subvolDim-2)^3
bool HOSegmentation::PDEUpdate()
{
  /// Update phi by single time step.
  ///Compute the higher and lower order terms.
  //Warning: h_subvolSpare, d_volPPtr have been used in initializing textures in loops. Do not use later them later in the loop or modifiy accordingly.
  debugPrint("Updating the PDE (level set and input volume)... ");
  float mu_factor = mu/(h*h);
#if UPDATE_HIGHER_ORDER_PDE_ON_CPU
  clock_t t0, t1;
  t0 = clock();
  //float divFD, divCS;
  float internal_energy, external_energy, phival, volval, delta_dirac;
  for( int _z = 0; _z< volSize.depth; _z++)
    for( int _y = 0; _y< volSize.height; _y++)
      for( int _x = 0; _x< volSize.width; _x++)
      {
        ///Compute higher order term (internal energy term) using Cubic splines
        //internal_energy = mu_factor*divergenceOfNormalizedGradient_FiniteDifference(_x, _y, _z);		
        internal_energy = mu_factor*divergenceOfNormalizedGradient_CubicSpline(_x, _y, _z);
/*	
         //-----------------------This code is to compare divergence computation here with the one from Zhang (see ../HLevelSet/KLevelSet_Recon2.cpp)--------------
         float partials[10];
         float phix=0.0, phiy=0.0, phiz=0.0, phixx=0.0, phiyy=0.0, phizz=0.0, phixy=0.0, phiyz=0.0, phizx=0.0;  
         EvaluateCubicSplineOrder2PartialsAtGridPoint(h_coeff, 1, 1, 1, volSize.width, volSize.height, volSize.depth, _x, _y, _z, partials);
         phix = partials[1];
         phiy = partials[2];
         phiz = partials[3];
         phixx = partials[4];
         phiyy = partials[7];
         phizz = partials[9];
         phixy = partials[5];
         phiyz = partials[8];
         phizx = partials[6];
         float coeff = phix*phix + phiy*phiy + phiz*phiz + TINY;
         coeff *=sqrtf(coeff);
         float divNormGrad = phix*phix*(phiyy + phizz) + phiy*phiy*(phixx + phizz) + phiz*phiz*(phixx + phiyy) - 2.0*(phix*phiy*phixy + phiy*phiz*phiyz + phiz*phix*phizx);
         divNormGrad /=coeff;	
         internal_energy = mu_factor*divNormGrad;  
         //------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	
*/	 
        volval = VOXEL_AT(h_vol, _x, _y, _z);
        external_energy = -nu - lambda1*(volval - c1)*(volval - c1) + lambda2*(volval - c2)*(volval - c2);	
	
        ///Compute external energy
        phival = VOXEL_AT(h_buffer[currentPhi], _x, _y, _z);
        delta_dirac = (epsilon/M_PI)/(epsilon*epsilon + phival*phival);
        VOXEL_AT(h_buffer[currentPhi], _x, _y, _z) = phival + delta_t*delta_dirac*(internal_energy + external_energy);
      }
  t1 = clock();
  milliseconds = (float)(t1-t0)*1000.0/CLOCKS_PER_SEC;    
#else
  const dim3 PDEblockSize(PDEBlockDim, PDEBlockDim, PDEBlockDim);
  const dim3 PDEgridSize(iDivUp(subvolDim, PDEBlockDim)*iDivUp(subvolDim, PDEBlockDim)*iDivUp(subvolDim, PDEBlockDim), 1);//linearized grid, since, we cannot have a 3D grid
  cudaExtent logicalGridSize = make_cudaExtent(iDivUp(subvolDim, PDEBlockDim), iDivUp(subvolDim, PDEBlockDim), iDivUp(subvolDim, PDEBlockDim));
  subvolIndicesExtents = make_cudaExtent(iDivUp(volSize.width-2, subvolDim-2), iDivUp(volSize.height-2, subvolDim-2), iDivUp(volSize.depth-2, subvolDim-2));
  cudaExtent copyvol = make_cudaExtent(subvolDim-2, subvolDim-2, subvolDim-2);
  cudaPos src_off = make_cudaPos(1, 1, 1);
  cudaPos dst_off = make_cudaPos(0, 0, 0);
  char _msg[256];
  float maxDist = sqrtf(float(volSize.width*volSize.width + volSize.height*volSize.height + volSize.depth*volSize.depth));

  //Since d_coeffArray has also been used while computing coefficients, we neet to reallocate memory and bind again to texCoeff.
  
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
  CUDA_SAFE_CALL(cudaMalloc3DArray(&d_coeffArray, &channelDesc, subvolSize));
  CHK_CUDA_ERR("\n");
  CUDA_SAFE_CALL(cudaBindTextureToArray(texCoeff, d_coeffArray, channelDesc)); //bind texture
  CUDA_SAFE_CALL(cudaBindTextureToArray(texPhi, d_phiArray, channelDesc)); //bind texture
  CUDA_SAFE_CALL(cudaBindTextureToArray(texVol, d_volArray, channelDesc)); //bind texture
  CHK_CUDA_ERR("\n");
  
  //We copy the boundary to the texture, so the extra volume uploaded must have the same value as the outside of the interface initialization.
  for(int i=0; i< subvolDim*subvolDim*subvolDim; i++)      
    *(h_subvolSpare + i) = -maxDist;	
  
  cudaExtent copyvol_upload = make_cudaExtent(0,0,0);
  cudaPos offset_upload = make_cudaPos(0, 0, 0);
  CUT_SAFE_CALL(cutStartTimer(timer));
  for(int _z = 0; _z< subvolIndicesExtents.depth; _z++)
    for(int _y = 0; _y< subvolIndicesExtents.height; _y++)
      for(int _x = 0; _x< subvolIndicesExtents.width; _x++)
      {
        adjustUploadPDESubvolSizeAndOffset(_x, _y, _z, copyvol_upload, offset_upload);
        //adjustUploadSubvolSize(_x, _y, _z, 2, copyvol_upload);
        
        //load texture values for coeff, phi and intensity values
        // Watchout for boundaries at the end. We don't want to read garbage values outside the allocated host volume!!
        if(zeroOutArray)
        {
          //zero out Arrays first so they do not have undesired values stored in them from previous load
          //N.B.: there seems to be no API to memset a cudaArray, so we do it indirectly
          cudaMemset3D(d_volPPtr, 0, make_cudaExtent(subvolDim*sizeof(float), subvolDim, subvolDim));// zero out the device mem
          copy3DMemToArray(d_volPPtr, d_volArray); 
          copy3DMemToArray(d_volPPtr, d_coeffArray); 
          //Fill the phi texture with -maxDist so that it doesn't have undesirable values when bundary volume parts are copied.
          copy3DHostToArray(h_subvolSpare, d_phiArray, subvolSize, subvolSize, make_cudaPos(0, 0, 0)); 
        }

        copy3DHostToArray(h_coeff, d_coeffArray, copyvol_upload, volSize, offset_upload);
        //copy3DHostToArray(h_coeff, d_coeffArray, copyvol_upload, volSize, make_cudaPos(_x*(subvolDim-4), _y*(subvolDim-4), _z*(subvolDim-4)));
        CHK_CUDA_ERR("->");
        copy3DHostToArray(h_buffer[currentPhi], d_phiArray, copyvol_upload, volSize, offset_upload);
        //copy3DHostToArray(h_buffer[currentPhi], d_phiArray, copyvol_upload, volSize, make_cudaPos(_x*(subvolDim-4), _y*(subvolDim-4), _z*(subvolDim-4)));
        CHK_CUDA_ERR("->");
        copy3DHostToArray(h_vol, d_volArray, copyvol_upload, volSize, offset_upload);
        //copy3DHostToArray(h_vol, d_volArray, copyvol_upload, volSize, make_cudaPos(_x*(subvolDim-4), _y*(subvolDim-4), _z*(subvolDim-4)));
        CHK_CUDA_ERR("->");

        //Update the PDE
        HigherOrder_PDE_update<<<PDEgridSize, PDEblockSize>>>(d_volPPtr, logicalGridSize, mu_factor, nu, make_float2(lambda1, lambda2), make_float2(c1, c2), epsilon, delta_t, subvolDim, PDEBlockDim);

        CUT_CHECK_ERROR("Kernel failed");
        sprintf(_msg,  "---%d-%d-%d", _x, _y, _z);
        CHK_CUDA_ERR(_msg);     	
        cudaThreadSynchronize();

        // Copy results back to h_phi and h_vol
        // Watchout for boundaries at the end. We don't want to write to volume outside allocated array!!
        adjustDnloadSubvolSize(_x, _y, _z, 2, copyvol);
        dst_off.x = 1+_x*(subvolDim-2); dst_off.y = 1+_y*(subvolDim-2); dst_off.z = 1+_z*(subvolDim-2);
        copy3DMemToHost(d_volPPtr, h_buffer[(currentPhi+1)%2], copyvol, volSize, src_off, dst_off); 
        CHK_CUDA_ERR("<-%s");
      }

  CUDA_SAFE_CALL(cudaUnbindTexture(texCoeff)); //unbind texture  
  CHK_CUDA_ERR("\n");    
  CUDA_SAFE_CALL(cudaFreeArray(d_coeffArray));//free already allocated mem
  CHK_CUDA_ERR("\n");    

  currentPhi = (currentPhi+1)%2;
  h_phi = h_buffer[currentPhi];

  CUT_SAFE_CALL(cutStopTimer(timer));
  milliseconds = cutGetAverageTimerValue(timer);
  CUT_SAFE_CALL(cutResetTimer(timer));

#endif
  debugPrint(" %f ms\n", milliseconds);  

  return true;
}

//Debug this routine for possible bugs.
void HOSegmentation::adjustUploadPDESubvolSizeAndOffset(int _x, int _y, int _z, cudaExtent &copyvol_upload, cudaPos &offset_upload) 
{
  cudaExtent _subvolDim = make_cudaExtent(subvolDim, subvolDim, subvolDim);
  //Handle volume size smaller than subvolume size!
  if(volSize.width < subvolDim) _subvolDim.width = volSize.width;
  if(volSize.height < subvolDim) _subvolDim.height = volSize.height;
  if(volSize.depth < subvolDim) _subvolDim.depth = volSize.depth;

  zeroOutArray = false;
  copyvol_upload.width = _subvolDim.width;
  copyvol_upload.height = _subvolDim.height;
  copyvol_upload.depth = _subvolDim.depth;

  offset_upload.x = _x*(_subvolDim.width - 2) - 1;
  offset_upload.y = _y*(_subvolDim.height - 2) - 1;
  offset_upload.z = _z*(_subvolDim.depth - 2) - 1;
  
  if(_x == 0)
  {
    copyvol_upload.width--; //1 zero value voxel to the left should be excluded
    offset_upload.x = 0;
    zeroOutArray = true;
  }
  else if(_x == (subvolIndicesExtents.width-1))
  {
    copyvol_upload.width = ((volSize.width) % (_subvolDim.width-2)) + 1; //1 voxel border to the left should be included
    //if(_x == 0 ) copyvol_upload.width--; //In case subvolIndicesExtents.width = 1
    zeroOutArray = true;
  }
  if(_y == 0)
  {
    copyvol_upload.height--; //1 zero value voxel to the left should be excluded
    offset_upload.y = 0;
    zeroOutArray = true;
  }  
  else if(_y == (subvolIndicesExtents.height-1))
  {
    copyvol_upload.height = ((volSize.height) % (_subvolDim.height-2)) + 1;
    //if(_y == 0 ) copyvol_upload.height--;     
    zeroOutArray = true;
  }
  if(_z == 0)
  {
    copyvol_upload.depth--; //1 zero value voxel to the left should be excluded
    offset_upload.z = 0;
    zeroOutArray = true;
  }
  else if(_z == (subvolIndicesExtents.depth-1))
  {
    copyvol_upload.depth =( (volSize.depth) % (_subvolDim.depth-2)) + 1; 
    //if(_z == 0 ) copyvol_upload.depth--;    
    zeroOutArray = true;
  }
}

void HOSegmentation::copy3DMemToArray(cudaPitchedPtr _src, cudaArray *_dst, cudaExtent copy_extent)
{
  cudaMemcpy3DParms copyParams = {0};
  copyParams.srcPtr = 	_src;
  copyParams.dstArray = _dst;
  copyParams.kind = cudaMemcpyDeviceToDevice;
  copyParams.extent = copy_extent;

  CUDA_SAFE_CALL(cudaMemcpy3D(&copyParams));
  CUT_CHECK_ERROR("Mem -> Array Memcpy failed\n");
}

void HOSegmentation::copy3DArrayToHost(cudaArray *_src, float *_dst, cudaExtent copy_extent, cudaExtent dst_extent, cudaPos dst_offset)
{
  cudaMemcpy3DParms copyParams = {0};
  float *h_target = _dst + dst_offset.x + dst_offset.y*dst_extent.width + dst_offset.z*dst_extent.width*dst_extent.height;//For some reason, using copyParams.dstPos doesn't give correct results, so we set the offset here.
  copyParams.dstPtr = make_cudaPitchedPtr((void*)h_target, dst_extent.width*sizeof(float), dst_extent.width, dst_extent.height);
  copyParams.srcArray = _src;
  copyParams.kind = cudaMemcpyDeviceToHost;
  copyParams.extent = copy_extent;

  CUDA_SAFE_CALL(cudaMemcpy3D(&copyParams));
  CUT_CHECK_ERROR("Host -> Array Memcpy failed\n");
}

void HOSegmentation::copy3DArrayToMem(cudaArray *_src, cudaPitchedPtr _dst, cudaExtent copy_extent)
{
  cudaMemcpy3DParms copyParams = {0};
  copyParams.dstPtr = _dst;
  copyParams.srcArray = _src;
  copyParams.kind = cudaMemcpyDeviceToDevice;
  copyParams.extent = copy_extent;

  CUDA_SAFE_CALL(cudaMemcpy3D(&copyParams));
  CUT_CHECK_ERROR("Host -> Array Memcpy failed\n");
}

/*-----------------------------------------------------------------------------*/
void  HOSegmentation::EvaluateCubicSplineOrder2PartialsAtGridPoint(float *c, float dx, float dy, float dz, int nx, int ny, int nz, int u, int v, int w, float *partials)
//              float   *c,    /* the spline  coefficients                            */
//              float   dx     /* spacing in x direction                              */
//              float   dy     /* spacing in y direction                              */
//              float   dz     /* spacing in z direction                              */
//              int     nx,    /* number of samples or coefficients in x direction    */
//              int     ny,    /* number of samples or coefficients in y direction    */
//              int     nz,    /* number of samples or coefficients in z direction    */
//              float   *partial /* partial derivatives                               */
{
float c27[27], indx, indy, indz;
int  i;
// Local variables
float  TensorF[27],   TensorFx[27],  TensorFy[27],  TensorFz[27],  TensorFxx[27], TensorFxy[27], TensorFxz[27], TensorFyy[27], TensorFyz[27], TensorFzz[27];
ComputeTensorXYZ(TensorF, TensorFx, TensorFy, TensorFz, TensorFxx, TensorFxy, TensorFxz, TensorFyy, TensorFyz, TensorFzz);
	      
Take_27_Coefficients(c, nx, ny, nz, u, v, w, c27);

for (i = 0; i < 10; i++) {
   partials[i] = 0.0;
}

for (i = 0; i < 27; i++) {
   partials[0] = partials[0] + c27[i]*TensorF[i];

   partials[1] = partials[1] + c27[i]*TensorFx[i];
   partials[2] = partials[2] + c27[i]*TensorFy[i];
   partials[3] = partials[3] + c27[i]*TensorFz[i];

   partials[4] = partials[4] + c27[i]*TensorFxx[i];
   partials[5] = partials[5] + c27[i]*TensorFxy[i];
   partials[6] = partials[6] + c27[i]*TensorFxz[i];

   partials[7] = partials[7] + c27[i]*TensorFyy[i];
   partials[8] = partials[8] + c27[i]*TensorFyz[i];

   partials[9] = partials[9] + c27[i]*TensorFzz[i];

   //printf("i = %d   Coeff = %f,     Fxx = %f\n", i, c27[i], TensorFxx[i]);
}      

indx = 1.0/dx;
indy = 1.0/dy;
indz = 1.0/dz;

//printf("indx = %f\n", partials[4]);

partials[1] = partials[1] * indx;
partials[2] = partials[2] * indy;
partials[3] = partials[3] * indz;

partials[4] = partials[4] * indx*indx;
partials[5] = partials[5] * indx*indy;
partials[6] = partials[6] * indx*indz;

partials[7] = partials[7] * indy*indy;
partials[8] = partials[8] * indy*indz;

partials[9] = partials[9] * indz*indz;

//printf("indx = %f\n", partials[4]);

}

void HOSegmentation::Take_27_Coefficients(float *c, int nx, int ny, int nz, int u, int v, int w, float *c27)
//              float   *c,    /* the spline  coefficients                            */
//              int     nx,    /* number of samples or coefficients in x direction    */
//              int     ny,    /* number of samples or coefficients in y direction    */
//              int     nz,    /* number of samples or coefficients in z direction    */
//              float   *c27   /* 27 coefficients                                     */
{
if (( u  > 0  && u < nx-1) &&
    ( v  > 0  && v < ny-1) &&
    ( w  > 0  && w < nz-1) ) {

    c27[0] = TakeACoefficient_Fast(c, nx, ny, nz, u-1,v-1,w-1);
    c27[1] = TakeACoefficient_Fast(c, nx, ny, nz, u-1,v-1,w);
    c27[2] = TakeACoefficient_Fast(c, nx, ny, nz, u-1,v-1,w+1);

    c27[3] = TakeACoefficient_Fast(c, nx, ny, nz, u-1,v,w-1);
    c27[4] = TakeACoefficient_Fast(c, nx, ny, nz, u-1,v,w);
    c27[5] = TakeACoefficient_Fast(c, nx, ny, nz, u-1,v,w+1);
                      
    c27[6] = TakeACoefficient_Fast(c, nx, ny, nz, u-1,v+1,w-1);
    c27[7] = TakeACoefficient_Fast(c, nx, ny, nz, u-1,v+1,w);
    c27[8] = TakeACoefficient_Fast(c, nx, ny, nz, u-1,v+1,w+1);

    c27[9]  = TakeACoefficient_Fast(c, nx, ny, nz, u,v-1,w-1);
    c27[10] = TakeACoefficient_Fast(c, nx, ny, nz, u,v-1,w);
    c27[11] = TakeACoefficient_Fast(c, nx, ny, nz, u,v-1,w+1);

    c27[12] = TakeACoefficient_Fast(c, nx, ny, nz, u,v,w-1);
    c27[13] = TakeACoefficient_Fast(c, nx, ny, nz, u,v,w);
    c27[14] = TakeACoefficient_Fast(c, nx, ny, nz, u,v,w+1);
                      
    c27[15] = TakeACoefficient_Fast(c, nx, ny, nz, u,v+1,w-1);
    c27[16] = TakeACoefficient_Fast(c, nx, ny, nz, u,v+1,w);
    c27[17] = TakeACoefficient_Fast(c, nx, ny, nz, u,v+1,w+1);

    c27[18] = TakeACoefficient_Fast(c, nx, ny, nz, u+1,v-1,w-1);
    c27[19] = TakeACoefficient_Fast(c, nx, ny, nz, u+1,v-1,w);
    c27[20] = TakeACoefficient_Fast(c, nx, ny, nz, u+1,v-1,w+1);

    c27[21] = TakeACoefficient_Fast(c, nx, ny, nz, u+1,v,w-1);
    c27[22] = TakeACoefficient_Fast(c, nx, ny, nz, u+1,v,w);
    c27[23] = TakeACoefficient_Fast(c, nx, ny, nz, u+1,v,w+1);
                      
    c27[24] = TakeACoefficient_Fast(c, nx, ny, nz, u+1,v+1,w-1);
    c27[25] = TakeACoefficient_Fast(c, nx, ny, nz, u+1,v+1,w);
    c27[26] = TakeACoefficient_Fast(c, nx, ny, nz, u+1,v+1,w+1);

    return;
}

    c27[0] = TakeACoefficient_Slow(c, nx, ny, nz, u-1,v-1,w-1);
    c27[1] = TakeACoefficient_Slow(c, nx, ny, nz, u-1,v-1,w);
    c27[2] = TakeACoefficient_Slow(c, nx, ny, nz, u-1,v-1,w+1);

    c27[3] = TakeACoefficient_Slow(c, nx, ny, nz, u-1,v,w-1);
    c27[4] = TakeACoefficient_Slow(c, nx, ny, nz, u-1,v,w);
    c27[5] = TakeACoefficient_Slow(c, nx, ny, nz, u-1,v,w+1);
                      
    c27[6] = TakeACoefficient_Slow(c, nx, ny, nz, u-1,v+1,w-1);
    c27[7] = TakeACoefficient_Slow(c, nx, ny, nz, u-1,v+1,w);
    c27[8] = TakeACoefficient_Slow(c, nx, ny, nz, u-1,v+1,w+1);

    c27[9]  = TakeACoefficient_Slow(c, nx, ny, nz, u,v-1,w-1);
    c27[10] = TakeACoefficient_Slow(c, nx, ny, nz, u,v-1,w);
    c27[11] = TakeACoefficient_Slow(c, nx, ny, nz, u,v-1,w+1);

    c27[12] = TakeACoefficient_Slow(c, nx, ny, nz, u,v,w-1);
    c27[13] = TakeACoefficient_Slow(c, nx, ny, nz, u,v,w);
    c27[14] = TakeACoefficient_Slow(c, nx, ny, nz, u,v,w+1);
                      
    c27[15] = TakeACoefficient_Slow(c, nx, ny, nz, u,v+1,w-1);
    c27[16] = TakeACoefficient_Slow(c, nx, ny, nz, u,v+1,w);
    c27[17] = TakeACoefficient_Slow(c, nx, ny, nz, u,v+1,w+1);

    c27[18] = TakeACoefficient_Slow(c, nx, ny, nz, u+1,v-1,w-1);
    c27[19] = TakeACoefficient_Slow(c, nx, ny, nz, u+1,v-1,w);
    c27[20] = TakeACoefficient_Slow(c, nx, ny, nz, u+1,v-1,w+1);

    c27[21] = TakeACoefficient_Slow(c, nx, ny, nz, u+1,v,w-1);
    c27[22] = TakeACoefficient_Slow(c, nx, ny, nz, u+1,v,w);
    c27[23] = TakeACoefficient_Slow(c, nx, ny, nz, u+1,v,w+1);
                      
    c27[24] = TakeACoefficient_Slow(c, nx, ny, nz, u+1,v+1,w-1);
    c27[25] = TakeACoefficient_Slow(c, nx, ny, nz, u+1,v+1,w);
    c27[26] = TakeACoefficient_Slow(c, nx, ny, nz, u+1,v+1,w+1);
}
float HOSegmentation::TakeACoefficient_Fast(float *c, int nx, int ny, int nz, int u, int v, int w)
{
    //return(c[(u*ny + v)*nz + w]);
    return(c[(w*ny + v)*nx + u]);
    
}


/*-----------------------------------------------------------------------------*/
float HOSegmentation::TakeACoefficient_Slow(float *c, int nx, int ny, int nz, int u, int v, int w)
{
float result;

result = 0.0;
if (( u  >= 0  && u < nx) &&
    ( v  >= 0  && v < ny) &&
    ( w  >= 0  && w < nz) ) {

    //result = c[(u*ny + v)*nz + w];
    result = c[(w*ny + v)*nx + u];
    
}

return(result);
}

void HOSegmentation::ComputeTensorXYZ(float *TensorF, float *TensorFx, float *TensorFy, float *TensorFz, float *TensorFxx, float *TensorFxy, float *TensorFxz, float *TensorFyy, float *TensorFyz, float *TensorFzz)
{
float f[3], fx[3], fxx[3];

f[0] = 1.0/6.0;
f[1] = 2.0/3.0;
f[2] = f[0];

fx[0] = -0.5;
fx[1] = 0.0;
fx[2] = 0.5;

fxx[0] = 1.0;
fxx[1] = -2.0;
fxx[2] = 1.0;

Tensor_333(f, f, f, TensorF);

Tensor_333(fx, f, f, TensorFx);
Tensor_333(f, fx, f, TensorFy);
Tensor_333(f, f, fx, TensorFz);

Tensor_333(fxx,f, f, TensorFxx);
Tensor_333(fx, fx,f, TensorFxy);
Tensor_333(fx, f, fx,TensorFxz);

Tensor_333(f, fxx,f,  TensorFyy);
Tensor_333(f, fx, fx, TensorFyz);
Tensor_333(f, f,  fxx,TensorFzz);
}

void   HOSegmentation::Tensor_333(float *xx, float *yy, float *zz, float *result)
{
int i, j, k, l;

l = 0;
for (i = 0; i < 3; i++) {
   for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
          result[l] = xx[i]*yy[j]*zz[k];
          l = l + 1;
      }
   }
}  
}
/*-----------------------------------------------------------------------------*/

// The cubic basis evaluated at point x
//Assuming x \in [-2, 2] range. returns 0 if outside.
float HOSegmentation::evalCubic(float x)
{
  float x2 = x*x;
  float x3 = x2*x;
  //float fx = (x<-2.0)? 0.0: ((x<=-1)? powf(2.0 + x, 3.0)/6.0: ((x<= 0.0)? (4.0 - 6.0*x2 - 3.0*x3)/6.0 : ((x<= 1.0)? (4.0 - 6.0*x2 + 3.0*x3)/6.0 : ((x<=2.0)?  powf(2.0 - x, 3.0)/6.0 : 0.0))));
  float fx;
  if(x >=-2.0 && x<-1.0)
    fx = powf(2.0 + x, 3.0)/6.0;
  else if (x>=-1.0 && x<0.0)  
    fx = (4.0 - 6.0*x2 - 3.0*x3)/6.0;
  else if (x>= 0.0 && x<1.0)
    fx = (4.0 - 6.0*x2 + 3.0*x3)/6.0;
  else if (x>=1.0 && x<2.0)
    fx = powf(2.0 - x, 3.0)/6.0;
  else
    fx = 0.0;
          
  return fx;
}

// First differntial of the cubic basis above
float HOSegmentation::evalCubic_Dx(float x)
{
  float x2 = x*x;
  //float fx = (x<-2.0)? 0.0: ((x<=-1)? powf(2.0 + x, 2.0)/2.0: ((x<= 0.0)? (-4.0*x - 3.0*x2)/2.0 : ((x<= 1.0)? (- 4.0*x + 3.0*x2)/2.0 : ((x<=2.0)? -powf(2.0 - x, 2.0)/2.0 : 0.0))));
  float fx;
  if(x >=-2.0 && x<-1.0)
    fx = powf(2.0 + x, 2.0)/2.0;
  else if (x>=-1.0 && x<0.0)  
    fx = (-4.0*x - 3.0*x2)/2.0;
  else if (x>= 0.0 && x<1.0)
    fx = (- 4.0*x + 3.0*x2)/2.0;
  else if (x>=1.0 && x<2.0)
    fx = -powf(2.0 - x, 2.0)/2.0;
  else
    fx = 0.0;  
  return fx;
}

// Second differential of the cubic basis
float HOSegmentation::evalCubic_Dxx(float x)
{
  //float fx = (x<-2.0)? 0.0: ((x<=-1)? 2.0 + x : ((x<= 0.0)? -2.0 - 3.0*x : ((x<= 1.0)? -2.0 + 3.0*x : ((x<=2.0)?  2.0 - x : 0.0))));
  float fx;  
  if(x >=-2.0 && x<-1.0)
    fx = 2.0 + x;
  else if (x>=-1.0 && x<0.0)  
    fx = -2.0 - 3.0*x;
  else if (x>= 0.0 && x<1.0)
    fx = -2.0 + 3.0*x;
  else if (x>=1.0 && x<2.0)
    fx = 2.0 - x;
  else
    fx = 0.0;    
  return fx;
}

float HOSegmentation::evalTriCubic(float x, float y, float z)
{
  return evalCubic(x)*evalCubic(y)*evalCubic(z);
}

float  HOSegmentation::evalTriCubic_Dx(float x, float y, float z)
{
  return evalCubic_Dx(x)*evalCubic(y)*evalCubic(z);
}

float HOSegmentation::evalTriCubic_Dy(float x, float y, float z)
{
  return evalTriCubic_Dx(y, z, x);
}
  
float HOSegmentation::evalTriCubic_Dz(float x, float y, float z)
{
  return  evalTriCubic_Dx(z, x, y);
}

float HOSegmentation::evalTriCubic_Dxx(float x, float y, float z)
{
  return evalCubic_Dxx(x)*evalCubic(y)*evalCubic(z);
}

float HOSegmentation::evalTriCubic_Dyy(float x, float y, float z) 
{
  return evalTriCubic_Dxx(y, z, x);
}

float HOSegmentation::evalTriCubic_Dzz(float x, float y, float z)
{
  return evalTriCubic_Dxx(z, x, y);
}

float HOSegmentation::evalTriCubic_Dxy(float x, float y, float z)
{
  return evalCubic_Dx(x)*evalCubic_Dx(y)*evalCubic(z);
}

float HOSegmentation::evalTriCubic_Dyz(float x, float y, float z)
{
  return evalTriCubic_Dxy(y, z, x);
}

float HOSegmentation::evalTriCubic_Dzx(float x, float y, float z)
{
  return evalTriCubic_Dxy(z, x, y);
}

bool HOSegmentation::initInterfaceMulti()
 {
    return MSLevelSet::initInterfaceMulti();
 }
 
bool HOSegmentation::computeSDTEikonal()
{
    return MSLevelSet::computeSDTEikonal();
}

bool HOSegmentation::computeAverageIntensities()
{
  return MSLevelSet::computeAverageIntensities();
}

bool HOSegmentation::normalizeVolume()
{
  return MSLevelSet::normalizeVolume();
}
