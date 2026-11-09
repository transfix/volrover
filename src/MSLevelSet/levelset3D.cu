//============================================================================
/* Copyright (c) The Technical University of Denmark
 * Author: Ojaswa Sharma
 * E-mail: os@imm.dtu.dk
 * File: levelset3D.cu
 * Main levelset implementation and host functions to the CUDA device functions
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

#include <MSLevelSet/levelset3D.h>
#include <MSLevelSet/private/levelset3D_global.cu.h>
#include <MSLevelSet/private/init_interface_kernels.cu.h>
#include <MSLevelSet/private/signed_distance_kernels.cu.h>
#include <MSLevelSet/private/average_intensity_kernels.cu.h>
#include <MSLevelSet/private/pde_update_kernels.cu.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <cutil_inline.h>

using namespace MumfordShahLevelSet;

double MSLevelSet::getTime()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (double)tv.tv_sec + (double)tv.tv_usec / 1000000.0;
}

MSLevelSet::MSLevelSet(bool doInitCUDA)
{
  subvolDim = 128;
  cudaDevice = 0;
  h_vol = NULL; // The full volume
  h_phi = NULL; // The phi function 
  h_buffer[0] = NULL;
  h_buffer[1] = NULL;
  currentPhi = 0;
  h_subvolSpare = NULL; //Spare subvolume on the host for ntermediate computations. N.B: Size: subvolDim^3
  d_phiArray = 0; // subvolume bound to texture - used for texture reads in phi
  d_volArray = 0; // subvolume bound to texture - used for texture reads in intensity
  tPitch = 0;
  timer = 0;
  milliseconds = 0.0;
  mu = 0.0005f*255.0f*255.0f;
  h = 1.0f;
  nu = 0.0f;
  lambda1 = 1.0f;
  lambda2 = 1.0f;
  delta_t = 0.75f;
  nIter = 10; DTWidth = 10;
  converged = false;
  epsilon = 1.0f;
  iter_num = 0;
  BBoxOffset = 5;
  DoInitializeCUDA = doInitCUDA;
  multi_init_r = 10;
  multi_init_dr = 5;
  multi_init_s = 5;
  reinit_niter = 8;
}

int MSLevelSet::iDivUp(int a, int b)
{
  return ((a % b) != 0)? (a / b + 1): (a / b);
}

unsigned long MSLevelSet::inKB(unsigned long bytes)
{ return bytes/1024; }

unsigned long MSLevelSet::inMB(unsigned long bytes)
{ return bytes/(1024*1024); }

void MSLevelSet::printStats(unsigned long freemem, unsigned long total)
{
  debugPrint("GPU mem: total: %lu MB, free: %lu MB, %.1f%%\n", inMB(total), inMB(freemem), 100.0*(double)freemem/(double)total);
}

void MSLevelSet::printMemInfo()
{
  size_t freemem, total;
  //int gpuCount, i;
  CUresult res;
  CUdevice dev;
  CUcontext ctx;

  //cuInit(0);
  //cuDeviceGetCount(&gpuCount);
  //printf("Detected %d GPU\n",gpuCount);

  //for (i=0; i<gpuCount; i++)
  //{
  //cuDeviceGet(&dev,i);
  cuDeviceGet(&dev,cudaDevice);
  cuCtxCreate(&ctx, 0, dev);
  res = cuMemGetInfo(&freemem, &total);
  if(res != CUDA_SUCCESS)
      debugPrint("Error! cuMemGetInfo failed! (status = %x)\n", res);
 //debugPrint("^^^^ Device: %d\n",i);
  printStats(freemem, total);
  cuCtxDetach(ctx);
  //}
}

// Host - host mem, Mem - device mem, Array - device array (texture bound)
void MSLevelSet::copy3DHostToArray(float *_src, cudaArray *_dst, cudaExtent copy_extent, cudaExtent src_extent, cudaPos src_offset)
{
  cudaMemcpy3DParms copyParams = {0};
  float *h_source = _src + src_offset.x + src_offset.y*volSize.width + src_offset.z*volSize.width*volSize.height;
  copyParams.srcPtr = make_cudaPitchedPtr((void*)h_source, src_extent.width*sizeof(float), src_extent.width, src_extent.height);
  copyParams.dstArray = _dst;
  copyParams.kind = cudaMemcpyHostToDevice;
  //copyParams.extent = make_cudaExtent(subvolDim, subvolDim, subvolDim);
  copyParams.extent = copy_extent;

  CUDA_SAFE_CALL(cudaMemcpy3D(&copyParams));
  CUT_CHECK_ERROR("Host -> Array Memcpy failed\n");
}

void MSLevelSet::copy3DHostToMem(float *_src, cudaPitchedPtr _dst, cudaExtent copy_extent, cudaExtent src_extent, cudaPos src_offset)
{
  cudaMemcpy3DParms copyParams = {0};
  float *h_source = _src + src_offset.x + src_offset.y*volSize.width + src_offset.z*volSize.width*volSize.height;
  copyParams.srcPtr = make_cudaPitchedPtr((void*)h_source, src_extent.width*sizeof(float), src_extent.width, src_extent.height);
  copyParams.dstPtr = _dst;
  copyParams.kind = cudaMemcpyHostToDevice;
  //copyParams.extent = make_cudaExtent(subvolDim, subvolDim, subvolDim);
  copyParams.extent =  make_cudaExtent(copy_extent.width*sizeof(float), copy_extent.height, copy_extent.depth);

  CUDA_SAFE_CALL(cudaMemcpy3D(&copyParams));
  CUT_CHECK_ERROR("Host -> Mem Memcpy failed\n");
}

void MSLevelSet::copy3DMemToArray(cudaPitchedPtr _src, cudaArray *_dst)
{
  cudaMemcpy3DParms copyParams = {0};
  copyParams.srcPtr = 	_src;
  copyParams.dstArray = _dst;
  copyParams.kind = cudaMemcpyDeviceToDevice;
  copyParams.extent = make_cudaExtent(subvolDim, subvolDim, subvolDim);

  CUDA_SAFE_CALL(cudaMemcpy3D(&copyParams));
  CUT_CHECK_ERROR("Mem -> Array Memcpy failed\n");
}

void MSLevelSet::copy3DMemToHost(cudaPitchedPtr _src, float *_dst, cudaExtent copy_extent, cudaExtent dst_extent, cudaPos src_offset, cudaPos dst_offset)
{
  cudaMemcpy3DParms copyParams = {0};
  copyParams.srcPtr = _src;
  float *h_target = _dst + dst_offset.x + dst_offset.y*dst_extent.width + dst_offset.z*dst_extent.width*dst_extent.height;//For some reason, using copyParams.dstPos doesn't give correct results, so we set the offset here.
  copyParams.dstPtr = make_cudaPitchedPtr((void*)h_target, dst_extent.width*sizeof(float), dst_extent.width, dst_extent.height);
  copyParams.kind = cudaMemcpyDeviceToHost;
  copyParams.extent = make_cudaExtent(copy_extent.width*sizeof(float), copy_extent.height, copy_extent.depth);
  copyParams.srcPos = make_cudaPos(src_offset.x*sizeof(float), src_offset.y, src_offset.z); // We want to copy copy_extent sized volume starting at (x_off, y_off, z_off).

  CUDA_SAFE_CALL(cudaMemcpy3D(&copyParams));
  CUT_CHECK_ERROR("Mem -> Host Memcpy failed\n");
}

void MSLevelSet::copy3DArrayToHost(cudaArray *_src, float *_dst, cudaExtent copy_extent, cudaExtent dst_extent, cudaPos src_offset, cudaPos dst_offset)
{
  cudaMemcpy3DParms copyParams = {0};
  copyParams.srcArray = _src;
  float *h_target = _dst + dst_offset.x + dst_offset.y*dst_extent.width + dst_offset.z*dst_extent.width*dst_extent.height;//For some reason, using copyParams.dstPos doesn't give correct results, so we set the offset here.
  copyParams.dstPtr = make_cudaPitchedPtr((void*)h_target, dst_extent.width*sizeof(float), dst_extent.width, dst_extent.height);
  copyParams.kind = cudaMemcpyDeviceToHost;
  copyParams.extent = make_cudaExtent(copy_extent.width, copy_extent.height, copy_extent.depth);
  copyParams.srcPos = make_cudaPos(src_offset.x, src_offset.y, src_offset.z); // We want to copy copy_extent sized volume starting at (x_off, y_off, z_off).

  CUDA_SAFE_CALL(cudaMemcpy3D(&copyParams));
  CUT_CHECK_ERROR("Mem -> Host Memcpy failed\n");
}

bool MSLevelSet::initCuda()
{
  //create 3D array
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();

  //printMemInfo();
  debugPrint("Allocating CUDA arrays...\n");
  //Allocate mem and bind phiArray to 3D texture
  cutilSafeCall(cudaMalloc3DArray(&d_phiArray, &channelDesc, subvolSize)); 
  CHK_CUDA_ERR("\n");
  texPhi.normalized = false;
  texPhi.filterMode = cudaFilterModePoint;
  texPhi.addressMode[0] = cudaAddressModeWrap; 
  texPhi.addressMode[1] = cudaAddressModeWrap; 
  texPhi.addressMode[2] = cudaAddressModeWrap; 
  CUDA_SAFE_CALL(cudaBindTextureToArray(texPhi, d_phiArray, channelDesc));
  CHK_CUDA_ERR("\n");

  //Allocate mem and bind phiArray to 3D texture
  CUDA_SAFE_CALL(cudaMalloc3DArray(&d_volArray, &channelDesc, subvolSize));
  CHK_CUDA_ERR("\n");
  texVol.normalized = false;
  texVol.filterMode = cudaFilterModePoint;
  texVol.addressMode[0] = cudaAddressModeWrap; 
  texVol.addressMode[1] = cudaAddressModeWrap; 
  texVol.addressMode[2] = cudaAddressModeWrap; 
  CUDA_SAFE_CALL(cudaBindTextureToArray(texVol, d_volArray, channelDesc));
  CHK_CUDA_ERR("\n");

  //Allocate memory for 3D device memory for kernel output
  cudaExtent pitchedVolSize = make_cudaExtent(subvolDim*sizeof(float), subvolDim, subvolDim); 
  CUDA_SAFE_CALL(cudaMalloc3D(&d_volPPtr, pitchedVolSize));
  CHK_CUDA_ERR("\n");

  //Allocate memory for spare subvolumes
  CUDA_SAFE_CALL(cudaMalloc3D(&d_spare1PPtr, pitchedVolSize));
  CHK_CUDA_ERR("\n");
  CUDA_SAFE_CALL(cudaMalloc3D(&d_spare2PPtr, pitchedVolSize));
  CHK_CUDA_ERR("\n");
  CUDA_SAFE_CALL(cudaMalloc3D(&d_spare3PPtr, pitchedVolSize));
  CHK_CUDA_ERR("\n");
  
  printMemInfo();

  return true;
}

void MSLevelSet::writeRawData(float *data, char* filename)
{
  //cudaExtent size = make_cudaExtent(subvolDim-2, subvolDim-2, subvolDim-2);
  cudaExtent size = subvolSize;
  
  FILE *fp = fopen(filename, "w");
  if(!fp) {
    debugPrint("Error opening file '%s'\n", filename);
    return;
  }

  for(int _z = 0; _z<size.depth; _z++)
  {
    for(int _y = 0; _y<size.height; _y++)
    {
      for(int _x = 0; _x<size.width; _x++)
          fprintf(fp, "%+07.3f ", *(data + _x + _y*size.width + _z*size.width*size.height));
	//fprintf(fp, "%c", *(data + _x + _y*size.width + _z*size.width*size.height) > 0.0 ? '+': '-');
      fprintf(fp, "\n");
    }
    fprintf(fp, "--(%d)--\n", _z+1);
  }
  fclose(fp);
}

void MSLevelSet::writeSlices(float *data, char* filename_prefix, int nslices)
{
  int delta_slices = iDivUp(volSize.depth, nslices);

  for(int _z = 0; _z<volSize.depth; _z+=delta_slices)
  {
    char fname[256];
    sprintf(fname, "./out/%s-%d.raw", filename_prefix, _z);
    FILE *fp = fopen(fname, "w");
    if(!fp) {
      debugPrint("Error opening file '%s'\n", fname);
      continue;
    }

    for(int _y = 0; _y<volSize.height; _y++)
    {
      for(int _x = 0; _x<volSize.width; _x++)
        fprintf(fp, "%6.3f ", *(data + _x + _y*volSize.width + _z*volSize.width*volSize.height));
      fprintf(fp, "\n");
    }
    fprintf(fp, "----\n");
    fclose(fp);
  }
}

void MSLevelSet::adjustDnloadSubvolSize(int _x, int _y, int _z, unsigned int offset, cudaExtent &copyvol_dnload) //Download to Host. offset is the shared border around subvolume.
{
  copyvol_dnload.width = subvolDim - offset;
  copyvol_dnload.height = subvolDim - offset;
  copyvol_dnload.depth = subvolDim - offset;
  if(_x == (subvolIndicesExtents.width-1))
  {
    copyvol_dnload.width = (volSize.width-offset) % (subvolDim-offset); 
    if(copyvol_dnload.width == 0) copyvol_dnload.width = subvolDim - offset;
  }
  if(_y == (subvolIndicesExtents.height-1))
  {
    copyvol_dnload.height = (volSize.height-offset) % (subvolDim-offset); 
    if(copyvol_dnload.height == 0) copyvol_dnload.height = subvolDim - offset;
  }
  if(_z == (subvolIndicesExtents.depth-1))
  {
    copyvol_dnload.depth = (volSize.depth-offset) % (subvolDim-offset); 
    if(copyvol_dnload.depth == 0) copyvol_dnload.depth = subvolDim - offset;
  }
}

void MSLevelSet::adjustUploadSubvolSize(int _x, int _y, int _z, unsigned int offset, cudaExtent &copyvol_upload) //upload to GPU. offset is the shared border around the subvolume.
{
  zeroOutArray = false;
  copyvol_upload.width = subvolDim;
  copyvol_upload.height = subvolDim;
  copyvol_upload.depth = subvolDim;

  if(_x == (subvolIndicesExtents.width-1))
  {
    copyvol_upload.width = (volSize.width-offset) % (subvolDim-offset); 
    if(copyvol_upload.width == 0)
      copyvol_upload.width = subvolDim;
    else {
      zeroOutArray = true;
      copyvol_upload.width += offset/2; //We need extra offset/2 voxel layer only on one side.
    }
  }
  if(_y == (subvolIndicesExtents.height-1))
  {
    copyvol_upload.height = (volSize.height-offset) % (subvolDim-offset); 
    if(copyvol_upload.height == 0)
      copyvol_upload.height = subvolDim;
    else {
      zeroOutArray = true;
      copyvol_upload.height += offset/2; //We need extra offset/2 voxel layer only on one side.
    }
  }
  if(_z == (subvolIndicesExtents.depth-1))
  {
    copyvol_upload.depth = (volSize.depth-offset) % (subvolDim-offset); 
    if(copyvol_upload.depth == 0)
      copyvol_upload.depth = subvolDim;
    else {
      zeroOutArray = true;
      copyvol_upload.depth += offset/2; //We need extra offset/2 voxel layer only on one side.
    }
  }
}

bool MSLevelSet::initInterfaceMulti()
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

  const dim3 SDTblockSize(SDTBlockDim, SDTBlockDim, SDTBlockDim);
  const dim3 SDTgridSize(iDivUp(subvolDim, SDTBlockDim)*iDivUp(subvolDim, SDTBlockDim)*iDivUp(subvolDim, SDTBlockDim), 1);//linearized grid, since, we cannot have a 3D grid
  cudaExtent logicalGridSize = make_cudaExtent(iDivUp(subvolDim, SDTBlockDim), iDivUp(subvolDim, SDTBlockDim), iDivUp(subvolDim, SDTBlockDim));
  char _msg[256];
  //The 1 voxel border is not initialized by our cuda kernel. Do it on the CPU :p. Outside is -ve
  for(int _x=0;_x<volSize.width; _x++) for(int _y=0; _y<volSize.height; _y++) { *(h_phi+ _x + _y*volSize.width) = maxDist; *(h_phi + _x + _y*volSize.width + (volSize.depth-1)*volSize.width*volSize.height) = maxDist; }
  for(int _x=0;_x<volSize.width; _x++) for(int _z=0; _z<volSize.depth; _z++) { *(h_phi+ _x + _z*volSize.width*volSize.height) = maxDist; *(h_phi + _x + (volSize.height-1)*volSize.width + _z*volSize.width*volSize.height) = maxDist; }
  for(int _z=0;_z<volSize.depth; _z++) for(int _y=0; _y<volSize.height; _y++) { *(h_phi + _y*volSize.width + _z*volSize.width*volSize.height) = maxDist; *(h_phi + volSize.width - 1 + _y*volSize.width + _z*volSize.width*volSize.height) = maxDist; }

  CUT_SAFE_CALL(cutStartTimer(timer));
  cudaExtent subvolIdx = make_cudaExtent(0, 0, 0);
  cudaExtent copyvol = make_cudaExtent(subvolDim-2, subvolDim-2, subvolDim-2);
  cudaPos src_off = make_cudaPos(1, 1, 1);
  cudaPos dst_off = make_cudaPos(0, 0, 0);
  float3 offset = make_float3(l_offx, l_offy, l_offz);
  float2 radii = make_float2(rad0, rad1);
  float3 upper_limit = make_float3(float(l_nx*l_box), float(l_ny*l_box), float(l_nz*l_box));


  for(int _z = 0; _z< subvolIndicesExtents.depth; _z++)
    for(int _y = 0; _y< subvolIndicesExtents.height; _y++)
      for(int _x = 0; _x< subvolIndicesExtents.width; _x++)
      {
        subvolIdx.width = _x; subvolIdx.height = _y; subvolIdx.depth = _z;
        d_init_interface_multi_smooth<<<SDTgridSize, SDTblockSize>>>(d_volPPtr, logicalGridSize, subvolIdx, maxDist, subvolDim, SDTBlockDim, offset, l_boxf, radius_reach, radii, upper_limit);	
        CUT_CHECK_ERROR("Kernel failed");
        sprintf(_msg,  "---%d-%d-%d", _x, _y, _z);
        CHK_CUDA_ERR(_msg);
        cudaThreadSynchronize();

        //copy the results to main volume
        //Watchout for boundaries at the end. We don't want to write to volume outside allocated array!!
        adjustDnloadSubvolSize(_x, _y, _z, 2, copyvol);
        dst_off.x = 1 + _x*(subvolDim-2); dst_off.y = 1 + _y*(subvolDim-2); dst_off.z = 1 + _z*(subvolDim-2);
        copy3DMemToHost(d_volPPtr, h_phi, copyvol, volSize, src_off, dst_off); 
        CHK_CUDA_ERR("->");
      }

  CUT_SAFE_CALL(cutStopTimer(timer));
  milliseconds = cutGetAverageTimerValue(timer);
  CUT_SAFE_CALL(cutResetTimer(timer));
  debugPrint(" %f ms\n", milliseconds);

  return true;
}

bool MSLevelSet::computeSDTEikonal()
{



  //For some reason the texture is to be bound again to the array.
  //Needed by MultiPhaseSeg - os/2010-10-10
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
  texPhi.normalized = false;
  texPhi.filterMode = cudaFilterModePoint;
  texPhi.addressMode[0] = cudaAddressModeWrap; 
  texPhi.addressMode[1] = cudaAddressModeWrap; 
  texPhi.addressMode[2] = cudaAddressModeWrap; 
  CUDA_SAFE_CALL(cudaBindTextureToArray(texPhi, d_phiArray, channelDesc));

  const dim3 SDTblockSize(SDTBlockDim, SDTBlockDim, SDTBlockDim);
  const dim3 SDTgridSize(iDivUp(subvolDim, SDTBlockDim)*iDivUp(subvolDim, SDTBlockDim)*iDivUp(subvolDim, SDTBlockDim), 1);//linearized grid, since, we cannot have a 3D grid
  cudaExtent logicalGridSize = make_cudaExtent(iDivUp(subvolDim, SDTBlockDim), iDivUp(subvolDim, SDTBlockDim), iDivUp(subvolDim, SDTBlockDim));
  cudaExtent subvolIndicesNewExtents = make_cudaExtent(iDivUp(volSize.width-2, subvolDim-2), iDivUp(volSize.height-2, subvolDim-2), iDivUp(volSize.depth-2, subvolDim-2));
  debugPrint("Computing signed distance field (Eikonal PDE)... ");
  cudaExtent copyvol = make_cudaExtent(subvolDim-2, subvolDim-2, subvolDim-2);
  cudaExtent copyvol_upload = make_cudaExtent(subvolDim, subvolDim, subvolDim);
  cudaPos src_off = make_cudaPos(1, 1, 1);
  cudaPos dst_off = make_cudaPos(0, 0, 0);

  float maxDist = sqrtf(float(volSize.width*volSize.width + volSize.height*volSize.height + volSize.depth*volSize.depth));
  if((reinit_niter%2) == 0) reinit_niter++; //So that the final written buffer is 1 for phi

  //We copy the boundary to the texture, so the extra volume uploaded must have the same value as the outside of the interface initialization.
  for (int _z = 0; _z<subvolDim; _z++)
    for (int _y = 0; _y<subvolDim; _y++)
      for (int _x = 0; _x<subvolDim; _x++)
        *(h_subvolSpare + _x + _y*subvolDim + _z*subvolDim*subvolDim) = -maxDist;	

  CUT_SAFE_CALL(cutStartTimer(timer));
  for(int j = 0; j < reinit_niter; j++) {
    debugPrint("*");
    for(int _z = 0; _z< subvolIndicesNewExtents.depth; _z++)
      for(int _y = 0; _y< subvolIndicesNewExtents.height; _y++)
        for(int _x = 0; _x< subvolIndicesNewExtents.width; _x++)
        {
          adjustUploadSubvolSize(_x, _y, _z, 2, copyvol_upload);
          //Copy the subvolume to texture
          // Watchout for boundaries at the end. We don't want to read garbage values outside the allocated host volume!!
          if(zeroOutArray) {
            //Fill the texture with -maxDist so that it doesn't have undesirable values when bundary volume parts are copied.
            copy3DHostToArray(h_subvolSpare, d_phiArray, subvolSize, subvolSize, make_cudaPos(0, 0, 0)); 
          }
          copy3DHostToArray(h_buffer[currentPhi], d_phiArray, copyvol_upload, volSize, make_cudaPos(_x*(subvolDim-2), _y*(subvolDim-2), _z*(subvolDim-2)));
          CHK_CUDA_ERR("->");
          d_signed_distance_eikonal<<<SDTgridSize, SDTblockSize>>>(d_volPPtr, logicalGridSize, subvolDim, SDTBlockDim);
          CUT_CHECK_ERROR("Kernel failed");
          CHK_CUDA_ERR("->");
          cudaThreadSynchronize();

          //Copy results back to host mem from device mem
          //Watchout for boundaries at the end. We don't want to write to volume outside allocated array!!
          adjustDnloadSubvolSize(_x, _y, _z, 2, copyvol);
          dst_off.x = 1 + _x*(subvolDim-2); dst_off.y = 1 + _y*(subvolDim-2); dst_off.z = 1 + _z*(subvolDim-2);
          copy3DMemToHost(d_volPPtr, h_buffer[(currentPhi+1)%2], copyvol, volSize, src_off, dst_off); 
          CHK_CUDA_ERR("%s<-\n");
        }
    currentPhi = (currentPhi+1)%2;
  }
  h_phi = h_buffer[currentPhi];
  CUT_SAFE_CALL(cutStopTimer(timer));
  milliseconds = cutGetAverageTimerValue(timer);
  debugPrint(" %f ms\n", milliseconds);
  CUT_SAFE_CALL(cutResetTimer(timer));


  return true;
}

//TODO: Corect CUDA kernel for average intensities and it's correspondng function.
bool MSLevelSet::computeAverageIntensities()
{
  ///Compute average intensity value inside and outside the volume
  //Warning: d_spare3PPtr has been used in initializing textures in loops. Do not use later them later in the loop or modifiy accordingly.
  debugPrint("Computing average intensities... ");
  float c1_d = 0.0, c1_n = 0.0, c2_d = 0.0, c2_n = 0.0; 
  unsigned int nvoxels = datainfo.n_input[0]*datainfo.n_input[1]*datainfo.n_input[2];

#if COMPUTE_AVERAGE_ON_CPU
  clock_t t0, t1;
  t0 = clock();
  float l_phi = 0.0;
  float l_I = 0.0;
  float l_H = 0.0;
  h_phi = h_buffer[currentPhi];


  for(int z = 0; z< (volSize.depth); z++)
    for(int y = 0; y< (volSize.height); y++)
      for(int x = 0; x< (volSize.width); x++)
      {
        l_phi = *(h_phi + x + y*volSize.width + z*volSize.width*volSize.height);
        l_I = *(h_vol + x + y*volSize.width + z*volSize.width*volSize.height);
        l_H = 0.5*(1.0 + (2.0/PI)*atanf(l_phi/epsilon));
        c1_d += l_H;	
        c1_n += l_H*l_I;
        c2_n += (1.0 - l_H)*l_I;
     	if(isnan(c2_n) || isnan(c1_n))
	    debugPrint("after l-H=%f, l_I=%f,  c1_d=%f, c1n=%f, c2n=%f\n", l_H, l_I, c1_d, c1_n, c2_n); 
      }
  t1 = clock();
  milliseconds = (float)(t1-t0)*1000.0/CLOCKS_PER_SEC;
#else
  //debugPrint("---%f, %f, %f---\n", c1_d, c1_n, c2_n);
  //c1_d = 0.0; c1_n = 0.0; c2_n = 0.0;
  const dim3 accumblockSize(avgBlockDim, avgBlockDim, avgBlockDim);
  const dim3 accumgridSize(iDivUp(subvolDim, avgBlockDim)*iDivUp(subvolDim, avgBlockDim)*iDivUp(subvolDim, avgBlockDim), 1);//linearized grid, since, we cannot have a 3D grid
  cudaExtent logicalGridSize = make_cudaExtent(iDivUp(subvolDim, avgBlockDim), iDivUp(subvolDim, avgBlockDim), iDivUp(subvolDim, avgBlockDim));
  cudaExtent copyvol_upload = make_cudaExtent(subvolDim, subvolDim, subvolDim);
  cudaExtent subvolIndicesNewExtents = make_cudaExtent(iDivUp(volSize.width, subvolDim), iDivUp(volSize.height, subvolDim), iDivUp(volSize.depth, subvolDim));
  char _msg[256];
  float maxDist = (float)DTWidth + 1.0;
  
  //We copy the boundary to the texture, so the extra volume uploaded must have the same value as the outside of the interface initialization.
  for (int _z = 0; _z<subvolDim; _z++)
    for (int _y = 0; _y<subvolDim; _y++)
      for (int _x = 0; _x<subvolDim; _x++)
	*(h_subvolSpare + _x + _y*subvolDim + _z*subvolDim*subvolDim) = -maxDist;	
  
  //Zero out the accumulators first
  cudaMemset3D(d_volPPtr, 0, make_cudaExtent(subvolDim*sizeof(float), subvolDim, subvolDim));// zero out the device mem
  cudaMemset3D(d_spare1PPtr, 0, make_cudaExtent(subvolDim*sizeof(float), subvolDim, subvolDim));// zero out the device mem
  cudaMemset3D(d_spare2PPtr, 0, make_cudaExtent(subvolDim*sizeof(float), subvolDim, subvolDim));// zero out the device mem

  CUT_SAFE_CALL(cutStartTimer(timer));
  for(int _z = 0; _z< subvolIndicesNewExtents.depth; _z++)
    for(int _y = 0; _y< subvolIndicesNewExtents.height; _y++)
      for(int _x = 0; _x< subvolIndicesNewExtents.width; _x++)
      {
        adjustUploadSubvolSize(_x, _y, _z, 0, copyvol_upload);
        //load texture values for phi and intensity values
        // Watchout for boundaries at the end. We don't want to read garbage values outside the allocated host volume!!
        if(zeroOutArray)
        {
          //zero out Array first so they do not have undesired values stored in them from previous load in the extraneous region.
          //N.B.: there seems to be no API to memset a cudaArray, so we do it indirectly
          cudaMemset3D(d_spare3PPtr, 0, make_cudaExtent(subvolDim*sizeof(float), subvolDim, subvolDim));// zero out the device mem
          //copy3DMemToArray(d_spare3PPtr, d_phiArray); 
          copy3DMemToArray(d_spare3PPtr, d_volArray); 
          copy3DHostToArray(h_subvolSpare, d_phiArray, subvolSize, subvolSize, make_cudaPos(0, 0, 0)); //<-- Use ths instead of copy3DMemtoArray since we need -maxDist instead zeros in H!!!
        }

        copy3DHostToArray(h_phi, d_phiArray, copyvol_upload, volSize, make_cudaPos(_x*subvolDim, _y*subvolDim, _z*subvolDim));
        CHK_CUDA_ERR("->");
        copy3DHostToArray(h_vol, d_volArray, copyvol_upload, volSize, make_cudaPos(_x*subvolDim, _y*subvolDim, _z*subvolDim));
        CHK_CUDA_ERR("->");
        //Accumulate H_in, I*H_in, and I*H_out to d_volPPtr, d_spare1PPtr, and d_spare2PPtr respectively
        accumulate_average_intensities<<<accumgridSize, accumblockSize>>>(d_volPPtr, d_spare1PPtr, d_spare2PPtr, logicalGridSize, epsilon, (float)nvoxels, subvolDim, avgBlockDim);
        CUT_CHECK_ERROR("Kernel failed");
        sprintf(_msg,  "---%d-%d-%d", _x, _y, _z);
        CHK_CUDA_ERR(_msg);
        cudaThreadSynchronize();
      }

  //Sum up the accumulators to get the final values of the average intensities. NB: This is done on the CPU (subvolDim^3 loop)
  cudaExtent copyvol = make_cudaExtent(subvolDim, subvolDim, subvolDim);
  cudaExtent dstextent = make_cudaExtent(subvolDim, subvolDim, subvolDim);
  cudaPos src_off = make_cudaPos(0, 0, 0);
  cudaPos dst_off = make_cudaPos(0, 0, 0);

  copy3DMemToHost(d_volPPtr, h_subvolSpare, copyvol, dstextent, src_off, dst_off); 
  CHK_CUDA_ERR("<-%s");
  for(int _z=0; _z<subvolDim; _z++) 
    for(int _y=0; _y<subvolDim; _y++) 
      for(int _x=0; _x<subvolDim; _x++) 
	c1_d += *(h_subvolSpare + _x + _y*subvolDim + _z*subvolDim*subvolDim);	

  copy3DMemToHost(d_spare1PPtr, h_subvolSpare, copyvol, dstextent, src_off, dst_off); 
  CHK_CUDA_ERR("<-%s");
  for(int _z=0; _z<subvolDim; _z++) 
    for(int _y=0; _y<subvolDim; _y++) 
      for(int _x=0; _x<subvolDim; _x++) 
	c1_n += *(h_subvolSpare + _x + _y*subvolDim + _z*subvolDim*subvolDim);	

  copy3DMemToHost(d_spare2PPtr, h_subvolSpare, copyvol, dstextent, src_off, dst_off); 
  CHK_CUDA_ERR("<-%s");
  for(int _z=0; _z<subvolDim; _z++) 
    for(int _y=0; _y<subvolDim; _y++) 
      for(int _x=0; _x<subvolDim; _x++) 
	c2_n += *(h_subvolSpare + _x + _y*subvolDim + _z*subvolDim*subvolDim);	

  CUT_SAFE_CALL(cutStopTimer(timer));
  milliseconds = cutGetAverageTimerValue(timer);
  CUT_SAFE_CALL(cutResetTimer(timer));
#endif

  debugPrint("---%f, %f, %f---\n", c1_d, c1_n, c2_n);
  //Compute the averages...
  c2_d = (float)nvoxels - c1_d;
  if(c1_d == 0.0f)
    c1 = 0.0f;
  else
    c1 = c1_n/c1_d;

  if(c2_d == 0.0f)
    c2 = 0.0f;
  else
    c2 = c2_n/c2_d;

  debugPrint("<(%d), %16.10f, %16.10f>", nvoxels, c1, c2);
  debugPrint(" %f ms\n", milliseconds);
  
  return true;
}

bool MSLevelSet::PDEUpdate()
{
  ///Compute external and internal energy and update phi and vol
  //Warning: h_subvolSpare, d_volPPtr have been used in initializing textures in loops. Do not use later them later in the loop or modifiy accordingly.
  debugPrint("Updating the PDE (level set and input volume)... ");
  float mu_factor = mu/(h*h);
  h_phi = h_buffer[currentPhi];
  double t0, t1;
  float grad_left[3], grad_right[3];
  float l_phi, l_I;
  float grad_norm, delta_dirac;
  float curvature_energy, external_energy;
  unsigned long i = 0;
  unsigned long nxy = volSize.width*volSize.height;
  unsigned long nelem  = nxy*volSize.depth;
  unsigned int x, y, z;
#if UPDATE_PDE_ON_CPU
  t0 = getTime();
#ifdef _OPENMP
#pragma omp parallel for private(i, x, y, z, l_phi, curvature_energy, grad_left, grad_norm, grad_right, l_I, external_energy, delta_dirac)
#endif
  for(i=0; i< nelem; i++) {
    z = i/nxy;
    y = (i - z*nxy)/volSize.width; 
    x = (i - z*nxy - y*volSize.width);
    l_phi = VOXEL_AT(h_phi, x, y, z);
    // Parabolic term
    curvature_energy = 0.;
    if(x > 1 && x < (volSize.width-2)) {
      grad_left[0] = (l_phi - VOXEL_AT(h_phi, x - 2, y, z))/2.0;
      grad_right[0] = (VOXEL_AT(h_phi, x + 2, y, z) - l_phi)/2.0;
    }

    if(y > 1 && y < (volSize.height-2)) {
      grad_left[1] = (l_phi - VOXEL_AT(h_phi, x, y - 2, z))/2.0;
      grad_right[1] = (VOXEL_AT(h_phi, x, y + 2, z) - l_phi)/2.0;
    } 

    if(z > 1 && z < (volSize.depth-2)) {
      grad_left[2] = (l_phi - VOXEL_AT(h_phi, x, y, z - 2))/2.0;
      grad_right[2] = (VOXEL_AT(h_phi, x, y, z + 2) - l_phi)/2.0;
    }
 
    grad_norm = sqrt(grad_left[0]*grad_left[0] + grad_left[1]*grad_left[1] + grad_left[2]*grad_left[2] + TINY);
    grad_left[0] /= grad_norm; grad_left[1] /= grad_norm; grad_left[2] /= grad_norm;

    grad_norm = sqrt(grad_right[0]*grad_right[0] + grad_right[1]*grad_right[1] + grad_right[2]*grad_right[2] + TINY);
    grad_right[0] /= grad_norm; grad_right[1] /= grad_norm; grad_right[2] /= grad_norm;

    curvature_energy = (grad_right[0] - grad_left[0] + grad_right[1] - grad_left[1] + grad_right[2] - grad_left[2])/2.0;
    curvature_energy *= mu_factor;

    //External energy
    l_I = VOXEL_AT(h_vol, x, y, z);
    external_energy = -nu - lambda1*(l_I - c1)*(l_I - c1) + lambda2*(l_I - c2)*(l_I - c2);

    //PDE update
    delta_dirac = (epsilon/M_PI)/(epsilon*epsilon + l_phi*l_phi);
    VOXEL_AT(h_buffer[(currentPhi+1)%2], x, y, z) = l_phi + delta_t*delta_dirac*(curvature_energy + external_energy);
  }
  t1 = getTime();
  milliseconds = (t1 - t0)*1000.0;
#else
  const dim3 PDEblockSize(PDEBlockDim, PDEBlockDim, PDEBlockDim);
  const dim3 PDEgridSize(iDivUp(subvolDim, PDEBlockDim)*iDivUp(subvolDim, PDEBlockDim)*iDivUp(subvolDim, PDEBlockDim), 1);//linearized grid, since, we cannot have a 3D grid
  cudaExtent logicalGridSize = make_cudaExtent(iDivUp(subvolDim, PDEBlockDim), iDivUp(subvolDim, PDEBlockDim), iDivUp(subvolDim, PDEBlockDim));
  cudaExtent subvolIndicesNewExtents = make_cudaExtent(iDivUp(volSize.width-4, subvolDim-4), iDivUp(volSize.height-4, subvolDim-4), iDivUp(volSize.depth-4, subvolDim-4));
  cudaExtent copyvol = make_cudaExtent(subvolDim-4, subvolDim-4, subvolDim-4);
  cudaPos src_off = make_cudaPos(2, 2, 2);
  cudaPos dst_off = make_cudaPos(0, 0, 0);
  char _msg[256];
  float maxDist = sqrtf(float(volSize.width*volSize.width + volSize.height*volSize.height + volSize.depth*volSize.depth));

  //We copy the boundary to the texture, so the extra volume uploaded must have the same value as the outside of the interface initialization.
  for (int _z = 0; _z<subvolDim; _z++)
    for (int _y = 0; _y<subvolDim; _y++)
      for (int _x = 0; _x<subvolDim; _x++)
        *(h_subvolSpare + _x + _y*subvolDim + _z*subvolDim*subvolDim) = -maxDist;	

  cudaExtent copyvol_upload = make_cudaExtent(0,0,0);
  CUT_SAFE_CALL(cutStartTimer(timer));
  for(int _z = 0; _z< subvolIndicesNewExtents.depth; _z++)
    for(int _y = 0; _y< subvolIndicesNewExtents.height; _y++)
      for(int _x = 0; _x< subvolIndicesNewExtents.width; _x++)
      {
        adjustUploadSubvolSize(_x, _y, _z, 4, copyvol_upload);
        //load texture values for phi and intensity values
        // Watchout for boundaries at the end. We don't want to read garbage values outside the allocated host volume!!
        if(zeroOutArray)
        {
          //zero out Arrays first so they do not have undesired values stored in them from previous load
          //N.B.: there seems to be no API to memset a cudaArray, so we do it indirectly
          cudaMemset3D(d_volPPtr, 0, make_cudaExtent(subvolDim*sizeof(float), subvolDim, subvolDim));// zero out the device mem
          copy3DMemToArray(d_volPPtr, d_volArray); 
          //Fill the phi texture with -maxDist so that it doesn't have undesirable values when bundary volume parts are copied.
          copy3DHostToArray(h_subvolSpare, d_phiArray, subvolSize, subvolSize, make_cudaPos(0, 0, 0)); 
        }

        copy3DHostToArray(h_phi, d_phiArray, copyvol_upload, volSize, make_cudaPos(_x*(subvolDim-4), _y*(subvolDim-4), _z*(subvolDim-4)));
        CHK_CUDA_ERR("->");
        copy3DHostToArray(h_vol, d_volArray, copyvol_upload, volSize, make_cudaPos(_x*(subvolDim-4), _y*(subvolDim-4), _z*(subvolDim-4)));
        CHK_CUDA_ERR("->");

        //Update the PDE
        PDE_update<<<PDEgridSize, PDEblockSize>>>(d_volPPtr, logicalGridSize, mu_factor, nu, make_float2(lambda1, lambda2), make_float2(c1, c2), epsilon, delta_t, subvolDim, PDEBlockDim);
        CUT_CHECK_ERROR("Kernel failed");
        sprintf(_msg,  "---%d-%d-%d", _x, _y, _z);
        CHK_CUDA_ERR(_msg);     	
        cudaThreadSynchronize();

        // Copy results back to h_phi
        // Watchout for boundaries at the end. We don't want to write to volume outside allocated array!!
        adjustDnloadSubvolSize(_x, _y, _z, 4, copyvol);
        dst_off.x = 2+ _x*(subvolDim-4); dst_off.y = 2 + _y*(subvolDim-4); dst_off.z = 2 + _z*(subvolDim-4);
        copy3DMemToHost(d_volPPtr, h_buffer[(currentPhi+1)%2], copyvol, volSize, src_off, dst_off); //N.B.: Results written to h_buffer+1 and NOT to h_phi
        CHK_CUDA_ERR("<-%s");
      }
  CUT_SAFE_CALL(cutStopTimer(timer));
  milliseconds = cutGetAverageTimerValue(timer);
  CUT_SAFE_CALL(cutResetTimer(timer));
#endif
  //The computation above is correct for the volume of size [2, Nx-1]x[2, Ny-1]x[2, Nz-1]. Boundry correction is done as follows
  //Compute correct PDE evolution on the boundary of the full volume. It is a good idea to do it hear for computational speed! (fewer "if" conditions in kernels)
  debugPrint(" Boundary correction...");
  t0 = getTime();
#ifdef _OPENMP
#pragma omp parallel for private(i, x, y, z, l_phi, curvature_energy, grad_left, grad_norm, grad_right, l_I, external_energy, delta_dirac)
#endif
  for(i=0; i< nelem; i++) {
    z = i/nxy;
    y = (i - z*nxy)/volSize.width; 
    x = (i - z*nxy - y*volSize.width);
    if(x > 1 && x < (volSize.width-2) && y > 1 && y < (volSize.height-2) && z > 1 && z < (volSize.depth-2)) continue; //Early reject
    l_phi = VOXEL_AT(h_phi, x, y, z);
    // Parabolic term
    curvature_energy = 0.;
    if (x == 1) {
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

    if (y == 1) {
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

    if (z == 1) {
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

    curvature_energy = (grad_right[0] - grad_left[0] + grad_right[1] - grad_left[1] + grad_right[2] - grad_left[2])/2.0;
    curvature_energy *= mu_factor;

    //External energy
    l_I = VOXEL_AT(h_vol, x, y, z);
    external_energy = -nu - lambda1*(l_I - c1)*(l_I - c1) + lambda2*(l_I - c2)*(l_I - c2);

    //PDE update
    delta_dirac = (epsilon/M_PI)/(epsilon*epsilon + l_phi*l_phi);
    VOXEL_AT(h_buffer[(currentPhi+1)%2], x, y, z) = l_phi + delta_t*delta_dirac*(curvature_energy + external_energy);
  }
  //Point h_phi to correct phi array
  currentPhi = (currentPhi+1)%2;
  h_phi = h_buffer[currentPhi];
  t1 = getTime();
  milliseconds += (t1 - t0)*1000.0;

  debugPrint(" %f ms\n", milliseconds);
  return true;
}

bool MSLevelSet::solve()
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
    if(!PDEUpdate()) return false; //update phi and vol via level set PDE's
  }
  debugPrint("Solver done.\n");
  return true;
}

int MSLevelSet::solverMain() 
{
  if(DoInitializeCUDA)
  {
  cuInit(0); //Any driver API can also be called now on...
  cudaSetDevice(cutGetMaxGflopsDeviceId());
  }

  timer = 0;
  cutilCheckError(cutCreateTimer(&timer));
  cutilCheckError(cutResetTimer(timer));
  
  // set the volume size and number of subvolumes
  volSize = make_cudaExtent(datainfo.n[0], datainfo.n[1], datainfo.n[2]);
  subvolIndicesExtents = make_cudaExtent(iDivUp(volSize.width-2, subvolDim-2), iDivUp(volSize.height-2, subvolDim-2), iDivUp(volSize.depth-2, subvolDim-2));
  subvolSize = make_cudaExtent(subvolDim, subvolDim, subvolDim);

  //TODO: h_vol must be in the range [0..255]. If not, scale the intensities.

  // Allocate two buffers for host memory for phi (Signed Distance function)
  h_buffer[0] = (float*)malloc(volSize.width*volSize.height*volSize.depth*sizeof(float));
  if(h_buffer[0] == NULL)
  {
    debugPrint("Failed to allocate %lu bytes of memory\n", volSize.width*volSize.height*volSize.depth*sizeof(float));
    freeHostMem();
    return false;
  } 
  h_buffer[1] = (float*)malloc(volSize.width*volSize.height*volSize.depth*sizeof(float));
  if(h_buffer[1] == NULL)
  {
    debugPrint("Failed to allocate %lu bytes of memory\n", volSize.width*volSize.height*volSize.depth*sizeof(float));
    freeHostMem();
    return false;
  } 
  h_phi = h_buffer[0];
  currentPhi = 0;

  //Allocate memory for the spare subvolume on the host
  h_subvolSpare = (float*)malloc((subvolDim)*(subvolDim)*(subvolDim)*sizeof(float));
  if(h_subvolSpare == NULL)
  {
    debugPrint("Failed to allocate %lu bytes of memory\n", (subvolDim)*(subvolDim)*(subvolDim)*sizeof(float));
    freeHostMem();
    return false;
  } 

  // Allocate cuda array and global mem for subvol. Bind texture to cuda array
  initCuda();

  debugPrint("Total volume: %zux %zux %zu\n", volSize.width, volSize.height, volSize.depth);
  debugPrint("Total input volume: %ux%ux%u\n", datainfo.n_input[0], datainfo.n_input[1], datainfo.n_input[2]);
  debugPrint("Total subvolumes: %zux%zux%zu\n", subvolIndicesExtents.width, subvolIndicesExtents.height, subvolIndicesExtents.depth);

  //scale the input volume
  normalizeVolume();

  initInterfaceMulti();// Initialize level set interface for the phi function
  computeSDTEikonal();

  //writeRawData(h_subvol, "raw.out");
  //writeSlices(h_phi, "res", 20);

  if(nIter > 0)
  {
    debugPrint("Solving the PDE...\n");
    for (int i = 0; i<nIter; i++)
      if(!solve()) return false;
   }

  //writeRawData(h_phi, "raw.out");
  //float offset = (float)(datainfo.n_input[0] + datainfo.n_input[1] + datainfo.n_input[2]);
  //writevtk("phi.vtk", &datainfo, h_phi, offset); // pads the input volume to match the integral subvolume partitions in the volume

  memcpy(h_vol, h_buffer[currentPhi], volSize.width*volSize.height*volSize.depth*sizeof(float));
  cleanUp();
  return true;
}

void MSLevelSet::cleanUp()
{
  freeHostMem();
  CUDA_SAFE_CALL(cudaFreeArray(d_phiArray));
  CUDA_SAFE_CALL(cudaFreeArray(d_volArray));
  CUDA_SAFE_CALL(cudaUnbindTexture(texPhi)); //unbind texture
  CUDA_SAFE_CALL(cudaUnbindTexture(texVol)); //unbind texture
  CUDA_SAFE_CALL(cudaFree(d_volPPtr.ptr));
  CUDA_SAFE_CALL(cudaFree(d_spare1PPtr.ptr));
  CUDA_SAFE_CALL(cudaFree(d_spare2PPtr.ptr));
  CUDA_SAFE_CALL(cudaFree(d_spare3PPtr.ptr));

  CUT_SAFE_CALL(cutDeleteTimer(timer));
}

void MSLevelSet::freeHostMem()
{
  //if(h_vol) free(h_vol); //<--needed by volRover for display, isosurf extraction, etc.
  if(h_buffer[0]) {free(h_buffer[0]); h_buffer[0] = NULL;}
  if(h_buffer[1]) {free(h_buffer[1]); h_buffer[1] = NULL;}
  if(h_subvolSpare) { free(h_subvolSpare); h_subvolSpare = NULL;}
}

void MSLevelSet::setParameters(int _width, int _height, int _depth, MSLevelSetParams *MSLSParams)
{
  datainfo.n[0] = datainfo.n_input[0] = _width;
  datainfo.n[1] = datainfo.n_input[1] = _height;
  datainfo.n[2] = datainfo.n_input[2] = _depth;
  
  //extract MS LevelSet parameters
  lambda1 = MSLSParams->lambda1; 
  lambda2 = MSLSParams->lambda2 ;
  mu = MSLSParams->mu; 
  nu = MSLSParams->nu;
  delta_t = MSLSParams->deltaT; 
  epsilon = MSLSParams->epsilon;
  nIter = MSLSParams->nIter; 
  DTWidth = MSLSParams->DTWidth; 
  subvolDim = MSLSParams->subvolDim;
  PDEBlockDim = MSLSParams->PDEBlockDim; 
  avgBlockDim = MSLSParams->avgBlockDim;
  SDTBlockDim = MSLSParams->SDTBlockDim; 
  superEllipsoidPower = MSLSParams->superEllipsoidPower;
  BBoxOffset = MSLSParams->BBoxOffset;
  init_interface_method = MSLSParams->init_interface_method;
  volval_min = MSLSParams->volval_min;
  volval_max = MSLSParams->volval_max;
  
  //pretty print all the parameters to console!
  debugPrint("---Level set parameters---\n");
  debugPrint("Lambda1: %f\n", lambda1);
  debugPrint("Lambda2: %f\n", lambda2);
  debugPrint("Mu: %f\n", mu);
  debugPrint("Nu: %f\n", nu);
  debugPrint("Delta T: %f\n", delta_t);
  debugPrint("Epsilon: %f\n", epsilon);
  debugPrint("Max solver iterations: %d\n", nIter);
  debugPrint("DT band width: %d\n", DTWidth);
  debugPrint("BBox offset for initialization: %d\n", BBoxOffset);
  debugPrint("SuperEllipsoidal power: %f\n", superEllipsoidPower);
  debugPrint("Use BBox or Super-Ellipsoid for interface initialization?: %s\n", (init_interface_method==BBOX)?"BBox":(init_interface_method==SUPER_ELLIPSOID)?"Super-Ellipsoid":"Unspecified!");
  debugPrint("---CUDA parameters---\n");
  debugPrint("Sub volume dimension: %d^3\n", subvolDim);
  debugPrint("SDT kernel block dimension: %d\n", SDTBlockDim);
  debugPrint("Average kernel block dimension: %d\n", avgBlockDim);
  debugPrint("PDE kernel block dimension: %d\n\n", PDEBlockDim);
}

bool MSLevelSet::runSolver(float *vol, int _width, int _height, int _depth, MSLevelSetParams *MSLSParams,
               		   void (*evolutionCallback)(const float *vol, int dimx, int dimy, int dimz))
{
  h_vol = vol;
  if(!h_vol) return false;
  
  setParameters (_width, _height, _depth, MSLSParams);

  int result = solverMain();
  return result;
}

bool MSLevelSet::normalizeVolume()//To lie in the range [0 255]
{
  debugPrint("Normalizing input volume from [%f, %f] to [0 255]... ", volval_min, volval_max);
  unsigned long nelem = volSize.width*volSize.height*volSize.depth;
  float l_factor = 255.0/(volval_max - volval_min);
  double t0, t1;
  t0 = getTime();
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
  for(unsigned long i=0; i< nelem; i++) {
    h_vol[i] = (h_vol[i]  - volval_min)*l_factor;
  }
  t1 = getTime();
  milliseconds = t1 - t0;
  debugPrint(" %f ms\n", milliseconds);
  return true;
}
