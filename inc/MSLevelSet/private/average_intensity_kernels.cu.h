//============================================================================
/* Copyright (c) The Technical University of Denmark
 * Author: Ojaswa Sharma
 * E-mail: os@imm.dtu.dk
 * File: average_intensity_kerenels.cu
 * CUDA kernels to compute average intensities inside and outside the level set interface.
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
#ifndef _AVERAGE_INTENSITY_KERNELS_
#define _AVERAGE_INTENSITY_KERNELS_

//Accumulate values of the Heaviside function (H_in), Intensity * H_in,  and Intensity * (1 - H_in) in accum_H_in, accum_IH_in,  and accum_IH_out.
__global__ void accumulate_average_intensities(cudaPitchedPtr accum_H_in, cudaPitchedPtr accum_IH_in, cudaPitchedPtr accum_IH_out, cudaExtent logicalGridSize, float epsilon, float nvoxels, int subvolDim, int avgBlockDim)
{
  int __x, __y, __z;
  unsigned int pitchz;
  pitchz = logicalGridSize.width*logicalGridSize.height;
  __z = blockIdx.x/pitchz;
  __y = (blockIdx.x - pitchz*__z)/logicalGridSize.width;
  __x = blockIdx.x - pitchz*__z - logicalGridSize.width*__y;

  //compute coordinates local (within subvolume)
  __x = __umul24(avgBlockDim, __x) + threadIdx.x;
  __y = __umul24(avgBlockDim, __y) + threadIdx.y;
  __z = __umul24(avgBlockDim, __z) + threadIdx.z;

  if(__x < 0 || __x > (subvolDim-1) || __y <0 || __y > (subvolDim-1) || __z < 0 || __z > (subvolDim-1)) return;
  
  volatile float x = (float)__x; 
  volatile float y = (float)__y; 
  volatile float z = (float)__z; 

  float phi_val = tex3D(texPhi, x, y, z);
  float H_in_val = 0.5*(1.0 + (2.0/PI)*atanf(phi_val/epsilon));  
  float IH_in_val = H_in_val*tex3D(texVol, x, y, z);
  float IH_out_val = (1.0 - H_in_val)*tex3D(texVol, x, y, z);

  float *row = (float*)((char*)accum_H_in.ptr + (__z*accum_H_in.ysize + __y)*accum_H_in.pitch);//output voxel: row[__x]
  row[__x] += H_in_val;
  
  row = (float*)((char*)accum_IH_in.ptr + (__z*accum_IH_in.ysize + __y)*accum_IH_in.pitch);//output voxel: row[__x]
  row[__x] += IH_in_val;

  row = (float*)((char*)accum_IH_out.ptr + (__z*accum_IH_out.ysize + __y)*accum_IH_out.pitch);//output voxel: row[__x]
  row[__x] += IH_out_val;
}

#endif
