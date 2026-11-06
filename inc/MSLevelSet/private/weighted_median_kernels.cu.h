//============================================================================
/* Copyright (c) The Technical University of Denmark
 * Author: Ojaswa Sharma
 * E-mail: os@imm.dtu.dk
 * File: weighted_median_kerenels.cu
 * CUDA kernels for weighted median computation.
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
#ifndef _WEIGHTED_MEDIAN_KERNELS_
#define _WEIGHTED_MEDIAN_KERNELS_

//#include <cutil_math.h>

__device__ float d_medianGuess;
__global__ void d_weighted_median(cudaExtent logicalGridSize, cudaPitchedPtr accum_lowerWeightSum, cudaPitchedPtr accum_upperWeightSum, float epsilon, int subvolDim, int medianBlockDim)
{
  int __x, __y, __z;
  unsigned int pitchz;
  pitchz = logicalGridSize.width*logicalGridSize.height;
  __z = blockIdx.x/pitchz;
  __y = (blockIdx.x - pitchz*__z)/logicalGridSize.width;
  __x = blockIdx.x - pitchz*__z - logicalGridSize.width*__y;

  //compute coordinates local (within subvolume)
  __x = __umul24(medianBlockDim, __x) + threadIdx.x;
  __y = __umul24(medianBlockDim, __y) + threadIdx.y;
  __z = __umul24(medianBlockDim, __z) + threadIdx.z;

  if(__x < 0 || __x > (subvolDim-1) || __y < 0 || __y > (subvolDim-1) || __z < 0 || __z > (subvolDim-1)) return;
  
  volatile float x = (float)__x; 
  volatile float y = (float)__y; 
  volatile float z = (float)__z; 

  volatile float weight = tex3D(texPhi, x, y, z);
  //weight = 0.5*(1.0 + (2.0/M_PI)*atanf(weight/epsilon));
  weight = 0.5 + 0.318309886*atanf(weight/epsilon);

  volatile float val = tex3D(texVol, x, y, z);
  float *row;
  if(val < d_medianGuess)
    row = (float*)((char*)accum_lowerWeightSum.ptr + (__z*accum_lowerWeightSum.ysize + __y)*accum_lowerWeightSum.pitch);//output voxel: row[__x]
  else
    row = (float*)((char*)accum_upperWeightSum.ptr + (__z*accum_upperWeightSum.ysize + __y)*accum_upperWeightSum.pitch);//output voxel: row[__x]

  row[__x] += weight; 
}
#endif
