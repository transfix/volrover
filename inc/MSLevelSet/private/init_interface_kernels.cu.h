//============================================================================
/* Copyright (c) The Technical University of Denmark
 * Author: Ojaswa Sharma
 * E-mail: os@imm.dtu.dk
 * File: signed_distance_kerenels.cu
 * CUDA kernels for smooth interface initialization
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
#ifndef _INIT_INTERFACE_KERNELS_
#define _INIT_INTERFACE_KERNELS_

#include <cutil_math.h>
typedef unsigned int uint;

#define MIN(x,y) (x) < (y) ? (x) : (y)
#define MAX(x,y) (x) > (y) ? (x) : (y)
#define EPS 0.01

__global__ void d_init_interface_multi_smooth(cudaPitchedPtr d_volPPtr, cudaExtent logicalGridSize, cudaExtent subvolIndex, float max_dist, int subvolDim, int SDTBlockDim, float3 offset, float box, float radius_reach, float2 radii, float3 upper_limit)
{
  unsigned int _w, _h, _d;
  unsigned int pitchz = logicalGridSize.width*logicalGridSize.height;

  _d = blockIdx.x/pitchz;
  _h = (blockIdx.x - (_d*pitchz))/logicalGridSize.width;
  _w = blockIdx.x - (_d*pitchz) - (_h*logicalGridSize.width);

  //compute coordinates local (wthin subvolume) and global (for the large volume) coordinates
  unsigned int l_x = SDTBlockDim*_w + threadIdx.x;
  unsigned int l_y = SDTBlockDim*_h + threadIdx.y;
  unsigned int l_z = SDTBlockDim*_d + threadIdx.z;

  if(l_x < 1 || l_x > (subvolDim-2) || l_y < 1 || l_y > (subvolDim-2) || l_z < 1 || l_z > (subvolDim-2)) return;

  // Reuse vars _w, _h, _d to store global coordinates of the voxel
  _w = l_x + __umul24(subvolDim-2, subvolIndex.width);
  _h = l_y + __umul24(subvolDim-2, subvolIndex.height);
  _d = l_z + __umul24(subvolDim-2, subvolIndex.depth);

  float3 coord, cen;
  float *row = (float*)((char*)d_volPPtr.ptr + (l_z*d_volPPtr.ysize + l_y)*d_volPPtr.pitch);
  coord.x = float(_w) - offset.x/2.;
  coord.y = float(_h) - offset.y/2.;
  coord.z = float(_d) - offset.z/2.;
  if(coord.x < 0. || coord.x > upper_limit.x || coord.y < 0. || coord.y > upper_limit.y || coord.z < 0. || coord.z > upper_limit.z) {
    row[l_x] = max_dist;
    return;
  }
  cen.x = floorf(coord.x/box)*box + radius_reach;
  cen.y = floorf(coord.y/box)*box + radius_reach;
  cen.z = floorf(coord.z/box)*box + radius_reach;
  float dist = sqrtf((coord.x - cen.x)*(coord.x - cen.x) + (coord.y - cen.y)*(coord.y - cen.y) + (coord.z - cen.z)*(coord.z - cen.z));
  row[l_x] = (dist <= radii.y)?(radii.x - dist):max_dist;  
}

#endif
