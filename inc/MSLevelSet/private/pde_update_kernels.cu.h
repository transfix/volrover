//============================================================================
/* Copyright (c) The Technical University of Denmark
 * Author: Ojaswa Sharma
 * E-mail: os@imm.dtu.dk
 * File: pde_update_kerenels.cu
 * CUDA kernels for PDE update.
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
#ifndef _PDE_UPDATE_KERNELS_
#define _PDE_UPDATE_KERNELS_

// This kernel processes a (subvolDim-4)^3 volume instead of (subvolDim-2)^3 since the curvature computation requires that two voxels are present at either end of the border.
__global__ void PDE_update(cudaPitchedPtr phi_out, cudaExtent logicalGridSize, float mu_factor, float nu, float2 lambda, float2 c, float epsilon, float delta_t, int subvolDim, int PDEBlockDim)
{
  unsigned int __x, __y, __z;
  unsigned int pitchz;
  pitchz = logicalGridSize.width*logicalGridSize.height;
  __z = blockIdx.x/pitchz;
  __y = (blockIdx.x - pitchz*__z)/logicalGridSize.width;
  __x = blockIdx.x - pitchz*__z - logicalGridSize.width*__y;

  //compute coordinates local (within subvolume)
  __x = __umul24(PDEBlockDim, __x) + threadIdx.x;
  __y = __umul24(PDEBlockDim, __y) + threadIdx.y;
  __z = __umul24(PDEBlockDim, __z) + threadIdx.z;

  if(__x < 2 || __x > (subvolDim-3) || __y < 2 || __y > (subvolDim-3) || __z < 2 || __z > (subvolDim-3)) return;
  
  volatile float x = (float)__x; 
  volatile float y = (float)__y; 
  volatile float z = (float)__z; 

  /// Compute internal energy (curvature energy)
  //1. compute gradient.
  volatile float grad_left[3];
  volatile float grad_right[3];
  volatile float l_phi = tex3D(texPhi, x, y, z);
  float grad_norm;  
  grad_left[0] = (l_phi - tex3D(texPhi, x-2, y, z))/2.0;
  grad_left[1] = (l_phi - tex3D(texPhi, x, y-2, z))/2.0;
  grad_left[2] = (l_phi - tex3D(texPhi, x, y, z-2))/2.0;
  grad_right[0] = (tex3D(texPhi, x+2, y, z) - l_phi)/2.0;
  grad_right[1] = (tex3D(texPhi, x, y+2, z) - l_phi)/2.0;
  grad_right[2] = (tex3D(texPhi, x, y, z+2) - l_phi)/2.0;
  grad_norm = sqrt(grad_left[0]*grad_left[0] + grad_left[1]*grad_left[1] + grad_left[2]*grad_left[2] + TINY);
  grad_left[0] /= grad_norm; grad_left[1] /= grad_norm; grad_left[2] /= grad_norm;
  grad_norm = sqrt(grad_right[0]*grad_right[0] + grad_right[1]*grad_right[1] + grad_right[2]*grad_right[2] + TINY);
  grad_right[0] /= grad_norm; grad_right[1] /= grad_norm; grad_right[2] /= grad_norm;
  
  //2. compute curvature energy
  float curvature_energy;
  curvature_energy = (grad_right[0] - grad_left[0] + grad_right[1] - grad_left[1] + grad_right[2] - grad_left[2])/2.0;
  curvature_energy *= mu_factor;

  /// Compute extrenal energy (using Mumford-Shah functional)
  float external_energy;
  float l_I = tex3D(texVol, x, y, z);
  external_energy = -nu - lambda.x*(l_I - c.x)*(l_I - c.x) + lambda.y*(l_I - c.y)*(l_I - c.y);

  ///PDEupdate
  //1. Update phi
  float delta_dirac = (epsilon/M_PI)/(epsilon*epsilon + l_phi*l_phi);  
  float *row = (float*)((char*)phi_out.ptr + (__z*phi_out.ysize + __y)*phi_out.pitch);//output voxel: row[__x]
  row[__x] = l_phi + delta_t*delta_dirac*(curvature_energy + external_energy);
}

#endif
