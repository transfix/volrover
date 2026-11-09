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
__global__ void MultiPhase_PDE_update(cudaPitchedPtr phi_out, unsigned int nImplicits, int3 texPhi_partition, unsigned int i, cudaExtent logicalGridSize, float mu_factor,  float epsilon, float delta_t, int subvolDim, int PDEBlockDim)
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
  
  /// Compute internal energy (curvature energy)
  //1. compute gradient.
  volatile float grad_left[3];
  volatile float grad_right[3];
  volatile float x = (float)__x + (i%texPhi_partition.x)*subvolDim;
  volatile float y = (float)__y + ((i/texPhi_partition.x)%texPhi_partition.y)*subvolDim;
  volatile float z = (float)__z + (i/(texPhi_partition.x*texPhi_partition.y))*subvolDim; 
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
  float l_intEnergy;
  l_intEnergy = (grad_right[0] - grad_left[0] + grad_right[1] - grad_left[1] + grad_right[2] - grad_left[2])/2.0;
  l_intEnergy *= mu_factor;

  //3. Compute  external energy
  volatile float l_vol, l_H;
  volatile float l_extEnergy, l_Hprod;
  volatile float l_p1, l_p0;
  volatile unsigned int l_terms = 1<<(nImplicits-1);  
  ///->value substitution: 32//volatile unsigned int l_bsz = sizeof(nImplicits) <<4;
  volatile unsigned int l_index0, l_index1;
  volatile unsigned int l_bit;
  volatile float l_delta_dirac;

  l_extEnergy = 0.0;
  l_vol = tex3D(texVol, x, y, z);
  for(unsigned int k=0; k<l_terms; k++) {
    l_index0 = ((k>>i)<<(i+1)) | ((i==0)?0:((k<<(32-i))>>(32-i)));//insert 0 at ith place
    l_index1 = ((k>>i)<<(i+1)) | (0x1<<i) | ((i==0)?0:((k<<(32-i))>>(32-i)));//insert 1 at ith place
    l_p0 = l_vol - tex1Dfetch(texAverageValues, l_index0);
    l_p1 = l_vol - tex1Dfetch(texAverageValues, l_index1);	    
    l_p0 *= l_p0;
    l_p1 *= l_p1;
    l_Hprod = 1.0;
    for (unsigned int j=0; j<nImplicits; j++)//bit position
    {
      if(j==i) continue; //There is no multiplier for H_i
      l_bit = (l_index0>>j)&0x1;//Get the jth bit
      x = (float)__x + (j%texPhi_partition.x)*subvolDim;
      y = (float)__y + ((j/texPhi_partition.x)%texPhi_partition.y)*subvolDim;
      z = (float)__z + (j/(texPhi_partition.x*texPhi_partition.y))*subvolDim; 
      l_H = 0.5*(1.0 + (2.0/M_PI)*atanf(tex3D(texPhi, x, y, z)/epsilon)); 
      l_Hprod *= (float)l_bit + (1.0 - 2.0*float(l_bit))*l_H;
    }
    l_extEnergy += (l_p1 - l_p0)*l_Hprod;
  }

  ///PDEupdate
  //1. Update phi
  l_delta_dirac = (epsilon/M_PI)/(epsilon*epsilon + l_phi*l_phi);  
  float *row = (float*)((char*)phi_out.ptr + (__z*phi_out.ysize + __y)*phi_out.pitch);//output voxel: row[__x]
  row[__x] = l_phi + delta_t*l_delta_dirac*(l_intEnergy + l_extEnergy);
}

#endif
