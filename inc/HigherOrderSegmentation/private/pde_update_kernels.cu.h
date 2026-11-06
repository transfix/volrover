//============================================================================
/* Copyright (c) The Technical University of Denmark
 * Developed at the Computational Visualization Center (CVC), The University of Texas at Austin
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

#define MIN(x,y) (x) < (y) ? (x) : (y)
#define MAX(x,y) (x) > (y) ? (x) : (y)

//Constant array to store derivatives (fx, fy, fz, fxx, fyy, fzz, fxy, fyz, fzx) at 27 points. Compuation is done on host side.
__device__ __constant__ float d_cubicDerivatives[27*9];

//Compute divergence of normalized gradient of the cubic spline function phi
//Since the function is evalated only at grid point, we need just the three neighbors in each x, y and z
// N.B.: Make sure that texCoeff texture has an extra 1 voxel border around. If near boundary of volume, add an extra border
// surrounding the volume with zeroes. This is needed to avoid 'if conditions' here. Further, at the boundaries, there should be no
// coefficient c_{-1} and c_{n} (1-D example). This is taken care of if a 1-voxel cover is layered around the volume and texture.
__device__ float divergenceOfNormalizedGradient(float x, float y, float z)
{
  float phix=0.0, phiy=0.0, phiz=0.0, phixx=0.0, phiyy=0.0, phizz=0.0, phixy=0.0, phiyz=0.0, phizx=0.0;
  float coeff = 0.0;
  for(int k=0; k<3; k++)
    for(int j=0; j<3; j++)
      for(int i=0; i<3; i++)
      {
        coeff = tex3D(texCoeff, x+(float)i -1.0, y+(float)j - 1.0, z+(float)k - 1.0);
        phix += coeff*d_cubicDerivatives[(i+3*j+9*k)*9];
        phiy += coeff*d_cubicDerivatives[(i+3*j+9*k)*9 + 1];
        phiz += coeff*d_cubicDerivatives[(i+3*j+9*k)*9 + 2];
        phixx += coeff*d_cubicDerivatives[(i+3*j+9*k)*9 + 3];
        phiyy += coeff*d_cubicDerivatives[(i+3*j+9*k)*9 + 4];
        phizz += coeff*d_cubicDerivatives[(i+3*j+9*k)*9 + 5];
        phixy += coeff*d_cubicDerivatives[(i+3*j+9*k)*9 + 6];
        phiyz += coeff*d_cubicDerivatives[(i+3*j+9*k)*9 + 7];
        phizx += coeff*d_cubicDerivatives[(i+3*j+9*k)*9 + 8];
      }
  
  //Reuse variable coeff to store ||grad(phi)||^2
  coeff = phix*phix + phiy*phiy + phiz*phiz + TINY;
  //coeff *=sqrtf(coeff);
  //float divNormGrad = phix*phix*(phiyy + phizz) + phiy*phiy*(phixx + phizz) + phiz*phiz*(phixx + phiyy) - 2.0*(phix*phiy*phixy + phiy*phiz*phiyz + phiz*phix*phizx);
  //divNormGrad /=coeff;
  
  float divNormGrad = phix*phix*phixx + phiy*phiy*phiyy + phiz*phiz*phizz;
  divNormGrad += 2.0*(phix*phiy*phixy + phiy*phiz*phiyz + phiz*phix*phizx);
  divNormGrad /= coeff;
  divNormGrad = (phixx + phiyy + phizz) - divNormGrad;
  divNormGrad /=sqrt(coeff);
  
  return divNormGrad;
}

__global__ void HigherOrder_PDE_update(cudaPitchedPtr phi_out, cudaExtent logicalGridSize, float mu_factor, float nu, float2 lambda, float2 c, float epsilon, float delta_t, int subvolDim, int PDEBlockDim)
{
  //Currently supports Mumford-Shah functional. More to be included :)
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

  if(__x < 1 || __x > (subvolDim-2) || __y < 1 || __y > (subvolDim-2) || __z < 1 || __z > (subvolDim-2)) return;
  
  volatile float x = (float)__x; 
  volatile float y = (float)__y; 
  volatile float z = (float)__z; 
  
  ///Compute higher order term (internal energy term) using Cubic splines
  float internal_energy = mu_factor*divergenceOfNormalizedGradient(x, y, z);
  
  ///Compute extrenal energy
  float external_energy;
  float volval = tex3D(texVol, x, y, z);
  external_energy = -nu - lambda.x*(volval - c.x)*(volval - c.x) + lambda.y*(volval - c.y)*(volval - c.y);

  ///PDEupdate
  //Update phi
  float phival = tex3D(texPhi, x, y, z);
  float delta_dirac = (epsilon/M_PI)/(epsilon*epsilon + phival*phival);  
  
  float *row = (float*)((char*)phi_out.ptr + (__z*phi_out.ysize + __y)*phi_out.pitch);//output voxel: row[__x]
  row[__x] = phival + delta_t*delta_dirac*(internal_energy + external_energy);
}
#endif
