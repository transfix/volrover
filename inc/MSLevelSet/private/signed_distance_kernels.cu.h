//============================================================================
/* Copyright (c) The Technical University of Denmark
 * Author: Ojaswa Sharma
 * E-mail: os@imm.dtu.dk
 * File: signed_distance_kerenels.cu
 * CUDA kernels for signed distance computation.
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
#ifndef _SIGNED_DISTANCE_KERNELS_
#define _SIGNED_DISTANCE_KERNELS_

#include <cutil_math.h>
typedef unsigned int uint;

#define MIN(x,y) (x) < (y) ? (x) : (y)
#define MAX(x,y) (x) > (y) ? (x) : (y)
#define EPS 0.01

__global__ void d_signed_distance_eikonal(cudaPitchedPtr d_volPPtr, cudaExtent logicalGridSize, int subvolDim, int SDTBlockDim)
{
  //An upwinding scheme to solve the Eikonal equation for redistancing.
  unsigned int __x, __y, __z;
  unsigned int pitchz;
  pitchz = logicalGridSize.width*logicalGridSize.height;
  __z = blockIdx.x/pitchz;
  __y = (blockIdx.x - pitchz*__z)/logicalGridSize.width;
  __x = blockIdx.x - pitchz*__z - logicalGridSize.width*__y;

  //compute coordinates local (within subvolume)
  __x = __umul24(SDTBlockDim, __x) + threadIdx.x;
  __y = __umul24(SDTBlockDim, __y) + threadIdx.y;
  __z = __umul24(SDTBlockDim, __z) + threadIdx.z;

  if(__x < 1 || __x > (subvolDim-2) || __y < 1 || __y > (subvolDim-2) || __z < 1 || __z > (subvolDim-2)) return;
  volatile float x = (float)__x; 
  volatile float y = (float)__y; 
  volatile float z = (float)__z; 
  
  //float reinit_delta_t = 0.2;
  float3 Dphi_p, Dphi_m;
  float a, b, c, d, e, f;
  float _sgn, _phi;
  float G;

  _phi = tex3D(texPhi, x, y, z);
  Dphi_m.x = _phi - tex3D(texPhi, x - 1., y, z);
  Dphi_p.x = tex3D(texPhi, x + 1., y, z) - _phi;
  Dphi_m.y = _phi - tex3D(texPhi, x, y - 1., z);
  Dphi_p.y = tex3D(texPhi, x, y + 1., z) - _phi;
  Dphi_m.z = _phi - tex3D(texPhi, x, y, z - 1.);
  Dphi_p.z = tex3D(texPhi, x, y, z + 1.) - _phi;
  a = (tex3D(texPhi, x + 1., y, z) - tex3D(texPhi, x - 1., y, z))/2.;
  b = (tex3D(texPhi, x, y + 1., z) - tex3D(texPhi, x, y - 1., z))/2.;
  c = (tex3D(texPhi, x, y, z + 1.) - tex3D(texPhi, x, y, z - 1.))/2.;
  G = a*a + b*b + c*c;
  _sgn = _phi/sqrtf(_phi*_phi + G); // Note G is squared magnitude of gradient!
  if(_sgn > 0.) {
    a = fmaxf(Dphi_m.x, 0.); b = fminf(Dphi_p.x, 0.);
    c = fmaxf(Dphi_m.y, 0.); d = fminf(Dphi_p.y, 0.);
    e = fmaxf(Dphi_m.z, 0.); f = fminf(Dphi_p.z, 0.);
    G = sqrtf(fmaxf(a*a, b*b) + fmaxf(c*c, d*d) + fmaxf(e*e, f*f)) - 1.0;
  } else if(_sgn < 0.) {
    a = fminf(Dphi_m.x, 0.); b = fmaxf(Dphi_p.x, 0.);
    c = fminf(Dphi_m.y, 0.); d = fmaxf(Dphi_p.y, 0.);
    e = fminf(Dphi_m.z, 0.); f = fmaxf(Dphi_p.z, 0.);
    G = sqrtf(fmaxf(a*a, b*b) + fmaxf(c*c, d*d) + fmaxf(e*e, f*f)) - 1.0;
  }
  else
    G = 0.;

  float *row = (float*)((char*)d_volPPtr.ptr + (__z*d_volPPtr.ysize + __y)*d_volPPtr.pitch);//output voxel: row[__x]
  row[__x] = _phi - 0.1*_sgn*G;// 0.1 is delta T for eikonal eqn solution
  return;
}  

  volatile float texval;
  volatile float texval2;
#endif
