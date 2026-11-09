//============================================================================
/* Copyright (c) The Technical University of Denmark
 * Developed at the Computational Visualization Center (CVC), The University of Texas at Austin
 * Author: Ojaswa Sharma
 * E-mail: os@imm.dtu.dk
 * File: cubic_coefficients_kernels.cu
 * CUDA kernels for cubic coefficients computation for phi.
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
#ifndef _CUBIC_COEFFICIENTS_KERNELS_
#define _CUBIC_COEFFICIENTS_KERNELS_

__global__ void d_cubic_coefficients_1DZ(cudaPitchedPtr coeff_out, cudaExtent volSize)
{
  //compute image x and y coordinates for the current thread.
  unsigned int _x, _y;
  unsigned int rindex = blockIdx.x*blockDim.x + threadIdx.x; //running index
  _x = rindex%volSize.width;
  _y = rindex/volSize.width;
  if(_y >= volSize.height) return; //Out of bounds of the image.
  
  int _k;
  float c0_plus = 0.0;
  float z1 = sqrtf(3.0) - 2.0;
  float z1_2n = powf(z1, 2.0*((float)volSize.depth + 1.0));
  float z1_k = 0.0;
  for(_k=1; _k<=volSize.depth; _k++)
  {
    z1_k = powf(z1, (float)_k);
    c0_plus += tex3D(texCoeff, _x, _y, _k - 1)*(z1_k - z1_2n/z1_k);
  }
  c0_plus = -c0_plus/(1.0 - z1_2n);
  
  float ck_plus, ck_minus = 0.0;
  for(_k=volSize.depth; _k>=1; _k--)
  {
    //compute c0_plus on the fly for every value. There is not enough memory to store it.
    ck_plus = c0_plus;
    for(rindex =1; rindex<=_k; rindex++) //reuse rindex
      ck_plus = tex3D(texCoeff, _x, _y, rindex - 1) + z1*ck_plus;
    ck_minus = z1*(ck_minus - ck_plus);  
    *((float*)((char*)coeff_out.ptr + (_y + (_k - 1)*coeff_out.ysize)*coeff_out.pitch) + _x) =  6.0*ck_minus;
  }
}

__global__ void d_cubic_coefficients_1DZ_Fast(cudaPitchedPtr coeff_out, cudaExtent volSize)
{
  //compute image x and y coordinates for the current thread.
  unsigned int _x, _y;
  unsigned int rindex = blockIdx.x*blockDim.x + threadIdx.x; //running index
  _x = rindex%volSize.width;
  _y = rindex/volSize.width;
  if(_y >= volSize.height) return; //Out of bounds of the image.
  
  int _k;
  float c0_plus = 0.0;
  float z1 = sqrtf(3.0) - 2.0;
  float z1_2n = powf(z1, 2.0*((float)volSize.depth + 1.0));
  for(_k=1; _k<=volSize.depth; _k++)
    c0_plus = (tex3D(texCoeff, _x, _y, volSize.depth - _k) +  c0_plus)*z1;
  c0_plus = -c0_plus/(1.0 - z1_2n);
  
  float ck_plus, ck_minus = 0.0;
  for(_k=volSize.depth; _k>=1; _k--)
  {
    //compute c0_plus on the fly for every value. There is not enough memory to store it.
    ck_plus = c0_plus;
    for(rindex =1; rindex<=_k; rindex++) //reuse rindex
      ck_plus = tex3D(texCoeff, _x, _y, rindex - 1) + z1*ck_plus;
    ck_minus = z1*(ck_minus - ck_plus); 
    *((float*)((char*)coeff_out.ptr + (_y + (_k - 1)*coeff_out.ysize)*coeff_out.pitch) + _x) = 6.0*ck_minus;
  }
}

__global__ void d_cubic_coefficients_1DY(cudaPitchedPtr coeff_out, cudaExtent volSize)
{
  //compute image x and y coordinates for the current thread.
  unsigned int _x, _z;
  unsigned int rindex = blockIdx.x*blockDim.x + threadIdx.x; //running index
  _x = rindex%volSize.width;
  _z = rindex/volSize.width;
  if(_z >= volSize.depth) return; //Out of bounds of the image.
  
  int _k;
  float c0_plus = 0.0;
  float z1 = sqrtf(3.0) - 2.0;
  float z1_2n = powf(z1, 2.0*((float)volSize.height + 1.0));
  float z1_k = 0.0;
  for(_k=1; _k<=volSize.height; _k++)
  {
    z1_k = powf(z1, (float)_k);
    c0_plus += tex3D(texCoeff, _x, _k - 1, _z)*(z1_k - z1_2n/z1_k);
  }
  c0_plus = -c0_plus/(1.0 - z1_2n);
  
  float ck_plus, ck_minus = 0.0;
  for(_k=volSize.height; _k>=1; _k--)
  {
    //compute c0_plus on the fly for every value. There is not enough memory to store it.
    ck_plus = c0_plus;
    for(rindex =1; rindex<=_k; rindex++) //reuse rindex
      ck_plus = tex3D(texCoeff, _x, rindex - 1, _z) + z1*ck_plus;
    ck_minus = z1*(ck_minus - ck_plus);  
    *((float*)((char*)coeff_out.ptr + ((_k - 1) + _z*coeff_out.ysize)*coeff_out.pitch) + _x) = 6.0*ck_minus; 
  }
}

__global__ void d_cubic_coefficients_1DY_Fast(cudaPitchedPtr coeff_out, cudaExtent volSize)
{
  //compute image x and y coordinates for the current thread.
  unsigned int _x, _z;
  unsigned int rindex = blockIdx.x*blockDim.x + threadIdx.x; //running index
  _x = rindex%volSize.width;
  _z = rindex/volSize.width;
  if(_z >= volSize.depth) return; //Out of bounds of the image.
  
  int _k;
  float c0_plus = 0.0;
  float z1 = sqrtf(3.0) - 2.0;
  float z1_2n = powf(z1, 2.0*((float)volSize.height + 1.0));
  for(_k=1; _k<=volSize.height; _k++)
    c0_plus = (tex3D(texCoeff, _x, volSize.height - _k, _z) +  c0_plus)*z1;  
  c0_plus = -c0_plus/(1.0 - z1_2n);
  
  float ck_plus, ck_minus = 0.0;
  for(_k=volSize.height; _k>=1; _k--)
  {
    //compute c0_plus on the fly for every value. There is not enough memory to store it.
    ck_plus = c0_plus;
    for(rindex =1; rindex<=_k; rindex++) //reuse rindex
      ck_plus = tex3D(texCoeff, _x, rindex - 1, _z) + z1*ck_plus;
    ck_minus = z1*(ck_minus - ck_plus);  
    *((float*)((char*)coeff_out.ptr + ((_k - 1) + _z*coeff_out.ysize)*coeff_out.pitch) + _x) = 6.0*ck_minus; 
  }
}

__global__ void d_cubic_coefficients_1DX(cudaPitchedPtr coeff_out, cudaExtent volSize)
{
  //compute image x and y coordinates for the current thread.
  unsigned int _y, _z;
  unsigned int rindex = blockIdx.x*blockDim.x + threadIdx.x; //running index
  _y = rindex%volSize.height;
  _z = rindex/volSize.height;
  if(_z >= volSize.depth) return; //Out of bounds of the image.
  
  int _k;
  float c0_plus = 0.0;
  float z1 = sqrtf(3.0) - 2.0;
  float z1_2n = powf(z1, 2.0*((float)volSize.width + 1.0));
  float z1_k = 0.0;
  for(_k=1; _k<=volSize.width; _k++)
  {
    z1_k = powf(z1, (float)_k);
    c0_plus += tex3D(texCoeff, _k - 1, _y, _z)*(z1_k - z1_2n/z1_k);
  }
  c0_plus = -c0_plus/(1.0 - z1_2n);
  
  float ck_plus, ck_minus = 0.0;
  float *coeffs = (float*)((char*)coeff_out.ptr + (_z*coeff_out.ysize + _y)*coeff_out.pitch);
  for(_k=volSize.width; _k>=1; _k--)
  {
    //compute c0_plus on the fly for every value. There is not enough memory to store it.
    ck_plus = c0_plus;
    for(rindex =1; rindex<=_k; rindex++) //reuse rindex
      ck_plus = tex3D(texCoeff, rindex - 1, _y, _z) + z1*ck_plus;
    ck_minus = z1*(ck_minus - ck_plus);  
    coeffs[_k - 1] = 6.0*ck_minus;  
  }
}

__global__ void d_cubic_coefficients_1DX_Fast(cudaPitchedPtr coeff_out, cudaExtent volSize)
{
  //compute image x and y coordinates for the current thread.
  unsigned int _y, _z;
  unsigned int rindex = blockIdx.x*blockDim.x + threadIdx.x; //running index
  _y = rindex%volSize.height;
  _z = rindex/volSize.height;
  if(_z >= volSize.depth) return; //Out of bounds of the image.
  
  int _k;
  float c0_plus = 0.0;
  float z1 = sqrtf(3.0) - 2.0;
  float z1_2n = powf(z1, 2.0*((float)volSize.width + 1.0));
  for(_k=1; _k<=volSize.width; _k++)
    c0_plus = (tex3D(texCoeff, volSize.width - _k, _y, _z) +  c0_plus)*z1;  
  c0_plus = -c0_plus/(1.0 - z1_2n);
  
  float ck_plus, ck_minus = 0.0;
  float *coeffs = (float*)((char*)coeff_out.ptr + (_z*coeff_out.ysize + _y)*coeff_out.pitch);
  for(_k=volSize.width; _k>=1; _k--)
  {
    //compute c0_plus on the fly for every value. There is not enough memory to store it.
    ck_plus = c0_plus;
    for(rindex =1; rindex<=_k; rindex++) //reuse rindex
      ck_plus = tex3D(texCoeff, rindex - 1, _y, _z) + z1*ck_plus;
    ck_minus = z1*(ck_minus - ck_plus);  
    coeffs[_k - 1] = 6.0*ck_minus;  
  }
}
#endif