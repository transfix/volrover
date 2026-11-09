//============================================================================
/* Copyright (c) The Technical University of Denmark
 * Developed at the Computational Visualization Center (CVC), The University of Texas at Austin
 * Author: Ojaswa Sharma
 * E-mail: os@imm.dtu.dk
 * File: higherorderseg_global.cu
 * CUDA global texture declarations.
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
#ifndef _HIGHER_ORDER_SEG_GLOBAL_H_
#define _HIGHER_ORDER_SEG_GLOBAL_H_

#include <cutil_math.h>

texture<float, 3, cudaReadModeElementType> texPhi; //3D texture for phi
texture<float, 3, cudaReadModeElementType> texVol; //3D texture for intensity
texture<float, 3, cudaReadModeElementType> texCoeff; //3D texture for cubic coefficient computation
 
#endif
