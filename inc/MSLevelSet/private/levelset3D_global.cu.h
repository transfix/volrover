//============================================================================
/* Copyright (c) The Technical University of Denmark
 * Author: Ojaswa Sharma
 * E-mail: os@imm.dtu.dk
 * File: levelset3D_global.cu
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
#ifndef _LEVELSET3D_GLOBAL_H_
#define _LEVELSET3D_GLOBAL_H_

#include <cutil_math.h>

texture<float, 3, cudaReadModeElementType> texPhi; //3D texture for phi
texture<float, 3, cudaReadModeElementType> texVol; //3D texture for intensity

#endif
