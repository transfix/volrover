/*
  Copyright 2011 The University of Texas at Austin

	Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of MolSurf.

  MolSurf is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.


  MolSurf is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
#ifndef COMPUTE_H
#define COMPUTE_H

#include <Utility/utility.h>
#include <Contour/basic.h>

void sortVerts(float* &x0, float* &x1, float* &x2, float* &x3, float &v0, float& v1, float& v2, float& v3);
void triSurfIntegral(double* x1, double* x2, double* x3, float v1, float v2, float v3, float* fx, float* val, int nbucket, float min, float max, float scaling);
void triVolIntegral(double* x1, double* x2, double* x3, float v1, float v2, float v3, float* fx, float* val, float* cum, u_int nbucket, float min, float max, float scaling);
void tetSurfIntegral(float* x1, float* x2, float* x3, float* x4, float v1, float v2, float v3, float v4, float* fx, float* val, int nbucket, float min, float max, float scaling);
void tetVolIntegral(float* x1, float* x2, float* x3, float* x4, float v1, float v2, float v3, float v4, float* fx, float* val, float* cum, u_int nbucket, float min, float max, float scaling);

void intVolIntegral(float** p, float* u, float* v, float fx1[], float** val1,float fx2[], float** val2, int nbucket, float min1, float max1, float min2, float max2, float scaling);

#endif
