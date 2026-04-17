/*
  Copyright 2006-2007 The University of Texas at Austin

        Authors: Dr. Xu Guo Liang <xuguo@ices.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of HLevelSet.

  HLevelSet is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  HLevelSet is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

// HLevelSet.cpp: implementation of the HLevelSet class.
//
//////////////////////////////////////////////////////////////////////

#include <HLevelSet/HLevelSet.h>
//#include "SimpleVolumeData.h"
//#include "../SimpleVolumeDataIsocontourer.h"
//#include "Atom.h"
//#include "BlurMapsDataManager.h"
#include <HLevelSet/Misc.h>
//#include "Geometry.h"

#include <iostream>
#include <vector>
#include <math.h>

#include <boost/tuple/tuple.hpp>

static float   InitialAntiCausalCoefficient(float *, int, float);
static float   InitialCausalCoefficient(float *, int, float, float);
static void    ConvertToInterpolationCoefficients_Qu(float *, int, float *, int ,float);

static void    ConvertToInterpolationCoefficients_1D(float *, int, float);
static void    ConvertToInterpolationCoefficients_2D(float *, int, int, float);
static void    ConvertToInterpolationCoefficients_3D(float *, int, int, int, float);

static void    Tensor_333(float *, float *, float *, float *);
static void    Tensor_444(float *, float *, float *, float *);
static void    Take_27_Coefficients(float *, int, int, int, int, int, int, float *);
static void    Take_64_Coefficients(float *, int, int, int, int, int, int, float *);
static void    Evaluat_Four_Basis(float x, float *values);

void ReSamplingCubicSpline1OrderPartials(float *coeff, int nx, int ny, int nz, float *dxyz, float *minExtent, float *zeroToTwoPartialValue,  float *p) ;


static void    EvaluateCubicSplineOrder1PartialsAtGridPoint(float *, float, float, float, int, int, int, int, int, int, float*);
static void    EvaluateCubicSplineOrder2PartialsAtGridPoint(float *, float, float, float, int, int, int, int, int, int, float*);
static void    Divided_DifferenceOrder2PartialsAtGridPoint(float *, float, float, float, int, int, int, int, int, int, float *);

static void    InitialData_Sphere(float *, float, float, float, float, float, float, int, int, int);
static void    InitialData_Box(float *, float, float, float, float, float, float, int, int, int);
static void    InitialData_IrBox(float *, float, float, float, float, float, float, int, int, int);
static void    InitialData_Scatter(float *v, float *, float * init, int vnum, float minx, float maxx, 
                  float miny, float maxy, float minz, float maxz, int nx, int ny, int nz, float *i, float);
static void    Fast_Distance_Function(float *, float * init, int vnum, float minx, float maxx, 
                  float miny, float maxy, float minz, float maxz, int nx, int ny, int nz, float);
static void    Fast_Sweeping(float *, int nx, int ny, int nz, float dx);


static void   Constraint_MeanCurvatureFlow(float *func_h, float *fun, float *fun_bak, float *coeff, 
                          float weight, float dt, 
                          float dx, float dy, float dz, int nx, int ny, int nz, int t);
static void    MeanCurvatureFlow(float *, float *, float *, float, float, float, float, int, int, int, int);
static void    ComputeTensorXYZ();

static void    ReInitilazation(float *, float *, float, float, float, int, int, int);
static void    ReInitilazation_Upwind_Eno_Engquist(float *, float *, float *, float, float, float, int, int, int);
static void    ReInitilazation_Upwind_Eno_Godunov(float *, float *, float *, float, float, float, int, int, int, int);
static void    ReInitilazation_Upwind_Eno_Godunov_Xu(float *, float *, float *, float, float, float, int, int, int, int);
//void    Fast_Reintialization(float *, float *, float *, unsigned int* , float* , float* );


static void    DiviededDifferencing_3j(float *, float *, float *, float *, float);
static void    DiviededDifferencing_2j(float *, float *, float *, float);
static void    DiviededDifferencing_1j(float *, float *, float);

static void    Cubic_Eno_Interpolation(float *f0, float *u_plus, float *u_minus, float dx);
static void    Quadr_Eno_Interpolation(float *f0, float *u_plus, float *u_minus, float dx);
static void    Linear_Eno_Interpolation(float *f0, float *u_plus, float *u_minus, float dx);

static void    Get_Seven_Function_Values_X(float *coeff, int nx, int ny, int nz, int i, int j, int k, float *fx);
static void    Get_Seven_Function_Values_Y(float *coeff, int nx, int ny ,int nz, int i, int j, int k, float *fy);
static void    Get_Seven_Function_Values_Z(float *coeff, int nx, int ny, int nz, int i, int j, int k, float *fz);

static void    ReSamplingCubicSpline(float *coeff, int nx, int ny, int nz, float *funvalues, int Nx, int Ny, int Nz);
static void    Find_Box_Dim(int Max_dim, const std::vector<float> vertex_Positions,float* minExtent, float* maxExtent, unsigned int *dim);
static void    Find_Box_Dim(int Max_dim, const std::vector<float> vertex_Positions,float* minExtent, float* maxExtent, unsigned int *dim, float edgelength);
static void    Switch(unsigned int *dim, float *minExt, float *maxExt);
static void    One_Step_Sweeping(float *dSquare, int nx, int ny, int nz, float dx, int i, int j, int k);


static void    HLevel_set();

static float   TakeACoefficient_Fast(float *c, int nx, int ny, int nz, int u, int v, int w);
static float   TakeACoefficient_Slow(float *c, int nx, int ny, int nz, int u, int v, int w);
static float   EvaluateCubicSplineAtGridPoint(float *, int, int, int, int, int, int);
static float   EvaluateCubicSplineAtAnyGivenPoint(float *, int, int, int, float, float, float);
static float   Extreme_Positive(float a, float b);
static float   Extreme_Negative(float a, float b);
static float   Gradient_2(float *fun, float dx, float dy, float dz, int i, int j, int k, int nx, int ny, int nz);

static bool    Un_Treat_Boundary_Point(int i, int j, int k, float *signs,float *dis, float offset, int nx, int ny, int nz);
static void    ComputeDistanceFunction(float *distanceD,float *init,  int wid, int nnx, int nny, int nnz, 
                                float ddx, float ddy, float ddz, float xS, float yS, float zS);
static void ComputeDistanceFunction_To_Triangle(float *distanceD,  float *p1, float *p2, 
                             float * p3,  int wid, int nnx, int nny, int nnz, 
                             float ddx, float ddy, float ddz, 
                             float xS,  float yS,  float zS);


//float Max_Of_Two(float a, float b);
//float Min_Of_Two(float a, float b);
//float Max_Of_Three(float a,float b,float c);
//float Min_Of_Three(float a,float b,float c);
static float PointToLineDistance(float p[], float p1[], float p2[]);
static float PointToTriangleDistance(float p[], float p1[], float p2[], float p3[]);

// External variables
//extern std::vector<double> vertexPositions;
//extern int vNumber;

// arand: eliminating globals...
//extern float edgelength;
//extern int end , Max_dim;

// Local variables
static float  TensorF[27],   TensorFx[27],  TensorFy[27],  TensorFz[27],  TensorFxx[27],
              TensorFxy[27], TensorFxz[27], TensorFyy[27], TensorFyz[27], TensorFzz[27];

static float  Height0, Height;
static float  TotalTime, offset;



#define CVC_DBL_EPSILON 1.0e-9f
#define XUGUO_Z1  sqrt(3.0) - 2.0
#define OneSix3 1.0/216.0
#define OneSix2 2.0/108.0
#define OneSix1 4.0/54.0
#define OneSix0 8.0/27.0

// class Geometry;
using namespace HLevelSetNS;
using namespace PDBParser;

HLevelSet::HLevelSet()
{

}

HLevelSet::HLevelSet(float el, int e, int md)
{
    edgelength=el;
    end=e;
    Max_dim=md;
    
}



HLevelSet::~HLevelSet()
{

}

bool HLevelSet::computeFunction_Zhang(const std::vector<float> vertex_Positions, float* funcvalue, unsigned int* dim, float* minExt, float* maxExt, float edgelength, int end)
{
float dx, dy, dz;
float minx, maxx, miny, maxy, minz, maxz;
float dt;
float weight;
int   i, j, k, nx, ny, nz, t, m, i3,size,numbpts;//end;
float *vertPosition, *funcvalue_bak, *coefficent, *dfang;
float point[3];

vertPosition   =new float [vertex_Positions.size()];
size=vertex_Positions.size();
numbpts=size/3;
//printf("numbpts=%d",numbpts);getchar();
memcpy(vertPosition,&(vertex_Positions[0]), vertex_Positions.size()*sizeof(float));
//vertexPositions.clear();

//printf("xyzStart %f  %f  %f \n",xyzStart[0],xyzStart[1], xyzStart[2]);
//Switch(dim, minExt, maxExt);



nx = dim[0];
ny = dim[1];
nz = dim[2];

Dim[0] = nx;
Dim[1] = ny;
Dim[2] = nz;

minx = minExt[0];
maxx = maxExt[0];

miny = minExt[1];
maxy = maxExt[1]; 

minz = minExt[2];
maxz = maxExt[2];

printf("minxyz= %f, %f, %f  maxxyz = %f   %f   %f\n", minx, miny, minz, maxx, maxy, maxz); 
//getchar();
ComputeTensorXYZ();

dx = (maxx - minx)/(nx-1.0);
dy = (maxy - miny)/(ny-1.0);
dz = (maxz - minz)/(nz-1.0);

Dxyz[0] = dx;
Dxyz[1] = dy;
Dxyz[2] = dz;

// malloc memery
funcvalue_bak = new float [dim[0]*dim[1]*dim[2]];
coefficent    = new float [dim[0]*dim[1]*dim[2]];
dfang         = new float [dim[0]*dim[1]*dim[2]];
boundarysign  = new float [dim[0]*dim[1]*dim[2]];

Funcvalue_bak = new float [dim[0]*dim[1]*dim[2]];
Coefficent    = new float [dim[0]*dim[1]*dim[2]];
Funcvalue     = new float [dim[0]*dim[1]*dim[2]];
Dfang         = new float [dim[0]*dim[1]*dim[2]];

k = log(dx*dx*0.005/1.732)/log(0.26795);
//k = 8;
k = k + 2;
printf("k = %d\n", k);

//Height0 = (k-4)*dx;
Height = k*dx;
Height0 = Height/2.0;

printf("Height0, Height = %f, %f, dxyz = %f, %f %f,\n", Height0, Height,dx, dy, dz);

// Compute initial data
Fast_Distance_Function(dfang, vertPosition, numbpts, minx,  
                    maxx, miny, maxy, minz, maxz, nx, ny,  nz, edgelength);

int  load_data;


load_data = 0; 

//printf("Please change your load_data !\n");
//scanf("%d", &load_data);
//printf("load_data== %d\n", load_data);

if (load_data == 1) {
   float gamma=0.0;
  delete [] vertPosition;
   float * vertPosit;

   // load data
   FILE * file;
   //file = fopen("/home/xuguo/data/Raw_data/box1.raw", "r");
   file = fopen("/home/xuguo/data/data_load.raw", "r");
   //file = fopen("/home/xuguo/data/Tri_raw/icosa20480.raw", "r");
   int nV, nT;
   fscanf(file, "%d %d \n", &nV, &nT);

   vertPosit = new float [3*nV];
   for(i=0; i<nV; i++)   
	fscanf(file,"%f %f %f \n", vertPosit+3*i, vertPosit+3*i+1, vertPosit+3*i+2);
   
   int i0,i1,i2;
   for(i=0; i<nT; i++) 
   {
   fscanf(file, "%d %d %d \n", &i0, &i1, &i2);
   float distance0= Distance_Two_points3D(vertPosit+3*i0,vertPosit+3*i1);
   if(distance0>gamma) gamma= distance0;

   distance0=Distance_Two_points3D(vertPosit+3*i1,vertPosit+3*i2);
   if(distance0>gamma) gamma= distance0;

   distance0=Distance_Two_points3D(vertPosit+3*i0,vertPosit+3*i2);
   if(distance0>gamma) gamma= distance0;
   }

   fclose(file);
   Fast_Distance_Function(funcvalue_bak, vertPosit, nV, minx,  
                    maxx, miny, maxy, minz, maxz, nx, ny,  nz, gamma);

   InitialData_Scatter(funcvalue, funcvalue_bak, vertPosit, nV, minx,  
                    maxx, miny, maxy, minz, maxz, nx, ny,  nz, coefficent, gamma);

   delete [] vertPosit;

} else {
   InitialData_Scatter(funcvalue, dfang, vertPosition, numbpts, minx,  
                    maxx, miny, maxy, minz, maxz, nx, ny,  nz, coefficent, edgelength);
   delete [] vertPosition;
}


// for(i = 0; i < nx*ny*nz; i++)
//   if(funcvalue[i] < 0)  printf("funcvalue=%f ", funcvalue[i]);

/********************************************************************/ // Li add. for seperate the outer boundary points of HIV. 
/*
int *vertag;
int xi,yi,zi,ind,ind1;

vertag = (int *)malloc(numbpts*sizeof(int));
for (i = 0; i < numbpts; i++)
  vertag[i] = 0;


for(i = 0; i < numbpts; i++) {
   j = 3*i;
   point[0] = vertPosition[j  ];
   point[1] = vertPosition[j+1];
   point[2] = vertPosition[j+2];

   xi = (point[0]-minx)/dx;
   yi = (point[1]-miny)/dy;
   zi = (point[2]-minz)/dz;

   if( xi > nx/2.0 && yi > ny/2.0)  ind = ((xi-1)*ny+yi-1)*nz+zi;
   if(xi > nx/2.0 && yi < ny/2.0)   ind = ((xi-1)*ny+yi+1)*nz+zi;
   if( xi < nx/2.0 && yi > ny/2.0)  ind = ((xi+1)*ny+yi-1)*nz+zi;
   if( xi < nx/2.0 && yi < ny/2.0)  ind = ((xi+1)*ny+yi+1)*nz+zi;

  // ind1 = ((xi+1)*ny+yi)*nz+zi;  
   if(funcvalue[ind] > 0 ) {
     vertag[i] = 1;  
   }
}

 FILE * fp;
/* fp = fopen("/h1/liming/getOuterBounds/vertags.txt","w");
 for(i = 0; i < numbpts; i++)
   fprintf(fp,"%d\n", vertag[i]);
 fclose(fp);
 free(vertag);

 printf("\n tag finished.\n");
*/
/*******************************************************************/ // OVer.

// Taking a power p of distance function
float p;
p = 2;
for(i=0; i<nx*ny*nz; i++) {
   //funcvalue[i] = dfang[i]; // For Test distance function
   dfang[i] = pow(dfang[i], p);
}


//InitialData_Box(funcvalue, minx, maxx, miny, maxy, minz, maxz, nx, ny, nz);
//InitialData_IrBox(funcvalue, minx, maxx, miny, maxy, minz, maxz, nx, ny, nz);
//InitialData_Sphere(funcvalue, minx, maxx, miny, maxy, minz, maxz, nx, ny, nz);

printf("minx= %f, maxx=%f, maxy =%f  %f   %f   %f\n", minx, maxx, miny, maxy, minz,maxz); 

// Convert to spline function
for (i = 0; i < nx*ny*nz; i++) {
   coefficent[i] = funcvalue[i];
}
ConvertToInterpolationCoefficients_3D(coefficent, nx,ny,nz,CVC_DBL_EPSILON);
ConvertToInterpolationCoefficients_3D(dfang, nx,ny,nz,CVC_DBL_EPSILON);

// Evolution
TotalTime = 0.0;
//dt = 0.01;
// Iterarion times
weight=0.2;
weight=0.05*edgelength*edgelength;
weight=4*edgelength*edgelength;
//weight=0.0;
dt = dx*dx;
Dt = dt;
//dt=0.0;
//end = 350;
//end = 258;
//end = 1000;
//end = 4;
//end = 3000;
//end = 2;
 End = end;
 printf("end=%d ",end);
for (t = 0; t < end; t++) {

   printf("Interation for time = %d\n", t);
   if (t < 2)  MeanCurvatureFlow(funcvalue, funcvalue_bak, coefficent, dt, dx, dy, dz, nx, ny, nz, t);
   Constraint_MeanCurvatureFlow(dfang, funcvalue, funcvalue_bak,coefficent, weight, 
                              dt, dx,  dy,  dz,  nx,  ny,  nz,  t);
   
   for (i = 0; i < nx*ny*nz; i++) {
      coefficent[i] = funcvalue[i];
   }
   ConvertToInterpolationCoefficients_3D(coefficent, nx,ny,nz,CVC_DBL_EPSILON);


// Reinitilization
   //ReInitilazation(funcvalue, coefficent, dx, dy, dz, nx, ny, nz);
   //ReInitilazation_Upwind_Eno_Engquist(funcvalue, funcvalue_bak, coefficent, dx, dy, dz, nx, ny, nz);
   // Usually use the following one
   //if (t < end - 1) 
    ReInitilazation_Upwind_Eno_Godunov(funcvalue, funcvalue_bak, coefficent, dx, dy, dz, nx, ny, nz, t); 
   //ReInitilazation_Upwind_Eno_Godunov_Xu(funcvalue, funcvalue_bak, coefficent, dx, dy, dz, nx, ny, nz, t);
   //if (t == end - 1)
      //MeanCurvatureFlow(funcvalue, funcvalue_bak, coefficent, dt, dx, dy, dz, nx, ny, nz, t);
      //MeanCurvatureFlow(funcvalue, funcvalue_bak, coefficent, dt, dx, dy, dz, nx, ny, nz, t);
      //Fast_Reintialization(funcvalue, funcvalue_bak, coefficent, dim, minExt, maxExt);
      //MeanCurvatureFlow(funcvalue, funcvalue_bak, coefficent, dt, dx, dy, dz, nx, ny, nz, t);
      //MeanCurvatureFlow(funcvalue, funcvalue_bak, coefficent, dt, dx, dy, dz, nx, ny, nz, t);
}

  for (i = 0; i < nx*ny*nz; i++) {
      Coefficent[i] = coefficent[i];
      Funcvalue[i] = funcvalue[i];
      Funcvalue_bak[i] = funcvalue_bak[i];
      Dfang[i] = dfang[i];
  }


//  Constraint_MeanCurvatureFlow(dfang, funcvalue,  funcvalue_bak,coefficent,weight, dt, dx,  dy,  dz,  nx,  ny,  nz,  t);

//MeanCurvatureFlow(funcvalue, funcvalue_bak, coefficent, dt, dx, dy, dz, nx, ny, nz, t);
//MeanCurvatureFlow(funcvalue, funcvalue_bak, coefficent, dt, dx, dy, dz, nx, ny, nz, t);

float minf, maxf, w;

minf = 1000.0;
maxf = -1000.0;
for (i = 0; i < nx*ny*nz; i++) {
   if (funcvalue[i] > maxf) maxf = funcvalue[i];
   if (funcvalue[i] < minf) minf = funcvalue[i];
}

//printf("minf, maxf = %f, %f, Toltal Time = %f\n", minf, maxf, TotalTime);
printf("minf, maxf = %f, %f, Toltal Time = %f\n", minf, maxf, TotalTime);
printf("-offset = %f, 0 = %f\n", (-offset-minf)/(maxf - minf), -minf/(maxf - minf));
 
maxf = 1.0/(maxf - minf);
for (i = 0; i < nx*ny*nz; i++) {
   funcvalue[i] = (funcvalue[i] - minf)*maxf;
}

for (i = 0; i < nx*ny*nz; i++) {
   coefficent[i] = funcvalue[i];
}
ConvertToInterpolationCoefficients_3D(coefficent, nx,ny,nz,CVC_DBL_EPSILON);
printf("-offset1 = %f, 0 = %f\n", (-offset-minf)/(maxf - minf), -minf/(maxf - minf));

//ReSamplingCubicSpline(coefficent, nx, ny, nz, funcvalue, dim[0], dim[1], dim[2]);
/*
dim[0] = nx;
dim[1] = ny;
dim[2] = nz;
*/

/*
// Test convert 
ReSamplingCubicSpline(coefficent, nx, ny, nz, funcvalue_bak, nx, ny, nz);
maxf = -1000.0;
minf = 1000.0;
int  ii, j;
for (i = 0; i < nx*ny*nz; i++) {
   w = fabs(funcvalue[i]  - funcvalue_bak[i]);
   if (w > maxf) maxf = w;
   if (w < minf) minf = w;

   j = i/nx;
   k = i - j*nx;
   ii = j/nx;
   j = j - ii*nx;
   if (ii == 32 && j == 32) printf("i, j, k = %d, %d, %d, error = %e\n", ii, j, k, w); 
}
printf("Maximal error = %e, minimal error = %e\n", maxf, minf);
*/

// Data form transform from Xu form to TexMol form
int  index, index1;
for (i = 0; i < nx; i++) {
   for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
         index = (i*ny + j)*nz + k;
         index1 = (k*ny + j)*nx + i;
         funcvalue_bak[index1] = funcvalue[index];
      }
   }
}

printf("-offset2 = %f, 0 = %f\n", (-offset-minf)/(maxf - minf), -minf/(maxf - minf));

for (i = 0; i < nx*ny*nz; i++) {
   funcvalue[i]  =  funcvalue_bak[i];
}

printf("-offset3 = %f, 0 = %f\n", (-offset-minf)/(maxf - minf), -minf/(maxf - minf));
delete [] funcvalue_bak;
return true;
}


bool HLevelSet::computeFunction_Zhang_N(const std::vector<float> vertex_Positions, float* funcvalue, unsigned int* dim, float* minExt, float* maxExt, float& isovalue)
{
float dx, dy, dz;
float minx, maxx, miny, maxy, minz, maxz;
float dt;
float weight;
int   i, j, k, nx, ny, nz, t, m, i3,size,numbpts;//end;
float *vertPosition, *funcvalue_bak, *coefficent, *dfang;
float point[3];

vertPosition   =new float [vertex_Positions.size()];
size=vertex_Positions.size();
numbpts=size/3;
//printf("numbpts=%d",numbpts);getchar();
memcpy(vertPosition,&(vertex_Positions[0]), vertex_Positions.size()*sizeof(float));
//vertexPositions.clear();

//printf("xyzStart %f  %f  %f \n",xyzStart[0],xyzStart[1], xyzStart[2]);
//Switch(dim, minExt, maxExt);



nx = dim[0];
ny = dim[1];
nz = dim[2];


minx = minExt[0];
maxx = maxExt[0];

miny = minExt[1];
maxy = maxExt[1]; 

minz = minExt[2];
maxz = maxExt[2];

//printf("minxyz= %f, %f, %f  maxxyz = %f   %f   %f\n", minx, miny, minz, maxx, maxy, maxz); 
//getchar();
ComputeTensorXYZ();

dx = (maxx - minx)/(nx-1.0);
dy = (maxy - miny)/(ny-1.0);
dz = (maxz - minz)/(nz-1.0);

// malloc memery
funcvalue_bak = new float [dim[0]*dim[1]*dim[2]];
coefficent    = new float [dim[0]*dim[1]*dim[2]];
dfang         = new float [dim[0]*dim[1]*dim[2]];


k = log(dx*dx*0.005/1.732)/log(0.26795);
//k = 8;  //for brain data uncomment this line
k = k + 2;
k=8; // test for neuron

// printf("k = %d\n", k);

//Height0 = (k-4)*dx;
Height = k*dx;
Height0 = Height/2.0;

//printf("Height0, Height = %f, %f, dxyz = %f, %f %f,\n", Height0, Height,dx, dy, dz);

// Compute initial data
Fast_Distance_Function(dfang, vertPosition, numbpts, minx,  
                    maxx, miny, maxy, minz, maxz, nx, ny,  nz, edgelength);


InitialData_Scatter(funcvalue, dfang, vertPosition, numbpts, minx,  
                    maxx, miny, maxy, minz, maxz, nx, ny,  nz, coefficent, edgelength);

delete [] vertPosition;


// Taking a power p of distance function
float p;
p = 2;
for(i=0; i<nx*ny*nz; i++) {
   //funcvalue[i] = dfang[i]; // For Test distance function
   dfang[i] = pow(dfang[i], p);
}


//InitialData_Box(funcvalue, minx, maxx, miny, maxy, minz, maxz, nx, ny, nz);
//InitialData_IrBox(funcvalue, minx, maxx, miny, maxy, minz, maxz, nx, ny, nz);
//InitialData_Sphere(funcvalue, minx, maxx, miny, maxy, minz, maxz, nx, ny, nz);

//printf("minx= %f, maxx=%f, maxy =%f  %f   %f   %f\n", minx, maxx, miny, maxy, minz,maxz); 

// Convert to spline function
for (i = 0; i < nx*ny*nz; i++) {
   coefficent[i] = funcvalue[i];
}
ConvertToInterpolationCoefficients_3D(coefficent, nx,ny,nz,CVC_DBL_EPSILON);
ConvertToInterpolationCoefficients_3D(dfang, nx,ny,nz,CVC_DBL_EPSILON);

// Evolution
TotalTime = 0.0;
//dt = 0.01;
// Iterarion times
weight=0.2;
weight=0.05*edgelength*edgelength;
weight=0.01*edgelength*edgelength;
//weight=0.0;
dt = dx;
//dt = dx*dx;  //test for neuron




//printf("end=%d ",end);
for (t = 0; t < end; t++) {

 //  printf("Interation for time = %d\n", t);
   if (t < 2)  MeanCurvatureFlow(funcvalue, funcvalue_bak, coefficent, dt, dx, dy, dz, nx, ny, nz, t);
   Constraint_MeanCurvatureFlow(dfang, funcvalue, funcvalue_bak,coefficent, weight, 
                              dt, dx,  dy,  dz,  nx,  ny,  nz,  t);
   
   for (i = 0; i < nx*ny*nz; i++) {
      coefficent[i] = funcvalue[i];
   }
   ConvertToInterpolationCoefficients_3D(coefficent, nx,ny,nz,CVC_DBL_EPSILON);


// Reinitilization
   //ReInitilazation(funcvalue, coefficent, dx, dy, dz, nx, ny, nz);
   //ReInitilazation_Upwind_Eno_Engquist(funcvalue, funcvalue_bak, coefficent, dx, dy, dz, nx, ny, nz);
   // Usually use the following one
   if (t %3 ==0) 
    ReInitilazation_Upwind_Eno_Godunov(funcvalue, funcvalue_bak, coefficent, dx, dy, dz, nx, ny, nz, t); 
   //ReInitilazation_Upwind_Eno_Godunov_Xu(funcvalue, funcvalue_bak, coefficent, dx, dy, dz, nx, ny, nz, t);
   //if (t == end - 1)
      //MeanCurvatureFlow(funcvalue, funcvalue_bak, coefficent, dt, dx, dy, dz, nx, ny, nz, t);
      //MeanCurvatureFlow(funcvalue, funcvalue_bak, coefficent, dt, dx, dy, dz, nx, ny, nz, t);
      //Fast_Reintialization(funcvalue, funcvalue_bak, coefficent, dim, minExt, maxExt);
      //MeanCurvatureFlow(funcvalue, funcvalue_bak, coefficent, dt, dx, dy, dz, nx, ny, nz, t);
      //MeanCurvatureFlow(funcvalue, funcvalue_bak, coefficent, dt, dx, dy, dz, nx, ny, nz, t);
}



Constraint_MeanCurvatureFlow(dfang, funcvalue,  funcvalue_bak,coefficent,weight, dt, dx,  dy,  dz,  nx,  ny,  nz,  t);

//MeanCurvatureFlow(funcvalue, funcvalue_bak, coefficent, dt, dx, dy, dz, nx, ny, nz, t);
//MeanCurvatureFlow(funcvalue, funcvalue_bak, coefficent, dt, dx, dy, dz, nx, ny, nz, t);

float minf, maxf, w;

minf = 1000.0;
maxf = -1000.0;
for (i = 0; i < nx*ny*nz; i++) {
   if (funcvalue[i] > maxf) maxf = funcvalue[i];
   if (funcvalue[i] < minf) minf = funcvalue[i];
}

if (std::fabs(minf)<std::fabs(maxf)) offset *= 0.9;  //test for neuron

//printf("minf, maxf = %f, %f, Toltal Time = %f\n", minf, maxf, TotalTime);
//printf("minf, maxf = %f, %f, Toltal Time = %f\n", minf, maxf, TotalTime);
printf("-offset = %f, 0 = %f\n", (-offset-minf)/(maxf - minf), -minf/(maxf - minf));
isovalue=(-offset-minf)/(maxf - minf);


maxf = 1.0/(maxf - minf);
for (i = 0; i < nx*ny*nz; i++) {
   funcvalue[i] = (funcvalue[i] - minf)*maxf;
}

for (i = 0; i < nx*ny*nz; i++) {
   coefficent[i] = funcvalue[i];
}
ConvertToInterpolationCoefficients_3D(coefficent, nx,ny,nz,CVC_DBL_EPSILON);
//printf("-offset1 = %f, 0 = %f\n", (-offset-minf)/(maxf - minf), -minf/(maxf - minf));


// Data form transform from Xu form to TexMol form
int  index, index1;
for (i = 0; i < nx; i++) {
   for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
         index = (i*ny + j)*nz + k;
         index1 = (k*ny + j)*nx + i;
         funcvalue_bak[index1] = funcvalue[index];
      }
   }
}

//printf("-offset2 = %f, 0 = %f\n", (-offset-minf)/(maxf - minf), -minf/(maxf - minf));

for (i = 0; i < nx*ny*nz; i++) {
   funcvalue[i]  =  funcvalue_bak[i];
}


ConvertToInterpolationCoefficients_3D(funcvalue_bak, nx,ny,nz,CVC_DBL_EPSILON);


//printf("-offset3 = %f, 0 = %f\n", (-offset-minf)/(maxf - minf), -minf/(maxf - minf));
delete[] funcvalue_bak;
delete[] coefficent;
delete[] dfang;
return true;
}





/****************************************************************/
bool HLevelSet::computeFunction_Zhang_sdf(float* funcvalue, float edgelength, int end)
{
float dx, dy, dz;
float minx, maxx, miny, maxy, minz, maxz;
float dt;
float weight;
int   i, j, k, nx, ny, nz, t, m, i3,size,numbpts;//end;
float *vertPosition, *funcvalue_bak, *coefficent,*dfang;
float point[3];


nx = Dim[0];
ny = Dim[1];
nz = Dim[2];



//printf("minxyz= %f, %f, %f  maxxyz = %f   %f   %f\n", minx, miny, minz, maxx, maxy, maxz); 
//getchar();
//ComputeTensorXYZ();

dx = Dxyz[0];
dy = Dxyz[1];
dz = Dxyz[2];

printf("sdf dx=%f dy= %f dz=%f",dx,dy,dz);
// malloc memery
funcvalue_bak = new float [Dim[0]*Dim[1]*Dim[2]];
coefficent    = new float [Dim[0]*Dim[1]*Dim[2]];
dfang         = new float [Dim[0]*Dim[1]*Dim[2]];

k = log(dx*dx*0.005/1.732)/log(0.26795);
//k = 8;
k = k + 2;
printf("k = %d\n", k);

//Height0 = (k-4)*dx;
Height = k*dx;
Height0 = Height/2.0;

printf("Height0, Height = %f, %f, dxyz = %f, %f %f,\n", Height0, Height,dx, dy, dz);


// Evolution
TotalTime = 0.0;
//dt = 0.01;
// Iterarion times
weight=0.2;
weight=0.05*edgelength*edgelength;
weight=4*edgelength*edgelength;
//weight=0.0;
dt = Dt;

for (i = 0; i < nx*ny*nz; i++) {
  coefficent[i] = Coefficent[i];
  funcvalue[i] = Funcvalue[i];
  funcvalue_bak[i] = Funcvalue_bak[i];
  dfang[i] = Dfang[i];
}

 printf("end=%d ",end);
for (t = 0; t < end; t++) {

   printf("Interation for time = %d\n", t);
   if (End < 2)  MeanCurvatureFlow(funcvalue, funcvalue_bak, coefficent, dt, dx, dy, dz, nx, ny, nz, t);
   Constraint_MeanCurvatureFlow(dfang, funcvalue, funcvalue_bak,coefficent, weight, 
                              dt, dx,  dy,  dz,  nx,  ny,  nz,  t);
   
   for (i = 0; i < nx*ny*nz; i++) {
      coefficent[i] = funcvalue[i];
   }
   ConvertToInterpolationCoefficients_3D(coefficent, nx,ny,nz,CVC_DBL_EPSILON);
    ReInitilazation_Upwind_Eno_Godunov(funcvalue, funcvalue_bak, coefficent, dx, dy, dz, nx, ny, nz, t); 
}

for (i = 0; i < nx*ny*nz; i++) {
  Coefficent[i] = coefficent[i];
  Funcvalue[i] = funcvalue[i];
  Funcvalue_bak[i] = funcvalue_bak[i];
  Dfang[i] = dfang[i];
}

 End = End + end;
 printf("\n End = %d ",End);


float minf, maxf, w;

minf = 1000.0;
maxf = -1000.0;
for (i = 0; i < nx*ny*nz; i++) {
   if (funcvalue[i] > maxf) maxf = funcvalue[i];
   if (funcvalue[i] < minf) minf = funcvalue[i];
}

//printf("minf, maxf = %f, %f, Toltal Time = %f\n", minf, maxf, TotalTime);
printf("minf, maxf = %f, %f, Toltal Time = %f\n", minf, maxf, TotalTime);
printf("-offset = %f, 0 = %f\n", (-offset-minf)/(maxf - minf), -minf/(maxf - minf));
 
maxf = 1.0/(maxf - minf);
for (i = 0; i < nx*ny*nz; i++) {
   funcvalue[i] = (funcvalue[i] - minf)*maxf;
}

for (i = 0; i < nx*ny*nz; i++) {
   coefficent[i] = funcvalue[i];
}
ConvertToInterpolationCoefficients_3D(coefficent, nx,ny,nz,CVC_DBL_EPSILON);
printf("-offset1 = %f, 0 = %f\n", (-offset-minf)/(maxf - minf), -minf/(maxf - minf));


// Data form transform from Xu form to TexMol form
int  index, index1;
for (i = 0; i < nx; i++) {
   for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
         index = (i*ny + j)*nz + k;
         index1 = (k*ny + j)*nx + i;
         funcvalue_bak[index1] = funcvalue[index];
      }
   }
}

printf("-offset2 = %f, 0 = %f\n", (-offset-minf)/(maxf - minf), -minf/(maxf - minf));

for (i = 0; i < nx*ny*nz; i++) {
   funcvalue[i]  =  funcvalue_bak[i];
}

printf("-offset3 = %f, 0 = %f\n", (-offset-minf)/(maxf - minf), -minf/(maxf - minf));
delete [] funcvalue_bak;

return true;
}

//-----------------------------------------------------------------------------
/*void Fast_Reintialization(float *data, float *distance, float *coeff, unsigned int* dim, 
                          float* minExt, float* maxExt) 
{
   int   i, nx, ny, nz, numbtris, numbvert, m, c, c3;
   float isovalue, p1[3], p2[3], p3[3], dx, dy, dz;

   nx = dim[0];
   ny = dim[1];
   nz = dim[2];

   for(i=0; i < nx*ny*nz; i++) {
      coeff[i] = data[i];
   }

   dx = (maxExt[0] - minExt[0])/(nx - 1.0);
   dy = (maxExt[1] - minExt[1])/(ny - 1.0);
   dz = (maxExt[2] - minExt[2])/(nz - 1.0);

   SimpleVolumeData* sData = new SimpleVolumeData( dim );
   sData->setDimensions( dim );
   sData->setNumberOfVariables(1);
   sData->setData(0, data);
   sData->setType(0, SimpleVolumeData::FLOAT);
   sData->setName(0, "HOrderLevelSet");
   sData->setMinExtent(minExt);
   sData->setMaxExtent(maxExt);

   //isovalue = 0.5;
   isovalue = 0.0;
   // extract iso-surface
   printf("before extract iso-surface\n");
   Geometry* geometry = SimpleVolumeDataIsocontourer::getIsocontour(sData, isovalue);

   // remove sData
   //delete sData;

   numbvert = geometry->m_NumTriVerts;
   numbtris = geometry->m_NumTris;
   printf("NVert= %d, NTRi=%d in Fast_Reinitilization, after exact iso-surface\n", 
           numbvert, numbtris);

   // set initial values 
   for(i = 0; i < nx*ny*nz; i++) {
      distance[i]= 1000000.0;
   }

   // compute distance locally 
   m = 3;
   //m =Height/dx + 1;
   printf("Before initialization dxyz = %f, %f, %f,    m = %d\n", dx, dy, dz, m);


   int ii, jj, kk;
   for (c = 0; c < numbtris; c++) {
       c3 = c*3;
       ii = 3*geometry->m_Tris[c3];
       jj = 3*geometry->m_Tris[c3+1];
       kk = 3*geometry->m_Tris[c3+2];
       for (i = 0; i < 3; i++) {
          p1[i] = geometry->m_TriVerts[ii + i];
          p2[i] = geometry->m_TriVerts[jj + i];
          p3[i] = geometry->m_TriVerts[kk + i];
       }
       ComputeDistanceFunction_To_Triangle(distance, p1, p2, p3, m, nx, ny, nz, dx, dy, dz,
                               minExt[0], minExt[1], minExt[2]);
    }

   printf("Before fast sweeping dxyz = %f, %f, %f\n", dx, dy, dz);
   // set initial values
   for(i=0; i < nx*ny*nz; i++) {
      distance[i]= sqrt(distance[i]);
   }

   // Fast sweeping
   Fast_Sweeping(distance, nx, ny, nz, dx);

   printf("After fast sweeping \n");

   // set signs
    for(i=0; i < nx*ny*nz; i++) {
      if (coeff[i] >=  0.0) {
         data[i] = distance[i];
      }  else {
         data[i] = -distance[i];
      }

      // Cut 
      if (data[i] < -Height) data[i] = -Height;
      if (data[i] > Height) data[i] = Height;

      //data[i] = coeff[i];   //  For Test
      coeff[i] = data[i];
   }

   ConvertToInterpolationCoefficients_3D(coeff, nx,ny,nz,CVC_DBL_EPSILON);


   delete geometry;
}
*/
//-----------------------------------------------------------------------------
/*bool HLevelSet::getAtomListAndExtent( GroupOfAtoms* molecule, std::vector<Atom*> &atomList, float* minExtent, float* maxExtent )
{
	if( !molecule || !minExtent || !maxExtent ) return false;

	CollectionData* collectionData = 0;
	if( molecule->type == COLLECTION_TYPE ) collectionData = molecule->m_CollectionData;

	GroupOfAtoms::RADIUS_TYPE radiusType = GroupOfAtoms::VDW_RADIUS;

	BlurMapsDataManager::flattenGOA(molecule, atomList, collectionData, 0, 0, 0, radiusType, ATOM_TYPE, false  );

	minExtent[0] = minExtent[1] = minExtent[2] = 0.;
	maxExtent[0] = maxExtent[1] = maxExtent[2] = 0.;

	double probeRadius = 1.4;
	BlurMapsDataManager::getBoundingBox(atomList, minExtent, maxExtent, radiusType, probeRadius*2+1, 0);

	return true;
}
*/
//------------------------------------------------------------------------------
boost::tuple<bool,VolMagick::Volume> HLevelSet::getHigherOrderLevelSetSurface_Xu_Li( std::vector<float> vertex_Positions,unsigned int* dim , float edgelength, int end, int Max_dim )
{
	VolMagick::Volume result;

	// get the atom list from the molecule
	float minExt[3]={0.0,0.0,0.0};
	float maxExt[3]={0.0,0.0,0.0};
        float *data, dx,maxlevel=-1000000000.0,minlevel=100000000;
        float cx,cy,cz,minx,miny,minz,maxx,maxy,maxz;
        int i,j,c,size,numbpts;
        //int   Max_dim;
        
        maxlevel=-1000000000.0;
        minlevel=100000000.0;
        // Sizes
        //Max_dim =  256;
        //Max_dim =  400;
       // Max_dim =  150;
        //Max_dim =  100;
        //Max_dim =  400;


/*****************************/ // Li add for seperate the outer boundary points of HIV.
/*
size=vertex_Positions.size();
numbpts=size/3.0;


for (i=0; i<3; i++) {
    minExt[i] = vertex_Positions[i];
    maxExt[i] = vertex_Positions[i];

}


for (c=0; c<numbpts; c++) {
    for (i=0; i<3;i++) {
        if(minExt[i] > vertex_Positions[3*c+i]) minExt[i] = vertex_Positions[3*c+i];
        if(maxExt[i] < vertex_Positions[3*c+i]) maxExt[i] = vertex_Positions[3*c+i];
    }
}


minx = minExt[0];
maxx = maxExt[0];

miny = minExt[1];
maxy = maxExt[1];

minz = minExt[2];
maxz = maxExt[2];



//determine the center and radius   
cx = (maxx + minx)/2.0;
cy = (maxy + miny)/2.0;
cz = (maxz + minz)/2.0;

for(i = 0; i < numbpts; i++) {
    j = 3*i;
    vertex_Positions[j  ] = vertex_Positions[j  ] - cx;
    vertex_Positions[j+1] = vertex_Positions[j+1] - cy;
    vertex_Positions[j+2] = vertex_Positions[j+2] - cz;
}

printf("____min=%f %f %f max=%f %f %f \n centers=%f %f %f",minx,miny,minz,maxx,maxy,maxz,cx,cy,cz);
*/
/*************************************************************/ // over.

        Find_Box_Dim(Max_dim, vertex_Positions,minExt, maxExt, dim, edgelength);
        MinExt[0] = minExt[0]; MinExt[1] = minExt[1]; MinExt[2] = minExt[2];
        MaxExt[0] = maxExt[0]; MaxExt[1] = maxExt[1]; MaxExt[2] = maxExt[2];


        data = new float [dim[0]*dim[1]*dim[2]];
      
	// compute the function
	if( !computeFunction_Zhang(vertex_Positions,data, dim, minExt, maxExt,edgelength,end ) ) 
                             { delete []data; data = 0; return boost::make_tuple(false,result); }  

	
	result.voxelType(VolMagick::Float);
	result.dimension(VolMagick::Dimension(dim[0],dim[1],dim[2]));
	result.boundingBox(VolMagick::BoundingBox(minExt[0],minExt[1],minExt[2],
					          maxExt[0],maxExt[1],maxExt[2]));
        printf("minExt=%f %f %f maxExt=%f %f %f", minExt[0],minExt[1],minExt[2],
                                                  maxExt[0],maxExt[1],maxExt[2]);
       // getchar();
	for(VolMagick::uint64 i = 0; i < dim[0]; i++)
	  for(VolMagick::uint64 j = 0; j < dim[1]; j++)
            for(VolMagick::uint64 k = 0; k < dim[2]; k++) {
              result(i,j,k, data[k*dim[0]*dim[1]+j*dim[0]+i]);  
              if (data[k*dim[0]*dim[1]+j*dim[0]+i] > maxlevel) maxlevel = data[k*dim[0]*dim[1]+j*dim[0]+i]; //find the 0 levelset position of transfer function in volrover. 
              if (data[k*dim[0]*dim[1]+j*dim[0]+i] < minlevel) minlevel = data[k*dim[0]*dim[1]+j*dim[0]+i];

            } 

        printf("\n maxlevel=%f minlevel=%f \n",maxlevel,minlevel);   
	result.desc("HOrderLevelSet");
	return boost::make_tuple(true,result);
	
	// create volume data and return it.
/*	SimpleVolumeData* sData = new SimpleVolumeData( dim );
	sData->setDimensions( dim );
	sData->setNumberOfVariables(1);
	sData->setData(0, data);
	sData->setType(0, SimpleVolumeData::FLOAT);
	sData->setName(0, "HOrderLevelSet");
	sData->setMinExtent(minExt);
	sData->setMaxExtent(maxExt);

	return sData;
*/
}

boost::tuple<bool,VolMagick::Volume> HLevelSet::getHigherOrderLevelSetSurface_Xu_Li_N( std::vector<float> vertex_Positions,unsigned int* dim, VolMagick::BoundingBox& bb, float& isovalue)
{
	VolMagick::Volume result;

	// get the atom list from the molecule
	float minExt[3]={0.0,0.0,0.0};
	float maxExt[3]={0.0,0.0,0.0};
        float *data, dx,maxlevel=-1000000000.0,minlevel=100000000;
        float cx,cy,cz,minx,miny,minz,maxx,maxy,maxz;
        int i,j,c,size,numbpts;
        //int   Max_dim;
        
        maxlevel=-1000000000.0;
        minlevel=100000000.0;

        minExt[0]=bb.minx; maxExt[0]=bb.maxx;
        minExt[1]=bb.miny; maxExt[1]=bb.maxy;
        minExt[2]=bb.minz; maxExt[2]=bb.maxz;
        MinExt[0] = minExt[0]; MinExt[1] = minExt[1]; MinExt[2] = minExt[2]; 
        MaxExt[0] = maxExt[0]; MaxExt[1] = maxExt[1]; MaxExt[2] = maxExt[2];








        data = new float [dim[0]*dim[1]*dim[2]];
      
	// compute the function
	if( !computeFunction_Zhang_N(vertex_Positions,data, dim, minExt, maxExt,isovalue ) ) 
                             { delete []data; data = 0; return boost::make_tuple(false,result); }  

	
	result.voxelType(VolMagick::Float);
	result.dimension(VolMagick::Dimension(dim[0],dim[1],dim[2]));
	result.boundingBox(VolMagick::BoundingBox(minExt[0],minExt[1],minExt[2],
					          maxExt[0],maxExt[1],maxExt[2]));
//        printf("minExt=%f %f %f maxExt=%f %f %f", minExt[0],minExt[1],minExt[2],
  //                                                maxExt[0],maxExt[1],maxExt[2]);
       // getchar();
	for(VolMagick::uint64 i = 0; i < dim[0]; i++)
	  for(VolMagick::uint64 j = 0; j < dim[1]; j++)
            for(VolMagick::uint64 k = 0; k < dim[2]; k++) {
              result(i,j,k, data[k*dim[0]*dim[1]+j*dim[0]+i]);  
              if (data[k*dim[0]*dim[1]+j*dim[0]+i] > maxlevel) maxlevel = data[k*dim[0]*dim[1]+j*dim[0]+i]; //find the 0 levelset position of transfer function in volrover. 
              if (data[k*dim[0]*dim[1]+j*dim[0]+i] < minlevel) minlevel = data[k*dim[0]*dim[1]+j*dim[0]+i];

            } 

  //      printf("\n maxlevel=%f minlevel=%f \n",maxlevel,minlevel);   
	result.desc("HOrderLevelSet");
	delete[] data;
	return boost::make_tuple(true,result);
	
	// create volume data and return it.
/*	SimpleVolumeData* sData = new SimpleVolumeData( dim );
	sData->setDimensions( dim );
	sData->setNumberOfVariables(1);
	sData->setData(0, data);
	sData->setType(0, SimpleVolumeData::FLOAT);
	sData->setName(0, "HOrderLevelSet");
	sData->setMinExtent(minExt);
	sData->setMaxExtent(maxExt);

	return sData;
*/
}


boost::tuple<bool,VolMagick::Volume> HLevelSet::getHigherOrderLevelSetSurface_sdf(float edgelength, int end)
{
        VolMagick::Volume result;
        float minExt[3]={0.0,0.0,0.0};
        float maxExt[3]={0.0,0.0,0.0};
        float *data, dx,maxlevel=-1000000000.0,minlevel=100000000;
        float cx,cy,cz,minx,miny,minz,maxx,maxy,maxz;
        int i,j,c,size,numbpts,dim[3];
     
        dim[0] = Dim[0];dim[1] = Dim[1]; dim[2] = Dim[2];        
        minExt[0] = MinExt[0]; minExt[1] = MinExt[1]; minExt[2] = MinExt[2];
        maxExt[0] = MaxExt[0]; maxExt[1] = MaxExt[1]; maxExt[2] = MaxExt[2];


        printf("Dim=%d %d %d \n",Dim[0],Dim[1],Dim[2]);

        data = new float [Dim[0]*Dim[1]*Dim[2]];
        if( !computeFunction_Zhang_sdf(data,edgelength,end) )
                             { delete []data; data = 0; return boost::make_tuple(false,result); }



        result.voxelType(VolMagick::Float);
        result.dimension(VolMagick::Dimension(dim[0],dim[1],dim[2]));
        result.boundingBox(VolMagick::BoundingBox(minExt[0],minExt[1],minExt[2],
                                                  maxExt[0],maxExt[1],maxExt[2]));
        printf("minExt=%f %f %f maxExt=%f %f %f", minExt[0],minExt[1],minExt[2],
                                                  maxExt[0],maxExt[1],maxExt[2]);
        for(VolMagick::uint64 i = 0; i < dim[0]; i++)
          for(VolMagick::uint64 j = 0; j < dim[1]; j++)
            for(VolMagick::uint64 k = 0; k < dim[2]; k++) {
              result(i,j,k, data[k*dim[0]*dim[1]+j*dim[0]+i]);
            }

        result.desc("HOrderLevelSet");
        return boost::make_tuple(true,result);
}





bool HLevelSet::computeFunction_Xu_Li(float *vertPosition, int size, float *gridvalue, float* funcvalue, unsigned int* dim, float* minExt, float* maxExt,int *newobject,int numobject, float edgelength)
{
float dx, dy, dz;
float minx, maxx, miny, maxy, minz, maxz;
float dt;
float weight;
int   i, j, k, nx, ny, nz, t, m, i3,numbpts;//end;
float *funcvalue_bak, *coefficent, *dfang;
float *oneobject;
//vertPosition   =new float [vertex_Positions.size()];
//size=vertex_Positions.size();
//numbpts=size/3;
//printf("numbpts=%d",numbpts);getchar();
//memcpy(vertPosition,&(vertex_Positions[0]), vertex_Positions.size()*sizeof(float));
//vertexPositions.clear();

//printf("xyzStart %f  %f  %f \n",xyzStart[0],xyzStart[1], xyzStart[2]);
//Switch(dim, minExt, maxExt);



nx = dim[0];
ny = dim[1];
nz = dim[2];

minx = minExt[0];
maxx = maxExt[0];

miny = minExt[1];
maxy = maxExt[1]; 

minz = minExt[2];
maxz = maxExt[2];

printf("minxyz= %f, %f, %f  maxxyz = %f   %f   %f\n", minx, miny, minz, maxx, maxy, maxz); 
//getchar();
ComputeTensorXYZ();

dx = (maxx - minx)/(nx-1.0);
dy = (maxy - miny)/(ny-1.0);
dz = (maxz - minz)/(nz-1.0);

// malloc memery
funcvalue_bak = new float [dim[0]*dim[1]*dim[2]];
coefficent    = new float [dim[0]*dim[1]*dim[2]];
dfang         = new float [dim[0]*dim[1]*dim[2]];


k = log(dx*dx*0.005/1.732)/log(0.26795);
//k = 8;
k = k + 2;
printf("k = %d\n", k);

//Height0 = (k-4)*dx;
Height = k*dx;
Height0 = Height/2.0;

printf("Height0, Height = %f, %f, dxyz = %f, %f %f,\n", Height0, Height,dx, dy, dz);
edgelength = 1.0;
// Compute initial data
Fast_Distance_Function(dfang, vertPosition, size, minx,  
                    maxx, miny, maxy, minz, maxz, nx, ny,  nz, edgelength);




   InitialData_Scatter(funcvalue, dfang, vertPosition, size, minx,  
                    maxx, miny, maxy, minz, maxz, nx, ny,  nz, coefficent, edgelength);
   delete [] vertPosition;



   for(i = 0; i < nx*ny*nz; i++)
     if(funcvalue[i] < 0 ) gridvalue[i] = gridvalue[i] + 1;

   return true;
}



//------------------------------------------------------------------------------
/*SimpleVolumeData* HLevelSet::getHigherOrderLevelSetSurface_Zhang( unsigned int* dim )
{
	// get the atom list from the molecule
	std::vector<Atom*> atomList;
	float minExt[3];
	float maxExt[3];
        float *data, dx;
        int   Max_dim;
        
        // Sizes
        Max_dim =  256;
        Max_dim =  400;
        Max_dim =  150;
        Max_dim =  100;
        //Max_dim =  400;
        Find_Box_Dim(Max_dim, minExt, maxExt, dim);
        data = new float [dim[0]*dim[1]*dim[2]];

	// compute the function
	if( !computeFunction_sZhang(data, dim, minExt, maxExt ) ) 
                             { delete []data; data = 0; return 0; }

	// create volume data and return it.
	SimpleVolumeData* sData = new SimpleVolumeData( dim );
	sData->setDimensions( dim );
	sData->setNumberOfVariables(1);
	sData->setData(0, data);
	sData->setType(0, SimpleVolumeData::FLOAT);
	sData->setName(0, "HOrderLevelSet");
	sData->setMinExtent(minExt);
	sData->setMaxExtent(maxExt);

	return sData;
}
*/

/*-----------------------------------------------------------------------------*/
void  InitialData_Scatter(float *v, float *dSquare, float *init, int vnum, float minx, 
                          float maxx, float miny, float maxy, float minz, float maxz, 
                          int nx, int ny, int nz, float *signs, float edgelong)
//              float   *v,    /* the volume data of initial signed distance function */
//		float   *dSquare /* return the square of distance function            */
//  		float   *init   /* the scattered data coordinates		      */
//		int     vnum   /*  the given scattered data number  		      */
//              float   minx   /* left end-point of x diection                        */
//              float   maxx   /* right end-point of x diection                       */
//              float   miny   /* left end-point of y diection                        */
//              float   maxy   /* right end-point of y diection                       */
//              float   minz   /* left end-point of z diection                        */
//              float   maxz   /* right end-point of z diection                       */
//              int     nx,    /* number of points in x direction                     */
//              int     ny,    /* number of points in y direction                     */
//              int     nz,    /* number of points in z direction                     */
{
int    i, j, k, l, m, index, index1, p;
float  dx, dy, dz, x, y, z;
float  cx, cy, cz, radius, dis, point[3];

dx = (maxx - minx)/(nx - 1.0);
dy = (maxy - miny)/(ny - 1.0);
dz = (maxz - minz)/(nz - 1.0);

//printf("nx, ny, nz in InitialData_Scatter= %d, %d, %d\n", nx, ny, nz);
//determine the center and radius 
cx = (maxx + minx)/2.0;
cy = (maxy + miny)/2.0;
cz = (maxz + minz)/2.0;


float gridInterval=Max_Of_Three(maxx-minx, maxy-miny, maxz-minz);
float gridSize= Max_Of_Three(dx,dy,dz);


// compute signs
offset = edgelong;
offset = 0.4*edgelong;  // Test only
offset = 0.9*edgelong; 

if (offset < dx) offset = dx;  //for brain data comment this line

// set initial signs
for (i = 0; i < nx; i++) {
   for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
         index = (i*ny + j)*nz + k;
         signs[index] = 1.0;
         if (i > 0 && i < nx-1 && j > 0 && j < ny -1 && k > 0 && k < nz - 1) signs[index] = -1.0;
      }
   }
}
// search in z direction 
//printf("serarch in z direction\n");
for (i = 1; i < nx-1; i++) {
   for (j = 1; j < ny-1; j++) {
      for (k = 1; k < nz-1; k++) {
         index = (i*ny + j)*nz + k;
         if (Un_Treat_Boundary_Point(i, j, k,  signs, dSquare, offset, nx, ny, nz) ||
            signs[index] == 0.0 )   {
             signs[index] = 0.0; 
             break;
         }
         if (dSquare[index] > offset)  signs[index] = 1.0; 
      }

      for (k = nz-2; k > 0; k--) {
         index = (i*ny + j)*nz + k;
         if (Un_Treat_Boundary_Point(i, j, k,  signs, dSquare, offset, nx, ny, nz) ||
             signs[index] == 0.0 )   {
             signs[index] = 0.0; 
             break;
         }
         if (dSquare[index] > offset)  signs[index] = 1.0;
      }
   }
}

// search in y direction 
//printf("serarch in y direction\n");
for (i = 1; i < nx-1; i++) {
   for (k = 1; k < nz-1; k++) {
      for (j = 1; j < ny-1; j++) {
         index =  (i*ny + j)*nz + k;
         if (Un_Treat_Boundary_Point(i, j, k,  signs, dSquare, offset, nx, ny, nz) ||
             signs[index] == 0.0 )   {
             signs[index] = 0.0; 
             break;
         }
         if (dSquare[index] > offset)  signs[index] = 1.0;
      }

      for (j = ny-2; j > 0; j--) {
         index = (i*ny + j)*nz + k;
         if (Un_Treat_Boundary_Point(i, j, k,  signs, dSquare, offset, nx, ny, nz) ||
             signs[index] == 0.0 )   {
             signs[index] = 0.0;
             break;
         }
         if (dSquare[index] > offset)  signs[index] = 1.0;
      }
   }
}


// search in x direczztion
//printf("serarch in x direction\n");
for (k = 1; k < nz-1; k++) {
   for (j = 1; j < ny-1; j++) {
      for (i = 1; i < nx-1; i++) {
         index = (i*ny + j)*nz + k;
         if (Un_Treat_Boundary_Point(i, j, k,  signs, dSquare, offset, nx, ny, nz) ||
            signs[index] == 0.0 )   {
             signs[index] = 0.0;
             break;
         }
         if (dSquare[index] > offset)  signs[index] = 1.0;
      }

      for (i = nx-2; i > 0; i--) {
         index = (i*ny + j)*nz + k;
         if (Un_Treat_Boundary_Point(i, j, k,  signs, dSquare, offset, nx, ny, nz) ||
             signs[index] == 0.0 )   {
             signs[index] = 0.0;
             break;
         }
         if (dSquare[index] > offset)  signs[index] = 1.0;
      }
   }
}


// Compute remaining boundary points;
int  size, size2;

size = 0;
for (i = 0; i < nx; i++) {
   for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
         index = (i*ny + j)*nz + k;
         if (signs[index] == 0.0) {   // boundary point
             v[3*size]   = i;
             v[3*size+1] = j;
             v[3*size+2] = k;
             size = size + 1;
         }
      }
   }
}


//printf("Number of Boundary points = %d\n", size);
size2 = 0;

// Find remaining un-reconganized boundary points.
while (size > size2) {
   size2 = size;
   for (l = 0; l < size2; l++) {
   i = v[3*l] ;
   j = v[3*l+1];
   k = v[3*l+2];

   if (Un_Treat_Boundary_Point(i-1, j, k,  signs, dSquare, offset, nx, ny, nz)) {
      signs[((i-1)*ny + j)*nz + k] = 0.0;
      v[3*size] = i - 1;
      v[3*size+1] = j;
      v[3*size+2] = k;
      size = size + 1;
   } 
   
   if (Un_Treat_Boundary_Point(i+1, j, k,  signs, dSquare, offset, nx, ny, nz)) {
      signs[((i+1)*ny + j)*nz + k] = 0.0;
      v[3*size] = i + 1;
      v[3*size+1] = j;
      v[3*size+2] = k;
      size = size + 1;
   }

   if (Un_Treat_Boundary_Point(i, j-1, k,  signs, dSquare, offset, nx, ny, nz)) {
      signs[(i*ny + j-1)*nz + k] = 0.0;
      v[3*size] = i;
      v[3*size+1] = j - 1;
      v[3*size+2] = k;
      size = size + 1;
   }

   if (Un_Treat_Boundary_Point(i, j+1, k,  signs, dSquare, offset, nx, ny, nz)) {
      signs[(i*ny + j+1)*nz + k] = 0.0;
      v[3*size] = i;
      v[3*size+1] = j + 1;
      v[3*size+2] = k;
      size = size + 1;
   }

   if (Un_Treat_Boundary_Point(i, j, k-1, signs, dSquare, offset, nx, ny, nz)) {
      signs[(i*ny + j)*nz + k-1] = 0.0;
      v[3*size] = i;
      v[3*size+1] = j;
      v[3*size+2] = k - 1;
      size = size + 1;
   }

   if (Un_Treat_Boundary_Point(i, j, k+1, signs, dSquare, offset, nx, ny, nz)) {
      signs[(i*ny + j)*nz + k+1] = 0.0;
      v[3*size] = i;
      v[3*size+1] = j;
      v[3*size+2] = k + 1;
      size = size + 1;
   }
   }
//printf("size = %d after while loop \n", size);
}

size2 = 1;

// Find remaining un-reconganized external points.
while (size2 > 0) {
   size2 = 0;
   for (i = 1; i < nx - 1; i++) {
   for (j = 1; j < ny - 1; j++) {
   for (k = 1; k < nz - 1; k++) {
   if (signs[(i*ny + j)*nz + k] == 1.0) {

   index = ((i-1)*ny + j)*nz + k;
   if (dSquare[index] > offset && signs[index] == -1.0) {
      signs[index] = 1.0;
      size2 = size2 + 1;
   }

   index = ((i+1)*ny + j)*nz + k;
   if (dSquare[index] > offset && signs[index] == -1.0) {
      signs[index] = 1.0;
      size2 = size2 + 1;
   }

   index = (i*ny + j-1)*nz + k;
   if (dSquare[index] > offset && signs[index] == -1.0) {
      signs[index] = 1.0;
      size2 = size2 + 1;
   }

   index = (i*ny + j+1)*nz + k;
   if (dSquare[index] > offset && signs[index] == -1.0) {
      signs[index] = 1.0;
      size2 = size2 + 1;
   }

   index = (i*ny + j)*nz + k-1;
   if (dSquare[index] > offset && signs[index] == -1.0) {
      signs[index] = 1.0;
      size2 = size2 + 1;
   }

   index = (i*ny + j)*nz + k+1;
   if (dSquare[index] > offset && signs[index] == -1.0) {
      signs[index] = 1.0;
      size2 = size2 + 1;
   }

   }}}}
//printf("size2 = %d after while loop \n", size2);
}


// set initial values 
float w;
w = (Height + edgelong)*(Height + edgelong);
for(i=0; i<nx*ny*nz; i++) {
   v[i] = w;
}

// compute distancen
//printf("compute distance to boundary data\n");
//m = (Height + 2*edgelong)/dx + 1;
m = (Height + edgelong)/dx + 1;
m = (Height + 2.0*edgelong)/dx + 1;

//printf("m = %d\n", m);
for (i = 0; i < nx; i++) {
   for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
         index = (i*ny + j)*nz + k;
         if (signs[index] == 0.0) {
            //printf("Find boundary = index = %d\n", index);
            point[0] = minx + i*dx;
            point[1] = miny + j*dy;
            point[2] = minz + k*dz;
            ComputeDistanceFunction(v, point, m, nx, ny, nz, dx, dy, dz, minx, miny, minz);
         }
      }
   }
}

for(i=0; i<nx*ny*nz; i++) {
   v[i] = sqrt(v[i]);
}

// specify sign 
//printf("Specify signs\n");
for(i=0; i<nx*ny*nz; i++) {
   //v[i] = dSquare[i];
   //v[i] = dSquare[i]*signs[i];
   v[i] = v[i]*signs[i];
   //v[i]=v[i]+offset - dx;
//   v[i]=v[i]+offset;
   //v[i]=v[i];   // Test

   //dSquare[i] = dSquare[i]*dSquare[i];
   //dSquare[i] = pow(dSquare[i], p);

   if (v[i] > Height) v[i] = Height;
   if (v[i] < -Height) v[i] = -Height;

  // printf("v[i]= %d ", v[i]);
 
}
//getchar();

}





/*-----------------------------------------------------------------------------*/
void  Fast_Distance_Function(float *dSquare, float *init, int vnum, float minx, 
                          float maxx, float miny, float maxy, float minz, float maxz, 
                          int nx, int ny, int nz, float gamma)
//		float   *dSquare /* return the square of distance function            */
//  		float   *init   /* the scattered data coordinates		      */
//		int     vnum   /*  the given scattered data number  		      */
//              float   minx   /* left end-point of x diection                        */
//              float   maxx   /* right end-point of x diection                       */
//              float   miny   /* left end-point of y diection                        */
//              float   maxy   /* right end-point of y diection                       */
//              float   minz   /* left end-point of z diection                        */
//              float   maxz   /* right end-point of z diection                       */
//              int     nx,    /* number of points in x direction                     */
//              int     ny,    /* number of points in y direction                     */
//              int     nz,    /* number of points in z direction                     */
{
int    i, j, k, l, m, index, index1;
float  dx, dy, dz;
float  cx, cy, cz, point[3];

dx = (maxx - minx)/(nx - 1.0);
dy = (maxy - miny)/(ny - 1.0);
dz = (maxz - minz)/(nz - 1.0);

//printf("nx, ny, nz in InitialData_Scatter= %d, %d, %d\n", nx, ny, nz);
//determine the center and radius 
cx = (maxx + minx)/2.0;
cy = (maxy + miny)/2.0;
cz = (maxz + minz)/2.0;

// set initial values

float w;
w = 100000000.0;

for(i=0; i < nx*ny*nz; i++) {
   dSquare[i]= w;
}

// Compute distance to the scattered data
//printf("compute distance to sacttered data\n");
m = gamma/dx + 2;

//printf("m = %d in Fast_Distance_Function\n", m);
for(int vn=0; vn<vnum; vn++) {
   j = 3*vn;
   point[0] = init[j  ];
   point[1] = init[j+1];
   point[2] = init[j+2];
   ComputeDistanceFunction(dSquare, point, m, nx, ny, nz, dx, dy, dz, minx, miny, minz);
}

// set initial values
for(i=0; i < nx*ny*nz; i++) {
   dSquare[i]= sqrt(dSquare[i]);
}

Fast_Sweeping(dSquare, nx, ny, nz, dx);

}

/*-----------------------------------------------------------------------------*/
/* Fast sweeping, suppose initial data are ready                               */
/*-----------------------------------------------------------------------------*/
void Fast_Sweeping(float *dSquare, int nx, int ny, int nz, float dx)
// fast sweeping
{
int i, j, k;

for (i = 0; i < nx; i++) {
   for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
         One_Step_Sweeping(dSquare, nx, ny, nz, dx, i, j, k);
      }
   }
}

for (i = 0; i < nx; i++) {
   for (j = 0; j < ny; j++) {
      for (k = nz-1; k >= 0; k--) {
         One_Step_Sweeping(dSquare, nx, ny, nz, dx, i, j, k);
      }
   }
}

for (i = 0; i < nx; i++) {
   for (j = ny - 1; j >=0; j--) {
      for (k = 0; k < nz; k++) {
         One_Step_Sweeping(dSquare, nx, ny, nz, dx, i, j, k);
      }
   }
}

for (i = 0; i < nx; i++) {
   for (j = ny - 1; j >= 0; j--) {
      for (k = nz-1; k >= 0; k--) {
         One_Step_Sweeping(dSquare, nx, ny, nz, dx, i, j, k);
      }
   }
}


for (i = nx - 1; i >= 0; i--) {
   for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
         One_Step_Sweeping(dSquare, nx, ny, nz, dx, i, j, k);
      }
   }
}

for (i = nx - 1; i >=0; i--) {
   for (j = 0; j < ny; j++) {
      for (k = nz-1; k >= 0; k--) {
         One_Step_Sweeping(dSquare, nx, ny, nz, dx, i, j, k);
      }
   }
}

for (i = nx - 1; i >=0; i--) {
   for (j = ny - 1; j >= 0; j--) {
      for (k = 0; k < nz; k++) {
         One_Step_Sweeping(dSquare, nx, ny, nz, dx, i, j, k);
      }
   }
}

for (i = nx - 1; i >= 0; i--) {
   for (j = ny - 1; j >= 0; j--) {
      for (k = nz-1; k >= 0; k--) {
         One_Step_Sweeping(dSquare, nx, ny, nz, dx, i, j, k);
      }
   }
}

// one more time
for (i = 0; i < nx; i++) {
   for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
         One_Step_Sweeping(dSquare, nx, ny, nz, dx, i, j, k);
      }
   }
}

// Cut the function by Height
for(i=0; i<nx*ny*nz;i++)
{
   if(dSquare[i]>Height) dSquare[i]=Height;
}



}

/*-----------------------------------------------------------------------------------------*/
void One_Step_Sweeping(float *dSquare, int nx, int ny, int nz, float dx, int i, int j, int k)
{
float ddx, ddy, ddz, a, b, c, abc, dijk, h2dx, h2;

   h2   = dx*dx;
   h2dx = h2 + h2;

   if (i > 0 && i < nx - 1) {
      ddx = Min_Of_Two(dSquare[((i+1)*ny + j)*nz + k], dSquare[((i-1)*ny + j)*nz + k]);
   }  else {
      if(i==0) ddx= dSquare[(ny+j)*nz+k];
      if(i==nx-1)  ddx=dSquare[((nx-2)*ny+j)*nz+k];
   }

   if (j > 0 && j < ny - 1) {
      ddy = Min_Of_Two(dSquare[(i*ny + j+1)*nz + k],   dSquare[(i*ny + j-1)*nz + k]);
   }  else {
      if(j==0) ddy=dSquare[(i*ny+1)*nz+k];
      if(j==ny-1) ddy=dSquare[(i*ny+ny-2)*nz+k];
   }

   if (k > 0 && k < nz - 1) {
      ddz = Min_Of_Two(dSquare[(i*ny + j)*nz + k+1],   dSquare[(i*ny + j)*nz + k-1]);
   }  else {
      if(k==0) ddz=dSquare[(i*ny+j)*nz+1];
      if(k==nz-1) ddz=dSquare[(i*ny+j)*nz+nz-2];
   }

   a = Min_Of_Three(ddx,ddy,ddz);
   c = Max_Of_Three(ddx,ddy,ddz);
   abc = ddx + ddy + ddz;
   b = abc - a - c;

//printf("a, b, c %f, %f, %f\n", a, b, c);
   if (b - a >= dx) {
      dijk = a + dx;
   }  else {
      dijk = 0.5*(a + b + sqrt(h2dx - (b - a)*(b - a)));
      if (dijk >= c) {
         dijk = 0.3333333333*(abc + sqrt(3*(h2 - a*a - b*b - c*c) + abc*abc));
      }             
   }

   if (dSquare[(i*ny + j)*nz + k] > dijk) dSquare[(i*ny + j)*nz + k] = dijk;
   //dSquare[(i*ny + j)*nz + k] = dijk;
}








// New conversion code by Xuguo
/*-----------------------------------------------------------------------------*/
void ConvertToInterpolationCoefficients_1D(float *s, int DataLength, float EPSILON) 
//		float	*s,		/* input samples --> output coefficients */
//		int	DataLength,	/* number of samples or coefficients     */
//		float	EPSILON		/* admissible relative error             */

{
int   i, n, ni, ni1, K;
float sum, z1, w1, w2;

n = DataLength + 1;
z1 = sqrt(3.0) - 2.0;
K = log(EPSILON)/log(fabs(z1));
//printf("K = %i\n", K);

// compute initial value s(0)
sum = 0.0;
w2 = pow(z1, 2*n);
if (n < K) {
   for (i = 1; i < n; i++){
      w1 = pow(z1, i);
      sum = sum + s[i-1]*(w1 - w2/w1);
   }
} else {
   for (i = 1; i < n; i++){
      sum = (sum + s[n- i-1])*z1;
   }
}
sum = -sum/(1.0 - w2);


// compute c^{+}
n = DataLength;
s[0]  = s[0] + z1*sum;
for (i = 1; i < n; i++) {
   s[i]  = s[i] + z1*s[i-1];
   //printf("cp[%i] = %e, %f \n", i, cp[i], z1);
}

// compute c^- 
s[n-1] = -z1*s[n-1];
for (i = 1; i < n; i++) {
   ni = n - i; 
   ni1 = ni - 1;
   s[ni1]  = z1*(s[ni] - s[ni1]);
}

for (i = 0; i < n; i++) {
   s[i]  = 6.0*s[i];
}

}

/*-----------------------------------------------------------------------------*/
void Find_Box_Dim(int Max_dim, const std::vector<float> vertex_Positions,float* minExt, float* maxExt, unsigned int *dim, float edgelength)
{
int   i,j, m, c,size,numbpts,test;
float dx, xyzInterval[3];
float *vertPositions,ii,jj,kk,mindis,dis;
float cx,cy,cz,minx,miny,minz,maxx,maxy,maxz;
vertPositions   =new float [vertex_Positions.size()];

size=vertex_Positions.size();
numbpts=size/3.0;
//printf("numbpts=%d ",numbpts);
//getchar();
memcpy(vertPositions,&(vertex_Positions[0]), vertex_Positions.size()*sizeof(float));



 //Compute the Max and Min values of initial vertices
for (i=0; i<3; i++) {
    minExt[i] = vertPositions[i];
    maxExt[i] = vertPositions[i];

}

/*
edgelength=0.0;

for(i = 0 ; i < numbpts-1 ; i++) {

   ii = vertPositions[3*i+0];
   jj = vertPositions[3*i+1];
   kk = vertPositions[3*i+2];

   mindis=10000000; test=0;
   for(j = i + 1 ; j < numbpts ; j++ ) {
      dis = sqrt((ii-vertPositions[3*j+0])*(ii-vertPositions[3*j+0]) + 
                (jj-vertPositions[3*j+1])*(jj-vertPositions[3*j+1]) +
                (kk-vertPositions[3*j+2])*(kk-vertPositions[3*j+2]) );

     if(dis < mindis) mindis = dis;
     if(dis >= 10000000)  { printf("dis=%f i=%d j=%d ii=%f jj=%f kk=%f j1=%f j2=%f j3=%f",dis,i,j,ii,jj,kk,vertPositions[3*j+0],
                               vertPositions[3*j+1],vertPositions[3*j+2]); getchar();} 
     
     test=1;
  }
  if(mindis >= 10000000) printf("dis=%f i=%d j=%d ii=%f jj=%f kk=%f numbpts=%d mindis=%f test=%d", dis,i,j,ii,jj,kk,numbpts,mindis,test);  
  if (mindis > edgelength ) edgelength = mindis;

}

edgelength=.10;
*/
for (c=0; c<numbpts; c++) {
    for (int i=0; i<3;i++) {
        if(minExt[i] > vertPositions[3*c+i]) minExt[i] = vertPositions[3*c+i];
        if(maxExt[i] < vertPositions[3*c+i]) maxExt[i] = vertPositions[3*c+i];
    }
}
/*
 for(i = 0; i < 3; i++) {
   minExt[i] = minExt[i] - edgelength;
   maxExt[i] = maxExt[i] + edgelength; 
   
}
*/
printf("\n dd minExt=%f %f %f maxExt=%f %f %f", minExt[0],minExt[1],minExt[2],
                                                  maxExt[0],maxExt[1],maxExt[2]);




//Compute the exact interval widths in three directions
for (i=0; i<3; i++) {
    xyzInterval[i] = maxExt[i]-minExt[i];
}

// Maximal width among the three
float gridInterval;
gridInterval = Max_Of_Three(xyzInterval[0], xyzInterval[1], xyzInterval[2]);
float xyzStart[3], dd[3], scale[3];

dx = gridInterval/Max_dim;
m = 2*edgelength/dx + 3; 
if (m < 8) m = 8;
if (m > Max_dim/4) m = Max_dim/4;
// m = 8; 

/*
int m1;
m1 = -1;
while (m1 < 0) {
m1 = Min_Of_Three(dd[0] - 2.0*m, dd[1] - 2.0*m,dd[2] - 2.0*m); 
if (m1 < 0) m = m - 1;
}
m = m - 1;
*/
for (i=0; i < 3; i++) {
   dd[i] = xyzInterval[i]/dx;
   scale[i] = dd[i]/(dd[i] - 2.0*m);
}

//edgelength = 1.0;

//scale = 1.4;
// scale[0] = scale[1] = scale[2] = 1.4;

printf("\n edgelength=%f  scale = %f, %f, %f, m = %d\n", edgelength,scale[0], scale[1], scale[2], m);

//Compute the start point
for (int i=0;i<3;i++) {
    xyzStart[i]=(maxExt[i]+minExt[i])/2.0 - xyzInterval[i]*scale[i]/2.0;
}


//Compute the Box with minimal point minExt[i], maximal point maxExt[i]
for (int i=0;i<3;i++) {
   minExt[i]=xyzStart[i];
   xyzInterval[i] = xyzInterval[i]*scale[i];
   maxExt[i]=xyzStart[i]+xyzInterval[i];
}

printf("\n 1 minExt=%f %f %f maxExt=%f %f %f", minExt[0],minExt[1],minExt[2],
                                                  maxExt[0],maxExt[1],maxExt[2]);


gridInterval = Max_Of_Three(xyzInterval[0], xyzInterval[1], xyzInterval[2]);
dx  = gridInterval/(Max_dim - 1.0);

for (i=0; i < 3; i++) {
   dim[i] = (maxExt[i] - minExt[i])/dx + 0.5 + 1.0;
   maxExt[i] = minExt[i] + dx*(dim[i]- 1);
}

printf("\n 2 minExt=%f %f %f maxExt=%f %f %f", minExt[0],minExt[1],minExt[2],
                                                  maxExt[0],maxExt[1],maxExt[2]);

printf("dim = %d, %d, %d\n", dim[0], dim[1], dim[2]);
//Switch(dim, minExt, maxExt);


}
/*-----------------------------------------------------------------------------*/
void Switch(unsigned int *dim, float *minExt, float *maxExt)
{
int  i; 
float dx;

i = dim[0];
dim[0] = dim[2];
dim[2] = i;

dx = minExt[0];
minExt[0] = minExt[2];
minExt[2] = dx;

dx = maxExt[0];
maxExt[0] = maxExt[2];
maxExt[2] = dx;

}
/*-----------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------*/
void ConvertToInterpolationCoefficients_2D(float *s, int nx, int ny, float EPSILON) 
//		float	*s,     /* input samples --> output coefficients               */
//		int	nx,	/* number of samples or coefficients in x direction    */
//		int	ny,	/* number of samples or coefficients in y direction    */
//		float	EPSILON	/* admissible relative error                           */
{
float *d, z1;
int    i, l;

d = (float *) malloc (nx*sizeof (float)); 

z1 = sqrt(3.0) - 2.0;
// x-direction interpolation
for (l = 0; l < ny; l++) {
   for (i = 0; i < nx; i++) {
      d[i] = s[i*ny + l];
   }
   ConvertToInterpolationCoefficients_1D(d, nx,CVC_DBL_EPSILON);
   //ConvertToInterpolationCoefficients_Qu(d, nx,&z1, 1, CVC_DBL_EPSILON);


   for (i = 0; i < nx; i++) {
      s[i*ny + l] = d[i];
   }
}

// y-direction interpolation
for (i = 0; i < nx; i++) {
   ConvertToInterpolationCoefficients_1D(s+i*ny, ny,CVC_DBL_EPSILON);
   //ConvertToInterpolationCoefficients_Qu(s+i*ny, ny, &z1, 1, CVC_DBL_EPSILON);

}

free(d);
}


/*-----------------------------------------------------------------------------*/
void ConvertToInterpolationCoefficients_3D_Xu(float *s, int nx, int ny, 
                                                     int nz,   float EPSILON)
//              float   *s,    /* input samples --> output coefficients               */
//              int     nx,    /* number of samples or coefficients in x direction    */
//              int     ny,    /* number of samples or coefficients in y direction    */
//              int     nz,    /* number of samples or coefficients in z direction    */
//              float   EPSILON/* admissible relative error                           */
{
float *d, z1;
int    u, v, w, k, kk;

d = (float *) malloc (nx*sizeof (float));

k = ny*nz;
z1 = sqrt(3.0) - 2.0;

// x-direction interpolation
for (v = 0; v < ny; v++) {
   for (w = 0; w < nz; w++) {
      kk = v*nz + w;
      for (u = 0; u < nx; u++) {
         d[u] = s[u*k + kk];
      }
      ConvertToInterpolationCoefficients_1D(d, nx,CVC_DBL_EPSILON);
      //ConvertToInterpolationCoefficients_Qu(d, nx,&z1, 1, CVC_DBL_EPSILON);


      for (u = 0; u < nx; u++) {
         s[u*k + kk] = d[u];
      }
   }
}

for (u = 0; u < nx; u++) {
   ConvertToInterpolationCoefficients_2D(s+u*k, ny, nz, CVC_DBL_EPSILON);
}

free(d);
}


/*-----------------------------------------------------------------------------*/
void ConvertToInterpolationCoefficients_3D(float *s, int nx, int ny, int nz, float EPSILON)
{
  float *d = new float[nx*ny*nz];
  float sum;
    
  int i, j, k;
    
  float lbda[7]={0.0};
  lbda[1] = sqrt(3.0) - 2.0;
  lbda[2] = lbda[1]*lbda[1];
  lbda[3] = lbda[2]*lbda[1];
  lbda[4] = lbda[2]*lbda[2];
  lbda[5] = lbda[3]*lbda[2];
  lbda[6] = lbda[3]*lbda[3];

  for(i = 0; i < nx; i++)
  {
    for(j = 0; j < ny; j++)
    {
      for(k = 0; k < nz; k++)
      {
        sum = 0.0;
        for(int ii = 1; ii < 6; ii++)
        {
          if(k+ii < nx)     sum += s[(i*ny+j)*nz+k+ii]*lbda[ii];
          if(k-ii >=0 )     sum += s[(i*ny+j)*nz+k-ii]*lbda[ii];
        }
        d[(i*ny+j)*nz+k] = sqrt(3.0)*(s[(i*ny+j)*nz+k]+ sum);
      }
    }
  }

  for(i = 0; i < nx; i++)
  {
    for(j = 0; j < ny; j++)
    {
      for(k = 0; k < nz; k++)
      {
        sum = 0.0;
        for(int ii = 1; ii < 6; ii++)
        {
          if(j+ii < ny)     sum += d[(i*ny+j+ii)*nz+k]*lbda[ii];
          if(j-ii >=0 )     sum += d[(i*ny+j-ii)*nz+k]*lbda[ii];
        }
        s[(i*ny+j)*nz+k] = sqrt(3.0)*(d[(i*ny+j)*nz+k]+ sum);
      }
    }
  }

  for(i = 0; i < nx; i++)
  {
    for(j = 0; j < ny; j++)
    {
      for(k = 0; k < nz; k++)
      {
        sum = 0.0;
        for(int ii = 1; ii < 6; ii++)
        {
          if(i+ii < nz)     sum += s[((i+ii)*ny+j)*nz+k]*lbda[ii];
          if(i-ii >=0 )     sum += s[((i-ii)*ny+j)*nz+k]*lbda[ii];
        }
        d[(i*ny+j)*nz+k] = sqrt(3.0)*(s[(i*ny+j)*nz+k]+ sum);
      }
    }
  }

  for(i = 0; i < nx*ny*nz; i++)
  {
    s[i] = d[i];
  }
  delete []d;
}







/*-----------------------------------------------------------------------------*/
float   EvaluateCubicSplineAtGridPoint(float *c, int nx, int ny, int nz, int u, int v, int w)
//              float   *c,    /* the spline  coefficients                            */
//              int     nx,    /* number of samples or coefficients in x direction    */
//              int     ny,    /* number of samples or coefficients in y direction    */
//              int     nz,    /* number of samples or coefficients in z direction    */
{
float result;

if (( u  > 0  && u < nx-1) &&
    ( v  > 0  && v < ny-1) &&
    ( w  > 0  && w < nz-1) ) {

    result = OneSix3*(TakeACoefficient_Fast(c, nx, ny, nz, u-1,v-1,w-1) +
                      TakeACoefficient_Fast(c, nx, ny, nz, u-1,v-1,w+1) +
                      TakeACoefficient_Fast(c, nx, ny, nz, u-1,v+1,w+1) +
                      TakeACoefficient_Fast(c, nx, ny, nz, u-1,v+1,w-1) +
                      TakeACoefficient_Fast(c, nx, ny, nz, u+1,v-1,w-1) +
                      TakeACoefficient_Fast(c, nx, ny, nz, u+1,v-1,w+1) +
                      TakeACoefficient_Fast(c, nx, ny, nz, u+1,v+1,w+1) +
                      TakeACoefficient_Fast(c, nx, ny, nz, u+1,v+1,w-1)); 

   result = OneSix2*(TakeACoefficient_Fast(c, nx, ny, nz, u,v-1,w-1) +
                     TakeACoefficient_Fast(c, nx, ny, nz, u,v-1,w+1) +
                     TakeACoefficient_Fast(c, nx, ny, nz, u,v+1,w-1) +
                     TakeACoefficient_Fast(c, nx, ny, nz, u,v+1,w+1) 
                     +
                     TakeACoefficient_Fast(c, nx, ny, nz, u-1,v,w-1) +
                     TakeACoefficient_Fast(c, nx, ny, nz, u-1,v,w+1) +
                     TakeACoefficient_Fast(c, nx, ny, nz, u+1,v,w-1) +
                     TakeACoefficient_Fast(c, nx, ny, nz, u+1,v,w+1) 
                     +
                     TakeACoefficient_Fast(c, nx, ny, nz, u-1,v-1,w) +
                     TakeACoefficient_Fast(c, nx, ny, nz, u-1,v+1,w) +
                     TakeACoefficient_Fast(c, nx, ny, nz, u+1,v-1,w) +
                     TakeACoefficient_Fast(c, nx, ny, nz, u+1,v+1,w)) + result; 


   result = OneSix1*(TakeACoefficient_Fast(c, nx, ny, nz, u,v,w-1) +
                     TakeACoefficient_Fast(c, nx, ny, nz, u,v,w+1) +
                     TakeACoefficient_Fast(c, nx, ny, nz, u,v-1,w) +
                     TakeACoefficient_Fast(c, nx, ny, nz, u,v+1,w) +
                     TakeACoefficient_Fast(c, nx, ny, nz, u-1,v,w) +
                     TakeACoefficient_Fast(c, nx, ny, nz, u+1,v,w)) + result;

   result = OneSix0*TakeACoefficient_Fast(c, nx, ny, nz, u,v,w) + result;
   return(result);
} 

    result = OneSix3*(TakeACoefficient_Slow(c, nx, ny, nz, u-1,v-1,w-1) +
                      TakeACoefficient_Slow(c, nx, ny, nz, u-1,v-1,w+1) +
                      TakeACoefficient_Slow(c, nx, ny, nz, u-1,v+1,w+1) +
                      TakeACoefficient_Slow(c, nx, ny, nz, u-1,v+1,w-1) +
                      TakeACoefficient_Slow(c, nx, ny, nz, u+1,v-1,w-1) +
                      TakeACoefficient_Slow(c, nx, ny, nz, u+1,v-1,w+1) +
                      TakeACoefficient_Slow(c, nx, ny, nz, u+1,v+1,w+1) +
                      TakeACoefficient_Slow(c, nx, ny, nz, u+1,v+1,w-1)); 

   result = OneSix2*(TakeACoefficient_Slow(c, nx, ny, nz, u,v-1,w-1) +
                     TakeACoefficient_Slow(c, nx, ny, nz, u,v-1,w+1) +
                     TakeACoefficient_Slow(c, nx, ny, nz, u,v+1,w-1) +
                     TakeACoefficient_Slow(c, nx, ny, nz, u,v+1,w+1) 
                     +
                     TakeACoefficient_Slow(c, nx, ny, nz, u-1,v,w-1) +
                     TakeACoefficient_Slow(c, nx, ny, nz, u-1,v,w+1) +
                     TakeACoefficient_Slow(c, nx, ny, nz, u+1,v,w-1) +
                     TakeACoefficient_Slow(c, nx, ny, nz, u+1,v,w+1) 
                     +
                     TakeACoefficient_Slow(c, nx, ny, nz, u-1,v-1,w) +
                     TakeACoefficient_Slow(c, nx, ny, nz, u-1,v+1,w) +
                     TakeACoefficient_Slow(c, nx, ny, nz, u+1,v-1,w) +
                     TakeACoefficient_Slow(c, nx, ny, nz, u+1,v+1,w)) + result; 


   result = OneSix1*(TakeACoefficient_Slow(c, nx, ny, nz, u,v,w-1) +
                     TakeACoefficient_Slow(c, nx, ny, nz, u,v,w+1) +
                     TakeACoefficient_Slow(c, nx, ny, nz, u,v-1,w) +
                     TakeACoefficient_Slow(c, nx, ny, nz, u,v+1,w) +
                     TakeACoefficient_Slow(c, nx, ny, nz, u-1,v,w) +
                     TakeACoefficient_Slow(c, nx, ny, nz, u+1,v,w)) + result;

   result = OneSix0*TakeACoefficient_Slow(c, nx, ny, nz, u,v,w) + result;
   return(result);

}


/*-----------------------------------------------------------------------------*/
float TakeACoefficient_Fast(float *c, int nx, int ny, int nz, int u, int v, int w)
{
    return(c[(u*ny + v)*nz + w]);
}


/*-----------------------------------------------------------------------------*/
float TakeACoefficient_Slow(float *c, int nx, int ny, int nz, int u, int v, int w)
{
float result;

result = 0.0;
if (( u  >= 0  && u < nx) &&
    ( v  >= 0  && v < ny) &&
    ( w  >= 0  && w < nz) ) {

    result = c[(u*ny + v)*nz + w];
}

return(result);
}


/*-----------------------------------------------------------------------------*/
void   Tensor_333(float *xx, float *yy, float *zz, float *result)
{
int i, j, k, l;

l = 0;
for (i = 0; i < 3; i++) {
   for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
          result[l] = xx[i]*yy[j]*zz[k];
          l = l + 1;
      }
   }
}  

}

/*-----------------------------------------------------------------------------*/
void   Tensor_444(float *xx, float *yy, float *zz, float *result)
{
int i, j, k, l;
                                                                                                                 
l = 0;
for (i = 0; i < 4; i++) {
   for (j = 0; j < 4; j++) {
      for (k = 0; k < 4; k++) {
          result[l] = xx[i]*yy[j]*zz[k];
          l = l + 1;
      }
   }
}
                                                                                                                 
}
                                                                                                                 

   
/*-----------------------------------------------------------------------------*/
void Take_27_Coefficients(float *c, int nx, int ny, int nz, int u, int v, int w, float *c27)
//              float   *c,    /* the spline  coefficients                            */
//              int     nx,    /* number of samples or coefficients in x direction    */
//              int     ny,    /* number of samples or coefficients in y direction    */
//              int     nz,    /* number of samples or coefficients in z direction    */
//              float   *c27   /* 27 coefficients                                     */
{
if (( u  > 0  && u < nx-1) &&
    ( v  > 0  && v < ny-1) &&
    ( w  > 0  && w < nz-1) ) {

    c27[0] = TakeACoefficient_Fast(c, nx, ny, nz, u-1,v-1,w-1);
    c27[1] = TakeACoefficient_Fast(c, nx, ny, nz, u-1,v-1,w);
    c27[2] = TakeACoefficient_Fast(c, nx, ny, nz, u-1,v-1,w+1);

    c27[3] = TakeACoefficient_Fast(c, nx, ny, nz, u-1,v,w-1);
    c27[4] = TakeACoefficient_Fast(c, nx, ny, nz, u-1,v,w);
    c27[5] = TakeACoefficient_Fast(c, nx, ny, nz, u-1,v,w+1);
                      
    c27[6] = TakeACoefficient_Fast(c, nx, ny, nz, u-1,v+1,w-1);
    c27[7] = TakeACoefficient_Fast(c, nx, ny, nz, u-1,v+1,w);
    c27[8] = TakeACoefficient_Fast(c, nx, ny, nz, u-1,v+1,w+1);

    c27[9]  = TakeACoefficient_Fast(c, nx, ny, nz, u,v-1,w-1);
    c27[10] = TakeACoefficient_Fast(c, nx, ny, nz, u,v-1,w);
    c27[11] = TakeACoefficient_Fast(c, nx, ny, nz, u,v-1,w+1);

    c27[12] = TakeACoefficient_Fast(c, nx, ny, nz, u,v,w-1);
    c27[13] = TakeACoefficient_Fast(c, nx, ny, nz, u,v,w);
    c27[14] = TakeACoefficient_Fast(c, nx, ny, nz, u,v,w+1);
                      
    c27[15] = TakeACoefficient_Fast(c, nx, ny, nz, u,v+1,w-1);
    c27[16] = TakeACoefficient_Fast(c, nx, ny, nz, u,v+1,w);
    c27[17] = TakeACoefficient_Fast(c, nx, ny, nz, u,v+1,w+1);

    c27[18] = TakeACoefficient_Fast(c, nx, ny, nz, u+1,v-1,w-1);
    c27[19] = TakeACoefficient_Fast(c, nx, ny, nz, u+1,v-1,w);
    c27[20] = TakeACoefficient_Fast(c, nx, ny, nz, u+1,v-1,w+1);

    c27[21] = TakeACoefficient_Fast(c, nx, ny, nz, u+1,v,w-1);
    c27[22] = TakeACoefficient_Fast(c, nx, ny, nz, u+1,v,w);
    c27[23] = TakeACoefficient_Fast(c, nx, ny, nz, u+1,v,w+1);
                      
    c27[24] = TakeACoefficient_Fast(c, nx, ny, nz, u+1,v+1,w-1);
    c27[25] = TakeACoefficient_Fast(c, nx, ny, nz, u+1,v+1,w);
    c27[26] = TakeACoefficient_Fast(c, nx, ny, nz, u+1,v+1,w+1);

    return;
}

    c27[0] = TakeACoefficient_Slow(c, nx, ny, nz, u-1,v-1,w-1);
    c27[1] = TakeACoefficient_Slow(c, nx, ny, nz, u-1,v-1,w);
    c27[2] = TakeACoefficient_Slow(c, nx, ny, nz, u-1,v-1,w+1);

    c27[3] = TakeACoefficient_Slow(c, nx, ny, nz, u-1,v,w-1);
    c27[4] = TakeACoefficient_Slow(c, nx, ny, nz, u-1,v,w);
    c27[5] = TakeACoefficient_Slow(c, nx, ny, nz, u-1,v,w+1);
                      
    c27[6] = TakeACoefficient_Slow(c, nx, ny, nz, u-1,v+1,w-1);
    c27[7] = TakeACoefficient_Slow(c, nx, ny, nz, u-1,v+1,w);
    c27[8] = TakeACoefficient_Slow(c, nx, ny, nz, u-1,v+1,w+1);

    c27[9]  = TakeACoefficient_Slow(c, nx, ny, nz, u,v-1,w-1);
    c27[10] = TakeACoefficient_Slow(c, nx, ny, nz, u,v-1,w);
    c27[11] = TakeACoefficient_Slow(c, nx, ny, nz, u,v-1,w+1);

    c27[12] = TakeACoefficient_Slow(c, nx, ny, nz, u,v,w-1);
    c27[13] = TakeACoefficient_Slow(c, nx, ny, nz, u,v,w);
    c27[14] = TakeACoefficient_Slow(c, nx, ny, nz, u,v,w+1);
                      
    c27[15] = TakeACoefficient_Slow(c, nx, ny, nz, u,v+1,w-1);
    c27[16] = TakeACoefficient_Slow(c, nx, ny, nz, u,v+1,w);
    c27[17] = TakeACoefficient_Slow(c, nx, ny, nz, u,v+1,w+1);

    c27[18] = TakeACoefficient_Slow(c, nx, ny, nz, u+1,v-1,w-1);
    c27[19] = TakeACoefficient_Slow(c, nx, ny, nz, u+1,v-1,w);
    c27[20] = TakeACoefficient_Slow(c, nx, ny, nz, u+1,v-1,w+1);

    c27[21] = TakeACoefficient_Slow(c, nx, ny, nz, u+1,v,w-1);
    c27[22] = TakeACoefficient_Slow(c, nx, ny, nz, u+1,v,w);
    c27[23] = TakeACoefficient_Slow(c, nx, ny, nz, u+1,v,w+1);
                      
    c27[24] = TakeACoefficient_Slow(c, nx, ny, nz, u+1,v+1,w-1);
    c27[25] = TakeACoefficient_Slow(c, nx, ny, nz, u+1,v+1,w);
    c27[26] = TakeACoefficient_Slow(c, nx, ny, nz, u+1,v+1,w+1);
}

/*-----------------------------------------------------------------------------*/
void Take_64_Coefficients(float *c, int nx, int ny, int nz, int u, int v, int w, float *c64)
//              float   *c,    /* the spline  coefficients                            */
//              int     nx,    /* number of samples or coefficients in x direction    */
//              int     ny,    /* number of samples or coefficients in y direction    */
//              int     nz,    /* number of samples or coefficients in z direction    */
//              float   *c64   /* 64 coefficients                                     */
{
int i, j, k, l;


if (( u  > 0  && u < nx-2) &&
    ( v  > 0  && v < ny-2) &&
    ( w  > 0  && w < nz-2) ) {

    l = 0;
    for (i = 0; i < 4; i++) {
       for (j = 0; j < 4; j++) {
          for (k = 0; k < 4; k++) {
            c64[l] = TakeACoefficient_Fast(c, nx, ny, nz, u-1+i,v-1+j,w-1+k);
            l = l + 1;
          }
      }
   }
   return;
}

l = 0;
for (i = 0; i < 4; i++) {
   for (j = 0; j < 4; j++) {
       for (k = 0; k < 4; k++) {
         c64[l] = TakeACoefficient_Slow(c, nx, ny, nz, u-1+i,v-1+j,w-1+k);
         l = l + 1;
       }
   }
}

}


/*-----------------------------------------------------------------------------*/
void  EvaluateCubicSplineOrder2PartialsAtGridPoint(float *c, float dx, float dy, float dz,
                                             int nx, int ny, int nz, 
                                             int u, int v, int w, float *partials)
//              float   *c,    /* the spline  coefficients                            */
//              float   dx     /* spacing in x direction                              */
//              float   dy     /* spacing in y direction                              */
//              float   dz     /* spacing in z direction                              */
//              int     nx,    /* number of samples or coefficients in x direction    */
//              int     ny,    /* number of samples or coefficients in y direction    */
//              int     nz,    /* number of samples or coefficients in z direction    */
//              float   *partial /* partial derivatives                               */
{
float c27[27], indx, indy, indz;
int  i;

Take_27_Coefficients(c, nx, ny, nz, u, v, w, c27);

for (i = 0; i < 10; i++) {
   partials[i] = 0.0;
}

for (i = 0; i < 27; i++) {
   partials[0] = partials[0] + c27[i]*TensorF[i];

   partials[1] = partials[1] + c27[i]*TensorFx[i];
   partials[2] = partials[2] + c27[i]*TensorFy[i];
   partials[3] = partials[3] + c27[i]*TensorFz[i];

   partials[4] = partials[4] + c27[i]*TensorFxx[i];
   partials[5] = partials[5] + c27[i]*TensorFxy[i];
   partials[6] = partials[6] + c27[i]*TensorFxz[i];

   partials[7] = partials[7] + c27[i]*TensorFyy[i];
   partials[8] = partials[8] + c27[i]*TensorFyz[i];

   partials[9] = partials[9] + c27[i]*TensorFzz[i];

   //printf("i = %d   Coeff = %f,     Fxx = %f\n", i, c27[i], TensorFxx[i]);
}      

indx = 1.0/dx;
indy = 1.0/dy;
indz = 1.0/dz;

//printf("indx = %f\n", partials[4]);

partials[1] = partials[1] * indx;
partials[2] = partials[2] * indy;
partials[3] = partials[3] * indz;

partials[4] = partials[4] * indx*indx;
partials[5] = partials[5] * indx*indy;
partials[6] = partials[6] * indx*indz;

partials[7] = partials[7] * indy*indy;
partials[8] = partials[8] * indy*indz;

partials[9] = partials[9] * indz*indz;

//printf("indx = %f\n", partials[4]);

}

/*-----------------------------------------------------------------------------*/
void  EvaluateCubicSplineOrder1PartialsAtGridPoint(float *c, float dx, float dy, float dz,
                                             int nx, int ny, int nz,
                                             int u, int v, int w, float *partials)
//              float   *c,    /* the spline  coefficients                            */
//              float   dx     /* spacing in x direction                              */
//              float   dy     /* spacing in y direction                              */
//              float   dz     /* spacing in z direction                              */
//              int     nx,    /* number of samples or coefficients in x direction    */
//              int     ny,    /* number of samples or coefficients in y direction    */
//              int     nz,    /* number of samples or coefficients in z direction    */
//              float   *partial /* partial derivatives                               */
{
float c27[27], indx, indy, indz;
int  i;

Take_27_Coefficients(c, nx, ny, nz, u, v, w, c27);

for (i = 0; i < 4; i++) {
   partials[i] = 0.0;
}

for (i = 0; i < 27; i++) {
   partials[0] = partials[0] + c27[i]*TensorF[i];

   partials[1] = partials[1] + c27[i]*TensorFx[i];
   partials[2] = partials[2] + c27[i]*TensorFy[i];
   partials[3] = partials[3] + c27[i]*TensorFz[i];
}     

indx = 1.0/dx;
indy = 1.0/dy;
indz = 1.0/dz;

//printf("indx = %f, %f, %f\n", indx, indy, indz);

partials[1] = partials[1] * indx;
partials[2] = partials[2] * indy;
partials[3] = partials[3] * indz;

}




//------------------------------------------------------------------------------
// return true if it is a untreat boundary point
bool Un_Treat_Boundary_Point(int i, int j, int k, float *signs, float *dis, float offset, 
                             int nx, int ny, int nz)
{
bool result;
int index, index1, index2, index3, index4, index5, index6;

result = false;
if (i < 1 || i >= nx-1) {
//    printf("i is out off the range\n");
    return(result);
}
if (j < 1 || j >= ny-1) {
//    printf("j is out off the range\n");
    return(result);
}
if (k < 1 || k >= nz-1) {
 //   printf("k is out off the range\n");
    return(result);
}

result = false;
index = (i*ny + j)*nz + k;
if (signs[index] != 0.0) {

index1 = ((i+1)*ny + j)*nz + k;
index2 = ((i-1)*ny + j)*nz + k;
index3 = (i*ny + j+1)*nz + k;
index4 = (i*ny + j-1)*nz + k;
index5 = (i*ny + j)*nz + k + 1;
index6 = (i*ny + j)*nz + k - 1;

if ((dis[index] > offset) && (
   (dis[index1] <= offset) || (dis[index2] <= offset) ||
   (dis[index3] <= offset) || (dis[index4] <= offset) ||
   (dis[index5] <= offset) || (dis[index6] <= offset) ) ) result = true;

}

return(result);
}

//-----------------------------------------------------------------------------
void ComputeDistanceFunction(float *distanceD,  float *init,  int wid, 
                             int nnx, int nny, int nnz, float ddx, float ddy, 
                             float ddz, float xS, float yS, float zS)
{
int i, j, k, ii, jj, kk, m1,m2, m3, mm, nyz;
int k1, k2, j1, j2, i1, i2;
float x, y, z, cx, cy, cz, xx, yy, zz;
float distTemp;

// The given point 
cx = init[0];
cy = init[1];
cz = init[2];

// The nearest point indices
ii=(cx-xS)/ddx + 0.5;
jj=(cy-yS)/ddy + 0.5;
kk=(cz-zS)/ddz + 0.5;

nyz = nny*nnz;

// Range of k
k1 = kk-wid+1;
if (k1 < 0) k1 = 0;
k2 = kk+wid;
if (k2 > nnz) k2 = nnz;

// Range of j
j1 = jj-wid+1;
if (j1 < 0) j1 = 0;
j2 = jj+wid;
if (j2 > nny) j2 = nny;

// Range of i
i1 = ii-wid+1;
if (i1 < 0) i1 = 0;
i2 = ii+wid; 
if (i2 > nnx) i2 = nnx;

// printf("\n k1=%d k2=%d j1=%d j2=%d i1=%d i2=%d",k1,k2,j1,j2,i1,i2);

xx = xS - cx;
yy = yS - cy;	
zz = zS - cz;
// loop in a cube
for (k = k1; k < k2; k++) {
   z = zz + k*ddz;
   z = z*z;
   for (j = j1; j < j2; j++) {
      y = yy + j*ddy;
      y = y*y + z;
      m2=j*nnz + k;
      for(i = i1; i < i2; i++) {
         x = xx + i*ddx;
         x = x*x;
         mm = m2 + i*nyz;
         //distTemp= Distance_Two_points3D(position, init);
         distTemp= x + y;
         if (distanceD[mm] > distTemp) distanceD[mm] = distTemp;
      }
   }
} 

}

//-----------------------------------------------------------------------------
void ComputeDistanceFunction_To_Triangle(float *distanceD,  float *p1, float *p2, 
                             float *p3,  int wid, int nnx, int nny, int nnz, 
                             float ddx, float ddy, float ddz, 
                             float xS,  float yS,  float zS)
{
int i, j, k, ii1, jj1, kk1, ii2, jj2, kk2, ii3, jj3, kk3, m1,m2, m3, mm, nyz;
int k1, k2, j1, j2, i1, i2;
int mini, minj, mink, maxi, maxj, maxk;
float distTemp, p[3];

// The nearest point indices for the three points
ii1=(p1[0]-xS)/ddx + 0.5;
jj1=(p1[1]-yS)/ddy + 0.5;
kk1=(p1[2]-zS)/ddz + 0.5;

ii2=(p2[0]-xS)/ddx + 0.5;
jj2=(p2[1]-yS)/ddy + 0.5;
kk2=(p2[2]-zS)/ddz + 0.5;

ii3=(p3[0]-xS)/ddx + 0.5;
jj3=(p3[1]-yS)/ddy + 0.5;
kk3=(p3[2]-zS)/ddz + 0.5;

// Find minimal indices
mini = ii1;
minj = jj1;
mink = kk1;
if (mini > ii2) mini = ii2;
if (mini > ii3) mini = ii3;

if (minj > jj2) minj = jj2;
if (minj > jj3) minj = jj3;

if (mink > kk2) mink = kk2;
if (mink > kk3) mink = kk3;

// Find maximal indices
maxi = ii1;
maxj = jj1;
maxk = kk1;
if (maxi < ii2) maxi = ii2;
if (maxi < ii3) maxi = ii3;

if (maxj < jj2) maxj = jj2;
if (maxj < jj3) maxj = jj3;

if (maxk < kk2) maxk = kk2;
if (maxk < kk3) maxk = kk3;

// Range of k
k1 = mink-wid;
if (k1 < 0) k1 = 0;
k2 = maxk+wid;
if (k2 > nnz) k2 = nnz;

// Range of j
j1 = minj-wid;
if (j1 < 0) j1 = 0;
j2 = maxj+wid;
if (j2 > nny) j2 = nny;

// Range of i
i1 = mini-wid;
if (i1 < 0) i1 = 0;
i2 = maxi+wid; 
if (i2 > nnx) i2 = nnx;

nyz = nny*nnz;

// loop in a cube determined
for (k = k1; k < k2; k++) {
   p[2] = zS + k*ddz;
   for (j = j1; j < j2; j++) {
      p[1] = yS + j*ddy;
      m2=j*nnz + k;
      for(i = i1; i < i2; i++) {

         p[0] = xS + i*ddx;
         mm = m2 + i*nyz;
         distTemp= PointToTriangleDistance(p, p1, p2, p3);
         if (distanceD[mm] > distTemp) distanceD[mm] = distTemp;

      }
   }
} 

}


/*-----------------------------------------------------------------------------*/
void  InitialData_Sphere(float *v, float minx, float maxx, float miny, float maxy, 
                  float minz, float maxz, int nx, int ny, int nz)
//              float   *v,    /* the volume data of initial signed distance function */
//              float   minx   /* left end-point of x diection                        */
//              float   maxx   /* right end-point of x diection                       */
//              float   miny   /* left end-point of y diection                        */
//              float   maxy   /* right end-point of y diection                       */
//              float   minz   /* left end-point of z diection                        */
//              float   maxz   /* right end-point of z diection                       */
//              int     nx,    /* number of points in x direction                     */
//              int     ny,    /* number of points in y direction                     */
//              int     nz,    /* number of points in z direction                     */
{
int    i, j, k, l;
float  dx, dy, dz, x, y, z;
float  cx, cy, cz, radius, dis;

dx = (maxx - minx)/(nx - 1.0);
dy = (maxy - miny)/(ny - 1.0);
dz = (maxz - minz)/(nz - 1.0);

//determine the center and radius 
cx = (maxx + minx)/2.0;
cy = (maxy + miny)/2.0;
cz = (maxz + minz)/2.0;

radius = (maxx - minx)/2.0;
//Xu's .
/*
if (radius > (maxy - miny)/2.0) radius = (maxy - miny)/2.0;
if (radius > (maxz - minz)/2.0) radius = (maxz - minz)/2.0;
*/

// Ming's.
//if (radius < (maxy - miny)/2.0) radius = (maxy - miny)/2.0;
//if (radius < (maxz - minz)/2.0) radius = (maxz - minz)/2.0;


//radius = radius - 2*(dx + dy + dz);
radius = radius - (dx + dy + dz);
//radius = radius - 14*(dx + dy + dz); // Test
//radius = radius - dx;
printf("radius = %f\n", radius);

// compute function value at grid points
l = 0;
for (i = 0; i < nx; i++) {
   x = minx + i*dx;
   for (j = 0; j < ny; j++) {
      y = miny + j*dy;
      for (k = 0; k < nz; k++) {
         z = minz + k*dz;
         dis = sqrt((x-cx)*(x-cx) + (y-cy)*(y-cy)+(z-cz)*(z-cz));
         dis = dis - radius;

         // cut 
         if (dis < -Height) dis = -Height;
         if (dis > Height) dis = Height;

         v[l] = dis;
         //v[l] = dis;
         //printf("volume[%d, %d, %d] = %f\n", i,j,k , v[l]);
         l = l + 1;
      }
   }
}

}

/*-----------------------------------------------------------------------------*/
void  InitialData_Box(float *v, float minx, float maxx, float miny, float maxy,
                  float minz, float maxz, int nx, int ny, int nz)
//              float   *v,    /* the volume data of initial signed distance function */
//              float   minx   /* left end-point of x diection                        */
//              float   maxx   /* right end-point of x diection                       */
//              float   miny   /* left end-point of y diection                        */
//              float   maxy   /* right end-point of y diection                       */
//              float   minz   /* left end-point of z diection                        */
//              float   maxz   /* right end-point of z diection                       */
//              int     nx,    /* number of points in x direction                     */
//              int     ny,    /* number of points in y direction                     */
//              int     nz,    /* number of points in z direction                     */
{
int    i, j, k, l;
float  dx, dy, dz, x, y, z;
float  cx, cy, cz, radius, dis;

dx = (maxx - minx)/(nx - 1.0);
dy = (maxy - miny)/(ny - 1.0);
dz = (maxz - minz)/(nz - 1.0);

//determine the center and radius
cx = (maxx + minx)/2.0;
cy = (maxy + miny)/2.0;
cz = (maxz + minz)/2.0;

radius = (maxx - minx)/2.0;

if (radius > (maxy - miny)/2.0) radius = (maxy - miny)/2.0;
if (radius > (maxz - minz)/2.0) radius = (maxz - minz)/2.0;

radius = radius - 2*(dx + dy + dz);
//radius = radius - dx;
printf("radius = %f\n", radius);

// compute function value at grid points
l = 0;
for (i = 0; i < nx; i++) {
   x = fabs(minx + i*dx - cx);
   for (j = 0; j < ny; j++) {
      y = fabs(miny + j*dy - cy);
      for (k = 0; k < nz; k++) {
         z = fabs(minz + k*dz - cz);

         if (x >= radius && y >= radius && z >= radius) { 
             dis = sqrt((x-radius)*(x-radius) + (y-radius)*(y-radius)+(z-radius)*(z-radius));
         } else {

            dis = 0.0;
            if (x >= radius) dis = (x - radius)*(x - radius);
            if (y >= radius) dis = (y - radius)*(y - radius) + dis;
            if (z >= radius) dis = (z - radius)*(z - radius) + dis;
            dis = sqrt(dis);

            if (x >= radius && y <= radius && z <= radius) dis = x - radius; 
            if (x <= radius && y >= radius && z <= radius) dis = y - radius; 
            if (x <= radius && y <= radius && z >= radius) dis = z - radius; 

            if (x <= radius && y <= radius && z <= radius) {
               dis = radius - x;
               if (dis > radius - y) dis = radius - y; 
               if (dis > radius - z) dis = radius - z; 
               dis = -dis;
            }
         }

         // cut
         if (dis < -Height) dis = -Height;
         if (dis > Height) dis = Height;

         v[l] = dis;
         //v[l] = dis;
         //printf("volume[%d, %d, %d] = %f\n", i,j,k , v[l]);
         l = l + 1;
      }
   }
}

}


/*-----------------------------------------------------------------------------*/
void  InitialData_IrBox(float *v, float minx, float maxx, float miny, float maxy,
                  float minz, float maxz, int nx, int ny, int nz)
//              float   *v,    /* the volume data of initial signed distance function */
//              float   minx   /* left end-point of x diection                        */
//              float   maxx   /* right end-point of x diection                       */
//              float   miny   /* left end-point of y diection                        */
//              float   maxy   /* right end-point of y diection                       */
//              float   minz   /* left end-point of z diection                        */
//              float   maxz   /* right end-point of z diection                       */
//              int     nx,    /* number of points in x direction                     */
//              int     ny,    /* number of points in y direction                     */
//              int     nz,    /* number of points in z direction                     */
{
int    i, j, k, l;
float  dx, dy, dz, x, y, z;
float  cx, cy, cz, radius[3], dis;

dx = (maxx - minx)/(nx - 1.0);
dy = (maxy - miny)/(ny - 1.0);
dz = (maxz - minz)/(nz - 1.0);

//determine the center and radius
cx = (maxx + minx)/2.0;
cy = (maxy + miny)/2.0;
cz = (maxz + minz)/2.0;

radius[0] = (maxx - minx)/2.0;
radius[1] = (maxy - miny)/2.0;
radius[2] = (maxz - minz)/2.0;

radius[0] = radius[0] - 2*(dx + dy + dz);
radius[1] = radius[1] - 2*(dx + dy + dz);
radius[2] = radius[2] - 2*(dx + dy + dz);
//radius = radius - dx;
printf("radius[] = %f, %f, %f\n", radius[0], radius[1], radius[2]);

// compute function value at grid points
l = 0;
for (i = 0; i < nx; i++) {
   x = fabs(minx + i*dx - cx);
   for (j = 0; j < ny; j++) {
      y = fabs(miny + j*dy - cy);
      for (k = 0; k < nz; k++) {
         z = fabs(minz + k*dz - cz);

         if (x >= radius[0] && y >= radius[1] && z >= radius[2]) { 
             dis = sqrt((x-radius[0])*(x-radius[0]) 
                      + (y-radius[1])*(y-radius[1]) 
                      + (z-radius[2])*(z-radius[2]));
         } else {

            dis = 0.0;
            if (x >= radius[0]) dis = (x - radius[0])*(x - radius[0]);
            if (y >= radius[1]) dis = (y - radius[1])*(y - radius[1]) + dis;
            if (z >= radius[2]) dis = (z - radius[2])*(z - radius[2]) + dis;
            dis = sqrt(dis);

            if (x >= radius[0] && y <= radius[1] && z <= radius[2]) dis = x - radius[0]; 
            if (x <= radius[0] && y >= radius[1] && z <= radius[2]) dis = y - radius[1]; 
            if (x <= radius[0] && y <= radius[1] && z >= radius[2]) dis = z - radius[2]; 

            if (x <= radius[0] && y <= radius[1] && z <= radius[2]) {
               dis = radius[0] - x;
               if (dis > radius[1] - y) dis = radius[1] - y; 
               if (dis > radius[2] - z) dis = radius[2] - z; 
               dis = -dis;
            }
         }

         // cut
         if (dis < -Height) dis = -Height;
         if (dis > Height) dis = Height;

         v[l] = dis;
         //v[l] = dis;
         //printf("volume[%d, %d, %d] = %f\n", i,j,k , v[l]);
         l = l + 1;
      }
   }
}

}



/*-----------------------------------------------------------------------------*/
void    MeanCurvatureFlow(float *fun, float *fun_bak, float *coeff, float dt, 
                          float dx, float dy, float dz, int nx, int ny, int nz, int t)
//              float   *fun   /* the volume data of initial signed distance function */
//              float   *coeff /* spline coefficients                                 */
//              float   dt     /* time step size                                      */
//              float   dx     /* x direction spacing                                 */
//              float   dy     /* y direction spacing                                 */
//              float   dz     /* z direction spacing                                 */
//              int     nx,    /* number of points in x direction                     */
//              int     ny,    /* number of points in y direction                     */
//              int     nz,    /* number of points in z direction                     */
{
int   i, j, k, l, index;
float partials[10], grad1, grad2, H, w, ww, ww1, ww2;
float x, y, z, r, fvalue, maxf, maxt;
float radius, totalgrad, maxgrad;

radius = 5.0 - 6*dx;

ww1 =  Height - 3*Height0;
ww2 = (Height-Height0)*(Height-Height0)*(Height-Height0);
ww2 = 1.0/ww2;
maxf = -1000.0;

l = 0;
totalgrad = 0.0;
maxgrad = 0.0;
for (i = 0; i < nx; i++) {
   for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {

         index = (i*ny + j)*nz + k;
         fvalue = fun[index];
         fun_bak[index] = 0.0;

         //if (fabs(fvalue) < Height0) {
         if (fabs(fvalue) < Height) {
         //if (i == 10 & j == 25 & k == 25) {
         //EvaluateCubicSplineOrder2PartialsAtGridPoint(coeff, dx, dy, dz, nx, ny, nz, i, j, k, partials); 
         Divided_DifferenceOrder2PartialsAtGridPoint(fun, dx, dy, dz, nx, ny, nz, i, j, k, partials); 

         grad2 = partials[1]*partials[1] + partials[2]*partials[2] + partials[3]*partials[3];
         if (grad2 < 0.0001) grad2 = 0.0001;

         grad1 = sqrt(grad2);
if (fabs(fvalue) < 0.5*Height) {
   if (maxgrad < fabs(grad1 - 1.0)) maxgrad = fabs(grad1 - 1.0);
   totalgrad = totalgrad + grad1;
   //printf("grad = %f\n", grad1);
   l = l + 1;
}
         H = -partials[1]*(partials[4]*partials[1] + partials[5]*partials[2] + partials[6]*partials[3])
             -partials[2]*(partials[5]*partials[1] + partials[7]*partials[2] + partials[8]*partials[3])
             -partials[3]*(partials[6]*partials[1] + partials[8]*partials[2] + partials[9]*partials[3]);

/*
         H = -partials[1]*(partials[4] + partials[5] + partials[6]) 
             -partials[2]*(partials[5] + partials[7] + partials[8])
             -partials[3]*(partials[6] + partials[8] + partials[9]);
*/
         //H = (H/grad2 + partials[4] + partials[7] + partials[9])/grad1;
       H = H/grad2 + partials[4] + partials[7] + partials[9];

//	H= -(pow(partials[1],2)*(partials[7]+partials[9])-2*partials[1]*partials[2]*partials[5]
//	   +pow(partials[2],2)*(partials[4]+partials[9])-2*partials[1]*partials[3]*partials[6]
//	   +pow(partials[3],2)*(partials[4]+partials[7])-2*partials[2]*partials[3]*partials[8])/(grad2*grad1);

x = -5.0 + i*dx;
y = -5.0 + j*dy;
z = -5.0 + k*dz;
r = sqrt(x*x+ y*y + z*z);

         /*
         if (fabs(fvalue) < 0.01) { 
         printf("i, j, k = %d, %d, %d H = %f, Exact H = %f, %e, %f, %e\n", 
                                            i, j, k, H/2.0, 1.0/r, partials[0], grad1, r-radius-partials[0]);
         //printf("Dxyz = %f, %f, %f, %f,%f, %f, %f, %f,%f,%f\n", partials[0], partials[1],partials[2], partials[3],
         //                                  partials[4], partials[5],partials[6], partials[7],
         //                                  partials[8], partials[9]);
         }
         */

         w = H;
         ww = fabs(w); 
         //printf("ww = %f\n", ww);
         if (fabs(fvalue) <Height0) {
            if (ww > maxf) maxf = ww;
         }  else {
           //ww = (fabs(fvalue) - Height)*(fabs(fvalue) - Height)*(2*fabs(fvalue) + ww1)*ww2;
           //printf("ww = %f\n", ww);
           //w = ww*w;
         }
         fun_bak[index] = w;
         }
      }
   }
}

//printf("Total Grad = %f, Max Grad = %f\n", totalgrad/l, maxgrad);
maxt = dx/maxf;
if (maxt > dt) maxt = dt;
//printf("maxf = %f, dx__1 = %e, maxt = %e\n", maxf, dx, maxt);
//if (t == 1) maxt = 0.013401;

TotalTime = TotalTime + maxt;

for (i = 0; i < nx*ny*nz; i++) {
   fun[i] =  fun[i] + maxt*fun_bak[i];
   
   // cut the function by Height
   if (fun[i] < -Height) fun[i] = -Height;
   if (fun[i] > Height) fun[i] = Height;
}


}







 /*-------------------------------------------------------------------------------------*/

void   Constraint_MeanCurvatureFlow(float *func_h, float *fun, float *fun_bak, float *coeff, 
                          float weight, float dt, 
                          float dx, float dy, float dz, int nx, int ny, int nz, int t)
//              float   *fun   /* the volume data of initial signed distance function */
//              float   *coeff /* spline coefficients                                 */
//              float   dt     /* time step size                                      */
//              float   dx     /* x direction spacing                                 */
//              float   dy     /* y direction spacing                                 */
//              float   dz     /* z direction spacing                                 */
//              int     nx,    /* number of points in x direction                     */
//              int     ny,    /* number of points in y direction                     */
//              int     nz,    /* number of points in z direction                     */
{
int   i, j, k, index;
float partials[10], partials_h[10], grad1, grad2, H, w, ww, ww1, ww2;
float x, y, z, r, fvalue, maxf, maxt;
float radius, fx[7], fy[7], fz[7], Lphi, Hphi;
float a, b, c, d, e, f,temp;

radius = 5.0 - 6*dx;

ww1 =  Height - 3*Height0;
ww2 = (Height-Height0)*(Height-Height0)*(Height-Height0);
ww2 = 1.0/ww2;
maxf = -1000.0;

for (i = 0; i < nx; i++) {
   for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {

         index = (i*ny + j)*nz + k;
         fvalue = fun[index];
         fun_bak[index] = 0.0;
 
	//if (fabs(fvalue) < Height0) {
         if (fabs(fvalue) < Height) {
         //if (i == 10 & j == 25 & k == 25) {
         //EvaluateCubicSplineOrder2PartialsAtGridPoint(coeff, dx, dy, dz, nx, ny, nz, i, j, k, partials);
         Divided_DifferenceOrder2PartialsAtGridPoint(fun, dx, dy, dz, nx, ny, nz, i, j, k, partials); 

         grad2 = partials[1]*partials[1] + partials[2]*partials[2] + partials[3]*partials[3];
         if (grad2 < 0.0001) grad2 = 0.0001;

         grad1 = sqrt(grad2);

         H = -partials[1]*(partials[4]*partials[1] + partials[5]*partials[2] + partials[6]*partials[3])
             -partials[2]*(partials[5]*partials[1] + partials[7]*partials[2] + partials[8]*partials[3])
             -partials[3]*(partials[6]*partials[1] + partials[8]*partials[2] + partials[9]*partials[3]);

         H = H/grad2 + partials[4] + partials[7] + partials[9];
x = -5.0 + i*dx;
y = -5.0 + j*dy;
z = -5.0 + k*dz;
r = sqrt(x*x+ y*y + z*z);

         /*
         if (fabs(fvalue) < 0.01) { 
         printf("i, j, k = %d, %d, %d H = %f, Exact H = %f, %e, %f, %e\n", 
                                            i, j, k, H/2.0, 1.0/r, partials[0], grad1, r-radius-partials[0]);
         //printf("Dxyz = %f, %f, %f, %f,%f, %f, %f, %f,%f,%f\n", partials[0], partials[1],partials[2], partials[3],
         //                                  partials[4], partials[5],partials[6], partials[7],
         //                                  partials[8], partials[9]);
         }
         */

         Divided_DifferenceOrder2PartialsAtGridPoint(func_h, dx, dy, dz, nx, ny, nz, i, j, k, partials_h); 
         //EvaluateCubicSplineOrder2PartialsAtGridPoint(func_h, dx, dy, dz, nx, ny, nz, i, j, k, partials_h); 
         //w = (func_h[index] + weight)*H*grad1 
   //      w = (partials_h[0]*partials_h[0] + weight)*H*grad1 
   //            + 4*partials_h[0]*(partials_h[1]*partials[1] 
     //                           + partials_h[2]*partials[2] 
       //                         + partials_h[3]*partials[3]);


         // Eno scheme for partials
         Get_Seven_Function_Values_X(fun, nx, ny, nz, i, j, k, fx);
         Get_Seven_Function_Values_Y(fun, nx, ny, nz, i, j, k, fy);
         Get_Seven_Function_Values_Z(fun, nx, ny, nz, i, j, k, fz);

         //Cubic_Eno_Interpolation(fx, &b, &a, dx);
         Quadr_Eno_Interpolation(fx, &b, &a, dx);
         //Linear_Eno_Interpolation(fx, &b, &a, dx);

         //Cubic_Eno_Interpolation(fy, &d, &c, dy);
         Quadr_Eno_Interpolation(fy, &d, &c, dy);
         //Linear_Eno_Interpolation(fy, &d, &c, dy);

         //Cubic_Eno_Interpolation(fz, &f, &e, dz);
         Quadr_Eno_Interpolation(fz, &f, &e, dz);
         //Linear_Eno_Interpolation(fz, &f, &e, dz);
         //partials[1] = fx;
         //partials[2] = fy;
         //partials[3] = fz;
	
	float zero=0.0f;
	//EO scheme

	if(partials_h[1]>0) partials[1]=Max_Of_Two(a,zero)+Min_Of_Two(b,zero);
	else partials[1]=Min_Of_Two(a,zero)+Max_Of_Two(b,zero);

	if(partials_h[2]>0) partials[2]=Max_Of_Two(c,zero)+Min_Of_Two(d,zero);
	else partials[2]=Min_Of_Two(c,zero)+Max_Of_Two(d,zero);

	if(partials_h[3]>0) partials[3]=Max_Of_Two(e,zero)+Min_Of_Two(f,zero);
	else partials[3]=Min_Of_Two(e,zero)+Max_Of_Two(f,zero);

//         Lphi = (partials_h[0] + weight)*H; 


         Lphi = (partials_h[0] + weight)*H; 
         Hphi =  2*(partials_h[1]*partials[1] 
                  + partials_h[2]*partials[2] 
                  + partials_h[3]*partials[3]);
         //printf("Lphi, Hphi = %f, %f\n",Lphi, Hphi);
         w = Lphi + Hphi;
         ww = fabs(w); 
         //printf("ww = %f\n", ww);
         if (fabs(fvalue) <Height0) {
            if (ww > maxf) {
               maxf = ww;
               //printf("ij,k = %d, %d, %d, ww = %f, grad1 = %f, H = %f\n", i,j,k,
               //                         partials_h[1],partials_h[2],partials_h[3]);
            }
         }  else {
           //ww = (fabs(fvalue) - Height)*(fabs(fvalue) - Height)*(2*fabs(fvalue) + ww1)*ww2;
           //printf("ww = %f\n", ww);
           //w = ww*w;
         }
         fun_bak[index] = w;

         }
      }
   }
}

maxt = dx/maxf;
if (maxt > dt) maxt = dt;
 if (maxt<0.00001) maxt =0.00001;
//printf("maxf = %f, dx_2 = %f, maxt = %f\n", maxf, dx, maxt);

TotalTime = TotalTime + maxt;

// update fun
for (i = 0; i < nx*ny*nz; i++) {

   fun[i] =  fun[i] + maxt*fun_bak[i];
   
   // cut the function by Height
   if (fun[i] < -Height) fun[i] = -Height;
   if (fun[i] > Height) fun[i] = Height;
}

}

void ComputeTensorXYZ()
{
float f[3], fx[3], fxx[3];

f[0] = 1.0/6.0;
f[1] = 2.0/3.0;
f[2] = f[0];

fx[0] = -0.5;
fx[1] = 0.0;
fx[2] = 0.5;

fxx[0] = 1.0;
fxx[1] = -2.0;
fxx[2] = 1.0;

Tensor_333(f, f, f, TensorF);

Tensor_333(fx, f, f, TensorFx);
Tensor_333(f, fx, f, TensorFy);
Tensor_333(f, f, fx, TensorFz);

Tensor_333(fxx,f, f, TensorFxx);
Tensor_333(fx, fx,f, TensorFxy);
Tensor_333(fx, f, fx,TensorFxz);

Tensor_333(f, fxx,f,  TensorFyy);
Tensor_333(f, fx, fx, TensorFyz);
Tensor_333(f, f,  fxx,TensorFzz);
}


/*-----------------------------------------------------------------------------*/
void    ReInitilazation(float *fun, float *coeff, float dx, float dy, float dz,
                          int nx, int ny, int nz)
//              float   *fun   /* the volume data of initial signed distance function */
//              float   *coeff /* spline coefficients                                 */
//              float   dx     /* x direction spacing                                 */
//              float   dy     /* y direction spacing                                 */
//              float   dz     /* z direction spacing                                 */
//              int     nx,    /* number of points in x direction                     */
//              int     ny,    /* number of points in y direction                     */
//              int     nz,    /* number of points in z direction                     */
{
int   i, j, k, index, t;
float partials[4], grad1, grad2;
float x, y, z, r;
float dt, sd;


dt = 0.01;

for (t = 0; t < 1000; t++) {

printf("Initialization t = %d\n", t);
for (i = 0; i < nx; i++) {
   for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
         EvaluateCubicSplineOrder1PartialsAtGridPoint(coeff, dx, dy, dz, nx, ny, nz, i, j, k, partials);
         //printf("i , j, k = %d, %d, %d Dxyz = %f, %f, %f, %f\n", i, j, k,partials[0], partials[1],partials[2], partials[3]);

         grad2 = partials[1]*partials[1] + partials[2]*partials[2] + partials[3]*partials[3];
         grad1 = sqrt(grad2);
         sd = partials[0]/sqrt(partials[0]*partials[0]+grad2*dx*dx);

x = -5.0 + i*dx;
y = -5.0 + j*dy;
z = -5.0 + k*dz;
r = sqrt(x*x+ y*y + z*z);

         //if (r < 4.0)
         if (i == 44 && j== 31 && k ==  26)
         printf("t = %d, i , j, k = %d, %d, %d, grad1 =  %f  S(d) = %f\n", t, i, j, k, grad1, sd);
         index = (i*ny + j)*nz + k;
         fun[index] = fun[index] - dt*sd * (grad1 - 1.0);
      }
   }
}

for (i = 0; i < nx*ny*nz; i++) {
   coeff[i] = fun[i];
}
ConvertToInterpolationCoefficients_3D(coeff, nx,ny,nz,CVC_DBL_EPSILON);

}
}

/*-----------------------------------------------------------------------------*/
void    ReInitilazation_Upwind_Eno_Engquist(float *fun, float *fun_bak, float *coeff, float dx, float dy, float dz,
                          int nx, int ny, int nz)
//              float   *fun   /* the volume data of initial signed distance function */
//              float   *coeff /* spline coefficients                                 */
//              float   dx     /* x direction spacing                                 */
//              float   dy     /* y direction spacing                                 */
//              float   dz     /* z direction spacing                                 */
//              int     nx,    /* number of points in x direction                     */
//              int     ny,    /* number of points in y direction                     */
//              int     nz,    /* number of points in z direction                     */
{
int   i, j, k, index, t;
float fvalue, grad1, grad2, ww, ww1;
float x, y, z, r, fx[7], fy[7], fz[7];
float dt, sd, sdp, sdm;
float a,  b,  c,  d,  e,  f;
float ap, bp, cp, dp, ep, fp;
float am, bm, cm, dm, em, fm;


dt = 0.5*dx;
dt = 0.05;

for (t = 0; t < 100; t++) {
printf("Initialization t = %d\n", t);
for (i = 0; i < nx; i++) {
   for (j = 0; j < ny; j++) {
//for (i = 10; i < 11; i++) {
//   for (j = 10; j < 11; j++) {
      for (k = 0; k < nz; k++) {

         Get_Seven_Function_Values_X(fun, nx, ny, nz, i, j, k, fx);
         //Cubic_Eno_Interpolation(fx, &b, &a, dx);
         Quadr_Eno_Interpolation(fx, &b, &a, dx);
         //Linear_Eno_Interpolation(fx, &b, &a, dx);

         ap = a; 
         if (a < 0.0) ap = 0.0;
         bp = b; 
         if (b < 0.0) bp = 0.0;

         am = a;   
         if (a > 0.0) am = 0.0;
         bm = b;
         if (b > 0.0) bm = 0.0;

         Get_Seven_Function_Values_Y(fun, nx, ny, nz, i, j, k, fy);
         //Cubic_Eno_Interpolation(fy, &d, &c, dy);
         Quadr_Eno_Interpolation(fy, &d, &c, dy);
         //Linear_Eno_Interpolation(fy, &d, &c, dy);


         cp = c;   
         if (c < 0.0) cp = 0.0;
         dp = d;
         if (d < 0.0) dp = 0.0;

         cm = c;            
         if (c > 0.0) cm = 0.0;
         dm = d;
         if (d > 0.0) dm = 0.0;


         Get_Seven_Function_Values_Z(fun, nx, ny, nz, i, j, k, fz);
         //printf("%f,%f,%f,%f,%f,%f,%f\n", fz[0]-fz[1],fz[1]-fz[2],fz[2]-fz[3],fz[3]-fz[4],fz[4]-fz[5],fz[5]-fz[6],dx); 
         //Cubic_Eno_Interpolation(fz, &f, &e, dz);
         Quadr_Eno_Interpolation(fz, &f, &e, dz);
         //Linear_Eno_Interpolation(fz, &f, &e, dz);


         ep = e;   
         if (e < 0.0) ep = 0.0;
         fp = f;
         if (f < 0.0) fp = 0.0;

         em = e;            
         if (e > 0.0) em = 0.0;
         fm = f;
         if (f > 0.0) fm = 0.0;

//printf("i,j,k = %d, %d, %d a--f = %f, %f, %f, %f, %f, %f, grad = %f, %f \n", i, j, k, a, b, c, d, e, f,
//                           a*a + c*c + e*e, b*b + d*d + f*f);

         grad2 = 0.5*(a*a + c*c + e*e + b*b + d*d + f*f);
         grad1 = sqrt(grad2);
         index = (i*ny + j)*nz + k;

         fvalue = fun[index];
         sd = fvalue/sqrt(fvalue*fvalue + grad2*dx*dx);
         sdp = sd;
         if (sd < 0.0) sdp = 0.0;

         sdm = sd;
         if (sd > 0.0) sdm = 0.0;


x = -5.0 + i*dx;
y = -5.0 + j*dy;
z = -5.0 + k*dz;
r = sqrt(x*x+ y*y + z*z);

         //if (r < 0.4)
         //if (r  > 1.0)
         //if (i == 44 && j== 31 && k ==  26)
         //printf("t = %d, i , j, k = %d, %d, %d, grad1 =  %f  S(d) = %f\n", t, i, j, k, grad1, sd);
         if (sdp > 0.0) {
            ap = ap*ap;
            cp = cp*cp;
            ep = ep*ep;
            bm = bm*bm;
            dm = dm*dm;
            fm = fm*fm;

          
            // Engquist-Osher scheme
            ww = sqrt(ap + cp + ep + bm + dm + fm);
            // end of Engquist-Osher scheme
          
/*
            // Godunov scheme
            ww = ap;
            if (ww < bm) ww = bm;
            ww1 = cp;
            if (ww1 < dm) ww1 = dm;
            ww = ww + ww1;
            ww1 = ep;
            if (ww1 < fm) ww1 = fm;
            ww = ww + ww1;
            ww = sqrt(ww);
            // end Godunov scheme            
*/
 
            fun_bak[index] = fvalue - dt*sdp * (ww - 1.0);
         }
         if (sdm < 0.0) {
            am = am*am; 
            cm = cm*cm; 
            em = em*em;
            bp = bp*bp;
            dp = dp*dp;
            fp = fp*fp;

            
            // Engquist-Osher scheme
             ww = sqrt(am + cm + em + bp + dp + fp);
            // end Engquist-Osher scheme
            
/*
            // Godunov scheme
            ww = am;
            if (ww < bp) ww = bp;
            ww1 = cm;
            if (ww1 < dp) ww1 = dp;
            ww = ww + ww1;
            ww1 = em;
            if (ww1 < fp) ww1 = fp;
            ww = ww + ww1;
            ww = sqrt(ww);
            // end Godunov scheme            
*/

            fun_bak[index] = fvalue - dt*sdm * (ww - 1.0);
         }

         //if (fabs(fvalue)  < 0.3)
         //if (r  < 4.01 && r > 3.9)
         if (r  < 0.4)
         //if (i == nx/2 && j == ny/2 && k == nz/2)
         //if (i == 0 && j == 0 && k == 0)
         printf("t = %d, i , j, k = %d, %d, %d, grad1 =  %f  S(d) = %f, fun = %f,  ww = %f\n", t, i, j, k, grad1, sd, fvalue, ww);

      }
   }
}

for (i = 0; i < nx*ny*nz; i++) {
   fun[i] = fun_bak[i];
}

}

for (i = 0; i < nx*ny*nz; i++) {
   coeff[i] = fun[i];
}
ConvertToInterpolationCoefficients_3D(coeff, nx,ny,nz,CVC_DBL_EPSILON);
}


/*-----------------------------------------------------------------------------*/
void    ReInitilazation_Upwind_Eno_Godunov(float *fun, float *fun_bak, float *coeff, float dx, float dy, float dz,
                          int nx, int ny, int nz, int tt)
//              float   *fun   /* the volume data of initial signed distance function */
//              float   *coeff /* spline coefficients                                 */
//              float   dx     /* x direction spacing                                 */
//              float   dy     /* y direction spacing                                 */
//              float   dz     /* z direction spacing                                 */
//              int     nx,    /* number of points in x direction                     */
//              int     ny,    /* number of points in y direction                     */
//              int     nz,    /* number of points in z direction                     */
{
int   i, j, k, index, t;
float fvalue, fvalue_new, fvalue_old, grad1, grad2, ww, ww1;
float x, y, z, r, fx[7], fy[7], fz[7];
float dt, sd, sdp, sdm, maxt, maxf;
float a,  b,  c,  d,  e,  f, width;
int temp;

dt = 0.5*dx;
//dt = 0.05;
width = dx;

// coeff keep the initial value of \phi
for (i = 0; i < nx*ny*nz; i++) {
   coeff[i] = fun[i];
}


//for (t = 0; t < 10; t++) {
for (t = 0; t < 20; t++) {

maxf = -100000.0;
//printf("Initialization t = %d\n", t);
for (i = 0; i < nx; i++) {
   for (j = 0; j < ny; j++) {
//for (i = 10; i < 11; i++) {
//   for (j = 10; j < 11; j++) {
      for (k = 0; k < nz; k++) {

         index = (i*ny + j)*nz + k;
         fvalue = fun[index];
         fvalue_old = coeff[index];
         fun_bak[index] = 0.0;

         Get_Seven_Function_Values_X(fun, nx, ny, nz, i, j, k, fx);
         Get_Seven_Function_Values_Y(fun, nx, ny, nz, i, j, k, fy);
         Get_Seven_Function_Values_Z(fun, nx, ny, nz, i, j, k, fz);

         if (fabs(fx[2]) < Height || fabs(fx[4]) < Height ||
             fabs(fy[2]) < Height || fabs(fy[4]) < Height ||
             fabs(fz[2]) < Height || fabs(fz[4]) < Height ||
             fabs(fvalue) < Height) {

         //Cubic_Eno_Interpolation(fx, &b, &a, dx);
         Quadr_Eno_Interpolation(fx, &b, &a, dx);
         //Linear_Eno_Interpolation(fx, &b, &a, dx);

         //Cubic_Eno_Interpolation(fy, &d, &c, dy);
         Quadr_Eno_Interpolation(fy, &d, &c, dy);
         //Linear_Eno_Interpolation(fy, &d, &c, dy);

         //Cubic_Eno_Interpolation(fz, &f, &e, dz);
         Quadr_Eno_Interpolation(fz, &f, &e, dz);
         //Linear_Eno_Interpolation(fz, &f, &e, dz);

/*
if (i == nx/2 && j == ny/2 && k == nz/2)
printf("i,j,k = %d, %d, %d a--f = %f, %f, %f, %f, %f, %f, grad = %f, %f \n", i, j, k, a, b, c, d, e, f,
                           a*a + c*c + e*e, b*b + d*d + f*f);
*/

         //grad2 =  Gradient_2(fun, dx, dy, dz, i, j, k, nx, ny, nz);

         grad2 = 0.5*(a*a + c*c + e*e + b*b + d*d + f*f);
         grad1 = sqrt(grad2);

         sd = 1.0; 

         // Using updated \phi values
         
         fvalue = fvalue_old; // Using old phi values
         //sd = fvalue/Height; // Xu Test The result is bad
         //sd =  fvalue/sqrt(fvalue*fvalue + grad2*dx*dx); // Zhao Hong Kai's choice is bad
         //sd =  fvalue*(3*Height*Height- fvalue*fvalue)/(2*Height*Height*Height);  // The result is bad
         if (fvalue < 0.0) sd = -1.0;   // Xu's choice
         if (fabs(fvalue) < width) sd = fvalue*(3*width*width - fvalue*fvalue)/(2*width*width*width); // Xu's

         // combintation
         //sd =  0.5*(sd + fvalue/sqrt(fvalue*fvalue + grad2*dx*dx)); // Xu_zhao Result is good 
         //sd =  0.5*(sd + fvalue/Height);  // The result is good Xu_Xuf
         // Hermite interpolation
         //sd =  0.5*(sd + fvalue*(3*Height*Height- fvalue*fvalue)/(2*Height*Height*Height));  
         // Xu_Hermite The result is bad
         

         sdp = sd;
         if (sd < 0.0) sdp = 0.0;

         sdm = sd;
         if (sd > 0.0) sdm = 0.0;


x = -5.0 + i*dx;
y = -5.0 + j*dy;
z = -5.0 + k*dz;
r = sqrt(x*x+ y*y + z*z);

         //if (r < 0.4)
         //if (r  > 1.0)
         //if (i == 44 && j== 31 && k ==  26)
         //printf("t = %d, i , j, k = %d, %d, %d, grad1 =  %f  S(d) = %f\n", t, i, j, k, grad1, sd);

         if (sdp >= 0.0) {
 
            ww = sqrt(Extreme_Positive(a, b)+ Extreme_Positive(c,d)+ Extreme_Positive(e,f));
            fvalue_new = - dt*sdp * (ww - 1.0);   // Zhao's method
         }
         if (sdm < 0.0) {

            ww = sqrt(Extreme_Negative(a, b)+ Extreme_Negative(c,d)+ Extreme_Negative(e,f));
            fvalue_new =  - dt*sdm * (ww - 1.0);   // Zhao's method
         }

         fun_bak[index] = fvalue_new;
         if (fabs(fun_bak[index]) > maxf) maxf = fabs(fun_bak[index]);

         //fun_bak[index] = fvalue - dt*fvalue * (grad1 - 1.0)/grad1; // Xu's method

         //if (fabs(fvalue)  < 0.001*Height)
         //if (r  < 4.01 && r > 3.9)
         //if (r  < 0.4)
         //if (i == nx/2 && j == ny/2 && k == nz/2)
         //if (i == 0 && j == 0 && k == 0)
         //if (tt == 9 && t == 99)
         //printf("t = %d, i , j, k = %d, %d, %d, grad1 =  %f  S(d) = %f, fun = %f,  ww = %f, r=%f\n", 
         //        t, i, j, k, grad1, sd, fvalue, ww, r);

         } 
      }
   }
}

maxt = dx/maxf;
//if (maxt > t) maxt = t;    // Tem chage for Test
//printf("maxf = %f, dx = %f, maxt = %f\n", maxf, dx, maxt);


// Update function 
for (i = 0; i < nx*ny*nz; i++) {
  x = fun[i] + maxt*fun_bak[i];
   //if (fun[i]*x > 0.0 ){  // Added laterly Bu Guang Hua 
   //if (fun[i]*x > 0.0 ){  // Added laterly 
      fun[i] = x;
   //}  else {
   //   fun[i] = 0.0;
  // }

   // cut the function by Height
   if (fun[i] < -Height) fun[i] = -Height;
   if (fun[i] > Height) fun[i] = Height;
}

}


for (i = 0; i < nx*ny*nz; i++) {
   coeff[i] = fun[i];
}
ConvertToInterpolationCoefficients_3D(coeff, nx,ny,nz,CVC_DBL_EPSILON);
}

/*-----------------------------------------------------------------------------*/
void    ReInitilazation_Upwind_Eno_Godunov_Xu(float *fun, float *fun_bak, float *coeff, float dx, float dy, float dz,
                          int nx, int ny, int nz, int tt)
//              float   *fun   /* the volume data of initial signed distance function */
//              float   *coeff /* spline coefficients                                 */
//              float   dx     /* x direction spacing                                 */
//              float   dy     /* y direction spacing                                 */
//              float   dz     /* z direction spacing                                 */
//              int     nx,    /* number of points in x direction                     */
//              int     ny,    /* number of points in y direction                     */
//              int     nz,    /* number of points in z direction                     */
{
int   i, j, k, index, t;
float fvalue, fvalue_new, grad1, grad2, ww, ww1;
float x, y, z, r, fx[7], fy[7], fz[7];
float dt, sd, sdp, sdm, maxt, maxf;
float a,  b,  c,  d,  e,  f;
                                                                                                                                                             
                                                                                                                                                             
dt = 0.5*dx;
//dt = 0.05;
                                                                                                                                                             
//for (t = 0; t < 10; t++) {
for (t = 0; t < 20; t++) {
                                                                                                                                                             
maxf = -100000.0;
printf("Initialization t = %d\n", t);
for (i = 0; i < nx; i++) {
   for (j = 0; j < ny; j++) {
//for (i = 10; i < 11; i++) {
//   for (j = 10; j < 11; j++) {
      for (k = 0; k < nz; k++) {
                                                                                                                                                             
         index = (i*ny + j)*nz + k;
         fvalue = fun[index];
         fun_bak[index] = 0.0;
                                                                                                                                                             
         Get_Seven_Function_Values_X(fun, nx, ny, nz, i, j, k, fx);
         Get_Seven_Function_Values_Y(fun, nx, ny, nz, i, j, k, fy);
         Get_Seven_Function_Values_Z(fun, nx, ny, nz, i, j, k, fz);
   
                                                                                                                                                             
         if (fabs(fx[2]) < Height || fabs(fx[4]) < Height ||
             fabs(fy[2]) < Height || fabs(fy[4]) < Height ||
             fabs(fz[2]) < Height || fabs(fz[4]) < Height ||
             fabs(fvalue) < Height) {
                                                                                                                                                             
         //Cubic_Eno_Interpolation(fx, &b, &a, dx);
         Quadr_Eno_Interpolation(fx, &b, &a, dx);
         //Linear_Eno_Interpolation(fx, &b, &a, dx);
                                                                                                                                                             
         //Cubic_Eno_Interpolation(fy, &d, &c, dy);
         Quadr_Eno_Interpolation(fy, &d, &c, dy);
         //Linear_Eno_Interpolation(fy, &d, &c, dy);
                                                                                                                                                             
         //Cubic_Eno_Interpolation(fz, &f, &e, dz);
         Quadr_Eno_Interpolation(fz, &f, &e, dz);
         //Linear_Eno_Interpolation(fz, &f, &e, dz);
                                                                                                                                                             
/*
if (i == nx/2 && j == ny/2 && k == nz/2)
printf("i,j,k = %d, %d, %d a--f = %f, %f, %f, %f, %f, %f, grad = %f, %f \n", i, j, k, a, b, c, d, e, f,
                           a*a + c*c + e*e, b*b + d*d + f*f);
*/
                                                                                                                                                             
         //grad2 =  Gradient_2(fun, dx, dy, dz, i, j, k, nx, ny, nz);
                                                                                                                                                             
         grad2 = 0.5*(a*a + c*c + e*e + b*b + d*d + f*f);
         grad1 = sqrt(grad2);
                                                                                                                                                             
         sd = fvalue/sqrt(fvalue*fvalue + grad2*dx*dx);
/*
sd = 1.0;
if (fvalue < 0.0)
sd = -1.0;
*/
         sdp = sd;
         if (sd < 0.0) sdp = 0.0;
                                                                                                                                                             
         sdm = sd;
         if (sd > 0.0) sdm = 0.0;
                                                                                                                                                             
                                                                                                                                                             
x = -5.0 + i*dx;
y = -5.0 + j*dy;
z = -5.0 + k*dz;
r = sqrt(x*x+ y*y + z*z);
                                                                                                                                                          
         //if (r < 0.4)
         //if (r  > 1.0)
         //if (i == 44 && j== 31 && k ==  26)
         //printf("t = %d, i , j, k = %d, %d, %d, grad1 =  %f  S(d) = %f\n", t, i, j, k, grad1, sd);
                                                                                                                                                             
         if (sdp >= 0.0) {
            ww = Extreme_Positive(a, b)+ Extreme_Positive(c,d)+ Extreme_Positive(e,f);
            fvalue_new = - 0.7*dt*sdp * (ww - 1.0);   // Zhao's method
         }
         if (sdm < 0.0) {
            ww = Extreme_Negative(a, b)+ Extreme_Negative(c,d)+ Extreme_Negative(e,f);
            fvalue_new =  - 0.7*dt*sdm * (ww - 1.0);   // Zhao's method
         }
                                                                                                                                                             
         fun_bak[index] = fvalue_new;
         if (fabs(fun_bak[index]) > maxf) maxf = fabs(fun_bak[index]);
                                                                                                                                                             
         //fun_bak[index] = fvalue - dt*fvalue * (grad1 - 1.0)/grad1; // Xu's method
                                                                                                                                                             
         //if (fabs(fvalue)  < 0.001*Height)
         //if (r  < 4.01 && r > 3.9)
         //if (r  < 0.4)
         //if (i == nx/2 && j == ny/2 && k == nz/2)
         //if (i == 0 && j == 0 && k == 0)
         //if (tt == 9 && t == 99)
         //printf("t = %d, i , j, k = %d, %d, %d, grad1 =  %f  S(d) = %f, fun = %f,  ww = %f, r=%f\n",
         //        t, i, j, k, grad1, sd, fvalue, ww, r);
                                                                                                                                                             
         }
      }
   }
}
                                                                                                                                                             
                                                                                                                                                             
maxt = dx/maxf;
if (maxt > t) maxt = t;
//printf("maxf = %f, dx = %f, maxt = %f\n", maxf, dx, maxt);
                                                                                                                                                             
                                                                                                                                                             
for (i = 0; i < nx*ny*nz; i++) {
   fun[i] = fun[i] + maxt*fun_bak[i];
                                                                                                                                                             
   // cut the function by Height
   if (fun[i] < -Height) fun[i] = -Height;
   if (fun[i] > Height) fun[i] = Height;
}
                                                                                                                                                             
}
                                                                                                                                                             
for (i = 0; i < nx*ny*nz; i++) {
   coeff[i] = fun[i];
}
ConvertToInterpolationCoefficients_3D(coeff, nx,ny,nz,CVC_DBL_EPSILON);
}


/*-----------------------------------------------------------------------------*/
void    DiviededDifferencing_3j(float *f0, float *f1, float *f2, float *f3, float dx)
{
int i; 
float indx1, indx2, indx3;


indx1 = 1.0/dx;

// The first order divided difference
for (i = 0; i < 6; i++) {
   f1[i] = (f0[i+1] - f0[i])*indx1;
}

// The second order divided difference
indx2 = indx1/2.0;
for (i = 0; i < 5; i++) {
   f2[i] = (f1[i+1] - f1[i])*indx2;
}

// The third order divided difference
indx3 = indx1/3.0;
for (i = 0; i < 4; i++) {
   f3[i] = (f2[i+1] - f2[i])*indx3;
}

}

/*-----------------------------------------------------------------------------*/
void    DiviededDifferencing_2j(float *f0, float *f1, float *f2, float dx)
{
int i;
float indx1, indx2, indx3;


indx1 = 1.0/dx;

// The first order divided difference
for (i = 1; i < 5; i++) {
   f1[i] = (f0[i+1] - f0[i])*indx1;
}

// The second order divided difference
indx2 = indx1/2.0;
for (i = 1; i < 4; i++) {
   f2[i] = (f1[i+1] - f1[i])*indx2;
}

}


/*-----------------------------------------------------------------------------*/
void    DiviededDifferencing_1j(float *f0, float *f1, float dx)
{
int i;
float indx1, indx2, indx3;


indx1 = 1.0/dx;

// The first order divided difference
for (i = 2; i < 4; i++) {
   f1[i] = (f0[i+1] - f0[i])*indx1;
}

}



/*-----------------------------------------------------------------------------*/
void    Cubic_Eno_Interpolation(float *f0, float *u_plus, float *u_minus, float dx)
{
int    k;
float f1[6], f2[5], f3[4];
float x1, x2, x3;

DiviededDifferencing_3j(f0, f1, f2, f3, dx);
//printf("%f,%f,%f,%f,%f,%f,%f\n", f0[0],f0[1],f0[2],f0[3],f0[4],f0[5],f0[6]);

// Compute u^-
k = 2;
if (fabs(f2[k-1]) < fabs(f2[k])) k = k - 1;
if (fabs(f3[k-1]) < fabs(f3[k])) k = k - 1;

x1 = (3 - k)*dx;
x2 = (2 - k)*dx;
x3 = (1 - k)*dx;

*u_minus = f1[k] + f2[k]*(x1 + x2) + f3[k]*(x1*x2 + x1*x3 + x2*x3);
//printf("k- = %d, u_minus = %f\n", k, *u_minus);

// Compute u^+
k = 3;
if (fabs(f2[k-1]) < fabs(f2[k])) k = k - 1;
if (fabs(f3[k-1]) < fabs(f3[k])) k = k - 1;

x1 = (3 - k)*dx;
x2 = (2 - k)*dx;
x3 = (1 - k)*dx;

*u_plus = f1[k] + f2[k]*(x1 + x2) + f3[k]*(x1*x2 + x1*x3 + x2*x3);
//printf("k+ = %d, *u_plus = %f\n", k, *u_plus);
}

/*-----------------------------------------------------------------------------*/
void    Quadr_Eno_Interpolation(float *f0, float *u_plus, float *u_minus, float dx)
{
int    k;
float f1[6], f2[5];
float x1, x2, x3;

DiviededDifferencing_2j(f0, f1, f2, dx);
//printf("%f,%f,%f,%f,%f,%f,%f\n", f0[0],f0[1],f0[2],f0[3],f0[4],f0[5],f0[6]);

// Compute u^-
k = 2;
if (fabs(f2[k-1]) < fabs(f2[k])) k = k - 1;

x1 = (3 - k)*dx;
x2 = (2 - k)*dx;

*u_minus = f1[k] + f2[k]*(x1 + x2);
//printf("k- = %d, u_minus = %f\n", k, *u_minus);

// Compute u^+
k = 3;
if (fabs(f2[k-1]) < fabs(f2[k])) k = k - 1;

x1 = (3 - k)*dx;
x2 = (2 - k)*dx;

*u_plus = f1[k] + f2[k]*(x1 + x2);
//printf("k+ = %d, *u_plus = %f\n", k, *u_plus);
}

/*-----------------------------------------------------------------------------*/
void    Linear_Eno_Interpolation(float *f0, float *u_plus, float *u_minus, float dx)
{
int    k;
float f1[6];
float x1, x2, x3;

DiviededDifferencing_1j(f0, f1, dx);
//printf("%f,%f,%f,%f,%f,%f,%f\n", f0[0],f0[1],f0[2],f0[3],f0[4],f0[5],f0[6]);

// Compute u^-
k = 2;
*u_minus = f1[k];

//printf("k- = %d, u_minus = %f\n", k, *u_minus);

// Compute u^+
k = 3;
*u_plus = f1[k];
//printf("k+ = %d, *u_plus = %f\n", k, *u_plus);
}

/*-----------------------------------------------------------------------------*/
void    Get_Seven_Function_Values_X(float *coeff, int nx, int ny, int nz, int u, int v, int w, float *fx)
{
int   j, l, uj;

if (u >= 3 && u <= nx-4) {  
   for (j = 0; j < 7; j++) {
      l = u - 3 + j; 
      fx[j] = coeff[(l*ny + v)*nz + w];
   }
   return;
}

if (u == 2) {
   for (j = 1; j < 7; j++) {
      l = j -1;
      fx[j] = coeff[(l*ny + v)*nz + w];
   }
   fx[0] = 2*fx[1] - fx[2];
   return;
}

if (u == 1) {
   for (j = 2; j < 7; j++) {
      l = j-2;
      fx[j] = coeff[(l*ny + v)*nz + w];
   }

   fx[1] = 2*fx[2] - fx[3];
   fx[0] = 2*fx[1] - fx[2];
   return;
}


if (u == 0) {
   for (j = 3; j < 7; j++) {
      l = j -3;
      fx[j] = coeff[(l*ny + v)*nz + w];
   }

   fx[2] = 2*fx[3] - fx[4];
   fx[1] = 2*fx[2] - fx[3];
   fx[0] = 2*fx[1] - fx[2];
   return;
}

/*
if (u <= 2) {
   for (j = 3 - u; j < 7; j++) {
      l = u - 3 + j;
      fx[j] = coeff[(l*ny + v)*nz + w];
   }  

   for (j = 0; j < 3 - u; j++) {
      uj = u - j;
      fx[2 - uj] = 2*fx[3 -  uj] - fx[4 - uj];
   }
   return;
}
*/

if (u == nx - 3) {
   for (j = 0; j < 6; j++) {
      l = u - 3 + j;
      fx[j] = coeff[(l*ny + v)*nz + w];
   } 
   fx[6] = 2*fx[5] - fx[4];
   return;
}

if (u == nx - 2) {
   for (j = 0; j < 5; j++) {
      l = u - 3 + j;
      fx[j] = coeff[(l*ny + v)*nz + w];
   } 
   fx[5] = 2*fx[4] - fx[3];
   fx[6] = 2*fx[5] - fx[4];
   return;
}

if (u == nx - 1) {
   for (j = 0; j < 4; j++) {
      l = u - 3 + j;
      fx[j] = coeff[(l*ny + v)*nz + w];
   }
   fx[4] = 2*fx[3] - fx[2];
   fx[5] = 2*fx[4] - fx[3];
   fx[6] = 2*fx[5] - fx[4];
   return;
}


}

/*-----------------------------------------------------------------------------*/
void    Get_Seven_Function_Values_Y(float *coeff, int nx, int ny, int nz, int u, int v, int w, float *fx)
{
int   j, l, uj;

if (v >= 3 && v <= ny-4) {
   for (j = 0; j < 7; j++) {
      l = v - 3 + j;
      fx[j] = coeff[(u*ny + l)*nz + w];
   }
   return;
}

if (v == 2) {
   for (j = 1; j < 7; j++) {
      l = j -1;
      fx[j] = coeff[(u*ny + l)*nz + w];
   }
   fx[0] = 2*fx[1] - fx[2];
   return;
}


if (v == 1) {
   for (j = 2; j < 7; j++) {
      l = j-2;
      fx[j] = coeff[(u*ny + l)*nz + w];
   }

   fx[1] = 2*fx[2] - fx[3];
   fx[0] = 2*fx[1] - fx[2];
   return;
}

if (v == 0) {
   for (j = 3; j < 7; j++) {
      l = j -3;
      fx[j] = coeff[(u*ny + l)*nz + w];
   }

   fx[2] = 2*fx[3] - fx[4];
   fx[1] = 2*fx[2] - fx[3];
   fx[0] = 2*fx[1] - fx[2];
   return;
}

/*
if (v <= 2) {
   for (j = 3 - v; j < 7; j++) {
      l = v - 3 + j;
      fx[j] = coeff[(u*ny + l)*nz + w];
   } 

   for (j = 0; j < 3 - v; j++) {
      uj = v - j;
      fx[2 - uj] = 2*fx[3 -  uj] - fx[4 - uj];
   }
   return;
}
*/

if (v == ny - 3) {
   for (j = 0; j < 6; j++) {
      l = v - 3 + j;
      fx[j] = coeff[(u*ny + l)*nz + w];
   }
   fx[6] = 2*fx[5] - fx[4];
   return;
}

if (v == ny - 2) {
   for (j = 0; j < 5; j++) {
      l = v - 3 + j;
      fx[j] = coeff[(u*ny + l)*nz + w];
   }
   fx[5] = 2*fx[4] - fx[3];
   fx[6] = 2*fx[5] - fx[4];
   return;
}

if (v == ny - 1) {
   for (j = 0; j < 4; j++) {
      l = v - 3 + j;
      fx[j] = coeff[(u*ny + l)*nz + w];
   }
   fx[4] = 2*fx[3] - fx[2];
   fx[5] = 2*fx[4] - fx[3];
   fx[6] = 2*fx[5] - fx[4];
   return;
}
}


/*-----------------------------------------------------------------------------*/
void    Get_Seven_Function_Values_Z(float *coeff, int nx, int ny, int nz, int u, int v, int w, float *fx)
{
int   j, l, uj;

if (w >= 3 && w <= nz-4) {
   for (j = 0; j < 7; j++) {
      l = w - 3 + j;
      fx[j] = coeff[(u*ny + v)*nz + l];
   }
   return;
}

if (w == 2) {
   for (j = 1; j < 7; j++) {
      l = j -1;
      fx[j] = coeff[(u*ny + v)*nz + l];
   }
   fx[0] = 2*fx[1] - fx[2];
   return;
}


if (w == 1) {
   for (j = 2; j < 7; j++) {
      l = j-2;
      fx[j] = coeff[(u*ny + v)*nz + l];
   }

   fx[1] = 2*fx[2] - fx[3];
   fx[0] = 2*fx[1] - fx[2];
   return;
}

if (w == 0) {
   for (j = 3; j < 7; j++) {
      l = j -3;
      fx[j] = coeff[(u*ny + v)*nz + l];
   }

   fx[2] = 2*fx[3] - fx[4];
   fx[1] = 2*fx[2] - fx[3];
   fx[0] = 2*fx[1] - fx[2];
   return;
}



/*
if (w <= 2) {
   for (j = 3 - w; j < 7; j++) {
      l = w - 3 + j;
      fx[j] = coeff[(u*ny + v)*nz + l];
   }

   for (j = 0; j < 3 - w; j++) {
      uj = w - j;
      fx[2 - uj] = 2*fx[3 -  uj] - fx[4 - uj];
   }
   return;
}
*/

if (w == nz - 3) {
   for (j = 0; j < 6; j++) {
      l = w - 3 + j;
      fx[j] = coeff[(u*ny + v)*nz + l];
   }
   fx[6] = 2*fx[5] - fx[4];
   return;
}

if (w == nz - 2) {
   for (j = 0; j < 5; j++) {
      l = w - 3 + j;
      fx[j] = coeff[(u*ny + v)*nz + l];
   }
   fx[5] = 2*fx[4] - fx[3];
   fx[6] = 2*fx[5] - fx[4];
   return;
}

if (w == nz - 1) {
   for (j = 0; j < 4; j++) {
      l = w - 3 + j;
      fx[j] = coeff[(u*ny + v)*nz + l];
   }
   fx[4] = 2*fx[3] - fx[2];
   fx[5] = 2*fx[4] - fx[3];
   fx[6] = 2*fx[5] - fx[4];
   return;
}
}


/*-----------------------------------------------------------------------------*/
float   Extreme_Positive(float a, float b)
{
float aa, bb;

aa = a*a;
bb = b*b;

if (a <= b) {
   if (b < 0.0) return(bb);
   if (a > 0.0) return(aa);
   return(0.0);
}

if (aa > bb) return(aa);
return(bb);
}

/*-----------------------------------------------------------------------------*/
float   Extreme_Negative(float a, float b)
{
float aa, bb;

aa = a*a;
bb = b*b;

if (a <= b) {
   if (aa > bb) return(aa);
   return(bb);
}

if (a < 0.0) return(aa);
if (b > 0.0) return(bb);
return(0.0);
}

/*-----------------------------------------------------------------------------*/
float Gradient_2(float *fun, float dx, float dy, float dz, int i, int j, int k, 
               int nx, int ny, int nz)
{
float fx, fy, fz;


// f_ijk = fun[(i*ny + j)*nz + k];

if (i > 0 && i < nx- 1) fx = (fun[((i+1)*ny + j)*nz + k] -  fun[((i-1)*ny + j)*nz + k])/(dx+dx); 
if (i == 0)             fx = (fun[((1)*ny + j)*nz + k]  - fun[j*nz + k])/dx;
if (i == nx - 1)        fx = (fun[((nx-1)*ny + j)*nz + k] - fun[((nx-2)*ny + j)*nz + k])/dx;

if (j > 0 && j < ny- 1) fy = (fun[(i*ny + j+1)*nz + k] -  fun[(i*ny + j-1)*nz + k])/(dy+dy);
if (j == 0)             fy = (fun[(i*ny + 1)*nz + k] -  fun[(i*ny)*nz + k])/dy; 
if (j == ny - 1)        fy = (fun[(i*ny + ny-1)*nz + k] - fun[(i*ny + ny-2)*nz + k])/dy;

if (k > 0 && k < nz- 1) fz = (fun[(i*ny + j)*nz + k+1] -  fun[(i*ny + j)*nz + k-1])/(dz+dz);
if (k == 0)             fz = (fun[(i*ny + j)*nz + 1] -  fun[(i*ny + j)*nz])/dz;             
if (k == nz - 1)        fz = (fun[(i*ny + j)*nz + nz-1] - fun[(i*ny + j)*nz + nz-2])/dz;


return(fx*fx + fy*fy + fz*fz);
}


/*-----------------------------------------------------------------------------*/
void Divided_DifferenceOrder2PartialsAtGridPoint(float *c, float dx, float dy, float dz,
                                             int nx, int ny, int nz,
                                             int u, int v, int w, float *partials)
//              float   *c,    /* the spline  coefficients                            */
//              float   dx     /* spacing in x direction                              */
//              float   dy     /* spacing in y direction                              */
//              float   dz     /* spacing in z direction                              */
//              int     nx,    /* number of samples or coefficients in x direction    */
//              int     ny,    /* number of samples or coefficients in y direction    */
//              int     nz,    /* number of samples or coefficients in z direction    */
//              float   *partial /* partial derivatives                               */
{
float c27[27], indx, indy, indz;
int  i;

Take_27_Coefficients(c, nx, ny, nz, u, v, w, c27);


indx = 1.0/dx;
indy = 1.0/dy;
indz = 1.0/dz;

partials[0] = c27[13];

// fx, fy, fz
partials[1] = (c27[22] - c27[4 ])*indx*0.5;
partials[2] = (c27[16] - c27[10])*indy*0.5;
partials[3] = (c27[14] - c27[12])*indz*0.5;

// fxx, fxy, fxz
partials[4] = (c27[22] + c27[4] - c27[13] - c27[13])*indx*indx;
partials[5] = (c27[25] + c27[1] - c27[19] - c27[7 ])*indx*indy*0.25;
partials[6] = (c27[23] + c27[3] - c27[21] - c27[5 ])*indx*indz*0.25;

// fyy, fyz
partials[7] = (c27[16] + c27[10] - c27[13] - c27[13])*indy*indy;
partials[8] = (c27[17] + c27[9 ] - c27[15] - c27[11])*indy*indz*0.25;

// fzz
partials[9] = (c27[12] + c27[14] - c27[13] - c27[13])*indz*indz;
}

/*--------------------------------------------------------------------------*/
 void   ConvertToInterpolationCoefficients_Qu
                (
                        float   *c,             /* input samples --> output coefficients */
                        int     DataLength,     /* number of samples or coefficients */
                        float   *z,             /* poles */
                        int     NbPoles,        /* number of poles */
                        float   Tolerance       /* admissible relative error */
                )

{ /* begin ConvertToInterpolationCoefficients */

        float   Lambda = 1.0;
        int     n, k;

        /* special case required by mirror boundaries */
        if (DataLength == 1) {
                return;
        }
        /* compute the overall gain */
        for (k = 0; k < NbPoles; k++) {
                Lambda = Lambda * (1.0f - z[k]) * (1.0f - 1.0f / z[k]);
        }
        /* apply the gain */
        for (n = 0; n < DataLength; n++) {
                c[n] *= Lambda;
        }
        /* loop over all poles */
        for (k = 0; k < NbPoles; k++) {
                /* causal initialization */
                c[0] = InitialCausalCoefficient(c, DataLength, z[k], Tolerance);
                /* causal recursion */
                for (n = 1; n < DataLength; n++) {
                        c[n] += z[k] * c[n - 1];
//printf("c^+ Qin = %f\n", c[n]);
                }
                /* anticausal initialization */
                c[DataLength - 1] = InitialAntiCausalCoefficient(c, DataLength, z[k]);
                /* anticausal recursion */
                for (n = DataLength - 2; 0 <= n; n--) {
                        c[n] = z[k] * (c[n + 1] - c[n]);
                }
        }
} /* end ConvertToInterpolationCoefficients */

/*--------------------------------------------------------------------------*/
float   InitialCausalCoefficient
                (
                        float   *c,             /* coefficients */
                        int     DataLength,     /* number of coefficients */
                        float   z,                      /* actual pole */
                        float   Tolerance       /* admissible relative error */
                )

{ /* begin InitialCausalCoefficient */

        float   Sum, zn, z2n, iz;
        int     n, Horizon;

        /* this initialization corresponds to mirror boundaries */
        Horizon = DataLength;
        if (Tolerance > 0.0) {
                Horizon = (int)ceil(log(Tolerance) / log(fabs(z)));
        }
        if (Horizon < DataLength) {
                /* accelerated loop */
                zn = z;
                Sum = c[0];
                for (n = 1; n < Horizon; n++) {
                        Sum += zn * c[n];
                        zn *= z;
                }
                return(Sum);
        }
        else {
                /* full loop */
                zn = z;
                iz = 1.0f / z;
                z2n = (float)pow(z, (DataLength - 1));
                Sum = c[0] + z2n * c[DataLength - 1];
                z2n *= z2n * iz;
                for (n = 1; n <= DataLength - 2; n++) {
                        Sum += (zn + z2n) * c[n];
                        zn *= z;
                        z2n *= iz;
                }
                return(Sum / (1.0f - zn * zn));
        }
} /* end InitialCausalCoefficient */

/*--------------------------------------------------------------------------*/
float  InitialAntiCausalCoefficient
                (
                        float   *c,             /* coefficients */
                        int     DataLength,     /* number of samples or coefficients */
                        float   z                       /* actual pole */
                )

{ /* begin InitialAntiCausalCoefficient */

        /* this initialization corresponds to mirror boundaries */
        return ((z / (z * z - 1.0f)) * (z * c[DataLength - 2] + c[DataLength - 1]));
} /* end InitialAntiCausalCoefficient */


/*--------------------------------------------------------------------------*/
void    Evaluat_Four_Basis(float x, float *values)
{
float y;

values[0] = (1.0 - x)*(1.0 - x)*(1.0 - x)/6.0;
values[1] = 0.666666666 - x*x + 0.5*x*x*x; 
y = 1.0 - x;
values[2] = 0.666666666 - y*y + 0.5*y*y*y; 
y = 2.0 - x; 
values[3] = x*x*x/6.0;
}

/*--------------------------------------------------------------------------*/
void ReSamplingCubicSpline(float *coeff, int nx, int ny, int nz, float *funvalues, 
                           int Nx, int Ny, int Nz)
{
int   i, j, k, l, ix, jy, kz;
float Dx, Dy, Dz, x, y, z;
float values_x[4], values_y[4],values_z[4];
float c64[64], tensor[64], result;

Dx = (nx - 1.0)/(Nx - 1.0);
Dy = (ny - 1.0)/(Ny - 1.0);
Dz = (nz - 1.0)/(Nz - 1.0);

for (i = 0; i < Nx; i++) {
   x = i*Dx;
   ix = x;
   x = x - ix;
   for (j = 0; j < Ny; j++) {
      y = j*Dy;
      jy = y;
      y = y - jy;
      for (k = 0; k < Nz; k++) {
         z = k*Dz;
         kz = z;
         z = z - kz;
         Evaluat_Four_Basis(x, values_x);
         Evaluat_Four_Basis(y, values_y);
         Evaluat_Four_Basis(z, values_z);

         Tensor_444(values_x, values_y,values_z, tensor);
         Take_64_Coefficients(coeff, nx, ny, nz, ix, jy, kz, c64);
         
         result = 0.0;
         for (l = 0; l < 64; l++) {
            result = result + tensor[l]*c64[l];
         }
         funvalues[(i*Ny + j)*Nz + k] = result;
     }
   }
}

}  

/*--------------------------------------------------------------------------*/
float PointToLineDistance(float p[], float p1[], float p2[])
{
    int i;
    float t, tem1,tem2, d1,d2,d3, dis;
    float v1[3], v2[3];

    for (i=0; i<3; i++)  {
       v1[i]=p1[i]-p[i];
       v2[i]=p2[i]-p1[i];
    }

    tem1 = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
    tem2 = v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2];

    if(tem2 < 0.000000001) tem2=0.000000001;
    t=-tem1/tem2;
    
    if((t>0.0) && (t<1.0)) { 

        for(i = 0; i < 3; i++) {
	   v1[i] = v1[i]+t*v2[i];
        }
	dis=v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];

    }  else {
        dis=v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];
        
        for(i=0; i<3; i++) {
           v1[i]=p2[i]-p[i];
        }
 
        d2 = v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];
        if (dis > d2) dis = d2;        
    }
    return(dis);
}

/*--------------------------------------------------------------------------*/
float PointToTriangleDistance(float p[], float p1[], float p2[], float p3[])
{
    int i;
    float u, v, w;
    float t1[3], t2[3], t3[3];
    float a, b, c, d, e, d1,d2, d3, dis;
    float tem, tem1, tem2;
    
    for(i=0; i<3;i++) {
        t1[i]=p[i]  - p3[i];
	t2[i]=p1[i] - p3[i];
	t3[i]=p2[i] - p3[i];
    }
 
    a = 0.0;
    b = 0.0;
    c = 0.0;
    d = 0.0;
    e = 0.0;

    for(i=0; i<3; i++) {
       a = a + t1[i]*t2[i];
       b = b + t1[i]*t3[i];
       c = c + t2[i]*t2[i];
       d = d + t2[i]*t3[i];
       e = e + t3[i]*t3[i];
    }

    tem= c*e - d*d;
    tem1=a*e - b*d;
    tem2=b*c - a*d;

    if(tem<0.00000001) tem=0.00000001;
    u=tem1/tem;
    v=tem2/tem;
    w=1.0 - u - v;

    if((u>0) && (u<1) && (v>0) && (v<1) && (w>0)) {

       for(i=0; i<3; i++) t1[i]=u*t2[i] + v*t3[i] - t1[i];
       dis = t1[0]*t1[0] + t1[1]*t1[1] + t1[2]*t1[2];

    }  else {
       dis=PointToLineDistance(p, p1, p2);
       d2 =PointToLineDistance(p, p1, p3);
       d3 =PointToLineDistance(p, p2, p3);
       if (dis > d2) dis = d2;
       if (dis > d3) dis = d3;
    }
    return(dis);
}



/** Added by zqyork brado **/

/* --------------------------------------------------------------------------*/
void Evaluate_Four_Basis_First_Partial(float x, float *firstPartials)
{
  float y;
  firstPartials[0] = -0.5*(1.0-x)*(1.0-x);
  firstPartials[1] = 1.5*x*x - 2*x;
  y = 1.0 - x;
  firstPartials[2] = -1.5*y*y +2*y;
  firstPartials[3] = 0.5*x*x;
}


/*---------------------------------------------------------------------------*/
void Evaluate_Four_Basis_Second_Partial(float x, float *secondPartials)
{
  secondPartials[0] = 1.0 - x;
  secondPartials[1] = 3*x -2;
  secondPartials[2] = 1 - 3*x;
  secondPartials[3] = x;
}

void ReSamplingCubicSpline1OrderPartials(float *coeff, int nx, int ny, int nz, float *dxyz, float *minExtent, float *zeroToTwoPartialValue,  float *p) 
    /* p the point to be evaluated on                                          */
    /* coeff the spline coefficients                                           */
    /* firstParValue, 3 directional derivatives                         */
{
  int i, j, k, l;
  int u, v, w;
  float x, y, z;
  float funcValue_x[4], funcValue_y[4], funcValue_z[4];
  float f1Partial_x[4], f1Partial_y[4], f1Partial_z[4], f2Partial_xx[4], f2Partial_xy[4], f2Partial_xz[4], f2Partial_yy[4], f2Partial_yz[4], f2Partial_zz[4];
  float c64[64], tensor[64], tensorx[64], tensory[64], tensorz[64],tensorxx[64],tensorxy[64],tensorxz[64],tensoryy[64],tensoryz[64],tensorzz[64];
  float result[10];

  float dx = dxyz[2];
  float dy = dxyz[1];
  float dz = dxyz[0];

  float indx = 1.0/dx;
  float indy = 1.0/dy;
  float indz = 1.0/dz;


  u = (int)((p[0]-minExtent[0] )*indx);
  v = (int)((p[1]-minExtent[1] )*indy);
  w = (int)((p[2]-minExtent[2] )*indz);


 //  u = (int)((p[0] - minExtent[0])/dx);
 //  v = (int)((p[1] - minExtent[1])/dy);
 //  w = (int)((p[2] - minExtent[2])/dz);

  // cout<<" P coordinate: "<< p[0] << "  "<<p[1]<<" "<<p[2]<<"   u, v, w = "<< u<<" "<<v<<" "<<w;
  // cout<<"   dxyz="<<dx; 
  //
 
  float fu = (p[0] - minExtent[0])/dx;
  float fv = (p[1] - minExtent[1])/dy;
  float fw = (p[2] - minExtent[2])/dz;
 //printf("fu fv fw:%f %f %f\n",fu,fv,fw);

  if (fu - u > 0.999)   u = u+1;
  if (fv - v > 0.999)   v = v+1;
  if (fw - w > 0.999)   w = w+1;



  x = (p[0] - minExtent[0] - u*dx)/dx;
  y = (p[1] - minExtent[1] - v*dy)/dy;
  z = (p[2] - minExtent[2] - w*dz)/dz;

  if(fabs(x) < 0.0001) x = 0.0;
  if(fabs(y) < 0.0001) y = 0.0;
  if(fabs(z) < 0.0001) z = 0.0;
  
  if(x > 0.9999) x = 1.0;
  if(y > 0.9999) y = 1.0;
  if(z > 0.9999) z = 1.0;
   

  
  Evaluat_Four_Basis(x,funcValue_x);
  Evaluat_Four_Basis(y,funcValue_y);
  Evaluat_Four_Basis(z,funcValue_z);
 
  Evaluate_Four_Basis_First_Partial(x, f1Partial_x);
  Evaluate_Four_Basis_First_Partial(y, f1Partial_y);
  Evaluate_Four_Basis_First_Partial(z, f1Partial_z);

  Evaluate_Four_Basis_Second_Partial(x, f2Partial_xx);
  Evaluate_Four_Basis_Second_Partial(y, f2Partial_yy);
  Evaluate_Four_Basis_Second_Partial(z, f2Partial_zz);



  Tensor_444(funcValue_x, funcValue_y, funcValue_z, tensor);
  Tensor_444(f1Partial_x, funcValue_y, funcValue_z, tensorx);
  Tensor_444(funcValue_x, f1Partial_y, funcValue_z, tensory);
  Tensor_444(funcValue_x, funcValue_y, f1Partial_z, tensorz);

  Tensor_444(f2Partial_xx, funcValue_y, funcValue_z, tensorxx);
  Tensor_444(f1Partial_x, f1Partial_y, funcValue_z, tensorxy);
  Tensor_444(f1Partial_x, funcValue_y, f1Partial_z, tensorxz);
  Tensor_444(funcValue_x, f2Partial_yy, funcValue_z, tensoryy);
  Tensor_444(funcValue_x, f1Partial_y, f1Partial_z, tensoryz);
  Tensor_444(funcValue_x, funcValue_y, f2Partial_zz, tensorzz);



  
//   x = (p[0] - minExtent[0] - u*dx)/dx; 
//   y = (p[1] - minExtent[1] - v*dy)/dy;
 //  z = (p[2] - minExtent[2] - w*dz)/dz;
//    cout<<"x, y, z   "<< x<<",  " <<  y <<",   " << z<<endl;
  if (( u  > 0  && u < nx-2) &&
        ( v  > 0  && v < ny-2) &&
        ( w  > 0  && w < nz-2) ){
  
    Take_64_Coefficients(coeff, nx, ny, nz, w, v, u, c64);
         
    for(l = 0; l < 10; l++)
    {
      result[l] = 0.0;
    }
    
    for(l = 0; l < 64; l++)
    {
      result[0] += tensor[l]*c64[l];
      result[1] += tensorx[l]*c64[l]; //*indx;
      result[2] += tensory[l]*c64[l]; //*indy;
      result[3] += tensorz[l]*c64[l]; //*indz;

      result[4] += tensorxx[l]*c64[l]; //*indx*indx;
      result[5] += tensorxy[l]*c64[l]; //*indx*indy;
      result[6] += tensorxz[l]*c64[l]; //*indx*indz;
      result[7] += tensoryy[l]*c64[l]; //*indy*indy;
      result[8] += tensoryz[l]*c64[l]; //*indy*indz;
      result[9] += tensorzz[l]*c64[l]; //*indz*indz;
    }
 
    for(l = 0; l < 10; l++)
    {
      zeroToTwoPartialValue[l] = result[l];
    }
        }
}  

