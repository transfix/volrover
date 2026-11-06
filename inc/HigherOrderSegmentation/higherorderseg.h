//============================================================================
/* Copyright (c) The Technical University of Denmark
 * Developed at the Computational Visualization Center (CVC), The University of Texas at Austin
 * Author: Ojaswa Sharma
 * E-mail: os@imm.dtu.dk
 * File: higherorderseg.h
 * Segmentation using Higher-order mumford shah and other functionals.
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

#ifndef __HIGHER_ORDER_SEG__
#define __HIGHER_ORDER_SEG__

#include <cutil.h>
#include <cuda.h>
#include <driver_types.h>
#include <vector_types.h>
#include "../MSLevelSet/levelset3D.h"

using namespace MumfordShahLevelSet;

namespace HigherOrderSegmentation
{

  #define COMPUTE_COEFFICIENTS_ON_CPU 1
  #define UPDATE_HIGHER_ORDER_PDE_ON_CPU 0
  
  #define EPSILON 3 //Redefinition here. also defined in ../Curation/smax.cpp
  #define SAFE_VOXEL_AT(vol, x, y, z) (( ((x)>=0) && ((x)<volSize.width) && ((y)>=0) && ((y)<volSize.height) && ((z)>=0) && ((z)<volSize.depth) )? (*(vol + x + (y)*volSize.width + (z)*volSize.width*volSize.height)) : 0.0 )
  class HOSegmentation:public MSLevelSet
  {
    protected:    
    float *h_coeff;//cubic coefficients
    cudaArray *d_coeffArray; //For texture
    float h_cubicDerivatives[27*9];
    int CCblockDim;
    
    //member functions
    public:
    HOSegmentation(bool doInitCUDA = true);
   // Interface to volRover
   // Returns false on error. true otherwise.
   virtual bool runSolver(	float *vol,
                             int _width, int _height, int _depth, MSLevelSetParams *MSLSParams,
                             void (*evolutionCallback)(const float *vol, int dimx, int dimy, int dimz) = NULL);
    
    protected:
    virtual int solverMain();//entry point into the CUDA solver. 
    virtual  bool solve();
    virtual void cleanUp();
    virtual void freeHostMem();
    virtual bool initCuda();
    virtual bool PDEUpdate();


    virtual bool computeCubicCoefficients();
    bool getAvailbleGPUSliceSize(cudaExtent &gpuSlice);
    bool allocateSliceOnGPU(cudaArray*& d_coeffArray, cudaPitchedPtr &d_coeffPPtr, cudaExtent coeffSubvolSize);
    void adjustUploadPDESubvolSizeAndOffset(int _x, int _y, int _z, cudaExtent &copyvol_upload, cudaPos &offset_upload) ;
    bool getFreeGPUMem(unsigned int &free, unsigned int &total);
    float evaluateCubicSplineAtGridPoint(int x, int y, int z);
    float* accessX(float *arr, int y, int z, int k);
    float* accessY(float *arr, int z, int x, int k);
    float* accessZ(float *arr, int x, int y, int k);
    void cubicCoeff1D(int n, int _x1, int _x2, float z1_2n, float* c_plus, float K, float z1, float *readBuffer, float *writeBuffer, float* (HOSegmentation::*accessDim)(float*, int, int, int)) ;   
    using MSLevelSet::copy3DMemToArray;
    void copy3DMemToArray(cudaPitchedPtr _src, cudaArray *_dst, cudaExtent copy_extent);
    void copy3DArrayToHost(cudaArray *_src, float *_dst, cudaExtent copy_extent, cudaExtent dst_extent, cudaPos dst_offset);
    void copy3DArrayToMem(cudaArray *_src, cudaPitchedPtr _dst, cudaExtent copy_extent);
    float divergenceOfNormalizedGradient_FiniteDifference( int x,  int y,  int z);    
    float divergenceOfNormalizedGradient_CubicSpline( int x,  int y,  int z);
    float evaluateCubicSplineGenericDerivativeAtGridPoint(int x, int y, int z, float *xGridSplineValues, float *yGridSplineValues, float *zGridSplineValues);  
    
    //----------------------------------------------------------------------------------------------------------------------------------------------------------
    //Following functions are borrowed from ../HLevelSet/KLevelSet_Recon2.cpp to compare against the divergence computed here.
    void  EvaluateCubicSplineOrder2PartialsAtGridPoint(float *c, float dx, float dy, float dz, int nx, int ny, int nz, int u, int v, int w, float *partials);
    void Take_27_Coefficients(float *c, int nx, int ny, int nz, int u, int v, int w, float *c27);
    float TakeACoefficient_Fast(float *c, int nx, int ny, int nz, int u, int v, int w);
    float TakeACoefficient_Slow(float *c, int nx, int ny, int nz, int u, int v, int w);
    void ComputeTensorXYZ(float *TensorF, float *TensorFx, float *TensorFy, float *TensorFz, float *TensorFxx, float *TensorFxy, float *TensorFxz, float *TensorFyy, float *TensorFyz, float *TensorFzz);
    void   Tensor_333(float *xx, float *yy, float *zz, float *result);
    //----------------------------------------------------------------------------------------------------------------------------------------------------------
     
    float evalCubic(float x);
    float evalCubic_Dx(float x);
    float evalCubic_Dxx(float x);
    float evalTriCubic(float x, float y, float z);
    float evalTriCubic_Dx(float x, float y, float z);
    float evalTriCubic_Dy(float x, float y, float z);
    float evalTriCubic_Dz(float x, float y, float z);
    float evalTriCubic_Dxx(float x, float y, float z);
    float evalTriCubic_Dyy(float x, float y, float z);
    float evalTriCubic_Dzz(float x, float y, float z);
    float evalTriCubic_Dxy(float x, float y, float z);
    float evalTriCubic_Dyz(float x, float y, float z);
    float evalTriCubic_Dzx(float x, float y, float z);
    
    virtual bool initInterfaceMulti();
    virtual bool computeSDTEikonal();    
    virtual bool computeAverageIntensities();    
    virtual bool normalizeVolume();
  };
}
#endif
