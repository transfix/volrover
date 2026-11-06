//============================================================================
/* Copyright (c) The Technical University of Denmark
 * Developed at the Computational Visualization Center (CVC), The University of Texas at Austin
 * Author: Ojaswa Sharma
 * E-mail: os@imm.dtu.dk
 * File: multiphaseseg.h
 * Multi-phase segmentation using Higher-order mumford shah and other functionals.
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

#ifndef __MULTIPHASE_SEG__
#define __MULTIPHASE_SEG__

#include <cutil.h>
#include <cuda.h>
#include <driver_types.h>
#include <cutil_inline.h>
#include <cutil_math.h>
#include <vector_types.h>
#include "../HigherOrderSegmentation/higherorderseg.h"

using namespace HigherOrderSegmentation;

struct MPLevelSetParams : public MSLevelSetParams
{
  int nImplicits;
  char userSegFileName[512];
  //for a fine grain control we need to be set the subvoldim differently for different kernels
  int subvolDimSDT, subvolDimPDE, subvolDimAvg;
};

namespace MultiphaseSegmentation
{
  #define MAX_3D_TEX_SIZE 2048 // As of CUDA 2.1 (in all the three dimensions)
  #define INIT_MULTIPHASE_INTERFACE_ON_CPU 0
  #define COMPUTE_MULTPHASE_AVERAGE_ON_CPU 1
  #define UPDATE_MULTIPHASE_PDE_ON_CPU 0
  #define AVERAGE_INTENSITY_H_PRODUCT_THREASHOLD 0.000001f
  
  struct Factors
  {
    unsigned int x;
    unsigned int y;
    unsigned int z;
  };
  
  struct Coordinate3D
  {
    float x;
    float y;
    float z;
  };
  
  class MPSegmentation:public HOSegmentation
  {
    //member variables
    protected:    
    //multiphase variables
    unsigned int nImplicits;// Number of implcit functions. i.e., phi functions
    unsigned int nImplicitsX, nImplicitsY, nImplicitsZ;//optimal decomposition of nImplicits, such that  (nImplicitsX*nImplicitsY*nImplicitsZ - nImplicits) is minimized
    float **h_PHI; //vector PHI of nImplicits phi funtions.
    float **h_BUFFER[2];//spare arrays
    //No need to declare d_phiArray again here. It should just be fine to use parent member's variable. It is not allocated the same way anyway!
    //cudaArray *d_phiArray; // subvolume bound to texture - used for texture reads in phi. This is a concatenated array for all phis [phi1][phi2]...[phin]. 
    //Concatenation should take care of maximum 3D texture size along any dimension.
    float *c_avg;//average values for different classes
    float *d_c_avg;
    unsigned long nelem;
    int subvolDimSDT, subvolDimPDE, subvolDimAvg;

    //member functions
    public:
    MPSegmentation(bool doInitCUDA = true);    
   // Interface to volRover
   // Returns false on error. true otherwise.
   using MSLevelSet::runSolver;
   virtual bool runSolver(	float *vol, float **phi, int _width, int _height, int _depth, MPLevelSetParams *MPLSParams);
    protected:
    using MSLevelSet::solverMain;
    virtual int solverMain();//entry point into the CUDA solver.
    virtual  bool solve();

    bool allocateHost(float *vol, float **phi);
    bool freeHost();
    bool allocateDevice();
    bool freeDevice(int level = -1);
    virtual bool initCuda();

    virtual bool PDEUpdate();
    bool multiPhasePDEUpdate(bool onlyborder);
    bool computeOptimalFactors(unsigned int i, unsigned int n, unsigned int &a, unsigned int &b, unsigned int &c); //Solves an integer programming problem. See function documentation  
    
    virtual bool initInterfaceMulti();
    virtual bool computeSDTEikonal();
    bool computeAverageIntensities();
    using MSLevelSet::copy3DHostToArray;
    void copy3DHostToArray(float *_src, cudaArray *_dst, cudaExtent copy_extent, cudaExtent src_extent, cudaPos src_offset, cudaPos dst_offset);
    using MSLevelSet::setParameters;
    void setParameters(int _width, int _height, int _depth, MPLevelSetParams *MPLSParams);
    float divergenceOfNormalizedGradient_FiniteDifference( int x,  int y,  int z);
  };
}
#endif
