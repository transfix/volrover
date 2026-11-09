//============================================================================
/* Copyright (c) The Technical University of Denmark
 * Author: Ojaswa Sharma
 * E-mail: os@imm.dtu.dk
 * File: anisotropic_diffusion_kerenels.cu
 * Interface to libMSLevelSet.a
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

#ifndef __LEVELSET3D__
#define __LEVELSET3D__

#include <cutil.h>
#include <cuda.h>
#include <driver_types.h>
#include <vector_types.h>

// Other clases, newvolumemainwindow, and MSLSdialog box need this.
struct MSLevelSetParams
{
  //LevelSet Params
  float lambda1, lambda2, mu, nu, deltaT, epsilon;
  unsigned int nIter, DTWidth, medianIter;
  float medianTolerance;
  float superEllipsoidPower; 
  int BBoxOffset; 
  //CUDA params
  int subvolDim;
 //,subvolDimSDT,subvolDimAvg,subvolDimPDE;
  int PDEBlockDim, avgBlockDim, SDTBlockDim, medianBlockDim;
  int init_interface_method;
  int multi_init_r;
  int multi_init_dr;
  int multi_init_s;  
  int reint_niter;
  float volval_min;
  float volval_max;
};

//extern void qDebug(const char * format, ...);

namespace MumfordShahLevelSet
{
  #define TINY 0.0000001
  #define PI 3.14159265f

  #define CHK_CUDA_ERR(_str)    cuda_err = cudaGetLastError(); \
  if(cuda_err != cudaSuccess) { \
  fprintf(stderr,"%s: %s\n", _str, cudaGetErrorString(cuda_err)) ; \
  return(false); \
  }

#if __DEVICE_EMULATION__
#  define CUDA_DEVICE_INIT()
#else
#  define CUDA_DEVICE_INIT() {                                      \
    int deviceCount;                                                         \
    CUDA_SAFE_CALL_NO_SYNC(cudaGetDeviceCount(&deviceCount));                \
    if (deviceCount == 0) {                                                  \
        fprintf(stderr, "cuda error: no devices supporting CUDA.\n");       \
        return(false);                                                  \
    }                                                                        \
    cudaDevice = -1;                                                             \
    int cudaCapableDevices=0; \
    cudaDeviceProp deviceProp;                                               \
    for(int i = 0; i<deviceCount; i++) { \
      CUDA_SAFE_CALL_NO_SYNC(cudaGetDeviceProperties(&deviceProp, i));       \
      if (deviceProp.major < 1) {                                              \
        fprintf(stderr, "cuda error: device %d: %s does not support CUDA.\n", i, deviceProp.name);     \
        return(false);                                                  \
      }                                                                        \
      else { \
        cudaCapableDevices++; \
        if(cudaDevice == -1) cudaDevice = i; \
      } \
    } \
    fprintf(stderr, "%d/%d cuda capable devices available. Using device %d: %s\n", cudaCapableDevices, deviceCount, cudaDevice, deviceProp.name);       \
    CUDA_SAFE_CALL(cudaSetDevice(cudaDevice));                                      \
}
#endif

  #define VOXEL_AT(vol, x, y, z) (*(vol + x + (y)*volSize.width + (z)*volSize.width*volSize.height)) //Strictly for use with full volume
  #define SUB2IND(x, y, z) (x + (y)*(int)volSize.width + (z)*(int)volSize.width*(int)volSize.height) 
  #define CPUMIN(a, b) (a<b)?a:b
  #define CPUMAX(a, b) (a>b)?a:b
  #define COMPUTE_AVERAGE_ON_CPU 1
  #define UPDATE_PDE_ON_CPU 0
  #define INIT_INTERFACE_ON_CPU 0


  //#define debugPrint(format, ...) qDebug(format, ## __VA_ARGS__)
  #define debugPrint(format, ...) fprintf(stderr, format, ## __VA_ARGS__)
  
  struct DataInfo
  {
    unsigned int n[3]; // Number of elements in each dimension
    unsigned int n_input[3]; // NUmber of elements in the input volume
    unsigned long nTotal; // Total number of voxels = n[0]*n[1]*n[3]
    char dataType[256]; // data type of values. -ve value indicated signed data
    bool volumeCovered; //cover the input volume with a one voxel layer around - binary flag.
  };

  class MSLevelSet
  {
    //Member variables.
    protected:
    //CUDA variables
    int cudaDevice;
    cudaExtent volSize;// = make_cudaExtent(2*(subvolDim - 2) + 2, 3*(subvolDim - 2) + 2, 4*(subvolDim - 2) + 2); // Test size for a large volume that can encomapss 8x4x2 subvolumes
    cudaExtent subvolSize;// = make_cudaExtent(subvolDim, subvolDim, subvolDim);
    cudaExtent inputVolSize;
    cudaExtent subvolIndicesExtents;
    //const cudaExtent subvolIndicesExtents = make_cudaExtent(iDivUp(volSize.width-2, subvolDim-2), iDivUp(volSize.height-2, subvolDim-2), iDivUp(volSize.depth-2, subvolDim-2));
    DataInfo datainfo;// Stores information about the input volume

    float* h_vol; // The full volume
    float* h_phi; // The phi function 
    float* h_buffer[2];
    unsigned int currentPhi;
    unsigned int timer;
    float milliseconds;
    cudaError_t cuda_err;

    //CUDA parameters
    bool DoInitializeCUDA;
    int subvolDim;//, subvolDimSDT,subvolDimAvg, subvolDimPDE;
    int PDEBlockDim, avgBlockDim, SDTBlockDim;
    enum {BBOX, SUPER_ELLIPSOID};
    int init_interface_method;

    float* h_subvolSpare; //Spare subvolume on the host for ntermediate computations. N.B: Size: (subvolDim)^3
    cudaArray *d_phiArray; // subvolume bound to texture - used for texture reads in phi    
    cudaArray *d_volArray; // subvolume bound to texture - used for texture reads in intensitintensity
    cudaPitchedPtr d_volPPtr;// subvolume in device memory - used for updates in phi
    cudaPitchedPtr d_spare1PPtr, d_spare2PPtr, d_spare3PPtr; // 2 spare subvolumes to store ntermediae results of computation
    size_t tPitch;

    bool zeroOutArray;

    //Level set parameters
    float mu;
    float h;
    float nu;
    float lambda1;
    float lambda2;
    float delta_t;
    float c1, c2; 
    float old_0_c, old_1_c;
    int nIter, DTWidth;
    bool converged;
    float epsilon;
    int iter_num;
    int BBoxOffset;
    float superEllipsoidPower;
    int multi_init_r;
    int multi_init_dr;
    int multi_init_s;  
    int reinit_niter;
    float volval_min;
    float volval_max;
    
    //Member functions
    public:
    double getTime();
    MSLevelSet(bool doInitCUDA =true);
    // Interface to volRover
    // Returns false on error. true otherwise.
    virtual bool runSolver(	float *vol,
                   int _width, int _height, int _depth, MSLevelSetParams *MSLSParams,
                   void (*evolutionCallback)(const float *vol, int dimx, int dimy, int dimz) = NULL);
    
    protected:
    int iDivUp(int a, int b);
    unsigned long inKB(unsigned long bytes);
    unsigned long inMB(unsigned long bytes);
    void printStats(unsigned long free, unsigned long total);
    void printMemInfo();
    void copy3DHostToArray(float *_src, cudaArray *_dst, cudaExtent copy_extent, cudaExtent src_extent, cudaPos src_offset);
    void copy3DHostToMem(float *_src, cudaPitchedPtr _dst, cudaExtent copy_extent, cudaExtent src_extent, cudaPos src_offset);
    void copy3DMemToArray(cudaPitchedPtr _src, cudaArray *_dst);
    void copy3DMemToHost(cudaPitchedPtr _src, float *_dst, cudaExtent copy_extent, cudaExtent dst_extent, cudaPos src_offset, cudaPos dst_offset);
    void copy3DArrayToHost(cudaArray *_src, float *_dst, cudaExtent copy_extent, cudaExtent dst_extent, cudaPos src_offset, cudaPos dst_offset);
    void adjustDnloadSubvolSize(int _x, int _y, int _z, unsigned int offset, cudaExtent &copyvol); //Download to Host
    void adjustUploadSubvolSize(int _x, int _y, int _z, unsigned int offset, cudaExtent &copyvol_upload); //upload to GPU
    void writeRawData(float *data, char* filename);
    void writeSlices(float *data, char* filename_prefix, int nslices);
    
    virtual int solverMain();//entry point into the CUDA solver. 
    virtual  bool solve();
    virtual void cleanUp();
    virtual void freeHostMem();
    virtual bool initCuda();
    
    virtual bool normalizeVolume();//To lie in the range [0 255]
    virtual bool initInterfaceMulti();
    virtual bool computeSDTEikonal();
    virtual bool computeAverageIntensities();
    virtual bool PDEUpdate();
    void setParameters(int _width, int _height, int _depth, MSLevelSetParams *MSLSParams);
  };
}
#endif
