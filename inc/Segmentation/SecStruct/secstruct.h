/******************************************************************************
				Copyright   

This code is developed within the Computational Visualization Center at The 
University of Texas at Austin.

This code has been made available to you under the auspices of a Lesser General 
Public License (LGPL) (http://www.ices.utexas.edu/cvc/software/license.html) 
and terms that you have agreed to.

Upon accepting the LGPL, we request you agree to acknowledge the use of use of 
the code that results in any published work, including scientific papers, 
films, and videotapes by citing the following references:

C. Bajaj, Z. Yu, M. Auer
Volumetric Feature Extraction and Visualization of Tomographic Molecular Imaging
Journal of Structural Biology, Volume 144, Issues 1-2, October 2003, Pages 
132-143.

If you desire to use this code for a profit venture, or if you do not wish to 
accept LGPL, but desire usage of this code, please contact Chandrajit Bajaj 
(bajaj@ices.utexas.edu) at the Computational Visualization Center at The 
University of Texas at Austin for a different license.
******************************************************************************/


#ifndef SECSTRUCT_H
#define SECSTRUCT_H

#include <xmlrpc/XmlRpc.h>

namespace SecStruct {

typedef struct {
  float x;
  float y;
  float z;
}VECTOR;

typedef struct PNTSIndex PNTS;
struct PNTSIndex{
  float x;  
  float y;  
  float z; 
  PNTS *next;
};

typedef struct {
  float dx;
  float dy;
  float dz;
  float thickness;
  int size;
  PNTS *points;
}SEEDS;

void Diffuse(float *, int, int, int);
void GVF_Compute(int, int, int, float *, VECTOR*);
void Canny(int, int, int, float, float, float*);
void StructureTensor(int, int, int, VECTOR*,unsigned char*,float*,float*,
		     float,float,float,float,float,float,int*);
void read_data(int *, int *, int *, float **, float *, float *, const char *);
void write_rawiv_char(unsigned char*, FILE*);
void write_rawiv_short(unsigned short*, FILE*);
void write_rawiv_float(float*, FILE*);
void HelixHunter(int,int,int,float*,float*,float *,unsigned char *,
		 VECTOR *,float,float,float,FILE*,SEEDS*, int);
void SheetHunter(int,int,int,float*,float*,float *,unsigned char *,VECTOR *,
		 float,float,float,float,float,float,float,float,float,FILE*,SEEDS*,int,FILE*,FILE*);

};

int secondaryStructureDetection(XmlRpc::XmlRpcValue &params, XmlRpc::XmlRpcValue &result);

#endif
