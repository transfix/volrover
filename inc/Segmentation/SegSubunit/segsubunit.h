/*
 *  segsubunit.h
 *  
 *
 *  Created by Jose  Rivera on 6/28/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef SEGSUBUNIT_H
#define SEGSUBUNIT_H

#include <xmlrpc/XmlRpc.h>

namespace SegSubunit {

typedef struct {
  float x;
  float y;
  float z;
}VECTOR;

typedef struct {
  unsigned short x;
  unsigned short y;
  unsigned short z;
}INTVECT;

typedef struct {
  float sx;
  float sy;
  float sz;
  float ex;
  float ey;
  float ez;
}DB_VECTOR;

typedef struct CriticalPoint CPNT;
struct CriticalPoint{
  unsigned short x;
  unsigned short y;
  unsigned short z;
  CPNT *next;
};

typedef struct {
  float trans;
  float rotat;
  float angle;
}CVM;

void FindCriticalPoints(int,int,int,float*,CPNT**,float,int,int);
void LocalSymmetry(int, int, int, float*, float*, float*, CPNT*, int,int,DB_VECTOR*,
		   VECTOR*, DB_VECTOR*, int*,int*,int);
void LocalSymmetryRefine(int,int,int,float *,CPNT *,unsigned short *, VECTOR *,
			 int, DB_VECTOR *,FILE*,FILE*,int,int, float *, float *, float);
void swap_buffer(char *buffer, int count, int typesize);
void SubunitSegment(int, int, int, float*, float*, float *, unsigned short*,float, 
		    VECTOR *, DB_VECTOR *, int, int, int, int);
void SubunitSegmentRefine(int, int, int, float*, float*, float *, unsigned short*,float, 
			  VECTOR *, DB_VECTOR *, int, int, CVM*,int, int);
void CoVarianceRefine(int,int,int,float *,CPNT *,unsigned short *,
		      DB_VECTOR *,int,CVM*,int,FILE*);
void ComputeSubAverage(int, int, int, float *,unsigned short *,
		       DB_VECTOR *,int,int,CVM*,int,float *,float *,FILE *);
void MakeTransformMatrix(int, int, int,VECTOR *,DB_VECTOR *,int,int,
			 CVM*,float *, float *, FILE *);
void Write5fSubunit(int,int,int,float *,unsigned short *,float *,float *,FILE *);
void AsymSubunitSegment(int,int,int,float *, float*, float *,unsigned short *, 
			float, float, CPNT *,VECTOR *);
void read_data(int* , int*, int*, float**, float*,float*,char*,int,int);
void write_rawiv_char(unsigned char*, FILE*);
void write_rawiv_short(unsigned short*, FILE*);
void write_rawiv_float(float*, FILE*);
void write_rawv(float*,unsigned short *,int, int, float*, float*, int, FILE*);
void write_asym_rawv(float*,unsigned short *,float*, float*, FILE*);
//void write_asym_rawv(float*,unsigned short *, FILE*);

};

int virusSegSubunit(XmlRpc::XmlRpcValue &params, XmlRpc::XmlRpcValue &result);

#endif
