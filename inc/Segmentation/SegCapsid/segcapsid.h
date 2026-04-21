/*
 *  segcapsid.h
 *  
 *
 *  Created by Jose  Rivera on 6/28/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef SEGCAPSID_H
#define SEGCAPSID_H

#include <xmlrpc/XmlRpc.h>

namespace SegCapsid {

typedef struct {
	float x;
	float y;
	float z;
}VECTOR;

typedef struct {
	float sx;
	float sy;
	float sz;
	float ex;
	float ey;
	float ez;
}DB_VECTOR;

typedef struct {
	unsigned short x;
	unsigned short y;
	unsigned short z;
}INTVECT;

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

void Diffuse(int, int, int, float*);
void FindCriticalPoints(int,int,int,float*,CPNT**,float,int,int);
void GlobalSymmetry(int, int,int, float*, CPNT*, VECTOR*, FILE*);
void SimpleCapsidSegment(int, int,int, float*, float, int, int, int, VECTOR*);
void SingleCapsidSegment(int,int,int,float *,float, VECTOR *,
						 int, int, int,int, int, int);
void DoubleCapsidSegment(int,int,int,float,float *,unsigned short *, 
						 VECTOR *, float, float);
void CapsidSegmentScore(int,int,int,float *,float *,VECTOR *,DB_VECTOR *,
						CVM*,int,int, float,float);
void CapsidSegmentMarch(int,int,int,float,float *,float *,unsigned short *, 
						VECTOR *, DB_VECTOR *,CVM*,int,int,float, float);
void swap_buffer(char *buffer, int count, int typesize);
void read_data(int* , int*, int*, float**, float*,float*,char*,int,int);
void write_rawiv_char(unsigned char*, FILE*);
void write_rawiv_short(unsigned short*, FILE*);
void write_rawiv_float(float*, FILE*);

};

int virusSegCapsid(XmlRpc::XmlRpcValue &params, XmlRpc::XmlRpcValue &result);

#endif
