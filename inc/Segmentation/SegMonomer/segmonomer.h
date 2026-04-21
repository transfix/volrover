/*
 *  segmonomer.h
 *  
 *
 *  Created by Jose  Rivera on 6/28/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef SEGMONOMER_H
#define SEGMONOMER_H

#include <xmlrpc/XmlRpc.h>

namespace SegMonomer {

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

typedef struct CriticalPoint CPNT;
struct CriticalPoint{
  unsigned short x;
  unsigned short y;
  unsigned short z;
  CPNT *next;
};



void GVF_Compute(int, int, int, float *, float *, CPNT **);
void read_data(int* , int*, int*, float**, float*,float*,char*,int,int);
void write_rawiv_char(unsigned char*, FILE*);
void write_rawiv_short(unsigned short*, FILE*);
void write_rawiv_float(float*, FILE*);
void get_header(int *, int *, int *, float *, FILE *);
void MonomerSegment(int, int, int, float*, float*, float *, float *, CPNT *,
		    unsigned char *, DB_VECTOR *, int);
void get_average(int, int, int, float*, float*, float*, DB_VECTOR *, int);

};

int virusSegMonomer(XmlRpc::XmlRpcValue &params, XmlRpc::XmlRpcValue &result);

#endif
