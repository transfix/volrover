/*
 *  SegCapsid.cpp
 *  
 *
 *  Created by Jose  Rivera on 6/21/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>

#include <xmlrpc/XmlRpc.h>

#define _LITTLE_ENDIAN 1

#define IndexVect(i,j,k) ((k)*xdim*ydim + (j)*xdim + (i))
#define max2(x, y) ((x>y) ? (x):(y))
#define min2(x, y) ((x<y) ? (x):(y))
#define PolyNum       15
#define PIE           3.1415926536f
#define MAX_STRING    256

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

extern "C" {
	void Diffuse(int, int, int, float*);
	void FindCriticalPoints(int,int,int,float*,CPNT**,float,int,int);
	void GlobalSymmetry(int, int,int, float*, CPNT*, VECTOR*, FILE*);
	void SimpleCapsidSegment(int, int,int, float*, float, CPNT*, VECTOR*);
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

using namespace XmlRpc;

XmlRpcServer s;

class SegCapsid : public XmlRpcServerMethod
{
public:
	SegCapsid(XmlRpcServer* s) : XmlRpcServerMethod("SegCapsid", s) {}
	
	void execute(XmlRpcValue& params, XmlRpcValue& result)
	{
		int xdim,ydim,zdim;
		float *dataset;
		FILE *fp;
		char file_name[MAX_STRING];
		VECTOR *fivefold;
		CPNT *critical_list;
		float *span_tmp, *orig_tmp;
		float tlow,score;
		int i,j,k,type;
		float fx,fy,fz;
		float gx,gy,gz;
		time_t t1,t2;
		int sx1,sy1,sz1;
		int sx2,sy2,sz2;
		float *sym_score;
		unsigned short *_result;
		DB_VECTOR *sixfold,*threefold;
		float small_radius,large_radius;
		int numfold6,numfold5,numfold3;
		int numaxis6,numaxis3;
		CVM *coVarMatrix3,*coVarMatrix6;
		std::string filename;
		
		filename = std::string(params[0]);
		type = int(params[1]);
		tlow = double(params[2]);
		
		printf("begin reading rawiv.... \n");
		span_tmp = (float *)malloc(sizeof(float)*3);
		orig_tmp = (float *)malloc(sizeof(float)*3);
		(void)time(&t1); 
		read_data(&xdim,&ydim,&zdim,&dataset,span_tmp,orig_tmp,const_cast<char *>(filename.c_str()),0,1);
		printf("xdim: %d  ydim: %d  zdim: %d \n",xdim,ydim,zdim);
		(void)time(&t2); 
		printf("time to read dataset: %d seconds. \n\n",(int)(t2-t1));
		
		if (type == 0 || type == 1) 
		{
			printf("begin diffusion .... \n");
			(void)time(&t1); 
			Diffuse(xdim,ydim,zdim,dataset);
			(void)time(&t2); 
			printf("time to diffuse dataset: %d seconds. \n\n",(int)(t2-t1));
			
			printf("begin finding critical points....\n");
			(void)time(&t1); 
			FindCriticalPoints(xdim,ydim,zdim,dataset,&critical_list,tlow,3,2);
			(void)time(&t2); 
			printf("time to find critical points: %d seconds. \n\n",(int)(t2-t1));
			
			fivefold = (VECTOR*)malloc(sizeof(VECTOR)*12);
			printf("begin global symmetry detection ....\n");
			strcpy(file_name,filename.c_str());
			for(i = 0; i<MAX_STRING; i++) {
				if (file_name[i] == '.' && file_name[i+1] == 'r' &&
					file_name[i+2] == 'a' && file_name[i+3] == 'w')
					break;
			}
			file_name[i+12] = file_name[i+6];
			file_name[i] = '_';
			file_name[i+1] = '5';
			file_name[i+2] = 'f';
			file_name[i+3] = '_';
			file_name[i+4] = 'a';
			file_name[i+5] = 'x';
			file_name[i+6] = 'i';
			file_name[i+7] = 's';
			file_name[i+8] = '.';
			file_name[i+9] = 't';
			file_name[i+10] = 'x';
			file_name[i+11] = 't';
			if ((fp=fopen(file_name, "wb"))==NULL){
				printf("write error...\n");
				exit(0); 
			};    
			(void)time(&t1); 
			GlobalSymmetry(xdim,ydim,zdim, dataset, critical_list, fivefold, fp);
			(void)time(&t2); 
			printf("time to find global symmetry: %d seconds. \n\n",(int)(t2-t1));
		}
		else if (type == 2 || type == 3) 
		{
			strcpy(file_name,filename.c_str());
			for(i = 0; i<MAX_STRING; i++) {
				if (file_name[i] == '.' && file_name[i+1] == 'r' &&
					file_name[i+2] == 'a' && file_name[i+3] == 'w')
					break;
			}
			j = i;
			while (file_name[i] != '_')
				i--;
			file_name[i+12] = file_name[j+6];
			file_name[i] = '_';
			file_name[i+1] = '5';
			file_name[i+2] = 'f';
			file_name[i+3] = '_';
			file_name[i+4] = 'a';
			file_name[i+5] = 'x';
			file_name[i+6] = 'i';
			file_name[i+7] = 's';
			file_name[i+8] = '.';
			file_name[i+9] = 't';
			file_name[i+10] = 'x';
			file_name[i+11] = 't';
			if ((fp=fopen(file_name, "r"))==NULL){
				printf("write error...\n");
				exit(0); 
			};    
			fivefold = (VECTOR*)malloc(sizeof(VECTOR)*12);
			for (i = 0; i < 12; i++) {
				fscanf(fp, "%f %f %f \n",&fx,&fy,&fz);
				fivefold[i].x = fx;
				fivefold[i].y = fy;
				fivefold[i].z = fz;
			}
		}
		
		if (type == 0) 
		{
			printf("begin segmentation of capsid layer: type = 0 ....\n");
			(void)time(&t1); 
			SimpleCapsidSegment(xdim,ydim,zdim,dataset,tlow,critical_list,fivefold);
			(void)time(&t2); 
			printf("time to segment capsid layer: %d seconds. \n\n",(int)(t2-t1));
			strcpy(file_name,filename.c_str());
			for(i = 0; i<MAX_STRING; i++) {
				if (file_name[i] == '.' && file_name[i+1] == 'r' &&
					file_name[i+2] == 'a' && file_name[i+3] == 'w')
					break;
			}
			file_name[i+13] = file_name[i+6];
			file_name[i] = '_';
			file_name[i+1] = 'c';
			file_name[i+2] = 'a';
			file_name[i+3] = 'p';
			file_name[i+4] = 's';
			file_name[i+5] = 'i';
			file_name[i+6] = 'd';
			file_name[i+7] = '.';
			file_name[i+8] = 'r';
			file_name[i+9] = 'a';
			file_name[i+10] = 'w';
			file_name[i+11] = 'i';
			file_name[i+12] = 'v';
			if ((fp=fopen(file_name, "wb"))==NULL){
				printf("write error...\n");
				exit(0); 
			};    
			write_rawiv_float(dataset, fp);
		}
		else if (type == 1) 
		{
			printf("begin segmentation of capsid layer: type = 1 ....\n");
			sx1 = int(params[3]);
			sy1 = int(params[4]);
			sz1 = int(params[5]);
			sx2 = int(params[6]);
			sy2 = int(params[7]);
			sz2 = int(params[8]);
			(void)time(&t1); 
			SingleCapsidSegment(xdim,ydim,zdim, dataset, tlow, 
								fivefold, sx1,sy1,sz1, sx2,sy2,sz2);
			(void)time(&t2); 
			printf("time to segment capsid layer: %d seconds. \n\n",(int)(t2-t1));
			strcpy(file_name,filename.c_str());
			for(i = 0; i<MAX_STRING; i++) {
				if (file_name[i] == '.' && file_name[i+1] == 'r' &&
					file_name[i+2] == 'a' && file_name[i+3] == 'w')
					break;
			}
			file_name[i+13] = file_name[i+6];
			file_name[i] = '_';
			file_name[i+1] = 'c';
			file_name[i+2] = 'a';
			file_name[i+3] = 'p';
			file_name[i+4] = 's';
			file_name[i+5] = 'i';
			file_name[i+6] = 'd';
			file_name[i+7] = '.';
			file_name[i+8] = 'r';
			file_name[i+9] = 'a';
			file_name[i+10] = 'w';
			file_name[i+11] = 'i';
			file_name[i+12] = 'v';
			if ((fp=fopen(file_name, "wb"))==NULL){
				printf("write error...\n");
				exit(0); 
			};    
			write_rawiv_float(dataset, fp);
		}
		else if (type == 2) 
		{
			printf("begin segmentation of capsid layer: type = 2 ....\n");
			small_radius = double(params[6]);  
			large_radius = double(params[7]);  
			printf("begin double layer segmentation ....\n");
			_result = (unsigned short *)malloc(sizeof(unsigned short)*xdim*ydim*zdim);
			(void)time(&t1); 
			DoubleCapsidSegment(xdim,ydim,zdim, tlow, dataset,_result,fivefold,
								small_radius, large_radius);
			(void)time(&t2); 
			printf("time to segment double layers: %d seconds. \n\n",(int)(t2-t1));
			
			sym_score = (float *)malloc(sizeof(float)*xdim*ydim*zdim);
			printf("begin writing outer layer ....\n");
			for (k=0; k<zdim; k++)
				for (j=0; j<ydim; j++)
					for (i=0; i<xdim; i++) 
					{
						sym_score[IndexVect(i,j,k)] = 0;
						if(_result[IndexVect(i,j,k)] == 1 &&
						   dataset[IndexVect(i,j,k)] > tlow)
							sym_score[IndexVect(i,j,k)] = dataset[IndexVect(i,j,k)];
					}
						strcpy(file_name,filename.c_str());
			for(i = 0; i<MAX_STRING; i++) {
				if (file_name[i] == '.' && file_name[i+1] == 'r' &&
					file_name[i+2] == 'a' && file_name[i+3] == 'w')
					break;
			}
			j = i;
			while (file_name[i] != '_')
				i--;
			file_name[i+12] = file_name[j+6];
			file_name[i] = '_';
			file_name[i+1] = 'o';
			file_name[i+2] = 'u';
			file_name[i+3] = 't';
			file_name[i+4] = 'e';
			file_name[i+5] = 'r';
			file_name[i+6] = '.';
			file_name[i+7] = 'r';
			file_name[i+8] = 'a';
			file_name[i+9] = 'w';
			file_name[i+10] = 'i';
			file_name[i+11] = 'v';
			if ((fp=fopen(file_name, "wb"))==NULL){
				printf("write error...\n");
				exit(0); 
			}; 
			write_rawiv_float(sym_score,fp);
			printf("begin writing inner layer ....\n");
			for (k=0; k<zdim; k++)
				for (j=0; j<ydim; j++)
					for (i=0; i<xdim; i++) 
					{
						sym_score[IndexVect(i,j,k)] = 0;
						if (_result[IndexVect(i,j,k)] == 2 &&
							dataset[IndexVect(i,j,k)] > tlow)
							sym_score[IndexVect(i,j,k)] = dataset[IndexVect(i,j,k)];
					}
						strcpy(file_name,filename.c_str());
			for(i = 0; i<MAX_STRING; i++) {
				if (file_name[i] == '.' && file_name[i+1] == 'r' &&
					file_name[i+2] == 'a' && file_name[i+3] == 'w')
					break;
			}
			j = i;
			while (file_name[i] != '_')
				i--;
			file_name[i+12] = file_name[j+6];
			file_name[i] = '_';
			file_name[i+1] = 'i';
			file_name[i+2] = 'n';
			file_name[i+3] = 'n';
			file_name[i+4] = 'e';
			file_name[i+5] = 'r';
			file_name[i+6] = '.';
			file_name[i+7] = 'r';
			file_name[i+8] = 'a';
			file_name[i+9] = 'w';
			file_name[i+10] = 'i';
			file_name[i+11] = 'v';
			if ((fp=fopen(file_name, "wb"))==NULL){
				printf("write error...\n");
				exit(0); 
			}; 
			write_rawiv_float(sym_score,fp);
		}
		else if (type == 3) 
		{
			printf("begin segmentation of capsid layer: type = 3 ....\n");
			numfold3 = int(params[3]);
			numfold5 = int(params[4]);
			numfold6 = int(params[5]);
			
			strcpy(file_name,filename.c_str());
			for(i = 0; i<MAX_STRING; i++) {
				if (file_name[i] == '.' && file_name[i+1] == 'r' &&
					file_name[i+2] == 'a' && file_name[i+3] == 'w')
					break;
			}
			j = i;
			while (file_name[i] != '_')
				i--;
			file_name[i+18] = file_name[j+6];
			file_name[i] = '_';
			file_name[i+1] = 'o';
			file_name[i+2] = 'u';
			file_name[i+3] = 't';
			file_name[i+4] = 'e';
			file_name[i+5] = 'r';
			file_name[i+6] = '_';
			if (numfold3 > 0)       
				file_name[i+7] = '3';
			else if (numfold6 > 0)
				file_name[i+7] = '6';
			file_name[i+8] = 'f';
			file_name[i+9] = '_';
			file_name[i+10] = 'a';
			file_name[i+11] = 'x';
			file_name[i+12] = 'i';
			file_name[i+13] = 's';
			file_name[i+14] = '.';
			file_name[i+15] = 't';
			file_name[i+16] = 'x';
			file_name[i+17] = 't';
			if ((fp=fopen(file_name, "r"))==NULL){
				printf("read error 111...\n");
				exit(0); 
			};    
			if (numfold3 > 0) {   
				threefold = (DB_VECTOR*)malloc(sizeof(DB_VECTOR)*6000); 
				i = 0;
				while ((fscanf(fp, "%f %f %f %f %f %f %f \n",
							   &fx,&fy,&fz,&gx,&gy,&gz,&score)) != EOF)
				{
					threefold[i].sx = fx;
					threefold[i].sy = fy;
					threefold[i].sz = fz;
					threefold[i].ex = gx;
					threefold[i].ey = gy;
					threefold[i].ez = gz;
					i++;
				}
				fclose(fp);
				numaxis3 = i/60;
				numaxis6 = 0;
			}
			else if (numfold6 > 0) {   
				sixfold = (DB_VECTOR*)malloc(sizeof(DB_VECTOR)*6000); 
				i = 0;
				while ((fscanf(fp, "%f %f %f %f %f %f %f \n",
							   &fx,&fy,&fz,&gx,&gy,&gz,&score)) != EOF)
				{
					sixfold[i].sx = fx;
					sixfold[i].sy = fy;
					sixfold[i].sz = fz;
					sixfold[i].ex = gx;
					sixfold[i].ey = gy;
					sixfold[i].ez = gz;
					i++;
				}
				fclose(fp);
				numaxis6 = i/60;
				numaxis3 = 0;
			}
			strcpy(file_name,filename.c_str());
			for(i = 0; i<MAX_STRING; i++) {
				if (file_name[i] == '.' && file_name[i+1] == 'r' &&
					file_name[i+2] == 'a' && file_name[i+3] == 'w')
					break;
			}
			j = i;
			while (file_name[i] != '_')
				i--;
			file_name[i+27] = file_name[j+6];
			file_name[i] = '_';
			file_name[i+1] = 'o';
			file_name[i+2] = 'u';
			file_name[i+3] = 't';
			file_name[i+4] = 'e';
			file_name[i+5] = 'r';
			file_name[i+6] = '_';
			file_name[i+7] = 's';
			file_name[i+8] = 'i';
			file_name[i+9] = 'm';
			file_name[i+10] = 'i';
			file_name[i+11] = 'l';
			file_name[i+12] = 'a';
			file_name[i+13] = 'r';
			file_name[i+14] = 'i';
			file_name[i+15] = 't';
			file_name[i+16] = 'y';
			file_name[i+17] = '_';
			file_name[i+18] = 's';
			file_name[i+19] = 'c';
			file_name[i+20] = 'o';
			file_name[i+21] = 'r';
			file_name[i+22] = 'e';
			file_name[i+23] = '.';
			file_name[i+24] = 't';
			file_name[i+25] = 'x';
			file_name[i+26] = 't';
			if ((fp=fopen(file_name, "r"))==NULL){
				printf("read error 222...\n");
				exit(0); 
			};    
			if (numfold3 > 0) {   
				coVarMatrix3 = (CVM*)malloc(sizeof(CVM)*numaxis3*numaxis3); 
				for (j = 0; j < numaxis3; j++) {
					for (i = 0; i < numaxis3; i++) {
						fscanf(fp, "%f %f %f %f \n",&score,&fx,&fy,&fz);
						coVarMatrix3[j*numaxis3+i].trans = fx;
						coVarMatrix3[j*numaxis3+i].rotat = fy;
						coVarMatrix3[j*numaxis3+i].angle = fz;
					}
					fscanf(fp, "\n");
				}
				fclose(fp);
			}
			else if (numfold6 > 0) {   
				coVarMatrix6 = (CVM*)malloc(sizeof(CVM)*numaxis6*numaxis6); 
				for (j = 0; j < numaxis6; j++) {
					for (i = 0; i < numaxis6; i++) {
						fscanf(fp, "%f %f %f %f \n",&score,&fx,&fy,&fz);
						coVarMatrix6[j*numaxis6+i].trans = fx;
						coVarMatrix6[j*numaxis6+i].rotat = fy;
						coVarMatrix6[j*numaxis6+i].angle = fz;
					}
					fscanf(fp, "\n");
				}
				fclose(fp);
			}
			
			small_radius = double(params[6]);  
			large_radius = double(params[7]);  
			printf("begin double capsid symmetric score calculations ....\n");  
			sym_score = (float *)malloc(sizeof(float)*xdim*ydim*zdim);
			(void)time(&t1); 
			if (numfold3 > 0)
				CapsidSegmentScore(xdim,ydim,zdim,dataset,sym_score,fivefold,threefold,coVarMatrix3,
								   numaxis3,numfold3,small_radius,large_radius);
			else if (numfold6 > 0)
				CapsidSegmentScore(xdim,ydim,zdim,dataset,sym_score,fivefold,sixfold,coVarMatrix6,
								   numaxis6,numfold6,small_radius,large_radius);
			(void)time(&t2); 
			printf("time to compute symmetric scores: %d seconds. \n\n",(int)(t2-t1));
			
			printf("begin double layer segmentation ....\n");
			_result = (unsigned short *)malloc(sizeof(unsigned short)*xdim*ydim*zdim);
			(void)time(&t1); 
			if (numfold3 > 0)
				CapsidSegmentMarch(xdim,ydim,zdim, tlow, dataset,sym_score,_result,fivefold,threefold,
								   coVarMatrix3,numaxis3,numfold3,small_radius, large_radius);
			else if (numfold6 > 0)
				CapsidSegmentMarch(xdim,ydim,zdim, tlow, dataset,sym_score,_result,fivefold,sixfold,
								   coVarMatrix6,numaxis6,numfold6,small_radius, large_radius);
			(void)time(&t2); 
			printf("time to segment double layers: %d seconds. \n\n",(int)(t2-t1));
			
			printf("begin writing outer layer ....\n");
			for (k=0; k<zdim; k++)
				for (j=0; j<ydim; j++)
					for (i=0; i<xdim; i++) {
						sym_score[IndexVect(i,j,k)] = 0;
						if (_result[IndexVect(i,j,k)] == 1 &&
							dataset[IndexVect(i,j,k)] > tlow)
							sym_score[IndexVect(i,j,k)] = dataset[IndexVect(i,j,k)];
					}
						strcpy(file_name,filename.c_str());
			for(i = 0; i<MAX_STRING; i++) {
				if (file_name[i] == '.' && file_name[i+1] == 'r' &&
					file_name[i+2] == 'a' && file_name[i+3] == 'w')
					break;
			}
			j = i;
			while (file_name[i] != '_')
				i--;
			file_name[i+12] = file_name[j+6];
			file_name[i] = '_';
			file_name[i+1] = 'o';
			file_name[i+2] = 'u';
			file_name[i+3] = 't';
			file_name[i+4] = 'e';
			file_name[i+5] = 'r';
			file_name[i+6] = '.';
			file_name[i+7] = 'r';
			file_name[i+8] = 'a';
			file_name[i+9] = 'w';
			file_name[i+10] = 'i';
			file_name[i+11] = 'v';
			if ((fp=fopen(file_name, "wb"))==NULL){
				printf("write error...\n");
				exit(0); 
			}; 
			write_rawiv_float(sym_score,fp);
			printf("begin writing inner layer ....\n");
			for (k=0; k<zdim; k++)
				for (j=0; j<ydim; j++)
					for (i=0; i<xdim; i++) {
						sym_score[IndexVect(i,j,k)] = 0;
						if (_result[IndexVect(i,j,k)] == 2 &&
							dataset[IndexVect(i,j,k)] > tlow)
							sym_score[IndexVect(i,j,k)] = dataset[IndexVect(i,j,k)];
					}
						strcpy(file_name,filename.c_str());
			for(i = 0; i<MAX_STRING; i++) {
				if (file_name[i] == '.' && file_name[i+1] == 'r' &&
					file_name[i+2] == 'a' && file_name[i+3] == 'w')
					break;
			}
			j = i;
			while (file_name[i] != '_')
				i--;
			file_name[i+12] = file_name[j+6];
			file_name[i] = '_';
			file_name[i+1] = 'i';
			file_name[i+2] = 'n';
			file_name[i+3] = 'n';
			file_name[i+4] = 'e';
			file_name[i+5] = 'r';
			file_name[i+6] = '.';
			file_name[i+7] = 'r';
			file_name[i+8] = 'a';
			file_name[i+9] = 'w';
			file_name[i+10] = 'i';
			file_name[i+11] = 'v';
			if ((fp=fopen(file_name, "wb"))==NULL){
				printf("write error...\n");
				exit(0); 
			}; 
			write_rawiv_float(sym_score,fp);
		}
		
		result = "done";
	}
	
	std::string help() { return std::string("Segment virus capsid layer."); }
} SegCapsid(&s);

int main(int argc, char **argv)
{
	if(argc != 2)
	{
		std::cerr<< "Usage: " << argv[0] << " serverPort\n";
		return -1;
	}
	
	int port = atoi(argv[1]);
	XmlRpc::setVerbosity(5);
	s.bindAndListen(port);
	s.enableIntrospection(true);
	s.work(-1.0);
	
	return 0;
}
