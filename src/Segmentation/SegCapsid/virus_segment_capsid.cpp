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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <xmlrpc/XmlRpc.h>
#include <cmath>
#include <Segmentation/SegCapsid/segcapsid.h>

#define _LITTLE_ENDIAN 1

#define IndexVect(i,j,k) ((k)*xdim*ydim + (j)*xdim + (i))
#define max2(x, y) ((x>y) ? (x):(y))
#define min2(x, y) ((x<y) ? (x):(y))
#define PolyNum       15
#define PIE           3.1415926536f
#define MAX_STRING    256

using namespace XmlRpc;
using namespace SegCapsid;

static int v_num, t_num;
static VECTOR *vertex;
static INTVECT *triangle;
void print_5_fold_axis_raw(int xdim, int ydim, int zdim, VECTOR *fivefold, FILE *fp, FILE *fp2, float *span_tmp, float *orig_tmp);
void DrawLineCapsid(float sx, float sy, float sz, float ex, float ey, float ez, float radius);
void generate_5_fold_transform(FILE *fp2, FILE *fp);

int virusSegCapsid(XmlRpcValue &params, XmlRpcValue &result)
{
	int xdim,ydim,zdim;
	float *dataset;
	FILE *fp,*fp2;
	char file_name[MAX_STRING],file_nameA[MAX_STRING];
	VECTOR fivefold[12];
	CPNT *critical_list;
	float span_tmp[3], orig_tmp[3];
	float tlow,score;
	int i,j,k,type;
	float fx,fy,fz;
	float gx,gy,gz;
	time_t t1,t2;
	int sx1,sy1,sz1;
	int sx2,sy2,sz2;
	float *sym_score;
	unsigned short *_result;
	DB_VECTOR *sixfold=NULL,*threefold=NULL;
	float small_radius,large_radius;
	int numfold6,numfold5,numfold3;
	int numaxis6=0,numaxis3=0;
	CVM *coVarMatrix3=NULL,*coVarMatrix6=NULL;
	std::string filename;
	
	int fscanf_return = 0;
	
	filename = std::string(params[0]);
	type = int(params[1]);
	tlow = double(params[2]);
	
	printf("begin reading rawiv.... \n");
	//span_tmp = (float *)malloc(sizeof(float)*3);
	//orig_tmp = (float *)malloc(sizeof(float)*3);
	(void)time(&t1); 
	read_data(&xdim,&ydim,&zdim,&dataset,span_tmp,orig_tmp,const_cast<char *>(filename.c_str()),0,1);
	printf("xdim: %d  ydim: %d  zdim: %d \n",xdim,ydim,zdim);
	(void)time(&t2); 
	printf("time to read dataset: %d seconds. \n\n",(int)(t2-t1));
	
	if (type == 0 || type == 1) 
	{
	  if(type == 0 ? bool(params[6]) : bool(params[9]))
	    {
	      printf("begin diffusion .... \n");
	      (void)time(&t1); 
	      Diffuse(xdim,ydim,zdim,dataset);
	      (void)time(&t2); 
	      printf("time to diffuse dataset: %d seconds. \n\n",(int)(t2-t1));
	    }
	  else
	    printf("not running diffusion... \n\n");
	  
		
		printf("begin finding critical points....\n");
		(void)time(&t1); 
		FindCriticalPoints(xdim,ydim,zdim,dataset,&critical_list,tlow,3,2);
		(void)time(&t2); 
		printf("time to find critical points: %d seconds. \n\n",(int)(t2-t1));
		
		//fivefold = (VECTOR*)malloc(sizeof(VECTOR)*12);
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
			return 0;
		};    
		(void)time(&t1); 
		GlobalSymmetry(xdim,ydim,zdim, dataset, critical_list, fivefold, fp);
		(void)time(&t2); 
		printf("time to find global symmetry: %d seconds. \n\n",(int)(t2-t1));
//-Jesse-Adding code here
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
		file_name[i+9] = 'r';
		file_name[i+10] = 'a';
		file_name[i+11] = 'w';
		if ((fp=fopen(file_name, "wb"))==NULL){
			printf("write error...\n");
			return 0;
		};

		strcpy(file_nameA,filename.c_str());
		for(i = 0; i<MAX_STRING; i++) {
			if (file_nameA[i] == '.' && file_nameA[i+1] == 'r' &&
				file_nameA[i+2] == 'a' && file_nameA[i+3] == 'w')
				break;
		}
		file_nameA[i+17] = file_nameA[i+6];
		file_nameA[i] = '_';
		file_nameA[i+1] = '5';
		file_nameA[i+2] = 'f';
		file_nameA[i+3] = '_';
		file_nameA[i+4] = 'r';
		file_nameA[i+5] = 'e';
		file_nameA[i+6] = 'a';
		file_nameA[i+7] = 'l';
		file_nameA[i+8] = '_';
		file_nameA[i+9] = 'a';
		file_nameA[i+10] = 'x';
		file_nameA[i+11] = 'i';
		file_nameA[i+12] = 's';
		file_nameA[i+13] = '.';
		file_nameA[i+14] = 't';
		file_nameA[i+15] = 'x';
		file_nameA[i+16] = 't';
		if ((fp2=fopen(file_nameA, "wb"))==NULL){
			printf("write error...\n");
			return 0;
		};
		//For output the real axis, zq added FILE fp2.
		print_5_fold_axis_raw(xdim, ydim, zdim, fivefold, fp, fp2, span_tmp, orig_tmp);

		//Generate transformation matrix for 12 axes. 

		strcpy(file_name,filename.c_str());
		for(i = 0; i<MAX_STRING; i++) {
			if (file_name[i] == '.' && file_name[i+1] == 'r' &&
				file_name[i+2] == 'a' && file_name[i+3] == 'w')
				break;
		}
		file_name[i+14] = file_name[i+6];
		file_name[i] = '_';
		file_name[i+1] = '5';
		file_name[i+2] = 'f';
		file_name[i+3] = '_';
		file_name[i+4] = 'm';
		file_name[i+5] = 'a';
		file_name[i+6] = 't';
		file_name[i+7] = 'r';
		file_name[i+8] = 'i';
		file_name[i+9] = 'x';
		file_name[i+10] = '.';
		file_name[i+11] = 't';
		file_name[i+12] = 'x';
		file_name[i+13] = 't';
		if ((fp=fopen(file_name, "wb"))==NULL){
			printf("write error...\n");
			return 0;
		};

		if((fp2= fopen(file_nameA,"r"))==NULL){
			printf("read error of axis file ... \n");
			return 0;
		};

		generate_5_fold_transform(fp2, fp);

//-Jesse done adding code
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
			printf("read error...\n");
			return 0;
		};    
		//fivefold = (VECTOR*)malloc(sizeof(VECTOR)*12);
		for (i = 0; i < 12; i++) {
			fscanf_return = fscanf(fp, "%f %f %f \n",&fx,&fy,&fz);
			fivefold[i].x = fx;
			fivefold[i].y = fy;
			fivefold[i].z = fz;
		}
	}
	
	if (type == 0) 
	{
		printf("begin segmentation of capsid layer: type = 0 ....\n");
		sx1 = int(params[3]);
		sy1 = int(params[4]);
		sz1 = int(params[5]);
		(void)time(&t1); 
		SimpleCapsidSegment(xdim,ydim,zdim,dataset,tlow,sx1,sy1,sz1,fivefold);
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
			return 0; 
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
			return 0;
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
			return 0; 
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
			return 0;
		}; 
		write_rawiv_float(sym_score,fp);

		free(_result);
		free(sym_score);
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
			return 0; 
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
			return 0; 
		};    
		if (numfold3 > 0) {   
			coVarMatrix3 = (CVM*)malloc(sizeof(CVM)*numaxis3*numaxis3); 
			for (j = 0; j < numaxis3; j++) {
				for (i = 0; i < numaxis3; i++) {
					fscanf_return = fscanf(fp, "%f %f %f %f \n",&score,&fx,&fy,&fz);
					coVarMatrix3[j*numaxis3+i].trans = fx;
					coVarMatrix3[j*numaxis3+i].rotat = fy;
					coVarMatrix3[j*numaxis3+i].angle = fz;
				}
				fscanf_return = fscanf(fp, "\n");
			}
			fclose(fp);
		}
		else if (numfold6 > 0) {   
			coVarMatrix6 = (CVM*)malloc(sizeof(CVM)*numaxis6*numaxis6); 
			for (j = 0; j < numaxis6; j++) {
				for (i = 0; i < numaxis6; i++) {
					fscanf_return = fscanf(fp, "%f %f %f %f \n",&score,&fx,&fy,&fz);
					coVarMatrix6[j*numaxis6+i].trans = fx;
					coVarMatrix6[j*numaxis6+i].rotat = fy;
					coVarMatrix6[j*numaxis6+i].angle = fz;
				}
				fscanf_return = fscanf(fp, "\n");
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
			return 0; 
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
			return 0;
		}; 
		write_rawiv_float(sym_score,fp);

		if(numfold3 > 0)
		  {
		    free(threefold);
		    free(coVarMatrix3);
		  }
		else if(numfold6 > 0)
		  {
		    free(sixfold);
		    free(coVarMatrix6);
		  }

		free(sym_score);
		free(_result);
	}


	result = XmlRpcValue(true);
	return 1;
}

void print_5_fold_axis_raw(int xdim, int ydim, int zdim, VECTOR *fivefold, FILE *fp, FILE *fp2, float *span_tmp, float *orig_tmp){
	vertex = (VECTOR*)malloc(sizeof(VECTOR)*5000*PolyNum);
	triangle = (INTVECT*)malloc(sizeof(INTVECT)*5000*PolyNum);
	v_num = 0;
	t_num = 0;
	int n;

	for (n=0; n<12; n++) {
		DrawLineCapsid(fivefold[n].x, fivefold[n].y, fivefold[n].z,(float)(xdim/2),(float)(ydim/2),(float)(zdim/2),3.f);    
	}
	fprintf(fp, "%d %d\n",v_num,t_num);
	for (n = 0; n < v_num; n++){
		fprintf(fp, "%f %f %f\n",vertex[n].x*span_tmp[0]+orig_tmp[0], 
			vertex[n].y*span_tmp[1]+orig_tmp[1],vertex[n].z*span_tmp[2]+orig_tmp[2]);
	}
	for (n = 0; n < t_num; n++){
		fprintf(fp, "%d %d %d\n",triangle[n].x,triangle[n].y,triangle[n].z);
	}
	fclose(fp);
	
	float temp[3];
	for(int i=0; i<3; i++)
		temp[i]=0;

	for(n =1; n <v_num; n = n+2)
	{
	   temp[0] += vertex[n].x*span_tmp[0]+orig_tmp[0];
	   temp[1] += vertex[n].y*span_tmp[1]+orig_tmp[1];
	   temp[2] += vertex[n].z*span_tmp[2]+orig_tmp[2];

	   if((n+1)%30 ==0)
	   {
	   		fprintf(fp2, "%f %f %f\n", temp[0]/15.0, temp[1]/15.0, temp[2]/15.0);
			for(int i=0; i<3; i++)
			temp[i] =0;
	   }

	}
	fclose(fp2);

	free(triangle);
	free(vertex);

}

template <class T>
void MatrixMutipleVector(T* A, T* B, T* v, T* u)
{
    for(int l=0; l<3; l++)
    {
    u[l] = A[3*l]*v[0]+A[3*l+1]*v[1]+A[3*l+2]*v[2]+B[l];
    }
}

template <class T>
void MatrixMultiplyMatrix(T* A, T* B, T* C)
{
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			C[i][j] = 0.0;
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			for(int l = 0; l < 3; l++)
				C[i][j] += A[i][l] * B[l][j];
	
}


#define EPSILON 0.000001

#define CROSS(dest, v1, v2){                 \
          dest[0] = v1[1] * v2[2] - v1[2] * v2[1]; \
          dest[1] = v1[2] * v2[0] - v1[0] * v2[2]; \
          dest[2] = v1[0] * v2[1] - v1[1] * v2[0];}

#define DOT(v1, v2) (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2])

#define SUB(dest, v1, v2){       \
          dest[0] = v1[0] - v2[0]; \
          dest[1] = v1[1] - v2[1]; \
          dest[2] = v1[2] - v2[2];}

/*
 * A function for creating a rotation matrix that rotates a vector called
 * "from" into another vector called "to".
 * Input : from[3], to[3] which both must be *normalized* non-zero vectors
 * Output: mtx[3][3] -- a 3x3 matrix in colum-major form
 * Authors: Tomas Möller, John Hughes
 *          "Efficiently Building a Matrix to Rotate One Vector to Another"
 *          Journal of Graphics Tools, 4(4):1-4, 1999
 */
template <class T>
void fromToRotation(T from[3], T to[3], T mtx[3][3]) {
  T v[3];
  T e, h, f;

  CROSS(v, from, to);
  e = DOT(from, to);
  f = (e < 0)? -e:e;
  if (f > 1.0 - EPSILON)     /* "from" and "to"-vector almost parallel */
  {
    T u[3], v[3]; /* temporary storage vectors */
    T x[3];       /* vector most nearly orthogonal to "from" */
    T c1, c2, c3; /* coefficients for later use */
    int i, j;

    x[0] = (from[0] > 0.0)? from[0] : -from[0];
    x[1] = (from[1] > 0.0)? from[1] : -from[1];
    x[2] = (from[2] > 0.0)? from[2] : -from[2];

    if (x[0] < x[1])
    {
      if (x[0] < x[2])
      {
        x[0] = 1.0; x[1] = x[2] = 0.0;
      }
      else
      {
        x[2] = 1.0; x[0] = x[1] = 0.0;
      }
    }
    else
    {
      if (x[1] < x[2])
      {
        x[1] = 1.0; x[0] = x[2] = 0.0;
      }
      else
      {
        x[2] = 1.0; x[0] = x[1] = 0.0;
      }
    }

    u[0] = x[0] - from[0]; u[1] = x[1] - from[1]; u[2] = x[2] - from[2];
    v[0] = x[0] - to[0];   v[1] = x[1] - to[1];   v[2] = x[2] - to[2];

    c1 = 2.0 / DOT(u, u);
    c2 = 2.0 / DOT(v, v);
    c3 = c1 * c2  * DOT(u, v);

    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        mtx[i][j] =  - c1 * u[i] * u[j]
                     - c2 * v[i] * v[j]
                     + c3 * v[i] * u[j];
      }
      mtx[i][i] += 1.0;
    }
  }
  else  /* the most common case, unless "from"="to", or "from"=-"to" */
  {
#if 0
    /* unoptimized version - a good compiler will optimize this. */
    /* h = (1.0 - e)/DOT(v, v); old code */
    h = 1.0/(1.0 + e);      /* optimization by Gottfried Chen */
    mtx[0][0] = e + h * v[0] * v[0];
    mtx[0][1] = h * v[0] * v[1] - v[2];
    mtx[0][2] = h * v[0] * v[2] + v[1];

    mtx[1][0] = h * v[0] * v[1] + v[2];
    mtx[1][1] = e + h * v[1] * v[1];
    mtx[1][2] = h * v[1] * v[2] - v[0];

    mtx[2][0] = h * v[0] * v[2] - v[1];
    mtx[2][1] = h * v[1] * v[2] + v[0];
    mtx[2][2] = e + h * v[2] * v[2];
#else
    /* ...otherwise use this hand optimized version (9 mults less) */
    T hvx, hvz, hvxy, hvxz, hvyz;
    /* h = (1.0 - e)/DOT(v, v); old code */
    h = 1.0/(1.0 + e);      /* optimization by Gottfried Chen */
    hvx = h * v[0];
    hvz = h * v[2];
    hvxy = hvx * v[1];
    hvxz = hvx * v[2];
    hvyz = hvz * v[1];
    mtx[0][0] = e + hvx * v[0];
    mtx[0][1] = hvxy - v[2];
    mtx[0][2] = hvxz + v[1];

    mtx[1][0] = hvxy + v[2];
    mtx[1][1] = e + h * v[1] * v[1];
    mtx[1][2] = hvyz - v[0];

    mtx[2][0] = hvxz - v[1];
    mtx[2][1] = hvyz + v[0];
    mtx[2][2] = e + hvz * v[2];
#endif
  }
}

template <class T>
T InnerProduct(T u[3], T v[3])
{
	return (u[0]*v[0]+u[1]*v[1] + u[2]*v[2]);
}

template <class T>
T len(T u[3])
{
	return sqrt(InnerProduct(u,u));
}	

template <class T>
void normalize(T* u)
{	
	T leng = len(u);
	for(int i=0; i<3; i++)
	u[i] = u[i]/leng; 
}

// The following rotation matrix can be found in "Geometric Modeling
// and Quantative Visualization of Virus Ultrastructure" by Bajaj.
// The commented lines are for right hand side rotation.

template <class T>
void defineRotationMatrix(T matrix[3][3], T alpha, T* v)
{
	T c, s, t;
	c = cos(alpha);
	s = sin(alpha);
	t = 1-c;
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)
		matrix[i][j] = 0.0;
	matrix[0][0] = v[0]*v[0]*t + c;
	matrix[0][1] = v[0]*v[1]*t - v[2]*s;
	//matrix[0][1] = v[0]*v[1]*t + v[2]*s;
	matrix[0][2] = v[0]*v[2]*t + v[1]*s;
	//matrix[0][2] = v[0]*v[2]*t - v[1]*s;
	matrix[1][0] = v[0]*v[1]*t + v[2]*s;
	//matrix[1][0] = v[0]*v[1]*t - v[2]*s;
	matrix[1][1] = v[1]*v[1]*t + c;
	matrix[1][2] = v[1]*v[2]*t - v[0]*s;
	//matrix[1][2] = v[1]*v[2]*t + v[0]*s;
	matrix[2][0] = v[0]*v[2]*t - v[1]*s;
	//matrix[2][0] = v[0]*v[2]*t + v[1]*s;
	matrix[2][1] = v[1]*v[2]*t + v[0]*s;
	//matrix[2][1] = v[1]*v[2]*t - v[0]*s;
	matrix[2][2] = v[2]*v[2]*t + c;
}


void generate_5_fold_transform(FILE *fp2, FILE *fp)
{
	float center[3];
	float axis[12][3];
	float x, y, z;
	for(int i = 0; i < 3; i++)
	   center[i] = 0.0;

	for(int i = 0;i < 12; i++)
	{
		fscanf(fp2, "%f %f %f \n", &x, &y,&z);
		axis[i][0] = x;
		axis[i][1] = y;
		axis[i][2] = z;
		for(int j=0; j<3; j++)
			center[j]+= axis[i][j];
	}
	
	for(int j = 0; j < 3; j++)
		center[j] /= 12.0;
	
	for(int i = 0; i < 12; i++)
		for(int j = 0; j < 3; j++)
			axis[i][j] -= center[j];

	for(int i = 0; i < 12; i++)
		normalize(axis[i]);

	float mat[3][3], matrix[3][3], matrixTmp[3][3];
	float Matrix[9], Trans[3], MatrixN[9];
	float NewTrans[3];

	//From Zeyun's segmentation, the first axis is corresponding to the segmented 5-fold subunit.

	for(int i = 0; i < 12; i++)
	{
		if(InnerProduct(axis[0], axis[i])> -0.5 && InnerProduct(axis[0], axis[i]) <0.5)
			defineRotationMatrix(matrix, (float)PIE, axis[0]);
		else if (InnerProduct(axis[0],axis[i]) > 0.5)
			defineRotationMatrix(matrix, (float)0.0, axis[0]); 
		else 
			defineRotationMatrix(matrix, (float)(PIE*1.0/2.0), axis[0]);



		fromToRotation(axis[0], axis[i], mat);

		MatrixMultiplyMatrix(mat, matrix, matrixTmp);

		fprintf(fp, "%f %f %f \n", matrixTmp[0][0], matrixTmp[0][1], matrixTmp[0][2]);
		fprintf(fp, "%f %f %f \n", matrixTmp[1][0], matrixTmp[1][1], matrixTmp[1][2]);
		fprintf(fp, "%f %f %f \n", matrixTmp[2][0], matrixTmp[2][1], matrixTmp[2][2]);


		Matrix[0] = 1.0 - matrixTmp[0][0];
    	Matrix[1] =  -1.0 * matrixTmp[0][1];
    	Matrix[2] = -1.0 * matrixTmp[0][2];
		Matrix[3] = -1.0 * matrixTmp[1][0];
    	Matrix[4] = 1.0 - matrixTmp[1][1];
    	Matrix[5] = -1.0 * matrixTmp[1][2];
    	Matrix[6] = -1.0 * matrixTmp[2][0];
    	Matrix[7] =  -1.0 * matrixTmp[2][1];
    	Matrix[8] = 1.0 - matrixTmp[2][2];

		Trans[0] = 0; 
		Trans[1] = 0; 
		Trans[2] = 0; 

		MatrixMutipleVector(Matrix, Trans, center, NewTrans);

		fprintf(fp, "%f %f %f \n\n", NewTrans[0], NewTrans[1], NewTrans[2]);
	}


	
	fclose(fp);
	fclose(fp2);	

}





void DrawLineCapsid(float sx, float sy, float sz, float ex, float ey, float ez, float radius){
  float x,y,z;
  float xx,yy,zz;
  float xxx,yyy,zzz;
  float a[3][3],b[3][3];
  float theta, phi;
  int m;
  

  theta = (float)atan2(sy-ey,sx-ex);
  phi = (float)atan2(sz-ez, sqrt((sx-ex)*(sx-ex)+(sy-ey)*(sy-ey)));
  
  a[0][0] = (float)(cos(0.5*PIE-phi)*cos(theta));
  a[0][1] = (float)(cos(0.5*PIE-phi)*sin(theta));
  a[0][2] = (float)-sin(0.5*PIE-phi);
  a[1][0] = (float)-sin(theta);
  a[1][1] = (float)cos(theta);
  a[1][2] = 0;
  a[2][0] = (float)(sin(0.5*PIE-phi)*cos(theta));
  a[2][1] = (float)(sin(0.5*PIE-phi)*sin(theta));
  a[2][2] = (float)cos(0.5*PIE-phi);

  b[0][0] = (float)(cos(0.5*PIE-phi)*cos(theta));
  b[0][1] = (float)-sin(theta); 
  b[0][2] = (float)(sin(0.5*PIE-phi)*cos(theta)); 
  b[1][0] = (float)(cos(0.5*PIE-phi)*sin(theta));
  b[1][1] = (float)cos(theta);
  b[1][2] = (float)(sin(0.5*PIE-phi)*sin(theta));
  b[2][0] = (float)-sin(0.5*PIE-phi);
  b[2][1] = 0;
  b[2][2] = (float)cos(0.5*PIE-phi);

  
  xx = (float)sqrt((sy-ey)*(sy-ey)+(sx-ex)*(sx-ex));
  if (xx == 0) {
    x = radius+ex;
    y = ey;
    z = ez;
  }
  else {
    x = radius*(ey-sy)/xx+ex;
    y = radius*(sx-ex)/xx+ey;
    z = ez;
  }
  
  vertex[v_num].x = x;
  vertex[v_num].y = y;
  vertex[v_num].z = z;
  vertex[v_num+1].x = x+sx-ex;
  vertex[v_num+1].y = y+sy-ey;
  vertex[v_num+1].z = z+sz-ez;
    
  
  x = x-ex;
  y = y-ey;
  z = z-ez;
  
  xx = a[0][0]*x+a[0][1]*y+a[0][2]*z;
  yy = a[1][0]*x+a[1][1]*y+a[1][2]*z;
  zz = a[2][0]*x+a[2][1]*y+a[2][2]*z;
  
  for (m = 1; m < PolyNum; m++) {
    x = (float)(cos(2*PIE*(float)(m)/(float)(PolyNum))*xx - 
      sin(2*PIE*(float)(m)/(float)(PolyNum))*yy);
    y = (float)(sin(2*PIE*(float)(m)/(float)(PolyNum))*xx + 
      cos(2*PIE*(float)(m)/(float)(PolyNum))*yy);
    z = zz;
    
    xxx = b[0][0]*x+b[0][1]*y+b[0][2]*z+ex;
    yyy = b[1][0]*x+b[1][1]*y+b[1][2]*z+ey;
    zzz = b[2][0]*x+b[2][1]*y+b[2][2]*z+ez;
    
    vertex[v_num+2*m].x = xxx;
    vertex[v_num+2*m].y = yyy;
    vertex[v_num+2*m].z = zzz;
    vertex[v_num+2*m+1].x = xxx+sx-ex;
    vertex[v_num+2*m+1].y = yyy+sy-ey;
    vertex[v_num+2*m+1].z = zzz+sz-ez;
    
  }


  for (m = 0; m < PolyNum-1; m++) {
    triangle[t_num+2*m].x = v_num+2*m;
    triangle[t_num+2*m].y = v_num+2*m+1;
    triangle[t_num+2*m].z = v_num+2*m+2;
    triangle[t_num+2*m+1].x = v_num+2*m+1;
    triangle[t_num+2*m+1].y = v_num+2*m+2;
    triangle[t_num+2*m+1].z = v_num+2*m+3;
  }

  triangle[t_num+2*PolyNum-2].x = v_num+2*PolyNum-2;
  triangle[t_num+2*PolyNum-2].y = v_num+2*PolyNum-1;
  triangle[t_num+2*PolyNum-2].z = v_num;
  triangle[t_num+2*PolyNum-1].x = v_num+2*PolyNum-1;
  triangle[t_num+2*PolyNum-1].y = v_num;
  triangle[t_num+2*PolyNum-1].z = v_num+1;

  v_num += 2*PolyNum;
  t_num += 2*PolyNum;

}
