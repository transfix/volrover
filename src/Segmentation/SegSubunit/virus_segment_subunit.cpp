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
#include <math.h>
#include <Segmentation/SegSubunit/segsubunit.h>

#define _LITTLE_ENDIAN 1

#define IndexVect(i,j,k) ((k)*xdim*ydim + (j)*xdim + (i))
#define max2(x, y) ((x>y) ? (x):(y))
#define min2(x, y) ((x<y) ? (x):(y))
#define PolyNum       15
#define PIE           3.1415926536f
#define MAX_STRING    256

using namespace XmlRpc;
using namespace SegSubunit;

static int v_num, t_num;
static VECTOR *vertex;
static INTVECT *triangle;
void print_local_axis_raw(int xdim, int ydim, int zdim, DB_VECTOR *sym_axes, int num_axes, FILE *fp, float *span_tmp, float *orig_tmp);
void DrawLine_subunit(float sx, float sy, float sz, float ex, float ey, float ez, float radius);

int virusSegSubunit(XmlRpcValue &params, XmlRpcValue &result){
	using namespace SegSubunit;
	int xdim,ydim,zdim;
	float *dataset;
	unsigned short *_result;
	FILE *fp,*fp1;
	char file_name[MAX_STRING];
	VECTOR *fivefold;
	DB_VECTOR *sixfold=NULL,*threefold=NULL;
	CPNT *critical_list;
	float *span_tmp, *orig_tmp, *max_tmp;
	float tlow;
	int numfold6,numfold5,numfold3;
	int numaxis6,numaxis3;
	int h_num, k_num;
	int i,j;
	float fx,fy,fz;
	CVM *coVarMatrix3,*coVarMatrix6;
	int init_radius;
	std::string filename;
	
	vertex = (VECTOR*)malloc(sizeof(VECTOR)*5000*PolyNum);
	triangle = (INTVECT*)malloc(sizeof(INTVECT)*5000*PolyNum);
	int fscanf_return = 0;
	
	time_t t1,t2;
	
	filename = std::string(params[0]);
	
	printf("begin reading rawiv.... \n");
	span_tmp = (float *)malloc(sizeof(float)*3);
	orig_tmp = (float *)malloc(sizeof(float)*3);
	max_tmp = (float *)malloc(sizeof(float)*3);
	(void)time(&t1); 
	read_data(&xdim,&ydim,&zdim,&dataset,span_tmp,orig_tmp,const_cast<char *>(filename.c_str()),0,1);
	printf("xdim: %d  ydim: %d  zdim: %d \n",xdim,ydim,zdim);
	printf("span: %f, %f, %f\n", span_tmp[0], span_tmp[1], span_tmp[2]);
	printf("orig: %f, %f, %f\n", orig_tmp[0], orig_tmp[1], orig_tmp[2]); 
	max_tmp[0] = orig_tmp[0]+(xdim-1)*span_tmp[0];
	max_tmp[1] = orig_tmp[1]+(ydim-1)*span_tmp[1];
	max_tmp[2] = orig_tmp[2]+(zdim-1)*span_tmp[2];

	(void)time(&t2); 
	printf("time to read dataset: %d seconds. \n\n",(int)(t2-t1));
	
	
	tlow = 50.0f;  
	h_num = int(params[1]);
	k_num = int(params[2]);
	printf("begin finding critical points....\n");
	(void)time(&t1); 
	FindCriticalPoints(xdim,ydim,zdim,dataset,&critical_list,tlow,h_num,k_num);
	(void)time(&t2); 
	printf("time to find critical points: %d seconds. \n\n",(int)(t2-t1));
	
	
	fivefold = (VECTOR*)malloc(sizeof(VECTOR)*12);
	printf("reading global symmetry from file ....\n");
	
	strcpy(file_name,filename.c_str());
	for(i = 0; i<MAX_STRING; i++) {
		if (file_name[i] == '.' && file_name[i+1] == 'r' &&	file_name[i+2] == 'a' && file_name[i+3] == 'w'){
			break;
		}
	}
	j = i;
	while (file_name[i] != '_'){
		i--;
	}
	file_name[i+17] = file_name[j+6];
	file_name[i] = '_';
	file_name[i+1] = '5';
	file_name[i+2] = 'f';
	file_name[i+3] = '_';
	file_name[i+4] = 'r';
	file_name[i+5] = 'e';
	file_name[i+6] = 'a';
	file_name[i+7] = 'l';
	file_name[i+8] = '_';
	file_name[i+9] = 'a';
	file_name[i+10] = 'x';
	file_name[i+11] = 'i';
	file_name[i+12] = 's';
	file_name[i+13] = '.';
	file_name[i+14] = 't';
	file_name[i+15] = 'x';
	file_name[i+16] = 't';
	if ((fp=fopen(file_name, "r"))==NULL){
		printf("read error...\n");
		return 0; 
	}
	for(i=0; i<12; i++) {
		fscanf_return = fscanf(fp, "%f %f %f \n",&fx,&fy,&fz);
		fivefold[i].x = fx;
		fivefold[i].y = fy;
		fivefold[i].z = fz;
	}
	fclose(fp);
	
	numfold3 = int(params[3]);
	numfold5 = int(params[4]);
	numfold6 = int(params[5]);
	init_radius = int(params[6]);
	
	if ((numfold3 > 0 && numfold6 == 0) || (numfold6 > 0 && numfold3 == 0)) {
		
		numaxis3 = numfold3;
		numaxis6 = numfold6;
		if (numfold3 > 0){
			threefold = (DB_VECTOR*)malloc(sizeof(DB_VECTOR)*6000); 
		}
		if (numfold6 > 0){
			sixfold = (DB_VECTOR*)malloc(sizeof(DB_VECTOR)*6000); 
		}
		
		printf("begin initial local symmetry detection ....\n");
		(void)time(&t1); 
		LocalSymmetry(xdim,ydim,zdim,orig_tmp, span_tmp, dataset, critical_list,h_num,k_num,threefold,fivefold,sixfold, &numaxis3,&numaxis6,numfold5);
		(void)time(&t2); 
		printf("time to detect initial local symmetry: %d seconds. \n\n",(int)(t2-t1));
		
		//-Jesse-Adding code here
		strcpy(file_name,filename.c_str());
		for(i = 0; i<MAX_STRING; i++) {
			if (file_name[i] == '.' && file_name[i+1] == 'r' &&
				file_name[i+2] == 'a' && file_name[i+3] == 'w')
				break;
		}
		file_name[i+23] = file_name[i+6];
		file_name[i] = '_';
		file_name[i+1] = '3';
		file_name[i+2] = 'f';
		file_name[i+3] = '_';
		file_name[i+4] = 'a';
		file_name[i+5] = 'x';
		file_name[i+6] = 'e';
		file_name[i+7] = 's';
		file_name[i+8] = '_';
		file_name[i+9] = 'p';
		file_name[i+10] = 'r';
		file_name[i+11] = 'e';
		file_name[i+12] = '_';
		file_name[i+13] = 'r';
		file_name[i+14] = 'e';
		file_name[i+15] = 'f';
		file_name[i+16] = 'i';
		file_name[i+17] = 'n';
		file_name[i+18] = 'e';
		file_name[i+19] = '.';
		file_name[i+20] = 'r';
		file_name[i+21] = 'a';
		file_name[i+22] = 'w';
		
		if(numfold3>0)
		{
			if ((fp=fopen(file_name, "w"))==NULL){
				printf("write error...\n");
			return 0;
			}
			print_local_axis_raw(xdim, ydim, zdim, threefold, numaxis3, fp, span_tmp, orig_tmp);
		}
		//now print 6 fold axes
	
		if(numfold6>0)
		{
			file_name[i+1] = '6';
			if ((fp=fopen(file_name, "w"))==NULL){
				printf("write error...\n");
			return 0;
			}
			print_local_axis_raw(xdim, ydim, zdim, sixfold, numaxis6, fp, span_tmp, orig_tmp);
		}
		//-Jesse done adding code
		
		
		printf("begin initial subunit segmentation ....\n");
		_result = (unsigned short *)malloc(sizeof(unsigned short)*xdim*ydim*zdim);
		(void)time(&t1); 
		if (numfold3 > 0){
			SubunitSegment(xdim,ydim,zdim,orig_tmp, span_tmp, dataset, _result,tlow, fivefold, threefold, numfold3, numaxis3, numfold5, init_radius);
		}
		else if (numfold6 > 0){
			SubunitSegment(xdim,ydim,zdim,orig_tmp, span_tmp, dataset, _result,tlow, fivefold, sixfold, numfold6, numaxis6, numfold5, init_radius);
		}
		(void)time(&t2); 
		printf("time for initial subunit segmentation: %d seconds. \n\n",(int)(t2-t1));
		
		
		printf("begin local symmetry refinement ....\n");
		(void)time(&t1); 
		strcpy(file_name,filename.c_str());
		for(i = 0; i<MAX_STRING; i++) {
			if (file_name[i] == '.' && file_name[i+1] == 'r' &&
				file_name[i+2] == 'a' && file_name[i+3] == 'w')
				break;
		}
		file_name[i+18] = file_name[i+6];
		file_name[i] = '_';
		file_name[i+1] = 's';
		file_name[i+2] = 'y';
		file_name[i+3] = 'm';
		file_name[i+4] = 'm';
		file_name[i+5] = 'e';
		file_name[i+6] = 't';
		file_name[i+7] = 'r';
		file_name[i+8] = 'y';
		file_name[i+9] = '_';
		file_name[i+10] = 'a';
		file_name[i+11] = 'x';
		file_name[i+12] = 'i';
		file_name[i+13] = 's';
		file_name[i+14] = '.';
		file_name[i+15] = 'r';
		file_name[i+16] = 'a';
		file_name[i+17] = 'w';
		if ((fp1=fopen(file_name, "w"))==NULL){
			printf("write error...\n");
			return 0; 
		} 
		if (numfold3 > 0) {
			strcpy(file_name,filename.c_str());
			for(i = 0; i<MAX_STRING; i++) {
				if (file_name[i] == '.' && file_name[i+1] == 'r' &&
					file_name[i+2] == 'a' && file_name[i+3] == 'w')
					break;
			}
			file_name[i+12] = file_name[i+6];
			file_name[i] = '_';
			file_name[i+1] = '3';
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
			if ((fp=fopen(file_name, "w"))==NULL){
				printf("write error...\n");
				return 0; 
			}
			LocalSymmetryRefine(xdim,ydim,zdim, dataset, critical_list,_result,fivefold, numfold5, threefold, fp, fp1, numfold3,numaxis3, span_tmp, orig_tmp, tlow);
			//-Jesse-Adding code here
			strcpy(file_name,filename.c_str());
			for(i = 0; i<MAX_STRING; i++) {
				if (file_name[i] == '.' && file_name[i+1] == 'r' &&
					file_name[i+2] == 'a' && file_name[i+3] == 'w')
					break;
			}
			file_name[i+12] = file_name[i+6];
			file_name[i] = '_';
			file_name[i+1] = '3';
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
			}
			print_local_axis_raw(xdim, ydim, zdim, threefold, numaxis3, fp, span_tmp, orig_tmp);
			//-Jesse done adding code
		}
		else if (numfold6 > 0) {
			strcpy(file_name,filename.c_str());
			for(i = 0; i<MAX_STRING; i++) {
				if (file_name[i] == '.' && file_name[i+1] == 'r' &&
					file_name[i+2] == 'a' && file_name[i+3] == 'w')
					break;
			}
			file_name[i+12] = file_name[i+6];
			file_name[i] = '_';
			file_name[i+1] = '6';
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
			if ((fp=fopen(file_name, "w"))==NULL){
				printf("write error...\n");
				return 0; 
			};
			LocalSymmetryRefine(xdim,ydim,zdim, dataset, critical_list,_result,fivefold, numfold5, sixfold, fp, fp1, numfold6,numaxis6, span_tmp, orig_tmp, tlow);
			//-Jesse-Adding code here
			strcpy(file_name,filename.c_str());
			for(i = 0; i<MAX_STRING; i++) {
				if (file_name[i] == '.' && file_name[i+1] == 'r' &&
					file_name[i+2] == 'a' && file_name[i+3] == 'w')
					break;
			}
			file_name[i+12] = file_name[i+6];
			file_name[i] = '_';
			file_name[i+1] = '6';
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
			if ((fp=fopen(file_name, "w"))==NULL){
				printf("write error...\n");
				return 0;
			}
			print_local_axis_raw(xdim, ydim, zdim, sixfold, numaxis6, fp, span_tmp, orig_tmp);
			//-Jesse done adding code
		}
		(void)time(&t2); 
		printf("time for local symmetry refinement: %d seconds. \n\n",(int)(t2-t1));
		
		printf("begin covariance matrix calculation ....\n");
		coVarMatrix3 = (CVM*)malloc(sizeof(CVM)*numaxis3*numaxis3); 
		coVarMatrix6 = (CVM*)malloc(sizeof(CVM)*numaxis6*numaxis6); 
		strcpy(file_name,filename.c_str());
		for(i = 0; i<MAX_STRING; i++) {
			if (file_name[i] == '.' && file_name[i+1] == 'r' &&
				file_name[i+2] == 'a' && file_name[i+3] == 'w')
				break;
		}
		file_name[i+21] = file_name[i+6];
		file_name[i] = '_';
		file_name[i+1] = 's';
		file_name[i+2] = 'i';
		file_name[i+3] = 'm';
		file_name[i+4] = 'i';
		file_name[i+5] = 'l';
		file_name[i+6] = 'a';
		file_name[i+7] = 'r';
		file_name[i+8] = 'i';
		file_name[i+9] = 't';
		file_name[i+10] = 'y';
		file_name[i+11] = '_';
		file_name[i+12] = 's';
		file_name[i+13] = 'c';
		file_name[i+14] = 'o';
		file_name[i+15] = 'r';
		file_name[i+16] = 'e';
		file_name[i+17] = '.';
		file_name[i+18] = 't';
		file_name[i+19] = 'x';
		file_name[i+20] = 't';
		if ((fp=fopen(file_name, "w"))==NULL){
			printf("write error...\n");
			return 0; 
		}
		(void)time(&t1); 
		if (numfold3 > 0){
			CoVarianceRefine(xdim,ydim,zdim, dataset, critical_list,_result,threefold, numaxis3, coVarMatrix3,numfold5, fp);
		}
		else if (numfold6 > 0){
			CoVarianceRefine(xdim,ydim,zdim, dataset, critical_list,_result,sixfold,numaxis6, coVarMatrix6,numfold5, fp);
		}
		(void)time(&t2); 
		printf("time for covariance matrix calculation: %d seconds. \n\n",(int)(t2-t1));
		
		
		printf("begin subunit segmentation refinement ....\n");
		(void)time(&t1); 
		if (numfold3 > 0){
			SubunitSegmentRefine(xdim,ydim,zdim, orig_tmp, span_tmp, dataset, _result,tlow, fivefold, threefold, numfold3, numaxis3, coVarMatrix3,numfold5,init_radius);
		}
		else if (numfold6 > 0){
			SubunitSegmentRefine(xdim,ydim,zdim,orig_tmp, span_tmp, dataset, _result,tlow, fivefold, sixfold, numfold6, numaxis6, coVarMatrix6,numfold5,init_radius);
		}
		(void)time(&t2); 
		printf("time for subunit segmentation refinement: %d seconds. \n\n",(int)(t2-t1));
		
		
		(void)time(&t1); 
		printf("begin writing subunit indexing ....\n");
		strcpy(file_name,filename.c_str());
		for(i = 0; i<MAX_STRING; i++) {
			if (file_name[i] == '.' && file_name[i+1] == 'r' &&	file_name[i+2] == 'a' && file_name[i+3] == 'w'){
				break;
			}
		}
		file_name[i+12] = file_name[i+6];
		file_name[i] = '_';
		file_name[i+1] = 'i';
		file_name[i+2] = 'n';
		file_name[i+3] = 'd';
		file_name[i+4] = 'e';
		file_name[i+5] = 'x';
		file_name[i+6] = '.';
		file_name[i+7] = 'r';
		file_name[i+8] = 'a';
		file_name[i+9] = 'w';
		file_name[i+10] = 'i';
		file_name[i+11] = 'v';
		if ((fp=fopen(file_name, "wb"))==NULL){
			printf("write error...\n");
			return 0; 
		}
		write_rawiv_short(_result,fp);
		
		printf("begin writing segmentation results ....\n");
		strcpy(file_name,filename.c_str());
		for(i = 0; i<MAX_STRING; i++) {
			if (file_name[i] == '.' && file_name[i+1] == 'r' && file_name[i+2] == 'a' && file_name[i+3] == 'w'){
				break;
			}
		}
		file_name[i+9] = file_name[i+6];
		file_name[i] = '_';
		file_name[i+1] = 's';
		file_name[i+2] = 'e';
		file_name[i+3] = 'g';
		file_name[i+4] = '.';
		file_name[i+5] = 'r';
		file_name[i+6] = 'a';
		file_name[i+7] = 'w';
		file_name[i+8] = 'v';
		if ((fp=fopen(file_name, "wb"))==NULL){
			printf("write error...\n");
			return 0; 
		}
		if (numfold3 > 0){
			write_rawv(dataset,_result,numaxis3,numfold5,orig_tmp, span_tmp, 0,fp);
		}
		else if (numfold6 > 0){
			write_rawv(dataset,_result,numaxis6,numfold5,orig_tmp, span_tmp, 1,fp);
		}

		if (numfold3 > 0){
			printf("begin computing and writing the 3-fold averaged subunit....\n");
			strcpy(file_name,filename.c_str());
			for(i = 0; i<MAX_STRING; i++) {
				if (file_name[i] == '.' && file_name[i+1] == 'r' && file_name[i+2] == 'a' && file_name[i+3] == 'w'){
					break;
				}
			}
			file_name[i+16] = file_name[i+6];
			file_name[i] = '_';
			file_name[i+1] = '3';
			file_name[i+2] = 'f';
			file_name[i+3] = '_';
			file_name[i+4] = 's';
			file_name[i+5] = 'u';
			file_name[i+6] = 'b';
			file_name[i+7] = 'a';
			file_name[i+8] = 'v';
			file_name[i+9] = 'g';
			file_name[i+10] = '.';
			file_name[i+11] = 'r';
			file_name[i+12] = 'a';
			file_name[i+13] = 'w';
			file_name[i+14] = 'i';
			file_name[i+15] = 'v';
			if ((fp=fopen(file_name, "wb"))==NULL){
				printf("write error...\n");
				return 0; 
			}
			ComputeSubAverage(xdim,ydim,zdim,dataset,_result,threefold,numfold3,numaxis3, coVarMatrix3,numfold5,span_tmp,orig_tmp,fp);
		}
		else if (numfold6 > 0) {
			printf("begin computing and writing the 6-fold averaged subunit....\n");
			strcpy(file_name,filename.c_str());
			for(i = 0; i<MAX_STRING; i++) {
				if (file_name[i] == '.' && file_name[i+1] == 'r' && file_name[i+2] == 'a' && file_name[i+3] == 'w'){
					break;
				}
			}
			file_name[i+16] = file_name[i+6];
			file_name[i] = '_';
			file_name[i+1] = '6';
			file_name[i+2] = 'f';
			file_name[i+3] = '_';
			file_name[i+4] = 's';
			file_name[i+5] = 'u';
			file_name[i+6] = 'b';
			file_name[i+7] = 'a';
			file_name[i+8] = 'v';
			file_name[i+9] = 'g';
			file_name[i+10] = '.';
			file_name[i+11] = 'r';
			file_name[i+12] = 'a';
			file_name[i+13] = 'w';
			file_name[i+14] = 'i';
			file_name[i+15] = 'v';
			if ((fp=fopen(file_name, "wb"))==NULL){
				printf("write error...\n");
				return 0; 
			}
			ComputeSubAverage(xdim,ydim,zdim,dataset,_result,sixfold,numfold6,numaxis6, coVarMatrix6,numfold5,span_tmp,orig_tmp,fp);
		}
		
		if (numfold5 > 0){
			printf("begin writing the 5-fold subunit....\n");
			strcpy(file_name,filename.c_str());
			for(i = 0; i<MAX_STRING; i++) {
				if (file_name[i] == '.' && file_name[i+1] == 'r' && file_name[i+2] == 'a' && file_name[i+3] == 'w'){
					break;
				}
			}
			file_name[i+17] = file_name[i+6];
			file_name[i] = '_';
			file_name[i+1] = '5';
			file_name[i+2] = 'f';
			file_name[i+3] = '_';
			file_name[i+4] = 's';
			file_name[i+5] = 'u';
			file_name[i+6] = 'b';
			file_name[i+7] = 'u';
			file_name[i+8] = 'n';
			file_name[i+9] = 'i';
			file_name[i+10] = 't';
			file_name[i+11] = '.';
			file_name[i+12] = 'r';
			file_name[i+13] = 'a';
			file_name[i+14] = 'w';
			file_name[i+15] = 'i';
			file_name[i+16] = 'v';
			if ((fp=fopen(file_name, "wb"))==NULL){
				printf("write error...\n");
				return 0; 
			}
			Write5fSubunit(xdim,ydim,zdim,dataset,_result,span_tmp,orig_tmp,fp);
		}
		
		printf("begin making transformation matrices....\n");
		strcpy(file_name,filename.c_str());
		for(i = 0; i<MAX_STRING; i++) {
			if (file_name[i] == '.' && file_name[i+1] == 'r' && file_name[i+2] == 'a' && file_name[i+3] == 'w'){
				break;
			}
		}
		file_name[i+14] = file_name[i+6];
		file_name[i] = '_';
		file_name[i+1] = '6';
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
		}
		if (numfold3 > 0){
			MakeTransformMatrix(xdim,ydim,zdim,fivefold,threefold,numfold3,numaxis3,coVarMatrix3,span_tmp,orig_tmp,fp);
		}
		else if (numfold6 > 0){
			MakeTransformMatrix(xdim,ydim,zdim,fivefold,sixfold,numfold6,numaxis6,coVarMatrix6,span_tmp,orig_tmp,fp);
		}
		(void)time(&t2); 
		printf("time for writing results: %d seconds. \n\n",(int)(t2-t1));
		
	}
	else if (numfold5 > 0){
		_result = (unsigned short *)malloc(sizeof(unsigned short)*xdim*ydim*zdim);
		(void)time(&t1); 
		AsymSubunitSegment(xdim,ydim,zdim,orig_tmp, span_tmp, dataset,_result,tlow,200,critical_list,fivefold);
		(void)time(&t2); 
		printf("time for asymmetric subunit segmentation: %d seconds. \n\n",(int)(t2-t1));
		
		(void)time(&t1); 
		printf("begin writing subunit indexing ....\n");
		strcpy(file_name,filename.c_str());
		for(i = 0; i<MAX_STRING; i++){
			if (file_name[i] == '.' && file_name[i+1] == 'r' && file_name[i+2] == 'a' && file_name[i+3] == 'w'){
				break;
			}
		}
		file_name[i+12] = file_name[i+6];
		file_name[i] = '_';
		file_name[i+1] = 'i';
		file_name[i+2] = 'n';
		file_name[i+3] = 'd';
		file_name[i+4] = 'e';
		file_name[i+5] = 'x';
		file_name[i+6] = '.';
		file_name[i+7] = 'r';
		file_name[i+8] = 'a';
		file_name[i+9] = 'w';
		file_name[i+10] = 'i';
		file_name[i+11] = 'v';
		if ((fp=fopen(file_name, "wb"))==NULL){
			printf("write error...\n");
			return 0; 
		}
		write_rawiv_short(_result,fp);
		
		printf("begin writing segmentation results ....\n");
		strcpy(file_name,filename.c_str());
		for(i = 0; i<MAX_STRING; i++) {
			if (file_name[i] == '.' && file_name[i+1] == 'r' && file_name[i+2] == 'a' && file_name[i+3] == 'w'){
				break;
			}
		}
		file_name[i+9] = file_name[i+6];
		file_name[i] = '_';
		file_name[i+1] = 's';
		file_name[i+2] = 'e';
		file_name[i+3] = 'g';
		file_name[i+4] = '.';
		file_name[i+5] = 'r';
		file_name[i+6] = 'a';
		file_name[i+7] = 'w';
		file_name[i+8] = 'v';
		if ((fp=fopen(file_name, "wb"))==NULL){
			printf("write error...\n");
			return 0; 
		}
//		write_asym_rawv(dataset,_result,fp);
		write_asym_rawv(dataset,_result, span_tmp, orig_tmp, fp);
		
		printf("begin writing the 5-fold subunit....\n");
		strcpy(file_name,filename.c_str());
		for(i = 0; i<MAX_STRING; i++) {
			if (file_name[i] == '.' && file_name[i+1] == 'r' && file_name[i+2] == 'a' && file_name[i+3] == 'w'){
				break;
			}
		}
		file_name[i+19] = file_name[i+6];
		file_name[i] = '_';
		file_name[i+1] = 'a';
		file_name[i+2] = 's';
		file_name[i+3] = 'y';
		file_name[i+4] = 'm';
		file_name[i+5] = '_';
		file_name[i+6] = 's';
		file_name[i+7] = 'u';
		file_name[i+8] = 'b';
		file_name[i+9] = 'u';
		file_name[i+10] = 'n';
		file_name[i+11] = 'i';
		file_name[i+12] = 't';
		file_name[i+13] = '.';
		file_name[i+14] = 'r';
		file_name[i+15] = 'a';
		file_name[i+16] = 'w';
		file_name[i+17] = 'i';
		file_name[i+18] = 'v';
		if ((fp=fopen(file_name, "wb"))==NULL){
			printf("write error...\n");
			return 0; 
		}
		Write5fSubunit(xdim,ydim,zdim,dataset,_result,span_tmp,orig_tmp,fp);
		(void)time(&t2); 
		printf("time for writing results: %d seconds. \n\n",(int)(t2-t1));
		
	}
	else {
		printf("Failed: one of <numfold3, numfold5,numfold6> must positive integer !!!\n");
		return 0;
	}
	
	free(triangle);
	free(vertex);
	result = XmlRpcValue(true);
	return 1;
}

void print_local_axis_raw(int xdim, int ydim, int zdim, DB_VECTOR *sym_axes, int num_axes, FILE *fp, float *span_tmp, float *orig_tmp){
	v_num = 0;
	t_num = 0;
	int n;
	
	for (n=0; n<num_axes; n++) {
	  	DrawLine_subunit(sym_axes[n].sx, sym_axes[n].sy, sym_axes[n].sz, sym_axes[n].ex, sym_axes[n].ey, sym_axes[n].ez, 3.0f);
	}
	fprintf(fp, "%d %d\n",v_num,t_num);
	for (n = 0; n < v_num; n++){
//		fprintf(fp, "%f %f %f\n",vertex[n].x*span_tmp[0]+orig_tmp[0], 
//			vertex[n].y*span_tmp[1]+orig_tmp[1],vertex[n].z*span_tmp[2]+orig_tmp[2]);
		fprintf(fp, "%f %f %f\n",vertex[n].x, vertex[n].y,vertex[n].z);


	}
	for (n = 0; n < t_num; n++){
		fprintf(fp, "%d %d %d\n",triangle[n].x,triangle[n].y,triangle[n].z);
	}
	fclose(fp);
}

void DrawLine_subunit(float sx, float sy, float sz, float ex, float ey, float ez, float radius){
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
