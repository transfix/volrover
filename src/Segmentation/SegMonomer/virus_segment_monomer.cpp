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
#include <Segmentation/SegMonomer/segmonomer.h>

#define _LITTLE_ENDIAN 1

#define IndexVect(i,j,k) ((k)*xdim*ydim + (j)*xdim + (i))
#define max2(x, y) ((x>y) ? (x):(y))
#define min2(x, y) ((x<y) ? (x):(y))
#define PIE           3.1415926536f
#define MAX_STRING    256

using namespace XmlRpc;
using namespace SegMonomer;

int virusSegMonomer(XmlRpcValue &params, XmlRpcValue &result)
{
	int xdim,ydim,zdim;
	int xdim_big,ydim_big,zdim_big;
	float *dataset;
	float *span_tmp, *orig_tmp, *orig_big;
	float *edge_mag;
	FILE *fp;
	char file_name[MAX_STRING];
	CPNT *critical_list;
	unsigned char *_result;
	DB_VECTOR *symmetry_axis;
	int foldnum;
	float shift_x,shift_y,shift_z;
	float fx,fy,fz;
	float gx,gy,gz;
	float score;
	int i,j,k;
	int num;
	time_t t1,t2;
	std::string filename;
	
	int fscanf_return = 0;
	
	filename = std::string(params[0]);
	
	printf("begin reading rawiv.... \n");
	span_tmp = (float *)malloc(sizeof(float)*3);
	orig_tmp = (float *)malloc(sizeof(float)*3);
	orig_big = (float *)malloc(sizeof(float)*3);
	(void)time(&t1); 

	strcpy(file_name,filename.c_str());
	for(i = 0; i<MAX_STRING; i++) {
		if (file_name[i] == '.' && file_name[i+1] == 'r' &&
			file_name[i+2] == 'a' && file_name[i+3] == 'w')
			break;
	}
	j = i;
	while (file_name[i] != '_')
		i--;
	i--;
	while (file_name[i] != '_')
		i--;
	file_name[i+6] = file_name[j+6];
	file_name[i] = '.';
	file_name[i+1] = 'r';
	file_name[i+2] = 'a';
	file_name[i+3] = 'w';
	file_name[i+4] = 'i';
	file_name[i+5] = 'v';

	if ((fp=fopen(file_name, "r"))==NULL){
		printf("read error...\n");
		return 0; 
	};  
	get_header(&xdim_big,&ydim_big,&zdim_big,orig_big,fp);

	read_data(&xdim,&ydim,&zdim,&dataset,span_tmp,orig_tmp,const_cast<char *>(filename.c_str()),0,1);
	printf("xdim: %d  ydim: %d  zdim: %d \n",xdim,ydim,zdim);
	(void)time(&t2); 
	printf("time to read dataset: %d seconds. \n\n",(int)(t2-t1));
	
	symmetry_axis = (DB_VECTOR*)malloc(sizeof(DB_VECTOR));
	shift_x = (orig_tmp[0]-orig_big[0])/span_tmp[0];
	shift_y = (orig_tmp[1]-orig_big[1])/span_tmp[1];
	shift_z = (orig_tmp[2]-orig_big[2])/span_tmp[2];
	
	strcpy(file_name,filename.c_str());
	for(i = 0; i<MAX_STRING; i++) {
		if (file_name[i] == '.' && file_name[i+1] == 'r' &&
			file_name[i+2] == 'a' && file_name[i+3] == 'w')
			break;
	}
	j = i;
	while (file_name[i] != '_')
		i--;
	i--;
	while (file_name[i] != '_')
		i--;
	if (file_name[i+1] == '5') {
		i--;
		while (file_name[i] != '_')
			i--;
	//	file_name[i+12] = file_name[j+6];
//		file_name[i] = '_';
//		file_name[i+1] = '5';
//		file_name[i+2] = 'f';
//		file_name[i+3] = '_';
///		file_name[i+4] = 'a';
//		file_name[i+5] = 'x';
//		file_name[i+6] = 'i';
//		file_name[i+7] = 's';
//		file_name[i+8] = '.';
//		file_name[i+9] = 't';
//		file_name[i+10] = 'x';
//		file_name[i+11] = 't';
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
		};    
		fscanf_return = fscanf(fp, "%f %f %f \n",&fx,&fy,&fz);
		symmetry_axis->sx = fx; //-shift_x;
		symmetry_axis->sy = fy; //-shift_y;
		symmetry_axis->sz = fz; // -shift_z;
		symmetry_axis->ex = (orig_big[0]+ xdim_big*span_tmp[0])/2;  // (float)xdim_big/2.0f-shift_x;
		symmetry_axis->ey = (orig_big[1]+ ydim_big*span_tmp[1])/2;  //(float)ydim_big/2.0f-shift_y;
		symmetry_axis->ez = (orig_big[2]+ zdim_big*span_tmp[2])/2;  //(float)zdim_big/2.0f-shift_z;
		fclose(fp);
	}
	else {
		file_name[i+12] = file_name[j+6];
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
		fscanf_return = fscanf(fp, "%f %f %f %f %f %f %f\n",&fx,&fy,&fz,&gx,&gy,&gz,&score);
		symmetry_axis->sx = fx; //-shift_x;
		symmetry_axis->sy = fy; //-shift_y;
		symmetry_axis->sz = fz; //-shift_z;
		symmetry_axis->ex = gx; //-shift_x;
		symmetry_axis->ey = gy; //-shift_y;
		symmetry_axis->ez = gz; //-shift_z;
		fclose(fp);
	}
	
	foldnum = int(params[1]);
	printf("Compute averaging map ....\n");
	(void)time(&t1); 
	get_average(xdim,ydim,zdim, dataset,orig_tmp, span_tmp, symmetry_axis,foldnum);
	(void)time(&t2); 
	printf("time to compute averaging map : %d seconds. \n\n",(int)(t2-t1));
	
	
	edge_mag = (float*)malloc(sizeof(float)*xdim*ydim*zdim);
	printf("begin GVF computation....\n");
	(void)time(&t1); 
	GVF_Compute(xdim,ydim,zdim, dataset, edge_mag, &critical_list);
	(void)time(&t2); 
	printf("time to find critical points: %d seconds. \n\n",(int)(t2-t1));
	
	
	_result = (unsigned char*)malloc(sizeof(unsigned char)*xdim*ydim*zdim);
	printf("begin segmentation ....\n");
	(void)time(&t1); 
	MonomerSegment(xdim,ydim,zdim, orig_tmp, span_tmp, dataset, edge_mag, 
				   critical_list,_result,symmetry_axis,foldnum);
	(void)time(&t2); 
	printf("time to segment monomers: %d seconds. \n\n",(int)(t2-t1));
	
	
	printf("begin writing the monomers ....\n");
	for (num = 1; num <= foldnum; num++) {
		strcpy(file_name,filename.c_str());
		for(i = 0; i<MAX_STRING; i++) {
			if (file_name[i] == '.' && file_name[i+1] == 'r' &&
				file_name[i+2] == 'a' && file_name[i+3] == 'w')
				break;
		}
		file_name[i+17] = file_name[i+6];
		file_name[i] = '_';
		file_name[i+1] = 'm';
		file_name[i+2] = 'o';
		file_name[i+3] = 'n';
		file_name[i+4] = 'o';
		file_name[i+5] = 'm';
		file_name[i+6] = 'e';
		file_name[i+7] = 'r';
		file_name[i+8] = '_';
		file_name[i+9] = 48+num/10;
		file_name[i+10] = 48+num-10*(num/10);
		file_name[i+11] = '.';
		file_name[i+12] = 'r';
		file_name[i+13] = 'a';
		file_name[i+14] = 'w';
		file_name[i+15] = 'i';
		file_name[i+16] = 'v';
		if ((fp=fopen(file_name, "wb"))==NULL){
			printf("write error...\n");
			return 0; 
		};
		for (k=0; k<zdim; k++)
			for (j=0; j<ydim; j++) 
				for (i=0; i<xdim; i++) {
					if (_result[IndexVect(i,j,k)] == num)  
						edge_mag[IndexVect(i,j,k)] = dataset[IndexVect(i,j,k)];
					else
						edge_mag[IndexVect(i,j,k)] = 0;
				}
					write_rawiv_float(edge_mag,fp);
	}


  result = XmlRpcValue(true);
  return 1;
}
