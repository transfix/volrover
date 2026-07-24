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
#include <math.h>
#include <time.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <xmlrpc/XmlRpc.h>
#include <Segmentation/SecStruct/secstruct.h>

#include <Segmentation/SecStruct/fit_cylinder.h>

#define IndexVect(i,j,k) ((k)*XDIM*YDIM + (j)*XDIM + (i))
#define MAX_STRING    256

using namespace XmlRpc;
using namespace SecStruct;

int secondaryStructureDetection(XmlRpc::XmlRpcValue &params, XmlRpc::XmlRpcValue &result)
{
  char file_name[MAX_STRING];
  VECTOR* velocity;
  unsigned char *result_;
  float tlow;
  float helixwidth,min_hratio,max_hratio,min_hlength;
  float sheetwidth,min_sratio,max_sratio,sheet_extend;
  int XDIM, YDIM, ZDIM;
  float *dataset;
  float *span_tmp, *orig_tmp;
  FILE *fp,*fp1,*fp2;
  int i;
  time_t t1,t2;
  SEEDS *AllHelix_list; 
  int HSeedNum;
  
  std::string filename, newfilename;
  
  filename = std::string(params[0]);
  
#if 0
  if (argc != 10){
    printf("Usage: HelixSheetHunter <input> <helixwidth> <min_helixwidth_ratio> <max_helixwidth_ratio> <min_helixlength> \n");
    printf("                        <sheetwidth> <min_sheetwidth_ratio> <max_sheetwidth_ratio> <sheet_extend>  \n");
    printf("       <input>: Rawiv file \n");
    printf("       <helixwidth> : the thickness of typical helices (in pixels) \n");
    printf("       <min_helixwidth_ratio> : the low ratio of thickness (0~1) \n");
    printf("       <max_helixwidth_ratio> : the high fofthickness (> 1) \n");
    printf("       <min_helixlength> : the shortest helices (in pixels) \n");
    printf("       <sheetwidth> : the thickness of typical sheets (in pixels) \n");
    printf("       <min_sheetwidth_ratio> : the low ratio of thickness (0~1) \n");
    printf("       <max_sheetwidth_ratio> : the high ratio of thickness (> 1) \n");
    printf("       <sheet_extend> : the extension ratio of sheets (1~2) \n");
    exit(0);              
  }
#endif

  (void)time(&t1); 
  
  printf("Loading dataset... \n");
  span_tmp = (float *)malloc(sizeof(float)*3);
  orig_tmp = (float *)malloc(sizeof(float)*3);
  read_data(&XDIM,&YDIM,&ZDIM,&dataset,span_tmp,orig_tmp,filename.c_str());
  printf("Dataset loaded\n");
  

  /*
  printf("Begin Diffusion....\n");
  Diffuse(dataset,XDIM,YDIM,ZDIM);
  */
  
  velocity = (VECTOR*)malloc(sizeof(VECTOR)*XDIM*YDIM*ZDIM);
  printf("Begin GVF computation....\n");
  GVF_Compute(XDIM,YDIM,ZDIM,dataset,velocity);

  result_ = (unsigned char*)malloc(sizeof(unsigned char)*XDIM*YDIM*ZDIM);
  printf("begin detecting features ....\n");
  helixwidth = double(params[1]);
  min_hratio = double(params[2]);
  max_hratio = double(params[3]);
  min_hlength = double(params[4]);
  sheetwidth = double(params[5]);
  min_sratio = double(params[6]);
  max_sratio = double(params[7]);
  sheet_extend = double(params[8]);
  tlow =  double(params[9]); //if not selected tlow =-1; 
  
  StructureTensor(XDIM,YDIM,ZDIM,velocity,result_,dataset,&tlow,helixwidth,min_hratio,
		  max_hratio, sheetwidth,min_sratio,max_sratio, &HSeedNum);
  
  (void)time(&t2); 
  printf("time used for preprocessing: %d seconds. \n\n",(int)(t2-t1));
 
  size_t pos;

  pos = filename.find('.');

  filename = filename.erase(pos);

  (void)time(&t1); 
  #ifdef VRML_OUTPUT
  newfilename = filename + "_skeleton.wrl";
  #else
  newfilename = filename + "_skeleton.rawc";
  #endif

  if ((fp=fopen(newfilename.c_str(), "w"))==NULL){
    printf("write error 1...\n");
    result = int(0);
    return 0;
  };
  AllHelix_list = (SEEDS*)malloc(sizeof(SEEDS)*HSeedNum);
  HelixHunter(XDIM,YDIM,ZDIM,span_tmp,orig_tmp, dataset, result_, 
	      velocity,tlow,helixwidth,min_hlength,fp,AllHelix_list,HSeedNum);
  
  (void)time(&t2); 
  printf("time used for helix-hunter: %d seconds. \n\n",(int)(t2-t1));
  
  (void)time(&t1); 

  #ifdef VRML_OUTPUT
  newfilename = filename + "_sheet.wrl";
  #else
   newfilename = filename + "_sheet.rawc";
  #endif
  if ((fp=fopen(newfilename.c_str(), "w"))==NULL){
    printf("write error 2...\n");
    result = int(0);
	return 0;
  };  

  #ifdef VRML_OUTPUT
  newfilename = filename + "_helix_s.wrl";
  #else
  newfilename = filename + "_helix_s.rawc";
  #endif
  if ((fp1=fopen(newfilename.c_str(), "w"))==NULL){
    printf("write error 3...\n");
    result = int(0);
	return 0; 
  };  

  #ifdef VRML_OUTPUT
  newfilename = filename + "_helix_c.wrl";
  #else
  newfilename = filename + "_helix_c.rawc";
  #endif
  if ((fp2=fopen(newfilename.c_str(), "w"))==NULL){
    printf("write error 4...\n");
    result = int(0);
	return 0; 
  };

  SheetHunter(XDIM,YDIM,ZDIM,span_tmp,orig_tmp, dataset, result_, velocity, tlow,helixwidth,min_hratio,
	      max_hratio,min_hlength,sheetwidth,min_sratio,max_sratio,sheet_extend,fp,AllHelix_list,HSeedNum,fp1,fp2);
  
  (void)time(&t2); 
  printf("time used for sheet-hunter: %d seconds. \n\n",(int)(t2-t1));
  
 
  free(dataset);
  free(result_);
  free(velocity);
  
 #if 0
  std::vector<Point> pts;
  pts.push_back(Point(0.0,0.0,0.0));
  pts.push_back(Point(10.0,0.0,0.0));
  Vector vec(1.0,0.0,0.0);
  std::vector<std::vector<Point> > cyl_Vs(fit_cylinder(std::vector<Cylinder>(1,Cylinder(pts,vec))));
   
  for(std::vector<std::vector<Point> >::iterator cyl_Vi = cyl_Vs.begin();
      cyl_Vi != cyl_Vs.end();
      cyl_Vi++)
    {
      std::vector<Point> cyl_V;

      cyl_V = *cyl_Vi;

      char file_prefix[100] = "cyl";
      char op_fname[100];
      char extn[10];
      extn[0] = '_'; extn[1] = '0' + i/10; extn[2] = '0' + i%10; extn[3] = '\0';
      strcpy(op_fname, file_prefix);
      strcat(op_fname, extn);
      strcat(op_fname, ".cyl");
      std::cerr << "file : " << op_fname << std::endl;

      std::ofstream fout;
      fout.open(op_fname);
      if(! fout)
      {
         std::cerr << "Error in opening output file " << std::endl;
         return 0;
      }

      fout << (int)cyl_V.size() << " " << (int)cyl_V.size() << std::endl;
      for(int j = 0; j < (int)cyl_V.size(); j ++)
         fout << cyl_V[j] << std::endl;
      int N = (int)cyl_V.size()/2;
      for(int j = 0; j < N; j ++)
      {
         int v1 = j, v2 = j+N, v3 = (j+1)%N; 
         fout << v1 << " " << v2 << " " << v3 << std::endl;
         v1 = v3, v3 = v2+1;
         v3 = (v3 == 2*N)? N : v3;
         fout << v1 << " " << v2 << " " << v3 << std::endl;
      }
      fout.close();
    }
      #endif

  result = XmlRpcValue(true);
  return(1);
}

