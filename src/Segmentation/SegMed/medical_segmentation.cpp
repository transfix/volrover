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
#include <string>
#include <xmlrpc/XmlRpc.h>
#include <Segmentation/SegMed/segmed.h>

#define IndexVect(i,j,k) ((k)*xdim*ydim + (j)*xdim + (i))
#define max(x, y) ((x>y) ? (x):(y))
#define min(x, y) ((x<y) ? (x):(y))

using namespace XmlRpc;
using namespace SegMed;
using namespace std;

int medicalSegmentation(XmlRpcValue &params, XmlRpcValue &result)
{
  int xdim,ydim,zdim;
  float *image;
  float *smth_img;
  float *edge_mag;
  FILE *fp;
//  FILE *fp2;
  CPNT *critical_list;
  float tlow;
//  float thigh;
  unsigned char *_result;
  string filename;
  string output_filename;

  /*
  if (argc != 4){
    printf("Usage: CCVskel <input_filename> <output_rawiv> <tlow>\n");
    printf("       <input_filename>:   RAWIV file \n");
    printf("       <output_rawiv>:  RAWIV file \n");
    printf("       <tlow>: thresholds for segemntation (0-255) \n");
    exit(0);              
  }
  */
  
  printf("begin reading rawiv.... \n");
  filename = string(params[0]);
  read_data(&xdim,&ydim,&zdim,&image,filename.c_str());

  smth_img = image;
  /*
  smth_img = (float*)malloc(sizeof(float)*xdim*ydim*zdim);
  printf("begin diffusion .... \n");
  Diffuse(xdim,ydim,zdim, image, smth_img);
  */
  
  tlow = double(params[1]);  
  edge_mag = (float*)malloc(sizeof(float)*xdim*ydim*zdim);
  printf("begin GVF computation....\n");
  GVF_Compute(xdim,ydim,zdim, smth_img, edge_mag, &critical_list,tlow);
  

  _result = (unsigned char*)malloc(sizeof(unsigned char)*xdim*ydim*zdim);
  printf("begin segmentation ....\n");
  Segment(xdim,ydim,zdim, smth_img,edge_mag,_result,critical_list,tlow);
  
  output_filename = filename + ".seg.rawiv";
  if ((fp=fopen(output_filename.c_str(), "w"))==NULL){
    printf("read error...\n");
    return 0;
  }
  printf("Writing data....\n");
  write_data(xdim,ydim,zdim, smth_img,_result, fp);

  free(edge_mag);
  free(_result);

  result = XmlRpcValue(true);
  return 1;
}
