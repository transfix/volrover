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
#include <string.h>
#include <iostream>
#include <string>
#include <xmlrpc/XmlRpc.h>
#include <Segmentation/GenSeg/genseg.h>

using namespace XmlRpc;
using namespace GenSeg;
using namespace std;

int generalSegmentation(XmlRpcValue& params, XmlRpcValue& result)
{
  int xdim,ydim,zdim;
  float *image;
  
  SDPNT **seed_list,*cur;
  int index = 5,i,j;
	
  const char *filename = string(params[0]).c_str();
  int tlow = int(params[1]);
  int thigh = int(params[2]);
  int num_classes = int(params[3]);
  const char *output_filename = string(params[4]).c_str();
	
  printf("filename: %s\n",filename);
  printf("output_filename_prefix: %s\n",output_filename);
  printf("tlow = %d, thigh = %d\n",tlow,thigh);
  printf("num_classes = %d\n",num_classes);
	
  // initialize the list of seed points
  seed_list = new SDPNT*[num_classes];
  memset(seed_list,0,sizeof(SDPNT*)*num_classes);
  for(i=0; i<num_classes; i++)
    {
      int num_points = int(params[index]);
      printf("Class %d: Num points: %d, ",i,num_points);
      index++;
	 
      if(num_points>0) cur = seed_list[i] = new SDPNT;
      for(j=0; j<num_points; j++)
	{
	  cur->x = int(params[index+0]);
	  cur->y = int(params[index+1]);
	  cur->z = int(params[index+2]);
	  printf("(%d,%d,%d) ",cur->x,cur->y,cur->z);
	  index+=3;
	  if(j+1<num_points)
	    {
	      cur->next = new SDPNT;
	      cur = cur->next;
	    }
	  else
	    {
	      cur->next = NULL;
	    }
	}
      printf("\n");
    }

  //execute segmentation
  printf("begin reading rawiv.... \n");
  read_data(&xdim,&ydim,&zdim,&image,filename);
  printf("begin segmentation ....\n");
  Segment(xdim,ydim,zdim, image, tlow, thigh, seed_list,num_classes,output_filename);
	

  // cleanup	
  for(i=0; i<num_classes; i++)
    {
      cur = seed_list[i];
      while(cur != NULL)
	{
	  SDPNT *tmp = cur;
	  cur = cur->next;
	  delete tmp;
	}
    }
  delete [] seed_list;
    
  result = int(true);

  return 1;
}
