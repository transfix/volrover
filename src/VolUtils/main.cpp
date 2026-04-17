/*
  Copyright 2005-2011 The University of Texas at Austin

        Authors: Joe Rivera <transfix@ices.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of VolUtils.

  VolUtils is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  VolUtils is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

/* $Id: main.cpp 4742 2011-10-21 22:09:44Z transfix $ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <fcntl.h>
#ifndef _WIN32
#include <unistd.h>
#endif
#include <search.h>
#include <VolUtils/MappedRawIVFile.h>
#include <VolUtils/MappedRawVFile.h>

typedef struct
{
  float min[3];
  float max[3];
  unsigned int numVerts;
  unsigned int numCells;
  unsigned int dim[3];
  float origin[3];
  float span[3];
} RawIVHeader;

typedef struct
{
  unsigned char r;
  unsigned char g;
  unsigned char b;
} color_t;

MappedVolumeFile *vol = NULL;
FILE *outvol = NULL;

static void cleanup()
{
  printf("Cleaning up...\n");
  if(outvol) fclose(outvol);
  if(vol) delete vol;
}

static int color_cmp(const void *a, const void *b)
{
  if(((color_t *)a)->r == ((color_t *)b)->r &&
     ((color_t *)a)->g == ((color_t *)b)->g &&
     ((color_t *)a)->b == ((color_t *)b)->b) return 0;
  return 1;
}

int main(int argc, char **argv)
{
  unsigned long long i,j,k;
  unsigned int red_var,green_var,blue_var,alpha_var,timestep;
  color_t colors[257]; /* one extra for error handling */
  size_t num_colors = 0;

  atexit(cleanup);

  if(argc != 2 && argc != 8)
    {
      printf("Every color in the input RawVRGBA file (RawV with 4 variables: red, green, blue, and alpha)\n"
	     "will be mapped to some voxel value range within [0-255] in the output rawiv file (unsigned char volume).\n"
	     "The number of values in each color's range is dependent on the number of different voxel colors in the\n"
	     "RawVRGBA volume. The input voxel's alpha value will determine which value within each color's selected\n"
	     "range the output voxel will take.  If the input voxel's alpha value is 0, then the value in the output volume\n"
	     "will be 0. Due to the fact that the output volume voxels are of unsigned char\n"
	     "type, the maximum number of different colors that can be represented in the output volume are 256.  If\n"
	     "there are more than 256 colors, those colors will cause their corresponding output voxel to take the value\n"
	     "0.\n\n");
      printf("Usage: %s <input rawv with >=4 variables> <red var> <green var> <blue var> <alpha var> <timestep> <output rawiv>\n\n",argv[0]);
      printf("Example: %s heart.rawv 0 1 2 3 0 heart.rawiv\n"
	     "Most RawVRGBA files have the RGBA variables in order, thus the above example should work in most cases.\n"
	     "To produce a rendering similar to the RGBA rendering of the input RawV file, make sure to map each voxel\n"
	     "value range to the color used in the RawV file.",argv[0]);
      return 0;
    }

  vol = new MappedRawVFile(argv[1],true,true);
  if(!vol->isValid())
    {
      printf("Error loading %s!\n",argv[1]);
      return 1;
    }

  printf("File: %s\n",argv[1]);
  printf("Num Vars: %d\n",vol->numVariables());
  printf("Vars: ");
  for(i=0; i<vol->numVariables(); i++) printf("%s ",vol->get(i,0)->name());
  printf("\n");
  printf("Num Timesteps: %d\n",vol->numTimesteps());
  printf("Dimensions: %d x %d x %d\n",vol->XDim(),vol->YDim(),vol->ZDim());
  printf("Span: %lf x %lf x %lf\n",vol->XSpan(),vol->YSpan(),vol->ZSpan());
  printf("TSpan: %lf\n",vol->TSpan());

  if(argc == 2) { return 0; } /* only need to print out volume info */

  red_var = atoi(argv[2]);
  green_var = atoi(argv[3]);
  blue_var = atoi(argv[4]);
  alpha_var = atoi(argv[5]);
  timestep = atoi(argv[6]);

  outvol = fopen(argv[7],"wb+");
  if(outvol == NULL)
    {
      char err_str[512];
      sprintf(err_str,"Error opening %s",argv[7]);
      perror(err_str);
      return 1;
    }

  unsigned long long len = vol->XDim()*vol->YDim()*vol->ZDim()*sizeof(unsigned char)+68;

  /* create the header for the new file */
  RawIVHeader header;
  header.min[0] = 0.0; header.min[1] = 0.0; header.min[2] = 0.0;
  header.max[0] = (vol->XDim()-1)*vol->XSpan();
  header.max[1] = (vol->YDim()-1)*vol->YSpan();
  header.max[2] = (vol->ZDim()-1)*vol->ZSpan();
  header.numVerts = vol->XDim()*vol->YDim()*vol->ZDim();
  header.numCells = (vol->XDim()-1)*(vol->YDim()-1)*(vol->ZDim()-1);
  header.dim[0] = vol->XDim(); header.dim[1] = vol->YDim(); header.dim[2] = vol->ZDim();
  header.origin[0] = 0.0; header.origin[1] = 0.0; header.origin[2] = 0.0;
  header.span[0] = vol->XSpan(); header.span[1] = vol->YSpan(); header.span[2] = vol->ZSpan();
	
  if(!big_endian())
    {
      for(i=0; i<3; i++) SWAP_32(&(header.min[i]));
      for(i=0; i<3; i++) SWAP_32(&(header.max[i]));
      SWAP_32(&(header.numVerts));
      SWAP_32(&(header.numCells));
      for(i=0; i<3; i++) SWAP_32(&(header.dim[i]));
      for(i=0; i<3; i++) SWAP_32(&(header.origin[i]));
      for(i=0; i<3; i++) SWAP_32(&(header.span[i]));
    }

  fwrite(&header,sizeof(RawIVHeader),1,outvol);
	
  /* get a slice at a time because it's quicker */
  unsigned char *slice[5]; /* [0-3] == rgba slices, [4] == output slice */
  for(i=0; i<5; i++) slice[i] = (unsigned char *)malloc(vol->XDim()*vol->YDim()*sizeof(unsigned char));
	
  for(k=0; k<vol->ZDim(); k++)
    {
      /* get the colors for this slice (ignoring alpha for now) */
      vol->get(red_var,timestep)->getMapped(0,0,k,vol->XDim(),vol->YDim(),1,slice[0]);
      vol->get(green_var,timestep)->getMapped(0,0,k,vol->XDim(),vol->YDim(),1,slice[1]);
      vol->get(blue_var,timestep)->getMapped(0,0,k,vol->XDim(),vol->YDim(),1,slice[2]);
	  
      /* check each colored voxel and add all new colors to the list to determine output voxel ranges */
      color_t color;
      for(i=0; i<vol->XDim(); i++)
	for(j=0; j<vol->YDim(); j++)
	  {
	    if(num_colors > 256)
	      {
		printf("Warning: more than 256 colors! Any color not in the list will result in a zero output voxel.\n");
		num_colors = 256;
		goto docalc;
	      }
		
	    color.r = slice[0][i+j*vol->XDim()];
	    color.g = slice[1][i+j*vol->XDim()];
	    color.b = slice[2][i+j*vol->XDim()];
		  
	    if(lsearch(&color,colors,&num_colors,sizeof(color_t),color_cmp) == NULL)
	      {
		printf("Error in lsearch()!\n");
		return 1;
	      }
	  }
		
      fprintf(stderr,"Determining color list... %5.2f %%   \r",(((float)k)/((float)((int)(vol->ZDim()-1))))*100.0);
    }
  printf("\n");
	
 docalc:;

  printf("Number of colors: %d\n",num_colors);
  printf("Colors: ");
  for(i=0; i<num_colors; i++)
    printf("(%d,%d,%d) ",colors[i].r,colors[i].g,colors[i].b);
  printf("\n");
  unsigned int range_size = 256/num_colors; /* range size == the whole space divided by the number of found colors */
  printf("Range size: %d\n",range_size);
	
  /* now write the output volume */
  for(k=0; k<vol->ZDim(); k++)
    {
      /* get the colors for this slice */
      vol->get(red_var,timestep)->getMapped(0,0,k,vol->XDim(),vol->YDim(),1,slice[0]);
      vol->get(green_var,timestep)->getMapped(0,0,k,vol->XDim(),vol->YDim(),1,slice[1]);
      vol->get(blue_var,timestep)->getMapped(0,0,k,vol->XDim(),vol->YDim(),1,slice[2]);
      vol->get(alpha_var,timestep)->getMapped(0,0,k,vol->XDim(),vol->YDim(),1,slice[3]);
	  
      /* lookup each color to determine it's output voxel value */
      /* check each colored voxel and add all new colors to the list to determine output voxel ranges */
      color_t color,*cur;
      unsigned int index, min, max;
      for(i=0; i<vol->XDim(); i++)
	for(j=0; j<vol->YDim(); j++)
	  {
	    color.r = slice[0][i+j*vol->XDim()];
	    color.g = slice[1][i+j*vol->XDim()];
	    color.b = slice[2][i+j*vol->XDim()];
		  
	    cur = (color_t *)lfind(&color,colors,&num_colors,sizeof(color_t),color_cmp);
	    if(cur == NULL)
	      {
		slice[4][i+j*vol->XDim()] = 0;
		continue;
	      }
	    index = ((unsigned int)(cur - colors)); /* determine the color's index */
	    min = index*range_size; /* find the start of this color's range */
	    max = min+range_size-1; /* find the end of this color's range */ /* Note: due to the discreet nature of unsigned char,
										we may not use the entire available 256 voxel values.
									     */
	    /* now use the color's alpha value to determine where on the range the output voxel is */
	    slice[4][i+j*vol->XDim()] = slice[3][i+j*vol->XDim()] == 0 ? 0 : (unsigned char)(min + float(range_size-1)*(float(slice[3][i+j*vol->XDim()])/255.0));
	  }
	
      fwrite(slice[4],sizeof(unsigned char),vol->XDim()*vol->YDim(),outvol);
	
      fprintf(stderr,"Writing output volume... %5.2f %%   \r",(((float)k)/((float)((int)(vol->ZDim()-1))))*100.0);
    }
  printf("\n");
	
  for(i=0; i<5; i++) free(slice[i]);
	
  return 0;
}
