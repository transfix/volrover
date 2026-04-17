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

#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

#include <sys/types.h>
#include <sys/stat.h>
#ifndef _WIN32
#include <unistd.h>
#endif
#include <errno.h>
#include <math.h>

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple_io.hpp>

#include <VolMagick/VolMagick.h>
#include <VolMagick/VolumeCache.h>
#include "VolMagick/endians.h"


#define _LITTLE_ENDIAN 1


void swap_buffer(char *buffer, int count, int typesize);

int main(int argc, char **argv)
{
  if(argc < 3)
    {
      std:: cerr << 
	"Usage: " << argv[0] << 

	
	"  number of volumes, first volume, second volume,..., last volume,  output volume. \n";

      return 1;
    }

  try
    {

	unsigned int xdim,ydim,zdim;
	int i,j,k;
  
	unsigned int	MagicNumW=0xBAADBEEF;
	unsigned int	NumTimeStep, NumVariable;
  	float		MinXYZT[4], MaxXYZT[4];
  	unsigned char	VariableType[100];
  	char		*VariableName[100];
	unsigned char c_unchar;


	int flag=0;

	size_t fwrite_return = 0;

	  int Classes = atoi(argv[1]);
	  FILE *fp = fopen(argv[argc-1],"w");


  	 std::vector<float> r(Classes), g(Classes), b(Classes);

      VolMagick::Volume inputVol;

  //    VolMagick::Volume inputVol2;
      std::vector <VolMagick::Volume> inputVols;

	  for(i=0; i< Classes; i++)
	  {
	  	VolMagick::readVolumeFile(inputVol, argv[i+2]);
		inputVols.push_back(inputVol);
	  }
      
      VolMagick::Volume outputVol;

           

      VolMagick::VolumeFileInfo volinfo1;
      volinfo1.read(argv[2]);
	//  std::cout << volinfo1.filename() << ":" <<std::endl;

	  NumTimeStep = 1;
	  NumVariable = 4;
	  MinXYZT[0] = volinfo1.boundingBox().minx;
	  MinXYZT[1] = volinfo1.boundingBox().miny;
	  MinXYZT[2] = volinfo1.boundingBox().minz;
	  MinXYZT[3] = 0.0;
	  MaxXYZT[0] = volinfo1.boundingBox().maxx;
	  MaxXYZT[1] = volinfo1.boundingBox().maxy;
	  MaxXYZT[2] = volinfo1.boundingBox().maxz;
	  MaxXYZT[3] = 1.0;

      VariableType[0] = 1;
	  VariableName[0] = (char*)malloc(sizeof(char)*64);
	  strcpy (VariableName[0], "red");
	  VariableType[1] = 1;
	  VariableName[1] = (char*)malloc(sizeof(char)*64);
	  strcpy (VariableName[1], "green");
	  VariableType[2] = 1;
	  VariableName[2] = (char*)malloc(sizeof(char)*64);
	  strcpy (VariableName[2], "blue");
	  VariableType[3] = 1;
	  VariableName[3] = (char*)malloc(sizeof(char)*64);
	  strcpy (VariableName[3], "alpha");

	  xdim = volinfo1.XDim();
	  ydim = volinfo1.YDim();
	  zdim = volinfo1.ZDim();

#ifdef _LITTLE_ENDIAN
	  swap_buffer((char *)&MagicNumW, 1, sizeof(unsigned int));
	  swap_buffer((char *)&xdim, 1, sizeof(unsigned int));
	  swap_buffer((char *)&ydim, 1, sizeof(unsigned int));
	  swap_buffer((char *)&zdim, 1, sizeof(unsigned int));
	  swap_buffer((char *)&NumTimeStep, 1, sizeof(unsigned int));
	  swap_buffer((char *)&NumVariable, 1, sizeof(unsigned int));
	  swap_buffer((char *)MinXYZT, 4, sizeof(float));
	  swap_buffer((char *)MaxXYZT, 4, sizeof(float));
#endif 
  
	  fwrite_return = fwrite(&MagicNumW, sizeof(unsigned int), 1, fp);
	  fwrite_return = fwrite(&xdim, sizeof(unsigned int), 1, fp);
	  fwrite_return = fwrite(&ydim, sizeof(unsigned int), 1, fp);
	  fwrite_return = fwrite(&zdim, sizeof(unsigned int), 1, fp);
	  fwrite_return = fwrite(&NumTimeStep, sizeof(unsigned int), 1, fp);
	  fwrite_return = fwrite(&NumVariable, sizeof(unsigned int), 1, fp);
	  fwrite_return = fwrite(MinXYZT, sizeof(float), 4, fp);
	  fwrite_return = fwrite(MaxXYZT, sizeof(float), 4, fp);

#ifdef _LITTLE_ENDIAN
  	  swap_buffer((char *)&MagicNumW, 1, sizeof(unsigned int));
	  swap_buffer((char *)&xdim, 1, sizeof(unsigned int));
	  swap_buffer((char *)&ydim, 1, sizeof(unsigned int));
	  swap_buffer((char *)&zdim, 1, sizeof(unsigned int));
	  swap_buffer((char *)&NumTimeStep, 1, sizeof(unsigned int));
	  swap_buffer((char *)&NumVariable, 1, sizeof(unsigned int));
	  swap_buffer((char *)MinXYZT, 4, sizeof(float));
	  swap_buffer((char *)MaxXYZT, 4, sizeof(float));
#endif 

     for (int m=0; m<4; m++) {
       fwrite_return = fwrite(&VariableType[m], sizeof(unsigned char), 1, fp);
 	   fwrite_return = fwrite(VariableName[m], sizeof(unsigned char), 64, fp);
     }
  
	int num = (int)(Classes/7);
	int remainder = Classes%7;
	std::map<unsigned int, std::vector<float> >  colors;
	std::map<unsigned int, std::vector<float> >::iterator it;

	float red, green, blue;
   if(Classes%7==0) 
   {
   		num--;
		remainder=7;
	}
   	
    int t=0;

	for(i = num; i >= 0; i --)
	{
		
	float invnum = (float)((i+1.0)/(num+1)*255);

	std::vector<float> color;

	if(i>0)
	{

	 for(j=0; j<7; j++)
	 {
	  t++;
	  	color.clear();

	 	switch(j){
		case 0:	
		red = invnum;
		green = 0;
		blue = 0;
		break;

		case 1:
		red = 0;
		green = invnum;
		blue = 0;
		break;

		case 2:
		red = 0;
		green = 0;
		blue = invnum;
		break;

		case 3:
		red = invnum;
		green = 0;
		blue = invnum;
		break;

		case 4:
		red = invnum;
		green = invnum;
		blue = 0;
		break;

		case 5:
		red = 0;
		green = invnum;
		blue = invnum;
		break;

		case 6:
		red = invnum;
		green = invnum;
		blue = invnum;
		break;
	}
	color.push_back(red);
	color.push_back(green);
	color.push_back(blue);
	colors[t]=color;
	}
   }
   else{
   for(j = 0; j<remainder; j++)
   {
   	t++;
	color.clear();
	if(j==0)
	 {
		red = invnum;
		green = 0;
		blue = 0;
	 }
	 else if (j==1)
	 {
		red = 0;
		green = invnum;
		blue = 0;
	 }
	 else if(j== 2)
	 {
		red = 0;
		green = 0;
		blue = invnum;
	 }
	 else if(j==3)
	 {
		red = invnum;
		green = 0;
		blue = invnum;
	}
	else if(j == 4)
	{
	 	red = invnum;
		green = invnum;
		blue = 0;
	}
	else if (j==5)
	{
		red = 0;
		green = invnum;
		blue = invnum;
	}
	else if (j==6)
	{
		red = invnum;
		green = invnum;
		blue = invnum;
		
	}
	
	color.push_back(red);
	color.push_back(green);
	color.push_back(blue);
	colors[t]=color;
	}
  }
  }
  i=0;
  for(it=colors.begin(); it != colors.end(); it++)
  {
  	r[i] = (*it).second.at(0);
  	g[i] = (*it).second.at(1);
  	b[i] = (*it).second.at(2);
//	std::cout<<"i= " << i<< "  r,g,b = " << r[i] <<", " << g[i] <<", " << b[i] << std::endl;
    i++;
  }


	  	for(k=0;k<zdim; k++)
	  	   for(j=0;j<ydim; j++)
	  	      for(i=0;i<xdim; i++)
			  {
				  flag = 0;
			  	  for(int classes=0; classes<Classes; classes++)
	 			 {

			  		 float temp = inputVols[classes](i,j,k);
					 if (temp > inputVols[classes].min())
					 {
					 	c_unchar = (unsigned char) r[classes];
					 	fwrite_return = fwrite(&c_unchar, sizeof(unsigned char), 1, fp);
						flag = 1;
						break;
					 }
				  }
				  if(flag ==0)
				  {
					 	c_unchar = 0;
					 	fwrite_return = fwrite(&c_unchar, sizeof(unsigned char), 1, fp);
				  }
	 		 }



	  	for(k=0;k<zdim; k++)
	  	   for(j=0;j<ydim; j++)
	  	      for(i=0;i<xdim; i++)
			  {
			     flag = 0;
			  	  for(int classes=0; classes<Classes; classes++)
	 			 {

			  		 float temp = inputVols[classes](i,j,k);
					 if (temp > inputVols[classes].min())
					 {
					 	c_unchar = (unsigned char) g[classes];
					 	fwrite_return = fwrite(&c_unchar, sizeof(unsigned char), 1, fp);
						flag = 1;
						break;
					 }
				 }
				 if(flag ==0)
				 {
					 	c_unchar = 0;
					 	fwrite_return = fwrite(&c_unchar, sizeof(unsigned char), 1, fp);
			 	 }
	 		 }




	  	for(k=0;k<zdim; k++)
	  	   for(j=0;j<ydim; j++)
	  	      for(i=0;i<xdim; i++)
			  {
			     flag = 0;
			  	  for(int classes=0; classes<Classes; classes++)
	 			 {

			  		 float temp = inputVols[classes](i,j,k);
					 if (temp > inputVols[classes].min())
					 {
					 	c_unchar = (unsigned char) b[classes];
					 	fwrite_return = fwrite(&c_unchar, sizeof(unsigned char), 1, fp);
						flag = 1;
						break;
					 }
				 }
				 if(flag == 0)
				 {
					 	c_unchar = 0;
					 	fwrite_return = fwrite(&c_unchar, sizeof(unsigned char), 1, fp);
				 }
	 		 }


	  	for(k=0;k<zdim; k++)
	  	   for(j=0;j<ydim; j++)
	  	      for(i=0;i<xdim; i++)
			  {
			     flag = 0;
			  	 for(int classes=0; classes<Classes; classes++)
	 			 {

			  		 float temp = inputVols[classes](i,j,k);
					 if (temp > inputVols[classes].min())
					 {
					    c_unchar = (unsigned char) temp;
					 	fwrite_return = fwrite(&c_unchar, sizeof(unsigned char), 1, fp);
						flag = 1;
						break;
					 }
				 }
				 if(flag == 0)
				 {
					 	c_unchar =  0;
					 	fwrite_return = fwrite(&c_unchar, sizeof(unsigned char), 1, fp);
				 }
			 	
	 		 }

	fclose(fp);
	std::cout<<"done" << std::endl;

     }

  catch(VolMagick::Exception &e)
    {
      std:: cerr << e.what() << std::endl;
    }
  catch(std::exception &e)
    {
      std::cerr << e.what() << std::endl;
    }

  return 0;
}



void swap_buffer(char *buffer, int count, int typesize)
{
  char sbuf[4];
  int i;
  int temp = 1;
  unsigned char* chartempf = (unsigned char*) &temp;
  if(chartempf[0] > '\0') {
  
	// swapping isn't necessary on single byte data
	if (typesize == 1)
		return;
  
  
	for (i=0; i < count; i++)
    {
		memcpy(sbuf, buffer+(i*typesize), typesize);
      
		switch (typesize)
		{
			case 2:
			{
				buffer[i*typesize] = sbuf[1];
				buffer[i*typesize+1] = sbuf[0];
				break;
			}
			case 4:
			{
				buffer[i*typesize] = sbuf[3];
				buffer[i*typesize+1] = sbuf[2];
				buffer[i*typesize+2] = sbuf[1];
				buffer[i*typesize+3] = sbuf[0];
				break;
			}
			default:
				break;
		}
    }

  }
}

