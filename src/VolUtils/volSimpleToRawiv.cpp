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
#include <VolMagick/endians.h>



using namespace std;
int main(int argc, char **argv)
{
  if(argc < 3)
    {
      std:: cerr << 
	"Usage: " << argv[0] << 

	
	" <first volume>  <output simple volume> \n";

      return 1;
    }

  try
    {
      
	   string ss=string(argv[1]);
	   FILE* fp;
	   if ((fp= fopen(ss.c_str(), "rb"))==NULL)
	   {
	   		cout<<"read error ... \n";
			exit(0);
	   }

	  float xmin, ymin, zmin;
	  float xmax, ymax, zmax;
	  fread(&xmin, sizeof(float), 1, fp);
	  fread(&ymin, sizeof(float), 1, fp);
	  fread(&zmin, sizeof(float), 1, fp);

	  unsigned int dim[3];
	  fread(dim,sizeof(unsigned int),3, fp);

	  float span[3], origin[3];
	  fread(span, sizeof(float), 3, fp);

	  origin[0] = xmin;
	  origin[1] = ymin;
	  origin[2] = zmin;

	  xmax = xmin + span[0]*(dim[0]-1);
	  ymax = ymin + span[1]*(dim[1]-1);
	  zmax = zmin + span[2]*(dim[2]-1);

      unsigned int numverts = dim[0]*dim[1]*dim[2];
	  unsigned int numcells = (dim[0]-1)*(dim[1]-1)*(dim[2]-1);

	  unsigned int d0, d1, d2;
	  d0 = dim[0];
	  d1 = dim[1];
	  d2 = dim[2];

      if(!big_endian())
	  {
	  	SWAP_32(&xmin);
		SWAP_32(&ymin);
		SWAP_32(&zmin);

		SWAP_32(&xmax);
		SWAP_32(&ymax);
		SWAP_32(&zmax);

		SWAP_32(&numverts);
		SWAP_32(&numcells);
		
		for(int i=0; i<3; i++) 
		{
			SWAP_32(&(dim[i]));
			SWAP_32(&(origin[i]));
			SWAP_32(&(span[i]));
		}
	  } 
	  string out = string(argv[2]);
	  FILE* fout = fopen(out.c_str(), "wb");

	  fwrite(&xmin, sizeof(float), 1, fout); 
	  fwrite(&ymin, sizeof(float), 1, fout); 
	  fwrite(&zmin, sizeof(float), 1, fout); 
	  fwrite(&xmax, sizeof(float), 1, fout); 
	  fwrite(&ymax, sizeof(float), 1, fout); 
	  fwrite(&zmax, sizeof(float), 1, fout); 

	  fwrite(&numverts, sizeof(unsigned int), 1, fout);
	  fwrite(&numcells, sizeof(unsigned int), 1, fout);

	  for(int i=0; i<3; i++)
	  {	
	  	fwrite(&dim[i], sizeof(unsigned int), 1, fout);
	  }


	  for(int i=0; i<3; i++)
	  {	
	  	fwrite(&origin[i], sizeof(float), 1, fout);
	  }

	  for(int i=0; i<3; i++)
	  {	
	  	fwrite(&span[i], sizeof(float), 1, fout);
	  }

      
		
	float tmp, tmp1;

	 std::cout<<d0 <<" " << d1 <<" "<<d2 << std::endl;
     for( int ix = 0; ix<d0*d1*d2; ix++)
		    {
				fread(&tmp, sizeof(float), 1, fp);
		 //	std::cout<<tmp<< std::endl;
			tmp1 = tmp;
		if(!big_endian())		SWAP_32(&tmp1);
		//	if(ix+jy+kz<100)	std::cout<<tmp<<" "<< std::endl;
				fwrite(&tmp1, sizeof(float), 1, fout);
			}
    

	fclose(fp);
	fclose(fout);

    std::cout<<"done!"<<std::endl;


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
