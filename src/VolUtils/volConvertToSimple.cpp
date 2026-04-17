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
      VolMagick::Volume inputVol;


      
      
      VolMagick::readVolumeFile(inputVol,argv[1]); ///first argument is input volume
      
	
	  VolMagick::Volume  outputVol;

	
      VolMagick::VolumeFileInfo volinfo1;
      volinfo1.read(argv[1]);
      std::cout << volinfo1.filename() << ":" <<std::endl;

      
      std::cout<<"minVol1 , maxVol1: "<<volinfo1.min()<<" "<<volinfo1.max()<<std::endl;;
 
	  float span[3];
	  span[0]=inputVol.XSpan(); 
	  span[1]=inputVol.YSpan(); 
	  span[2]=inputVol.ZSpan();
	
	  string ss=string(argv[2]);
	 
	  FILE* fp = fopen(ss.c_str(), "wb");
	  if(!fp){
	  	cout<<"Cannot write. "<<argv[2] << endl;
		return 0;
	  }
	
	  int xd = inputVol.XDim();
	  int yd = inputVol.YDim();
	  int zd = inputVol.ZDim();
	  float xmin = inputVol.XMin();
	  float ymin = inputVol.YMin();
	  float zmin = inputVol.ZMin();

	  float Minimum = inputVol.min();

	  fwrite(&xmin, sizeof(float), 1, fp);
	  fwrite(&ymin, sizeof(float), 1, fp);
	  fwrite(&zmin, sizeof(float), 1, fp);

		
	  fwrite(&xd, sizeof(int), 1, fp);
	  fwrite(&yd, sizeof(int), 1, fp);
	  fwrite(&zd, sizeof(int), 1, fp);

	  fwrite(&span[0], sizeof(float), 1, fp);
	  fwrite(&span[1], sizeof(float), 1, fp);
	  fwrite(&span[2], sizeof(float), 1, fp);


	float tmp;

	  for( int kz = 0; kz<inputVol.ZDim(); kz++)
	   for( int jy = 0; jy<inputVol.YDim(); jy++)
	     for( int ix = 0; ix<inputVol.XDim(); ix++)
		    {
				tmp = inputVol(kz,jy,ix)-Minimum;
		//	if(ix+jy+kz<100)	std::cout<<tmp<<" "<< std::endl;
				fwrite(&tmp, sizeof(float), 1, fp);
			}
    

	fclose(fp);

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
