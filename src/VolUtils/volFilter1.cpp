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




int main(int argc, char **argv)
{
  if(argc < 5)
    {
      std:: cerr << 
	"Usage: " << argv[0] <<"  inputfile   t_low   t_high   outputFile \n";
      return 1;
    }

  try
    {
      VolMagick::Volume inputVol;

      VolMagick::Volume inputVol2;
      
      VolMagick::Volume outputVol;

      VolMagick::readVolumeFile(inputVol,argv[1]); ///first argument is input volume
      
      

      VolMagick::VolumeFileInfo volinfo1;
      volinfo1.read(argv[1]);
      std::cout << volinfo1.filename() << ":" <<std::endl;

      
      outputVol.voxelType(inputVol.voxelType());
      outputVol.dimension(inputVol.dimension());
      outputVol.boundingBox(inputVol.boundingBox());
      
      
      std::cout<<"voxeltype "<<inputVol.voxelType()<<std::endl;
      	
     float tlow = atof(argv[2]);
     float thigh = atof(argv[3]); 
     
	 if(tlow < 0)
	 {
	 	tlow = 0;
		std::cout<<"tlow should be bigger than 0, set to 0. "<<std::endl;
	}
	 if(thigh>255) 
	 {
	 	thigh = 255;
		std::cout<<"thigh should be smaller than 255, set to 255. " <<std::endl;
	 }

      for( int kz = 0; kz<inputVol.ZDim(); kz++)
	 {
	    for( int jy = 0; jy<inputVol.YDim(); jy++)
	      for( int ix = 0; ix<inputVol.XDim(); ix++)
		{
		  float temp = inputVol(ix, jy, kz);
		  float temp1 = 255.0*(temp-volinfo1.min())/(volinfo1.max()-volinfo1.min());
		  if( (temp1 >= tlow) && (temp1 < thigh))
		  	outputVol(ix, jy, kz, temp);
		  else outputVol(ix, jy, kz, volinfo1.min());

		}
	}	
      

      VolMagick::createVolumeFile(outputVol, argv[argc-1]);


	std::cout<<" done!" <<std::endl;


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
