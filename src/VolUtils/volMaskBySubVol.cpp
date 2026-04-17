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
#include <CVC/BoundingBox.h>


template <class T>
T mini(T a, T b)
{
	if(a<=b) return a;
	else return b;
}

template <class T>
T maxi(T a, T b)
{
	if(a>=b) return a;
	else return b;
}


using namespace std;
int main(int argc, char **argv)
{
  if(argc < 6)
    {
      std:: cerr << 
	"Usage: " << argv[0] << 

	
	"  <first volume>  <second volume> < threshold> <output volume left> <output volume masked> \n";
	std::cerr<<"second volume (son) is subset of first volume (mother), threshold is in [0,255], output masked is the volume in first volume with > threshold in the second volume. output volume left is the first volume - output volume masked. " << std::endl;

      return 1;
    }

  try
    {
      VolMagick::Volume inputVol;

      VolMagick::Volume inputVol2;

      
      
      VolMagick::readVolumeFile(inputVol,argv[1]); ///first argument is input volume
      
      VolMagick::readVolumeFile(inputVol2, argv[2]); /// second argument is mask volume
      
	  float threshold = atof(argv[3]);


	  VolMagick::Volume  outputLeft;
	  VolMagick::Volume  outputMasked;

      if(!inputVol2.boundingBox().isWithin(inputVol.boundingBox()) )
	  	throw VolMagick::SubVolumeOutOfBounds("The mask volume bounding box must be within the mother volume's bounding box.");



      outputLeft.voxelType(inputVol.voxelType());
      outputLeft.dimension(inputVol.dimension());
      outputLeft.boundingBox(inputVol.boundingBox());
      
      outputMasked.voxelType(inputVol.voxelType());
      outputMasked.dimension(inputVol.dimension());
      outputMasked.boundingBox(inputVol.boundingBox());
    
   

	  float x, y, z;
	  float fvalue;
	  float minvalue = (float)inputVol2.min();
	  float maxvalue = (float)inputVol2.max();

 	  for( int kz = 0; kz<inputVol.ZDim(); kz++)
	  {
	  	for( int jy = 0; jy<inputVol.YDim(); jy++)
	    {
			for( int ix = 0; ix<inputVol.XDim(); ix++)
		   	{
		     	z = inputVol.ZMin() + kz*inputVol.ZSpan();
			 	y = inputVol.YMin() + jy*inputVol.YSpan();
				x = inputVol.XMin() + ix*inputVol.XSpan();
				if((z >= (float)inputVol2.ZMin()) && (z <= (float)inputVol2.ZMax()) && ( y >=(float)inputVol2.YMin()) && (y <= (float)inputVol2.YMax()) &&(x >= (float)inputVol2.XMin()) && (x <= (float)inputVol2.XMax()) )
				{
				  fvalue = 255.0*(inputVol2.interpolate(x,y,z)-minvalue)/(maxvalue-minvalue);
				  if(fvalue > threshold)
				  {
				  	outputMasked(ix,jy,kz, inputVol(ix,jy,kz));
					outputLeft(ix,jy,kz, inputVol.min());
				  }	  
				  else
				  {
				  	outputMasked(ix,jy,kz, inputVol.min());
					outputLeft(ix,jy,kz, inputVol(ix,jy,kz));
				  }
				}
				else
				{
					outputMasked(ix,jy,kz, inputVol.min());
					outputLeft(ix,jy,kz, inputVol(ix,jy,kz));
				}
		     }		   
		  }	
		}



	VolMagick::createVolumeFile(outputLeft, argv[4]);
	VolMagick::createVolumeFile(outputMasked, argv[5]);

  
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
