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
  if(argc < 4)
    {
      std:: cerr << 
	"Usage: " << argv[0] << 

	
	"<first volume>  <second volume> <output volume> \n";
	std::cerr<<"second volume is subset of first volume" << std::endl;

      return 1;
    }

  try
    {
      VolMagick::Volume inputVol;

      VolMagick::Volume inputVol2;

      
      
      VolMagick::readVolumeFile(inputVol,argv[1]); ///first argument is input volume
      
      VolMagick::readVolumeFile(inputVol2, argv[2]); /// second argument is mask volume
      
	
	  VolMagick::Volume  outputVol;


      if((inputVol2.XSpan()- inputVol.XSpan()> 0.0001) || (inputVol2.YSpan()- inputVol.YSpan()>0.0001) || (inputVol2.ZSpan()- inputVol.ZSpan()>0.0001)
	  || inputVol2.XDim()>inputVol.XDim() || inputVol2.YDim() > inputVol.YDim() || inputVol2.ZDim()> inputVol.ZDim())
	  {
	  	std::cerr<<" inputVol2 is not a subset of inputVol1, you cannot use this program, sorry." << std::endl;
		std::cerr<<" inputVol2: Dims, spans: " << inputVol2.XDim() <<" " << inputVol2.XSpan() << " "<< inputVol2.YDim() <<" " << inputVol2.YSpan() << " "
<< inputVol2.ZDim() <<" " << inputVol2.ZSpan() << std::endl;
		std::cerr<<" inputVol: Dims, spans: " << inputVol.XDim() <<" " << inputVol.XSpan() << " "<< inputVol.YDim() <<" " << inputVol.YSpan() << " "
<< inputVol.ZDim() <<" " << inputVol.ZSpan() << std::endl;
	

		return 1;
	  }


      outputVol.voxelType(inputVol.voxelType());
      outputVol.dimension(inputVol.dimension());
      outputVol.boundingBox(inputVol.boundingBox());
      
      
      
	  int i, j, k;

	  int shiftx = (int)((inputVol2.XMin()-inputVol.XMin())/inputVol.XSpan()+0.5);
	  int shifty = (int)((inputVol2.YMin()-inputVol.YMin())/inputVol.YSpan()+0.5);
	  int shiftz = (int)((inputVol2.ZMin()-inputVol.ZMin())/inputVol.ZSpan()+0.5);
	  

	  float x, y, z;

 	  for( int kz = 0; kz<inputVol.ZDim(); kz++)
	   for( int jy = 0; jy<inputVol.YDim(); jy++)
	     for( int ix = 0; ix<inputVol.XDim(); ix++)
		   { 
		   		i = ix- shiftx;
				j = jy - shifty;
				k = kz - shiftz;
			
			if(i>=0 && j >=0 && k >=0 && i< inputVol2.XDim() && j < inputVol2.YDim() && k < inputVol2.ZDim() && (inputVol2(i,j,k)+inputVol2(maxi(i-1,0),j,k)
	+ inputVol2(mini(i+1, (int)(inputVol2.XDim())-1), j, k)  + inputVol2(i, maxi(j-1,0), k)+inputVol2(i, mini(j+1, (int)(inputVol2.YDim())-1), k) + inputVol2(i,j,maxi(k-1,0))+inputVol2(i,j,mini(k+1, (int)(inputVol2.ZDim())-1))) > inputVol2.min())
			{
				outputVol(ix, jy, kz, inputVol(ix, jy,kz));
			}
			else  outputVol(ix, jy, kz, inputVol.min());
		}	

/*			 if(inputVol2(ix, jy, kz)>= inputVol2.min())
			 {
//		     	i = ix + shiftx;
//				j = jy + shifty;
//				k = kz + shiftz;
			    
				x = inputVol2.XMin() + ix*inputVol2.XSpan();
				y = inputVol2.YMin() + jy*inputVol2.YSpan();
				z = inputVol2.ZMin() + kz*inputVol2.ZSpan();

			

				outputVol(ix, jy, kz, inputVol.interpolate(x,y,z));
			}
		//	else outputVol(ix, jy, kz, inputVol.min()); */
	//	 }


	VolMagick::createVolumeFile(outputVol, argv[3]);
  
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
