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
  if(argc < 4)
    {
      std:: cerr << 
	"Usage: " << argv[0] << 

	
	"<first volume>  <second volume> <output volume>  \n";

      return 1;
    }

  try
    {
      VolMagick::Volume inputVol;

      VolMagick::Volume inputVol2;

      
      
      VolMagick::readVolumeFile(inputVol,argv[1]); ///first argument is input volume
      
      VolMagick::readVolumeFile(inputVol2, argv[2]); /// second argument is mask volume
      
	 // int number = atoi(argv[4]);
	
	  VolMagick::Volume  outputVol;

	
      VolMagick::VolumeFileInfo volinfo1;
      volinfo1.read(argv[1]);
      std::cout << volinfo1.filename() << ":" <<std::endl;

      VolMagick::VolumeFileInfo volinfo2;
      volinfo2.read(argv[2]);
      std::cout << volinfo2.filename() << ":" <<std::endl;
      
      std::cout<<"minVol1 , maxVol1: "<<volinfo1.min()<<" "<<volinfo1.max()<<std::endl;;

      std::cout<<"minVol2 , maxVol2: "<<volinfo2.min()<<" "<<volinfo2.max()<<std::endl;

      
      outputVol.voxelType(inputVol.voxelType());
      outputVol.dimension(inputVol.dimension());
      outputVol.boundingBox(inputVol.boundingBox());
      
      

      
      
	//for (int num = 1; num < number; num ++)
    {
	  //std::cout<<"\xd"<<(float)(num)/(float)(number-1)*100.00<<"% is finished";
	  for( int kz = 0; kz<inputVol.ZDim(); kz++)
	   for( int jy = 0; jy<inputVol.YDim(); jy++)
	     for( int ix = 0; ix<inputVol.XDim(); ix++)
		   { 
		   float temp= inputVol2(ix,jy,kz)/(inputVol2(ix,jy,kz)-inputVol(ix, jy,kz));
		 //  if ( temp<-10.0 || temp> 1.0) temp = 0.0;
//		   	if (fabs(inputVol2(ix,jy,kz))>10 || fabs(inputVol(ix,jy,kz))>10) 
//			{	
//				outputVol(ix,jy,kz, -2); 
//				continue;
//			}
		    outputVol(ix,jy,kz, temp  );
			}
			stringstream s;
	 		s<<argv[3]<< ".rawiv";
	      VolMagick::createVolumeFile(outputVol, s.str());
    }

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
