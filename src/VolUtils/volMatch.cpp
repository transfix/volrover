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
#include <vector>
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

	
	"  first volume  second volume  output volume.  \n  first volume matches the size of second volume. first is rawv  \n";

      return 1;
    }

  try
    {
      vector<VolMagick::Volume> inputVol;

      VolMagick::Volume inputVol2;
      
      vector<VolMagick::Volume> outputVol;

      VolMagick::readVolumeFile(inputVol,argv[1]); ///first argument is input volume
      
      VolMagick::readVolumeFile(inputVol2, argv[2]); /// second argument is mask volume
      

      VolMagick::VolumeFileInfo volinfo1;
      volinfo1.read(argv[1]);
      std::cout << volinfo1.filename() << ":" <<std::endl;

      VolMagick::VolumeFileInfo volinfo2;
      volinfo2.read(argv[2]);

      std::cout << volinfo2.filename() << ":" <<std::endl;
      
      std::cout<<"minVol1 , maxVol1: "<<volinfo1.min()<<" "<<volinfo1.max()<<std::endl;;

      std::cout<<"minVol2 , maxVol2: "<<volinfo2.min()<<" "<<volinfo2.max()<<std::endl;

      VolMagick::Volume tmpVol;

      for(int i=0; i< inputVol.size(); i++)
	  {	
	  	tmpVol = inputVol[i];
		tmpVol.boundingBox(inputVol2.boundingBox());
	  	
		outputVol.push_back(tmpVol);
	 }
/*      outputVol.voxelType(inputVol.voxelType());
      outputVol.dimension(inputVol.dimension());
      outputVol.boundingBox(inputVol.boundingBox());
      
      
      std::cout<<"voxeltype "<<inputVol.voxelType()<<std::endl;
      	

      for(int i=0; i<inputVol2.size(); i++)
      for( int kz = 0; kz<inputVol2.ZDim(); kz++)
	    for( int jy = 0; jy<inputVol2.YDim(); jy++)
	      for( int ix = 0; ix<inputVol.XDim(); ix++)
		  		outputVol(ix, jy, kz, inputVol[i](ix, jy, kz));
					
		  

      VolMagick::createVolumeFile(outputVol, argv[3]);
*/
		VolMagick::writeVolumeFile(outputVol,argv[3]);



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
