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

/* $Id$ */

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

	
	"  <first volume>  <output volume>.   \n";

      return 1;
    }

  try
    {
      VolMagick::Volume inputVol;

      
      VolMagick::Volume outputVol;

      VolMagick::readVolumeFile(inputVol,argv[1]); ///first argument is input volume
      
      VolMagick::VolumeFileInfo volinfo1;
      volinfo1.read(argv[1]);
      std::cout << volinfo1.filename() << ":" <<std::endl;

       
      std::cout<<"minVol1 , maxVol1: "<<volinfo1.min()<<" "<<volinfo1.max()<<std::endl;;

  
      VolMagick::BoundingBox bbox;
	  bbox.minx = inputVol.XMin();
	  bbox.maxx = inputVol.XMax();
	  bbox.miny = inputVol.YMin();
	  bbox.maxy = inputVol.YMax();
	  bbox.minz = 2.0*inputVol.ZMin() - inputVol.ZMax();
	  bbox.maxz = inputVol.ZMax();

	  VolMagick::Dimension dim;
	  dim.xdim = inputVol.XDim();
	  dim.ydim = inputVol.YDim();
	  dim.zdim = (int)((bbox.maxz -bbox.minz)/inputVol.ZSpan())+1;


//	 cout<<bbox.minz <<" " << bbox.maxz<<" "<< bbox.maxy <<endl;
//	 cout<<dim.zdim <<" " << dim.ydim << endl;
      
      outputVol.voxelType(inputVol.voxelType());
      outputVol.dimension(dim);
      outputVol.boundingBox(bbox);

	
      
   
      
      for( int kz = 0; kz<outputVol.ZDim()/2; kz++)
	   for( int jy = 0; jy<inputVol.YDim(); jy++)
	      for( int ix = 0; ix<inputVol.XDim(); ix++)
  				outputVol(ix, jy, kz, inputVol(inputVol.XDim()-1-ix, inputVol.YDim()-1-jy, inputVol.ZDim()-1-kz));

      for( int kz=outputVol.ZDim()/2; kz < outputVol.ZDim(); kz++)
	   for( int jy = 0; jy<inputVol.YDim(); jy++)
	      for( int ix = 0; ix<inputVol.XDim(); ix++)
  				outputVol(ix, jy, kz, inputVol(ix, jy, kz+1-inputVol.ZDim()));
	  
	  
					
		  
			  
		 
	      

      VolMagick::createVolumeFile(outputVol, argv[2]);


      cout<<"done!"<<endl;

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
