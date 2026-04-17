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
  if(argc < 6)
    {
      std:: cerr << 
	"Usage: " << argv[0] << 

	
	" <input volume>  <shift xmin> <shift ymin> <shift zmin>  <output volume> \n";

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
	
	  VolMagick::BoundingBox bbox;
	  bbox.minx = atof(argv[2])+ volinfo1.XMin();
	  bbox.miny = atof(argv[3])+ volinfo1.YMin();
	  bbox.minz = atof(argv[4])+ volinfo1.ZMin();
	  bbox.maxx = atof(argv[2])+ volinfo1.XMax();
	  bbox.maxy = atof(argv[3])+ volinfo1.YMax();
	  bbox.maxz = atof(argv[4])+ volinfo1.ZMax();
	
/*	  VolMagick::Dimension dim;
	  dim.xdim = (int) ((bbox.maxx-bbox.minx)/ span[0])+1;
	  dim.ydim = (int) ((bbox.maxy-bbox.miny)/ span[1])+1;
	  dim.zdim = (int) ((bbox.maxz-bbox.minz)/ span[2])+1;

*/
 	
      outputVol.voxelType(inputVol.voxelType());
      outputVol.dimension(inputVol.dimension());
      outputVol.boundingBox(bbox);
      
 /*     
	 float temp;
	 int i, j, k;
     int xsh = (int)((bbox.minx - inputVol.XMin())/span[0]);
     int ysh = (int)((bbox.miny - inputVol.YMin())/span[1]);
     int zsh = (int)((bbox.minz - inputVol.ZMin())/span[2]);
	 
*/
	  for( int kz = 0; kz<outputVol.ZDim(); kz++)
	   for( int jy = 0; jy<outputVol.YDim(); jy++)
	     for( int ix = 0; ix<outputVol.XDim(); ix++)
		    {
		/*		i = ix + xsh;
				j = jy + ysh;
				k = kz + zsh;
		//		cout <<"ijk "<< i <<" " << j <<" " <<k <<endl;
				if(i<0 || i >= inputVol.XDim()|| j<0 || j>= inputVol.YDim() || k <0 || k>=inputVol.ZDim()) outputVol(ix, jy, kz, 0.0);
				else */
				outputVol(ix,jy,kz, inputVol(ix,jy,kz)  );
			}
	//		stringstream s;
	//		s<<argv[3]<< num <<".rawiv";
	      VolMagick::createVolumeFile(outputVol, argv[5]);
    

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
