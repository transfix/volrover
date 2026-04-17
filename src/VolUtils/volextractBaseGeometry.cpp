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

/* $Id: volinv.cpp 1481 2010-03-08 00:19:37Z transfix $ */

#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>


#include <sys/types.h>
#include <sys/stat.h>
#ifndef _WIN32
#include <unistd.h>
#endif
#include <errno.h>
#include <math.h>

#include <VolMagick/VolMagick.h>
#include <VolMagick/VolumeCache.h>
#include <VolMagick/endians.h>

#include <cvcraw_geometry/cvcraw_geometry.h>


#include <fstream>

#define Pi 3.1415926

using namespace std;

typedef cvcraw_geometry::geometry_t geometry;
//typedef cvcraw_geometry::geometry_t::

class VolMagickOpStatus : public VolMagick::VoxelOperationStatusMessenger
{
public:
  void start(const VolMagick::Voxels *vox, Operation op, VolMagick::uint64 numSteps) const
  {
    _numSteps = numSteps;
  }

  void step(const VolMagick::Voxels *vox, Operation op, VolMagick::uint64 curStep) const
  {
    const char *opStrings[] = { "CalculatingMinMax", "CalculatingMin", "CalculatingMax",
				"SubvolumeExtraction", "Fill", "Map", "Resize", "Composite",
				"BilateralFilter", "ContrastEnhancement"};

    fprintf(stderr,"%s: %5.2f %%\r",opStrings[op],(((float)curStep)/((float)((int)(_numSteps-1))))*100.0);
  }

  void end(const VolMagick::Voxels *vox, Operation op) const
  {
    printf("\n");
  }

private:
  mutable VolMagick::uint64 _numSteps;
};

int main(int argc, char **argv)
{
  if(argc != 5)
    {
      cerr << "Usage: " << argv[0] << " <input volume file> <input geometry file> <output volume file> <int thresh>" << endl;
      return 1;
    }

  try
    {
      VolMagickOpStatus status;
      VolMagick::setDefaultMessenger(&status);

      VolMagick::Volume inputVol;
	
      VolMagick::readVolumeFile(inputVol, argv[1]);
//	  inputVol2info.read(argv[2]);
//	  fstream fs;
//	  fs.open(argv[3]);

//	  int seeda, seedb, seedc;

     VolMagick::Volume outputVol;
	  
	 outputVol.voxelType(inputVol.voxelType());
	 outputVol.dimension(inputVol.dimension());
	 outputVol.boundingBox(inputVol.boundingBox());


	 

	geometry geom;
	cvcraw_geometry::read(geom, string(argv[2]));


	
	for(int i=0; i <inputVol.XDim(); i++)
		  for(int j=0; j< inputVol.YDim(); j++)
			    for(int k=0; k<inputVol.ZDim(); k++)
				outputVol(i,j,k,inputVol.min());


	float x, y, z;

	float minx = inputVol.XMin();
	float miny = inputVol.YMin();
	float minz = inputVol.ZMin();

	float xspan = inputVol.XSpan();
	float yspan = inputVol.YSpan();
	float zspan = inputVol.ZSpan();


	int nx1, nx2, ny1, ny2, nz1, nz2;

	int thres=atoi(argv[4]);

	for(int i=0; i< geom.tris.size(); i++)
	{
	  for(int j=0; j< 3; j ++)
	  {
		x = geom.points[geom.tris[i][j]][0];
		y = geom.points[geom.tris[i][j]][1];
		z = geom.points[geom.tris[i][j]][2];

		nx1 = (int)floor((x-minx)/xspan);
		nx2 = (int)ceil((x-minx)/xspan);
		ny1 = (int)floor((y-miny)/yspan);
		ny2 = (int)ceil((y-miny)/yspan);
		nz1 = (int)floor((z-minz)/zspan);
		nz2 = (int)ceil((z-minz)/zspan);

		for(int a = nx1-thres; a <= nx2 + thres; a ++)
			for(int b = ny1 - thres; b <= ny2 + thres; b ++)
				for(int c = nz1- thres; c <= nz2 + thres; c++)
				{
					if(!(a<0) && !(a>=inputVol.XDim()) 
					 && !(b<0) && !(b>=inputVol.YDim()) 
					 && !(c<0) && !(c>=inputVol.ZDim()) ) 
						outputVol(a,b,c, inputVol(a,b,c));
				}
		}
	}

	VolMagick::createVolumeFile(outputVol, argv[3]);

	cout<<"done !" <<endl;

	}
  catch(VolMagick::Exception &e)
    {
      cerr << e.what() << endl;
    }
  catch(std::exception &e)
    {
      cerr << e.what() << endl;
    }

  return 0;
}
