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

/* $Id: volresize.cpp 1481 2010-03-08 00:19:37Z transfix $ */

#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#ifndef _WIN32
#include <unistd.h>
#endif
#include <errno.h>
#include <math.h>

#include <boost/cstdint.hpp>

#include <VolMagick/VolMagick.h>
#include <VolMagick/VolumeCache.h>
#include <VolMagick/endians.h>

#include <fstream>

typedef boost::uint32_t uint;

using namespace std;

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
  if(argc < 6)
    {
      cerr << "Usage: " << argv[0] << " <input volume file> <output volume file> <resample factor X dim> <resample factor Y dim> <resample factor Z dim> " << endl;
      return 1;
    }

  try
    {


      VolMagick::VolumeFileInfo volinfo;

      volinfo.read(argv[1]);

      VolMagick::Dimension targetDim( floor(volinfo.XDim()/atoi(argv[3])) +1 ,
					 floor(volinfo.YDim()/atoi(argv[4])) +1,
			  		  floor(volinfo.ZDim()/atoi(argv[5])) +1);


	std::cout<<"Read input Vol\n";

	uint xfact = atoi(argv[3]);
	uint yfact = atoi(argv[4]);
	uint zfact = atoi(argv[5]);

     VolMagick::BoundingBox bbox;
	bbox.minx = 0;
	bbox.maxx = targetDim.xdim-1;
	bbox.miny = 0;
	bbox.maxy = targetDim.ydim-1;
	bbox.minz = 0;
	bbox.maxz = targetDim.zdim-1;


	VolMagick::Volume inputVol;

	VolMagick::readVolumeFile(inputVol, argv[1]);
	
	VolMagick::Volume outputVol;

	outputVol.boundingBox(bbox);

	outputVol.dimension(targetDim);

	outputVol.voxelType(inputVol.voxelType());



      VolMagick::createVolumeFile(argv[2],
				  bbox,
				  targetDim,
				  volinfo.voxelTypes(),
				  volinfo.numVariables(),
				  volinfo.numTimesteps(),
				  volinfo.TMin(),volinfo.TMax());

	std::cout<<"Created outputFile\n";

	for(uint i =0, x=0; i<inputVol.XDim(); i=i+atoi(argv[3]), x++)
	 for(uint j=0, y=0; j<inputVol.YDim(); j=j+atoi(argv[4]), y++)
	  for(uint k=0, z=0; k<inputVol.ZDim(); k=k+atoi(argv[5]), z++)
		{
	
		

		outputVol(x,y,z, inputVol(i,j,k));


		}


			

		std::cout<<"writing outputFile\n";

	    writeVolumeFile(outputVol,argv[2]);
	
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
