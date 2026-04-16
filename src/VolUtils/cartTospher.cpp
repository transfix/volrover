/*
  Copyright 2005-2011 The University of Texas at Austin

        Authors: Joe Rivera <transfix@ices.utexas.edu>
                 Ajay Gopinath <ajay@ices.utexas.edu>
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

/** Ajay Gopinath:
 * convert a volume file from Cartesian to spherical coordinate system
 */

#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

#include <boost/cstdint.hpp>

#include <VolMagick/VolMagick.h>
#include <VolMagick/VolumeCache.h>
#include <VolMagick/endians.h>

#define PI 3.14159265

using namespace std;

// sys/types.h is being included on Linux from one of the above headers...
#ifdef _WIN32
typedef unsigned int uint;
#endif

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

typedef boost::tuple<double, double, double> Color;

int main(int argc, char **argv)
{
 if(argc < 3)
    {
      cerr << "Usage: inputfile, outputfile \n";


      return 1;
    }

  try
    {
      VolMagick::Volume inputVol;

      VolMagick::readVolumeFile(inputVol,argv[1]); ///first argument is input volume



        VolMagick::Volume spher;

        VolMagick::Dimension dim;

	std::cout<<"Cart x: "<<inputVol.XDim()<<" y: "<<inputVol.YDim()<<" z: "<<inputVol.ZDim()<<"\n";


	dim.xdim =sqrt( (inputVol.XDim()+1.0)*(inputVol.XDim()+1.0)/4.0 + (inputVol.YDim()+1.0)*(inputVol.YDim()+1)/4 + (inputVol.ZDim()
			+1.0)*(inputVol.ZDim()+1.0)/4.0 );

         dim.ydim = 360;

         dim.zdim = 180;

        std::cout<<"spher x: "<<dim.xdim <<" y: "<<dim.ydim<<" z: "<<dim.zdim<<"\n";

        std::cout<<"inputVol x:"<<inputVol.XDim()<<" y: "<<inputVol.YDim()<<" z: "<<inputVol.ZDim()<<"\n";

        VolMagick::BoundingBox outputBB;

         outputBB.minx = 0;

         outputBB.maxx = dim.xdim - 1;

         outputBB.maxy = dim.ydim - 1;

         outputBB.miny = 0;

         outputBB.minz = 0;

         outputBB.maxz = dim.zdim - 1;

         VolMagick::VoxelType vox = VolMagick::Float;

        spher.voxelType(vox);

        spher.dimension(dim);

        spher.boundingBox(outputBB);

        std::cout<<"enter conversion zone\n";

        float max = -10000000;

        for(uint i = 0; i< spher.XDim(); i++)
         for(uint j = 0; j<spher.YDim(); j++)
          for(uint k=0; k<spher.ZDim(); k++)

                {

                int x = (int) (i*cos(j*PI/180)*sin(k*PI/180) + (inputVol.XDim()+1)/2);

                int y = (int) (i*sin(j*PI/180)*sin(k*PI/180) + (inputVol.YDim()+1)/2);

                int z = (int) (i*cos(k*PI/180)  + (inputVol.ZDim()+1)/2);

		if(x<inputVol.XDim() && x >0 && y < inputVol.YDim() && y>0 && z< inputVol.ZDim() && z>0)
                {
                try
                        {

                        spher(i,j,k, inputVol(x,y,z) * i*i * sin(k*PI/180) ); // Need the r^2 sin(theta) term while performing integration. Volume integral in polar coordinates
                        if(spher(i,j,k) > max)
                                max = spher(i,j,k);
                        }

                  catch(VolMagick::Exception &e)
                    {
                      cerr << e.what() << endl;

                        std::cout<<x<<" "<<y<<" "<<z<<"  ...\n  ";

                        std::cout<<i<<" "<<j<<" "<<k<<"    ";

                        break;

                    }
                 }



                }


        std::cout<<"done cart to sphere\n";


	
      VolMagick::createVolumeFile(argv[2], spher);

      VolMagick::writeVolumeFile(spher, argv[2]);




 return 0;
 
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
