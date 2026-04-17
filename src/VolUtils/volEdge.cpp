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

/* AJAY GOPINATH: implementing the canny edge detector in 3D
 * Laplacian of Gaussian operator will be used.
*/


#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

#define _USE_MATH_DEFINES
#include <cmath>

#include <VolMagick/VolMagick.h>
#include <VolMagick/VolumeCache.h>
#include <VolMagick/endians.h>

typedef unsigned int uint;

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

typedef boost::tuple<double, double, double> Color;

int main(int argc, char **argv)
{
    typedef VolMagick::uint64 uint;

  if(argc < 4)
    {
      cerr << "Usage: inputfile, outputfile, sigma \n";


      return 1;
    }

  try
    {
      VolMagick::Volume inputVol;



      VolMagick::Volume outputVol;   


      VolMagick::Volume volLoG;


      VolMagick::readVolumeFile(inputVol,argv[1]); ///first argument is input volume


      outputVol.voxelType(inputVol.voxelType());
      outputVol.dimension(inputVol.dimension());
      outputVol.boundingBox(inputVol.boundingBox());

      volLoG.voxelType(inputVol.voxelType());
      volLoG.dimension(inputVol.dimension());
      volLoG.boundingBox(inputVol.boundingBox());

  
     
     float sigma = atof(argv[3]);


 //  Creating the LoG kernel operator

	double LoG[27]; 
	
	int count = 0;

//     for(int i = -2; i<3; i++)
//	 for(int j = -2; j<3; j++)
//		for(int k = -2; k<3; k++)

	{
		int i = -2;
		int j = -2;
        float term1 = 1/(pow(2*M_PI, 1.5) * pow(sigma, 5));

		int k = -2;


		double term2 = (i*i + j*j + k*k)/(sigma*sigma) - 3;

		double term3 = exp( (-i*i-j*j-k*k)/(2*sigma*sigma) );

		LoG[count] = term1*term2*term3;


		std::cout<<term1<<" "<<term2<<" "<<term3<<"\n";

		std::cout<<LoG[count]<<" ";
	
		count++;	


	}


 // Convolving LoG operator with input volume

	for(uint x=0; x<volLoG.XDim(); x++)
	 for(uint y=0; y<volLoG.YDim(); y++)
	  for(uint z=0; z<volLoG.ZDim(); z++)

	{

	float val = 0;

	uint count = 0;

		 for(int i = -2; i<3; i++)
                   for(int j = -2; j<3; j++)
	       	    for(int k = -2; k<3; k++)
			{

			if( (x+i)>0 && (x+i)<inputVol.XDim() && (y+j)>0 && (y+j)<inputVol.YDim() && (z+k)>0 &&(z+k)<inputVol.ZDim() )
	
				val = val + inputVol(x+i, y+j, z+k)*LoG[count];
	
			count++;
			

			}


	volLoG(x,y,z, val);
	
	
	
	}



/// Non-maximum supression

	

      VolMagick::createVolumeFile("./volLoG.rawiv", volLoG);

      VolMagick::writeVolumeFile(volLoG, "./volLoG.rawiv");



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
