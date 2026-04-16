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

#include <boost/cstdint.hpp>

#include <VolMagick/VolMagick.h>

#include <omp.h>
#include <time.h>

using namespace std;

#ifdef _WIN32
typedef unsigned int uint;
#endif

void computeLocalStats(const VolMagick::Volume currVol, uint *currentLoc, float &mean, float &var, const VolMagick::Volume gaussVol, int gaussianSize);

void computeCovar(const VolMagick::Volume inputVol, const VolMagick::Volume inputVol2, uint *currLoc, float &covar, float mean1, float mean2, const VolMagick::Volume gaussVol, int gaussianSize);


int main(int argc, char **argv)
{

if(argc<4)

  {
	cerr<< "Usage: inputvolume1, inputvolume2, outputSimilarityfile"<<"\n Given two input volumes, a structure similarity metric is generated based on the SSIM paper by Wang et al. 04 IEEE Transactions on Image Proc \n";
	return 1;	
 }

  try
    {
      VolMagick::Volume inputVol;

	  VolMagick::Volume inputVol2;
   
      VolMagick::readVolumeFile(inputVol,argv[1]); ///first argument is first input volume
	  
	  VolMagick::readVolumeFile(inputVol2,argv[2]); ///second argument is second input volume

	std::cout<<"Volume read \n Dimensions :"<<inputVol.XDim()<< " " << inputVol.YDim()<< " " <<inputVol.ZDim()<< " \n";


      VolMagick::VolumeFileInfo volinfo;

      volinfo.read(argv[1]);


      VolMagick::createVolumeFile(argv[3],
                                  volinfo.boundingBox(),
                                  volinfo.dimension(),
                                  volinfo.voxelTypes(),
                                  volinfo.numVariables(),
                                  volinfo.numTimesteps(),
                                  volinfo.TMin(),volinfo.TMax());



	VolMagick::Volume outputVol;   

      outputVol.voxelType(inputVol.voxelType());

      outputVol.dimension(inputVol.dimension());

      outputVol.boundingBox(inputVol.boundingBox());



     ///check dimensions are equal
	 
	 if(inputVol.XDim() != inputVol2.XDim() || inputVol.YDim() != inputVol2.YDim() || inputVol.ZDim() != inputVol2.ZDim())
		{
			std::cerr<<"Error, input volume sizes should be the same";
			return 1;
		}
		

	int numProcs = omp_get_num_procs();
	std::cout<<"Num of processors: "<< numProcs<<"\n";

	int startTime;
	startTime = clock();
	
	///Parameter constants as described in Wang04:
	//"Throughout this paper, the SSIM measure
	//uses the following parameter settings:K1 = 0.01 and K2 = 0.03
	//. These values are somewhat arbitrary, but we find that in
	//our current experiments, the performance of the SSIM index algorithm
	//is fairly insensitive to variations of these values."

	float K1 = 0.01;
	
	float K2 = 0.01;

	float L = 255; // Intensity range
	
	float C1 = (K1*L)*(K1*L);
	
	float C2 = (K2*L)*(K2*L);

	std::cout<<"Parameters C1: "<<C1<<"  C2: "<<C2<<"\n";


	/// 3D Gaussian kernel is modeled as a volmagic volume type
	//

  	VolMagick::Volume gaussVol;
	
	int gaussianSize = 11;
	
	float sigma  = 1;
	
	VolMagick::Dimension dim;

      dim.xdim = gaussianSize;
      dim.ydim = gaussianSize;
      dim.zdim = gaussianSize;


      VolMagick::BoundingBox bbox;

      bbox.minx = 0;
      bbox.maxx = dim.xdim;

      bbox.miny = 0;
      bbox.maxy = dim.ydim;	

      bbox.minz = 0;
      bbox.maxz = dim.zdim;

//      VolMagick::VoxelType gaussVoxType = VolMagick::VoxelType::Float;

      gaussVol.voxelType(VolMagick::Float);
      gaussVol.dimension(dim);
      gaussVol.boundingBox(bbox);

	float gaussSum = 0;

	for (uint i = 0; i< gaussVol.XDim(); i++)
	 for(uint j = 0; j<gaussVol.YDim(); j++)
	  for(uint k = 0; k<gaussVol.ZDim(); k++)
	   
	    {
	float PI = 3.14159;

	float denom = 1/(2*sqrt(PI)*sigma);
	
	float x = i - (gaussianSize - 1)/2;
	float y = j - (gaussianSize - 1)/2;
	float z = k - (gaussianSize - 1)/2;

	float numer = exp( (-x*x -y*y - z*z)/(sigma*sigma) );
 	
	float gaussVal = numer/denom;

	gaussSum = gaussSum + gaussVal;

	gaussVol(i,j,k, gaussVal);

	    }

	//Normalizing the gaussian window
	
		
	  for (uint i = 0; i< gaussVol.XDim(); i++)
         for(uint j = 0; j<gaussVol.YDim(); j++)
          for(uint k = 0; k<gaussVol.ZDim(); k++)

 	{

	gaussVol(i,j,k, gaussVol(i,j,k)/gaussSum);
	}

	float sum = 0;

	uint voxcount = 0;

#pragma omp parallel for schedule(static, inputVol.XDim()/numProcs)	      
		  
    for(boost::int64_t i = 0; i<inputVol.XDim(); i++)
	{ std::cerr<<i<<"..";
	
		for(uint j = 0; j<inputVol.YDim(); j++)
			for(uint k = 0; k<inputVol.ZDim(); k++)
			
			{
			
			uint currLoc[3];
			
			currLoc[0] = i; currLoc[1] = j; currLoc[2] = k;
			
			
			float mean1, mean2;
			
			float var1, var2;
	//		std::cout<<"currloc "<<i<<" "<<j<<" "<<k<<" ";
			
			computeLocalStats(inputVol, currLoc, mean1, var1, gaussVol, gaussianSize);
			
			computeLocalStats(inputVol2, currLoc, mean2, var2, gaussVol, gaussianSize);
			
			float covar;
			
			computeCovar(inputVol, inputVol2, currLoc, covar, mean1, mean2, gaussVol, gaussianSize);
			
			float SSIM = (2*mean1*mean2 + C1)*(2*covar + C2) / ( (mean1*mean1 + mean2*mean2 + C1)*( var1 + var2 + C2));

			sum = sum + SSIM;

			voxcount = voxcount + 1;

			outputVol(i,j,k, SSIM);
			
			
			}
	}
	
	int endTime;
	
	endTime = clock();



	std::cout<<"Time elapsed: "<<(endTime-startTime)/(1000);
	std::cout<<"\n";

	
	
	VolMagick::writeVolumeFile(outputVol, argv[3]);

	std::cout<<"MSSIM index: "<<sum/voxcount;
      
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


void computeLocalStats(const VolMagick::Volume currVol, uint *currentLoc, float &mean, float &var, const VolMagick::Volume gaussVol, int gaussianSize)
{

	float sum;
	
//	std::cout<<"Enter compute loc..";
//	std::cout<<currVol.XDim();
	
	sum = 0;

	int voxcount = 0;

	int window_size = (gaussianSize-1)/2;


	for(int i = -window_size; i <= window_size; i++)
		for(int j = -window_size; j <= window_size; j++)
			for(int k = -window_size; k <= window_size; k++)
			
			{
			int x = currentLoc[0] + i; int y = currentLoc[1] + j; int z = currentLoc[2] + k;
			
			if(x >= 0 && x < currVol.XDim()
			&& y >= 0 && y < currVol.YDim()
			&& z >= 0 && z < currVol.ZDim())
				{
				voxcount++;
				
				int gx = i + (gaussianSize - 1)/2;
			        int gy = j + (gaussianSize - 1)/2;
			        int gz = k + (gaussianSize - 1)/2;

							
				sum = sum +  gaussVol(gx,gy,gz)*currVol(x,y,z);
				}
					
			}
			
	mean = sum/voxcount;

	float sum_sd = 0;

		for(int i = -window_size; i <= window_size; i++)
	                for(int j = -window_size; j <= window_size; j++)
        	                for(int k = -window_size; k <= window_size; k++)

			{
			int x = currentLoc[0] + i; int y = currentLoc[1] + j; int z = currentLoc[2] + k;
			
			if(x >= 0 && x < currVol.XDim()
			&& y >= 0 && y < currVol.YDim()
			&& z >= 0 && z < currVol.ZDim())
		
			   {
			
				int gx = i + (gaussianSize - 1)/2;
                                int gy = j + (gaussianSize - 1)/2;
                                int gz = k + (gaussianSize - 1)/2;


				sum_sd = sum_sd +  gaussVol(gx,gy,gz)*(currVol(x,y,z) - mean)*(currVol(x,y,z) - mean);
			

			   }		
			}

	if(voxcount==0)
		std::cerr<<"ERROR!, voxel count is zero here:"<<currentLoc[0]<<" "<<currentLoc[1]<<" "<< currentLoc[2] <<" \n";

	var = sum_sd/voxcount;

}



void computeCovar(const VolMagick::Volume inputVol, const VolMagick::Volume inputVol2, uint *currentLoc, float &covar, float mean1, float mean2, const VolMagick::Volume gaussVol, int gaussianSize)

{

float sum_sd = 0;

int voxcount = 0;

 int window_size = (gaussianSize-1)/2;


        for(int i = -window_size; i <= window_size; i++)
                for(int j = -window_size; j <= window_size; j++)
                        for(int k = -window_size; k <= window_size; k++)
			
			{
			int x = currentLoc[0] + i; int y = currentLoc[1] + j; int z = currentLoc[2] + k;
			
			if(x >= 0 && x < inputVol.XDim()
			&& y >= 0 && y < inputVol.YDim()
			&& z >= 0 && z < inputVol.ZDim())
				{
				
				voxcount++;
		
	 			 int gx = i + (gaussianSize - 1)/2;
                                int gy = j + (gaussianSize - 1)/2;
                                int gz = k + (gaussianSize - 1)/2;

		
				sum_sd = sum_sd + gaussVol(gx, gy, gz)*(inputVol(x,y,z) - mean1)*(inputVol2(x,y,z) - mean2);

				}
						
			}

	covar = sum_sd/voxcount;


}

