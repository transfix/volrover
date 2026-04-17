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
#include <limits>

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

#if ! defined (__APPLE__)
#include <omp.h>
#endif
#include <time.h>

//#include <complex.h> 
#include <fftw3.h>

//VolMagick::Volume Filter(VolMagick::Volume inputVol);
VolMagick::Volume Filterfft(VolMagick::Volume inputVol);

using namespace std;


int main(int argc, char **argv)
{
  if(argc<5)
    {
      cerr <<
 
	"Usage: inputfile (projections), projectionAngles, number of slices, outputfile (reconstruction) filter(1=yes, 0=no)" << endl;


      return 1;
    }


  try
    {
           
      VolMagick::Volume inputVol;

      VolMagick::Volume outputVol;   

      VolMagick::Volume startPtVol;

     VolMagick::Volume endPtVol;

      VolMagick::Volume rotatedVol;

      VolMagick::readVolumeFile(inputVol,argv[1]); ///first argument is input volume
    
      
	 
      static float PI = 3.14159265;
 
      std::cout<<"vol read\n";



#if ! defined (__APPLE__)
       int numProcs = omp_get_num_procs();
#else 
       int numProcs = 0;
#endif
      std::cout<<"Num of processors: "<< numProcs<<"\n";

      uint numAngles = 0;
	
      FILE* fp = fopen(argv[2],"r");
     
      std::cout<<"Counting angles \n";
	
      float tmp;	
	

	while(fscanf(fp, "%f", &tmp) != EOF)
		numAngles++;

	fclose(fp);


   // initialize projection output volume

//	 outputVol.dimension(VolMagick::Dimension(inputVol.XDim(), inputVol.YDim(), numAngles));

	 VolMagick::Dimension dim;

	 dim.xdim = inputVol.XDim();

	 dim.ydim = inputVol.YDim();

	 dim.zdim = atoi(argv[3]);
	
	 std::cout<<"Number of slices in the reconstructed volume: "<< dim.zdim <<"\n";
	
	 VolMagick::VoxelType vox = VolMagick::Float;

	 outputVol.voxelType(vox);

         VolMagick::BoundingBox outputBB;

//change this to let user provide boundingbox 

	 outputBB.minx = 0;

	 outputBB.maxx = dim.xdim - 1;
	 
	 outputBB.miny = 0;	

	 outputBB.maxy = dim.ydim - 1;

	 outputBB.minz = 0;

	 outputBB.maxz = dim.zdim - 1;

	 outputVol.boundingBox(outputBB);

	 outputVol.dimension(dim);
	
		
	// check
	std::cout<<outputVol.XDim()<<" "<<outputVol.YDim()<<" "<<outputVol.ZDim()<<"\n";
	
	std::cout<<outputVol.boundingBox().maxx<<" "<<outputVol.boundingBox().maxy<<" "<<outputVol.boundingBox().maxz<<"\n";

	std::cout<<outputVol.boundingBox().minx<<" "<<outputVol.boundingBox().miny<<" "<<outputVol.boundingBox().minz<<"\n";

	std::cout<<outputVol.voxelType()<<"\n";
	
//	
//
//	VolMagick::createVolumeFile(argv[3], outputBB, dim, inputVol.voxelType());

        
	std::cout<<"Total number of angles: "<<numAngles<<"\n";

    float *angles = new float[numAngles];
	
	FILE* fp1 = fopen(argv[2], "r");
	
	for(uint i=0; i<numAngles; i++) {

	float tempangle;
	
	fscanf(fp1, "%f", &tempangle);

	tempangle = tempangle*PI/180;

	angles[i] = tempangle;

	}

	fclose(fp1);

//        int startTime;
//        startTime = clock();
	time_t t0, t1;
	
	t0 = time(NULL);	

	int filtering = atoi(argv[5]);

	std::cout<<"Filtering? (yes =1, no =0) "<< filtering<<"\n";

	if(filtering==1)
	{

	std::cout<<"Begin filtering\n";

	inputVol = Filterfft(inputVol);

       VolMagick::createVolumeFile("./filtered.rawiv", inputVol);


       VolMagick::writeVolumeFile(inputVol, "./filtered.rawiv");

	std::cout<<"End filtering\n";

	}

//return 0;

	float delta = 0.1;

	#pragma omp parallel for schedule(static, numAngles/numProcs)           
	for(uint x=0; x<outputVol.XDim(); x++)
	{
	 
	std::cerr<<x<<"..";
	
	    for(uint y=0; y<outputVol.YDim(); y++)
		for(uint z=0; z<outputVol.ZDim(); z++)
              {
		
		float intensity_update = 0;

		uint countintensity = 1;
			

		for(uint angle = 0; angle<numAngles; angle++)
		{
//		uint angle = 23;
		   float tempangle = angles[angle];
//			float tempangle = 0.5;				
		

		float x_shift = (float)x - (float)outputVol.XDim()/2;

		float z_shift = (float)z - (float)outputVol.ZDim()/2;

	  
		float xcoordinate = x_shift*cos(tempangle) + z_shift*sin(tempangle) + (float)outputVol.XDim()/2;
	// Linear interpolation

		if(floor(xcoordinate) > 0 && ceil(xcoordinate) < inputVol.XDim())
			  {		
			
		 float intensity;

		//std::cerr<<xcoordinate<<"..";

			if(floor(xcoordinate) != ceil(xcoordinate))
			 intensity = inputVol(floor(xcoordinate), y, angle) + ( inputVol(ceil(xcoordinate), y, angle) - inputVol(floor(xcoordinate), y, angle) )*(xcoordinate - floor(xcoordinate))/( ceil(xcoordinate) - floor(xcoordinate)); 

			else
				intensity = inputVol(xcoordinate, y, angle);

			 if(intensity)
				countintensity++;

		        intensity_update = intensity_update + intensity;


  
			   }	
	
        		}

//			outputVol(x,y,z, outputVol(x,y,z)/numAngles);
			outputVol(x,y,z, outputVol(x,y,z) + intensity_update/countintensity);			
		
		}
	}
    
	
		
	
//        int endTime;

 //       endTime = clock();
		
	t1 = time(NULL);
	
	std::cout<<"Time elapsed: "<<(long)(t1-t0);

//        std::cout<<"Time elapsed: "<<(endTime-startTime)/(1000);
        std::cout<<"\n";

	

      	VolMagick::createVolumeFile(argv[4], outputVol);

	std::cout<<"Created outputfile\n";


         
     VolMagick::writeVolumeFile(outputVol, argv[4]);

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


VolMagick::Volume Filterfft(VolMagick::Volume inputVol)
{

uint n = inputVol.XDim()*inputVol.YDim();

//fftw_complex *spatial = fftw_malloc(sizeof((fftw_complex)*n));

//fftw_complex *freq    = fftw_malloc(sizeof((fftw_complex)*n));


        VolMagick::Volume outputVol;

         VolMagick::VoxelType vox = VolMagick::Float;


         outputVol.boundingBox(inputVol.boundingBox());

         outputVol.dimension(inputVol.dimension());

         outputVol.voxelType(vox);

#if defined (__APPLE__)
	float min = 1000000000;
	float max = -100000000;
#else
	float min = 100000000000;
	float max = -10000000000;
#endif

for(int k=0;k<inputVol.ZDim();k++)
	{		


	fftw_complex *spatial, *freq;

	spatial = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * n );

	freq = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * n);


	int count = 0;

	for(int i=0;i<inputVol.XDim(); i++)
		for(int j=0; j<inputVol.YDim();j++)
		{
		spatial[count][0] = inputVol(i,j,k);
		spatial[count][1] = 0; // imaginary i
		count++;

		if(inputVol(i,j,k) < min)
			min = inputVol(i,j,k);
		
		if(inputVol(i,j,k) > max)
			max = inputVol(i,j,k);
			

		}

	 fftw_plan forward = fftw_plan_dft_2d(inputVol.XDim(), inputVol.YDim(), spatial, freq, FFTW_FORWARD, FFTW_ESTIMATE);
         
	fftw_execute(forward);
		
//Filtering

	float max =sqrt( (inputVol.XDim()/2)*(inputVol.XDim()/2) + (inputVol.YDim()/2)*(inputVol.YDim()/2) );

        for(int i=0;i<inputVol.XDim(); i++)
                for(int j=0; j<inputVol.YDim();j++)
                {
		
		float slope = (1 - 1/max)/sqrt( inputVol.XDim()/2 * inputVol.XDim()/2 + inputVol.YDim()/2 * inputVol.YDim()/2);
		
		float x = sqrt( (inputVol.XDim()/2 - i) * (inputVol.XDim()/2 - i) + 
				(inputVol.YDim()/2 - j) * (inputVol.YDim()/2 - j) ) ;

		
                freq[j+inputVol.YDim()*i][0] = freq[j+inputVol.YDim()*i][0] *4*( 1 - slope*x);

	//now the imaginary part
	//
	        freq[j+inputVol.YDim()*i][1] = freq[j+inputVol.YDim()*i][1] *4*( 1 - slope*x);


//	if(x>15)
//	outputVol(i,j,k, 2* (1 - slope*x));
	
// The fourier space is multiplied by a ramp function. In 2D, the ramp function peaks to 1 at the center and falls off to zero at the corenrs
// The minimum value of the ramp function is 1/max while the max value is 1. Every point in between is interpolated
		}


         fftw_plan inverse = fftw_plan_dft_2d(inputVol.XDim(), inputVol.YDim(), freq, spatial, FFTW_BACKWARD, FFTW_ESTIMATE);


	fftw_execute(inverse);

		
        for(int i=0;i<inputVol.XDim(); i++)
                for(int j=0; j<inputVol.YDim();j++)
		
			outputVol(i,j,k, spatial[j + inputVol.YDim()*i][0] )  ;
//			outputVol(i,j,k, (freq[j + inputVol.YDim()*i][0] ) )  ;


	free(spatial);
	free(freq);

	}

	std::cout<<"min, max: "<<min<<" "<<max<<"\n";

	float newmin = 1000000000;
	float newmax = -100000000;

	        for(int i=0;i<inputVol.XDim(); i++)
                for(int j=0; j<inputVol.YDim();j++)
		for(int k = 0; k<inputVol.ZDim(); k++) 
		{
			if(outputVol(i,j,k) < newmin)
				newmin = outputVol(i,j,k);

			if(outputVol(i,j,k) > newmax)
				newmax = outputVol(i,j,k);
		}



	std::cout<<"newmin, newmax: "<<newmin<<" "<<newmax<<"\n";

                for(int i=0;i<inputVol.XDim(); i++)
                for(int j=0; j<inputVol.YDim();j++)
                for(int k = 0; k<inputVol.ZDim(); k++)

			outputVol(i,j,k, min + (outputVol(i,j,k) - newmin)*max/(newmax - newmin) );

//			outputVol(i,j,k, (outputVol(i,j,k) - newmin)/(newmax-newmin) );


	return outputVol;

}
