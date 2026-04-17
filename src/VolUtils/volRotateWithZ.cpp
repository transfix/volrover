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

#include <fstream>

#define Pi 3.1415926

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
  if(argc != 4)
    {
      cerr << "Usage: " << argv[0] << " <input volume file> <Rotate Num> <output volume file>" << endl;
      return 1;
    }

  try
    {
      VolMagickOpStatus status;
      VolMagick::setDefaultMessenger(&status);

      VolMagick::Volume inputVol;
	  VolMagick::VolumeFileInfo inputVol2info;

      VolMagick::readVolumeFile(inputVol, argv[1]);
//	  inputVol2info.read(argv[2]);
//	  fstream fs;
//	  fs.open(argv[3]);

//	  int seeda, seedb, seedc;
	  int Num = atoi(argv[2]);

     VolMagick::Volume outputVol;
	  
	 outputVol.voxelType(inputVol.voxelType());
	 outputVol.dimension(inputVol.dimension());
	 outputVol.boundingBox(inputVol.boundingBox());

//	 float input2Xmin= inputVol2info.XMin();
//	 float input2Ymin = inputVol2info.YMin();
//	 float input2Zmin = inputVol2info.ZMin();



	 float cx = (inputVol.XMin()+inputVol.XMax())*0.5+inputVol.XSpan()/2.0;
	 float cy = (inputVol.YMin() + inputVol.YMax())*0.5+inputVol.YSpan()/2.0;
	 
	 int ni, nj, nk;
	 float x, y;
    
//	  vector <int>  seed[Num];
  
	 float x0, y0, z0;
	 
//	 while(fs>> seeda >> seedb >>seedc)
//	 {
//	    x0 = inputVol2info.XMin() + seeda * inputVol2info.XSpan();
//		y0 = inputVol2info.YMin() + seedb * inputVol2info.YSpan();
//	    z0 = inputVol2info.ZMin() + seedc * inputVol2info.ZSpan();
//		seed[0].push_back(seeda);
//		seed[0].push_back(seedb);
//		seed[0].push_back(seedc);

		 for(int num=0; num<Num; num ++)
		 {	
		 	for(int i=0; i <inputVol.XDim(); i++)
			  for(int j=0; j< inputVol.YDim(); j++)
			    for(int k=0; k<inputVol.ZDim(); k++)
				outputVol(i,j,k,0.0);


		    for(int i=0; i <inputVol.XDim(); i++)
			  for(int j=0; j< inputVol.YDim(); j++)
			    for(int k=0; k<inputVol.ZDim(); k++)
				{
				 if(inputVol(i,j,k)!=0)
				 {
 				   x0 = inputVol.XMin() + i * inputVol.XSpan();
	   		       y0 = inputVol.YMin() + j * inputVol.YSpan();
				   z0 = inputVol.ZMin() + k * inputVol.ZSpan();
				
			     

		//		 cout<<"x0, y0,z0 = "<< x0<<" " << y0 <<" "<<z0<<endl;
				   float  nx = cos(2.0*Pi*(float)num/(float)Num)*(x0-cx) - sin(2.0*Pi*(float)num/(float)Num)*(y0-cy) + cx;
			       float  ny = sin(2.0*Pi*(float)num/(float)Num)*(x0-cx) + cos(2.0*Pi*(float)num/(float)Num)*(y0-cy) + cy;
				   int ni0 = (int)((nx-inputVol.XMin())/inputVol.XSpan()+0.5);
	  			   int nj0 = (int)((ny-inputVol.YMin())/inputVol.YSpan()+0.5);
                   int nk0 = (int)((z0-inputVol.ZMin())/inputVol.ZSpan()+0.5);

				   for(ni = ni0-1; ni <= ni0+1; ni++)
				     for(nj = nj0-1; nj <= nj0+1; nj++)
					   for(nk = nk0-1; nk <= nk0+1; nk++)
					   {
				    	   if(ni>=inputVol.XDim()||nj >= inputVol.YDim()|| nk >= inputVol.ZDim()|| ni<0 || nj <0 || nk <0) continue;
					       else
			    		   {
						  	 float mapx = ni*inputVol.XSpan() + inputVol.XMin();
							 float mapy = nj*inputVol.YSpan() + inputVol.YMin();
							 float mapz = nk*inputVol.ZSpan() + inputVol.ZMin();

	
							 float origx = cos(-2.0*Pi*(float)num/(float)Num)*(mapx-cx) - sin(-2.0*Pi*(float)num/(float)Num)*(mapy-cy) + cx;
							 float origy = sin(-2.0*Pi*(float)num/(float)Num)*(mapx-cx) + cos(-2.0*Pi*(float)num/(float)Num)*(mapy-cy) + cy;
							 float origz = mapz;

				//	cout<<"new i,j,k = "<< ni <<" " << nj <<" "<<nk<<endl;
			//		cout<<"origx,y,z = "<< origx <<" "<< origy <<" " <<z0 <<endl;
							 if(!(origx>inputVol.XMax()|| origx<inputVol.XMin() || origy>inputVol.YMax()|| origy<inputVol.YMin() || origz<inputVol.ZMin() ||
							  origz>inputVol.ZMax())) 
							 outputVol(ni, nj, nk, inputVol.interpolate(origx, origy, origz));
						    } 
					   }
        		  }	
				}
			stringstream ss;
  			ss<<argv[3]<<num<<".rawiv";
			string file=ss.str();
			VolMagick::createVolumeFile(outputVol, file.c_str());
	  }
	 

	// }
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
