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

void readUntilNewline(ifstream &fin){
	char c='a';
	while(c!='\n'){
	fin.get(c);
    }
}


int main(int argc, char **argv)
{
  if(argc != 6)
    {
      cerr << "Usage: " << argv[0] << " <input Original volume file> <input cutvolume> <input seed file ><Rotate Num> <output seed file>" << endl;
      return 1;
    }

  try
    {
      VolMagickOpStatus status;
      VolMagick::setDefaultMessenger(&status);

      VolMagick::Volume inputVol;
	  VolMagick::VolumeFileInfo inputVol2info;

      VolMagick::readVolumeFile(inputVol, argv[1]);
	  inputVol2info.read(argv[2]);
	//  fstream fs;
//	  fs.open(argv[3]);

	  ifstream fs;
	  fs.open(argv[3], ios::in);
	  readUntilNewline(fs);
	  readUntilNewline(fs);
	  readUntilNewline(fs);

	  string s1, s2;

	  int seeda, seedb, seedc;
	  int Num = atoi(argv[4]);


	 float cx = (inputVol.XMin()+inputVol.XMax())*0.5;
	 float cy = (inputVol.YMin() + inputVol.YMax())*0.5;
	 
	 int ni, nj, nk;
	 float x, y;
    
	  std::vector<std::vector<int>> seed(Num);
  
	 float x0, y0, z0;
	 
	while(fs>>s1>>seedb>>seedc>>s2){
  		seeda = atoi(s1.erase(0,7).c_str());
	//	cout<<seeda<< " "<< seedb<<" " <<seedc<<" " <<endl;
	    x0 = inputVol2info.XMin() + seeda * inputVol2info.XSpan();
		y0 = inputVol2info.YMin() + seedb * inputVol2info.YSpan();
	    z0 = inputVol2info.ZMin() + seedc * inputVol2info.ZSpan();


		 for(int num=0; num<Num; num ++)
		 {	
	 	 	float  nx = cos(2.0*Pi*(float)num/(float)Num)*(x0-cx) - sin(2.0*Pi*(float)num/(float)Num)*(y0-cy) + cx;
			float  ny = sin(2.0*Pi*(float)num/(float)Num)*(x0-cx) + cos(2.0*Pi*(float)num/(float)Num)*(y0-cy) + cy;
			ni = (int)((nx-inputVol.XMin())/inputVol.XSpan()+0.5);
			nj = (int)((ny-inputVol.YMin())/inputVol.YSpan()+0.5);
            nk = (int)((z0-inputVol.ZMin())/inputVol.ZSpan()+0.5);
			if(ni>=inputVol.XDim()||nj >= inputVol.YDim()|| nk >= inputVol.ZDim()|| ni<0 || nj <0 || nk <0) continue;
			else
			{
			   seed[num].push_back(ni);
			   seed[num].push_back(nj);
			   seed[num].push_back(nk);
			}
		 }
	  }
	  fs.close();

	  fstream fs1;
	  fs1.open(argv[5], ios::out);
	  fs1<<"<!DOCTYPE pointclassdoc>"<<endl;
	  fs1<<"<pointclassdoc>" <<endl;

	  string str[13]={"ffff00", "ff0000","00ff00","ff00ff","0000ff","ff5500","336699", "00ffff", "c0c0c0","800000", "800080", "808000", "008080"};
	  
	  if(Num>13) cout<<"You need add more colors." <<endl;
	  for(int num=0; num< Num; num++)
	  {
	  	 fs1<<" <pointclass timestep=\"0\" name=\"Class "<<num<<"\" color="<<"\"#"<<str[num]<<"\" variable=\"0\" >"<<endl;
		 for(int j = 0; j < seed[num].size()/3; j++)
		  fs1<<"  <point>" <<seed[num][3*j] <<" "<< seed[num][3*j+1] <<" "<<seed[num][3*j+2]<<"</point>"<<endl;
		 fs1<<" </pointclass>"<<endl;
	  }
	  fs1<<"</pointclassdoc>"<<endl;

	  fs1.close();

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
