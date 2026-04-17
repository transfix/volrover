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





template <class T>
T InnerProduct(T u[3], T v[3])
{
	return (u[0]*v[0]+u[1]*v[1] + u[2]*v[2]);
}

template <class T>
T len(T u[3])
{
	return sqrt(InnerProduct(u,u));
}	

template <class T>
void normalize(T* u)
{	
	T leng = len(u);
	for(int i=0; i<3; i++)
	u[i] = u[i]/leng; 
}

template <class T>
void defineRotationMatrix(T matrix[3][3], T alpha, T* v)
{
	T c, s, t;
	c = cos(alpha);
	s = sin(alpha);
	t = 1-c;
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)
		matrix[i][j] = 0.0;
	matrix[0][0] = v[0]*v[0]*t + c;
	matrix[0][1] = v[0]*v[1]*t - v[2]*s;
	//matrix[0][1] = v[0]*v[1]*t + v[2]*s;
	matrix[0][2] = v[0]*v[2]*t + v[1]*s;
	//matrix[0][2] = v[0]*v[2]*t - v[1]*s;
	matrix[1][0] = v[0]*v[1]*t + v[2]*s;
	//matrix[1][0] = v[0]*v[1]*t - v[2]*s;
	matrix[1][1] = v[1]*v[1]*t + c;
	matrix[1][2] = v[1]*v[2]*t - v[0]*s;
	//matrix[1][2] = v[1]*v[2]*t + v[0]*s;
	matrix[2][0] = v[0]*v[2]*t - v[1]*s;
	//matrix[2][0] = v[0]*v[2]*t + v[1]*s;
	matrix[2][1] = v[1]*v[2]*t + v[0]*s;
	//matrix[2][1] = v[1]*v[2]*t - v[0]*s;
	matrix[2][2] = v[2]*v[2]*t + c;
}

template <class T>
void MatrixMultipleVector(T A[3][3], T* v, T* u)
{
	for(int i=0; i<3; i++)
	{
		u[i] = A[i][0]*v[0] + A[i][1]*v[1]+A[i][2]*v[2];
	}
}



int main(int argc, char **argv)
{
  if(argc != 6)
    {
      cerr << "Usage: " << argv[0] << " <input Original volume file>  <input seed file ><input real axis> <Rotate Num> <output seed file>" << endl;
      return 1;
    }

  try
    {
      VolMagickOpStatus status;
      VolMagick::setDefaultMessenger(&status);

      VolMagick::Volume inputVol;

      VolMagick::readVolumeFile(inputVol, argv[1]);

	  ifstream fs;
	  fs.open(argv[2], ios::in);
	  readUntilNewline(fs);
	  readUntilNewline(fs);
	  readUntilNewline(fs);


	  ifstream faxis;
	  faxis.open(argv[3], ios::in);

	  float a, b, c;
	  float v[3], center[3];

	  //The first line is the axis
	  faxis>>v[0]>>v[1]>>v[2];

	 faxis.seekg (0, ios::beg);

	  for(int i=0; i<3; i++)
	  	center[i]=0.0;
	  int axisNum=0;
	  while(faxis>> a>> b>>c){
	  	center[0] += a;
		center[1] += b;
		center[2] += b;
	  	axisNum++;
	  }
	  faxis.close();
	
	  for(int i=0; i<3; i++)
	  {
	  	center[i]/=(float)axisNum;
		v[i]-= center[i];
	  }


		
	  string s1, s2;


	  normalize(v);

	  float Matrix[3][3];


	  int seeda, seedb, seedc;
	  int Num = atoi(argv[4]);


	 int ni, nj, nk;
    
	  vector <int>  seed[Num];
  
	 float x0, y0, z0;
	 
	while(fs>>s1>>seedb>>seedc>>s2){
  		seeda = atoi(s1.erase(0,7).c_str());
	//	cout<<seeda<< " "<< seedb<<" " <<seedc<<" " <<endl;
	    x0 = inputVol.XMin() + seeda * inputVol.XSpan();
		y0 = inputVol.YMin() + seedb * inputVol.YSpan();
	    z0 = inputVol.ZMin() + seedc * inputVol.ZSpan();

		float inputseed[3];
		inputseed[0] = x0-center[0];
		inputseed[1] = y0-center[1];
		inputseed[2] = z0-center[2];

		 for(int num=0; num<Num; num ++)
		 {	
		    float alpha = 2.0*Pi*num/Num;
			defineRotationMatrix(Matrix, alpha, v);
			float u[3];
			MatrixMultipleVector(Matrix, inputseed, u);

			ni = (int)((u[0]+center[0]-inputVol.XMin())/inputVol.XSpan()+0.5);
			nj = (int)((u[1]+center[1]-inputVol.YMin())/inputVol.YSpan()+0.5);
            nk = (int)((u[2]+center[2]-inputVol.ZMin())/inputVol.ZSpan()+0.5);
			cout<<ni <<" " <<nj<<" "<<nk <<endl;
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
