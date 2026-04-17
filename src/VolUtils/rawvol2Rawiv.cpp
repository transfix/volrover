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
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>



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
  if(argc < 3)
    {
      std:: cerr << 
	"Usage: " << argv[0] << 

	
	" <int text file>  <output volume> \n";

      return 1;
    }

  try
    {

	 
    ifstream input(argv[1]);
    string line;
  
    string filename;
	string type;
	float span[3];
	

	unsigned int dim[3];

     

	  while(input.good())
	  {
	 	getline (input, line);
		if(strncmp(line.c_str(),"data = ", 7) == 0)
			filename= line.substr(7);
		if(strncmp(line.c_str(), "type = ",7) ==0)
		    type = line.substr(7);
		if(strncmp(line.c_str(), "dimension = ", 12) == 0)
		{
			string s(line.substr(12));
			istringstream iss(s);
			iss >> dim[0] >> dim[1] >> dim[2];
		}
		 if(strncmp(line.c_str(), "ratio = ", 8) ==0)
		 {
		 	string s(line.substr(8));
			istringstream iss(s);
			iss >> span[0] >> span[1] >> span[2];
		 }
	 }

        
	  input.close();

	  
	  float xmax, xmin, ymax, ymin, zmax, zmin;
	  xmin = 0;
	  xmax = (float)dim[0]*span[0];
	  ymin = 0;
	  ymax = (float)dim[1]*span[1];
	  zmin = 0;
	  zmax = (float)dim[2]*span[2];


	  float origin[3];

	  origin[0] = xmin;
	  origin[1] = ymin;
	  origin[2] = zmin;


      unsigned int numverts = dim[0]*dim[1]*dim[2];
	  unsigned int numcells = (dim[0]-1)*(dim[1]-1)*(dim[2]-1);


	  unsigned int d0, d1, d2;
	  d0 = dim[0];
	  d1 = dim[1];
	  d2 = dim[2];


      if(!big_endian())
	  {
	  	SWAP_32(&xmin);
		SWAP_32(&ymin);
		SWAP_32(&zmin);

		SWAP_32(&xmax);
		SWAP_32(&ymax);
		SWAP_32(&zmax);

		SWAP_32(&numverts);
		SWAP_32(&numcells);
		
		for(int i=0; i<3; i++) 
		{
			SWAP_32(&(dim[i]));
			SWAP_32(&(origin[i]));
			SWAP_32(&(span[i]));
		}
	  } 

	   FILE* fout = fopen(argv[2], "wb");
	  fwrite(&xmin, sizeof(float), 1, fout); 
	  fwrite(&ymin, sizeof(float), 1, fout); 
	  fwrite(&zmin, sizeof(float), 1, fout); 
	  fwrite(&xmax, sizeof(float), 1, fout); 
	  fwrite(&ymax, sizeof(float), 1, fout); 
	  fwrite(&zmax, sizeof(float), 1, fout); 

	  fwrite(&numverts, sizeof(unsigned int), 1, fout);
	  fwrite(&numcells, sizeof(unsigned int), 1, fout);

	  for(int i=0; i<3; i++)
	  {	
	  	fwrite(&dim[i], sizeof(unsigned int), 1, fout);
	  }


	  for(int i=0; i<3; i++)
	  {	
	  	fwrite(&origin[i], sizeof(float), 1, fout);
	  }

	  for(int i=0; i<3; i++)
	  {	
	  	fwrite(&span[i], sizeof(float), 1, fout);
	  }




       FILE* fp=fopen(filename.c_str(), "rb");
	  if(!fp) 
	  {
	  	 cout<<"File " << filename <<" has problem." << endl;
		 return 1;
	  }
	
	
      if(type=="FLOAT"){

  	 float tmp, tmp1;

     for( int ix = 0; ix<d0*d1*d2; ix++)
		    {
				fread(&tmp, sizeof(float), 1, fp);
		 //	std::cout<<tmp<< std::endl;
			tmp1 = tmp;
		if(!big_endian())		SWAP_32(&tmp1);
		//	if(ix+jy+kz<100)	std::cout<<tmp<<" "<< std::endl;
				fwrite(&tmp1, sizeof(float), 1, fout);
			}
    }
   else if(type == "UCHAR") 
   {
 	 char tmp, tmp1;

	 std::cout<<d0 <<" " << d1 <<" "<<d2 << std::endl;
     for( int ix = 0; ix<d0*d1*d2; ix++)
		    {
				fread(&tmp, sizeof(char), 1, fp);
		 //	std::cout<<tmp<< std::endl;
			tmp1 = tmp;
		if(!big_endian())		SWAP_16(&tmp1);
		//	if(ix+jy+kz<100)	std::cout<<tmp<<" "<< std::endl;
				fwrite(&tmp1, sizeof(char), 1, fout);
			}
    }
   else if(type == "USHORT") 
   {
 	 short tmp, tmp1;

	 std::cout<<d0 <<" " << d1 <<" "<<d2 << std::endl;
     for( int ix = 0; ix<d0*d1*d2; ix++)
		    {
				fread(&tmp, sizeof(short), 1, fp);
		 //	std::cout<<tmp<< std::endl;
			tmp1 = tmp;
		if(!big_endian())		SWAP_16(&tmp1);
		//	if(ix+jy+kz<100)	std::cout<<tmp<<" "<< std::endl;
				fwrite(&tmp1, sizeof(short), 1, fout);
			}
    }

    else throw VolMagick::UnsupportedVolumeFileType("Unknown volume type");


	fclose(fp);
	fclose(fout);

     cout<< "done! "<<endl;
      /*
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
	
	  string ss=string(argv[2]);
	 
	  FILE* fp = fopen(ss.c_str(), "wb");
	  if(!fp){
	  	cout<<"Cannot write. "<<argv[2] << endl;
		return 0;
	  }
	
	  int xd = inputVol.XDim();
	  int yd = inputVol.YDim();
	  int zd = inputVol.ZDim();
	  float xmin = inputVol.XMin();
	  float ymin = inputVol.YMin();
	  float zmin = inputVol.ZMin();

	  float Minimum = inputVol.min();

	  fwrite(&xmin, sizeof(float), 1, fp);
	  fwrite(&ymin, sizeof(float), 1, fp);
	  fwrite(&zmin, sizeof(float), 1, fp);

		
	  fwrite(&xd, sizeof(int), 1, fp);
	  fwrite(&yd, sizeof(int), 1, fp);
	  fwrite(&zd, sizeof(int), 1, fp);

	  fwrite(&span[0], sizeof(float), 1, fp);
	  fwrite(&span[1], sizeof(float), 1, fp);
	  fwrite(&span[2], sizeof(float), 1, fp);


	float tmp;

	  for( int kz = 0; kz<inputVol.ZDim(); kz++)
	   for( int jy = 0; jy<inputVol.YDim(); jy++)
	     for( int ix = 0; ix<inputVol.XDim(); ix++)
		    {
				tmp = inputVol(kz,jy,ix)-Minimum;
		//	if(ix+jy+kz<100)	std::cout<<tmp<<" "<< std::endl;
				fwrite(&tmp, sizeof(float), 1, fp);
			}
    

	fclose(fp);

    std::cout<<"done!"<<std::endl;
 
   */

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
