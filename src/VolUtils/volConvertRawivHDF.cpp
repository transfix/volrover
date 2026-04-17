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


#include <H5IM.h>
//#include <pal_rgb.h>
#include <stdlib.h>

#define DATASETNAME "Converted Volume from .rawiv" 
#define RANK   3

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
  if(argc < 3)
    {
      cerr << 
	"Usage: inputvile, outputvile"<< endl;

      return 1;
    }

  VolMagick::Volume inputVol;


  VolMagick::readVolumeFile(inputVol,argv[1]); ///first argument is input volume

      
    hid_t       file, dataset;         /* file and dataset handles */
    hid_t       datatype, dataspace;   /* handles */
    hsize_t     dimsf[3];              /* dataset dimensions */
    herr_t      status;           


    switch(inputVol.voxelType()) 
      {
      case 0:
	{
	  datatype = H5Tcopy(H5T_NATIVE_CHAR);
	  std::cout<<"allocating data\n";                  
	  //unsigned char         data[inputVol.XDim()][inputVol.YDim()][inputVol.ZDim()];          /* data to write */
	  unsigned char*   data = new unsigned char[inputVol.XDim()*inputVol.YDim()*inputVol.ZDim()];          /* data to write */

	  std::cout<<"done allocating\n";
          for(unsigned int i=0; i<inputVol.XDim(); i++)
	    for(unsigned int j=0; j<inputVol.YDim(); j++)
	      for(unsigned int k=0; k<inputVol.ZDim(); k++)

		{
		  data[i*inputVol.XDim()*inputVol.YDim() + j*inputVol.YDim() + k] = inputVol(i,j,k);
		}


	  std::cout<<"file created here \n";
    
	  file = H5Fcreate(argv[2], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	  dimsf[0] = inputVol.XDim();
	  dimsf[1] = inputVol.YDim();
	  dimsf[2] = inputVol.ZDim();

	  dataspace = H5Screate_simple(RANK, dimsf, NULL); 
	
 
	  std::cout<<"Voxel Type: "<<inputVol.voxelType()<<std::endl;
    
 

	  status = H5Tset_order(datatype, H5T_ORDER_LE);


	  std::cout<<"File being created\n"; 
	  dataset = H5Dcreate(file, DATASETNAME, datatype, dataspace, H5P_DEFAULT, H5P_DATASET_CREATE_DEFAULT, H5P_DATASET_ACCESS_DEFAULT);
	  std::cout<<"done file creation\n";
	  status = H5Dwrite(dataset, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	  //       status = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    
	  std::cout<<"data written out \n";

	  /*
	   *      * Close/release resources.
	   *           */
	  H5Sclose(dataspace);
	  H5Tclose(datatype);
	  H5Dclose(dataset);
	  H5Fclose(file);
	}
	   
	break;
       
      case 1:
	{
	  datatype = H5Tcopy(H5T_NATIVE_SHORT);
	
	  std::cout<<"allocating data\n";
        
	  short*   data = new short[inputVol.XDim()*inputVol.YDim()*inputVol.ZDim()];          /* data to write */

	  std::cout<<"done allocating\n";
          for(unsigned int i=0; i<inputVol.XDim(); i++)
	    for(unsigned int j=0; j<inputVol.YDim(); j++)
	      for(unsigned int k=0; k<inputVol.ZDim(); k++)

		{
		  data[i*inputVol.XDim()*inputVol.YDim() + j*inputVol.YDim() + k] = inputVol(i,j,k);
		}

	   
	  std::cout<<"file created here \n";
    
	  file = H5Fcreate(argv[2], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	  dimsf[0] = inputVol.XDim();
	  dimsf[1] = inputVol.YDim();
	  dimsf[2] = inputVol.ZDim();

	  dataspace = H5Screate_simple(RANK, dimsf, NULL); 
	
 
	  std::cout<<"Voxel Type: "<<inputVol.voxelType()<<std::endl;
    
 

	  status = H5Tset_order(datatype, H5T_ORDER_LE);


	  std::cout<<"File being created\n"; 
	  dataset = H5Dcreate(file, DATASETNAME, datatype, dataspace, H5P_DEFAULT, H5P_DATASET_CREATE_DEFAULT, H5P_DATASET_ACCESS_DEFAULT);
	  std::cout<<"done file creation\n";
	  status = H5Dwrite(dataset, H5T_NATIVE_SHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	  //       status = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    
	  std::cout<<"data written out \n";

	  /*
	   *      * Close/release resources.
	   *           */
	  H5Sclose(dataspace);
	  H5Tclose(datatype);
	  H5Dclose(dataset);
	  H5Fclose(file);
	}
	break;
       
      case 2:
	{
	  datatype = H5Tcopy(H5T_NATIVE_INT);;

	  std::cout<<"allocating data\n";

	  int*   data = new int[inputVol.XDim()*inputVol.YDim()*inputVol.ZDim()];          /* data to write */

	  std::cout<<"done allocating\n";
          for(unsigned int i=0; i<inputVol.XDim(); i++)
	    for(unsigned int j=0; j<inputVol.YDim(); j++)
	      for(unsigned int k=0; k<inputVol.ZDim(); k++)

		{
		  data[i*inputVol.XDim()*inputVol.YDim() + j*inputVol.YDim() + k] = inputVol(i,j,k);
		}

	
	   
	  std::cout<<"file created here \n";
    
	  file = H5Fcreate(argv[2], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	  dimsf[0] = inputVol.XDim();
	  dimsf[1] = inputVol.YDim();
	  dimsf[2] = inputVol.ZDim();

	  dataspace = H5Screate_simple(RANK, dimsf, NULL); 
	
 
	  std::cout<<"Voxel Type: "<<inputVol.voxelType()<<std::endl;
    
 

	  status = H5Tset_order(datatype, H5T_ORDER_LE);


	  std::cout<<"File being created\n"; 
	  dataset = H5Dcreate(file, DATASETNAME, datatype, dataspace, H5P_DEFAULT, H5P_DATASET_CREATE_DEFAULT, H5P_DATASET_ACCESS_DEFAULT);
	  std::cout<<"done file creation\n";
	  status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	  //       status = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    
	  std::cout<<"data written out \n";

	  /*
	   *      * Close/release resources.
	   *           */
	  H5Sclose(dataspace);
	  H5Tclose(datatype);
	  H5Dclose(dataset);
	  H5Fclose(file);
	}

	break;

      case 3:
	
	{
	  datatype = H5Tcopy(H5T_NATIVE_FLOAT);
	  std::cout<<"allocating data\n";
	  float*   data = new float[inputVol.XDim()*inputVol.YDim()*inputVol.ZDim()];          /* data to write */

	  std::cout<<"done allocating\n";
	  for(unsigned int i=0; i<inputVol.XDim(); i++)
	    for(unsigned int j=0; j<inputVol.YDim(); j++)
	      for(unsigned int k=0; k<inputVol.ZDim(); k++)
	  
		{
		  data[i*inputVol.XDim()*inputVol.YDim() + j*inputVol.YDim() + k] = inputVol(i,j,k);	      
		}
	   
	  std::cout<<"file created here \n";
    
	  file = H5Fcreate(argv[2], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	  dimsf[0] = inputVol.XDim();
	  dimsf[1] = inputVol.YDim();
	  dimsf[2] = inputVol.ZDim();

	  dataspace = H5Screate_simple(RANK, dimsf, NULL); 
	
 
	  std::cout<<"Voxel Type: "<<inputVol.voxelType()<<std::endl;
    
 

	  status = H5Tset_order(datatype, H5T_ORDER_LE);


	  std::cout<<"File being created\n"; 
	  dataset = H5Dcreate(file, DATASETNAME, datatype, dataspace, H5P_DEFAULT, H5P_DATASET_CREATE_DEFAULT, H5P_DATASET_ACCESS_DEFAULT);
	  std::cout<<"done file creation\n";
	  status = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	  //       status = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    
	  std::cout<<"data written out \n";

	  /*
	   *      * Close/release resources.
	   *           */
	  H5Sclose(dataspace);
	  H5Tclose(datatype);
	  H5Dclose(dataset);
	  H5Fclose(file);

	}
	break;

      case 4:
 	{
	  datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
	  std::cout<<"allocating data\n";
	  double*   data = new double[inputVol.XDim()*inputVol.YDim()*inputVol.ZDim()];          /* data to write */

	  std::cout<<"done allocating\n";
          for(unsigned int i=0; i<inputVol.XDim(); i++)
	    for(unsigned int j=0; j<inputVol.YDim(); j++)
	      for(unsigned int k=0; k<inputVol.ZDim(); k++)

		{
		  data[i*inputVol.XDim()*inputVol.YDim() + j*inputVol.YDim() + k] = inputVol(i,j,k);
		}

	   
	  std::cout<<"file created here \n";
    
	  file = H5Fcreate(argv[2], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	  dimsf[0] = inputVol.XDim();
	  dimsf[1] = inputVol.YDim();
	  dimsf[2] = inputVol.ZDim();

	  dataspace = H5Screate_simple(RANK, dimsf, NULL); 
	
 
	  std::cout<<"Voxel Type: "<<inputVol.voxelType()<<std::endl;
    
 

	  status = H5Tset_order(datatype, H5T_ORDER_LE);


	  std::cout<<"File being created\n"; 
	  dataset = H5Dcreate(file, DATASETNAME, datatype, dataspace, H5P_DEFAULT, H5P_DATASET_CREATE_DEFAULT, H5P_DATASET_ACCESS_DEFAULT);
	  std::cout<<"done file creation\n";
	  status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	  //       status = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    
	  std::cout<<"data written out \n";

	  /*
	   *      * Close/release resources.
	   *           */
	  H5Sclose(dataspace);
	  H5Tclose(datatype);
	  H5Dclose(dataset);
	  H5Fclose(file);

	}
	break;


      case 5:
	std::cout<<"\n Unsupported conversion"<<endl;
	return 0;
    
      default:
	break;
      }

    /*

    for(unsigned int i=0; i<inputVol.XDim(); i++)
    for(unsigned int j=0; j<inputVol.YDim(); j++)
    for(unsigned int k=0; k<inputVol.ZDim(); k++)
	  
    {
    data[i][j][k] = inputVol(i,j,k);	      
    }
    
 
    std::cout<<"file created here \n";
    
    file = H5Fcreate(argv[2], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    dimsf[0] = inputVol.XDim();
    dimsf[1] = inputVol.YDim();
    dimsf[2] = inputVol.ZDim();

    dataspace = H5Screate_simple(RANK, dimsf, NULL); 
	
 
    std::cout<<"Voxel Type: "<<inputVol.voxelType()<<std::endl;
    
 

    status = H5Tset_order(datatype, H5T_ORDER_LE);


    std::cout<<"File being created\n"; 
    dataset = H5Dcreate(file, DATASETNAME, datatype, dataspace, H5P_DEFAULT, H5P_DATASET_CREATE_DEFAULT, H5P_DATASET_ACCESS_DEFAULT);
    std::cout<<"done file creation\n";

    status = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    //       status = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    
    std::cout<<"data written out \n";


    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Dclose(dataset);
    H5Fclose(file);
    */
 
    return 0;
  

}
