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



using namespace std;

int main(int argc, char **argv)
{
  if(argc < 3)
    {
      std:: cerr << 
	"Usage: " << argv[0] << 

	
	"  <first volume>  <output volume> [ bool 0/1 ].   \n"; 
	cerr<<"bool 1 keeps the original volume, bool 0 only the reflection, default 0."  << endl;

      return 1;
    }

  try
    {
      VolMagick::Volume inputVol;

      
      VolMagick::Volume outputVol;

      VolMagick::readVolumeFile(inputVol,argv[1]); ///first argument is input volume
      
      VolMagick::VolumeFileInfo volinfo1;
      volinfo1.read(argv[1]);
      std::cout << volinfo1.filename() << ":" <<std::endl;

       
      std::cout<<"minVol1 , maxVol1: "<<volinfo1.min()<<" "<<volinfo1.max()<<std::endl;;

  
      VolMagick::BoundingBox bbox;
	  bbox.minx = inputVol.XMin();
	  bbox.maxx = inputVol.XMax();
	  bbox.miny = inputVol.YMin();
	  bbox.maxy = inputVol.YMax();
	  bbox.minz = inputVol.ZMin();
	  bbox.maxz = inputVol.ZMax();

	  VolMagick::Dimension dim;
	  dim.xdim = inputVol.XDim();
	  dim.ydim = inputVol.YDim();
	  dim.zdim = inputVol.ZDim();


//	 cout<<bbox.minz <<" " << bbox.maxz<<" "<< bbox.maxy <<endl;
//	 cout<<dim.zdim <<" " << dim.ydim << endl;
      
      outputVol.voxelType(inputVol.voxelType());
      outputVol.dimension(dim);
      outputVol.boundingBox(bbox);

	
      
   
      //Works for GroEL4.2
      for( int kz = 0; kz<outputVol.ZDim(); kz++)
	   for( int jy = 0; jy<inputVol.YDim(); jy++)
	      for( int ix = 0; ix<inputVol.XDim(); ix++)
		       {
			        float temp=0.0;
					if(inputVol.ZDim() >= kz && inputVol.ZDim() - kz < inputVol.ZDim()
					&& inputVol.YDim() >= jy && inputVol.YDim() - jy < inputVol.YDim() )
					temp = inputVol(ix, inputVol.YDim()-jy, inputVol.ZDim()-kz);
	  				outputVol(ix, jy, kz, temp);
			   }	
	 
	 if(argc==4 && atoi(argv[3])== 1)
	 {
      for( int kz = 0; kz<inputVol.ZDim(); kz++)
	   for( int jy = 0; jy<inputVol.YDim(); jy++)
	      for( int ix = 0; ix<inputVol.XDim(); ix++)
		       {
			   		outputVol(ix, jy, kz, outputVol(ix,jy,kz)+inputVol(ix,jy,kz));
			   }		  
	 }
		  
			  
		 
	      

      VolMagick::createVolumeFile(outputVol, argv[2]);


      cout<<"done!"<<endl;

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
