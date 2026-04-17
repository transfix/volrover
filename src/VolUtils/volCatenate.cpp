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



/*using namespace std;

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
*/
int main(int argc, char **argv)
{
  if(argc < 3)
    {
      std:: cerr << 
	"Usage: " << argv[0] << 

	
	"  <int catVolsNum>  <first volume> [ second volume ...] <output volume> [dimX dimY dimZ]. \n";

      return 1;
    }

  try
    {
  
  //    VolMagick::Volume inputVol2;

      
      VolMagick::Volume outputVol;

	  int volNums = atoi(argv[1]);

	  VolMagick::Dimension dim;

//	  VolMagick::VolumeFileInfo volinfo[volNums];
	  VolMagick::Volume inputVol[volNums];


	  for(int i=0; i < volNums; i++)
	  {
		VolMagick::readVolumeFile(inputVol[i],argv[2+i]); 
	  }

//Set initial bounding box
	  VolMagick::BoundingBox bbox;
	  bbox.minx = inputVol[0].XMin();
	  bbox.maxx = inputVol[0].XMax();
	  bbox.miny = inputVol[0].YMin();
	  bbox.maxy = inputVol[0].YMax();
	  bbox.minz = inputVol[0].ZMin();
	  bbox.maxz = inputVol[0].ZMax();

	  for(int i=1; i < volNums; i++)
	  {
	  	  bbox.minx= bbox.minx>inputVol[i].XMin()? inputVol[i].XMin():bbox.minx;
	  	  bbox.maxx= bbox.maxx<inputVol[i].XMax()? inputVol[i].XMax():bbox.maxx;
	  	  bbox.miny= bbox.miny>inputVol[i].YMin()? inputVol[i].YMin():bbox.miny;
	  	  bbox.maxy= bbox.maxy<inputVol[i].YMax()? inputVol[i].YMax():bbox.maxy;
	  	  bbox.minz= bbox.minz>inputVol[i].ZMin()? inputVol[i].ZMin():bbox.minz;
	  	  bbox.maxz= bbox.maxz<inputVol[i].ZMax()? inputVol[i].ZMax():bbox.maxz;
	  }


  	 float span[3];

//Set dimensions

	  if(argc - 6 == volNums) 
	  {
	  	dim.xdim  = atoi(argv[argc-3]);
		dim.ydim  = atoi(argv[argc-2]);
		dim.zdim  = atoi(argv[argc-1]);

		span[0] =(float) (bbox.maxx - bbox.minx)/(dim.xdim-1);
		span[1] =(float) (bbox.maxy - bbox.miny)/(dim.ydim-1);
		span[2] =(float) (bbox.maxz - bbox.minz)/(dim.zdim-1);

	  }
	  else
	  {	 
		span[0] = inputVol[0].XSpan();
	  	span[1] = inputVol[0].YSpan();
	  	span[2] = inputVol[0].ZSpan();
		for(int i=1; i< volNums; i++)
		{
			span[0] = span[0]>inputVol[i].XSpan()?inputVol[i].XSpan():span[0];
			span[1] = span[1]>inputVol[i].YSpan()?inputVol[i].YSpan():span[1];
			span[2] = span[2]>inputVol[i].ZSpan()?inputVol[i].ZSpan():span[2];
		}
		dim.xdim = (int)((bbox.maxx - bbox.minx)/span[0])+1;
		dim.ydim = (int)((bbox.maxy - bbox.miny)/span[1])+1;
		dim.zdim = (int)((bbox.maxz - bbox.minz)/span[2])+1;
	  }


	  bbox.maxx = bbox.minx + (dim.xdim-1) * span[0];
	  bbox.maxy = bbox.miny + (dim.ydim-1) * span[1];
	  bbox.maxz = bbox.minz + (dim.zdim-1) * span[2];
	  
	  outputVol.boundingBox(bbox);
	  outputVol.dimension(dim);


      float volmin = inputVol[0].min();
	  for(int i=0; i < volNums; i++)
	  {
	  	volmin = volmin > inputVol[i].min()? inputVol[i].min():volmin;
	   }

	  outputVol.voxelType(inputVol[0].voxelType());
      
	  for(int kz = 0; kz < outputVol.ZDim(); kz ++)
	    for(int jy = 0; jy < outputVol.YDim(); jy ++)
		   for(int ix = 0; ix < outputVol.XDim(); ix ++)
		     outputVol(ix, jy, kz, volmin);

	  float x,y,z;
	  
	  for(int kz = 0; kz < outputVol.ZDim(); kz ++)
	    for(int jy = 0; jy < outputVol.YDim(); jy ++)
		   for(int ix = 0; ix < outputVol.XDim(); ix ++)
		   {
		   		x = bbox.minx + ix * span[0];
				y = bbox.miny + jy * span[1];
				z = bbox.minz + kz * span[2];

				float temp = 0.0;
				int t=0;
				for(int i=0; i< volNums; i++)
				{
					if(x>inputVol[i].XMin() && x<inputVol[i].XMax() 
						&& y>inputVol[i].YMin() && y<inputVol[i].YMax() 
						&& z>inputVol[i].ZMin() && z<inputVol[i].ZMax())
					{
						temp += inputVol[i].interpolate(x,y,z);
						t++;
					}
				}
				if(t>0) outputVol(ix,jy,kz, (float)(temp/t));

		   }


 

      VolMagick::createVolumeFile(outputVol, argv[2+volNums]);





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
