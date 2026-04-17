/*
  Copyright 2006 The University of Texas at Austin

        Authors: Sangmin Park <smpark@ices.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of PEDetection.

  PEDetection is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  PEDetection is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#ifndef FILE_MEMBRANE_H
#define FILE_MEMBRANE_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#ifndef _WIN32
#include <unistd.h>
#endif
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stream.h>


template <class _DataType> 
class cMembraneSeg
{
	protected:
		_DataType		*Data_mT;
		float			MinData_mf, MaxData_mf;
		int				Width_mi, Height_mi, Depth_mi;
		int				WtimesH_mi, WHD_mi;
		double			ZeroSurfaceValue_md;
		char			*TargetName_mc;
		unsigned char	*ClassData_muc;


	public:
		cMembraneSeg();
		~cMembraneSeg();

	public:
		void setData(_DataType *Data, float Min, float Max);
		void setWHD(int W, int H, int D);
		void setZeroSurface(double ZeroSurface);
		void setFileName(char *FileName);
		_DataType* TopSurfaceSeg();
		_DataType* DownSurfaceSeg();
		_DataType* Class(int ClassNum);
		
		
	public:
		int	Index(int X, int Y, int Z);


};

#endif


