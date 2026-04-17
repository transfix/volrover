/*
  Copyright 2007-2011 The University of Texas at Austin

        Authors: Joe Rivera <transfix@ices.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of VolMagick.

  VolMagick is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  VolMagick is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

/* $Id: Voxels.cpp 4742 2011-10-21 22:09:44Z transfix $ */

#include <VolMagick/Voxels.h>
#include <VolMagick/VoxelOperationStatusMessenger.h>
#include <VolMagick/CompositeFunction.h>
#include <VolMagick/Utility.h>

#include <CVC/App.h>

#include <boost/current_function.hpp>

namespace VolMagick
{
  Voxels::Voxels(const Dimension& d, VoxelType vt) 
    : _dimension(d), _voxelType(vt), _minIsSet(false), _maxIsSet(false), _vosm(vosmDefault),
      _histogramSize(0), _histogramDirty(true)
  {
    try
      {
	_voxels.reset(new unsigned char[XDim()*YDim()*ZDim()*voxelSize()]);
	memset(_voxels.get(),0,XDim()*YDim()*ZDim()*voxelSize());
      }
    catch(std::bad_alloc& e)
      {
	throw MemoryAllocationError("Could not allocate memory for voxels!");
      }
  }

  Voxels::Voxels(const void *v, const Dimension& d, VoxelType vt)
    : _dimension(d), _voxelType(vt), _minIsSet(false), _maxIsSet(false), _vosm(vosmDefault),
      _histogramSize(0), _histogramDirty(true)
  {
    try
      {
	_voxels.reset(new unsigned char[XDim()*YDim()*ZDim()*voxelSize()]);
	memcpy(_voxels.get(),v,XDim()*YDim()*ZDim()*VoxelTypeSizes[_voxelType]);
      }
    catch(std::bad_alloc& e)
      {
	throw MemoryAllocationError("Could not allocate memory for voxels!");
      }
  }

  Voxels::Voxels(const Voxels& v)
    : _dimension(v.dimension()), 
      _voxelType(v.voxelType()), _minIsSet(false), _maxIsSet(false), _vosm(vosmDefault),
      _histogramSize(0), _histogramDirty(true)
  {
    _voxels = v._voxels;
    if(v.minIsSet() && v.maxIsSet())
      {
	min(v.min());
	max(v.max());
      }
  }

  Voxels::~Voxels() { }

  // ---------------------
  // Voxels::dimension
  // ---------------------
  // Purpose:
  //   Changes the dimensions of this Voxels dataset.
  // ---- Change History ----
  // ??/??/2007 -- Joe R. -- Initial implementation.
  // 08/26/2011 -- Joe R. -- Added voxels argument
  void Voxels::dimension(const Dimension& d, boost::shared_array<unsigned char> voxels)
  {
    CVC::ThreadInfo ti(BOOST_CURRENT_FUNCTION);

    if(d.isNull()) throw NullDimension("Null volume dimension.");

    if(dimension() == d && !voxels) return;

    //If we have voxels initialized for us, just use those.  Else make one for ourselves
    //that is the right size for dimension d.
    if(voxels)
      {
        _dimension = d;
        _voxels = voxels;
      }
    else
      {
        Voxels bak(*this); //backup voxels into bak

        //allocate for the new dimension
        try
          {
            //in case this throws...
            boost::shared_array<unsigned char> tmp(new unsigned char[d.xdim*d.ydim*d.zdim*voxelSize()]);
            _voxels = tmp;
          }
        catch(std::bad_alloc& e)
          {
            throw MemoryAllocationError("Could not allocate memory for voxels!");
          }
    
        _dimension = d;
        memset(_voxels.get(),0,XDim()*YDim()*ZDim()*voxelSize());
        
        //copy the voxels back
        for(uint64 k = 0; k < ZDim() && k < bak.ZDim(); k++)
          for(uint64 j = 0; j < YDim() && j < bak.YDim(); j++)
            for(uint64 i = 0; i < XDim() && i < bak.XDim(); i++)
              (*this)(i,j,k, bak(i,j,k));
      }
  }
  
  void Voxels::voxelType(VoxelType vt)
  {
    CVC::ThreadInfo ti(BOOST_CURRENT_FUNCTION);

    if(voxelType() == vt) return;

    Voxels bak(*this); // backup voxels into bak

    //allocate for the new voxel size
    try
      {
	//in case this throws...
	boost::shared_array<unsigned char> tmp(new unsigned char[XDim()*YDim()*ZDim()*VoxelTypeSizes[vt]]);
	_voxels = tmp;
      }
    catch(std::bad_alloc& e)
      {
	throw MemoryAllocationError("Could not allocate memory for voxels!");
      }

    _voxelType = vt;
    memset(_voxels.get(),0,XDim()*YDim()*ZDim()*voxelSize());

    //copy the voxels back
    uint64 len = XDim()*YDim()*ZDim();
    for(uint64 i = 0; i<len; i++)
      (*this)(i,bak(i));
  }

  void Voxels::calcMinMax() const
  {
    CVC::ThreadInfo ti(BOOST_CURRENT_FUNCTION);

    if(_vosm) _vosm->start(this, VoxelOperationStatusMessenger::CalculatingMinMax, ZDim());
    double val;
    size_t len = XDim()*YDim()*ZDim(), i, count=0;
    size_t slice_len = XDim()*YDim();
    if(len == 0) return;
    val = (*this)(0);
    _min = _max = val;

#if 0
    for(i=0; i<len; i++)
      {
	val = (*this)(i);
	if(val < _min) _min = val;
	if(val > _max) _max = val;

	if(_vosm && (i % (XDim()*YDim())) == 0)
	  _vosm->step(this, VoxelOperationStatusMessenger::CalculatingMinMax, count++);
      }
#endif

    switch(voxelType())
      {
      case CVC::UChar:
	{
	  unsigned char v;
	  unsigned char uchar_min = (unsigned char)(_min);
	  unsigned char uchar_max = (unsigned char)(_max);
	  for(i=0; i<len; i++)
	    {
	      v = *((unsigned char *)(_voxels.get()+i*sizeof(unsigned char)));
	      if(v < uchar_min) uchar_min = v;
	      if(v > uchar_max) uchar_max = v;
	      if((i % slice_len) == 0)
                {
                  cvcapp.threadProgress(float(count)/float(ZDim()));
                  if(_vosm) _vosm->step(this, VoxelOperationStatusMessenger::CalculatingMinMax, count);
                  count++;
                }
	    }
	  _min = double(uchar_min);
	  _max = double(uchar_max);
	  break;
	}
      case CVC::UShort:
	{
	  unsigned short v;
	  unsigned short ushort_min = (unsigned short)(_min);
	  unsigned short ushort_max = (unsigned short)(_max);
	  for(i=0; i<len; i++)
	    {
	      v = *((unsigned short *)(_voxels.get()+i*sizeof(unsigned short)));
	      if(v < ushort_min) ushort_min = v;
	      if(v > ushort_max) ushort_max = v;
	      if((i % slice_len) == 0)
                {
                  cvcapp.threadProgress(float(count)/float(ZDim()));
                  if(_vosm) _vosm->step(this, VoxelOperationStatusMessenger::CalculatingMinMax, count);
                  count++;
                }
	    }
	  _min = double(ushort_min);
	  _max = double(ushort_max);
	  break;
	}
      case CVC::UInt:
	{
	  unsigned int v;
	  unsigned int uint_min = (unsigned int)(_min);
	  unsigned int uint_max = (unsigned int)(_max);
	  for(i=0; i<len; i++)
	    {
	      v = *((unsigned int *)(_voxels.get()+i*sizeof(unsigned int)));
	      if(v < uint_min) uint_min = v;
	      if(v > uint_max) uint_max = v;
	      if((i % slice_len) == 0)
                {
                  cvcapp.threadProgress(float(count)/float(ZDim()));
                  if(_vosm) _vosm->step(this, VoxelOperationStatusMessenger::CalculatingMinMax, count);
                  count++;
                }
	    }
	  _min = double(uint_min);
	  _max = double(uint_max);
	  break;
	}
      case CVC::Float:
	{
	  float v;
	  float float_min = (float)(_min);
	  float float_max = (float)(_max);
	  for(i=0; i<len; i++)
	    {
	      v = *((float *)(_voxels.get()+i*sizeof(float)));
	      if(v < float_min) float_min = v;
	      if(v > float_max) float_max = v;
	      if((i % slice_len) == 0)
                {
                  cvcapp.threadProgress(float(count)/float(ZDim()));
                  if(_vosm) _vosm->step(this, VoxelOperationStatusMessenger::CalculatingMinMax, count);
                  count++;
                }
	    }
	  _min = double(float_min);
	  _max = double(float_max);
	  break;
	}
      case CVC::Double:
	{
	  double v;
	  double double_min = (double)(_min);
	  double double_max = (double)(_max);
	  for(i=0; i<len; i++)
	    {
	      v = *((double *)(_voxels.get()+i*sizeof(double)));
	      if(v < double_min) double_min = v;
	      if(v > double_max) double_max = v;
	      if((i % slice_len) == 0)
                {
                  cvcapp.threadProgress(float(count)/float(ZDim()));
                  if(_vosm) _vosm->step(this, VoxelOperationStatusMessenger::CalculatingMinMax, count);
                  count++;
                }
	    }
	  _min = double(double_min);
	  _max = double(double_max);
	  break;
	}
      case CVC::UInt64:
	{
	  uint64 v;
	  uint64 uint64_min = (uint64)(_min);
	  uint64 uint64_max = (uint64)(_max);
	  for(i=0; i<len; i++)
	    {
	      v = *((uint64 *)(_voxels.get()+i*sizeof(uint64)));
	      if(v < uint64_min) uint64_min = v;
	      if(v > uint64_max) uint64_max = v;
	      if((i % slice_len) == 0)
                {
                  cvcapp.threadProgress(float(count)/float(ZDim()));
                  if(_vosm) _vosm->step(this, VoxelOperationStatusMessenger::CalculatingMinMax, count);
                  count++;
                }
	    }
	  _min = double(uint64_min);
	  _max = double(uint64_max);
	  break;
	}
      }

    _minIsSet = _maxIsSet = true;
    cvcapp.threadProgress(1.0f);
    if(_vosm) _vosm->end(this, VoxelOperationStatusMessenger::CalculatingMinMax);
  }

  void Voxels::calcHistogram(uint64 size) const
  {
    CVC::ThreadInfo ti(BOOST_CURRENT_FUNCTION);

    if(!_histogramDirty && _histogramSize == size) return;

    _histogramSize = size;
    _histogram.reset(new uint64[size]);
    memset(_histogram.get(),0,sizeof(uint64)*size);

    if(_vosm) _vosm->start(this, VoxelOperationStatusMessenger::CalculatingHistogram, ZDim());

    for(uint64 k = 0; k<ZDim(); k++)
      {
	for(uint64 j = 0; j<YDim(); j++)
	  for(uint64 i = 0; i<XDim(); i++)
	    {
	      uint64 offset = 
		uint64((((*this)(i,j,k) - min())/(max() - min())) * double(size-1));
	      _histogram[offset]++;
	    }
        cvcapp.threadProgress(float(k)/float(ZDim()));
	if(_vosm) _vosm->step(this, VoxelOperationStatusMessenger::CalculatingHistogram, k);
      }

    cvcapp.threadProgress(1.0f);
    if(_vosm) _vosm->end(this, VoxelOperationStatusMessenger::CalculatingHistogram);

    _histogramDirty = false;
  }

  double Voxels::min(uint64 off_x, uint64 off_y, uint64 off_z,
		     const Dimension& dim) const
  {
    CVC::ThreadInfo ti(BOOST_CURRENT_FUNCTION);

    if(_vosm) _vosm->start(this, VoxelOperationStatusMessenger::CalculatingMin, dim[2]);

    double val;
    uint64 i,j,k;
    val = (*this)(0,0,0);
    for(k=0; k<dim[2]; k++)
      {
	for(j=0; j<dim[1]; j++)
	  for(i=0; i<dim[0]; i++)
	    if(val > (*this)(i+off_x,j+off_y,k+off_z))
	      val = (*this)(i+off_x,j+off_y,k+off_z);
        cvcapp.threadProgress(float(k)/float(dim[2]));
	if(_vosm) _vosm->step(this, VoxelOperationStatusMessenger::CalculatingMin, k);
      }
    
    cvcapp.threadProgress(1.0f);
    if(_vosm) _vosm->end(this, VoxelOperationStatusMessenger::CalculatingMin);
    return val;
  }

  double Voxels::max(uint64 off_x, uint64 off_y, uint64 off_z,
		     const Dimension& dim) const
  {
    CVC::ThreadInfo ti(BOOST_CURRENT_FUNCTION);

    if(_vosm) _vosm->start(this, VoxelOperationStatusMessenger::CalculatingMax, dim[2]);

    double val;
    uint64 i,j,k;
    val = (*this)(0,0,0);
    for(k=0; k<dim[2]; k++)
      {
	for(j=0; j<dim[1]; j++)
	  for(i=0; i<dim[0]; i++)
	    if(val < (*this)(i+off_x,j+off_y,k+off_z))
	      val = (*this)(i+off_x,j+off_y,k+off_z);
        cvcapp.threadProgress(float(k)/float(dim[2]));        
	if(_vosm) _vosm->step(this, VoxelOperationStatusMessenger::CalculatingMax, k);
      }
    
    cvcapp.threadProgress(1.0f);
    if(_vosm) _vosm->end(this, VoxelOperationStatusMessenger::CalculatingMax);
    return val;
  }
  
  Voxels& Voxels::copy(const Voxels& vox)
  {
    if(this == &vox)
      return *this;

    //voxelType(vox.voxelType());
    //dimension(vox.dimension());
    //memcpy(_voxels,*vox,XDim()*YDim()*ZDim()*VoxelTypeSizes[voxelType()]);
    _voxelType = vox._voxelType;
    _dimension = vox._dimension;
    _voxels = vox._voxels;
    if(vox.minIsSet() && vox.maxIsSet())
      {
	min(vox.min());
	max(vox.max());
      }
    else
      unsetMinMax();

    _histogram = vox._histogram;
    _histogramSize = vox._histogramSize;
    _histogramDirty = vox._histogramDirty;

    return *this;
  }

  Voxels& Voxels::sub(uint64 off_x, uint64 off_y, uint64 off_z,
		      const Dimension& subvoldim)
  {
    CVC::ThreadInfo ti(BOOST_CURRENT_FUNCTION);

    if(off_x+subvoldim[0]-1 >= dimension()[0] || 
       off_y+subvoldim[1]-1 >= dimension()[1] || 
       off_z+subvoldim[2]-1 >= dimension()[2])
      throw IndexOutOfBounds("Subvolume offset and/or dimension is out of bounds");

    if(_vosm) _vosm->start(this, VoxelOperationStatusMessenger::SubvolumeExtraction, subvoldim[2]);

    Voxels tmp(*this); // back up this object into tmp

    dimension(subvoldim); // change this object's dimension to the subvolume dimension

    //copy the subvolume voxels
    for(uint64 k=0; k<dimension()[2]; k++)
      {
	for(uint64 j=0; j<dimension()[1]; j++)
	  for(uint64 i=0; i<dimension()[0]; i++)
	    (*this)(i,j,k,tmp(i+off_x,j+off_y,k+off_z));
        cvcapp.threadProgress(float(k)/float(dimension()[2]));
	if(_vosm) _vosm->step(this, VoxelOperationStatusMessenger::SubvolumeExtraction, k);
      }

    cvcapp.threadProgress(1.0f);
    if(_vosm) _vosm->end(this, VoxelOperationStatusMessenger::SubvolumeExtraction);
    return *this;
  }

  Voxels& Voxels::fill(double val)
  {
    return fillsub(0,0,0,dimension(),val);
  }

  Voxels& Voxels::fillsub(uint64 off_x, uint64 off_y, uint64 off_z,
			  const Dimension& subvoldim, double val)
  {
    CVC::ThreadInfo ti(BOOST_CURRENT_FUNCTION);

    if(off_x+subvoldim[0]-1 >= dimension()[0] || 
       off_y+subvoldim[1]-1 >= dimension()[1] || 
       off_z+subvoldim[2]-1 >= dimension()[2])
      throw IndexOutOfBounds("Subvolume offset and/or dimension is out of bounds");

    if(_vosm) _vosm->start(this, VoxelOperationStatusMessenger::Fill, subvoldim[2]);

    for(uint64 k=0; k<subvoldim[2]; k++)
      {
	for(uint64 j=0; j<subvoldim[1]; j++)
	  for(uint64 i=0; i<subvoldim[0]; i++)
	    (*this)(i+off_x,j+off_y,k+off_z,val);
        cvcapp.threadProgress(float(k)/float(subvoldim[2]));
	if(_vosm) _vosm->step(this, VoxelOperationStatusMessenger::Fill, k);
      }

    cvcapp.threadProgress(1.0f);
    if(_vosm) _vosm->end(this, VoxelOperationStatusMessenger::Fill);
    return *this;
  }

  Voxels& Voxels::map(double min_, double max_)
  {
    CVC::ThreadInfo ti(BOOST_CURRENT_FUNCTION);

    if(_vosm) _vosm->start(this, VoxelOperationStatusMessenger::Map, ZDim());
    uint64 len = XDim()*YDim()*ZDim(), count=0;
    for(uint64 i=0; i<len; i++)
      {
	(*this)(i,min_ + (((*this)(i) - min())/(max() - min()))*(max_ - min_));
	if((i % (XDim()*YDim())) == 0)
          {
            cvcapp.threadProgress(float(count)/float(ZDim()));
            if(_vosm) _vosm->step(this, VoxelOperationStatusMessenger::Map, count);
            count++;
          }
      }
    min(min_); max(max_); // set the new min and max
    if(_vosm) _vosm->end(this, VoxelOperationStatusMessenger::Map);
    return *this;
  }

  Voxels& Voxels::resize(const Dimension& newdim)
  {
    CVC::ThreadInfo ti(BOOST_CURRENT_FUNCTION);

    double inSpaceX, inSpaceY, inSpaceZ;
    double val[8];
    uint64 resXIndex = 0, resYIndex = 0, resZIndex = 0;
    uint64 ValIndex[8];
    double xPosition = 0, yPosition = 0, zPosition = 0;
    double xRes = 0, yRes = 0, zRes = 0;
    uint64 i,j,k;
    double x,y,z;

    if(newdim.isNull()) throw NullDimension("Null voxels dimension.");

    if(dimension() == newdim) return *this; //nothing needs to be done

    if(_vosm) _vosm->start(this, VoxelOperationStatusMessenger::Resize, newdim[2]);

    Voxels newvox(newdim,voxelType());

    //we require a dimension of at least 2^3 
    if(newdim < Dimension(2,2,2)) 
      {
	//resize this object as if it was 2x2x2
	resize(Dimension(2,2,2));

	//copy it into newvox
	newvox.copy(*this);

	//change this object's dimension to the real dimension (destroying voxel values, hence the backup)
	dimension(newdim);

	for(k=0; k<ZDim(); k++)
	  for(j=0; j<YDim(); j++)
	    for(i=0; i<XDim(); i++)
	      (*this)(i,j,k,newvox(i,j,k));

	return *this;
      }

    // inSpace calculation
    inSpaceX = (double)(dimension()[0]-1)/(newdim[0]-1);
    inSpaceY = (double)(dimension()[1]-1)/(newdim[1]-1);
    inSpaceZ = (double)(dimension()[2]-1)/(newdim[2]-1);

    for(k = 0; k < newvox.ZDim(); k++)
      {
	z = double(k)*inSpaceZ;
	resZIndex = uint64(z);
	zPosition = z - uint64(z);
	zRes = 1;
	
	for(j = 0; j < newvox.YDim(); j++)
	  {
	    y = double(j)*inSpaceY;
	    resYIndex = uint64(y);
	    yPosition = y - uint64(y);
	    yRes =  1;

	    for(i = 0; i < newvox.XDim(); i++)
	      {
		x = double(i)*inSpaceX;
		resXIndex = uint64(x);
		xPosition = x - uint64(x);
		xRes = 1;

		// find index to get eight voxel values
		ValIndex[0] = resZIndex*dimension()[0]*dimension()[1] + resYIndex*dimension()[0] + resXIndex;
		ValIndex[1] = ValIndex[0] + 1;
		ValIndex[2] = resZIndex*dimension()[0]*dimension()[1] + (resYIndex+1)*dimension()[0] + resXIndex;
		ValIndex[3] = ValIndex[2] + 1;
		ValIndex[4] = (resZIndex+1)*dimension()[0]*dimension()[1] + resYIndex*dimension()[0] + resXIndex;
		ValIndex[5] = ValIndex[4] + 1;
		ValIndex[6] = (resZIndex+1)*dimension()[0]*dimension()[1] + (resYIndex+1)*dimension()[0] + resXIndex;
		ValIndex[7] = ValIndex[6] + 1;

		if(resXIndex>=dimension()[0]-1)
		  {
		    ValIndex[1] = ValIndex[0];
		    ValIndex[3] = ValIndex[2];
		    ValIndex[5] = ValIndex[4];
		    ValIndex[7] = ValIndex[6];
		  }
		if(resYIndex>=dimension()[1]-1)
		  {
		    ValIndex[2] = ValIndex[0];
		    ValIndex[3] = ValIndex[1];
		    ValIndex[6] = ValIndex[4];
		    ValIndex[7] = ValIndex[5];
		  }
		if(resZIndex>=dimension()[2]-1) 
		  {
		    ValIndex[4] = ValIndex[0];
		    ValIndex[5] = ValIndex[1];
		    ValIndex[6] = ValIndex[2];
		    ValIndex[7] = ValIndex[3];
		  }

		for(int Index = 0; Index < 8; Index++) 
		  val[Index] = (*this)(ValIndex[Index]);
		  
		newvox(i,j,k,
		       getTriVal(val, xPosition, yPosition, zPosition, xRes, yRes, zRes));
	      }
	  }

        cvcapp.threadProgress(float(k)/float(newvox.ZDim()));
	if(_vosm) _vosm->step(this, VoxelOperationStatusMessenger::Resize, k);
      }

    copy(newvox); //make this into a copy of the interpolated voxels

    cvcapp.threadProgress(1.0f);
    if(_vosm) _vosm->end(this, VoxelOperationStatusMessenger::Resize);

    return *this;
  }

  Voxels& Voxels::composite(const Voxels& compVox, int64 off_x, int64 off_y, int64 off_z, const CompositeFunction& func)
  {
    CVC::ThreadInfo ti(BOOST_CURRENT_FUNCTION);

    uint64 i,j,k;

    if(_vosm) _vosm->start(this, VoxelOperationStatusMessenger::Composite, compVox.ZDim());

    for(k=0; k<compVox.ZDim(); k++)
      {
	for(j=0; j<compVox.YDim(); j++)
	  for(i=0; i<compVox.XDim(); i++)
	    if((int64(i)+off_x >= 0) && (int64(i)+off_x < int64(XDim())) &&
	       (int64(j)+off_y >= 0) && (int64(j)+off_y < int64(YDim())) &&
	       (int64(k)+off_z >= 0) && (int64(k)+off_z < int64(ZDim())))
	      (*this)(int64(i) + off_x, int64(j) + off_y, int64(k) + off_z,
		      func(compVox,i,j,k,
			   *this, int64(i) + off_x, int64(j) + off_y, int64(k) + off_z));
        cvcapp.threadProgress(float(k)/float(compVox.ZDim()));
	if(_vosm) _vosm->step(this, VoxelOperationStatusMessenger::Composite, k);
      }

    cvcapp.threadProgress(1.0f);
    if(_vosm) _vosm->end(this, VoxelOperationStatusMessenger::Composite);
    
    return *this;
  }
}

