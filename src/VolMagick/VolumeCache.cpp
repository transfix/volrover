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

/* $Id: VolumeCache.cpp 4742 2011-10-21 22:09:44Z transfix $ */

#define VOLMAGICK_VOLUMECACHE_BASIC
#ifdef DEBUG
#pragma message("Remove this file ASAP.  It is not needed now that we have BoundingBox based VolMagick::readVolumeFile()")
#endif

#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/array.hpp>

#ifndef __WINDOWS__
#ifndef VOLMAGICK_VOLUMECACHE_BASIC
#include <openssl/md5.h>
#endif
#endif

#include <VolMagick/VolMagick.h>
#include <VolMagick/VolumeCache.h>



namespace
{
  void append_string(const std::string& filename,const std::string& line)
  {
    using namespace std;
    using namespace boost;

    ofstream outfile(filename.c_str(),ofstream::app);
    if(!outfile)
      throw VolMagick::VolumeCacheDirectoryFileError(
	   str(format("Error appending cache entry to cache directory file %1%") % filename)
	   ); 
    outfile << line;
  }
}

namespace VolMagick
{
  void VolumeCache::add(const VolumeFileInfo& vfi)
  {
    if(!vfi.isSet())
      throw ReadError("Cannot add an uninitialized VolumeFileInfo object to the cache.");
    
    /*
      Make sure that it has the same bounding box and the same number of variables/timesteps
      as the rest of the cache volumes. Else throw an exception.
    */
    for(std::set<VolumeFileInfo, dimcmp>::iterator cur = _volCacheInfo.begin();
	cur != _volCacheInfo.end(); cur++)
      {
	if((*cur).boundingBox() != vfi.boundingBox())
	  throw VolumePropertiesMismatch("All VolumeFileInfo objects in a VolumeCache object "
					 "must have the same bounding box.");
	if((*cur).numVariables() != vfi.numVariables())
	  throw VolumePropertiesMismatch("All VolumeFileInfo objects in a VolumeCache object "
					 "must have the same number of variables.");
	if((*cur).numTimesteps() != vfi.numTimesteps())
	  throw VolumePropertiesMismatch("All VolumeFileInfo objects in a VolumeCache object "
					 "must have the same number of timesteps.");
      }

    if(_volCacheInfo.size() == 0)
      {
	_globalMin.resize(
	  boost::extents[vfi.numVariables()][vfi.numTimesteps()]
	  );
	_globalMax.resize(
	  boost::extents[vfi.numVariables()][vfi.numTimesteps()]
	  );
			  
	for(unsigned int i = 0; i < vfi.numVariables(); i++)
	  for(unsigned int j = 0; j < vfi.numTimesteps(); j++)
	    {
	      _globalMin[i][j] = vfi.min(i,j);
	      _globalMax[i][j] = vfi.max(i,j);
	    }
      }
    else
      {
	for(unsigned int i = 0; i < vfi.numVariables(); i++)
	  for(unsigned int j = 0; j < vfi.numTimesteps(); j++)
	    {
	      if(_globalMin[i][j] > vfi.min(i,j))
		_globalMin[i][j] = vfi.min(i,j);
	      if(_globalMax[i][j] < vfi.max(i,j))
		_globalMax[i][j] = vfi.max(i,j);
	    }
      }
    
    _volCacheInfo.insert(vfi);
  }

  Volume VolumeCache::get(const BoundingBox& requestRegion, unsigned int var, unsigned int time) const
  {
#ifdef VOLMAGICK_VOLUMECACHE_BASIC
    Volume vol;
    if(size()==0) return vol;
    VolumeFileInfo cacheFile = *(_volCacheInfo.begin());
    VolMagick::readVolumeFile(vol,
			      cacheFile.filename(),
			      var, time, requestRegion);
    return vol;
#else
    Dimension curDim, closestDim;
    VolumeFileInfo closestCacheFile;
    Volume retVol;
    
    if(!requestRegion.isWithin(boundingBox()))
      throw SubVolumeOutOfBounds("Request region bounding box must be within the bounding box defined by the cache.");
    
    if(size()>0)
      {
	//initialize the closest cache file info to be the first volume in the cache
	closestCacheFile = *(_volCacheInfo.begin());
	/*
	closestDim = Dimension(uint64(((requestRegion.maxx-requestRegion.minx)/closestCacheFile.XSpan())+1.0),
			       uint64(((requestRegion.maxy-requestRegion.miny)/closestCacheFile.YSpan())+1.0),
			       uint64(((requestRegion.maxz-requestRegion.minz)/closestCacheFile.ZSpan())+1.0));
	*/
	/*
	closestDim = Dimension(ceil((requestRegion.maxx-requestRegion.minx)/closestCacheFile.XSpan())+1,
			       ceil((requestRegion.maxy-requestRegion.miny)/closestCacheFile.YSpan())+1,
			       ceil((requestRegion.maxz-requestRegion.minz)/closestCacheFile.ZSpan())+1);
	*/
	closestDim = Dimension(((requestRegion.maxx-requestRegion.minx)/closestCacheFile.XSpan())+1,
			       ((requestRegion.maxy-requestRegion.miny)/closestCacheFile.YSpan())+1,
			       ((requestRegion.maxz-requestRegion.minz)/closestCacheFile.ZSpan())+1);
	
	for(std::set<VolumeFileInfo, dimcmp>::iterator cur = _volCacheInfo.begin();
	    cur != _volCacheInfo.end(); cur++)
	  {
	    /*
	    curDim = Dimension(uint64(((requestRegion.maxx-requestRegion.minx)/(*cur).XSpan())+1.0),
			       uint64(((requestRegion.maxy-requestRegion.miny)/(*cur).YSpan())+1.0),
			       uint64(((requestRegion.maxz-requestRegion.minz)/(*cur).ZSpan())+1.0));
	    */
	    /*
	    curDim = Dimension(ceil((requestRegion.maxx-requestRegion.minx)/cur->XSpan())+1,
			       ceil((requestRegion.maxy-requestRegion.miny)/cur->YSpan())+1,
			       ceil((requestRegion.maxz-requestRegion.minz)/cur->ZSpan())+1);
	    */
	    curDim = Dimension(((requestRegion.maxx-requestRegion.minx)/cur->XSpan())+1,
			       ((requestRegion.maxy-requestRegion.miny)/cur->YSpan())+1,
			       ((requestRegion.maxz-requestRegion.minz)/cur->ZSpan())+1);
	    
	    //get the closest cache file to the _maxDimension
	    if(abs(int64(curDim.size()) - int64(_maxDimension.size())) <
	       abs(int64(closestDim.size()) - int64(_maxDimension.size())))
	      {
		closestCacheFile = *cur;
		closestDim = curDim;
	      }
	  }
	
	//extract the smallest subvolume within voxel boundaries
	readVolumeFile(retVol,
		       closestCacheFile.filename(),
		       var, time,
		       /*floor*/((requestRegion.minx-boundingBox().minx)/closestCacheFile.XSpan()),
		       /*floor*/((requestRegion.miny-boundingBox().miny)/closestCacheFile.YSpan()),
		       /*floor*/((requestRegion.minz-boundingBox().minz)/closestCacheFile.ZSpan()),
		       closestDim);
	
	//sometimes the request region is slightly larger than retVol (off by a voxel)
	if(!requestRegion.isWithin(retVol.boundingBox()))
	  retVol.boundingBox(requestRegion);

	//extract the exact subvolume and return a volume exactly _maxDimension in size
	if(_maxDimension.size() < closestDim.size())
	  retVol.sub(requestRegion,_maxDimension);

	//force the min and max
	retVol.min(_globalMin[var][time]);
	retVol.max(_globalMax[var][time]);
      }
    
    return retVol;
#endif
  }

  VolumeCache buildCache(const VolumeFileInfo& vfi, const std::string& dir)
  {
#ifdef VOLMAGICK_VOLUMECACHE_BASIC
    VolumeCache cache;
    return cache;
#else
    using namespace std;
    using namespace boost;
    using namespace boost::algorithm;
    using namespace boost::filesystem;
    VolumeCache cache;

    //just disable this for now on windows due to a lack of MD5
#ifndef __WINDOWS__
    //generate an MD5 hash of the filename for our cache filename prefix
    unsigned char hash_buf[16];
    vector<unsigned char> local_filename(vfi.filename().size());
    memcpy(&(local_filename[0]),vfi.filename().c_str(),vfi.filename().size());
    MD5(&(local_filename[0]),
	local_filename.size(),hash_buf);
    string hashed_filename("");
    for(int i = 0; i < 16; i++)
      hashed_filename += str(format("%x") % int(hash_buf[i]));
    
    //hashed_filename += ".rawiv";
    //cout << "hashed filename: " << hashed_filename << endl;

    //read the cache directory file if it exists
    regex expression("^(.*):(.*)$");
    cmatch what;
    map<string, unsigned long long> cache_directory;
    path dir_path(dir);

    if(!exists(dir_path))
      {
	try
	  {
	    create_directory(dir_path);
	  }
	catch(std::exception& e)
	  {
	    throw VolumeCacheDirectoryFileError(
		str(format("Could not create new cache directory %1%")
		    % dir_path.external_file_string())
		);
	  }
      }

    if(exists(dir_path / "directory"))
      {
	ifstream infile(string((dir_path / "directory").external_file_string()).c_str());
	unsigned int linenum = 1;
	if(!infile)
	  throw VolumeCacheDirectoryFileError("Could not open directory file!");
	while(!infile.eof())
	  {
	    string line;
	    getline(infile,line);
	    if(!infile)
	      throw VolumeCacheDirectoryFileError(
		  str(format("Could not read directory file at line %1%") % linenum)
		  );
	    if(regex_match(line.c_str(),what,expression))
	      {
		string file_path(what[1]);
		unsigned long long modified_time = lexical_cast<unsigned long long>(what[2]);
		if(cache_directory.find(file_path) !=
		   cache_directory.end())
		  throw VolumeCacheDirectoryFileError(
		      str(format("Duplicate entry in directory file at line %1%") % linenum)
		      );
		cache_directory[ file_path ] = modified_time;
	      }
	    else
	      throw VolumeCacheDirectoryFileError(
		  str(format("Invalid cache directory file entry at line %1%") % linenum)
		  );
	    linenum++;
	  }
      }
    
    //build cache for vfi if not already present or needs updating
    if(cache_directory.find(vfi.filename()) == cache_directory.end() ||
       cache_directory[vfi.filename()] < static_cast<unsigned long long>(last_write_time(vfi.filename())))
      {
	array<VolMagick::uint64,3> dims = {{ vfi.XDim(), vfi.YDim(), vfi.ZDim() }};
	VolMagick::uint64 maxdim = *max_element(dims.begin(),dims.end());
	
	const VolMagick::uint64 max_chunk_dim_size = 128;
	const VolMagick::Dimension max_chunk_dim(max_chunk_dim_size,
						 max_chunk_dim_size,
						 max_chunk_dim_size);
	
	double xdim_ratio = ceil(max(double(vfi.XDim())/double(max_chunk_dim_size),1.0));
	double ydim_ratio = ceil(max(double(vfi.YDim())/double(max_chunk_dim_size),1.0));
	double zdim_ratio = ceil(max(double(vfi.ZDim())/double(max_chunk_dim_size),1.0));
	double x_interval = (vfi.XMax() - vfi.XMin())/xdim_ratio;
	double y_interval = (vfi.YMax() - vfi.YMin())/ydim_ratio;
	double z_interval = (vfi.ZMax() - vfi.ZMin())/zdim_ratio;

	for(unsigned int var = 0; var < vfi.numVariables(); var++)
	  for(unsigned int t = 0; t < vfi.numTimesteps(); t++)
	    {
	      for(VolMagick::uint64 dim = 2 << 3;
		  dim < maxdim;
		  dim <<= 1)
		{
		  VolMagick::Volume vol;
		  VolMagick::Dimension cur_dim(min(dim,vfi.XDim()),
					       min(dim,vfi.YDim()),
					       min(dim,vfi.ZDim()));
		  
		  string outfilename = str(format("%1%/%2%.%3%.%4%.%5%.rawiv")
					   % dir
					   % hashed_filename
					   % var
					   % t
					   % dim);
		  
		  VolMagick::createVolumeFile(outfilename,
					      vfi.boundingBox(),
					      cur_dim,
					      std::vector<VolMagick::VoxelType>(1,VolMagick::UChar));
		  for(double z_split = vfi.ZMin();
		      z_split < vfi.ZMax();
		      z_split += z_interval)
		    for(double y_split = vfi.YMin();
			y_split < vfi.YMax();
			y_split += y_interval)
		      for(double x_split = vfi.XMin();
			  x_split < vfi.XMax();
			  x_split += x_interval)
			{
			  VolMagick::readVolumeFile(vol,vfi.filename(),
						    var,t,
						    VolMagick::BoundingBox(x_split,y_split,z_split,
									   x_split+x_interval,
									   y_split+y_interval,
									   z_split+z_interval));
			  vol.min(vfi.min());
			  vol.max(vfi.max());
			  vol.map(0.0,255.0);
			  vol.voxelType(VolMagick::UChar);
			  VolMagick::writeVolumeFile(vol,outfilename,
						     var,t,
						     VolMagick::BoundingBox(x_split,y_split,z_split,
									    x_split+x_interval,
									    y_split+y_interval,
									    z_split+z_interval));
			}
		} 
	    }

	cache_directory[ vfi.filename() ] = static_cast<unsigned long long>(time(NULL));
	append_string((dir_path / "directory").external_file_string(),
		      str(format("%1%:%2%\n")
			  % vfi.filename()
			  % cache_directory[vfi.filename()]));
      }

    //now load everything into cache and return it
    regex filename_exp("(.*)\\.rawiv");
    directory_iterator end_itr; // default construction yields past-the-end
    for ( directory_iterator itr( dir_path );
	  itr != end_itr;
	  ++itr )
      {
	string name(itr->leaf());
	if(regex_match(name.c_str(),what,expression))
	  {
	    string fileprefix(what[1]);
	    vector<string> split_prefix;
	    split(split_prefix,fileprefix,is_any_of("."));
	    if(split_prefix[0] == hashed_filename)
	      cache.add(VolumeFileInfo(itr->path().external_file_string()));
	  }
      }
#else
    cache.add(vfi);
#endif

    return cache;
#endif
  }
};
