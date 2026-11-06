/*
  Copyright 2007-2011 The University of Texas at Austin
  
	Authors: Jose Rivera <transfix@ices.utexas.edu>
	Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of Volume Rover.

  Volume Rover is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  Volume Rover is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

/* $Id: HDF5_IO.cpp 5536 2012-05-07 22:12:01Z arand $ */

#include <VolMagick/HDF5_IO.h>

#include <VolMagick/VolMagick.h>
#include <VolMagick/endians.h>

#include <CVC/upToPowerOfTwo.h>

#if defined(WIN32)
#include <cpp/H5PredType.h>
#else
#include <H5PredType.h>
#endif

#include <boost/current_function.hpp>
#include <boost/format.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/scoped_array.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/array.hpp>

#include <iostream>
#include <algorithm>
#include <limits>
#include <cmath>

namespace
{
  // --------------
  // BuildHierarchy
  // --------------
  // Purpose:
  //   Thread for building the multi-res hierarchy for a volume dataset
  // ---- Change History ----
  // 09/09/2011 -- Joe R. -- Initial implementation.
  // 09/11/2011 -- Joe R. -- Finishes hierarchy calculation without throwing an exception
  // 09/18/2011 -- Joe R. -- Using property volmagick.hdf5_io.buildhierarchy.current to report
  //                         info about the hierarchy dataset currently being processed.  This
  //                         gives clients a chance to update themselves when more data is
  //                         available.
  // 09/30/2011 -- Joe R. -- Counting the number of steps for useful thread progress reporting.
  //                         Also added chunk_size property.
  // 10/09/2011 -- Joe R. -- Using startThread
  class BuildHierarchy
  {
  public:
    BuildHierarchy(const std::string& threadKey,
                   const std::string& hdf5_filename,
                   const std::string& hdf5_volumeDataSet)
      : _threadKey(threadKey),
        _hdf5_filename(hdf5_filename),
        _hdf5_volumeDataSet(hdf5_volumeDataSet) {}

    BuildHierarchy(const BuildHierarchy& t)
      : _threadKey(t._threadKey),
        _hdf5_filename(t._hdf5_filename),
        _hdf5_volumeDataSet(t._hdf5_volumeDataSet) {}

    BuildHierarchy& operator=(const BuildHierarchy& t)
    {
      _threadKey = t._threadKey;
      _hdf5_filename = t._hdf5_filename;
      _hdf5_volumeDataSet = t._hdf5_volumeDataSet;
    }

    //lazy way to count the number of steps
    VolMagick::uint64 countNumSteps(const VolMagick::Dimension& fullDim, 
                                    const VolMagick::BoundingBox& bbox)
    {
      using namespace VolMagick;

      Dimension prevDim(fullDim);

      uint64 steps = 0;
      while(1)
        {
          VolMagick::uint64 maxdim = std::max(prevDim.xdim,std::max(prevDim.ydim,prevDim.zdim));
          maxdim = CVC::upToPowerOfTwo(maxdim) >> 1; //power of 2 less than maxdim
          Dimension curDim(maxdim,maxdim,maxdim);
          for(int i = 0; i < 3; i++)
            curDim[i] = std::min(curDim[i],prevDim[i]);
          if(curDim.size()==1) break; //we're done if the dims hit 1
          
          {
            Dimension targetDim = curDim;
            const uint64 maxdim_size = 
              cvcapp.properties<VolMagick::uint64>("volmagick.hdf5_io.buildhierarchy.chunk_size");
            boost::array<double,3> theSize =
              {
                maxdim_size*bbox.XSpan(fullDim),
                maxdim_size*bbox.YSpan(fullDim),
                maxdim_size*bbox.ZSpan(fullDim)
              };
            for(double off_z = bbox.minz;
                off_z < bbox.maxz;
                off_z += theSize[2])
              for(double off_y = bbox.miny;
                  off_y < bbox.maxy;
                  off_y += theSize[1])
                for(double off_x = bbox.minx;
                    off_x < bbox.maxx;
                    off_x += theSize[0])
                  steps++;
          }

          prevDim = curDim;
        }
      return steps;
    }

    void operator()()
    {
      using namespace VolMagick;
      using namespace CVC::HDF5_Utilities;
      using namespace boost;

      CVC::ThreadFeedback feedback;

      //read/write 128^3 chunks by default
      if(!cvcapp.hasProperty("volmagick.hdf5_io.buildhierarchy.chunk_size"))
        cvcapp.properties("volmagick.hdf5_io.buildhierarchy.chunk_size",uint64(128));

      //Sleep for a second before beginning so we don't thrash about
      //if writeVolumeFile is called several times in succession.
      {
        CVC::ThreadInfo ti("sleeping");
        boost::xtime xt;
        boost::xtime_get( &xt, boost::TIME_UTC_ );
        xt.sec++;
        boost::thread::sleep( xt );
      }

      try
        {
          Dimension fullDim = getObjectDimension(_hdf5_filename,_hdf5_volumeDataSet);
          BoundingBox bbox = getObjectBoundingBox(_hdf5_filename,_hdf5_volumeDataSet);
          Dimension prevDim(fullDim);

          uint64 numSteps = countNumSteps(fullDim,bbox);
          uint64 steps = 0;

          while(1)
            {
              uint64 maxdim = std::max(prevDim.xdim,std::max(prevDim.ydim,prevDim.zdim));
              maxdim = CVC::upToPowerOfTwo(maxdim) >> 1; //power of 2 less than maxdim
              Dimension curDim(maxdim,maxdim,maxdim);
              for(int i = 0; i < 3; i++)
                curDim[i] = std::min(curDim[i],prevDim[i]);
              if(curDim.size()==1) break; //we're done if the dims hit 1

              std::string hier_volume_name =
                str(format("%1%_%2%x%3%x%4%") 
                    % _hdf5_volumeDataSet
                    % curDim[0] % curDim[1] % curDim[2]);

              int isDirty = 0;
              getAttribute(_hdf5_filename, hier_volume_name, "dirty", isDirty);
              if(isDirty)
                {
                  cvcapp.log(1,str(format("%1% :: computing %2%\n")
                                   % BOOST_CURRENT_FUNCTION
                                   % hier_volume_name));

                  {
                    Dimension targetDim = curDim;
                    const uint64 maxdim_size = 
                      cvcapp.properties<uint64>("volmagick.hdf5_io.buildhierarchy.chunk_size");
                    boost::array<double,3> theSize =
                      {
                        maxdim_size*bbox.XSpan(fullDim),
                        maxdim_size*bbox.YSpan(fullDim),
                        maxdim_size*bbox.ZSpan(fullDim)
                      };
                    for(double off_z = bbox.minz;
                        off_z < bbox.maxz;
                        off_z += theSize[2])
                      {
                        for(double off_y = bbox.miny;
                            off_y < bbox.maxy;
                            off_y += theSize[1])
                          for(double off_x = bbox.minx;
                              off_x < bbox.maxx;
                              off_x += theSize[0])
                            {
                              VolMagick::Volume vol;
                              VolMagick::BoundingBox subvolbox(
                                off_x,off_y,off_z,
                                std::min(off_x+theSize[0],bbox.maxx),
                                std::min(off_y+theSize[1],bbox.maxy),
                                std::min(off_z+theSize[2],bbox.maxz)
                              );
                              VolMagick::readVolumeFile(
                                vol,
                                _hdf5_filename + "|" + _hdf5_volumeDataSet,
                                0,0,subvolbox
                              );
                              VolMagick::writeVolumeFile(
                                vol,
                                _hdf5_filename + "|" + hier_volume_name,
                                0,0,subvolbox
                              );
                              cvcapp.threadProgress(float(++steps)/float(numSteps));
                            }
                      }
                  }

                  //Done, now mark this one clean.
                  setAttribute(_hdf5_filename, hier_volume_name, "dirty", 0);

                  cvcapp.properties("volmagick.hdf5_io.buildhierarchy.latest",
                                    _hdf5_filename + "|" + hier_volume_name);

                }
              else
                {
                  cvcapp.log(1,str(format("%1% :: %2% not dirty\n")
                                   % BOOST_CURRENT_FUNCTION
                                   % hier_volume_name));
                }

              prevDim = curDim;
            }

          cvcapp.threadProgress(1.0f);
        }
      catch(boost::thread_interrupted&)
        {
          cvcapp.log(6,str(format("%1% :: thread %2% interrupted\n")
                           % BOOST_CURRENT_FUNCTION
                           % cvcapp.threadKey()));
        }
      catch(VolMagick::Exception& e)
        {
          cvcapp.log(1,str(format("%1% :: ERROR :: %2%\n")
                           % BOOST_CURRENT_FUNCTION
                           % e.what()));
        }
    }

    static void start(const std::string& threadKey,
                      const std::string& hdf5_filename,
                      const std::string& hdf5_volumeDataSet)
    {
      cvcapp.startThread(threadKey,
                         BuildHierarchy(threadKey,
                                        hdf5_filename,
                                        hdf5_volumeDataSet));
    }

  protected:
    std::string _threadKey;
    std::string _hdf5_filename;
    std::string _hdf5_volumeDataSet;
  };
}

namespace VolMagick
{
  // ----------------
  // HDF5_IO::HDF5_IO
  // ----------------
  // Purpose:
  //   Initializes the extension list and id.
  // ---- Change History ----
  // 12/04/2009 -- Joe R. -- Initial implementation.
  // 01/04/2009 -- Joe R. -- Adding maxdim arg used with the bounding
  //                         box version of readVolumeFile
  // 09/17/2011 -- Joe R. -- Maxdim is now on the property map.
  // 09/30/2011 -- Joe R. -- Checking that the maxdim property doesn't exist
  //                         before setting it.
  HDF5_IO::HDF5_IO()
    : _id("HDF5_IO : v1.0")
  {
    _extensions.push_back(".h5");
    _extensions.push_back(".hdf5");
    _extensions.push_back(".hdf");
    _extensions.push_back(".cvc");

    if(!cvcapp.hasProperty("volmagick.hdf5_io.maxdim"))
      cvcapp.properties("volmagick.hdf5_io.maxdim","128,128,128");
  }

  // -----------
  // HDF5_IO::id
  // -----------
  // Purpose:
  //   Returns a string that identifies this VolumeFile_IO object.  This should
  //   be unique, but is freeform.
  // ---- Change History ----
  // 12/04/2009 -- Joe R. -- Initial implementation.
  const std::string& HDF5_IO::id() const
  {
    return _id;
  }

  const VolumeFile_IO::ExtensionList& HDF5_IO::extensions() const
  {
    return _extensions;
  }

  // --------------------------
  // HDF5_IO::getVolumeFileInfo
  // --------------------------
  // Purpose:
  //   Returns information about the volume file
  // ---- Change History ----
  // 12/04/2009 -- Joe R. -- Initial implementation.
  // 07/24/2011 -- Joe R. -- Using CVC::HDF5_Utilities.
  // 07/30/2011 -- Joe R. -- Fixed some issues with the move to HDF5_Utilities
  // 09/02/2011 -- Joe R. -- Forgot to copy filename to data.
  // 09/09/2011 -- Joe R. -- Adding support for ungrouped, lone datasets to make
  //                         multi-res hierarchy thread code simpler.
  void HDF5_IO::getVolumeFileInfo(VolumeFileInfo::Data& data,
				  const std::string& filename) const
  {
    using namespace CVC::HDF5_Utilities;
    using namespace boost;

    CVC::ThreadInfo ti(BOOST_CURRENT_FUNCTION);

    std::string actualFileName;
    std::string objectName;

    boost::tie(actualFileName, objectName) =
      splitRawFilename(filename);

    data._filename = filename;

    //check if it is an old style volmagick cvc-hdf5 file
    bool oldVolMagick = false;
    try
      {
        uint64 vm_version;
        getAttribute(actualFileName,objectName,"VolMagick_version",vm_version);
        oldVolMagick = true;
      }
    catch(HDF5Exception &)
      {
        std::string version;
        getAttribute(actualFileName,objectName,"libcvc_version",version);
        oldVolMagick = false;
      }

    uint64 numVariables, numTimesteps;
    if(!oldVolMagick)
      {
        cvcapp.log(6,str(format("%s :: using new (Aug2011) VolMagick cvc-hdf5 format\n")
                         % BOOST_CURRENT_FUNCTION));

        data._boundingBox = getObjectBoundingBox(actualFileName,objectName);
        data._dimension   = getObjectDimension(actualFileName,objectName);

        getAttribute(actualFileName,objectName,"numVariables",numVariables);
        getAttribute(actualFileName,objectName,"numTimesteps",numTimesteps);
        data._numVariables = numVariables;
        data._numTimesteps = numTimesteps;
        
        getAttribute(actualFileName,objectName,"min_time",data._tmin);
        getAttribute(actualFileName,objectName,"max_time",data._tmax);
      }
    else
      {
        cvcapp.log(5,str(format("%s :: using old VolMagick cvc-hdf5 format\n")
                         % BOOST_CURRENT_FUNCTION));

        //some older files will have attributes named like this...
        getAttribute(actualFileName,objectName,"XMin",data._boundingBox.minx);
        getAttribute(actualFileName,objectName,"YMin",data._boundingBox.miny);
        getAttribute(actualFileName,objectName,"ZMin",data._boundingBox.minz);
        getAttribute(actualFileName,objectName,"XMax",data._boundingBox.maxx);
        getAttribute(actualFileName,objectName,"YMax",data._boundingBox.maxy);
        getAttribute(actualFileName,objectName,"ZMax",data._boundingBox.maxz);

        getAttribute(actualFileName,objectName,"XDim",data._dimension.xdim);
        getAttribute(actualFileName,objectName,"YDim",data._dimension.ydim);
        getAttribute(actualFileName,objectName,"ZDim",data._dimension.zdim);

        getAttribute(actualFileName,objectName,"numVariables",numVariables);
        getAttribute(actualFileName,objectName,"numTimesteps",numTimesteps);
        data._numVariables = numVariables;
        data._numTimesteps = numTimesteps;

        getAttribute(actualFileName,objectName,"min_time",data._tmin);
        getAttribute(actualFileName,objectName,"max_time",data._tmax);
      }

    data._minIsSet.resize(numVariables);
    data._min.resize(numVariables);
    data._maxIsSet.resize(numVariables);
    data._max.resize(numVariables);
    data._voxelTypes.resize(numVariables);
    data._names.resize(numVariables);

    //The current HDF5_IO implementation prefers storing 3D datasets in groups.
    //Check if the object is a group.  If not, read it as a 1 var 1 timestep lone dataset.
    if(isGroup(actualFileName,objectName))
      {
        for(unsigned int i = 0; i < data._numVariables; i++)
          {
            data._minIsSet[i].resize(numTimesteps);
            data._min[i].resize(numTimesteps);
            data._maxIsSet[i].resize(numTimesteps);
            data._max[i].resize(numTimesteps);
	    
            for(unsigned int j = 0; j < data._numTimesteps; j++)
              {
                std::string volume_name =
                  boost::str(
                             boost::format("%1%/%2%:%3%:%4%") %
                             objectName %
                             DEFAULT_VOLUME_NAME %
                             i % j
                             );

                data._min[i][j] = getDataSetMinimum(actualFileName,volume_name);
                data._minIsSet[i][j]=true;
                data._max[i][j] = getDataSetMaximum(actualFileName,volume_name);
                data._maxIsSet[i][j]=true;
                data._names[i] = getDataSetInfo(actualFileName,volume_name);
                uint64 voxelType;
                if(oldVolMagick)
                  getAttribute(actualFileName,volume_name,"voxelType",voxelType);
                else
                  getAttribute(actualFileName,volume_name,"dataType",voxelType);
                data._voxelTypes[i]=VoxelType(voxelType);
              }
          }
      }
    else
      {
        if(numVariables > 1 || numTimesteps > 1)
          cvcapp.log(1,str(format("%s :: WARNING - HDF5_IO doesn't yet support lone datasets "
                                  "with more than one variable or timestep!\n")
                           % BOOST_CURRENT_FUNCTION));

        std::string volume_name = objectName;
        unsigned int i = 0, j = 0;

        numVariables = 1; numTimesteps = 1;
        data._minIsSet.resize(numVariables);
        data._min.resize(numVariables);
        data._maxIsSet.resize(numVariables);
        data._max.resize(numVariables);
        data._voxelTypes.resize(numVariables);
        data._names.resize(numVariables);
        data._minIsSet[i].resize(numTimesteps);
        data._min[i].resize(numTimesteps);
        data._maxIsSet[i].resize(numTimesteps);
        data._max[i].resize(numTimesteps);

        data._min[i][j] = getDataSetMinimum(actualFileName,volume_name);
        data._minIsSet[i][j]=true;
        data._max[i][j] = getDataSetMaximum(actualFileName,volume_name);
        data._maxIsSet[i][j]=true;
        data._names[i] = getDataSetInfo(actualFileName,volume_name);
        uint64 voxelType;
        if(oldVolMagick)
          getAttribute(actualFileName,volume_name,"voxelType",voxelType);
        else
          getAttribute(actualFileName,volume_name,"dataType",voxelType);
        data._voxelTypes[i]=VoxelType(voxelType);
      }
  }

  // -----------------------
  // HDF5_IO::readVolumeFile
  // -----------------------
  // Purpose:
  //   Reads from a volume file
  // ---- Change History ----
  // 12/04/2009 -- Joe R. -- Initial implementation.
  // 08/05/2011 -- Joe R. -- Using HDF5 Utilities now.
  // 09/09/2011 -- Joe R. -- Adding support for ungrouped, lone datasets to make
  //                         multi-res hierarchy thread code simpler.
  void HDF5_IO::readVolumeFile(Volume& vol,
			       const std::string& filename, 
			       unsigned int var, unsigned int time,
			       uint64 off_x, uint64 off_y, uint64 off_z,
			       const Dimension& subvoldim) const
  {
    using namespace H5;
    using namespace CVC::HDF5_Utilities;
    using namespace boost;

    CVC::ThreadInfo ti(BOOST_CURRENT_FUNCTION);

    std::string actualFileName;
    std::string objectName;

    if(subvoldim.isNull())
      throw HDF5Exception("Null subvoldim");

    boost::tie(actualFileName, objectName) =
      splitRawFilename(filename);

    std::string volume_name;

    if(isGroup(actualFileName,objectName))
      {
        //Name of actual dataset.  The group name in 'objectName' can contain several
        //instances of the same dataset at various resolutions.  This function assumes
        //you want the highest resolution dataset, which uses the following naming convention.
        volume_name =
          boost::str(
            boost::format("%1%/%2%:%3%:%4%") %
            objectName %
            DEFAULT_VOLUME_NAME %
            var % time
          );
      }
    else //Ungrouped dataset
      {
        if(var > 0 || time > 0)
          cvcapp.log(1,str(format("%s :: WARNING - HDF5_IO doesn't yet support lone datasets "
                                  "with more than one variable or timestep!\n")
                           % BOOST_CURRENT_FUNCTION));
        var = 0; time = 0;
        volume_name = objectName;
      }
    
    VolumeFileInfo vfi(filename);
    BoundingBox boundingBox = vfi.boundingBox();
    Dimension dimension = vfi.dimension();
    
    if(off_x + subvoldim[0] > dimension[0] ||
       off_y + subvoldim[1] > dimension[1] ||
       off_z + subvoldim[2] > dimension[2])
      throw InvalidHDF5File("Dimension and offset out of bounds");

    //calculate the subvolume boundingbox for the requested region
    BoundingBox subbbox
      ( 
       boundingBox.minx + 
       (boundingBox.maxx - boundingBox.minx)*
       (off_x/(dimension.xdim-1)),
       boundingBox.miny + 
       (boundingBox.maxy - boundingBox.miny)*
       (off_y/(dimension.ydim-1)),
       boundingBox.minz + 
       (boundingBox.maxz - boundingBox.minz)*
       (off_z/(dimension.zdim-1)),
       boundingBox.minx + 
       (boundingBox.maxx - boundingBox.minx)*
       ((off_x+subvoldim.xdim)/(dimension.xdim-1)),
       boundingBox.miny + 
       (boundingBox.maxy - boundingBox.miny)*
       ((off_y+subvoldim.ydim)/(dimension.ydim-1)),
       boundingBox.minz + 
       (boundingBox.maxz - boundingBox.minz)*
       ((off_z+subvoldim.zdim)/(dimension.zdim-1))
       );
    
    vol.voxelType(vfi.voxelTypes(var));
    vol.desc(vfi.name(var));
    vol.boundingBox(subbbox);
    vol.dimension(subvoldim);

    //set the vol's min/max if it is to be the same size as the 
    //volume in the file.
    if(off_x == 0 && off_y == 0 && off_z == 0 &&
       vfi.dimension() == subvoldim)
      {
        vol.min(getDataSetMinimum(actualFileName,volume_name));
        vol.max(getDataSetMaximum(actualFileName,volume_name));
      }
    
    readDataSet(actualFileName,volume_name,
                off_x,off_y,off_z,
                subvoldim,
                vol.voxelType(),
                *vol);
  }

  // -----------------------------
  // HDF5_IO::readVolumeFile
  // -----------------------------
  // Purpose:
  //   Same as above except uses a bounding box for specifying the
  //   subvol.  Uses maxdim to define a stride to use when reading
  //   for subsampling.
  // ---- Change History ----
  // 01/04/2010 -- Joe R. -- Initial implementation.
  // 08/26/2011 -- Joe R. -- Using HDF5 Utilities now.
  // 09/09/2011 -- Joe R. -- Adding support for ungrouped, lone datasets to make
  //                         multi-res hierarchy thread code simpler.
  // 09/17/2011 -- Joe R. -- Picking out the closest dimension to the maxdim in
  //                         the hierarchy.
  void HDF5_IO::readVolumeFile(Volume& vol, 
			       const std::string& filename, 
			       unsigned int var,
			       unsigned int time,
			       const BoundingBox& subvolbox) const
  {
    using namespace CVC::HDF5_Utilities;
    using namespace boost;

    CVC::ThreadInfo ti(BOOST_CURRENT_FUNCTION);

    std::string actualFileName;
    std::string objectName;

    boost::tie(actualFileName, objectName) =
      splitRawFilename(filename);

    std::string volume_name;
    if(isGroup(actualFileName,objectName))
      {
        //The group name in 'objectName' can contain several instances of the same dataset
        //at various resolutions.

        //get the maximum dimensions to extract
        std::vector<uint64> maxdim_vec = cvcapp.listProperty<uint64>("volmagick.hdf5_io.maxdim");
        while(maxdim_vec.size() < 3)
          maxdim_vec.push_back(128);
        CVC::Dimension maxdim(maxdim_vec);

        // ---
        // Find the dataset with size closest to the maxdim
        // ---
        
        //filter out other variables and timesteps
        std::string filter = str(format("%1%:%2%") % var % time);

        std::vector<std::string> hierarchy_objects = 
          getChildObjects(actualFileName, objectName, filter);
        if(hierarchy_objects.empty())
          throw HDF5Exception(str(format("%s :: no child objects!")
                                  % BOOST_CURRENT_FUNCTION));
        
        CVC::Dimension dim = 
          getDataSetDimensionForBoundingBox(actualFileName,
                                            objectName+"/"+hierarchy_objects[0],
                                            subvolbox);
        std::string hierarchy_object = hierarchy_objects[0];

        //find the first non dirty dataset
        for (const auto& obj : hierarchy_objects) {
            int isDirty = 0;
            try
              {
                getAttribute(actualFileName,
                             objectName+"/"+obj,
                             "dirty",isDirty);
              }
            catch(std::exception&){}
            if(!isDirty)
              {
                hierarchy_object = obj;
                break;
              }
          }

        //now select a dataset
        for (const auto& obj : hierarchy_objects) {
            cvcapp.log(3,str(format("%s: %s\n")
                             % BOOST_CURRENT_FUNCTION
                             % obj));

            //if this object is dirty, lets skip it for now.
            int isDirty = 0;
            try
              {
                getAttribute(actualFileName,
                             objectName+"/"+obj,
                             "dirty",isDirty);
              }
            catch(std::exception&){}

            CVC::Dimension newdim = 
              getDataSetDimensionForBoundingBox(actualFileName,
                                                objectName+"/"+obj,
                                                subvolbox);
            if(!isDirty &&
               std::abs(int64(newdim.size()) - int64(maxdim.size())) <
               std::abs(int64(dim.size()) - int64(maxdim.size())))
              {
                dim = newdim;
                hierarchy_object = obj;
              }
          }

        cvcapp.log(2,str(format("%s: selected object %s\n")
                         % BOOST_CURRENT_FUNCTION
                         % hierarchy_object));

        volume_name = objectName+"/"+hierarchy_object;
      }
    else //Ungrouped dataset
      {
        if(var > 0 || time > 0)
          cvcapp.log(1,str(format("%s :: WARNING - HDF5_IO doesn't yet support lone datasets "
                                  "with more than one variable or timestep!\n")
                           % BOOST_CURRENT_FUNCTION));
        var = 0; time = 0;
        volume_name = objectName;
      }

    VolumeFileInfo vfi(filename);

    vol.voxelType(vfi.voxelTypes(var));
    vol.desc(vfi.name(var));
    vol.boundingBox(subvolbox);

    //set the vol's min/max if it is to be the same size as the 
    //volume in the file.
    if(vfi.boundingBox() == subvolbox)
      {
        vol.min(getDataSetMinimum(actualFileName,volume_name));
        vol.max(getDataSetMaximum(actualFileName,volume_name));
      }

    boost::shared_array<unsigned char> data;
    Dimension dim;
    
    boost::tie(data,dim) = 
      readDataSet(actualFileName, volume_name,
                  subvolbox, vol.voxelType());
    
    vol.dimension(dim,data);
  }

  // -------------------------
  // HDF5_IO::createVolumeFile
  // -------------------------
  // Purpose:
  //   Creates an empty volume file to be later filled in by writeVolumeFile
  // ---- Change History ----
  // 12/04/2009 -- Joe R. -- Initial implementation.
  // 12/28/2009 -- Joe R. -- HDF5 file schema implementation
  // 08/26/2011 -- Joe R. -- Using HDF5 Utilities now.
  // 09/02/2011 -- Joe R. -- Calling new createHDF5File function.
  // 09/08/2011 -- Joe R. -- Only creating a file if none exists, else just re-using it.
  // 09/30/2011 -- Joe R. -- Added objectType attribute.
  void HDF5_IO::createVolumeFile(const std::string& filename,
				 const BoundingBox& boundingBox,
				 const Dimension& dimension,
				 const std::vector<VoxelType>& voxelTypes,
				 unsigned int numVariables, unsigned int numTimesteps,
				 double min_time, double max_time) const
  {
    using namespace CVC::HDF5_Utilities;
    namespace fs = boost::filesystem;
    using boost::filesystem::path;

    CVC::ThreadInfo ti(BOOST_CURRENT_FUNCTION);

    std::string actualFileName;
    std::string objectName;

    boost::tie(actualFileName, objectName) =
      splitRawFilename(filename);

    if(voxelTypes.empty())
      throw InvalidHDF5File("VoxelTypes array not large enough!");

    //create it if it doesn't exist, else just reuse the file.
    fs::path full_path(actualFileName);
    if(!fs::exists(full_path))
       createHDF5File(actualFileName);

    createGroup(actualFileName, objectName, true);
    setAttribute(actualFileName, objectName, "objectType", "VolMagick::Volume");
    setAttribute(actualFileName, objectName, "libcvc_version", VOLMAGICK_VERSION_STRING);
    setObjectBoundingBox(actualFileName, objectName, boundingBox);
    setObjectDimension(actualFileName, objectName, dimension);
    setAttribute(actualFileName, objectName, "dataTypes", voxelTypes.size(), &(voxelTypes[0]));
    setAttribute(actualFileName, objectName, "numVariables", numVariables);
    setAttribute(actualFileName, objectName, "numTimesteps", numTimesteps);
    setAttribute(actualFileName, objectName, "min_time", min_time);
    setAttribute(actualFileName, objectName, "max_time", max_time);
    
    unsigned int steps = numVariables*numTimesteps;
    unsigned int cur_step = 0;
    for(unsigned int var = 0; var < numVariables; var++)
      for(unsigned int time = 0; time < numTimesteps; time++)
        {
          cvcapp.threadProgress(float(cur_step)/float(steps));
          
          //Name of actual dataset.  The group name in 'objectName' can contain several
          //instances of the same dataset at various resolutions.  This function just makes space
          //for the highest resolution dataset.  writeVolumeFile should trigger updates to the
          //hierarchy via a seperate thread for use by the BoundingBox based readVolumeFile above.
          std::string volume_name =
            boost::str(
                boost::format("%1%/%2%:%3%:%4%") %
                objectName %
                DEFAULT_VOLUME_NAME %
                var % time
              );

          if(voxelTypes.size() <= var)
            throw InvalidHDF5File("VoxelTypes array not large enough!");

          createDataSet(actualFileName, volume_name, 
                        boundingBox, dimension, voxelTypes[var]);
        }

    cvcapp.threadProgress(1.0);
  }

  // ------------------------
  // HDF5_IO::writeVolumeFile
  // ------------------------
  // Purpose:
  //   Writes the volume contained in wvol to the specified volume file. Should create
  //   a volume file if the filename provided doesn't exist.  Else it will simply
  //   write data to the existing file.  A common user error arises when you try to
  //   write over an existing volume file using this function for unrelated volumes.
  //   If what you desire is to overwrite an existing volume file, first run
  //   createVolumeFile to replace the volume file.
  // ---- Change History ----
  // 12/04/2009 -- Joe R. -- Initial implementation.
  // 08/28/2011 -- Joe R. -- Using HDF5 Utilities now.
  // 09/09/2011 -- Joe R. -- Adding support for ungrouped, lone datasets to make
  //                         multi-res hierarchy thread code simpler.
  void HDF5_IO::writeVolumeFile(const Volume& wvol, 
				const std::string& filename,
				unsigned int var, unsigned int time,
				uint64 off_x, uint64 off_y, uint64 off_z) const
  {
    using namespace CVC::HDF5_Utilities;
    using namespace boost;

    CVC::ThreadInfo ti(BOOST_CURRENT_FUNCTION);

    std::string actualFileName;
    std::string objectName;

    boost::tie(actualFileName, objectName) =
      splitRawFilename(filename);

    std::string volume_name;
    bool doBuildHierarchy = false;
    if(isGroup(actualFileName,objectName))
      {
        //Name of actual dataset.  The group name in 'objectName' can contain several
        //instances of the same dataset at various resolutions.  This function assumes
        //you want the highest resolution dataset, which uses the following naming convention.
        volume_name =
          str(
            format("%1%/%2%:%3%:%4%") %
            objectName %
            DEFAULT_VOLUME_NAME %
            var % time
          );

        doBuildHierarchy = true;

        //If we have a thread running already computing the hierarchy, stop it!
        std::string threadKey(volume_name + " hierarchy_thread");
        if(cvcapp.hasThread(threadKey))
          {
            CVC::ThreadPtr t = cvcapp.threads(threadKey);
            t->interrupt(); //initiate thread quit
            t->join(); //wait for it to quit
          }

        //Now mark all of the hierarchy datasets dirty, creating them if they don't exist.
        //The BuildHierarchy thread fills these in with proper data later.
        Dimension prevDim(getObjectDimension(actualFileName,objectName));
        while(1)
          {
            uint64 maxdim = std::max(prevDim.xdim,std::max(prevDim.ydim,prevDim.zdim));
            maxdim = CVC::upToPowerOfTwo(maxdim) >> 1; //power of 2 less than maxdim
            Dimension curDim(maxdim,maxdim,maxdim);
            for(int i = 0; i < 3; i++)
              curDim[i] = std::min(curDim[i],prevDim[i]);
            if(curDim.size()==1) break; //we're done if the dims hit zero

            std::string hier_volume_name =
              str(format("%1%_%2%x%3%x%4%") 
                  % volume_name
                  % curDim[0] % curDim[1] % curDim[2]);

            cvcapp.log(10,str(format("%1% :: marking %2% dirty\n")
                              % BOOST_CURRENT_FUNCTION
                              % hier_volume_name));

            if(!isDataSet(actualFileName,hier_volume_name))
              createVolumeDataSet(actualFileName,hier_volume_name,
                                  getObjectBoundingBox(actualFileName,volume_name),
                                  curDim, wvol.voxelType());
            setAttribute(actualFileName, hier_volume_name, "dirty", 1);

            prevDim = curDim;
          }
      }
    else //Ungrouped dataset
      {
        if(var > 0 || time > 0)
          cvcapp.log(1,str(format("%s :: WARNING - HDF5_IO doesn't yet support lone datasets "
                                  "with more than one variable or timestep!\n")
                           % BOOST_CURRENT_FUNCTION));
        var = 0; time = 0;
        volume_name = objectName;
      }

    writeDataSet(actualFileName, volume_name, 
                 off_x, off_y, off_z,
                 wvol.dimension(),
                 wvol.voxelType(),
                 *wvol,
                 wvol.min(), wvol.max());
    
    setAttribute(actualFileName, volume_name, "info", wvol.desc());

    if(doBuildHierarchy)
      {
        std::string threadKey("build_hierarchy_" + volume_name);
        BuildHierarchy::start(threadKey,
                              actualFileName,
                              volume_name);
      }
  }

  // -------------------------------
  // VolumeFile_IO::writeBoundingBox
  // -------------------------------
  // Purpose:
  //   Writes the specified bounding box to the file.
  // ---- Change History ----
  // 04/06/2012 -- Joe R. -- Initial implementation.
  void HDF5_IO::writeBoundingBox(const BoundingBox& bbox, const std::string& filename) const
  {
    using namespace CVC::HDF5_Utilities;
    using namespace boost;

    CVC::ThreadInfo ti(BOOST_CURRENT_FUNCTION);

    std::string actualFileName;
    std::string objectName;

    boost::tie(actualFileName, objectName) =
      splitRawFilename(filename);
    
    if(isGroup(actualFileName,objectName))
      {
        std::vector<std::string> children = getChildObjects(actualFileName,
                                                            objectName);
        for (const auto& val : children)
          setObjectBoundingBox(actualFileName, objectName + "/" + val, bbox);
      }

    setObjectBoundingBox(actualFileName, objectName, bbox);
  }

  // ----------------------------
  // HDF5_IO::createVolumeDataSet
  // ----------------------------
  // Purpose:
  //   Creates a volume dataset without a group.  Used for building the multi-res
  //   hierarchy.
  // ---- Change History ----
  // 09/09/2011 -- Joe R. -- Initial implementation.
  // 09/30/2011 -- Joe R. -- Added objectType attribute.
  void HDF5_IO::createVolumeDataSet(const std::string& hdf5_filename,
                                    const std::string& volumeDataSet,
                                    const CVC::BoundingBox& boundingBox,
                                    const CVC::Dimension& dimension,
                                    VolMagick::VoxelType voxelType)
  {
    using namespace CVC::HDF5_Utilities;

    createDataSet(hdf5_filename, volumeDataSet,
                  boundingBox, dimension, voxelType, true);
    setAttribute(hdf5_filename, volumeDataSet, 
                 "objectType", "VolMagick::Volume");
    setAttribute(hdf5_filename, volumeDataSet, 
                 "libcvc_version", VOLMAGICK_VERSION_STRING);
    setObjectBoundingBox(hdf5_filename, volumeDataSet, boundingBox);
    setObjectDimension(hdf5_filename, volumeDataSet, dimension);
    setAttribute(hdf5_filename, volumeDataSet, "dataTypes", 1, &voxelType);
    setAttribute(hdf5_filename, volumeDataSet, "numVariables", 1);
    setAttribute(hdf5_filename, volumeDataSet, "numTimesteps", 1);
    setAttribute(hdf5_filename, volumeDataSet, "min_time", 0.0);
    setAttribute(hdf5_filename, volumeDataSet, "max_time", 0.0);
  }
}
