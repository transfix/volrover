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

#ifndef __HDF5_IO_H__
#define __HDF5_IO_H__

#include <VolMagick/VolumeFile_IO.h>
#include <VolMagick/Exceptions.h>
#include <VolMagick/Types.h>
#include <VolMagick/Dimension.h>

#include <cvc_compat.h>
#if defined (WIN32)
#include <cpp/H5Cpp.h>
#else 
#include <H5Cpp.h>
#endif

namespace VolMagick
{
  //VOLMAGICK_DEF_EXCEPTION(InvalidHDF5File);
  //VOLMAGICK_DEF_EXCEPTION(HDF5Exception);
  typedef CVC::InvalidHDF5File InvalidHDF5File;
  typedef CVC::HDF5Exception HDF5Exception;

  // -------
  // HDF5_IO
  // -------
  // Purpose: 
  //   Provides HDF5 file support.
  // ---- Change History ----
  // 12/04/2009 -- Joe R. -- Initial implementation.
  struct HDF5_IO : public VolumeFile_IO
  {
    /*
     * CVC hdf5 schema -- Joe R. -- 01/04/2010
     *  / - root
     *  |_ cvc/ - root cvc data hierarchy
     *     |_ geometry/ - placeholder, define fully later
     *     |_ transfer_functions/ - placeholder, define fully later
     *     |_ volumes/
     *        |_ <volume name> - volume group containing image or volume data
     *           |               Default volume name is 'volume'.  This group has the following attribs:
     *           |               - VolMagick_version (uint64)
     *           |               - XMin, YMin, ZMin,
     *           |                 XMax, YMax, ZMax  (double) - bounding box
     *           |               - XDim, YDim, ZDim (uint64) - volume dimensions
     *           |               - voxelTypes (uint64 array) - the type of each variable
     *           |               - numVariables (uint64)
     *           |               - numTimesteps (uint64)
     *           |               - min_time (double)
     *           |               - max_time (double)
     *           |_ <volume name>:<variable (int)>:<timestep (int)> - dataset for a volume.  each variable
     *                                                                of each timestep has it's own dataset.
     *                                                                Each volume dataset has the following
     *                                                                attributes:
     *                                                                - min (double) - min voxel value
     *                                                                - max (double) - max voxel value
     *                                                                - name (string) - variable name
     *                                                                - voxelType (uint64) - type of this dataset
     *  
     */
    static const hsize_t VOLUME_ATTRIBUTE_STRING_MAXLEN = 255;

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
    HDF5_IO();

    // -----------
    // HDF5_IO::id
    // -----------
    // Purpose:
    //   Returns a string that identifies this VolumeFile_IO object.  This should
    //   be unique, but is freeform.
    // ---- Change History ----
    // 12/04/2009 -- Joe R. -- Initial implementation.
    virtual const std::string& id() const;

    // -------------------
    // HDF5_IO::extensions
    // -------------------
    // Purpose:
    //   Returns a list of extensions that this VolumeFile_IO object supports.
    // ---- Change History ----
    // 12/04/2009 -- Joe R. -- Initial implementation.
    virtual const ExtensionList& extensions() const;

    // --------------------------
    // HDF5_IO::getVolumeFileInfo
    // --------------------------
    // Purpose:
    //   Writes to a structure containing all info that VolMagick needs
    //   from a volume file.
    // ---- Change History ----
    // 12/04/2009 -- Joe R. -- Initial implementation.
    virtual void getVolumeFileInfo(VolumeFileInfo::Data& data,
				   const std::string& filename) const;

    // -----------------------
    // HDF5_IO::readVolumeFile
    // -----------------------
    // Purpose:
    //   Writes to a Volume object after reading from a volume file.
    // ---- Change History ----
    // 12/04/2009 -- Joe R. -- Initial implementation.
    virtual void readVolumeFile(Volume& vol,
				const std::string& filename, 
				unsigned int var, unsigned int time,
				uint64 off_x, uint64 off_y, uint64 off_z,
				const Dimension& subvoldim) const;

    // -----------------------------
    // VolumeFile_IO::readVolumeFile
    // -----------------------------
    // Purpose:
    //   Same as above except uses a bounding box for specifying the
    //   subvol.  Uses maxdim to define a stride to use when reading
    //   for subsampling.
    // ---- Change History ----
    // 01/04/2010 -- Joe R. -- Initial implementation.
    virtual void readVolumeFile(Volume& vol, 
				const std::string& filename, 
				unsigned int var,
				unsigned int time,
				const BoundingBox& subvolbox) const;

    // -------------------------
    // HDF5_IO::createVolumeFile
    // -------------------------
    // Purpose:
    //   Creates an empty volume file to be later filled in by writeVolumeFile
    // ---- Change History ----
    // 12/04/2009 -- Joe R. -- Initial implementation.
    virtual void createVolumeFile(const std::string& filename,
				  const BoundingBox& boundingBox,
				  const Dimension& dimension,
				  const std::vector<VoxelType>& voxelTypes,
				  unsigned int numVariables, unsigned int numTimesteps,
				  double min_time, double max_time) const;

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
    virtual void writeVolumeFile(const Volume& wvol, 
				 const std::string& filename,
				 unsigned int var, unsigned int time,
				 uint64 off_x, uint64 off_y, uint64 off_z) const;

    // -------------------------------
    // VolumeFile_IO::writeBoundingBox
    // -------------------------------
    // Purpose:
    //   Writes the specified bounding box to the file.
    // ---- Change History ----
    // 04/06/2012 -- Joe R. -- Initial implementation.
    virtual void writeBoundingBox(const BoundingBox& bbox, const std::string& filename) const;

    // ----------------------------
    // HDF5_IO::createVolumeDataSet
    // ----------------------------
    // Purpose:
    //   Creates a volume dataset without a group.  Used for building the multi-res
    //   hierarchy.  TODO: this should probably be moved to HDF5_Utilities
    // ---- Change History ----
    // 09/09/2011 -- Joe R. -- Initial implementation.
    static void createVolumeDataSet(const std::string& hdf5_filename,
                                    const std::string& volumeDataSet,
                                    const CVC::BoundingBox& boundingBox,
                                    const CVC::Dimension& dimension,
                                    VolMagick::VoxelType voxelType);

  protected:
    std::string _id;
    ExtensionList _extensions;
  };
}

#endif
