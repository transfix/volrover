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

/* $Id: VolumeFile_IO.cpp 5355 2012-04-06 22:16:56Z transfix $ */

#include <VolMagick/VolumeFile_IO.h>
#include <VolMagick/Exceptions.h>

// **** VolMagick File I/O classes
#include <VolMagick/Null_IO.h>
#include <VolMagick/RawIV_IO.h>
#include <VolMagick/RawV_IO.h>
#include <VolMagick/MRC_IO.h>
#include <VolMagick/Spider_IO.h>
#include <VolMagick/VTK_IO.h>
#ifdef VOLMAGICK_USING_IMOD_MRC
#include <VolMagick/IMOD_MRC_IO.h>
#endif
#ifdef VOLMAGICK_USING_VOLMAGICK_INR
#include <VolMagick/INR_IO.h>
#endif
#ifdef VOLMAGICK_USING_HDF5
#include <VolMagick/HDF5_IO.h>
#endif
#ifdef VOLMAGICK_USING_QTIMAGE
#include <VolMagick/QtImage_IO.h>
#endif

#include <VolMagick/Utility.h>

#include <CVC/App.h>

#include <boost/foreach.hpp>

namespace VolMagick
{
  // 09/05/2011 -- Joe R. -- moved the following 2 const vars from HDF5_IO
  //The group where volumes are stored in files that support hierarchical data structures
  const char* VolumeFile_IO::CVC_VOLUME_GROUP = "/cvc/volumes";
  
  //The name to use when no volume name is specified
  const char* VolumeFile_IO::DEFAULT_VOLUME_NAME = "volume";

  //A regex to extract a filename extension
  const char* VolumeFile_IO::FILE_EXTENSION_EXPR = "^(.*)(\\.\\S*)$";

  // -------------------------------
  // VolumeFile_IO::splitRawFilename
  // -------------------------------
  // Purpose:
  //   Splits a filename into an actual file name and an object
  //   path.
  // ---- Change History ----
  // 07/24/2011 -- Joe R. -- Initial implementation.
  // 08/05/2011 -- Joe R. -- Using a '|' instead of ':'
  // 09/05/2011 -- Joe R. -- Moved here from HDF5_IO
  boost::tuple<
    std::string, /* actual file name */
    std::string  /* hdf5 object name */
  > VolumeFile_IO::splitRawFilename(const std::string& filename)
  {
    typedef std::vector<std::string> split_vector_type;
    split_vector_type names;
    boost::algorithm::split( names, filename, boost::algorithm::is_any_of("|") );
    names.resize(2);

    //Use default name if no object name was specified
    if(names[1].empty())
      names[1] = 
        std::string(CVC_VOLUME_GROUP) + "/" + 
        std::string(DEFAULT_VOLUME_NAME);

    return boost::make_tuple(names[0],names[1]);
  }

  // -----------------------------
  // VolumeFile_IO::readVolumeFile
  // -----------------------------
  // Purpose:
  //   Same as above except uses a bounding box for specifying the
  //   subvol.  A default implementation is provided. 
  // ---- Change History ----
  // 01/03/2010 -- Joe R. -- Initial implementation.
  // 05/11/2010 -- Joe R. -- Fixing off-by-one indexing problem.
  void VolumeFile_IO::readVolumeFile(Volume& vol, 
				     const std::string& filename, 
				     unsigned int var,
				     unsigned int time,
				     const BoundingBox& subvolbox) const
  {
    VolumeFileInfo volinfo(filename);
    if(!subvolbox.isWithin(volinfo.boundingBox()))
      throw SubVolumeOutOfBounds("The subvolume bounding box must be within the file's bounding box.");
    uint64 off_x = uint64((subvolbox.minx - volinfo.XMin())/volinfo.XSpan());
    uint64 off_y = uint64((subvolbox.miny - volinfo.YMin())/volinfo.YSpan());
    uint64 off_z = uint64((subvolbox.minz - volinfo.ZMin())/volinfo.ZSpan());
    Dimension dim;
    dim[0] = uint64((subvolbox.maxx - subvolbox.minx)/volinfo.XSpan())+1;
    dim[1] = uint64((subvolbox.maxy - subvolbox.miny)/volinfo.YSpan())+1;
    dim[2] = uint64((subvolbox.maxz - subvolbox.minz)/volinfo.ZSpan())+1;
    for(int i = 0; i < 3; i++) if(dim[i] == 0) dim[i]=1;
    if(dim[0] + off_x > volinfo.XDim()) dim[0] = volinfo.XDim() - off_x;
    if(dim[1] + off_y > volinfo.YDim()) dim[1] = volinfo.YDim() - off_y;
    if(dim[2] + off_z > volinfo.ZDim()) dim[2] = volinfo.ZDim() - off_z;
    readVolumeFile(vol,filename,var,time,off_x,off_y,off_z,dim);
    //vol.sub(subvolbox,dim); //get a subvolume that is exactly the size of subvolbox
    //just force the bounding box for now.. this might lead to aliasing errors
    vol.boundingBox(subvolbox);
  }
  
  // -------------------------------
  // VolumeFile_IO::writeBoundingBox
  // -------------------------------
  // Purpose:
  //   Writes the specified bounding box to the file.  The default implementation is slow
  //   because it has to read the entire file.  This can be sped up on an individual file type
  //   basis.
  // ---- Change History ----
  // 04/06/2012 -- Joe R. -- Initial implementation.
  void VolumeFile_IO::writeBoundingBox(const BoundingBox& bbox, const std::string& filename) const
  {
    std::vector<Volume> vols;
    VolumeFileInfo vfi(filename);
    vfi.boundingBox(bbox);
    VolMagick::readVolumeFile(vols,filename);
    for (auto& vol : vols)
      vol.boundingBox(bbox);
    VolMagick::createVolumeFile(filename,vfi); //TODO: don't overwrite existing file until temp file write is complete
    VolMagick::writeVolumeFile(vols,filename);
  }

  // --------------------------
  // VolumeFile_IO::handlersMap
  // --------------------------
  // Purpose:
  //   Static initialization of handler map.  Clients use the HandlerMap
  //   to add themselves to the collection of objects that are to be used
  //   to perform volume file i/o operations.
  // ---- Change History ----
  // 11/13/2009 -- Joe R. -- Initially implemented.
  VolumeFile_IO::HandlerMap& VolumeFile_IO::handlerMap()
  {
    //It's ok to leak: http://www.parashift.com/c++-faq-lite/ctors.html#faq-10.15
    static HandlerMap* p = initializeMap();
    return *p;
  }

  // ----------------------------
  // VolumeFile_IO::insertHandler
  // ----------------------------
  // Purpose:
  //   Convenence function for adding objects to the map.
  // ---- Change History ----
  // 11/13/2009 -- Joe R. -- Initially implemented.
  // 11/20/2009 -- Joe R. -- Removed extension arg.  Now using what the object provides.
  void VolumeFile_IO::insertHandler(const Ptr& vfio)
  {
    insertHandler(handlerMap(),vfio);
  }

  // ----------------------------
  // VolumeFile_IO::removeHandler
  // ----------------------------
  // Purpose:
  //   Convenence function for removing objects from the map.
  // ---- Change History ----
  // 11/13/2009 -- Joe R. -- Initially implemented.
  void VolumeFile_IO::removeHandler(const Ptr& vfio)
  {
    for(HandlerMap::iterator i = handlerMap().begin();
	i != handlerMap().end();
	i++)
      {
	Handlers handlers;
	for(Handlers::iterator j = i->second.begin();
	    j != i->second.end();
	    j++)
	  {
	    if(*j != vfio)
	      handlers.push_back(*j);
	  }
	i->second = handlers;
      }
  }

  // ----------------------------
  // VolumeFile_IO::removeHandler
  // ----------------------------
  // Purpose:
  //   Convenence function for removing objects from the map.
  // ---- Change History ----
  // 11/13/2009 -- Joe R. -- Initially implemented.
  void VolumeFile_IO::removeHandler(const std::string& id)
  {
    for(HandlerMap::iterator i = handlerMap().begin();
	i != handlerMap().end();
	i++)
      {
	Handlers handlers;
	for(Handlers::iterator j = i->second.begin();
	    j != i->second.end();
	    j++)
	  {
	    if((*j)->id() != id)
	      handlers.push_back(*j);
	  }
	i->second = handlers;
      }    
  }

  // ----------------------------
  // VolumeFile_IO::getExtensions
  // ----------------------------
  // Purpose:
  //   Returns the list of supported file extensions.
  // ---- Change History ----
  // 09/18/2011 -- Joe R. -- Initial implementation.    
  std::vector<std::string> VolumeFile_IO::getExtensions()
  {
    std::vector<std::string> ret;
    for (auto& i : handlerMap()) {
        ret.push_back(i.first);
      }
    return ret;
  }

  // ----------------------------
  // VolumeFile_IO::initializeMap
  // ----------------------------
  // Purpose:
  //   Adds the standard VolumeFile_IO objects to a new HandlerMap object
  // ---- Change History ----
  // 11/20/2009 -- Joe R. -- Initially implemented.
  VolumeFile_IO::HandlerMap *VolumeFile_IO::initializeMap()
  {
    HandlerMap *map = new HandlerMap;
    HandlerMap &ref = *map;

    insertHandler(ref, Ptr(new Null_IO));
    insertHandler(ref, Ptr(new RawIV_IO));
    insertHandler(ref, Ptr(new RawIV_IO(true)));
    insertHandler(ref, Ptr(new RawV_IO));
    insertHandler(ref, Ptr(new MRC_IO));
    insertHandler(ref, Ptr(new Spider_IO));
    insertHandler(ref, Ptr(new VTK_IO));

#ifdef VOLMAGICK_USING_IMOD_MRC
    insertHandler(ref, Ptr(new IMOD_MRC_IO));
#endif

#ifdef VOLMAGICK_USING_VOLMAGICK_INR
    insertHandler(ref, Ptr(new INR_IO));
#endif

#ifdef VOLMAGICK_USING_HDF5
    insertHandler(ref, Ptr(new HDF5_IO));
#endif

#ifdef VOLMAGICK_USING_QTIMAGE
    insertHandler(ref, Ptr(new QtImage_IO));
#endif
    
    return map;
  }

  // ----------------------------
  // VolumeFile_IO::insertHandler
  // ----------------------------
  // Purpose:
  //   Convenence function for adding objects to the specified map.
  // ---- Change History ----
  // 11/20/2009 -- Joe R. -- Initially implemented.
  void VolumeFile_IO::insertHandler(HandlerMap& hm,
				   const Ptr& vfio)
  {
    if(!vfio) return;
    for(VolumeFile_IO::ExtensionList::const_iterator i =   
	  vfio->extensions().begin();                      
	i != vfio->extensions().end();                     
	i++)
      {
#if 0
	std::cerr<<BOOST_CURRENT_FUNCTION<<": inserting handler '" 
		 << vfio->id() << "' for extension '" << *i << "'."<<std::endl;
#endif
	hm[*i].push_back(vfio);
      }
  }

  void readVolumeFile(Volume& vol, 
		      const std::string& filename,
		      unsigned int var, unsigned int time)
  {
    VolumeFileInfo volinfo(filename);
    readVolumeFile(vol,filename,var,time,0,0,0,volinfo.dimension());
  }

  // --------------
  // readVolumeFile
  // --------------
  // Purpose:
  //   The main readVolumeFile function.  Refers to the handler map to choose
  //   an appropriate IO object for reading the requested volume file.
  // ---- Change History ----
  // ??/??/2007 -- Joe R. -- Initially implemented.
  // 11/13/2009 -- Joe R. -- Re-implemented using VolumeFile_IO handler map.
  // 12/28/2009 -- Joe R. -- Collecting exception error strings
  // 09/05/2011 -- Joe R. -- Using splitRawFilename to extract real filename
  //                         if the provided filename is a file|obj tuple.
  void readVolumeFile(Volume& vol,
		      const std::string& filename, 
		      unsigned int var, unsigned int time,
		      uint64 off_x, uint64 off_y, uint64 off_z,
		      const Dimension& subvoldim)
  {
    vol.unsetMinMax();

    std::string errors;
    boost::smatch what;

    std::string actualFileName;
    std::string objectName;

    boost::tie(actualFileName, objectName) =
      VolumeFile_IO::splitRawFilename(filename);

    const boost::regex file_extension(VolumeFile_IO::FILE_EXTENSION_EXPR);
    if(boost::regex_match(actualFileName, what, file_extension))
      {
	if(VolumeFile_IO::handlerMap()[what[2]].empty())
	  throw UnsupportedVolumeFileType(std::string(BOOST_CURRENT_FUNCTION) + 
					  std::string(": Cannot read ") + filename);
	VolumeFile_IO::Handlers& handlers = VolumeFile_IO::handlerMap()[what[2]];
	//use the first handler that succeds
	for(VolumeFile_IO::Handlers::iterator i = handlers.begin();
	    i != handlers.end();
	    i++)
	  try
	    {
	      if(*i)
		{
		  (*i)->readVolumeFile(vol,filename,var,time,
				       off_x,off_y,off_z,subvoldim);
		  return;
		}
	    }
	  catch(VolMagick::Exception& e)
	    {
	      errors += std::string(" :: ") + e.what();
	    }
      }
    throw UnsupportedVolumeFileType(
      boost::str(
	boost::format("%1% : Cannot read '%2%'%3%") % 
	BOOST_CURRENT_FUNCTION %
	filename %
	errors
      )
    );
  }

  // --------------
  // readVolumeFile
  // --------------
  // Purpose:
  //    Same as above except it uses a bounding box.
  //   
  // ---- Change History ----
  // ??/??/2009 -- Joe R. -- Initially implemented.
  // 01/03/2010 -- Joe R. -- Re-implemented using VolumeFile_IO handler map.
  // 09/05/2011 -- Joe R. -- Using splitRawFilename to extract real filename
  //                         if the provided filename is a file|obj tuple.
  void readVolumeFile(Volume& vol,
		      const std::string& filename, 
		      unsigned int var, unsigned int time,
		      const BoundingBox& subvolbox)
  {
    vol.unsetMinMax();

    std::string errors;
    boost::smatch what;

    std::string actualFileName;
    std::string objectName;

    boost::tie(actualFileName, objectName) =
      VolumeFile_IO::splitRawFilename(filename);

    const boost::regex file_extension(VolumeFile_IO::FILE_EXTENSION_EXPR);
    if(boost::regex_match(actualFileName, what, file_extension))
      {
	if(VolumeFile_IO::handlerMap()[what[2]].empty())
	  throw UnsupportedVolumeFileType(std::string(BOOST_CURRENT_FUNCTION) + 
					  std::string(": Cannot read ") + filename);
	VolumeFile_IO::Handlers& handlers = VolumeFile_IO::handlerMap()[what[2]];
	//use the first handler that succeds
	for(VolumeFile_IO::Handlers::iterator i = handlers.begin();
	    i != handlers.end();
	    i++)
	  try
	    {
	      if(*i)
		{
		  (*i)->readVolumeFile(vol,filename,var,time,
				       subvolbox);
		  return;
		}
	    }
	  catch(VolMagick::Exception& e)
	    {
	      errors += std::string(" :::: ") + e.what();
	    }
      }
    throw UnsupportedVolumeFileType(
      boost::str(
	boost::format("%1% : Cannot read '%2%'%3%") % 
	BOOST_CURRENT_FUNCTION %
	filename %
	errors
      )
    );
  }

  void readVolumeFile(std::vector<Volume>& vols,
		      const std::string& filename)
  {
    VolumeFileInfo volinfo(filename);
    Volume vol;
    vols.clear();
    for(unsigned int var=0; var<volinfo.numVariables(); var++)
      for(unsigned int time=0; time<volinfo.numTimesteps(); time++)
	{
	  readVolumeFile(vol,filename,var,time);
	  vols.push_back(vol);
	}
  }

  // ---------------
  // writeVolumeFile
  // ---------------
  // Purpose:
  //   The main writeVolumeFile function.  Refers to the handler map to choose
  //   an appropriate IO object for writing to the requested volume file.
  // ---- Change History ----
  // ??/??/2007 -- Joe R. -- Initially implemented.
  // 11/13/2009 -- Joe R. -- Re-implemented using VolumeFile_IO handler map.
  // 12/28/2009 -- Joe R. -- Collecting exception error strings
  // 09/05/2011 -- Joe R. -- Using splitRawFilename to extract real filename
  //                         if the provided filename is a file|obj tuple.
  void writeVolumeFile(const Volume& vol, 
		       const std::string& filename,
		       unsigned int var, unsigned int time,
		       uint64 off_x, uint64 off_y, uint64 off_z)
  {
    std::string errors;
    boost::smatch what;

    std::string actualFileName;
    std::string objectName;

    boost::tie(actualFileName, objectName) =
      VolumeFile_IO::splitRawFilename(filename);

    const boost::regex file_extension(VolumeFile_IO::FILE_EXTENSION_EXPR);
    if(boost::regex_match(actualFileName, what, file_extension))
      {
	if(VolumeFile_IO::handlerMap()[what[2]].empty())
	  throw UnsupportedVolumeFileType(std::string(BOOST_CURRENT_FUNCTION) + 
					  std::string(": Cannot write ") + filename);
	VolumeFile_IO::Handlers& handlers = VolumeFile_IO::handlerMap()[what[2]];
	//use the first handler that succeds
	for(VolumeFile_IO::Handlers::iterator i = handlers.begin();
	    i != handlers.end();
	    i++)
	  try
	    {
	      if(*i)
		{
		  (*i)->writeVolumeFile(vol,filename,var,time,off_x,off_y,off_z);
		  return;
		}
	    }
	  catch(VolMagick::Exception& e)
	    {
	      errors += std::string(" :: ") + e.what();
	    }
      }
    throw UnsupportedVolumeFileType(
      boost::str(
	boost::format("%1% : Cannot read '%2%'%3%") % 
	BOOST_CURRENT_FUNCTION %
	filename %
	errors
      )
    );
  }

  // ---------------
  // writeVolumeFile
  // ---------------
  // Purpose:
  //   Writes a volume to the specified filename, using the specified bounding box
  //   as the target to fill in the write volume.  Interpolates the input volume
  //   appropriately.
  // ---- Change History ----
  // ??/??/2010 -- Joe R. -- I think this came about sometime in 2010.  It's a nasty
  //                         and slow function but I think it does the trick.
  // 09/10/2011 -- Joe R. -- Adding thread progress feedback.
  // 09/11/2011 -- Joe R. -- Fixing an indexing bug.
  void writeVolumeFile(const Volume& vol, 
		       const std::string& filename,
		       unsigned int var, unsigned int time,
		       const BoundingBox& subvolbox)
  {
    CVC::ThreadInfo ti(BOOST_CURRENT_FUNCTION);

    Volume localvol(vol);
    VolumeFileInfo volinfo(filename);
    if(!subvolbox.isWithin(volinfo.boundingBox()))
      throw SubVolumeOutOfBounds("The subvolume bounding box must be within the file's bounding box.");

    //set up a subvolume to represent the chunk we are writing to in the file
    uint64 off_x = uint64((subvolbox.minx - volinfo.XMin())/volinfo.XSpan());
    uint64 off_y = uint64((subvolbox.miny - volinfo.YMin())/volinfo.YSpan());
    uint64 off_z = uint64((subvolbox.minz - volinfo.ZMin())/volinfo.ZSpan());
    Dimension dim;
    dim[0] = uint64((subvolbox.maxx - subvolbox.minx)/volinfo.XSpan()+1);
    dim[1] = uint64((subvolbox.maxy - subvolbox.miny)/volinfo.YSpan()+1);
    dim[2] = uint64((subvolbox.maxz - subvolbox.minz)/volinfo.ZSpan()+1);
    if(dim[0] + off_x > volinfo.XDim()) dim[0] = volinfo.XDim() <= 1 ? 1 : volinfo.XDim() - off_x;
    if(dim[1] + off_y > volinfo.YDim()) dim[1] = volinfo.YDim() <= 1 ? 1 : volinfo.YDim() - off_y;
    if(dim[2] + off_z > volinfo.ZDim()) dim[2] = volinfo.ZDim() <= 1 ? 1 : volinfo.ZDim() - off_z;

    Volume subvol(dim,volinfo.voxelTypes(var),
		  VolMagick::BoundingBox(volinfo.XMin() + off_x*volinfo.XSpan(),
					 volinfo.YMin() + off_y*volinfo.YSpan(),
					 volinfo.ZMin() + off_z*volinfo.ZSpan(),
					 volinfo.XMin() + (off_x+dim[0])*volinfo.XSpan(),
					 volinfo.YMin() + (off_y+dim[1])*volinfo.YSpan(),
					 volinfo.ZMin() + (off_z+dim[2])*volinfo.ZSpan()));
    
    //reset the bounding box for easy interpolation
    localvol.boundingBox(VolMagick::BoundingBox(0.0,0.0,0.0,
						1.0,1.0,1.0));

    double xmax = subvol.XDim() > 1 ? double(subvol.XDim()-1) : 1.0;
    double ymax = subvol.YDim() > 1 ? double(subvol.YDim()-1) : 1.0;
    double zmax = subvol.ZDim() > 1 ? double(subvol.ZDim()-1) : 1.0;
    for(VolMagick::uint64 k = 0; k < subvol.ZDim(); k++)
      {
        for(VolMagick::uint64 j = 0; j < subvol.YDim(); j++)
          for(VolMagick::uint64 i = 0; i < subvol.XDim(); i++)
            {
              subvol(i,j,k, localvol.interpolate(double(i)/xmax,
                                                 double(j)/ymax,
                                                 double(k)/zmax));
            }
        cvcapp.threadProgress(float(k)/float(subvol.ZDim()));
      }
    
    cvcapp.threadProgress(1.0);
    writeVolumeFile(subvol,filename,var,time,off_x,off_y,off_z);
  }

  // ---------------
  // writeVolumeFile
  // ---------------
  // Purpose:
  //   Writes a collection of volumes to the specified file if the file type supports
  //   such an operation.
  // ---- Change History ----
  // ??/??/2007 -- Joe R. -- Initial implementation.
  void writeVolumeFile(const std::vector<Volume>& vols,
		       const std::string& filename)
  {
    uint64 i;

    if(vols.size() == 0) return;
    
    //create types vector
    std::vector<VoxelType> voxelTypes;
    for(i=0; i<vols.size(); i++) voxelTypes.push_back(vols[i].voxelType());

    //create the file and write the volume info
    createVolumeFile(filename,vols[0].boundingBox(),vols[0].dimension(),voxelTypes,vols.size());
    for(i=0; i<vols.size(); i++)
      writeVolumeFile(vols[i],filename,i);
  }

  // ----------------
  // createVolumeFile
  // ----------------
  // Purpose:
  //   The main createVolumeFile function.  Refers to the handler map to choose
  //   an appropriate IO object for creating the requested volume file.
  // ---- Change History ----
  // ??/??/2007 -- Joe R. -- Initially implemented.
  // 11/13/2009 -- Joe R. -- Re-implemented using VolumeFile_IO handler map.
  // 12/28/2009 -- Joe R. -- Collecting exception error strings
  // 09/05/2011 -- Joe R. -- Using splitRawFilename to extract real filename
  //                         if the provided filename is a file|obj tuple.
  void createVolumeFile(const std::string& filename,
			const BoundingBox& boundingBox,
			const Dimension& dimension,
			const std::vector<VoxelType>& voxelTypes,
			unsigned int numVariables, unsigned int numTimesteps,
			double min_time, double max_time)
  {
    std::string errors;
    boost::smatch what;

    std::string actualFileName;
    std::string objectName;

    boost::tie(actualFileName, objectName) =
      VolumeFile_IO::splitRawFilename(filename);

    const boost::regex file_extension(VolumeFile_IO::FILE_EXTENSION_EXPR);
    if(boost::regex_match(actualFileName, what, file_extension))
      {
	if(VolumeFile_IO::handlerMap()[what[2]].empty())
	  throw UnsupportedVolumeFileType(std::string(BOOST_CURRENT_FUNCTION) + 
					  std::string(": Cannot create ") + filename);
	VolumeFile_IO::Handlers& handlers = VolumeFile_IO::handlerMap()[what[2]];
	//use the first handler that succeds
	for(VolumeFile_IO::Handlers::iterator i = handlers.begin();
	    i != handlers.end();
	    i++)
	  try
	    {
	      if(*i)
		{
		  (*i)->createVolumeFile(filename,boundingBox,dimension,
					 voxelTypes,numVariables,numTimesteps,
					 min_time,max_time);
		  return;
		}
	    }
	  catch(VolMagick::Exception& e)
	    {
	      errors += std::string(" :: ") + e.what();
	    }
      }
    throw UnsupportedVolumeFileType(
      boost::str(
	boost::format("%1% : Cannot read '%2%'%3%") % 
	BOOST_CURRENT_FUNCTION %
	filename %
	errors
      )
    );
  }

  // ---------------
  // readBoundingBox
  // ---------------
  // Purpose:
  //  Returns a volume file's bounding box.
  // ---- Change History ----
  // 04/06/2012 -- Joe R. -- Initially implemented.
  BoundingBox readBoundingBox(const std::string& filename)
  {
    return VolumeFileInfo(filename).boundingBox();
  }

  // ----------------
  // writeBoundingBox
  // ----------------
  // Purpose:
  //  Changes a volume file's bounding box.
  // ---- Change History ----
  // 04/06/2012 -- Joe R. -- Initially implemented.
  void writeBoundingBox(const BoundingBox& bbox, const std::string& filename)
                        
  {
    std::string errors;
    boost::smatch what;

    std::string actualFileName;
    std::string objectName;

    boost::tie(actualFileName, objectName) =
      VolumeFile_IO::splitRawFilename(filename);

    const boost::regex file_extension(VolumeFile_IO::FILE_EXTENSION_EXPR);
    if(boost::regex_match(actualFileName, what, file_extension))
      {
	if(VolumeFile_IO::handlerMap()[what[2]].empty())
	  throw UnsupportedVolumeFileType(std::string(BOOST_CURRENT_FUNCTION) + 
					  std::string(": Cannot write ") + filename);
	VolumeFile_IO::Handlers& handlers = VolumeFile_IO::handlerMap()[what[2]];
	//use the first handler that succeds
	for(VolumeFile_IO::Handlers::iterator i = handlers.begin();
	    i != handlers.end();
	    i++)
	  try
	    {
	      if(*i)
		{
		  (*i)->writeBoundingBox(bbox,filename);
		  return;
		}
	    }
	  catch(VolMagick::Exception& e)
	    {
	      errors += std::string(" :: ") + e.what();
	    }
      }
    throw UnsupportedVolumeFileType(
      boost::str(
	boost::format("%1% : Cannot write '%2%'%3%") % 
	BOOST_CURRENT_FUNCTION %
	filename %
	errors
      )
    );
  }

}
