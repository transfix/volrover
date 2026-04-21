/*
  Copyright 2008-2011 The University of Texas at Austin

        Authors: Joe Rivera <transfix@ices.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of VolumeRover.

  VolumeRover is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  VolumeRover is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

/* $Id: main.cpp 5888 2012-07-22 22:57:49Z deukhyun $ */

#include <VolumeRover2/CVCMainWindow.h>
#include <VolumeRover2/VolumeInterface.h>
#include <VolumeRover2/VolumeMemoryInterface.h>
#include <VolumeRover2/GeometryInterface.h>
#include <VolumeRover2/Viewers.h>
#include <VolumeRover2/VolumeViewer2.h>
#include <VolumeRover2/Application.h>
#include <VolumeRover2/teebuf.h>
#include <VolumeRover2/windowbuf.h>

#include <QApplication>
#include <QFileInfo>
#include <QDir>
#include <QDateTime>
#include <QRegularExpression>

#include <VolMagick/VolMagick.h>
#include <cvcraw_geometry/cvcgeom.h>
#include <cvcraw_geometry/io.h>
#include <cvcraw_geometry/cvcraw_geometry.h>

#ifdef USING_TILING
#include <ContourTiler/reader_ser.h>
#include <ContourTiler/Contour.h>
#include <ContourTiler/Slice.h>
#include <cvcraw_geometry/contours.h>
#include <VolumeRover2/ContoursInterface.h>
#endif

#include <CVC/App.h>
#include <CVC/BoundingBox.h>
#include <CVC/Dimension.h>

#ifndef CVC_HDF5_DISABLED
#include <CVC/HDF5_Utilities.h>
#endif

#include <CVC/log4cplus_compat.h>

#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <boost/current_function.hpp>
#include <boost/foreach.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/scoped_ptr.hpp>

#include <iostream>
#include <cstdlib>

namespace
{
  bool load_volumeFileInfo(const std::string& filepath);


  // -------------
  // BatchThread
  // -------------
  // Purpose:
  //   Running predetermined commands/image capture
  //   from the command line
  // ---- Change History ----
  // 02/15/2012 -- arand - initial implementation
  class BatchThread
  {
  public:
    
    BatchThread(std::string vpmFile, std::string pngFile) 
      : _vpmFile(vpmFile), _pngFile(pngFile) {}

    //BatchThread(const BatchThread& t) {}

    //BatchThread& operator=(const BatchThread& t)
    // {
    //  return *this;
    //}

    void operator()()
    {

      CVC::ThreadFeedback feedback(BOOST_CURRENT_FUNCTION);


#ifdef WIN32
      Sleep(7);
#else
      sleep(7);
#endif

      std::cout << _vpmFile << " " << _pngFile << std::endl;
      cvcapp.readPropertyMap(_vpmFile);

      if (cvcapp.hasProperty("file_0_fullPath")) {
	cvcapp.readData(cvcapp.properties("file_0_fullPath"));
      }

      if (cvcapp.hasProperty("file_1_fullPath")) {
	cvcapp.readData(cvcapp.properties("file_1_fullPath"));
      }
      
      //if (cvcapp.hasProperty("file_2_fullPath")) {
      //cvcapp.readData(cvcapp.properties("file_2_fullPath"));
      //}

#ifdef WIN32
      Sleep(7);
#else
      sleep(7);
#endif
      cvcapp.readPropertyMap(_vpmFile);


#ifdef WIN32
      Sleep(3);
#else
      sleep(3);
#endif

      cvcapp.properties("thumbnail.snapshot_filename",_pngFile);

#ifdef WIN32
      Sleep(1);
#else
      sleep(1);
#endif

      cvcapp.properties("thumbnail.take_snapshot","true");

#ifdef WIN32
      Sleep(1);
#else
      sleep(1);
#endif
      
      cvcapp.properties("exitVolumeRover","true");
    }

  private:
    std::string _vpmFile;
    std::string _pngFile;
    
  };




  // -------------
  // CachingThread
  // -------------
  // Purpose:
  //   Copying to the cache a volume being read.
  // ---- Change History ----
  // 09/18/2011 -- Joe R. -- Initial implementation.
  class CachingThread
  {
  public:
    
    CachingThread(const std::string threadKey,
                  const std::string filepath)
      : _threadKey(threadKey),
        _filepath(filepath) {}

    CachingThread(const CachingThread& t)
      : _threadKey(t._threadKey),
        _filepath(t._filepath) {}

    CachingThread& operator=(const CachingThread& t)
    {
      _threadKey = t._threadKey;
      _filepath  = t._filepath;
      return *this;
    }

    static void start(const std::string threadKey,
                      const std::string filepath)
    {
      std::string uniqueThreadKey = threadKey;

      //For this thread, instead of starting a unique parallel thread,
      //we want to make sure this is the only one running.  So lets
      //try to stop the existing running thread with this key and wait for
      //it to actually stop before starting a new one.
      if(cvcapp.hasThread(threadKey))
        {
          CVC::ThreadPtr t = cvcapp.threads(threadKey);
          t->interrupt(); //initiate thread quit
          t->join(); //wait for it to quit
        }

      cvcapp.threads(
        uniqueThreadKey,
        CVC::ThreadPtr(
          new boost::thread(
            CachingThread(
               uniqueThreadKey,
               filepath
            )
          )
        )
      );
    }

    void operator()()
    {
#ifndef CVC_HDF5_DISABLED
      using namespace std;
      using namespace boost;
      using namespace CVC::HDF5_Utilities;

      CVC::ThreadFeedback feedback(BOOST_CURRENT_FUNCTION);

      cvcapp.log(2,str(format("%s :: caching %s\n")
                       % BOOST_CURRENT_FUNCTION
                       % _filepath));

      string cache_filename = cvcapp.properties("system.cache_file");
      vector<string> volumeobjs;
      try
        {
          volumeobjs = getChildObjects(cache_filename,"/cvc/volumes");
        }
      catch(std::exception&) {}
      string filename_to_open;
      string cache_object;
      bool copyToCache = true;
      for (auto& val : volumeobjs) {
          try
            {
              string origin;
              getAttribute(cache_filename,"/cvc/volumes/"+val,"origin",origin);
              CVC::uint64 modtime;
              getAttribute(cache_filename,"/cvc/volumes/"+val,"lastchange",modtime);

              QFileInfo fileinfo(QString::fromStdString(origin));
              CVC::uint64 new_modtime = fileinfo.lastModified().toSecsSinceEpoch();

              if(_filepath == origin && new_modtime <= modtime)
                {
                  filename_to_open = cache_filename + "|";
                  cache_object = "/cvc/volumes/" + val;
                  filename_to_open += cache_object;
                  copyToCache = false;
                  break;
                }
              else
                copyToCache = true;
            }
          catch(std::exception&)
            {
              copyToCache = true;
            }
        }

      try
        {
          if(copyToCache)
            {
              //extract just the filename from the path to use as the key
              QFileInfo fileinfo(QString::fromStdString(_filepath));
              std::string filename = fileinfo.baseName().toStdString();

              std::string cache_object = "/cvc/volumes/" + filename;

              //make sure the key is unique
              int i = 0;
              while(objectExists(cache_filename,cache_object))
                cache_object = "/cvc/volumes/" + filename + "_" + 
                  boost::lexical_cast<std::string>(i++);
           
              filename_to_open = cache_filename + "|" + cache_object;

              //convert to cache file hdf5 object.  This step will perform the hierarchy build.
              VolMagick::volconvert(cvcapp, _filepath,filename_to_open);
              setAttribute(cache_filename,cache_object,"origin",_filepath);
              CVC::uint64 modtime = 
                QFileInfo(QString::fromStdString(_filepath)).lastModified().toSecsSinceEpoch();
              setAttribute(cache_filename,cache_object,"lastchange",modtime);
              load_volumeFileInfo(filename_to_open);
            }
          else
            {
              load_volumeFileInfo(filename_to_open);
            }
        }
      catch(std::exception& e)
        {
          cvcapp.log(1,
            boost::str(
	      boost::format("%s: std::exception: %s\n")
	      % BOOST_CURRENT_FUNCTION
              % e.what()));
        }
#else
      load_volumeFileInfo(_filepath); //if no HDF5, just load directly w/o caching.
#endif
    }

  private:
    std::string _threadKey;
    std::string _filepath;
  };

  // 09/23/2011 - transfix - fixed loading so that this function doesn't try to cache geometry!
  bool load_volumeFileInfo(const std::string& filepath)
  {
    static log4cplus::Logger logger = FUNCTION_LOGGER;

    try
      {
        using namespace boost;

        std::string actualPath;
        std::string objectName;
        boost::tie(actualPath,
                   objectName)
          = VolMagick::VolumeFile_IO::splitRawFilename(filepath);

        QFileInfo fileInfo;
        fileInfo.setFile(QString::fromStdString(actualPath));
        std::string ext = fileInfo.suffix().toStdString();
        ext = "." + ext;

        //print out the list of supported extensions
        std::vector<std::string> extensions = 
          VolMagick::VolumeFile_IO::getExtensions();
        std::vector<std::string> filtered_extensions;
        bool foundExt = false;
        for (auto& val : extensions)
          {
            LOG4CPLUS_TRACE(logger, "supported extension " << val);
            // cvcapp.log(3,str(format("%s :: supported extension %s\n") 
            //                  % BOOST_CURRENT_FUNCTION
            //                  % val));
            if(val != ".cvc") //load cvc files outright
              filtered_extensions.push_back(val);

            //determine whether VolMagick recognizes this file extension
            if(val == ext)
              foundExt = true;
          }

        //if VolMagick doesn't recognize it, just return without trying to load it
        if(!foundExt)
          return false;

        //determine whether we need to launch the caching thread
        bool launchCachingThread = false;

        if(cvcapp.properties("system.caching_enabled") == "true")
          for (auto& val : filtered_extensions)
            {
              if(val == ext)
                {
                  launchCachingThread = true;
                  break;
                }
            }
        
        if(launchCachingThread)
          CachingThread::start(str(format("opening %s") % filepath),
                               filepath);
        else
          {
            cvcapp.log(2,str(format("%s :: opening %s directly\n") 
                             % BOOST_CURRENT_FUNCTION
                             % filepath));

            //extract just the filename from the path to use as the key
            std::string filename = fileInfo.completeBaseName().toStdString();

            //make sure the key is unique
            CVC::DataMap map = cvcapp.data();
            std::string key = filename;
            int i = 0;
            while(map.find(key)!=map.end())
              key = filename + "_" + boost::lexical_cast<std::string>(i++);
            
            //Try loading the VFI and inserting it into the datamap
            cvcapp.data(key,
                        VolMagick::VolumeFileInfo(filepath));

            //Set the viewers property for rendering the newly loaded VFI
            cvcapp.properties("viewers.vfi",key);
          }
      }
    catch(std::exception &e)
      {
        cvcapp.log(1,
          boost::str(
	    boost::format("%s: std::exception: %s\n")
	    % BOOST_CURRENT_FUNCTION
            % e.what()));
        
        return false;
      }

    return true;
  }

  bool load_cvcgeom(const std::string& filepath)
  {
    using namespace boost;

    try
      {
	QFileInfo fileInfo(QString::fromStdString(filepath));
        //extract just the filename from the path to use as the key
        std::string suffix = fileInfo.suffix().toStdString(); 
	if (!boost::iequals(suffix, "raw") &&
	    !boost::iequals(suffix, "rawc") &&
	    !boost::iequals(suffix, "rawn") &&
	    !boost::iequals(suffix, "rawnc") &&
	    !boost::iequals(suffix, "off") )
	  return false;

        //extract just the filename from the path to use as the key
        std::string filename = 
          QFileInfo(QString::fromStdString(filepath))
          .baseName()
          .toStdString();

        //make sure the key is unique
        CVC::DataMap map = cvcapp.data();
        std::string key = filename;
        int i = 0;
        while(map.find(key)!=map.end())
          key = filename + "_" + boost::lexical_cast<std::string>(i++);

	using namespace std;
        cvcapp.log(2,str(format("%s :: load start\n") 
                         % BOOST_CURRENT_FUNCTION));
        cvcraw_geometry::cvcgeom_t cvcgeom(filepath);
        cvcapp.log(2,str(format("%s :: load done\n") 
                         % BOOST_CURRENT_FUNCTION));
        //copy it to the data map
        cvcapp.data(key, cvcgeom);
	
	// arand, 5-3-2011: updated with listPropertyAppend
	cvcapp.listPropertyAppend("zoomed.geometries",key);
	cvcapp.listPropertyAppend("thumbnail.geometries",key);
      }
    catch(std::runtime_error&)
      {
        return false;
      }
    std::cout << "Leave load function" << std::endl;
    return true;
  }

  bool load_image(const std::string& filepath)
  {
    using namespace boost;

    try
      {
	QFileInfo fileInfo(QString::fromStdString(filepath));
        //extract just the filename from the path to use as the key
        std::string suffix = fileInfo.suffix().toStdString(); 
	if (!boost::iequals(suffix, "png") &&
	    !boost::iequals(suffix, "jpg") &&
	    !boost::iequals(suffix, "jpeg") &&
	    !boost::iequals(suffix, "bmp") &&
	    !boost::iequals(suffix, "xbm") &&
	    !boost::iequals(suffix, "xpm") &&
	    !boost::iequals(suffix, "pnm") &&
	    !boost::iequals(suffix, "mng") &&
	    !boost::iequals(suffix, "gif") )
	{
	  std::cout << "in load_image(no matching file format): " << filepath << std::endl;
	  return false;
	}

        //extract just the filename from the path to use as the key
        std::string filename = 
          QFileInfo(QString::fromStdString(filepath))
          .baseName()
          .toStdString();

	// generate vrtmp directory and put temp img files in there.
	// this part must be fixed to use volume in memory without storing the file
	if(!std::system("mkdir vrtmp"))
        cvcapp.log(2,str(format("%s :: failing to generate \n") 
                         % BOOST_CURRENT_FUNCTION));

        std::string obase = "vrtmp/tmpimg";
        std::string ofile = obase;

        //make sure the key is unique
        CVC::DataMap map = cvcapp.data();
        std::string key = filename;
        int i = 0;
        while(map.find(key)!=map.end()) {
          key = filename + "_" + boost::lexical_cast<std::string>(i++);
          ofile = obase + "_" + boost::lexical_cast<std::string>(i++);
        }
	ofile += ".rawiv";

	using namespace std;
        cvcapp.log(2,str(format("%s :: load start\n") 
                         % BOOST_CURRENT_FUNCTION));
        // load image
        QImage img(QString::fromStdString(filepath));
        if( img.isNull() ) {
           cvcapp.log(2,str(format("%s :: image loading failure!\n") 
                            % BOOST_CURRENT_FUNCTION));
           return false;
        }

	// convert the image to volume in (w,h,1) resolution
        int w = img.width(), h = img.height();
 
	VolMagick::Volume vol(VolMagick::Dimension(w,h,1),
                              VolMagick::UChar,
                              VolMagick::BoundingBox(0.0,0.0,0.0,w-1.0,h-1.0,1.0));

	for(int i = 0; i < w; i++)
	   for(int j = 0; j < h; j++)
              vol(i,j,0, qGray(img.pixel(i,j)));

	// write volume for test
	VolMagick::createVolumeFile(cvcapp, ofile.c_str(),
				  vol.boundingBox(),
				  vol.dimension(),
				  std::vector<VolMagick::VoxelType>(1, vol.voxelType()));

	VolMagick::writeVolumeFile(cvcapp, vol,ofile);

        cvcapp.log(2,str(format("%s :: load done\n") 
                         % BOOST_CURRENT_FUNCTION));

	VolMagick::VolumeFileInfo vfi(cvcapp, ofile);

        //copy it to the data map
        cvcapp.data(key, vfi);

        //Set the viewers property for rendering the newly loaded VFI
        cvcapp.properties("viewers.vfi",key);
      }
    catch(std::runtime_error&)
      {
        return false;
      }
    std::cout << "Leave load function" << std::endl;
    return true;
  }

  bool load_2DMRC(const std::string& filepath)
  {
    using namespace boost;

    try
      {
	QFileInfo fileInfo(QString::fromStdString(filepath));
        //extract just the filename from the path to use as the key
        std::string suffix = fileInfo.suffix().toStdString(); 
	if (!boost::iequals(suffix, "mrc") &&
	    !boost::iequals(suffix, "map"))
	{
	  std::cout << "in load_2DMRC(no matching file format): " << filepath << std::endl;
	  return false;
	}

        //extract just the filename from the path to use as the key
        std::string filename = 
          QFileInfo(QString::fromStdString(filepath))
          .baseName()
          .toStdString();

	// generate vrtmp directory and put temp img files in there.
	// this part must be fixed to use volume in memory without storing the file
	if(!std::system("mkdir vrtmp"))
        cvcapp.log(2,str(format("%s :: failing to generate \n") 
                         % BOOST_CURRENT_FUNCTION));

        std::string obase = "vrtmp/tmpimg";
        std::string ofile = obase;

        //make sure the key is unique
        CVC::DataMap map = cvcapp.data();
        std::string key = filename;
        int i = 0;
        while(map.find(key)!=map.end()) {
          key = filename + "_" + boost::lexical_cast<std::string>(i++);
          ofile = obase + "_" + boost::lexical_cast<std::string>(i++);
        }
	ofile += ".rawiv";

	using namespace std;
        cvcapp.log(2,str(format("%s :: load start\n") 
                         % BOOST_CURRENT_FUNCTION));
        // load image
	VolMagick::Volume img;
        VolMagick::readVolumeFile(cvcapp, img,filepath.c_str());
	img.map(0.0,255.0);
	img.voxelType(VolMagick::UChar);
	if(img.ZDim() > 1) {
           cvcapp.log(2,str(format("%s :: trying to load 3D MRC, please load it as a volume \n") 
                         % BOOST_CURRENT_FUNCTION));
           return false;
        }

	// write volume for re-read : this must be updated by copying data without file I/O
	// currently, if we set volume into data map, it doesn't initialize bounding box correctly.
	VolMagick::createVolumeFile(cvcapp, ofile.c_str(),
				  VolMagick::BoundingBox(0.0,0.0,0.0,img.XDim()-1.0, img.YDim()-1.0, 1.0),
				  img.dimension(),
				  std::vector<VolMagick::VoxelType>(1, img.voxelType()));

	VolMagick::writeVolumeFile(cvcapp, img,ofile);

        cvcapp.log(2,str(format("%s :: load done\n") 
                         % BOOST_CURRENT_FUNCTION));

	VolMagick::VolumeFileInfo vfi(cvcapp, ofile);

        //copy it to the data map
        cvcapp.data(key, vfi);

        //Set the viewers property for rendering the newly loaded VFI
        cvcapp.properties("viewers.vfi",key);
      }
    catch(std::runtime_error&)
      {
        return false;
      }
    std::cout << "Leave load function" << std::endl;
    return true;
  }

#ifdef USING_TILING
  // void zSpacingChanged(const std::string& key)
  // {
  //   log4cplus::Logger logger = log4cplus::Logger::getInstance("VolumeRover2.zSpacingChanged");
  //   LOG4CPLUS_TRACE(logger, "key = " << key);
  //   if (key.find("zSpacing") != string::npos) {
  //     LOG4CPLUS_TRACE(logger, "found zSpacing " << key);
  //     double zspacing = cvcapp.properties<double>(key);
  //     LOG4CPLUS_TRACE(logger, "zSpacing = " << zspacing);
  //     string dataset = key.substr(0, key.rfind('.'));
  //     LOG4CPLUS_TRACE(logger, "dataset = " << zspacing);
  //     cvcraw_geometry::contours_t c = cvcapp.data<cvcraw_geometry::contours_t>(dataset);
  //     c.set_z_spacing(zspacing);
  //     cvcapp.data(dataset, c);
  //   }
  // }

  bool load_ser(const std::string& filepath)
  {
    using namespace boost;
    static log4cplus::Logger logger = FUNCTION_LOGGER;

    try
      {
	QFileInfo fileInfo(QString::fromStdString(filepath));
        //extract just the filename from the path to use as the key
        std::string basename = fileInfo.completeBaseName().toStdString();
        std::string suffix = fileInfo.suffix().toStdString(); // ser
        std::string dir = fileInfo.dir().absolutePath().toStdString();
        std::string dirbasename = QFileInfo(QString::fromStdString(dir + "/" + basename))
	  .filePath().toStdString();

	if (!boost::iequals(suffix, "ser")) return false;

        //make sure the key is unique
        CVC::DataMap map = cvcapp.data();
        std::string key = basename;
        int i = 0;
        while(map.find(key)!=map.end())
          key = basename + "_" + boost::lexical_cast<std::string>(i++);

	using namespace std;
        LOG4CPLUS_TRACE(logger, "load start");
        // cvcapp.log(2,str(format("%s :: load start\n") 
        //                  % BOOST_CURRENT_FUNCTION));

	// Iterate over directory and find all files of the form
	// dirbasename.[0-9]* and find the smallest and largest numbers.
	QDir serDir = fileInfo.dir();
        QRegularExpression rx(".*" + QString::fromStdString(basename) + "\\.[0-9]+$");
	QStringList files = serDir.entryList();
	int start = 999999;
	int end = -1;
	for (QStringList::const_iterator it = files.begin(); it != files.end(); ++it) {
          string filename = it->toStdString();
          if (rx.match(*it).hasMatch()) {
              QFileInfo fi(*it);
              string sf = fi.suffix().toStdString();
              try {
                int idx = boost::lexical_cast<int>(fi.suffix().toStdString());
                start = min(start, idx);
                end = max(end, idx);
              } catch (boost::bad_lexical_cast& e) {
                // Do nothing
              }
            }
	}
	if (end == -1) {
	  throw runtime_error("No series files found");
	}

	list<Contour_handle> contours;
	list<Contour_exception> exceptions;
	vector<string> empty_components;
	// const int start = 59;
	// const int end = 160;
	const double smoothing_factor = -1;
	list<string> components;
	list<string> components_skip;
	string fn = dirbasename;
	read_contours_ser(fn, back_inserter(contours), back_inserter(exceptions),
			  start, end, smoothing_factor,
			  components.begin(), components.end(),
			  components_skip.begin(), components_skip.end());

	for (list<Contour_exception>::const_iterator it = exceptions.begin(); it != exceptions.end(); ++it) {
	  // cvcapp.log(4, str(format("Error in reading contours: %s\n") % it->what()));
	  // std::cout << "Error in reading contours: " << it->what();
	  LOG4CPLUS_WARN(logger, "Error in reading contours: " << it->what());
	}

	int slice_begin = start;
	int slice_end = end + 1;
	vector<Slice> slices(slice_end - slice_begin);
	for (list<Contour_handle>::iterator it = contours.begin(); it != contours.end(); ++it) {
	  const Contour_handle contour = *it;
	  int z = (int) contour->polygon()[0].z();
	  string component = contour->info().object_name();
	  // LOG4CPLUS_TRACE(logger, "Read " << component << ": " << pp(contour->polygon()));
	  Slice& slice = slices[z - slice_begin];
	  slice.push_back(component, contour);
	}

	set<string> component_names;
	for (list<Contour_handle>::iterator it = contours.begin(); it != contours.end(); ++it) {
	  const Contour_handle contour = *it;
	  string component = contour->info().object_name();
	  component_names.insert(component);
	}

	// for (set<string>::const_iterator it = component_names.begin(); it != component_names.end(); ++it) {
	//   string name = *it;
	//   // m_Contours[m_Variable->currentItem()][m_Timestep->value()][name] = *i;
	//   // QListWidgetItem *tmp = _ui->m_Objects->findItem(name.c_str(),0);
	//   // if(tmp) delete tmp;
	//   // _ui->m_Objects->insertItem(new QListWidgetItem(QString(name.c_str()), _ui->m_Objects, 0));
	//   // _ui->m_Objects->addItem(QString(name.c_str()));
	//   _ui->m_tileObjectsList->addItem(QString(name.c_str()));
	// }

        // cvcraw_geometry::cvcgeom_t cvcgeom(filepath);

        LOG4CPLUS_TRACE(logger, "load done");
        // cvcapp.log(2,str(format("%s :: load done\n") 
        //                  % BOOST_CURRENT_FUNCTION));
        //copy it to the data map
        // cvcapp.data(key, cvcgeom);
	cvcraw_geometry::contours_t contours_data(list<string>(component_names.begin(), 
							       component_names.end()), slices, 
						  slice_begin, slice_end-1, key);
        cvcapp.data(key, contours_data);
        // cvcapp.properties(key+".zSpacing", boost::lexical_cast<string>(contours_data.z_spacing()));
	// cvcapp.dataChanged
	// cvcapp.dataChanged.connect(MapChangeSignal::slot_type(zSpacingChanged));
	// cvcapp.propertiesChanged.connect(zSpacingChanged);
	
	// arand, 5-3-2011: updated with listPropertyAppend
	cvcapp.listPropertyAppend("zoomed.geometries",key);
	cvcapp.listPropertyAppend("thumbnail.geometries",key);
      }
    catch(std::runtime_error&)
      {
        return false;
      }
    return true;
  }
#endif
}

int main(int argc, char **argv)
{
  using namespace std;

  // Initialize log4cplus
  std::ifstream testfile("log4cplus.properties");
  if (testfile) {
    testfile.close();
    log4cplus::PropertyConfigurator::doConfigure("log4cplus.properties");
  }
  else {
    log4cplus::BasicConfigurator::doConfigure();
  }

  // Set properties for use in QSettings
  QCoreApplication::setOrganizationName("CVC");
  QCoreApplication::setOrganizationDomain("cvcweb.ices.utexas.edu");
  QCoreApplication::setApplicationName("VolumeRover");

  bool batchmode = false;

  string propertymapfile = "";
  string outputimagefile = "";

  // arand: check for batch mode...
  if (argc > 3) {
    if (boost::iequals(argv[1], "-batch")) {
      // setup for batch mode here...      
      batchmode = true;
      propertymapfile = argv[2];
      outputimagefile = argv[3];
    }
  }



  //Install loaders for volumes and geometry
  cvcapp.dataReader(load_volumeFileInfo);
  cvcapp.dataReader(load_cvcgeom);
  cvcapp.dataReader(load_image);
  cvcapp.dataReader(load_2DMRC);
#ifdef USING_TILING
  cvcapp.dataReader(load_ser);
#endif
  //Initialize the Qt part of the application
  CVC::Application app(argc,argv);

  //Install data widgets for objects that will end up in the datamap so
  //users can manipulate them.
  CVC::CVCMainWindow::instance().
    insertDataWidget<VolMagick::VolumeFileInfo>(new VolumeInterface);
  CVC::CVCMainWindow::instance().
    insertDataWidget<cvcraw_geometry::cvcgeom_t>(new GeometryInterface);
#ifdef USING_TILING
  CVC::CVCMainWindow::instance().
    insertDataWidget<cvcraw_geometry::contours_t>(new ContoursInterface);
#endif
  CVC::CVCMainWindow::instance().
    insertDataWidget<VolMagick::Volume>(new VolumeMemoryInterface);

  //Add the VolumeViewer container widget to the tab widget
  CVC::Viewers *viewers = new CVC::Viewers;
  viewers->setObjectName("viewers");

  //set logging output properties
  cvcapp.properties("system.log_output","stdout");
  cvcapp.properties("system.log_verbosity",3); // 0 is quiet

#ifdef USING_VOLUMEGRIDROVER
  cvcapp.log(5, "set volume grid rover pointer\n");
  viewers->setVolumeGridRoverPtr( CVC::CVCMainWindow::instance().volumeGridRover() );
  cvcapp.log(5, "set volume grid rover pointer end\n");
#endif

  CVC::CVCMainWindow::instance().
    tabWidget()->insertTab(0,
                           viewers,
                           "Viewers");
  //Make sure it is visible by default
  CVC::CVCMainWindow::instance().
    tabWidget()->setCurrentIndex(0);

  CVC::CVCMainWindow::instance().show();

#if 0
  //VolumeViewer2 testing...
  boost::scoped_ptr<CVC::VolumeViewer2> vv(new CVC::VolumeViewer2);
  vv->show();
#endif

  //Cache file to use
  QString cachePath(QDir::homePath() + "/.VolumeRover-cache.cvc");
  cvcapp.properties("system.cache_file",cachePath.toStdString());

  //Turn caching off for now until the gap bug is fixed -- transfix 10/09/2011
  cvcapp.properties("system.caching_enabled", "false" /*"true"*/);

  //Tell the system about the types we will be using with this application
  //so it can print programmer readable type names.
  cvcapp.registerDataType(cvcraw_geometry::cvcgeom_t);
  cvcapp.registerDataType(VolMagick::Volume);
  cvcapp.registerDataType(VolMagick::VolumeFileInfo);

  // arand: this version of "batchmode" is somewhat hacked...
  if (batchmode) {
    cvcapp.startThread("batch_mode_thread", BatchThread(propertymapfile, outputimagefile));
  }

  // Redirect stdout to the output window
  streambuf* stdoutbuf = cout.rdbuf();
  streambuf* stderrbuf = cerr.rdbuf();
  CVC::windowbuf* wbuf = new CVC::windowbuf(&CVC::CVCMainWindow::instance());//.createWindowStreamBuffer();
  CVC::teebuf* tbuf = new CVC::teebuf(wbuf, stdoutbuf);

  cout.rdbuf(tbuf);
  // cerr.rdbuf(&wbuf);

  // cout.rdbuf(stdoutbuf);
  // cerr.rdbuf(stderrbuf);
    
  int retval = app.exec();
  
  //Need to do this to delete the instance before QApplication is destroyed
  CVC::CVCMainWindow::terminate();

  return retval;
}
