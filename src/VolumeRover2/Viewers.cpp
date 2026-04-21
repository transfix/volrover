/*
  Copyright 2011 The University of Texas at Austin

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

/* $Id$ */

#include <VolumeRover2/Viewers.h>

#ifdef USING_VOLUMEGRIDROVER
#include <VolumeGridRover/VolumeGridRover.h>
#endif

#ifdef ISOCONTOURING_WITH_LBIE
#include <LBIE/LBIE_Mesher.h>
#endif

#ifdef ISOCONTOURING_WITH_FASTCONTOURING
#include <FastContouring/FastContouring.h>
#endif

#include <VolMagick/VolMagick.h>

#include <CVC/App.h>
#include <CVC/CVCEvent.h>
#include <CVC/Exception.h>
#include <CVC/log4cplus_compat.h>

#include <QTimer>
#include <QMessageBox>
#include <QSplitter>

#include <boost/lexical_cast.hpp>
#include <boost/thread/thread.hpp>
#include <boost/format.hpp>
#include <boost/bind/bind.hpp>

#include "ui_VolumeViewerPage.h"
#include "ui_VolumeViewerPageManipulators.h"

#include <map>
#include <set>
#include <iterator>

using namespace boost::placeholders;

namespace
{
#ifdef ISOCONTOURING_WITH_LBIE
  //utility function to convert from LBIE geoframe to cvcraw geometry
  // 07/29/2011 - transfix - added some thread feedback
  cvcraw_geometry::cvcgeom_t convert(const LBIE::geoframe& geo)
  {
    using namespace std;
    CVC::ThreadInfo ti(BOOST_CURRENT_FUNCTION);

    const double step_inc = 1.0/6.0;
    cvcapp.threadProgress(step_inc * 0.0);

    cvcraw_geometry::cvcgeom_t ret_geom;
    ret_geom.points().resize(geo.verts.size());
    copy(geo.verts.begin(),
	 geo.verts.end(),
	 ret_geom.points().begin());
    cvcapp.threadProgress(step_inc * 1.0);
    ret_geom.normals().resize(geo.normals.size());
    copy(geo.normals.begin(),
	 geo.normals.end(),
	 ret_geom.normals().begin());
    cvcapp.threadProgress(step_inc * 2.0);
    ret_geom.colors().resize(geo.color.size());
    copy(geo.color.begin(),
	 geo.color.end(),
	 ret_geom.colors().begin());
    cvcapp.threadProgress(step_inc * 3.0);
    ret_geom.boundary().resize(geo.bound_sign.size());
    for(vector<unsigned int>::const_iterator j = geo.bound_sign.begin();
	j != geo.bound_sign.end();
	j++)
      ret_geom.boundary()[distance(j,geo.bound_sign.begin())] = *j;
    cvcapp.threadProgress(step_inc * 4.0);
    ret_geom.triangles().resize(geo.triangles.size());
    copy(geo.triangles.begin(),
	 geo.triangles.end(),
	 ret_geom.triangles().begin());
    cvcapp.threadProgress(step_inc * 5.0);
    ret_geom.quads().resize(geo.quads.size());
    copy(geo.quads.begin(),
	 geo.quads.end(),
	 ret_geom.quads().begin());
    cvcapp.threadProgress(step_inc * 6.0);
    return ret_geom;
  }
#endif

#ifdef ISOCONTOURING_WITH_FASTCONTOURING
  // 07/29/2011 - transfix - added some thread feedback
  cvcraw_geometry::cvcgeom_t convert(const FastContouring::TriSurf& geo)
  {
    using namespace std;
    CVC::ThreadInfo ti(BOOST_CURRENT_FUNCTION);

    const double step_inc = 1.0/4.0;
    cvcapp.threadProgress(step_inc * 0.0);

    cvcraw_geometry::cvcgeom_t ret_geom;
    ret_geom.points().resize(geo.verts.size()/3);
    memcpy(&(ret_geom.points()[0]),
	   &(geo.verts[0]),
	   geo.verts.size()*sizeof(double));
    cvcapp.threadProgress(step_inc * 1.0);
    ret_geom.normals().resize(geo.normals.size()/3);
    memcpy(&(ret_geom.normals()[0]),
	   &(geo.normals[0]),
	   geo.normals.size()*sizeof(double));
    cvcapp.threadProgress(step_inc * 2.0);
    ret_geom.colors().resize(geo.colors.size()/3);
    memcpy(&(ret_geom.colors()[0]),
	   &(geo.colors[0]),
	   geo.colors.size()*sizeof(double));
    cvcapp.threadProgress(step_inc * 3.0);
    ret_geom.triangles().resize(geo.tris.size()/3);
    
    // arand: this doesn't work
    //        ret_geom.triangles() isn't a flat array?
    //memcpy(&(ret_geom.triangles()[0]),
    //	   &(geo.tris[0]),
    //	   geo.tris.size()*sizeof(unsigned int));


    // arand: this works, but but may be slow?
    for (int i=0; i<geo.tris.size()/3; i++) {
      ret_geom.triangles()[i][0] = geo.tris[3*i];
      ret_geom.triangles()[i][1] = geo.tris[3*i+1];
      ret_geom.triangles()[i][2] = geo.tris[3*i+2];
      //cout << geo.tris[3*i] << " " << geo.tris[3*i+1] << " " << geo.tris[3*i+2] << " " << ret_geom.triangles()[i][0] << " " << ret_geom.triangles()[i][1] << " " << ret_geom.triangles()[i][2] << endl;
    }

    cvcapp.threadProgress(step_inc * 4.0);
    return ret_geom;
  }
#endif

  const double MIN_RANGE = 0.0;
  const double MAX_RANGE = 1.0;
  template <typename T> static inline T clamp(T val)
  {
    return std::max(MIN_RANGE,std::min(MAX_RANGE,val));
  }

  class IsocontouringThread
  {
  public:
    IsocontouringThread(const std::string& threadKey,
			const std::string& viewerCtiDataKey,
                        const std::string& volDataKey,
                        const std::string& outputGeometryDataKey)
      : _threadKey(threadKey),
	_viewerCtiDataKey(viewerCtiDataKey),
        _volDataKey(volDataKey),
        _outputGeometryDataKey(outputGeometryDataKey) {}

    IsocontouringThread(const IsocontouringThread& t)
      : _threadKey(t._threadKey),
        _viewerCtiDataKey(t._viewerCtiDataKey),
        _volDataKey(t._volDataKey),
        _outputGeometryDataKey(t._outputGeometryDataKey) {}

    IsocontouringThread& operator=(const IsocontouringThread& t)
    {
      _threadKey = t._threadKey;
      _viewerCtiDataKey = t._viewerCtiDataKey;
      _volDataKey = t._volDataKey;
      _outputGeometryDataKey = t._outputGeometryDataKey;
      return *this;
    }

    // 07/29/2011 - transfix - fixed the thread progress feedback
    // 09/16/2011 - transfix - restarting the thread if the volume changed
    void operator()()
    {      
      using namespace boost;
      using namespace boost::algorithm;
      using namespace CVCColorTable;

      CVC::ThreadFeedback feedback;

      std::string ctiKey = _viewerCtiDataKey;
      ColorTable::color_table_info cti = cvcapp.data<ColorTable::color_table_info>(ctiKey);
      ColorTable::isocontour_nodes nodes = cti.isocontourNodes();
      ColorTable::color_nodes c_nodes = cti.colorNodes();

      std::string objName = _threadKey;

      VolMagick::Volume vol(cvcapp.data<VolMagick::Volume>(_volDataKey));

      if(!nodes.empty())
        {
          cvcraw_geometry::cvcgeom_t iso_geom;

          if(!cvcapp.hasProperty(objName + ".method") || //if no prop set, default to LBIE
             cvcapp.properties(objName + ".method") == "LBIE")
            {
              cvcapp.log(4,str(format("%s :: isocontouring using LBIE\n")
                               % BOOST_CURRENT_FUNCTION));

#ifdef ISOCONTOURING_WITH_LBIE
              LBIE::Mesher def_thumb_mesher;
              //TODO: allow other methods and mesh types
              def_thumb_mesher.extractionMethod(LBIE::Mesher::FASTCONTOURING);
              def_thumb_mesher.setVolume(vol);
              _meshers.resize(nodes.size(),def_thumb_mesher);
              int node_idx = 0;
              const double node_inc = nodes.size() <= 1 ? 1.0 : 
                1.0/nodes.size();
              for(ColorTable::isocontour_nodes::const_iterator i = nodes.begin();
                  i != nodes.end();
                  i++)
                {
                  int step = 0;
                  const double step_inc = 1.0/3.0;
                  cvcapp.threadProgress(node_idx*node_inc + (step*step_inc)/nodes.size());

                  double isoval = vol.min()+(vol.max()-vol.min())*i->position;
                  cvcraw_geometry::cvcgeom_t cur_geom = 
                    //isocache(isoval,vol.boundingBox());
                    cvcraw_geometry::cvcgeom_t();

                  if(cur_geom.empty())
                    {
                      _meshers[node_idx].isovalue(isoval);
                      LBIE::geoframe &extracted = 
                        _meshers[node_idx].extractMesh();
                      cur_geom = convert(extracted);
                    }
                  step++;
                  cvcapp.threadProgress(node_idx*node_inc + (step*step_inc)/nodes.size());
	      
                  //calculate a color for the isocontour via the color nodes
                  double r = 1.0, g = 1.0, b = 1.0;
                  if(c_nodes.size() >= 2)
                    {
                      ColorTable::color_nodes::const_iterator low_itr;
                      ColorTable::color_nodes::const_iterator high_itr = 
                        c_nodes.lower_bound(color_node(i->position));
                      if(high_itr != c_nodes.end())
                        {
                          low_itr =
                            high_itr == c_nodes.begin() ? high_itr : prior(high_itr);
                        }
                      ColorTable::color_node high = *high_itr;
                      ColorTable::color_node low = *low_itr;
                      double interval_pos = 
                        (i->position - low.position)/
                        (high.position - low.position);
                      r = clamp(low.r + (high.r - low.r)*interval_pos);
                      g = clamp(low.g + (high.g - low.g)*interval_pos);
                      b = clamp(low.b + (high.b - low.b)*interval_pos);
                    }
	      
                  cvcraw_geometry::geometry_t::color_t color = {{ r,g,b }};
                  cur_geom.colors().resize(cur_geom.points().size());
                  fill(cur_geom.colors().begin(),
                       cur_geom.colors().end(),
                       color);
                  step++;
                  cvcapp.threadProgress(node_idx*node_inc + (step*step_inc)/nodes.size());

                  //write it back to the cache -- there seems to be a problem with
                  //the cache but isocontouring seems fast enough to not need it now
                  //isocache(isoval,vol.boundingBox(),cur_geom);
                  iso_geom.merge(cur_geom);

                  node_idx++;
                }
              cvcapp.threadProgress(1.0);
#else
              cvcapp.log(2,str(format("%s :: LBIE isocontouring disabled\n")
                               % BOOST_CURRENT_FUNCTION));
#endif
            }
          else if(cvcapp.properties(objName + ".method") == "FastContouring")
            {
              cvcapp.log(4,str(format("%s :: isocontouring using FastContouring\n")
                               % BOOST_CURRENT_FUNCTION));

#ifdef ISOCONTOURING_WITH_FASTCONTOURING
              using namespace std;
              int node_idx = 0;
              const double node_inc = nodes.size() <= 1 ? 1.0 : 
                1.0/nodes.size();
              for(ColorTable::isocontour_nodes::const_iterator i = nodes.begin();
                  i != nodes.end();
                  i++)
                {
                  int step = 0;
                  const double step_inc = 1.0/3.0;
                  cvcapp.threadProgress(node_idx*node_inc + (step*step_inc)/nodes.size());

                  double isoval = vol.min()+(vol.max()-vol.min())*i->position;
                  cvcraw_geometry::cvcgeom_t cur_geom = 
                    //isocache(isoval,vol.boundingBox());
                    cvcraw_geometry::cvcgeom_t();

                  //calculate a color for the isocontour via the color nodes
                  double r = 1.0, g = 1.0, b = 1.0;
                  if(c_nodes.size() >= 2)
                    {
                      ColorTable::color_nodes::const_iterator low_itr;
                      ColorTable::color_nodes::const_iterator high_itr = 
                        c_nodes.lower_bound(color_node(i->position));
                      if(high_itr != c_nodes.end())
                        {
                          low_itr =
                            high_itr == c_nodes.begin() ? high_itr : prior(high_itr);
                        }
                      ColorTable::color_node high = *high_itr;
                      ColorTable::color_node low = *low_itr;
                      double interval_pos = 
                        (i->position - low.position)/
                        (high.position - low.position);
                      r = clamp(low.r + (high.r - low.r)*interval_pos);
                      g = clamp(low.g + (high.g - low.g)*interval_pos);
                      b = clamp(low.b + (high.b - low.b)*interval_pos);
                    }
                  step++;
                  cvcapp.threadProgress(node_idx*node_inc + (step*step_inc)/nodes.size());

                  if(cur_geom.empty()) {
                    _contourExtractor.setVolume(vol); // arand: is this slow???
                    cur_geom =                        // joe: shouldn't be, copy-on-write...
                      convert(
                        _contourExtractor.extractContour(isoval,r,g,b)
                      );
                  } else
                    {
                      //just fill in the color
                      cvcraw_geometry::geometry_t::color_t color = {{ r,g,b }};
                      cur_geom.colors().resize(cur_geom.points().size());
                      fill(cur_geom.colors().begin(),
                           cur_geom.colors().end(),
                           color);
                    }
                  step++;
                  cvcapp.threadProgress(node_idx*node_inc + (step*step_inc)/nodes.size());
                  
                  //write it back to the cache -- there seems to be a problem with
                  //the cache but isocontouring seems fast enough to not need it now
                  //isocache(isoval,vol.boundingBox(),cur_geom);
                  iso_geom.merge(cur_geom);

                  node_idx++;
                }
              cvcapp.threadProgress(1.0);
#else
              cvcapp.log(2,str(format("%s :: FastContouring isocontouring disabled\n")
                               % BOOST_CURRENT_FUNCTION));
#endif
            }
          else
            {
              cvcapp.log(1,str(format("%s :: Unknown isocontouring method %s\n")
                               % BOOST_CURRENT_FUNCTION
                               % cvcapp.properties(objName + ".method")));
            }

	  //if the nodes changed while we were busy, start a new thread!
	  if(nodes != 
	     cvcapp.data<ColorTable::color_table_info>(ctiKey).isocontourNodes())
            IsocontouringThread::start(*this);

          //if the volume changed while we were busy, start a new thread!
          if(vol !=
             cvcapp.data<VolMagick::Volume>(_volDataKey))
            IsocontouringThread::start(*this);            
          
          cvcapp.data(_outputGeometryDataKey,cvcraw_geometry::cvcgeom_t(iso_geom));
        }
      else
	{
	  // arand: fix to delete isocontours...
	  cvcraw_geometry::cvcgeom_t iso_geom;  
	  cvcapp.data(_outputGeometryDataKey,cvcraw_geometry::cvcgeom_t(iso_geom));
	}
    }

    static void start(const std::string& threadKey,
		      const std::string& viewerCtiDataKey,
		      const std::string& volDataKey,
		      const std::string& outputGeometryDataKey)
    {
      //Make a unique key name to use by adding a number to the key
      std::string uniqueThreadKey = threadKey;
      unsigned int i = 0;
      while(cvcapp.hasThread(uniqueThreadKey))
        uniqueThreadKey = 
          threadKey + boost::lexical_cast<std::string>(i++);

      cvcapp.threads(
	uniqueThreadKey,
	CVC::ThreadPtr(
	  new boost::thread(
	    IsocontouringThread(
	      uniqueThreadKey,
	      viewerCtiDataKey,
	      volDataKey,
	      outputGeometryDataKey
	    )
	  )
	)
      );
    }

    static void start(const IsocontouringThread& it)
    {
      cvcapp.threads(
        it.threadKey(),
	CVC::ThreadPtr(new boost::thread(IsocontouringThread(it)))
      );
    }

    //The isocache is used to avoid re-computing isocontours for all
    //nodes' isovalues, even if only 1 changed.
    typedef std::map<
      VolMagick::BoundingBox,
      cvcraw_geometry::cvcgeom_t
    > isogeom_t;

    typedef std::map<
      double, //isovalue
      isogeom_t
    > isocache_t;

    //to manipulate the isocache
    static cvcraw_geometry::cvcgeom_t isocache(double isoval, 
                                               const VolMagick::BoundingBox& bbox)
    {
      boost::mutex::scoped_lock lock(_isocacheMutex);
      return _isocache[isoval][bbox];
    }

    static void isocache(double isoval,
                         const VolMagick::BoundingBox& bbox,
                         const cvcraw_geometry::cvcgeom_t& geom)
    {
      boost::mutex::scoped_lock lock(_isocacheMutex);
      _isocache[isoval][bbox] = geom;
    }

    static isocache_t isocache()
    {
      boost::mutex::scoped_lock lock(_isocacheMutex);
      return _isocache;
    }

    static void isocache(const isocache_t& isoc)
    {
      boost::mutex::scoped_lock lock(_isocacheMutex);
      _isocache = isoc;
    }

    std::string threadKey() const 
    {
      boost::mutex::scoped_lock lock(_settingsMutex);
      return _threadKey; 
    }
    
    std::string viewerCtiDataKey() const
    {
      boost::mutex::scoped_lock lock(_settingsMutex);
      return _viewerCtiDataKey;
    }

    std::string volDataKey() const
    {
      boost::mutex::scoped_lock lock(_settingsMutex);
      return _volDataKey;
    }

    std::string outputGeometryDataKey() const
    {
      boost::mutex::scoped_lock lock(_settingsMutex);
      return _outputGeometryDataKey;
    }

    IsocontouringThread& threadKey(const std::string tk)
    {
      boost::mutex::scoped_lock lock(_settingsMutex);
      _threadKey = tk;
      return *this;
    }

    IsocontouringThread& viewerCtiDataKey(const std::string dk)
    {
      boost::mutex::scoped_lock lock(_settingsMutex);
      _viewerCtiDataKey = dk;
      return *this;
    }

    IsocontouringThread& volDataKey(const std::string dk)
    {
      boost::mutex::scoped_lock lock(_settingsMutex);
      _volDataKey = dk;
      return *this;
    }

    IsocontouringThread& outputGeometryDataKey(const std::string dk)
    {
      boost::mutex::scoped_lock lock(_settingsMutex);
      _outputGeometryDataKey = dk;
      return *this;
    }

  private:
    mutable boost::mutex _settingsMutex;
    std::string _threadKey;
    std::string _viewerCtiDataKey;
    std::string _volDataKey;
    std::string _outputGeometryDataKey;

#ifdef ISOCONTOURING_WITH_LBIE
    //for mesh extraction, 1 mesher per isocontour bar
    std::vector<LBIE::Mesher> _meshers;
#endif

#ifdef ISOCONTOURING_WITH_FASTCONTOURING
    FastContouring::ContourExtractor _contourExtractor;
#endif

    //Cache already computed isovalues.
    static boost::mutex _isocacheMutex;
    static isocache_t   _isocache;
  };

  boost::mutex IsocontouringThread::_isocacheMutex;
  IsocontouringThread::isocache_t IsocontouringThread::_isocache;

  // ----------------
  // ReadVolumeThread
  // ----------------
  // Purpose:
  //   Performs volume reads to produce the volume data needed for
  //   the viewers to visualize.
  // ---- Change History ----
  // 09/16/2011 -- Joe R. -- Initial implementation.
  class ReadVolumeThread
  {
  public:
    CVC_DEF_EXCEPTION(ReadVolumeThreadException);

    ReadVolumeThread(const std::string threadKey,
                     const std::string vfiDataKey,    //Can point to a vfi or a Volume.
                     const std::string bboxDataKey,   //The subvolume bbox to extract.
                     const std::string outVolDataKey, //The extracted Volume will end up here.
                     unsigned int variable_index = 0, //Var and time indices, ignored if
                     unsigned int timestep_index = 0) //vfiDataKey points to a Volume.
      : _threadKey(threadKey),
        _vfiDataKey(vfiDataKey),
        _bboxDataKey(bboxDataKey),
        _outVolDataKey(outVolDataKey),
        _variable(variable_index),
        _timestep(timestep_index) {}

    ReadVolumeThread(const ReadVolumeThread& t)
      : _threadKey(t._threadKey),
        _vfiDataKey(t._vfiDataKey),
        _bboxDataKey(t._bboxDataKey),
        _outVolDataKey(t._outVolDataKey),
        _variable(t._variable),
        _timestep(t._timestep) {}
        
    ReadVolumeThread& operator=(const ReadVolumeThread& t)
    {
      _threadKey     = t._threadKey;
      _vfiDataKey    = t._vfiDataKey;
      _bboxDataKey   = t._bboxDataKey;
      _outVolDataKey = t._outVolDataKey;
      _variable      = t._variable;
      _timestep      = t._timestep;
      return *this;
    }

    static void start(const std::string threadKey,
                      const std::string vfiDataKey,
                      const std::string bboxDataKey,
                      const std::string outVolDataKey,
                      unsigned int variable_index = 0,
                      unsigned int timestep_index = 0)
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
            ReadVolumeThread(
               uniqueThreadKey,
               vfiDataKey,
               bboxDataKey,
               outVolDataKey,
               variable_index,
               timestep_index
            )
          )
        )
      );
    }

    void operator()()
    {
      using namespace std;
      using namespace boost;

      CVC::ThreadFeedback feedback(BOOST_CURRENT_FUNCTION);

      cvcapp.log(1,str(format("%s :: updating :: %s\n")
                       % BOOST_CURRENT_FUNCTION
                       % _threadKey));
      
      try
        {
          VolMagick::BoundingBox bbox =
            cvcapp.data<VolMagick::BoundingBox>(_bboxDataKey);

          //If the vfi is set to point to an actual VolumeFileInfo object...
          if(cvcapp.isData<VolMagick::VolumeFileInfo>(_vfiDataKey))
            {
              VolMagick::VolumeFileInfo vfi =
                cvcapp.data<VolMagick::VolumeFileInfo>(_vfiDataKey);

              if(_variable >= vfi.numVariables() && vfi.numVariables()>0)
                {
                  cvcapp.log(2,str(format("%s: WARNING: variable index greater than number of"
                                          " variables in volume file, setting to last index\n")
                                   % BOOST_CURRENT_FUNCTION));
                  _variable = vfi.numVariables()-1;
                }

              if(_timestep >= vfi.numTimesteps() && vfi.numTimesteps()>0)
                {
                  cvcapp.log(2,str(format("%s: WARNING: volume_timestep greater than number of"
                                          " timesteps in volume file, setting to last index\n")
                                   % BOOST_CURRENT_FUNCTION));
                  _timestep = vfi.numTimesteps()-1;
                }

              VolMagick::Volume vol;
              VolMagick::readVolumeFile(cvcapp, vol,
                                        vfi.filename(),
                                        _variable,
                                        _timestep,
                                        bbox);

              cvcapp.data(_outVolDataKey,vol);
            }
          //else if it is pointing to an in-core Volume
          else if(cvcapp.isData<VolMagick::Volume>(_vfiDataKey))
            {
              if(_vfiDataKey == _outVolDataKey)
                throw ReadVolumeThreadException("ERROR: specified vfi data key is the same as the"
                                                " output volume data key, aborting!\n");
                
              VolMagick::Volume vol = cvcapp.data<VolMagick::Volume>(_vfiDataKey);
              cvcapp.data(_outVolDataKey,vol.sub(bbox));
            }
        }
      catch(std::exception& e)
        {
          std::string msg = str(boost::format("%s :: Error: %s")
                                % BOOST_CURRENT_FUNCTION
                                % e.what());
          cvcapp.log(1,msg+"\n");
        }
    }

  private:
    std::string  _threadKey;
    std::string  _vfiDataKey;
    std::string  _bboxDataKey;
    std::string  _outVolDataKey;
    unsigned int _variable;
    unsigned int _timestep;
  };
}

namespace CVC_NAMESPACE
{
  // 12/02/2011 - transfix - added a QSplitter between the viewers and the color table.
  Viewers::Viewers(QWidget *parent,
                   Qt::WindowFlags flags) : 
    QWidget(parent,flags),
    _thumbnailVolumeDirty(true),
    _subVolumeDirty(true),
    _thumbnailViewer(NULL),
    _subvolumeViewer(NULL),
    _subvolumeRenderQualitySlider(NULL),
    _subvolumeNearClipPlaneSlider(NULL),
    _thumbnailRenderQualitySlider(NULL),
    _thumbnailNearClipPlaneSlider(NULL),
    _colorTable(NULL),
    _ui(NULL),
    _uiManip(NULL),
#ifdef USING_VOLUMEGRIDROVER
    _volumeGridRoverPtr(NULL),
#endif
    _defaultSceneSet(false),
    _thumbnailPostInitFinished(false),
    _subvolumePostInitFinished(false)
  {
    _ui = new Ui::VolumeViewerPage;
    _ui->setupUi(this);

    _oldObjectName = objectName();

    QGridLayout *viewersFrameLayout = new QGridLayout(_ui->_viewersFrame);
    QSplitter *vsplitter = new QSplitter(_ui->_viewersFrame);
    vsplitter->setOrientation(Qt::Vertical);
    viewersFrameLayout->addWidget(vsplitter);
    
    //QFrame *viewers = new QFrame;
    QWidget *viewers = new QWidget;
    QGridLayout *viewersLayout = new QGridLayout(viewers);
    vsplitter->addWidget(viewers);
    //QFrame *manipulators = new QFrame;
    QWidget *manipulators = new QWidget;
    vsplitter->addWidget(manipulators);

    _subvolumeViewer = new VolumeViewer(NULL,flags);
    _thumbnailViewer = new VolumeViewer(NULL,flags);
    QSplitter *splitter = new QSplitter(viewers);
    viewersLayout->addWidget(splitter);
    splitter->addWidget(_subvolumeViewer);
    splitter->addWidget(_thumbnailViewer);

    _uiManip = new Ui::VolumeViewerPageManipulators;
    _uiManip->setupUi(manipulators);

    _subvolumeRenderQualitySlider = _uiManip->_subvolumeRenderQualitySlider;
    _subvolumeNearClipPlaneSlider = _uiManip->_subvolumeNearClipPlaneSlider;
    _thumbnailRenderQualitySlider = _uiManip->_thumbnailRenderQualitySlider;
    _thumbnailNearClipPlaneSlider = _uiManip->_thumbnailNearClipPlaneSlider;

    // arand, 6-14-2011: setting the slider to start in the middle
    //        this should match the default in the properties map even though they are 
    //        set in different parts of the code
    _subvolumeRenderQualitySlider->setValue((_subvolumeRenderQualitySlider->maximum()-
						    _subvolumeRenderQualitySlider->minimum())/2.0 +
						   _subvolumeRenderQualitySlider->minimum());
    _thumbnailRenderQualitySlider->setValue((_thumbnailRenderQualitySlider->maximum()-
						    _thumbnailRenderQualitySlider->minimum())/2.0 +
						   _thumbnailRenderQualitySlider->minimum());

    _colorTable = new CVCColorTable::ColorTable(_uiManip->_colortableFrame);
    QGridLayout *colorTableFrameLayout = new QGridLayout(_uiManip->_colortableFrame);
    colorTableFrameLayout->addWidget(_colorTable,0,0);
    _colorTable->interactiveUpdates(true);
    
    connect(_colorTable,SIGNAL(changed()),SLOT(updateColorTable()));

    connect(_subvolumeRenderQualitySlider,
            SIGNAL(valueChanged(int)),
            SLOT(setSubVolumeQuality(int)));
    connect(_subvolumeRenderQualitySlider,
            SIGNAL(sliderMoved(int)),
            SLOT(setSubVolumeQuality(int)));
    connect(_subvolumeNearClipPlaneSlider,
            SIGNAL(valueChanged(int)),
            SLOT(setSubVolumeNearPlane(int)));
    connect(_subvolumeNearClipPlaneSlider,
            SIGNAL(sliderMoved(int)),
            SLOT(setSubVolumeNearPlane(int)));
    connect(_thumbnailRenderQualitySlider,
            SIGNAL(valueChanged(int)),
            SLOT(setThumbnailQuality(int)));
    connect(_thumbnailRenderQualitySlider,
            SIGNAL(sliderMoved(int)),
            SLOT(setThumbnailQuality(int)));
    connect(_thumbnailNearClipPlaneSlider,
            SIGNAL(valueChanged(int)),
            SLOT(setThumbnailNearPlane(int)));
    connect(_thumbnailNearClipPlaneSlider,
            SIGNAL(sliderMoved(int)),
            SLOT(setThumbnailNearPlane(int)));

    _timer = new QTimer(this);
    _timer->start(100);
    connect(_timer,SIGNAL(timeout()),SLOT(timeout()));

    //Watch the global property map
    _propertiesConnection.disconnect();
    _propertiesConnection = 
      cvcapp.propertiesChanged.connect(
          MapChangeSignal::slot_type(
            &Viewers::propertiesChanged, this, _1
          )
        );

    //now read the properties map and pick up any settings from there
    handlePropertiesChanged("all");

    //Watch the global data map
    _dataConnection.disconnect();
    _dataConnection = 
      cvcapp.dataChanged.connect(
          MapChangeSignal::slot_type(
            &Viewers::dataChanged, this, _1
          )
        );

    //read the whole datamap and load necessary data
    handleDataChanged("all");

    //set the viewers' state to a sane default after they initialize
    connect(_subvolumeViewer,SIGNAL(postInit()),SLOT(setDefaultSubVolumeViewerState()));
    connect(_thumbnailViewer,SIGNAL(postInit()),SLOT(setDefaultThumbnailViewerState()));

    //Make sure the following types are registered so they show up nicely in the data map widget
    cvcapp.registerDataType(CVCColorTable::ColorTable::color_table_info);
    cvcapp.registerDataType(boost::shared_array<unsigned char>);
    cvcapp.registerDataType(boost::shared_array<float>);
  }

  Viewers::~Viewers()
  {
    delete _ui;
    delete _uiManip;
  }

  void Viewers::setDefaultSubVolumeViewerState()
  {
    if(!_defaultSceneSet)
      {
        setDefaultScene();
        _defaultSceneSet = true;
      }

    //so it responds to zoomed.* property changes
    _subvolumeViewer->setObjectName("zoomed");

    PropertyMap properties;
    properties["zoomed.rendering_mode"] = "colormapped";
    properties["zoomed.shaded_rendering_enabled"] = "false";
    properties["zoomed.draw_bounding_box"] = "true";
    properties["zoomed.draw_subvolume_selector"] = "false";
    // arand, 6-14-2011: changed default below to 0.5
    properties["zoomed.volume_rendering_quality"] = "0.5"; //[0.0,1.0]
    properties["zoomed.volume_rendering_near_plane"] = "0.0";
    properties["zoomed.projection_mode"] = "perspective";
    properties["zoomed.draw_corner_axis"] = "true";
    properties["zoomed.draw_geometry"] = "true";
    properties["zoomed.draw_volumes"] = "true";
    properties["zoomed.clip_geometry"] = "true";
    properties["zoomed.draw_geometry_normals"] = "false";
    properties["zoomed.geometries"] = 
      cvcapp.properties(getObjectName("zoomed_isocontour_geometry"));
    properties["zoomed.geometry_line_width"] = "1.2";
    properties["zoomed.volumes"] = 
      cvcapp.properties(getObjectName("subvolumes"));
    properties["zoomed.background_color"] = "#000000";

    #ifndef SYNC_THUMBNAIL_WITH_MULTITILESERVER
    //multi-tile server related properties
    properties["zoomed.syncCamera_with_multiTileServer"] = "false";
    properties["zoomed.syncTransferFunc_with_multiTileServer"] = "false";
    properties["zoomed.syncShadedRender_with_multiTileServer"] = "false";
    properties["zoomed.syncRenderMode_with_multiTileServer"] = "false";
    properties["zoomed.interactiveMode_with_multiTileServer"] = "false";
    properties["zoomed.syncMode_with_multiTileServer"] = "0";
    properties["zoomed.syncInitial_multiTileServer"] = "0";
    #endif

    //stereo related properties
    properties["zoomed.io_distance"] = "0.062";
    properties["zoomed.physical_distance_to_screen"] = "2.0";
    properties["zoomed.physical_screen_width"] = "1.8";
    properties["zoomed.focus_distance"] = "1000.0";

    properties["zoomed.fov"] = boost::lexical_cast<std::string>(M_PI/4.0f);

    //viewers_transfer_function should be a boost::shared_array<float>
    //on the data map.
    properties["zoomed.transfer_function"] = 
      cvcapp.properties(getObjectName("transfer_function"));

    cvcapp.addProperties(properties);

    //doesnt work like I want it :(
    //_subvolumeViewer->normalizeOnVolume(false);

    //Son't reset the subvolume view whenever it's volume changes.
    //Since it changes frequently because of subvolume box manipulation
    //in the thumbnail viewer, if we were to allow it to show entire
    //scene, thus changing the camera, it would make the thumbnail view
    //change too.
    _subvolumeViewer->showEntireSceneOnNormalize(false);

    _subvolumePostInitFinished = true;
  }
  
  void Viewers::setDefaultThumbnailViewerState()
  {
    if(!_defaultSceneSet)
      {
        setDefaultScene();
        _defaultSceneSet = true;
      }

    //so it responds to thumbnail.* property changes
    _thumbnailViewer->setObjectName("thumbnail");
    
    PropertyMap properties;
    properties["thumbnail.rendering_mode"] = "colormapped";    
    properties["thumbnail.shaded_rendering_enabled"] = "false";
    properties["thumbnail.draw_bounding_box"] = "true";
    properties["thumbnail.draw_subvolume_selector"] = "true";
    // arand, 6-14-2011: changed default below to 0.5
    properties["thumbnail.volume_rendering_quality"] = "0.5"; //[0.0,1.0]
    properties["thumbnail.volume_rendering_near_plane"] = "0.0";
    properties["thumbnail.projection_mode"] = "perspective";
    properties["thumbnail.draw_corner_axis"] = "true";
    properties["thumbnail.draw_geometry"] = "true";
    properties["thumbnail.draw_volumes"] = "true";
    properties["thumbnail.clip_geometry"] = "true";
    properties["thumbnail.draw_geometry_normals"] = "false";
    properties["thumbnail.geometries"] = 
      cvcapp.properties(getObjectName("thumbnail_isocontour_geometry"));
    properties["thumbnail.geometry_line_width"] = "1.2";
    properties["thumbnail.volumes"] = 
      cvcapp.properties(getObjectName("volumes"));
    properties["thumbnail.background_color"] = "#000000";
    //When the subvolume box of the thumbnail viewer is manipulated, the thumbnail
    //viewer it will write its bounding box data to the datamap using the datakey 
    //set as the value of its property "subvolume_box_data".
    properties["thumbnail.subvolume_box_data"] = 
      cvcapp.properties(getObjectName("subvolume_bounding_box"));

#ifdef SYNC_THUMBNAIL_WITH_MULTITILESERVER
    //multi-tile server related properties
    properties["thumbnail.syncCamera_with_multiTileServer"] = "false";
    properties["thumbnail.syncTransferFunc_with_multiTileServer"] = "false";
    properties["thumbnail.syncShadedRender_with_multiTileServer"] = "false";
    properties["thumbnail.syncRenderMode_with_multiTileServer"] = "false";
    properties["thumbnail.interactiveMode_with_multiTileServer"] = "false";
    properties["thumbnail.syncMode_with_multiTileServer"] = "0";
    properties["thumbnail.syncInitial_multiTileServer"] = "0";
#endif

    //stereo related properties
    properties["thumbnail.io_distance"] = "0.062";
    properties["thumbnail.physical_distance_to_screen"] = "2.0";
    properties["thumbnail.physical_screen_width"] = "1.8";
    properties["thumbnail.focus_distance"] = "1000.0";

    properties["thumbnail.fov"] = boost::lexical_cast<std::string>(M_PI/4.0f);

    //viewers_transfer_function should be a boost::shared_array<float>
    //on the data map.
    properties["thumbnail.transfer_function"] = 
      cvcapp.properties(getObjectName("transfer_function"));

    cvcapp.addProperties(properties);

    //Makes the subvolume conform to the thumbnail viewer as volumes
    //are loaded into it.
    //TODO: not really needed now that the camera is synced with the property
    //map before each render
    //_thumbnailViewer->copyCameraOnNormalize(true);

    _thumbnailPostInitFinished = true;
  }

  void Viewers::setThumbnailQuality(int q)
  {
    cvcapp.properties(
       "thumbnail.volume_rendering_quality",
       float(q-_thumbnailRenderQualitySlider->minimum())/
       float(_thumbnailRenderQualitySlider->maximum()-
             _thumbnailRenderQualitySlider->minimum()));
  }

  void Viewers::setSubVolumeQuality(int q)
  {
    cvcapp.properties(
       "zoomed.volume_rendering_quality",
       float(q-_subvolumeRenderQualitySlider->minimum())/
       float(_subvolumeRenderQualitySlider->maximum()-
             _subvolumeRenderQualitySlider->minimum()));
  }

  void Viewers::setThumbnailNearPlane(int q)
  {
    cvcapp.properties(
       "thumbnail.volume_rendering_near_plane",
       float(q-_thumbnailNearClipPlaneSlider->minimum())/
       float(_thumbnailNearClipPlaneSlider->maximum()-
             _thumbnailNearClipPlaneSlider->minimum()));
  }

  void Viewers::setSubVolumeNearPlane(int q)
  {
    cvcapp.properties(
       "zoomed.volume_rendering_near_plane",
       float(q-_subvolumeNearClipPlaneSlider->minimum())/
       float(_subvolumeNearClipPlaneSlider->maximum()-
             _subvolumeNearClipPlaneSlider->minimum()));
  }

  void Viewers::setDefaultScene()
  {
    setDefaultOptions();

    std::string thumbnail_isocontour_geometry = 
      cvcapp.properties(getObjectName("thumbnail_isocontour_geometry"));
    std::string zoomed_isocontour_geometry = 
      cvcapp.properties(getObjectName("zoomed_isocontour_geometry"));

    //Make sure the data exist that viewers expect to have available
    if(!cvcapp.isData<cvcraw_geometry::cvcgeom_t>(thumbnail_isocontour_geometry))
      cvcapp.data(thumbnail_isocontour_geometry,cvcraw_geometry::cvcgeom_t());
    if(!cvcapp.isData<cvcraw_geometry::cvcgeom_t>(zoomed_isocontour_geometry))
      cvcapp.data(zoomed_isocontour_geometry,cvcraw_geometry::cvcgeom_t());
    ensureVolumeAvailability();
    ensureSubVolumeAvailability();

    CVC_NAMESPACE::PropertyMap properties;

    //geometry options
    properties[thumbnail_isocontour_geometry + ".visible"] = "true";
    properties[thumbnail_isocontour_geometry + ".render_mode"] = "triangles";
    properties[zoomed_isocontour_geometry + ".visible"] = "true";
    properties[zoomed_isocontour_geometry + ".render_mode"] = "triangles";

    cvcapp.addProperties(properties);

    _colorTable->info() = CVCColorTable::ColorTable::default_transfer_function();
    _colorTable->update(); //redraw the table since it isn't automatic.  Change this!
    updateColorTable(); //copy the default color table to the datamap
  }

  // 11/04/2011 - transfix - syncing viewers by default now that it's working well.
  void Viewers::setDefaultOptions()
  {
    //Set up initial Viewer properties
    {
      CVC_NAMESPACE::PropertyMap properties;

      //Viewers options for subvolume extraction and color table handling

      //the VolumeFileInfo for the volume(s) to show in the viewers
      properties[getObjectName("vfi")] = "none";

      //set the default rendering mode so we know what we need to extract
      properties[getObjectName("rendering_mode")] = "colormapped";

      //indices of the volumes to render in the VFI
      properties[getObjectName("single_volume_index")] = "0";
      properties[getObjectName("red_volume_index")]    = "0";
      properties[getObjectName("green_volume_index")]  = "1";
      properties[getObjectName("blue_volume_index")]   = "2";
      properties[getObjectName("alpha_volume_index")]  = "3";
      properties[getObjectName("volume_timestep")]     = "0";
      
      properties[getObjectName("subvolume_bounding_box")] = "thumbnail_subvolume_boundingbox";
      properties[getObjectName("thumbnail_bounding_box")] = "overall_bounding_box";

      properties[getObjectName("color_table_opacity_cubed")] = "true";
      //properties[getObjectName("color_table_opacity_cubed")] = "false";
      
      //After reading volume data from a volume via the VFI, the data
      //will have the following names.  The values are comma separated lists
      properties[getObjectName("volumes")] = "thumbnail_volume";
      properties[getObjectName("subvolumes")] = "zoomed_volume";

      //Setting these to a non-zero integer will trigger an update to the thumbnail
      //or subvolume volumes
      properties[getObjectName("thumbnail_dirty")] = "1";
      properties[getObjectName("subvolume_dirty")] = "1";

      //The transfer function data names
      properties[getObjectName("transfer_function")] = "viewers_transfer_function";
      properties[getObjectName("color_table_info")] = "viewers_color_table_info";
      properties[getObjectName("colortable_interactive_updates")] = "true";

      //The volumes to use for isocontouring
      properties[getObjectName("thumbnail_isocontouring_volume")] = "thumbnail_volume";
      properties[getObjectName("zoomed_isocontouring_volume")] = "zoomed_volume";
      
      //The isocontouring thread names
      properties[getObjectName("thumbnail_isocontouring_thread")] = "thumbnail_isocontouring";
      properties[getObjectName("zoomed_isocontouring_thread")] = "zoomed_isocontouring";

      //The isocontouring output geometry data names
      properties[getObjectName("thumbnail_isocontour_geometry")] = "thumbnail_isocontour";
      properties[getObjectName("zoomed_isocontour_geometry")] = "zoomed_isocontour";

      //If nonzero, sync viewer angles
      properties[getObjectName("sync_viewers")] = "1";

      //special property to have the subvolume fit the screen when syncing
      properties[getObjectName("subvolume_fit_screen")] = "1";

      cvcapp.addProperties(properties);
    }

    //Properties for objects that are used by the Viewers object
    {
      CVC_NAMESPACE::PropertyMap properties;

      //The methods to use for isocontouring
      std::string thumbnail_isocontouring_thread = 
        cvcapp.properties(getObjectName("thumbnail_isocontouring_thread"));
      std::string zoomed_isocontouring_thread = 
        cvcapp.properties(getObjectName("zoomed_isocontouring_thread"));
      properties[thumbnail_isocontouring_thread + ".method"] = "LBIE";
      properties[zoomed_isocontouring_thread + ".method"] = "LBIE";

      cvcapp.addProperties(properties);
    }
  }

  std::string Viewers::getObjectName(const std::string& property) const
  {
    std::string objName = objectName().toStdString();
    if(objName.empty()) objName = "Viewers";
    if(!property.empty()) objName = objName + "." + property;
    return objName;
  }

  void Viewers::markThumbnailDirty(bool flag)
  {
    _thumbnailVolumeDirty = flag;
    cvcapp.properties<int>(getObjectName("thumbnail_dirty"),flag?1:0);
  }

  void Viewers::markSubVolumeDirty(bool flag)
  {
    _subVolumeDirty = flag;
    cvcapp.properties<int>(getObjectName("subvolume_dirty"),flag?1:0);
  }

  // 09/17/2011 - transfix - added
  void Viewers::loadThumbnail()
  {
    using namespace std;
    using namespace boost;

    static log4cplus::Logger logger = FUNCTION_LOGGER;

    try
      {
        string vfiDataKey = cvcapp.properties(getObjectName("vfi"));
        //if it doesn't point to a valid data, then resetting the loaded data using
        //geometry bounding boxes
        if(!cvcapp.isData<VolMagick::VolumeFileInfo>(vfiDataKey) &&
           !cvcapp.isData<VolMagick::Volume>(vfiDataKey))
          {
            LOG4CPLUS_TRACE(logger, "no VolumeFileInfo to load, resetting loaded data");
	    // cvcapp.log(2,str(format("%s: no VolumeFileInfo to load, resetting loaded data\n")
	    //     	     % BOOST_CURRENT_FUNCTION));

	    vector<string> data_keys = cvcapp.listProperty(getObjectName("volumes"));

            vector<cvcraw_geometry::cvcgeom_t> geoms =
              _thumbnailViewer->getGeometriesFromDatamap();
            LOG4CPLUS_TRACE(logger, "geoms.size() == " << geoms.size());
            // cvcapp.log(3,str(format("%s :: geoms.size() == %s\n")
            //                  % BOOST_CURRENT_FUNCTION
            //                  % geoms.size()));

            if(geoms.empty())
              {
                for (const auto& data_key : data_keys)
                  cvcapp.data(data_key,VolMagick::Volume());
              }
            else
              {
                VolMagick::BoundingBox bbox;
                cvcraw_geometry::cvcgeom_t initial_geom;

                //find a non empty geom to supply an initial bbox
                for (const auto& geom : geoms)
                  if(!geom.empty())
                    {
                      initial_geom = geom;
                      break;
                    }

                cvcraw_geometry::cvcgeom_t::point_t minpt;
                cvcraw_geometry::cvcgeom_t::point_t maxpt;
                if(initial_geom.empty())
                  {
                    string data_key;
                    if(!data_keys.empty()) 
                      data_key = data_keys[0];
                    bbox = cvcapp.data<VolMagick::Volume>(data_key).boundingBox();
                  }
                else
                  {
                    minpt = initial_geom.min_point();
                    maxpt = initial_geom.max_point();
		    // arand: slightly enlarging
		    double eps0 = (maxpt[0]-minpt[0])/20.0;
		    double eps1 = (maxpt[1]-minpt[1])/20.0;
		    double eps2 = (maxpt[2]-minpt[2])/20.0;

                    bbox = VolMagick::BoundingBox(minpt[0]-eps0,minpt[1]-eps1,minpt[2]-eps2,
                                                  maxpt[0]+eps0,maxpt[1]+eps1,maxpt[2]+eps2);
                  }

                //build a bounding box that encompasses all bounding boxes of the geometries
                for (const auto& geom : geoms) {
                    if(geom.empty()) continue;
                    minpt = geom.min_point();
                    maxpt = geom.max_point();
                    VolMagick::BoundingBox geobox(minpt[0],minpt[1],minpt[2],
                                                  maxpt[0],maxpt[1],maxpt[2]);

                    LOG4CPLUS_TRACE(logger, str(boost::format("geobox: (%f,%f,%f) (%f,%f,%f)")
                                                % geobox.minx % geobox.miny % geobox.minz
                                                % geobox.maxx % geobox.maxy % geobox.maxz));
                    LOG4CPLUS_TRACE(logger, str(boost::format("bbox: (%f,%f,%f) (%f,%f,%f)")
                                                % bbox.minx % bbox.miny % bbox.minz
                                                % bbox.maxx % bbox.maxy % bbox.maxz));
                    // cvcapp.log(5,str(boost::format("geobox: (%f,%f,%f) (%f,%f,%f)")
                    //                  % geobox.minx % geobox.miny % geobox.minz
                    //                  % geobox.maxx % geobox.maxy % geobox.maxz));
                    // cvcapp.log(5,str(boost::format("bbox: (%f,%f,%f) (%f,%f,%f)")
                    //                  % bbox.minx % bbox.miny % bbox.minz
                    //                  % bbox.maxx % bbox.maxy % bbox.maxz));
       
                    bbox += geobox;
                  }
                    
                VolMagick::Volume newvol;
                newvol.boundingBox(bbox);
                for (const auto& data_key : data_keys)
                  cvcapp.data(data_key,newvol);
              }
          }
        else
          {
            unsigned int single_volume_index = 0;
            unsigned int red_volume_index = 0;
            unsigned int green_volume_index = 0;
            unsigned int blue_volume_index = 0;
            unsigned int alpha_volume_index = 0;
            unsigned int volume_timestep = 0;

            VolMagick::BoundingBox bbox;
            if(cvcapp.isData<VolMagick::VolumeFileInfo>(vfiDataKey))
              {
                VolMagick::VolumeFileInfo vfi =
                  cvcapp.data<VolMagick::VolumeFileInfo>(vfiDataKey);

                //initialize the indicies from the property map
                if(cvcapp.properties(getObjectName("rendering_mode")) == "colormapped")
                  {
                    single_volume_index =
                      cvcapp.properties<unsigned int>(getObjectName("single_volume_index"));
                    volume_timestep =
                      cvcapp.properties<unsigned int>(getObjectName("volume_timestep"));
            
                    if(single_volume_index >= vfi.numVariables() && vfi.numVariables()>0)
                      {
                        cvcapp.log(2,str(format("%s: WARNING: single_volume_index greater than number of variables in volume file, setting to last index\n")
                                         % BOOST_CURRENT_FUNCTION));
                        single_volume_index = vfi.numVariables()-1;
                      }
                    if(volume_timestep >= vfi.numTimesteps() && vfi.numTimesteps()>0)
                      {
                        cvcapp.log(2,str(format("%s: WARNING: volume_timestep greater than number of timesteps in volume file, setting to last index\n")
                                         % BOOST_CURRENT_FUNCTION));
                        volume_timestep = vfi.numTimesteps()-1;
                      }
                  }
                else if(cvcapp.properties(getObjectName("rendering_mode")) == "rgba" ||
                        cvcapp.properties(getObjectName("rendering_mode")) == "RGBA")
                  {
                    red_volume_index =
                      cvcapp.properties<unsigned int>(getObjectName("red_volume_index"));
                    green_volume_index =
                      cvcapp.properties<unsigned int>(getObjectName("green_volume_index"));
                    blue_volume_index =
                      cvcapp.properties<unsigned int>(getObjectName("blue_volume_index"));
                    alpha_volume_index =
                      cvcapp.properties<unsigned int>(getObjectName("alpha_volume_index"));
                    volume_timestep =
                      cvcapp.properties<unsigned int>(getObjectName("volume_timestep"));

                    if(red_volume_index >= vfi.numVariables() && vfi.numVariables()>0)
                      {
                        cvcapp.log(2,str(format("%s: WARNING: red_volume_index greater than number of variables in volume file, setting to last index\n")
                                         % BOOST_CURRENT_FUNCTION));
                        red_volume_index = vfi.numVariables()-1;
                      }
                    if(green_volume_index >= vfi.numVariables() && vfi.numVariables()>0)
                      {
                        cvcapp.log(2,str(format("%s: WARNING: green_volume_index greater than number of variables in volume file, setting to last index\n")
                                         % BOOST_CURRENT_FUNCTION));
                        green_volume_index = vfi.numVariables()-1;
                      }
                    if(blue_volume_index >= vfi.numVariables() && vfi.numVariables()>0)
                      {
                        cvcapp.log(2,str(format("%s: WARNING: blue_volume_index greater than number of variables in volume file, setting to last index\n")
                                         % BOOST_CURRENT_FUNCTION));
                        blue_volume_index = vfi.numVariables()-1;
                      }
                    if(alpha_volume_index >= vfi.numVariables() && vfi.numVariables()>0)
                      {
                        cvcapp.log(2,str(format("%s: WARNING: alpha_volume_index greater than number of variables in volume file, setting to last index\n")
                                         % BOOST_CURRENT_FUNCTION));
                        alpha_volume_index = vfi.numVariables()-1;
                      }
                    if(volume_timestep >= vfi.numTimesteps() && vfi.numTimesteps()>0)
                      {
                        cvcapp.log(2,str(format("%s: WARNING: volume_timestep greater than number of timesteps in volume file, setting to last index\n")
                                         % BOOST_CURRENT_FUNCTION));
                        volume_timestep = vfi.numTimesteps()-1;
                      }
                  }

                bbox = vfi.boundingBox();
              }
            else if(cvcapp.isData<VolMagick::Volume>(vfiDataKey))
              {
                VolMagick::Volume vol =
                  cvcapp.data<VolMagick::Volume>(vfiDataKey);
                bbox = vol.boundingBox();
              }

            //make sure we have the thumbnail_bounding_box property pointing to data
            cvcapp.data(cvcapp.properties(getObjectName("thumbnail_bounding_box")),
                        bbox);

            vector<string> volDataKeys = cvcapp.listProperty(getObjectName("volumes"));
            if(volDataKeys.empty())
              {
                volDataKeys.resize(1);
                volDataKeys[0] = "thumbnail_volume";
              }

            if(cvcapp.properties(getObjectName("rendering_mode")) == "colormapped")
              {
                ReadVolumeThread::start(std::string("ReadVolumeThread_")+volDataKeys[0],
                                        vfiDataKey,
                                        cvcapp.properties(getObjectName("thumbnail_bounding_box")),
                                        volDataKeys[0],
                                        single_volume_index,
                                        volume_timestep);
              }
            else if(cvcapp.properties(getObjectName("rendering_mode")) == "rgba" ||
                    cvcapp.properties(getObjectName("rendering_mode")) == "RGBA")
              {
                unsigned int vol_indices[4] =
                  {
                    red_volume_index,
                    green_volume_index,
                    blue_volume_index,
                    alpha_volume_index
                  };

                for(int i = 0; i < 4; i++)
                  ReadVolumeThread::start(std::string("ReadVolumeThread_")+volDataKeys[i],
                                          vfiDataKey,
                                          cvcapp.properties(getObjectName("thumbnail_bounding_box")),
                                          volDataKeys[i],
                                          vol_indices[i],
                                          volume_timestep);
              }         
          }
      }
    catch(std::exception& e)
      {
	std::string msg = str(boost::format("%s :: Error: %s")
			 % BOOST_CURRENT_FUNCTION
			 % e.what());
	cvcapp.log(1,msg+"\n");
	QMessageBox::critical(this,"Error",
			      QString("%1")
			      .arg(QString::fromStdString(msg)));
      }

    markThumbnailDirty(false);
  }

  // 09/17/2011 - transfix - added
  void Viewers::loadSubVolume()
  {
    using namespace std;
    using namespace boost;

    static log4cplus::Logger logger = FUNCTION_LOGGER;

    try
      {
        string vfiDataKey = cvcapp.properties(getObjectName("vfi"));
        string bboxDataKey = cvcapp.properties(getObjectName("subvolume_bounding_box"));
        //if it doesn't point to a valid data, then resetting the loaded data using
        //geometry bounding boxes
        if(!cvcapp.isData<VolMagick::VolumeFileInfo>(vfiDataKey) &&
           !cvcapp.isData<VolMagick::Volume>(vfiDataKey))
          {
            LOG4CPLUS_TRACE(logger, "no VolumeFileInfo to load, resetting loaded data");
	    // cvcapp.log(2,str(format("%s: no VolumeFileInfo to load, resetting loaded data\n")
	    //     	     % BOOST_CURRENT_FUNCTION));

            //use the thumbnail volume as a template
            vector<string> thumbnail_data_keys = cvcapp.listProperty(getObjectName("volumes"));
            vector<string> zoomed_data_keys    = cvcapp.listProperty(getObjectName("subvolumes"));

            VolMagick::BoundingBox bbox = cvcapp.data<VolMagick::BoundingBox>(bboxDataKey);

            VolMagick::Volume newvol;

            if(!thumbnail_data_keys.empty())
              newvol = cvcapp.data<VolMagick::Volume>(thumbnail_data_keys[0]);
                
            newvol.boundingBox(bbox);
            for (const auto& zoomed_data_key : zoomed_data_keys)
              cvcapp.data(zoomed_data_key,newvol);
          }
        else
          {
            unsigned int single_volume_index = 0;
            unsigned int red_volume_index = 0;
            unsigned int green_volume_index = 0;
            unsigned int blue_volume_index = 0;
            unsigned int alpha_volume_index = 0;
            unsigned int volume_timestep = 0;

            if(cvcapp.isData<VolMagick::VolumeFileInfo>(vfiDataKey))
              {
                VolMagick::VolumeFileInfo vfi =
                  cvcapp.data<VolMagick::VolumeFileInfo>(vfiDataKey);

                //initialize the indicies from the property map
                if(cvcapp.properties(getObjectName("rendering_mode")) == "colormapped")
                  {
                    single_volume_index =
                      cvcapp.properties<unsigned int>(getObjectName("single_volume_index"));
                    volume_timestep =
                      cvcapp.properties<unsigned int>(getObjectName("volume_timestep"));
            
                    if(single_volume_index >= vfi.numVariables() && vfi.numVariables()>0)
                      {
                        cvcapp.log(2,str(format("%s: WARNING: single_volume_index greater than number of variables in volume file, setting to last index\n")
                                         % BOOST_CURRENT_FUNCTION));
                        single_volume_index = vfi.numVariables()-1;
                      }
                    if(volume_timestep >= vfi.numTimesteps() && vfi.numTimesteps()>0)
                      {
                        cvcapp.log(2,str(format("%s: WARNING: volume_timestep greater than number of timesteps in volume file, setting to last index\n")
                                         % BOOST_CURRENT_FUNCTION));
                        volume_timestep = vfi.numTimesteps()-1;
                      }
                  }
                else if(cvcapp.properties(getObjectName("rendering_mode")) == "rgba" ||
                        cvcapp.properties(getObjectName("rendering_mode")) == "RGBA")
                  {
                    red_volume_index =
                      cvcapp.properties<unsigned int>(getObjectName("red_volume_index"));
                    green_volume_index =
                      cvcapp.properties<unsigned int>(getObjectName("green_volume_index"));
                    blue_volume_index =
                      cvcapp.properties<unsigned int>(getObjectName("blue_volume_index"));
                    alpha_volume_index =
                      cvcapp.properties<unsigned int>(getObjectName("alpha_volume_index"));
                    volume_timestep =
                      cvcapp.properties<unsigned int>(getObjectName("volume_timestep"));

                    if(red_volume_index >= vfi.numVariables() && vfi.numVariables()>0)
                      {
                        cvcapp.log(2,str(format("%s: WARNING: red_volume_index greater than number of variables in volume file, setting to last index\n")
                                         % BOOST_CURRENT_FUNCTION));
                        red_volume_index = vfi.numVariables()-1;
                      }
                    if(green_volume_index >= vfi.numVariables() && vfi.numVariables()>0)
                      {
                        cvcapp.log(2,str(format("%s: WARNING: green_volume_index greater than number of variables in volume file, setting to last index\n")
                                         % BOOST_CURRENT_FUNCTION));
                        green_volume_index = vfi.numVariables()-1;
                      }
                    if(blue_volume_index >= vfi.numVariables() && vfi.numVariables()>0)
                      {
                        cvcapp.log(2,str(format("%s: WARNING: blue_volume_index greater than number of variables in volume file, setting to last index\n")
                                         % BOOST_CURRENT_FUNCTION));
                        blue_volume_index = vfi.numVariables()-1;
                      }
                    if(alpha_volume_index >= vfi.numVariables() && vfi.numVariables()>0)
                      {
                        cvcapp.log(2,str(format("%s: WARNING: alpha_volume_index greater than number of variables in volume file, setting to last index\n")
                                         % BOOST_CURRENT_FUNCTION));
                        alpha_volume_index = vfi.numVariables()-1;
                      }
                    if(volume_timestep >= vfi.numTimesteps() && vfi.numTimesteps()>0)
                      {
                        cvcapp.log(2,str(format("%s: WARNING: volume_timestep greater than number of timesteps in volume file, setting to last index\n")
                                         % BOOST_CURRENT_FUNCTION));
                        volume_timestep = vfi.numTimesteps()-1;
                      }
                  }
              }

            vector<string> volDataKeys = cvcapp.listProperty(getObjectName("subvolumes"));
            if(volDataKeys.empty())
              {
                volDataKeys.resize(1);
                volDataKeys[0] = "zoomed_volume";
              }

            if(cvcapp.properties(getObjectName("rendering_mode")) == "colormapped")
              {
                ReadVolumeThread::start(std::string("ReadVolumeThread_")+volDataKeys[0],
                                        vfiDataKey,
                                        bboxDataKey,
                                        volDataKeys[0],
                                        single_volume_index,
                                        volume_timestep);
              }
            else if(cvcapp.properties(getObjectName("rendering_mode")) == "rgba" ||
                    cvcapp.properties(getObjectName("rendering_mode")) == "RGBA")
              {
                unsigned int vol_indices[4] =
                  {
                    red_volume_index,
                    green_volume_index,
                    blue_volume_index,
                    alpha_volume_index
                  };

                for(int i = 0; i < 4; i++)
                  ReadVolumeThread::start(std::string("ReadVolumeThread_")+volDataKeys[i],
                                          vfiDataKey,
                                          bboxDataKey,
                                          volDataKeys[i],
                                          vol_indices[i],
                                          volume_timestep);
              }         
          }
      }
    catch(std::exception& e)
      {
	std::string msg = str(boost::format("%s :: Error: %s")
			 % BOOST_CURRENT_FUNCTION
			 % e.what());
	cvcapp.log(1,msg+"\n");
	QMessageBox::critical(this,"Error",
			      QString("%1")
			      .arg(QString::fromStdString(msg)));
      }
    
    markSubVolumeDirty(false);
  }

#ifdef USING_VOLUMEGRIDROVER
  void Viewers::setVolumeGridRoverPtr( VolumeGridRover *ptr )
  {
    if( ptr == NULL )
      fprintf( stderr, "VolumeGridRover Pointer is NULL\n");
    else {
      _volumeGridRoverPtr = ptr;
    }
  }
#endif

  // 10/28/2011 - transfix - syncing viewers once they're finished initializing
  // 11/04/2011 - transfix - fixed a bug preventing syncing from being turned off
  void Viewers::timeout()
  {
    if(_thumbnailVolumeDirty)
      QCoreApplication::postEvent(this,
        new CVCEvent("loadThumbnail")
      );

    if(_subVolumeDirty)
      QCoreApplication::postEvent(this,
        new CVCEvent("loadSubVolume")
      );

    //if our name changed, reset the options
    if(_oldObjectName != objectName())
      {
	setDefaultOptions();

	//remove all the properties under the old name
	//TODO: not really comfortable with this, but this is the
	//result of the decision to have this object set it's own
	//default options.  Might revisit this later if it introduces
	//bugs.
	PropertyMap map = cvcapp.properties();
        for (const auto& val : map) {
	    using namespace std;
	    using namespace boost;
	    using namespace boost::algorithm;

	    vector<string> key_idents;
	    split(key_idents, val.first, is_any_of("."));
	    //check if the key was meant for this object
	    if(key_idents.size() == 2 && 
	       key_idents[0] == _oldObjectName.toStdString())
	      cvcapp.properties(val.first,"");
	  }

	_oldObjectName = objectName();
      }

    if(_thumbnailPostInitFinished &&
       _subvolumePostInitFinished &&
       !_thumbnailVolumeDirty &&
       !_subVolumeDirty &&
       cvcapp.properties<int>(getObjectName("sync_viewers")))
      {
        //sync zoomed to thumbnail viewer once both have been
        //fully initialized
        syncViewers("thumbnail.orientation");
        syncViewers("thumbnail.position");
      }
  }
  
  void Viewers::ensureVolumeAvailability()
  {
    std::vector<std::string> volDataKeys = cvcapp.listProperty(getObjectName("volumes"));
    for (const auto& volDataKey : volDataKeys)
      if(!cvcapp.isData<VolMagick::Volume>(volDataKey))
        cvcapp.data(volDataKey,VolMagick::Volume());
  }

  void Viewers::ensureSubVolumeAvailability()
  {
    std::vector<std::string> volDataKeys = cvcapp.listProperty(getObjectName("subvolumes"));
    for (const auto& volDataKey : volDataKeys)
      if(!cvcapp.isData<VolMagick::Volume>(volDataKey))
        cvcapp.data(volDataKey,VolMagick::Volume());
  }

  void Viewers::customEvent(QEvent *event)
  {
    CVCEvent *mwe = dynamic_cast<CVCEvent*>(event);
    if(!mwe) return;

    if(mwe->name == "handlePropertiesChanged")
      handlePropertiesChanged(boost::any_cast<std::string>(mwe->data));
    if(mwe->name == "handleDataChanged")
      handleDataChanged(boost::any_cast<std::string>(mwe->data));
    if(mwe->name == "loadThumbnail")
      loadThumbnail();
    if(mwe->name == "loadSubVolume")
      loadSubVolume();
    if(mwe->name == "saveImage")
      saveImage(boost::any_cast<int>(mwe->data));
    if(mwe->name == "syncViewers")
      syncViewers(boost::any_cast<std::string>(mwe->data));
  }

  void Viewers::propertiesChanged(const std::string& key)
  {
    QCoreApplication::postEvent(this,
      new CVCEvent("handlePropertiesChanged",key)
    );
  }

  void Viewers::dataChanged(const std::string& key)
  {
    QCoreApplication::postEvent(this,
      new CVCEvent("handleDataChanged",key)
    );
  }

  void Viewers::handlePropertiesChanged(const std::string& key)
  {
    using namespace std;
    using namespace boost;
    using namespace boost::algorithm;

    static log4cplus::Logger logger = FUNCTION_LOGGER;

    //if everything changed, iterate across all keys
    if(key == "all")
      {
        PropertyMap map = cvcapp.properties();
        for (const auto& val : map) {
            assert(val.first!="all");
            handlePropertiesChanged(val.first);
          }
        return;
      }
    
    vector<string> key_idents;
    split(key_idents, key, is_any_of("."));

    // arand, 2/15/2012
    // kind of a hack for automatically loading a vinay file from the properties map
    if(key_idents.size() == 1 && key_idents[0] == "transfer_function_fullPath")
     {
        updateColorTable();
      }

    //check if the key was meant for this object
    if(key_idents.size() == 2 && key_idents[0] == getObjectName())
      {
	cvcapp.log(4,str(boost::format("%s :: object %s, property %s\n")
			 % BOOST_CURRENT_FUNCTION
			 % key_idents[0]
			 % key_idents[1]).c_str());
        try
          {
            if(key_idents[1] == "vfi"                 ||
	       key_idents[1] == "single_volume_index" ||
	       key_idents[1] == "red_volume_index"    ||
	       key_idents[1] == "green_volume_index"  ||
	       key_idents[1] == "blue_volume_index"   ||
	       key_idents[1] == "alpha_volume_index"  ||
	       key_idents[1] == "volume_timestep"     ||
	       key_idents[1] == "rendering_mode"      ||
	       key_idents[1] == "shaded_rendering_enabled")
              {
                if(key_idents[1] == "rendering_mode")
                  {
                    //Change the rendering mode in the viewers
                    cvcapp.properties("thumbnail.rendering_mode",
                                      cvcapp.properties(getObjectName("rendering_mode")));
                    cvcapp.properties("zoomed.rendering_mode",
                                      cvcapp.properties(getObjectName("rendering_mode")));

                    //Set the volume names to use
                    if(cvcapp.properties(getObjectName("rendering_mode"))=="colormapped")
                      {
                        cvcapp.properties(getObjectName("volumes"),"thumbnail_volume");
                        cvcapp.properties(getObjectName("subvolumes"),"zoomed_volume");

                        //Attempt to make the the single_volume_index less than
                        //the number of variables in the volume file if it isn't already.
                        if(cvcapp.hasProperty(getObjectName("vfi")) &&
                           cvcapp.isData<VolMagick::VolumeFileInfo>(
                             cvcapp.properties(getObjectName("vfi"))
                           ))
                          {
                            VolMagick::VolumeFileInfo vfi =
                              cvcapp.data<VolMagick::VolumeFileInfo>(
                                cvcapp.properties(getObjectName("vfi"))
                              );
                            unsigned int single_volume_index =
                              cvcapp.properties<unsigned int>(getObjectName("single_volume_index"));
                            if(single_volume_index >= vfi.numVariables() && vfi.numVariables()>0)
                              single_volume_index=vfi.numVariables()-1;
                            cvcapp.properties(getObjectName("single_volume_index"),single_volume_index);
                          }
                      }
                    else if(cvcapp.properties(getObjectName("rendering_mode")) == "rgba" ||
                            cvcapp.properties(getObjectName("rendering_mode")) == "RGBA")
                      {
                        string separators = cvcapp.properties("system.list_separators");
                        //use the first separator in the list if not empty. Else use ,
                        string separator;
                        if(separators.empty())
                          separator = ",";
                        else
                          separator.assign(separators,0,1);

                        cvcapp.properties(getObjectName("volumes"),
                                          "thumbnail_volume_red"+separator+
                                          "thumbnail_volume_green"+separator+
                                          "thumbnail_volume_blue"+separator+
                                          "thumbnail_volume_alpha");
                        cvcapp.properties(getObjectName("subvolumes"),
                                          "zoomed_volume_red"+separator+
                                          "zoomed_volume_green"+separator+
                                          "zoomed_volume_blue"+separator+
                                          "zoomed_volume_alpha");
                        
                        //For now, lets reset the indices so we get what the user expects
                        cvcapp.properties(getObjectName("red_volume_index"),0);
                        cvcapp.properties(getObjectName("green_volume_index"),1);
                        cvcapp.properties(getObjectName("blue_volume_index"),2);
                        cvcapp.properties(getObjectName("alpha_volume_index"),3);

                        //Attempt to make the the rgba indices less than
                        //the number of variables in the volume file if they aren't already.
                        if(cvcapp.hasProperty(getObjectName("vfi")) &&
                           cvcapp.isData<VolMagick::VolumeFileInfo>(
                             cvcapp.properties(getObjectName("vfi"))
                           ))
                          {
                            VolMagick::VolumeFileInfo vfi =
                              cvcapp.data<VolMagick::VolumeFileInfo>(
                                cvcapp.properties(getObjectName("vfi"))
                              );
                            unsigned int red_volume_index =
                              cvcapp.properties<unsigned int>(getObjectName("red_volume_index"));
                            if(red_volume_index >= vfi.numVariables() && vfi.numVariables()>0)
                              red_volume_index=vfi.numVariables()-1;
                            cvcapp.properties(getObjectName("red_volume_index"),red_volume_index);
                            unsigned int green_volume_index =
                              cvcapp.properties<unsigned int>(getObjectName("green_volume_index"));
                            if(green_volume_index >= vfi.numVariables() && vfi.numVariables()>0)
                              green_volume_index=vfi.numVariables()-1;
                            cvcapp.properties(getObjectName("green_volume_index"),green_volume_index);
                            unsigned int blue_volume_index =
                              cvcapp.properties<unsigned int>(getObjectName("blue_volume_index"));
                            if(blue_volume_index >= vfi.numVariables() && vfi.numVariables()>0)
                              blue_volume_index=vfi.numVariables()-1;
                            cvcapp.properties(getObjectName("blue_volume_index"),blue_volume_index);
                            unsigned int alpha_volume_index =
                              cvcapp.properties<unsigned int>(getObjectName("alpha_volume_index"));
                            if(alpha_volume_index >= vfi.numVariables() && vfi.numVariables()>0)
                              alpha_volume_index=vfi.numVariables()-1;
                            cvcapp.properties(getObjectName("alpha_volume_index"),alpha_volume_index);
                          }
                      }
                    else
                      {
                        cvcapp.log(1,str(boost::format("%s :: Invalid rendering_mode %s\n")
                                         % BOOST_CURRENT_FUNCTION
                                         % cvcapp.properties(getObjectName("rendering_mode"))));
                      }
                  }

		//When the thumbnail volume is updated, the subvolume 
		//bounding box reset in the thumbnail viewer causes the
		//subvolume to get updated later
                // cvcapp.log(2,str(boost::format("%s :: %s changed, loading thumbnail\n")
                //                  % BOOST_CURRENT_FUNCTION
                //                  % key));
                LOG4CPLUS_TRACE(logger, key << " changed, loading thumbnail");
		markThumbnailDirty();
	      }
	    else if(key_idents[1] == "subvolume_bounding_box")
	      {
                // cvcapp.log(2,str(boost::format("%s :: %s changed, loading subvolume\n")
                //                  % BOOST_CURRENT_FUNCTION
                //                  % key));
                LOG4CPLUS_TRACE(logger, key << " changed, loading subvolume");
		markSubVolumeDirty();
	      }
	    else if(key_idents[1] == "colortable_interactive_updates")
	      {
		bool flag = cvcapp.properties(key) == "true";
		_colorTable->interactiveUpdates(flag);
	      }
            else if(key_idents[1] == "volumes")
              {
                //If any volumes are missing from the datamap, add some for now until they get
                //updated.  
                ensureVolumeAvailability();

                //Make sure the viewers have the list of volumes
                cvcapp.properties("thumbnail.volumes",
                                  cvcapp.properties(getObjectName("volumes")));
                
                //Set the isocontouring volumes to the last volume in the list of volumes.
                //This way, for colormapped it wil pick the single volume in the list
                //and for rgba it will pick the last volume in the list, the alpha volume.
                std::vector<std::string> volumeKeys = cvcapp.listProperty(getObjectName("volumes"));
                cvcapp.properties(getObjectName("thumbnail_isocontouring_volume"),
                                  !volumeKeys.empty() ? 
                                  *prior(volumeKeys.end()) : 
                                  "thumbnail_volume");

                //trigger a re-load of thumbnail data
                // cvcapp.log(2,str(boost::format("%s :: %s changed, loading thumbnail\n")
                //                  % BOOST_CURRENT_FUNCTION
                //                  % key));
                LOG4CPLUS_TRACE(logger, key << " changed, loading thumbnail");
                markThumbnailDirty();
              }
            else if(key_idents[1] == "subvolumes")
              {
                //If any volumes are missing from the datamap, add some for now until they get
                //updated.  
                ensureSubVolumeAvailability();

                //Make sure the viewers have the list of volumes
                cvcapp.properties("zoomed.volumes",
                                  cvcapp.properties(getObjectName("subvolumes")));

                //Set the isocontouring volumes to the last volume in the list of volumes.
                //This way, for colormapped it wil pick the single volume in the list
                //and for rgba it will pick the last volume in the list, the alpha volume.
                std::vector<std::string> subvolumeKeys = cvcapp.listProperty(getObjectName("subvolumes"));
                cvcapp.properties(getObjectName("zoomed_isocontouring_volume"),
                                  !subvolumeKeys.empty() ?
                                  *prior(subvolumeKeys.end()) : 
                                  "zoomed_volume");

                //trigger a re-load of subvolume data
                // cvcapp.log(2,str(boost::format("%s :: %s changed, loading subvolume\n")
                //                  % BOOST_CURRENT_FUNCTION
                //                  % key));
                LOG4CPLUS_TRACE(logger, key << " changed, loading subvolume");
                markSubVolumeDirty();
              }
            else if(key_idents[1] == "thumbnail_isocontouring_volume" ||
                    key_idents[1] == "thumbnail_isocontour_geometry")
              {
                string thumbnail_isocontouring_thread =
                  cvcapp.properties(getObjectName("thumbnail_isocontouring_thread"));
                if(!cvcapp.hasThread(thumbnail_isocontouring_thread))
                  IsocontouringThread::start(
                    thumbnail_isocontouring_thread,
                    cvcapp.properties(getObjectName("color_table_info")),
                    cvcapp.properties(getObjectName("thumbnail_isocontouring_volume")),
                    cvcapp.properties(getObjectName("thumbnail_isocontour_geometry"))
                  );
              }
            else if(key_idents[1] == "zoomed_isocontouring_volume" ||
                    key_idents[1] == "zoomed_isocontour_geometry")
              {
                string zoomed_isocontouring_thread =
                  cvcapp.properties(getObjectName("zoomed_isocontouring_thread"));
                if(!cvcapp.hasThread(zoomed_isocontouring_thread))	    
                  IsocontouringThread::start(
                    zoomed_isocontouring_thread,
                    cvcapp.properties(getObjectName("color_table_info")),
                    cvcapp.properties(getObjectName("zoomed_isocontouring_volume")),
                    cvcapp.properties(getObjectName("zoomed_isocontour_geometry"))
                  );
              }
            else if(key_idents[1] == "color_table_opacity_cubed")
              {
                updateColorTable();
              }
	  }
	catch(std::exception& e)
          {
            string msg = str(boost::format("%s :: Error: %s")
                             % BOOST_CURRENT_FUNCTION
                             % e.what());
            // cvcapp.log(1,msg+"\n");
            LOG4CPLUS_ERROR(logger, e.what());
            QMessageBox::critical(this,"Error",
                                  QString("Error setting value for key %1 : %2")
                                  .arg(QString::fromStdString(key))
                                  .arg(QString::fromStdString(msg)));
          }
      }

    //update the volumes if there isn't a vfi and the thumbnail or zoomed geometry list changed
    if((key == "zoomed.geometries" || key == "thumbnail.geometries") &&
       cvcapp.properties(getObjectName("vfi")) == "none")
      {
        // cvcapp.log(2,str(boost::format("%s :: %s changed, loading thumbnail\n")
        //                  % BOOST_CURRENT_FUNCTION
        //                  % key));
        LOG4CPLUS_TRACE(logger, key << " changed, loading thumbnail");
        markThumbnailDirty();
      }

    //update the volumes if the hdf5 hierarchy thread reports that it is finished computing
    //something that we have loaded.
    if(key == "volmagick.hdf5_io.buildhierarchy.latest")
      {
        string hierarchy_filename = cvcapp.properties("volmagick.hdf5_io.buildhierarchy.latest");
        vector<string> hier_parts;
        split(hier_parts, hierarchy_filename, is_any_of("|"));
        if(!hier_parts.empty())
          hierarchy_filename = hier_parts[0];

        string current_vfi = cvcapp.properties(getObjectName("vfi"));
        if(cvcapp.isData<VolMagick::VolumeFileInfo>(current_vfi))
          {
            VolMagick::VolumeFileInfo vfi =
              cvcapp.data<VolMagick::VolumeFileInfo>(current_vfi);

            string current_vfi_filename = vfi.filename();
            LOG4CPLUS_TRACE(logger, "vfi = " << current_vfi_filename << 
                            ", hierarchy_filename = " << hierarchy_filename);
            // cvcapp.log(4,str(boost::format("%s :: vfi = %s, hierarchy_filename = %s\n")
            //                  % BOOST_CURRENT_FUNCTION
            //                  % current_vfi_filename
            //                  % hierarchy_filename));

            if(current_vfi_filename.find(hierarchy_filename) != std::string::npos)
              {
                cvcapp.log(2,str(boost::format("%s :: new hierarchy info, reload\n")
                                 % BOOST_CURRENT_FUNCTION));
                markThumbnailDirty();
              }
          }
      }

    //sync viewers
    if(cvcapp.properties<int>(getObjectName("sync_viewers")) &&
       _thumbnailPostInitFinished &&
       _subvolumePostInitFinished)
      {
        //QCoreApplication::postEvent(this,new CVCEvent("syncViewers",key));
        syncViewers(key);
      }
  } // Viewers::handlePropertiesChanged

  // 09/17/2011 - transfix - using loadThumbnail and loadSubVolume
  void Viewers::handleDataChanged(const std::string& key)
  {
    using namespace std;
    using namespace boost;
    using namespace boost::algorithm;

    static log4cplus::Logger logger = FUNCTION_LOGGER;

    //if everything changed, iterate across all keys
    if(key == "all")
      {
        DataMap map = cvcapp.data();
        for (const auto& val : map) {
            assert(val.first!="all");
            handleDataChanged(val.first);
          }
        return;
      }

    try
      {
	//if the subvolume bounding box data we're watching has changed,
	//we need to re-extract
	if(cvcapp.properties(getObjectName("subvolume_bounding_box")) == key)
	  {
            // cvcapp.log(2,str(boost::format("%s :: %s changed, loading subvolume\n")
            //                  % BOOST_CURRENT_FUNCTION
            //                  % key));
            LOG4CPLUS_TRACE(logger, key << " changed, loading subvolume");
	    markSubVolumeDirty();
	  }

	//If the vfi we're watching has changed, reload the thumbnail.
        //this will trigger a reload of the subvolume as well.
	if(cvcapp.properties(getObjectName("vfi")) == key)
	  {
            //clear the isocache if we're loading a new volume file
            IsocontouringThread::isocache(IsocontouringThread::isocache_t());

            // cvcapp.log(2,str(boost::format("%s :: %s changed, loading thumbnail\n")
            //                  % BOOST_CURRENT_FUNCTION
            //                  % key));
            LOG4CPLUS_TRACE(logger, key << " changed, loading thumbnail");

	    markThumbnailDirty();
	  }

        if(cvcapp.properties(getObjectName("color_table_info")) == key)
          {
            using namespace CVCColorTable;

            ColorTable::color_table_info cti = 
              cvcapp.data<ColorTable::color_table_info>(key);
            ColorTable::isocontour_nodes nodes = cti.isocontourNodes();

            //Return right away if the nodes haven't actually changed.
            //There is nothing to do here in that case.
            if(_oldNodes == nodes) return;
            _oldNodes = nodes;

            //Check to make sure that the min & max of thumb and zoomed volumes agree
            VolMagick::Volume thumb_vol =
              cvcapp.data<VolMagick::Volume>(cvcapp.properties(getObjectName("thumbnail_isocontouring_volume")));
            VolMagick::Volume zoom_vol =
              cvcapp.data<VolMagick::Volume>(cvcapp.properties(getObjectName("zoomed_isocontouring_volume")));
            if(thumb_vol.min() == zoom_vol.min() &&
               thumb_vol.max() == zoom_vol.max())
              {
                VolMagick::Volume& vol = thumb_vol;

                //Clear out data in the isocache that we don't need
                set<double> isovals;
                for (const auto& node : nodes) {
                    double isoval = vol.min()+(vol.max()-vol.min())*node.position;
                    isovals.insert(isoval);
                  }
                IsocontouringThread::isocache_t isoc = IsocontouringThread::isocache();
                set<double> isocache_vals;
                for (const auto& iv : isoc)
                  isocache_vals.insert(iv.first);
                set<double> clear_vals;
                set_difference(isocache_vals.begin(), isocache_vals.end(),
                               isovals.begin(), isovals.end(),
                               inserter(clear_vals, clear_vals.begin()));
                for (const auto& val : clear_vals) {
                    cvcapp.log(5,str(format("Erasing isoval %1%\n")%val));
                    isoc.erase(val);
                  }
                IsocontouringThread::isocache(isoc);
              }
            else
              {
                cvcapp.log(1,str(format("%s :: Warning: thumbnail and zoomed volumes' min/maxes do not agree\n")
                                 % BOOST_CURRENT_FUNCTION));
                //clear the isocache if we're loading a new volume file
                IsocontouringThread::isocache(IsocontouringThread::isocache_t());
              }

            std::string thumbnail_isocontouring_thread = 
              cvcapp.properties(getObjectName("thumbnail_isocontouring_thread"));
            std::string zoomed_isocontouring_thread = 
              cvcapp.properties(getObjectName("zoomed_isocontouring_thread"));
	    if(!cvcapp.hasThread(thumbnail_isocontouring_thread))
	      IsocontouringThread::start(
                thumbnail_isocontouring_thread,
		key,
		cvcapp.properties(getObjectName("thumbnail_isocontouring_volume")),
		cvcapp.properties(getObjectName("thumbnail_isocontour_geometry"))
	      );

	    if(!cvcapp.hasThread(zoomed_isocontouring_thread))	    
	      IsocontouringThread::start(
                zoomed_isocontouring_thread,
		key,
		cvcapp.properties(getObjectName("zoomed_isocontouring_volume")),
		cvcapp.properties(getObjectName("zoomed_isocontour_geometry"))
	      );
          }

        //start the isocontouring thread if the thumbnail_isocontouring_volume data has changed
        if(cvcapp.properties(getObjectName("thumbnail_isocontouring_volume")) == key)
          {
            string thumbnail_isocontouring_thread = 
              cvcapp.properties(getObjectName("thumbnail_isocontouring_thread"));
            if(!cvcapp.hasThread(thumbnail_isocontouring_thread)) 
              IsocontouringThread::start(
                thumbnail_isocontouring_thread,
                cvcapp.properties(getObjectName("color_table_info")),
                cvcapp.properties(getObjectName("thumbnail_isocontouring_volume")),
                cvcapp.properties(getObjectName("thumbnail_isocontour_geometry"))
              );
          }

        //start the isocontouring thread if the zoomed_isocontouring_volume data has changed
        if(cvcapp.properties(getObjectName("zoomed_isocontouring_volume")) == key)
          {
            string zoomed_isocontouring_thread = 
              cvcapp.properties(getObjectName("zoomed_isocontouring_thread"));
            if(!cvcapp.hasThread(zoomed_isocontouring_thread)) 
              IsocontouringThread::start(
                zoomed_isocontouring_thread,
                cvcapp.properties(getObjectName("color_table_info")),
                cvcapp.properties(getObjectName("zoomed_isocontouring_volume")),
                cvcapp.properties(getObjectName("zoomed_isocontour_geometry"))
              );
          }


        //if one of the subvolumes changed, set the contour spectrum to it
        //TODO: add a property to set the contour volume explicitly
        //TODO: also do the contour spectrum/tree calculations in a separate thread
        vector<string> volDataKeys = cvcapp.listProperty(getObjectName("subvolumes"));
        for (const auto& volDataKey : volDataKeys)
          if(volDataKey == key)
            {
              _colorTable->setContourVolume(cvcapp.data<VolMagick::Volume>(volDataKey));
            }

      }
    catch(std::exception& e)
      {
	string msg = str(format("%s :: Error: %s")
			 % BOOST_CURRENT_FUNCTION
			 % e.what());
	cvcapp.log(1,msg+"\n");
	QMessageBox::critical(this,"Error",
			      QString("Error handling data key %1 : %2")
			      .arg(QString::fromStdString(key))
			      .arg(QString::fromStdString(msg)));
      }
  }

  void Viewers::updateColorTable()
  {

    if ((_colorTable->transfer_function_filename).size() > 0) {
      cvcapp.properties("transfer_function_fullPath",_colorTable->transfer_function_filename); 
    } else if (cvcapp.hasProperty("transfer_function_fullPath")) {
      _colorTable->read(cvcapp.properties("transfer_function_fullPath"));
    }

    std::string opacity_cubed =
      cvcapp.properties(getObjectName("color_table_opacity_cubed"));
    _colorTable->opacityCubed(opacity_cubed == "false" ? false : true);

#ifdef VOLUMEROVER2_FLOAT_COLORTABLE
    cvcapp.data(cvcapp.properties(getObjectName("transfer_function")),
		_colorTable->getFloatTable());
#else
    cvcapp.data(cvcapp.properties(getObjectName("transfer_function")),
		_colorTable->getTable());
#endif
    cvcapp.data(cvcapp.properties(getObjectName("color_table_info")),
                _colorTable->info());

#ifdef USING_VOLUMEGRIDROVER
    if( _volumeGridRoverPtr )
    {
      unsigned char *byte_map = _colorTable->getCharTable();
      if( byte_map ) 
        _volumeGridRoverPtr->setTransferFunction( _colorTable->getCharTable() );
    }
#endif
  }
  
  // arand, 6-8-2011
  // arguments are (1) thumbnail or zoomed, (2) filename, (3) filetype
  void Viewers::saveImage(int panel) {

    VolumeViewer * vv = _subvolumeViewer;
    if (panel==1) {
      vv = _thumbnailViewer;
    }

    vv->saveSnapshot(false,false);
    // arand: old style... using QGLViewer is better
    //QImage saveImg = vv->grabFrameBuffer();
    //saveImg.save(sii.filename.c_str(), sii.filetype.c_str());
  }

  //for syncronising their orientations
  // 10/10/2011 -- transfix -- initial implementation
  void Viewers::syncViewers(const std::string& key)
  {
    if(cvcapp.properties<int>(getObjectName("sync_viewers")))
      {
        if(key == "thumbnail.orientation" &&
           cvcapp.properties("zoomed.orientation")!=cvcapp.properties("thumbnail.orientation"))
          cvcapp.properties("zoomed.orientation",
                            cvcapp.properties("thumbnail.orientation"));

        if(key == "zoomed.orientation" &&
           cvcapp.properties("zoomed.orientation")!=cvcapp.properties("thumbnail.orientation"))
          cvcapp.properties("thumbnail.orientation",
                            cvcapp.properties("zoomed.orientation"));

        if(key == "thumbnail.position" &&
           cvcapp.properties("zoomed.position")!=cvcapp.properties("thumbnail.position"))
          cvcapp.properties("zoomed.position",
                            cvcapp.properties("thumbnail.position"));

        if(key == "zoomed.position" &&
           cvcapp.properties("zoomed.position")!=cvcapp.properties("thumbnail.position"))
          cvcapp.properties("thumbnail.position",
                            cvcapp.properties("zoomed.position"));

#if 0
        //only attempt to fit the subvolume to screen if the thumbnail viewer changed
        if(cvcapp.properties<int>(getObjectName("subvolume_fit_screen")) &&
           (key == "thumbnail.orientation" || key == "thumbnail.position"))
          {
            qglviewer::Camera cam(*_thumbnailViewer->camera());
            cam.setFOVToFitScene();
            cvcapp.properties("zoomed.fov",cam.fieldOfView());
          }
        else
          {
            //use the thumbnail's fov if not fitting to screen
            cvcapp.properties("zoomed.fov",
                              cvcapp.properties("thumbnail.fov"));
          }
#endif
      }
  }
}
