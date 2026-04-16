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

#include <algorithm>
#include <boost/array.hpp>
#include <boost/scoped_array.hpp>
#include <boost/shared_array.hpp>
#include <boost/current_function.hpp>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/utility.hpp>
#include <boost/bind/bind.hpp>
#include <cmath>
#include <vector>
#include <VolumeRover2/VolumeViewer.h>

using namespace boost::placeholders;
#include <CVC/CVCEvent.h>
#include <CVC/State.h>
#include <log4cplus/logger.h>
#include <CVC/log4cplus_compat.h>
#include <log4cplus/loggingmacros.h>

#include <QGLViewer/quaternion.h>

#ifdef USING_TILING
#include <cvcraw_geometry/contours.h>
#endif

#include <QApplication>
#include <QProgressDialog>
#include <QMouseEvent>
#include <QWheelEvent>
#include <QMessageBox>
#include <QTimer>

#include <iostream>
#include <sstream>
#include <VolumeRover2/MultiTileServerDialog.h>

//#define DISABLE_VBO_SUPPORT
//#define DEBUG_VOLUMEVIEWER

namespace CVC_NAMESPACE
{
  //=========================================================================
  // XmlRpc thread routines
  void* XmlRpcThread::clientThreadEntryPoint( void * pthis ) {
     XmlRpcThread *ptr = (XmlRpcThread*) pthis;
     ptr->run();
  }

  //#define FAIL_RETURN

  void XmlRpcThread::run(void) {
     //send initialization information
     XmlRpc::XmlRpcValue volArg, resArg, result;
     int ntilex, ntiley;
     m_VolumeViewerPtr->multiTileServerNtiles( &ntilex, &ntiley );
     resArg[0] = ntilex;
     resArg[1] = ntiley;
    
     int trial = 0, maxTrial = 10;
     while( !m_client->execute( "updateResolution", resArg, result ) && (trial++ < maxTrial) );

     if( trial >= maxTrial ) {
        fprintf( stderr, "[client %d] Error calling updateResolution\n", m_id );
        #ifdef FAIL_RETURN
        return;
        #endif
     }
     else
        fprintf( stderr, "[client %d] updateResolution done with trial[%d]\n", m_id, trial );

     #ifndef SYNC_THUMBNAIL_WITH_MULTITILESERVER
     XmlRpc::XmlRpcValue bboxArg;
     double bbox[6];
     bbox[0] = m_VolumeViewerPtr->boundingBox().XMax();
     bbox[1] = m_VolumeViewerPtr->boundingBox().YMax();
     bbox[2] = m_VolumeViewerPtr->boundingBox().ZMax();
     bbox[3] = m_VolumeViewerPtr->boundingBox().XMin();
     bbox[4] = m_VolumeViewerPtr->boundingBox().YMin();
     bbox[5] = m_VolumeViewerPtr->boundingBox().ZMin();
     for(int i=0; i<6; i++) {
        bboxArg[i] = bbox[i];
     }
     
     trial = 0;
     while( !m_client->execute( "updateSubBbox", bboxArg, result ) && (trial++ < maxTrial) );

     if( trial >= maxTrial ) {
        fprintf( stderr, "[client %d] Error calling updateSubBbox\n", m_id );
        #ifdef FAIL_RETURN
        return; 
        #endif
     }
     else
        fprintf( stderr, "[client %d] updateSubBbox done with trial[%d]\n", m_id, trial );
 
     #endif

     std::string volFile = m_VolumeViewerPtr->multiTileServerVolumeFile();
     volArg[0] = volFile;
     trial = 0;
     while( !m_client->execute( "updateVolume", volArg, result ) && (trial++ < maxTrial) );

     if( trial >= maxTrial ) {
        fprintf( stderr, "[client %d] Error calling updateVolume\n", m_id );
        #ifdef FAIL_RETURN
        return;
        #endif
     }
     else
        fprintf( stderr, "[client %d] updateVolume done with trial[%d]\n", m_id, trial );

     while(1) {
        if( m_VolumeViewerPtr->needToUpdateMultiTileServer( m_id ) ) {

            if(m_VolumeViewerPtr->syncCameraWithMultiTileServer()) {
               XmlRpc::XmlRpcValue arg, result;
               double mv[16];
               m_VolumeViewerPtr->camera()->getModelViewMatrix( mv );
               for( int i =0; i < 16; i++ ) {
                  arg[i] = mv[i];
                  //fprintf( stderr, "%lf\n", mv[i] );
               }

               //double p[16];
               //camera()->getProjectionMatrix( p );

               if( !m_client->execute( "updateView", arg, result ) )
                  fprintf( stderr, "[client %d] Error calling updateView\n", m_id );
            }
            if(m_VolumeViewerPtr->syncTransferFuncWithMultiTileServer()) {
               XmlRpc::XmlRpcValue arg, result;
               const unsigned char *table = m_VolumeViewerPtr->colorTable();
               for(int i = 0; i < 1024; i++ ) {
                  int val = table[i];
                  arg[i] = val;
               }

               if( !m_client->execute( "updateTransFunc", arg, result ) )
                  fprintf( stderr, "[client %d] Error calling updateTransFunc\n", m_id );
            }
            if(m_VolumeViewerPtr->syncShadedRenderWithMultiTileServer()) {
               XmlRpc::XmlRpcValue arg, result;
               int cflag = m_VolumeViewerPtr->shadedRenderingEnabled();
               if( m_savedShadedRenderFlag != cflag ) {
                  arg[0] = cflag;

                  if( !m_client->execute( "updateSRFlag", arg, result ) )
                     fprintf( stderr, "[client %d] Error calling updateSRFlag\n", m_id );
                  m_savedShadedRenderFlag = cflag;
               }
            }
            if(m_VolumeViewerPtr->syncRenderModeWithMultiTileServer()) {
               XmlRpc::XmlRpcValue arg, result;
               int cflag = (int)m_VolumeViewerPtr->volumeRenderingType();
               if( m_savedRenderModeFlag != cflag ) {
                  arg[0] = cflag;

                  if( !m_client->execute( "updateRMFlag", arg, result ) )
                     fprintf( stderr, "[client %d] Error calling updateRMFlag\n", m_id );
                  m_savedRenderModeFlag = cflag;
               }
            }

            m_VolumeViewerPtr->needToUpdateMultiTileServer( m_id, false );
        }
        if( m_VolumeViewerPtr->terminateMultiTileClient( m_id ) ) {
           m_VolumeViewerPtr->terminateMultiTileClient( m_id, false );
           fprintf( stderr, "[client %d] thread will terminated\n", m_id );
           break;
        }
     }
  }

  int XmlRpcThread::start( VolumeViewer *_ptr, std::string _host, int _port, int _id ) {

     m_VolumeViewerPtr = _ptr;
     if( _host != m_hostname || m_port != _port ) {
        if( m_client ) delete m_client;
        m_hostname = _host;
        m_port = _port;
        m_client = new XmlRpc::XmlRpcClient(m_hostname.data(), m_port);
     }
     m_id = _id;

     int code = pthread_create( &m_thread_id, NULL, XmlRpcThread::clientThreadEntryPoint, this );
     return code;
  }

  void XmlRpcThread::end(void) {
  }

  //=========================================================================
  void VolumeViewer::defaultConstructor()
  {
    //NO GL CALLS ALLOWED HERE!!!! do them in init()
    //I mean it, you will have a nasty crash in the GL driver otherwise! This is because
    //this function is called without the guarantee that the GL context yet exists.

    _hasBeenInitialized = false;

    _drawBoundingBox = false;

    _volumeRenderingType = ColorMapped;
        
    _drawable = false;
    _cmVolumeUploaded = false;
    _rgbaVolumeUploaded = false;
    _shadedRenderingEnabled = false;
    
    _cmVolume = VolMagick::Volume(VolMagick::Dimension(4,4,4),
				  UChar,
				  VolMagick::BoundingBox(-0.5,-0.5,-0.5,0.5,0.5,0.5));
    _rgbaVolumes = std::vector<VolMagick::Volume>(4,VolMagick::Volume(VolMagick::Dimension(4,4,4),
								      UChar,
								      VolMagick::BoundingBox(-0.5,-0.5,-0.5,
											     0.5,0.5,0.5)));

    _drawSubVolumeSelector = false;
    _cmSubVolumeSelector = VolMagick::BoundingBox(-0.25,-0.25,-0.25,0.25,0.25,0.25);
    _rgbaSubVolumeSelector = VolMagick::BoundingBox(-0.25,-0.25,-0.25,0.25,0.25,0.25);

    _selectedObj = -1;

    _drawCornerAxis = false;

    _drawGeometry = true;
    _drawVolumes = true;
    _usingVBO = false;
    _vboUpdated = false;
    _geometries.clear();
    
    _clipGeometryToBoundingBox = true;
    _drawGeometryNormals = false;

    //Watch the global property map
    _propertiesConnection.disconnect();
    _propertiesConnection = 
      cvcapp.propertiesChanged.connect(
          MapChangeSignal::slot_type(
            &VolumeViewer::propertiesChanged, this, _1
          )
        );

    //now read the properties map and pick up any settings from there
    handlePropertiesChanged("all");

    //watch the global data map
    _dataConnection.disconnect();
    _dataConnection =
      cvcapp.dataChanged.connect(
	  MapChangeSignal::slot_type(
	    &VolumeViewer::dataChanged, this, _1
	  )
	);

    handleDataChanged("all");

    //watch the thread map to print messages for the user
    _threadsConnection.disconnect();
    _threadsConnection = 
      cvcapp.threadsChanged.connect(
          MapChangeSignal::slot_type(
            &VolumeViewer::threadsChanged, this, _1
          )
        );

    //do initial copy so subvolume extraction can happen if needed
    copyBoundingBoxToDataMap();

    _timer = new QTimer(this);
    _timer->start(500);
    connect(_timer,SIGNAL(timeout()),SLOT(timeout()));

    _glUpdated = false;

    _syncCameraWithMultiTileServer = false;
    _syncTransferFuncWithMultiTileServer = false;
    _syncShadedRenderWithMultiTileServer = false;
    _interactiveSyncWithMultiTileServer = false;
#ifndef USE_XmlRpc
    m_socket = NULL;
#endif
    _multiTileServer_nServer = 0;
    _terminateMultiTileClient = NULL;
    _updateMultiTileServer = NULL;

    _multiTileServer_hostList = NULL;
    _multiTileServer_portList = NULL;

    // to test handling XmlRpc delay
    _mousePressed = false;

    _doDisplayMessage = false;

    //mouseMoved is just used to prevent too many hits to the property map when rotating
    _mouseMoved = false;

    _copyCameraOnNormalize = false;
    _normalizeOnVolume = true;
    _showEntireSceneOnNormalize = true;
    _doNormalizeScene = false;

    _stateConnection.disconnect();
    _stateConnection = cvcstate(stateName()).childChanged.connect(
          MapChangeSignal::slot_type(
            &VolumeViewer::stateChanged, this, _1
          )
        );

    //TODO: replace properties/data usage with cvcstate
    cvcstate(stateName("rendering_mode")).value("colormapped");
    cvcstate(stateName("shaded_rendering_enabled")).value("false");
    cvcstate(stateName("geometry_rendering_mode")).value("triangles");
    cvcstate(stateName("draw_bounding_box")).value("true");
    cvcstate(stateName("draw_subvolume_selector")).value("true");
    cvcstate(stateName("volume_rendering_quality")).value("0.5"); //[0.0,1.0]
    cvcstate(stateName("volume_rendering_near_plane")).value("0.0");
    cvcstate(stateName("projection_mode")).value("perspective");
    cvcstate(stateName("draw_corner_axis")).value("true");
    cvcstate(stateName("draw_geometry")).value("true");
    cvcstate(stateName("draw_volumes")).value("true");
    cvcstate(stateName("clip_geometry")).value("true");
    cvcstate(stateName("draw_geometry_normals")).value("false");
    cvcstate(stateName("geometry_line_width")).value("1.2");
    cvcstate(stateName("background_color")).value("#000000");

    cvcstate(stateName("syncCamera_with_multiTileServer")).value("false");
    cvcstate(stateName("syncTransferFunc_with_multiTileServer")).value("false");
    cvcstate(stateName("syncShadedRender_with_multiTileServer")).value("false");
    cvcstate(stateName("syncRenderMode_with_multiTileServer")).value("false");
    cvcstate(stateName("interactiveMode_with_multiTileServer")).value("false");
    cvcstate(stateName("syncMode_with_multiTileServer")).value("0");
    cvcstate(stateName("syncInitial_multiTileServer")).value("0");

    //stereo related properties
    cvcstate(stateName("io_distance")).value("0.062");
    cvcstate(stateName("physical_distance_to_screen")).value("2.0");
    cvcstate(stateName("physical_screen_width")).value("1.8");
    cvcstate(stateName("focus_distance")).value("1000.0");

    cvcstate(stateName("fov")).value(boost::lexical_cast<std::string>(M_PI/4.0f));

    _geomStateConnection.disconnect();
    _geomStateConnection = cvcstate("VolumeViewer::scene_geometry_t").childChanged.connect(
          MapChangeSignal::slot_type(
            &VolumeViewer::geomStateChanged, this, _1
          )
    );
  }

  VolumeViewer::~VolumeViewer()
  {
    _propertiesConnection.disconnect();
    _dataConnection.disconnect();
    _threadsConnection.disconnect();
    _stateConnection.disconnect();
    _geomStateConnection.disconnect();
  }

  void VolumeViewer::setDefaultColorTable()
  {
#ifdef VOLUMEROVER2_FLOAT_COLORTABLE
    float table[256*4];
    memset(table,0,sizeof(float)*256*4);

    //set the color table to greyscale for now
    for(int i=0; i<256; i++)
      {
	table[i*4+0] = float(i)/255.0f;
	table[i*4+1] = float(i)/255.0f;
	table[i*4+2] = float(i)/255.0f;
	table[i*4+3] = float(i)/255.0f;
      }
#else
    unsigned char table[256*4];
    memset(table,0,sizeof(unsigned char)*256*4);

    //set the color table to greyscale for now
    for(int i=0; i<256; i++)
      {
	table[i*4+0] = i;
	table[i*4+1] = i;
	table[i*4+2] = i;
	table[i*4+3] = i;
      }
#endif
    colorTable(table);
  }

  void VolumeViewer::colorMappedVolume(const VolMagick::Volume& v)
  {
    _cmVolume = v;

    
    _cmSubVolumeSelector = VolMagick::BoundingBox(v.XMin() + (v.XMax()-v.XMin())/4.0,
						  v.YMin() + (v.YMax()-v.YMin())/4.0,
						  v.ZMin() + (v.ZMax()-v.ZMin())/4.0,
						  v.XMax() - (v.XMax()-v.XMin())/4.0,
						  v.YMax() - (v.YMax()-v.YMin())/4.0,
						  v.ZMax() - (v.ZMax()-v.ZMin())/4.0);
    
    
    /*
    // arand: hack to get the entire volume in the subvolume
    _cmSubVolumeSelector = VolMagick::BoundingBox(v.XMin(),
						  v.YMin(),
						  v.ZMin(),
						  v.XMax(),
						  v.YMax(),
						  v.ZMax());
    */

    _cmVolumeUploaded = false;
    copyBoxToDatamapAndEmitChange();
    if(normalizeOnVolume()) normalizeScene();
    if(hasBeenInitialized()) scheduleUpdateGL();
  }

  // 10/09/2011 -- transfix -- forcing volumes to have the same dimension and bounding box
  void VolumeViewer::rgbaVolumes(const std::vector<VolMagick::Volume>& v)
  {
    static log4cplus::Logger logger = FUNCTION_LOGGER;

    std::vector<VolMagick::Volume> vols = v;
    if(vols.empty())
      {
        LOG4CPLUS_INFO(logger, "no rgba volumes, using default volume");
        // cvcapp.log(2,
        //   boost::str(
	//     boost::format("%s: no rgba volumes, using default volume\n")
	//     % BOOST_CURRENT_FUNCTION));
        vols.push_back(VolMagick::Volume());
      }

    //duplicate the last vol until we have the needed 4
    while(vols.size() < 4)
      vols.push_back(*boost::prior(vols.end()));

    //Make sure the volume dimensions and bounding boxes all match.  Force non-matching volume
    //dimensions to match the first volume.
    for (auto& vol : vols) {
        if(vol.dimension() != vols[0].dimension())
          {
            cvcapp.log(2,
              boost::str(
                boost::format("%1% :: WARNING: %2% dimension does not match the first listed, resizing\n")
                % BOOST_CURRENT_FUNCTION
                % vol.desc()));
            vol.resize(vols[0].dimension());
          }
        if(vol.boundingBox() != vols[0].boundingBox())
          {
            cvcapp.log(2,
              boost::str(
                boost::format("%1% :: WARNING: %2% bounding box does not match the first listed\n")
                % BOOST_CURRENT_FUNCTION
                % vol.desc()));
            vol.boundingBox(vols[0].boundingBox());
          }
      }

    assert(vols[0].dimension() == vols[1].dimension() &&
	   vols[0].dimension() == vols[2].dimension() &&
	   vols[0].dimension() == vols[3].dimension());
    assert(vols[0].boundingBox() == vols[1].boundingBox() &&
	   vols[0].boundingBox() == vols[2].boundingBox() &&
	   vols[0].boundingBox() == vols[3].boundingBox());
    
    _rgbaVolumes = vols;
    _rgbaSubVolumeSelector = VolMagick::BoundingBox(vols[0].XMin() + (vols[0].XMax()-vols[0].XMin())/4.0,
						    vols[0].YMin() + (vols[0].YMax()-vols[0].YMin())/4.0,
						    vols[0].ZMin() + (vols[0].ZMax()-vols[0].ZMin())/4.0,
						    vols[0].XMax() - (vols[0].XMax()-vols[0].XMin())/4.0,
						    vols[0].YMax() - (vols[0].YMax()-vols[0].YMin())/4.0,
						    vols[0].ZMax() - (vols[0].ZMax()-vols[0].ZMin())/4.0);
    _rgbaVolumeUploaded = false;
    copyBoxToDatamapAndEmitChange();
    if(normalizeOnVolume()) normalizeScene();
    if(hasBeenInitialized()) scheduleUpdateGL();
  }

  std::string VolumeViewer::getObjectName(const std::string& property) const
  {
    std::string objName = objectName().toStdString();
    if(objName.empty()) objName = "VolumeViewer";
    if(!property.empty()) objName = objName + "." + property;
    return objName;
  }

  void VolumeViewer::addGeometry(const cvcraw_geometry::cvcgeom_t& geom,
                                 const std::string& name,
                                 scene_geometry_t::render_mode_t mode,
                                 bool do_updategl)
  {   
    // arand: grab current rendering mode...
    using namespace std;
    if (cvcapp.hasProperty(getObjectName("geometry_rendering_mode"))) {
      string defaultmode = cvcapp.properties(getObjectName("geometry_rendering_mode"));
      if (mode == scene_geometry_t::TRIANGLES) {
	if (defaultmode.compare("wireframe")==0) {
	  mode = scene_geometry_t::TRIANGLE_WIREFRAME;
	} else if (defaultmode.compare("filled_wireframe")==0) {
	  mode = scene_geometry_t::TRIANGLE_FILLED_WIRE;
	} else if (defaultmode.compare("flat")==0) {
	  mode = scene_geometry_t::TRIANGLES_FLAT;
	} else if (defaultmode.compare("flat_wireframe")==0) {
	  mode = scene_geometry_t::TRIANGLE_FLAT_FILLED_WIRE;
	}
      }
      if (mode == scene_geometry_t::QUADS) {
	if (defaultmode.compare("wireframe")==0) {
	  mode = scene_geometry_t::QUAD_WIREFRAME;
	} else if (defaultmode.compare("filled_wireframe")==0) {
	  mode = scene_geometry_t::QUAD_FILLED_WIRE;
	}
      }
    }    

    //calculate normals if we don't have any
    if(geom.normals().size() != geom.points().size())
      {
        cvcgeom_t newgeom(geom);
        newgeom.calculate_surf_normals();
        scene_geometry_t newscene_geom(newgeom,name,mode);
        _geometries[name] = scene_geometry_ptr(new scene_geometry_t(newscene_geom));
      }
    else
      {
        _geometries[name] =
          scene_geometry_ptr(new scene_geometry_t(geom,name,mode));
      }

    _vboUpdated = false;
    _allocatedBuffersToReInit.push_back(name);

    //override with property settings if available
    scene_geometry_ptr scene_geom_ptr = _geometries[name];
    if(cvcapp.hasProperty(name+".render_mode"))
      scene_geom_ptr->render_mode = 
        scene_geometry_t::render_mode_enum(cvcapp.properties(name+".render_mode"));
    if(cvcapp.hasProperty(name+".visible"))
      scene_geom_ptr->visible = cvcapp.properties(name+".visible") != "false";

    //generate the render_geometry for this geometry if we're 
    //rendering TETRA or HEXA - Joe, 4-19-2011
    if(scene_geom_ptr->render_mode == scene_geometry_t::TETRA ||
       scene_geom_ptr->render_mode == scene_geometry_t::HEXA)
      {
        cvcapp.log(3,str(boost::format("%s :: generating wire interior for geometry '%s'\n")
                         % BOOST_CURRENT_FUNCTION
                         % name));
        scene_geom_ptr->render_geometry =
          scene_geom_ptr->geometry.generate_wire_interior();
      }

    if(do_updategl) scheduleUpdateGL();
  }

  void VolumeViewer::drawBoundingBox(bool draw)
  {
    _drawBoundingBox = draw;
    if(hasBeenInitialized()) scheduleUpdateGL();
  }

  void VolumeViewer::volumeRenderingType(VolumeRenderingType t)
  {
    if(hasBeenInitialized()) makeCurrent();
    _volumeRenderingType = t;
    //Since internally there are 2 subvol bounding boxes, one for each
    //rendering type, we need to make sure that we have the currently
    //relevant one copied to the map.
    copyBoundingBoxToDataMap();
    if(normalizeOnVolume()) normalizeScene();
    if(hasBeenInitialized()) scheduleUpdateGL();

    if( interactiveSyncWithMultiTileServer() && syncRenderModeWithMultiTileServer() )
       #ifdef USE_XmlRpc
       needToUpdateMultiTileServer( -1, true );
       #else
       syncRenderModeFlag();
       #endif
  }

#ifdef VOLUMEROVER2_FLOAT_COLORTABLE
  void VolumeViewer::colorTable(const float *table)
  {
    memcpy(_colorTable, table, sizeof(float)*256*4);
    uploadColorTable();
    if(hasBeenInitialized()) scheduleUpdateGL();
  }
#else
  void VolumeViewer::colorTable(const unsigned char *table)
  {
    memcpy(_colorTable, table, sizeof(unsigned char)*256*4);
    uploadColorTable();
    if(hasBeenInitialized()) scheduleUpdateGL();
  }
#endif

  void VolumeViewer::drawSubVolumeSelector(bool show)
  {
    _drawSubVolumeSelector = show;
    if(hasBeenInitialized()) scheduleUpdateGL();
  }

  void VolumeViewer::quality(double q)
  {
    if(hasBeenInitialized()) makeCurrent();
    _cmRenderer.setQuality(q);
    _rgbaRenderer.setQuality(q);
    uploadColorTable();
    if(hasBeenInitialized()) scheduleUpdateGL();
  }

  void VolumeViewer::nearPlane(double n)
  {
    if(hasBeenInitialized()) makeCurrent();
    _cmRenderer.setNearPlane(n);
    _rgbaRenderer.setNearPlane(n);
    if(hasBeenInitialized()) scheduleUpdateGL();
  }

  void VolumeViewer::drawCornerAxis(bool draw)
  {
    _drawCornerAxis = draw;
    if(hasBeenInitialized()) scheduleUpdateGL();
  }

  void VolumeViewer::drawGeometry(bool draw)
  {
    _drawGeometry = draw;
    if(hasBeenInitialized()) scheduleUpdateGL();
  }

  void VolumeViewer::drawVolumes(bool draw)
  {
    _drawVolumes = draw;
    if(hasBeenInitialized()) scheduleUpdateGL();
  }

  void VolumeViewer::clipGeometryToBoundingBox(bool draw)
  {
    _clipGeometryToBoundingBox = draw;
    if(hasBeenInitialized()) scheduleUpdateGL();
  }

  void VolumeViewer::drawGeometryNormals(bool draw)
  {
    _drawGeometryNormals = draw;
    if(hasBeenInitialized()) scheduleUpdateGL();
  }

  // TODO: this is a hack from old code
  static inline void normalize(float vin[3])
  {
    float div = sqrt(vin[0]*vin[0]+vin[1]*vin[1]+vin[2]*vin[2]);
    
    if (div != 0.0)
      {
	vin[0] /= div;
	vin[1] /= div;
	vin[2] /= div;
      }
  }

  // arand, 4-21-2011, implemented
  void VolumeViewer::setShadedRenderingMode(bool val) {
    using namespace std;
    if (val != _shadedRenderingEnabled) { // don't churn if possible
      _shadedRenderingEnabled = val;
      if (val) {
	_cmRenderer.enableShadedRendering();
	_rgbaRenderer.enableShadedRendering();
      }
      else {
	_cmRenderer.disableShadedRendering();
	_rgbaRenderer.disableShadedRendering();
      }
      
      // reload the volume data
      _cmVolumeUploaded = false;
      _rgbaVolumeUploaded = false;
      
      if( interactiveSyncWithMultiTileServer() && syncShadedRenderWithMultiTileServer() )
#ifdef USE_XmlRpc
        needToUpdateMultiTileServer( -1, true );
#else
      syncShadedRenderFlag();
#endif
    }
  }

  VolumeViewer::scene_geometry_t::render_mode_t 
  VolumeViewer::scene_geometry_t::render_mode_enum(const std::string& str)
  {
    std::map<std::string,scene_geometry_t::render_mode_t> mapping;
    mapping["points"] =               POINTS;
    mapping["lines"] =                LINES;
    mapping["triangles"] =            TRIANGLES;
    mapping["triangle_wireframe"] =   TRIANGLE_WIREFRAME;
    mapping["triangle_filled_wire"] = TRIANGLE_FILLED_WIRE;
    mapping["triangles_flat"] =       TRIANGLES_FLAT;
    mapping["triangle_flat_filled_wire"] = TRIANGLE_FLAT_FILLED_WIRE;
    mapping["quads"] =                QUADS;
    mapping["quad_wireframe"] =       QUAD_WIREFRAME;
    mapping["quad_filled_wire"] =     QUAD_FILLED_WIRE;
    mapping["tetra"] =                TETRA;
    mapping["hexa"] =                 HEXA;
    if(mapping.find(str)==mapping.end())
      {
        cvcapp.log(1,boost::str(boost::format("%s: unknown render_mode %s, defaulting to 'points'\n")
                                % BOOST_CURRENT_FUNCTION
                                % str));
        return POINTS;
      }
    return mapping[str];
  }

  VolumeViewer::scene_geometry_collection VolumeViewer::geometries() const
  {
    using namespace std;
    scene_geometry_collection copy(_geometries);
    for(scene_geometry_collection::iterator i = copy.begin();
	i != copy.end();
	i++)
      i->second.reset(new scene_geometry_t(i->second->geometry,
					   i->second->name,
					   i->second->render_mode));
    return copy;
  }

  void VolumeViewer::geometries(const VolumeViewer::scene_geometry_collection& coll)
  {
    //Every geometry in the new collection that has an existing geometry
    //of the same name needs to be marked for VBO reinitialization.
    for (const auto& val : coll)
      if(_geometries.find(val.first)!=_geometries.end())
        _allocatedBuffersToReInit.push_back(val.first);
    _vboUpdated = false;

    _geometries = coll;
    for(scene_geometry_collection::iterator i = _geometries.begin();
	i != _geometries.end();
	i++)
      i->second.reset(new scene_geometry_t(i->second->geometry,
                                           i->second->name,
                                           i->second->render_mode));
  }

  void VolumeViewer::reset()
  {
    if(hasBeenInitialized()) makeCurrent();
    defaultConstructor();
    init();
    if(hasBeenInitialized()) scheduleUpdateGL();
  }

  void VolumeViewer::scheduleUpdateGL()
  {
    _glUpdated = false;
    QCoreApplication::postEvent(this,new CVCEvent("updateGL"));
  }

  // 09/18/2011 - transfix - adding display message functionality.  Doesn't work though :(
  //                         Something wrong with drawText() called by displayMessage()
  // 11/04/2011 - transfix - doing delayed scene normalization
  void VolumeViewer::timeout()
  {
    static log4cplus::Logger logger = log4cplus::Logger::getInstance("VolumeViewer.timeout");

    //If our name changed, re-copy the subvolume bounding box
    //if needed, referring to the new property as a result of the name
    //change.
    if(_oldObjectName != objectName())
      {
	copyBoundingBoxToDataMap();

	_oldObjectName = objectName();
      }

    //Display a message if one is available
    if(_doDisplayMessage)
      {
        LOG4CPLUS_TRACE(logger, BOOST_CURRENT_FUNCTION << " :: displaying message: " <<
                        _displayMessage.toStdString());
        // cvcapp.log(3,boost::str(boost::format("%s :: displaying message: %s\n")
        //                         % BOOST_CURRENT_FUNCTION
        //                         % _displayMessage.toStdString()));
        displayMessage(_displayMessage);
        _doDisplayMessage = false;
      }

    //normalize scene if needed
    if(_doNormalizeScene)
      {
        doNormalizeScene();
        _doNormalizeScene = false;
      }
  }

// void myReplace(std::string& str, const std::string& oldStr, const std::string& newStr)
// {
//   size_t pos = 0;
//   while((pos = str.find(oldStr, pos)) != std::string::npos)
//   {
//      str.replace(pos, oldStr.length(), newStr);
//      pos += newStr.length();
//   }
// }

  void VolumeViewer::init()
  {
    // static log4cplus::Logger logger = log4cplus::Logger::getInstance("VolumeViewer.init");
    // static log4cplus::Logger logger = getLogger(BOOST_CURRENT_FUNCTION);
    static log4cplus::Logger logger = FUNCTION_LOGGER;
    // static log4cplus::Logger logger = getLogger("VolumeViewer.init");

    setMouseTracking(true);

    // arand: size of 10 seems to work for me... 4-22-2011
    // transfix: moved this code here from CVC::Viewers 07/29/2011
    setSelectRegionWidth(10);
    setSelectRegionHeight(10);
    setTextIsEnabled(true);

    setDefaultColorTable();

    //initialize glew before initializing the renderers!!
    GLenum err = glewInit();
    if(GLEW_OK != err)
      LOG4CPLUS_ERROR(logger, glewGetErrorString(err));
      // cvcapp.log(1,
      //   boost::str(
      //     boost::format("%s: Error: %s\n")
      //     % BOOST_CURRENT_FUNCTION
      //     % glewGetErrorString(err)));
    LOG4CPLUS_DEBUG(logger, "Using GLEW " << glewGetString(GLEW_VERSION));
    // cvcapp.log(1,
    //   boost::str(
    //     boost::format("%s: Status: Using GLEW %s\n")
    //     % BOOST_CURRENT_FUNCTION
    //     % glewGetString(GLEW_VERSION)));
    
    if(!_cmRenderer.initRenderer() ||
       !_rgbaRenderer.initRenderer())
      _drawable = false;
    else
      {
	uploadColorTable();
	normalizeScene();
      }

    _hasBeenInitialized = true;
    emit postInit();
  }
  
  void VolumeViewer::draw()
  {
    if(interactiveSyncWithMultiTileServer() && syncCameraWithMultiTileServer())
      {
#ifdef USE_XmlRpc
        needToUpdateMultiTileServer( -1, true );
#else
        syncViewInformation();
#endif
      }

    copyCameraToPropertyMap();

    if(drawBoundingBox()) doDrawBoundingBox();
    if(drawSubVolumeSelector()) doDrawSubVolumeBoundingBox();
    if(drawGeometry()) doDrawGeometry();

    if(drawVolumes())
      {
        //Something to note about the volume renderer: it always draws volumes with the center of the volume
        //at the origin.  So to correctly position the volume in the scene, we must translate and scale it
        //so it fits its designated bounding box.
        glPushMatrix();
        glTranslatef(sceneCenter().x,sceneCenter().y,sceneCenter().z);
        switch(volumeRenderingType())
          {
          case ColorMapped:
            {
              if(!_cmVolumeUploaded)
                uploadColorMappedVolume();
	  
              double maxd = std::max(_cmVolume.XMax()-_cmVolume.XMin(),
                                     std::max(_cmVolume.YMax()-_cmVolume.YMin(),
                                              _cmVolume.ZMax()-_cmVolume.ZMin()));
              double aspectx = (_cmVolume.XMax()-_cmVolume.XMin())/maxd;
              double aspecty = (_cmVolume.YMax()-_cmVolume.YMin())/maxd;
              double aspectz = (_cmVolume.ZMax()-_cmVolume.ZMin())/maxd;

              glScalef((_cmVolume.XMax()-_cmVolume.XMin())/aspectx,
                       (_cmVolume.YMax()-_cmVolume.YMin())/aspecty,
                       (_cmVolume.ZMax()-_cmVolume.ZMin())/aspectz);

              _cmRenderer.renderVolume();
            }
            break;
          case RGBA:
            {
              if(!_rgbaVolumeUploaded)
                uploadRGBAVolume();
	  
              double maxd = std::max(_rgbaVolumes[0].XMax()-_rgbaVolumes[0].XMin(),
                                     std::max(_rgbaVolumes[0].YMax()-_rgbaVolumes[0].YMin(),
                                              _rgbaVolumes[0].ZMax()-_rgbaVolumes[0].ZMin()));
              double aspectx = (_rgbaVolumes[0].XMax()-_rgbaVolumes[0].XMin())/maxd;
              double aspecty = (_rgbaVolumes[0].YMax()-_rgbaVolumes[0].YMin())/maxd;
              double aspectz = (_rgbaVolumes[0].ZMax()-_rgbaVolumes[0].ZMin())/maxd;

              glScalef((_rgbaVolumes[0].XMax()-_rgbaVolumes[0].XMin())/aspectx,
                       (_rgbaVolumes[0].YMax()-_rgbaVolumes[0].YMin())/aspecty,
                       (_rgbaVolumes[0].ZMax()-_rgbaVolumes[0].ZMin())/aspectz);

              _rgbaRenderer.renderVolume();
            }
            break;
          default: break;
          }
        
        glPopMatrix();
      }

    if(drawSubVolumeSelector()) doDrawSubVolumeSelector();
  }

  void VolumeViewer::syncShadedRenderFlag()
  {
#ifndef USE_XmlRpc
     if( multiTileServerInitialized() ) {
        char flag = (char)_shadedRenderingEnabled;
        if( !m_socket->_send( SHDREN, &flag, sizeof(char) ) )
           fprintf( stderr, "syncShadedRenderFlag fail\n");
     }
#endif
  }

  void VolumeViewer::syncRenderModeFlag()
  {
#ifndef USE_XmlRpc
     if( multiTileServerInitialized() ) {
        char flag = (char)_volumeRenderingType;
        if( !m_socket->_send( RENMODE, &flag, sizeof(char) ) )
           fprintf( stderr, "syncRenderModeFlag fail\n");
     }
#endif
  }

  void VolumeViewer::syncColorTable()
  {
#ifndef USE_XmlRpc
     if( multiTileServerInitialized() ) {
#ifdef VOLUMEROVER2_FLOAT_COLORTABLE
        if( !m_socket->_send( COLTBL, _colorTable, sizeof(float)*1024 ) )
           fprintf( stderr, "syncColorTable fail\n");
#else
        if( !m_socket->_send( COLTBL, _colorTable, sizeof(unsigned char)*1024 ) )
           fprintf( stderr, "syncColorTable fail\n");
#endif
     }
#endif
  }

  void VolumeViewer::syncViewInformation()
  {
#ifndef USE_XmlRpc
     if( multiTileServerInitialized() ) {
       double mv[16], p[16];
       camera()->getModelViewMatrix( mv );
       if( !m_socket->_send( VIEW, mv, sizeof(double) * 16 ) )
         fprintf( stderr, "syncViewInformation fail\n");
#ifdef DEBUG_VOLUMEVIEWER
       else
         fprintf( stderr, "syncViewInformation\n");
#endif

       /*
         camera()->getProjectionMatrix( p );
         m_socket->send( PROJ, p, sizeof(double) * 16 );
       */
     }
#endif
  }

  void VolumeViewer::syncCurrentWithMultiTileServer()
  {
#ifdef USE_XmlRpc
    needToUpdateMultiTileServer( -1, true );
#else
    if( syncCameraWithMultiTileServer() )
      syncViewInformation();
    
    if( syncTransferFuncWithMultiTileServer() )
      syncColorTable();
    
    if( syncShadedRenderWithMultiTileServer() )
      syncShadedRenderFlag();

    if( syncRenderModeWithMultiTileServer() )
      syncRenderModeFlag();
#endif
  }

  void VolumeViewer::initializeMultiTileServer()
  {
     multiTileServerInitialized( true );

     _multiTileServer_nServer = cvcapp.properties<int>("MultiTileServer.nServer");
     _multiTileServer_hostList = new std::string[ _multiTileServer_nServer ];
     _multiTileServer_portList = new int[ _multiTileServer_nServer ];

     for( int i = 0; i < _multiTileServer_nServer; i ++ )
     {
        char key[2048];
        sprintf( key, "MultiTileServer.host_%d", i );
        _multiTileServer_hostList[i] = cvcapp.properties( key );
        sprintf( key, "MultiTileServer.port_%d", i );
        _multiTileServer_portList[i] = cvcapp.properties<int>( key );
     }
     int nTiles[2] = { cvcapp.properties<int>( "MultiTileServer.nTilesX" ), cvcapp.properties<int>( "MultiTileServer.nTilesY" ) };
     multiTileServerNtiles( nTiles[0], nTiles[1] );
     multiTileServerVolumeFile( cvcapp.properties( "MultiTileServer.volFile" ) );

#ifdef USE_XmlRpc
     if( _xmlrpcClients.size() < _multiTileServer_nServer ) {

        if( _xmlrpcClients.size() > 0 ) {
           // first terminate previous client threads
           terminateMultiTileClient( -1, true );
           bool allTerminated = false;
           while( !allTerminated ) {
              allTerminated = true;
              for( int i = 0; i < _xmlrpcClients.size(); i++ )
                 allTerminated = !terminateMultiTileClient( i );
           }
        }

        _xmlrpcClients.clear();
        for( int i = 0; i < _multiTileServer_nServer; i ++ ) {
           XmlRpcThread* xmlrpcT = new XmlRpcThread();
           _xmlrpcClients.push_back( xmlrpcT );
        }
        if( _updateMultiTileServer ) delete[] _updateMultiTileServer;
        _updateMultiTileServer = new bool[_multiTileServer_nServer];

        if( _terminateMultiTileClient ) delete[] _terminateMultiTileClient;
        _terminateMultiTileClient = new bool[_multiTileServer_nServer];

        needToUpdateMultiTileServer( -1, false );
        terminateMultiTileClient( -1, false );        
     }
     else {
        terminateMultiTileClient( -1, true );        
        bool allTerminated = false;
        while( !allTerminated ) {
           allTerminated = true;
           for( int i = 0; i < _xmlrpcClients.size(); i++ )
              allTerminated = !terminateMultiTileClient( i );
        }
     }
        
     for( int i = 0; i < _multiTileServer_nServer; i ++ ) {
        _xmlrpcClients[i]->start( this, _multiTileServer_hostList[i], _multiTileServer_portList[i], i );
     }
#else
     if( m_socket == NULL ) {
        m_socket = new CVCSocketTCP( CLIENT, _multiTileServer_nServer, 1024 );
        if( m_socket == NULL )
           multiTileServerInitialized( false );
     }

     m_socket->setTargets( _multiTileServer_nServer, _multiTileServer_hostList, _multiTileServer_portList);

     m_socket->_connect();

     if( !m_socket->_send( RESOL, nTiles, sizeof(int)*2 ) ) {
       fprintf( stderr, "syncTileResol fail\n");
       multiTileServerInitialized( false );
     }
#ifdef DEBUG_VOLUMEVIEWER
     else
       fprintf( stderr, "syncTileResol\n");
#endif
     
     const char *volFile = multiTileServerVolumeFile().data();
     if( !m_socket->_send( LOADVOL, volFile, strlen( volFile ) ) ) {
       fprintf( stderr, "syncVolumeFileName fail\n");
       multiTileServerInitialized( false );
     }
#ifdef DEBUG_VOLUMEVIEWER
     else
       fprintf( stderr, "syncVolumeFile\n");
#endif
     
#endif //USE_XmlRpc
  }

  void VolumeViewer::drawWithNames()
  {
    if(drawSubVolumeSelector()) doDrawSubVolumeSelector(true);
  }

  void VolumeViewer::postDraw()
  {
    QGLViewer::postDraw();
    if(drawCornerAxis()) doDrawCornerAxis();
  }

  void VolumeViewer::endSelection(const QPoint&)
  {
    glFlush();

    GLint nbHits = glRenderMode(GL_RENDER);
    
    if(nbHits <= 0)
      setSelectedName(-1);
    else
      {
	for(int i=0; i<nbHits; i++)
	  {
	    //prefer to select the heads over shafts
	    switch((selectBuffer())[i*4+3])
	      {
	      case MaxXHead:
	      case MinXHead:
	      case MaxYHead:
	      case MinYHead:
	      case MaxZHead:
	      case MinZHead:
		setSelectedName((selectBuffer())[i*4+3]);
		return;
	      }
	       
	    setSelectedName((selectBuffer())[i*4+3]);
	  }
      }
  }

  void VolumeViewer::postSelection(const QPoint& point)
  {
    //qDebug("Selected: %d", selectedName());
    _selectedPoint = point;
    _selectedObj = selectedName();
  }

  QString VolumeViewer::helpString() const
  {
    return QString("<h1>Nothing here right now...</h1>");
  }

  void VolumeViewer::doDrawBoundingBox()
  {
    double minx, miny, minz;
    double maxx, maxy, maxz;
    float bgcolor[4];

    switch(volumeRenderingType())
      {
      case RGBA:
	minx = _rgbaVolumes[0].XMin();
	miny = _rgbaVolumes[0].YMin();
	minz = _rgbaVolumes[0].ZMin();
	maxx = _rgbaVolumes[0].XMax();
	maxy = _rgbaVolumes[0].YMax();
	maxz = _rgbaVolumes[0].ZMax();
	break;
      default:;
      case ColorMapped:
	minx = _cmVolume.XMin();
	miny = _cmVolume.YMin();
	minz = _cmVolume.ZMin();
	maxx = _cmVolume.XMax();
	maxy = _cmVolume.YMax();
	maxz = _cmVolume.ZMax();
	break;
      }

    glPushAttrib(GL_CURRENT_BIT | GL_ENABLE_BIT);
    glLineWidth(1.0);
    glEnable(GL_LINE_SMOOTH);
    glDisable(GL_LIGHTING);
    glGetFloatv(GL_COLOR_CLEAR_VALUE, bgcolor);
    glColor4f(1.0-bgcolor[0],1.0-bgcolor[1],1.0-bgcolor[2],1.0-bgcolor[3]);
    /* front face */
    glBegin(GL_LINE_LOOP);
    glVertex3d(minx,miny,minz);
    glVertex3d(maxx,miny,minz);
    glVertex3d(maxx,maxy,minz);
    glVertex3d(minx,maxy,minz);
    glEnd();
    /* back face */
    glBegin(GL_LINE_LOOP);
    glVertex3d(minx,miny,maxz);
    glVertex3d(maxx,miny,maxz);
    glVertex3d(maxx,maxy,maxz);
    glVertex3d(minx,maxy,maxz);
    glEnd();
    /* connecting lines */
    glBegin(GL_LINES);
    glVertex3d(minx,maxy,minz);
    glVertex3d(minx,maxy,maxz);
    glVertex3d(minx,miny,minz);
    glVertex3d(minx,miny,maxz);
    glVertex3d(maxx,maxy,minz);
    glVertex3d(maxx,maxy,maxz);
    glVertex3d(maxx,miny,minz);
    glVertex3d(maxx,miny,maxz);
    glEnd();
    glDisable(GL_LINE_SMOOTH);
    glPopAttrib();
  }

  void VolumeViewer::doDrawSubVolumeBoundingBox()
  {
    double minx, miny, minz;
    double maxx, maxy, maxz;
    float bgcolor[4];
    
    switch(volumeRenderingType())
      {
      case RGBA:
	minx = _rgbaSubVolumeSelector.minx;
	miny = _rgbaSubVolumeSelector.miny;
	minz = _rgbaSubVolumeSelector.minz;
	maxx = _rgbaSubVolumeSelector.maxx;
	maxy = _rgbaSubVolumeSelector.maxy;
	maxz = _rgbaSubVolumeSelector.maxz;

	break;
      default:;
      case ColorMapped:
	minx = _cmSubVolumeSelector.minx;
	miny = _cmSubVolumeSelector.miny;
	minz = _cmSubVolumeSelector.minz;
	maxx = _cmSubVolumeSelector.maxx;
	maxy = _cmSubVolumeSelector.maxy;
	maxz = _cmSubVolumeSelector.maxz;

	break;
      }

    glPushAttrib(GL_CURRENT_BIT | GL_ENABLE_BIT);
    glLineWidth(1.0);
    glEnable(GL_LINE_SMOOTH);
    glDisable(GL_LIGHTING);
    glGetFloatv(GL_COLOR_CLEAR_VALUE, bgcolor);
    glColor4f(1.0-bgcolor[0],1.0-bgcolor[1],1.0-bgcolor[2],1.0-bgcolor[3]);
    /* front face */
    glBegin(GL_LINE_LOOP);
    glVertex3d(minx,miny,minz);
    glVertex3d(maxx,miny,minz);
    glVertex3d(maxx,maxy,minz);
    glVertex3d(minx,maxy,minz);
    glEnd();
    /* back face */
    glBegin(GL_LINE_LOOP);
    glVertex3d(minx,miny,maxz);
    glVertex3d(maxx,miny,maxz);
    glVertex3d(maxx,maxy,maxz);
    glVertex3d(minx,maxy,maxz);
    glEnd();
    /* connecting lines */
    glBegin(GL_LINES);
    glVertex3d(minx,maxy,minz);
    glVertex3d(minx,maxy,maxz);
    glVertex3d(minx,miny,minz);
    glVertex3d(minx,miny,maxz);
    glVertex3d(maxx,maxy,minz);
    glVertex3d(maxx,maxy,maxz);
    glVertex3d(maxx,miny,minz);
    glVertex3d(maxx,miny,maxz);
    glEnd();
    glDisable(GL_LINE_SMOOTH);
    glPopAttrib();
  }

  void VolumeViewer::doDrawSubVolumeSelector(bool withNames)
  {
    using namespace qglviewer;

    double minx, miny, minz;
    double maxx, maxy, maxz;
    
    switch(volumeRenderingType())
      {
      case RGBA:
	minx = _rgbaSubVolumeSelector.minx;
	miny = _rgbaSubVolumeSelector.miny;
	minz = _rgbaSubVolumeSelector.minz;
	maxx = _rgbaSubVolumeSelector.maxx;
	maxy = _rgbaSubVolumeSelector.maxy;
	maxz = _rgbaSubVolumeSelector.maxz;

	break;
      default:;
      case ColorMapped:
	minx = _cmSubVolumeSelector.minx;
	miny = _cmSubVolumeSelector.miny;
	minz = _cmSubVolumeSelector.minz;
	maxx = _cmSubVolumeSelector.maxx;
	maxy = _cmSubVolumeSelector.maxy;
	maxz = _cmSubVolumeSelector.maxz;

	break;
      }

    double centerx, centery, centerz;
    centerx = (maxx-minx)/2+minx;
    centery = (maxy-miny)/2+miny;
    centerz = (maxz-minz)/2+minz;

    /*** draw the selector handle arrows ***/
    float pointSize = withNames ? 6.0 : 6.0;
    float lineWidth = withNames ? 4.0 : 2.0;

    glClear(GL_DEPTH_BUFFER_BIT); //show arrows on top of everything

    drawHandle(Vec(centerx,centery,centerz),
    	       Vec(maxx,centery,centerz),
    	       1.0,0.0,0.0,
    	       MaxXHead,MaxXShaft,
   	       pointSize,lineWidth);
    drawHandle(Vec(centerx,centery,centerz),
	       Vec(minx,centery,centerz),
	       1.0,0.0,0.0,
	       MinXHead,MinXShaft,
	       pointSize,lineWidth);

    drawHandle(Vec(centerx,centery,centerz),
	       Vec(centerx,maxy,centerz),
	       0.0,1.0,0.0,
	       MaxYHead,MaxYShaft,
	       pointSize,lineWidth);
    drawHandle(Vec(centerx,centery,centerz),
	       Vec(centerx,miny,centerz),
	       0.0,1.0,0.0,
	       MinYHead,MinYShaft,
	       pointSize,lineWidth);
    
    drawHandle(Vec(centerx,centery,centerz),
	       Vec(centerx,centery,maxz),
	       0.0,0.0,1.0,
	       MaxZHead,MaxZShaft,
	       pointSize,lineWidth);
    drawHandle(Vec(centerx,centery,centerz),
	       Vec(centerx,centery,minz),
	       0.0,0.0,1.0,
	       MinZHead,MinZShaft,
	       pointSize,lineWidth);
  }

  void VolumeViewer::doDrawCornerAxis()
  {
    GLint viewport[4];
    GLint scissor[4];

    // The viewport and the scissor are changed to fit the lower left
    // corner. Original values are saved.
    glGetIntegerv(GL_VIEWPORT, viewport);
    glGetIntegerv(GL_SCISSOR_BOX, scissor);

    // Axis viewport size, in pixels
    const GLint size = 50;
    glViewport(0,0,size,size);
    glScissor(0,0,size,size);

    // The Z-buffer is cleared to make the axis appear over the
    // original image.
    glClear(GL_DEPTH_BUFFER_BIT);

    // Tune for best line rendering
    glDisable(GL_LIGHTING);
    glLineWidth(1.0);

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(-1, 1, -1, 1, -1, 1);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glMultMatrixd(camera()->orientation().inverse().matrix());

    glBegin(GL_LINES);
    glColor3f(1.0, 0.0, 0.0);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(1.0, 0.0, 0.0);

    glColor3f(0.0, 1.0, 0.0);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 1.0, 0.0);
	
    glColor3f(0.0, 0.0, 1.0);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 0.0, 1.0);
    glEnd();

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();

    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();

    glEnable(GL_LIGHTING);

    // The viewport and the scissor are restored.
    glScissor(scissor[0],scissor[1],scissor[2],scissor[3]);
    glViewport(viewport[0],viewport[1],viewport[2],viewport[3]);
/*
    // CHA: upon request of Joe, trying to use width and height of window for setting viewport.
    // This fixes the problem in OSX lion
    // This issue has been fixed by updating the version of glew
    float viewportsx = 0, viewportsy = 0, viewportw = width(), viewporth = height();
    glScissor (viewportsx, viewportsy, viewportw, viewporth );
    glViewport(viewportsx, viewportsy, viewportw, viewporth );
*/
  }

  void VolumeViewer::doDrawGeometry()
  {
    updateVBO();

    glPushAttrib(GL_LIGHTING_BIT | GL_ENABLE_BIT | GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
    
    using namespace qglviewer;
    const Vec camPos = camera()->position();
    //const GLfloat pos[4] = {camPos[0]+70.0f,camPos[1]+50.0f,camPos[2]+100.0f,1.0};
    const GLfloat pos[4] = {camPos[0],camPos[1],camPos[2],1.0};


    // arand: old lights
    //GLfloat diffuseColor[] = {1.0f, 1.0f, 1.0f, 1.0f};
    //GLfloat specularColor[] = {1.0f, 1.0f, 1.0f, 1.0f};
    //GLfloat ambientColor[] = {0.0f, 0.0f, 0.0f, 1.0f};

    GLfloat diffuseColor[] = {0.90f, 0.90f, 0.90f, 1.0f};
    GLfloat specularColor[] = {0.60f, 0.60f, 0.60f, 1.0f};
    GLfloat ambientColor[] = {0.0f, 0.0f, 0.0f, 1.0f};

    // arand, 7-19-2011
    // I am not sure if this code setting up the lighting options 
    // actually does anything... I just tried to copy TexMol
    glLightfv(GL_LIGHT0, GL_POSITION, pos);      
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseColor);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specularColor);
    glLightfv(GL_LIGHT0, GL_AMBIENT, ambientColor);


    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    
    glEnable(GL_LIGHT0);
    //glEnable(GL_LIGHT1);
    glEnable(GL_NORMALIZE);
    
    //// arand: added to render both sides of the surface... 4-12-2011
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

    if(clipGeometryToBoundingBox())
      setClipPlanes();

    for(scene_geometry_collection::iterator i = _geometries.begin();
	i != _geometries.end();
	i++)
      {
	if(!i->second) continue;

	scene_geometry_t &scene_geom = *(i->second.get());
        if(!scene_geom.visible) continue;

        //Choose render_geometry if available.  Used for TETRA and HEXA rendering.
        //-Joe, 4-19-2011
        cvcgeom_t &render_geometry = 
          !scene_geom.render_geometry.empty() ?
          scene_geom.render_geometry : scene_geom.geometry;
	cvcgeom_t &geom = render_geometry;

	GLint params[2];
	//back up current setting
	glGetIntegerv(GL_POLYGON_MODE,params);


	glEnable(GL_LIGHTING);

	//make sure we have normals!
	if(geom.const_points().size() != geom.const_normals().size()) {
	  geom.calculate_surf_normals();
	}

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);

	if (i->second->render_mode != scene_geometry_t::TRIANGLES_FLAT &&
	    i->second->render_mode != scene_geometry_t::TRIANGLE_FLAT_FILLED_WIRE) {

	  if(usingVBO())
	    {
	      glBindBufferARB(GL_ARRAY_BUFFER_ARB,scene_geom.vboArrayBufferID);
	      glVertexPointer(3, GL_DOUBLE, 0, (GLvoid*)scene_geom.vboArrayOffsets[0]);
	      glNormalPointer(GL_DOUBLE, 0, (GLvoid*)scene_geom.vboArrayOffsets[1]);
	    }
	  else
	    {
	      glVertexPointer(3, GL_DOUBLE, 0, &(geom.const_points()[0]));
	      glNormalPointer(GL_DOUBLE, 0, &(geom.const_normals()[0]));
	    }	  
	}
	
 
	if(geom.const_colors().size() == geom.const_points().size())
	  {
	    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	    glEnable(GL_COLOR_MATERIAL);
	    glEnableClientState(GL_COLOR_ARRAY);
	    if(usingVBO())
	      glColorPointer(3, GL_DOUBLE, 0, (GLvoid*)scene_geom.vboArrayOffsets[2]);
	    else
	      glColorPointer(3, GL_DOUBLE, 0, &(geom.const_colors()[0])); 
	  } 
	else {
	  glColor3f(1.0,1.0,1.0);
	}
	   
	switch(i->second->render_mode)
	  {
	  case scene_geometry_t::POINTS:
	    {
	      glPushAttrib(GL_LIGHTING_BIT | GL_POINT_BIT);
	      glDisable(GL_LIGHTING);
	      glPointSize(2.0);
	      glEnable(GL_POINT_SMOOTH);

	      //no need for normals when point rendering
	      glDisableClientState(GL_NORMAL_ARRAY);

	      if(usingVBO())
			glDrawArrays(GL_POINTS, 0, scene_geom.vboVertSize/3);
	      else
	    	glDrawArrays(GL_POINTS, 0, geom.const_points().size());

	      glPopAttrib();
	    }
	    break;
	  case scene_geometry_t::LINES:
	    glLineWidth(getWireframeLineWidth());
		glDisable(GL_LIGHTING);
		glEnable(GL_LINE_SMOOTH);
		//glHint(GL_LINE_SMOOTH,GL_NICEST); // arand: this doesn't seem to improve anything
		glDisableClientState(GL_NORMAL_ARRAY);

 		glPolygonMode(GL_FRONT, GL_LINE);

	    if(usingVBO())
	      {

		glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,
				scene_geom.vboLineElementArrayBufferID);
		glDrawElements(GL_LINES, 
			       geom.const_lines().size()*2,
			       GL_UNSIGNED_INT, 0);
	      }
	    else
	      glDrawElements(GL_LINES, 
			     geom.const_lines().size()*2,
			     GL_UNSIGNED_INT, &(geom.const_lines()[0]));
	    break;
	  case scene_geometry_t::TRIANGLES:

	    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	    if(usingVBO())
	      {
		glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,
				scene_geom.vboTriElementArrayBufferID);
		glDrawElements(GL_TRIANGLES, 
			       geom.const_triangles().size()*3,
			       GL_UNSIGNED_INT, 0);
	      }
	    else
	      glDrawElements(GL_TRIANGLES, 
			     geom.const_triangles().size()*3,
			     GL_UNSIGNED_INT, &(geom.const_triangles()[0]));


	    if (_drawGeometryNormals) {
	      // draw normals...
	      glLineWidth(2.0*getWireframeLineWidth());
	      glBegin(GL_LINES);
	      for (int i=0; i< geom.points().size(); i++) {
		double len = sqrt(geom.const_normals()[i][0]*geom.const_normals()[i][0]+geom.const_normals()[i][1]*geom.const_normals()[i][1]+geom.const_normals()[i][2]*geom.const_normals()[i][2]);
		len *= 2.0;
		glVertex3d(geom.const_points()[i][0],geom.const_points()[i][1],geom.const_points()[i][2]);
		glVertex3d(geom.const_points()[i][0]+geom.const_normals()[i][0]/len,
			   geom.const_points()[i][1]+geom.const_normals()[i][1]/len,
			   geom.const_points()[i][2]+geom.const_normals()[i][2]/len);	     	      
	      }
	      glEnd();

	    }

	    break;

	    case scene_geometry_t::TRIANGLES_FLAT:

	      // arand, 9-12-2011: added "flat" mode. This is probably much slower
	      //                   but displays un-oriented triangulations much nicer

	      glBegin(GL_TRIANGLES);
	      for (int i=0; i<geom.triangles().size(); i++) {
		int t1,t2,t0;
		t0 = geom.const_triangles()[i][0];
		t1 = geom.const_triangles()[i][1];
		t2 = geom.const_triangles()[i][2];

		double nx,ny,nz;		
		double v1x,v2x,v1y,v2y,v1z,v2z;
		
		v1x = geom.const_points()[t1][0] - geom.const_points()[t0][0];
		v1y = geom.const_points()[t1][1] - geom.const_points()[t0][1];
		v1z = geom.const_points()[t1][2] - geom.const_points()[t0][2];
		v2x = geom.const_points()[t2][0] - geom.const_points()[t0][0];
		v2y = geom.const_points()[t2][1] - geom.const_points()[t0][1];
		v2z = geom.const_points()[t2][2] - geom.const_points()[t0][2];

		double lv1 = sqrt(v1x*v1x+v1y*v1y+v1z*v1z);
		double lv2 = sqrt(v2x*v2x+v2y*v2y+v2z*v2z);

		nx = (v1y*v2z-v1z*v2y)/(lv1*lv2);
		ny = (v1z*v2x-v1x*v2z)/(lv1*lv2);
		nz = (v1x*v2y-v1y*v2x)/(lv1*lv2);

		glNormal3d(nx,ny,nz);
		if (!geom.const_colors().empty())
		  glColor3d(geom.const_colors()[t0][0],geom.const_colors()[t0][1],geom.const_colors()[t0][2]);
		glVertex3d(geom.const_points()[t0][0],geom.const_points()[t0][1],geom.const_points()[t0][2]);
		glVertex3d(geom.const_points()[t1][0],geom.const_points()[t1][1],geom.const_points()[t1][2]);
		glVertex3d(geom.const_points()[t2][0],geom.const_points()[t2][1],geom.const_points()[t2][2]);
	      }
	      glEnd();	      


	    if (_drawGeometryNormals) {
	      // draw normals...
	      glLineWidth(2.0*getWireframeLineWidth());
	      glBegin(GL_LINES);
	      for (int i=0; i< geom.points().size(); i++) {
		double len = sqrt(geom.const_normals()[i][0]*geom.const_normals()[i][0]+geom.const_normals()[i][1]*geom.const_normals()[i][1]+geom.const_normals()[i][2]*geom.const_normals()[i][2]);
		len *= 2.0;
		glVertex3d(geom.const_points()[i][0],geom.const_points()[i][1],geom.const_points()[i][2]);
		glVertex3d(geom.const_points()[i][0]+geom.const_normals()[i][0]/len,
			   geom.const_points()[i][1]+geom.const_normals()[i][1]/len,
			   geom.const_points()[i][2]+geom.const_normals()[i][2]/len);	     	      
	      }
	      glEnd();

	    }

	    break;
	  case scene_geometry_t::TRIANGLE_WIREFRAME:
	    glLineWidth(getWireframeLineWidth());

	    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	    if(usingVBO())
	      {
		glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,
				scene_geom.vboTriElementArrayBufferID);
		glDrawElements(GL_TRIANGLES, 
			       geom.const_triangles().size()*3,
			       GL_UNSIGNED_INT, 0);
	      }
	    else
	      glDrawElements(GL_TRIANGLES, 
			     geom.const_triangles().size()*3,
			     GL_UNSIGNED_INT, &(geom.const_triangles()[0]));
	    break;
	  case scene_geometry_t::TRIANGLE_FILLED_WIRE:
	    glLineWidth(getWireframeLineWidth());

	    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	    glEnable(GL_POLYGON_OFFSET_FILL);
	    glPolygonOffset(1.0,1.0);
	    if(usingVBO())
	      {
		glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,
				scene_geom.vboTriElementArrayBufferID);
		glDrawElements(GL_TRIANGLES, 
			       geom.const_triangles().size()*3,
			       GL_UNSIGNED_INT, 0);
	      }
	    else
	      glDrawElements(GL_TRIANGLES, 
			     geom.const_triangles().size()*3,
			     GL_UNSIGNED_INT, &(geom.const_triangles()[0]));	      


	    if (_drawGeometryNormals) {
	      // draw normals...
	      glLineWidth(2.0*getWireframeLineWidth());
	      glBegin(GL_LINES);
	      for (int i=0; i< geom.points().size(); i++) {
		double len = sqrt(geom.const_normals()[i][0]*geom.const_normals()[i][0]+geom.const_normals()[i][1]*geom.const_normals()[i][1]+geom.const_normals()[i][2]*geom.const_normals()[i][2]);
		len *= 2.0;
		glVertex3d(geom.const_points()[i][0],geom.const_points()[i][1],geom.const_points()[i][2]);
		glVertex3d(geom.const_points()[i][0]+geom.const_normals()[i][0]/len,
			   geom.const_points()[i][1]+geom.const_normals()[i][1]/len,
			   geom.const_points()[i][2]+geom.const_normals()[i][2]/len);	     	      
	      }
	      glEnd();

	    }

	    glPolygonOffset(0.0,0.0);
	    glDisable(GL_POLYGON_OFFSET_FILL);
	    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	    glDisable(GL_LIGHTING);
	    glDisableClientState(GL_COLOR_ARRAY);
	    glColor3f(0.0,0.0,0.0); //black wireframe
	    if(usingVBO())
	      {
		glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,
				scene_geom.vboTriElementArrayBufferID);
		glDrawElements(GL_TRIANGLES, 
			       geom.const_triangles().size()*3,
			       GL_UNSIGNED_INT, 0);
	      }
	    else
	      glDrawElements(GL_TRIANGLES, 
			     geom.const_triangles().size()*3,
			     GL_UNSIGNED_INT, &(geom.const_triangles()[0]));

	    
	    break;	  
	    
	    
	  case scene_geometry_t::TRIANGLE_FLAT_FILLED_WIRE:
	    glLineWidth(getWireframeLineWidth());

	    // new hack mode...
	    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	    glEnable(GL_POLYGON_OFFSET_FILL);
	    glPolygonOffset(1.0,1.0);
	    
	    glBegin(GL_TRIANGLES);
	    for (int i=0; i<geom.triangles().size(); i++) {
	      int t1,t2,t0;
	      t0 = geom.const_triangles()[i][0];
	      t1 = geom.const_triangles()[i][1];
	      t2 = geom.const_triangles()[i][2];
	      
	      double nx,ny,nz;
	      double v1x,v2x,v1y,v2y,v1z,v2z;
	      
	      v1x = geom.const_points()[t1][0] - geom.const_points()[t0][0];
	      v1y = geom.const_points()[t1][1] - geom.const_points()[t0][1];
	      v1z = geom.const_points()[t1][2] - geom.const_points()[t0][2];
	      v2x = geom.const_points()[t2][0] - geom.const_points()[t0][0];
	      v2y = geom.const_points()[t2][1] - geom.const_points()[t0][1];
	      v2z = geom.const_points()[t2][2] - geom.const_points()[t0][2];
		
	      double lv1 = sqrt(v1x*v1x+v1y*v1y+v1z*v1z);
	      double lv2 = sqrt(v2x*v2x+v2y*v2y+v2z*v2z);
	      
	      nx = (v1y*v2z-v1z*v2y)/(lv1*lv2);
	      ny = (v1z*v2x-v1x*v2z)/(lv1*lv2);
	      nz = (v1x*v2y-v1y*v2x)/(lv1*lv2);
	      
	      glNormal3d(nx,ny,nz);
	      if (!geom.const_colors().empty())
		glColor3d(geom.const_colors()[t0][0],geom.const_colors()[t0][1],geom.const_colors()[t0][2]);
	      glVertex3d(geom.const_points()[t0][0],geom.const_points()[t0][1],geom.const_points()[t0][2]);
	      glVertex3d(geom.const_points()[t1][0],geom.const_points()[t1][1],geom.const_points()[t1][2]);
	      glVertex3d(geom.const_points()[t2][0],geom.const_points()[t2][1],geom.const_points()[t2][2]);
	    }
	    glEnd();

	    if (_drawGeometryNormals) {
	      // draw normals...
	      glLineWidth(2.0*getWireframeLineWidth());
	      glBegin(GL_LINES);
	      for (int i=0; i< geom.points().size(); i++) {
		double len = sqrt(geom.const_normals()[i][0]*geom.const_normals()[i][0]+geom.const_normals()[i][1]*geom.const_normals()[i][1]+geom.const_normals()[i][2]*geom.const_normals()[i][2]);
		len *= 2.0;
		glVertex3d(geom.const_points()[i][0],geom.const_points()[i][1],geom.const_points()[i][2]);
		glVertex3d(geom.const_points()[i][0]+geom.const_normals()[i][0]/len,
			   geom.const_points()[i][1]+geom.const_normals()[i][1]/len,
			   geom.const_points()[i][2]+geom.const_normals()[i][2]/len);	     	      
	      }
	      glEnd();
	      
	    }
	      
	    glPolygonOffset(0.0,0.0);
	    glDisable(GL_POLYGON_OFFSET_FILL);
	    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	    glDisable(GL_LIGHTING);
	    glDisableClientState(GL_COLOR_ARRAY);
	    glColor3f(0.0,0.0,0.0); //black wireframe

	    // draw wireframe...
	    glBegin(GL_TRIANGLES);
	    for (int i=0; i<geom.triangles().size(); i++) {
	      int t1,t2,t0;
	      t0 = geom.const_triangles()[i][0];
	      t1 = geom.const_triangles()[i][1];
	      t2 = geom.const_triangles()[i][2];
	      
	      double nx,ny,nz;
	      double v1x,v2x,v1y,v2y,v1z,v2z;
	      
	      v1x = geom.const_points()[t1][0] - geom.const_points()[t0][0];
	      v1y = geom.const_points()[t1][1] - geom.const_points()[t0][1];
	      v1z = geom.const_points()[t1][2] - geom.const_points()[t0][2];
	      v2x = geom.const_points()[t2][0] - geom.const_points()[t0][0];
	      v2y = geom.const_points()[t2][1] - geom.const_points()[t0][1];
	      v2z = geom.const_points()[t2][2] - geom.const_points()[t0][2];
		
	      double lv1 = sqrt(v1x*v1x+v1y*v1y+v1z*v1z);
	      double lv2 = sqrt(v2x*v2x+v2y*v2y+v2z*v2z);
	      
	      nx = (v1y*v2z-v1z*v2y)/(lv1*lv2);
	      ny = (v1z*v2x-v1x*v2z)/(lv1*lv2);
	      nz = (v1x*v2y-v1y*v2x)/(lv1*lv2);
	      
	      glNormal3d(nx,ny,nz);
	      glVertex3d(geom.const_points()[t0][0],geom.const_points()[t0][1],geom.const_points()[t0][2]);
	      glVertex3d(geom.const_points()[t1][0],geom.const_points()[t1][1],geom.const_points()[t1][2]);
	      glVertex3d(geom.const_points()[t2][0],geom.const_points()[t2][1],geom.const_points()[t2][2]);
	    }
	    glEnd();

	    break;

	  case scene_geometry_t::QUADS:
	    if(usingVBO())
	      {
		glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,
				scene_geom.vboQuadElementArrayBufferID);
		glDrawElements(GL_QUADS, 
			       geom.const_quads().size()*4,
			       GL_UNSIGNED_INT, 0);
	      }
	    else
	      glDrawElements(GL_QUADS, 
			     geom.const_quads().size()*4,
			     GL_UNSIGNED_INT, &(geom.const_quads()[0]));
	    break;
	  case scene_geometry_t::QUAD_WIREFRAME:
	    glLineWidth(getWireframeLineWidth());

	    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	    if(usingVBO())
	      {
		glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,
				scene_geom.vboQuadElementArrayBufferID);
		glDrawElements(GL_QUADS, 
			       geom.const_quads().size()*4,
			       GL_UNSIGNED_INT, 0);
	      }
	    else
	      glDrawElements(GL_QUADS, 
			     geom.const_quads().size()*4,
			     GL_UNSIGNED_INT, &(geom.const_quads()[0]));
	    break;
	  case scene_geometry_t::QUAD_FILLED_WIRE:
	    glLineWidth(getWireframeLineWidth());

	    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	    glEnable(GL_POLYGON_OFFSET_FILL);
	    glPolygonOffset(1.0,1.0);
	    if(usingVBO())
	      {
		glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,
				scene_geom.vboQuadElementArrayBufferID);
		glDrawElements(GL_QUADS, 
			       geom.const_quads().size()*4,
			       GL_UNSIGNED_INT, 0);
	      }
	    else
	      glDrawElements(GL_QUADS, 
			     geom.const_quads().size()*4,
			     GL_UNSIGNED_INT, &(geom.const_quads()[0]));
	    glPolygonOffset(0.0,0.0);
	    glDisable(GL_POLYGON_OFFSET_FILL);
	    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	    glDisable(GL_LIGHTING);
	    glDisableClientState(GL_COLOR_ARRAY);
	    glColor3f(0.0,0.0,0.0); //black wireframe
	    if(usingVBO())
	      {
		glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,
				scene_geom.vboQuadElementArrayBufferID);
		glDrawElements(GL_QUADS, 
			       geom.const_quads().size()*4,
			       GL_UNSIGNED_INT, 0);
	      }
	    else
	      glDrawElements(GL_QUADS, 
			     geom.const_quads().size()*4,
			     GL_UNSIGNED_INT, &(geom.const_quads()[0]));
	    break;

	  case scene_geometry_t::TETRA: //TRIANGLES and LINES combined
	    glLineWidth(getWireframeLineWidth());

	    if(usingVBO())
	      {
		glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,
				scene_geom.vboTriElementArrayBufferID);
		glDrawElements(GL_TRIANGLES, 
			       geom.const_triangles().size()*3,
			       GL_UNSIGNED_INT, 0);
	      }
	    else
	      glDrawElements(GL_TRIANGLES,
			     geom.const_triangles().size()*3,
			     GL_UNSIGNED_INT, &(geom.const_triangles()[0]));

	    if(usingVBO())
	      {
		glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,
				scene_geom.vboLineElementArrayBufferID);
		glDrawElements(GL_LINES, 
			       geom.const_lines().size()*2,
			       GL_UNSIGNED_INT, 0);
	      }
	    else
	      glDrawElements(GL_LINES, 
			     geom.const_lines().size()*2,
			     GL_UNSIGNED_INT, &(geom.const_lines()[0]));

	    break;

	  case scene_geometry_t::HEXA: //QUADS and LINES combined
	    glLineWidth(getWireframeLineWidth());

	    if(usingVBO())
	      {
		glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,
				scene_geom.vboQuadElementArrayBufferID);
		glDrawElements(GL_QUADS, 
			       geom.const_quads().size()*4,
			       GL_UNSIGNED_INT, 0);
	      }
	    else
	      glDrawElements(GL_QUADS, 
			     geom.const_quads().size()*4,
			     GL_UNSIGNED_INT, &(geom.const_quads()[0]));

	    if(usingVBO())
	      {
		glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,
				scene_geom.vboLineElementArrayBufferID);
		glDrawElements(GL_LINES, 
			       geom.const_lines().size()*2,
			       GL_UNSIGNED_INT, 0);
	      }
	    else
	      glDrawElements(GL_LINES, 
			     geom.const_lines().size()*2,
			     GL_UNSIGNED_INT, &(geom.const_lines()[0]));
	  }

	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_TEXTURE_COORD_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);

	if(usingVBO())
	  {
	    glBindBufferARB(GL_ARRAY_BUFFER_ARB,0);
	    glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,0);
	  }

	//restore previous setting for polygon mode
	glPolygonMode(GL_FRONT,params[0]);
	glPolygonMode(GL_BACK,params[1]);
      }
      
    if(clipGeometryToBoundingBox())
      disableClipPlanes();

    glPopAttrib();
  }

  static inline unsigned int upToPowerOfTwo(unsigned int value)
  {
    unsigned int c = 0;
    unsigned int v = value;
    
    // round down to nearest power of two
    while (v>1) {
      v = v>>1;
      c++;
    }
    
    // if that isn't exactly the original value
    if ((v<<c)!=value) {
      // return the next power of two
      return (v<<(c+1));
    }
    else {
      // return this power of two
      return (v<<c);
    }
  }

  inline double texCoordOfSample(double sample, 
				 int bufferWidth, 
				 int canvasWidth,
				 double bufferMin, 
				 double bufferMax)
  {
    // get buffer min and max in the texture's space
    double texBufferMin = 0.5 / (double)canvasWidth;
    double texBufferMax = ((double)bufferWidth - 0.5) / (double)canvasWidth;
    
    return (sample-bufferMin)/(bufferMax-bufferMin) * (texBufferMax-texBufferMin) + texBufferMin;
  }

  void VolumeViewer::uploadColorMappedVolume()
  {
    if(hasBeenInitialized()) makeCurrent();

    if(_cmVolume.voxelType() != UChar)
      {
	_cmVolume.map(0.0,255.0);
	_cmVolume.voxelType(UChar);
      }
    
    VolMagick::uint64 upload_dimx = upToPowerOfTwo(_cmVolume.XDim());
    VolMagick::uint64 upload_dimy = upToPowerOfTwo(_cmVolume.YDim());
    VolMagick::uint64 upload_dimz = upToPowerOfTwo(_cmVolume.ZDim());
    boost::scoped_array<unsigned char> upload_buf(new unsigned char[upload_dimx*upload_dimy*upload_dimz]);
    for(VolMagick::uint64 k = 0; k < _cmVolume.ZDim(); k++)
      for(VolMagick::uint64 j = 0; j < _cmVolume.YDim(); j++)
	for(VolMagick::uint64 i = 0; i < _cmVolume.XDim(); i++)
	  upload_buf[k*upload_dimx*upload_dimy+j*upload_dimx+i] = 
	    static_cast<unsigned char>(_cmVolume(i,j,k));

    _cmRenderer.uploadColorMappedData(upload_buf.get(),upload_dimx,upload_dimy,upload_dimz);


    if (_shadedRenderingEnabled)
      {
	_cmRenderer.calculateGradientsFromDensities(upload_buf.get(),upload_dimx,upload_dimy,upload_dimz);
      }


    _cmRenderer.setAspectRatio(fabs(_cmVolume.XMax()-_cmVolume.XMin()),
			       fabs(_cmVolume.YMax()-_cmVolume.YMin()),
			       fabs(_cmVolume.ZMax()-_cmVolume.ZMin()));
    _cmRenderer.
      setTextureSubCube(texCoordOfSample(_cmVolume.XMin(),
					 _cmVolume.XDim(),
					 upload_dimx,
					 _cmVolume.XMin(),
					 _cmVolume.XMax()),
			texCoordOfSample(_cmVolume.YMin(),
					 _cmVolume.YDim(),
					 upload_dimy,
					 _cmVolume.YMin(),
					 _cmVolume.YMax()),
			texCoordOfSample(_cmVolume.ZMin(),
					 _cmVolume.ZDim(),
					 upload_dimz,
					 _cmVolume.ZMin(),
					 _cmVolume.ZMax()),
			texCoordOfSample(_cmVolume.XMax(),
					 _cmVolume.XDim(),
					 upload_dimx,
					 _cmVolume.XMin(),
					 _cmVolume.XMax()),
			texCoordOfSample(_cmVolume.YMax(),
					 _cmVolume.YDim(),
					 upload_dimy,
					 _cmVolume.YMin(),
					 _cmVolume.YMax()),
			texCoordOfSample(_cmVolume.ZMax(),
					 _cmVolume.ZDim(),
					 upload_dimz,
					 _cmVolume.ZMin(),
					 _cmVolume.ZMax()));

    _cmVolumeUploaded = true;
  }

  void VolumeViewer::uploadRGBAVolume()
  {
    if(hasBeenInitialized()) makeCurrent();

    for(unsigned int i = 0; i<4; i++)
      {
	if(_rgbaVolumes[i].voxelType() != UChar)
	  {
	    _rgbaVolumes[i].map(0.0,255.0);
	    _rgbaVolumes[i].voxelType(UChar);
	  }
      }

    VolMagick::uint64 upload_dimx = upToPowerOfTwo(_rgbaVolumes[0].XDim());
    VolMagick::uint64 upload_dimy = upToPowerOfTwo(_rgbaVolumes[0].YDim());
    VolMagick::uint64 upload_dimz = upToPowerOfTwo(_rgbaVolumes[0].ZDim());
    boost::scoped_array<unsigned char> upload_buf(new unsigned char[upload_dimx*upload_dimy*upload_dimz*4]);

    //QProgressDialog progress("Creating interleaved RGBA 3D texture...","Abort",0,_rgbaVolumes[0].ZDim());
    //progress.setWindowModality(Qt::WindowModal);
    
    for(VolMagick::uint64 k = 0; k < _rgbaVolumes[0].ZDim(); k++)
      {
	//progress.setValue(k);
	//qApp->processEvents();
	for(VolMagick::uint64 j = 0; j < _rgbaVolumes[0].YDim(); j++)
	  for(VolMagick::uint64 i = 0; i < _rgbaVolumes[0].XDim(); i++)
	    for(int cols = 0; cols < 4; cols++) {
	      upload_buf[k*upload_dimx*upload_dimy*4+
			 j*upload_dimx*4+
			 i*4+
			 cols] = 
		static_cast<unsigned char>(_rgbaVolumes[cols](i,j,k));
            }
      }
    //progress.setValue(_rgbaVolumes[0].ZDim());
    
    _rgbaRenderer.uploadRGBAData(upload_buf.get(),upload_dimx,upload_dimy,upload_dimz);

    if (_shadedRenderingEnabled)
      {
	_rgbaRenderer.calculateGradientsFromDensities(upload_buf.get(),upload_dimx,upload_dimy,upload_dimz);
      }

    _rgbaRenderer.setAspectRatio(fabs(_rgbaVolumes[0].XMax()-_rgbaVolumes[0].XMin()),
				 fabs(_rgbaVolumes[0].YMax()-_rgbaVolumes[0].YMin()),
				 fabs(_rgbaVolumes[0].ZMax()-_rgbaVolumes[0].ZMin()));
    _rgbaRenderer.
      setTextureSubCube(texCoordOfSample(_rgbaVolumes[0].XMin(),
					 _rgbaVolumes[0].XDim(),
					 upload_dimx,
					 _rgbaVolumes[0].XMin(),
					 _rgbaVolumes[0].XMax()),
			texCoordOfSample(_rgbaVolumes[0].YMin(),
					 _rgbaVolumes[0].YDim(),
					 upload_dimy,
					 _rgbaVolumes[0].YMin(),
					 _rgbaVolumes[0].YMax()),
			texCoordOfSample(_rgbaVolumes[0].ZMin(),
					 _rgbaVolumes[0].ZDim(),
					 upload_dimz,
					 _rgbaVolumes[0].ZMin(),
					 _rgbaVolumes[0].ZMax()),
			texCoordOfSample(_rgbaVolumes[0].XMax(),
					 _rgbaVolumes[0].XDim(),
					 upload_dimx,
					 _rgbaVolumes[0].XMin(),
					 _rgbaVolumes[0].XMax()),
			texCoordOfSample(_rgbaVolumes[0].YMax(),
					 _rgbaVolumes[0].YDim(),
					 upload_dimy,
					 _rgbaVolumes[0].YMin(),
					 _rgbaVolumes[0].YMax()),
			texCoordOfSample(_rgbaVolumes[0].ZMax(),
					 _rgbaVolumes[0].ZDim(),
					 upload_dimz,
					 _rgbaVolumes[0].ZMin(),
					 _rgbaVolumes[0].ZMax()));

    _rgbaVolumeUploaded = true;
  }

  void VolumeViewer::uploadColorTable()
  {
    if(hasBeenInitialized()) makeCurrent();
    _cmRenderer.uploadColorMap(_colorTable);
  }

  void VolumeViewer::updateVBO()
  {
    using namespace std;
    using namespace boost;

    if(_vboUpdated) return;

    //If this version of opengl doesn't support VBOs, don't do anything.
    if(!glewIsSupported("GL_ARB_vertex_buffer_object") 
#ifdef DISABLE_VBO_SUPPORT
       || true
#endif
       )
      {
	static bool already_printed = false;
	
	if(!already_printed)
	  {
#ifdef DISABLE_VBO_SUPPORT
	    cvcapp.log(2,str(boost::format("%s: VBO support disabled...\n") % BOOST_CURRENT_FUNCTION));
#else
	    cvcapp.log(2,str(boost::format("%s: no VBO support detected...\n") % BOOST_CURRENT_FUNCTION));
#endif
	    already_printed = true;
	  }

	_vboUpdated = true;
	_usingVBO = false;
	return;
      }

    //clean up all previously allocated buffers for geometry that has been removed
    vector<string> removeBuffers;
    for (auto val : _allocatedBuffersPerGeom) {
        if(_geometries.find(val.first)==_geometries.end())
          {
            removeBuffers.push_back(val.first);
            cvcapp.log(4,str(boost::format("%s: viewer %s removing buffers for geometry %s\n") 
                             % BOOST_CURRENT_FUNCTION
                             % getObjectName()
                             % val.first));
            vector<unsigned int> &allocatedBuffers = val.second;
            for(vector<unsigned int>::iterator i = allocatedBuffers.begin();
                i != allocatedBuffers.end();
                i++)
              if(glIsBufferARB(*i))
                glDeleteBuffersARB(1,reinterpret_cast<GLuint*>(&(*i)));
          }
      }

    //clean up VBO info for geometry that is going to be reinitialized
    for (const auto& name : _allocatedBuffersToReInit) {
        removeBuffers.push_back(name);
        cvcapp.log(4,str(boost::format("%s: viewer %s re-generating buffers for geometry %s\n") 
                         % BOOST_CURRENT_FUNCTION
                         % getObjectName()
                         % name));
        vector<unsigned int> &allocatedBuffers =
          _allocatedBuffersPerGeom[name];
        for(vector<unsigned int>::iterator i = allocatedBuffers.begin();
            i != allocatedBuffers.end();
            i++)
          if(glIsBufferARB(*i))
            glDeleteBuffersARB(1,reinterpret_cast<GLuint*>(&(*i)));
      }
    _allocatedBuffersToReInit.clear();

    //now remove the allocated buffers entry
    for (const auto& name : removeBuffers)
      _allocatedBuffersPerGeom.erase(name);

    for(scene_geometry_collection::iterator i = _geometries.begin();
	i != _geometries.end();
	i++)
      {
	if(!i->second) continue;

        //if we have already updated the vbo for this geometry previously, skip it
        if(_allocatedBuffersPerGeom.find(i->first)!=_allocatedBuffersPerGeom.end())
          continue;

        cvcapp.log(4,str(boost::format("%s: viewer %s generating buffers for geometry %s\n") 
                         % BOOST_CURRENT_FUNCTION
                         % getObjectName()
                         % i->first));

	scene_geometry_t &scene_geom = *(i->second.get());
        cvcgeom_t &render_geometry = 
          !scene_geom.render_geometry.empty() ?
          scene_geom.render_geometry : scene_geom.geometry;
	cvcgeom_t &geom = render_geometry;
	
	//make sure we have normals!
	if(geom.const_points().size() != geom.const_normals().size())
	  geom.calculate_surf_normals();
	
	//get new buffer ids
	glGenBuffersARB(1,
          reinterpret_cast<GLuint*>(&scene_geom.vboArrayBufferID));
	glGenBuffersARB(1,
          reinterpret_cast<GLuint*>(&scene_geom.vboLineElementArrayBufferID));
	glGenBuffersARB(1,
          reinterpret_cast<GLuint*>(&scene_geom.vboTriElementArrayBufferID));
	glGenBuffersARB(1,
          reinterpret_cast<GLuint*>(&scene_geom.vboQuadElementArrayBufferID));
	_allocatedBuffersPerGeom[i->first]
          .push_back(scene_geom.vboArrayBufferID);
	_allocatedBuffersPerGeom[i->first]
          .push_back(scene_geom.vboLineElementArrayBufferID);
	_allocatedBuffersPerGeom[i->first]
          .push_back(scene_geom.vboTriElementArrayBufferID);
	_allocatedBuffersPerGeom[i->first]
          .push_back(scene_geom.vboQuadElementArrayBufferID);
	
	unsigned int bufsize =
	  geom.const_points().size()*3 +
	  geom.const_normals().size()*3 +
	  geom.const_colors().size()*3;
	scoped_array<double>
	  buf(new double[bufsize]);
	  
	for(unsigned int i = 0; i < geom.const_points().size(); i++)
	  for(unsigned int j = 0; j < 3; j++)
	    buf[i*3+j] = geom.const_points()[i][j];
	  
	for(unsigned int i = 0, buf_i = geom.const_points().size();
	    i < geom.const_normals().size(); 
	    i++, buf_i++)
	  for(unsigned int j = 0; j < 3; j++)
	    buf[buf_i*3+j] = geom.const_normals()[i][j];
	  
	for(unsigned int i = 0, buf_i = geom.const_points().size()+geom.const_normals().size();
	    i < geom.const_colors().size(); 
	    i++, buf_i++)
	   for(unsigned int j = 0; j < 3; j++)
	   	buf[buf_i*3+j] = (double)geom.const_colors()[i][j];

	scene_geom.vboVertSize = geom.const_points().size()*3;
	scene_geom.vboArrayOffsets[0] = 0;
	scene_geom.vboArrayOffsets[1] = geom.const_points().size()*3*sizeof(double);
	scene_geom.vboArrayOffsets[2] = (geom.const_points().size()+geom.const_normals().size())*3*sizeof(double);
	  
	//upload vertex info
	glBindBufferARB(GL_ARRAY_BUFFER_ARB,scene_geom.vboArrayBufferID);
	glBufferDataARB(GL_ARRAY_BUFFER_ARB,
			bufsize*sizeof(double),
			buf.get(),
			GL_STATIC_DRAW_ARB);
	  
	//upload element info
	bufsize = geom.const_lines().size()*2;
	scoped_array<int> element_buf(new int[bufsize]);
	for(unsigned int i = 0; i < geom.const_lines().size(); i++)
	  for(unsigned int j = 0; j < 2; j++)
	    element_buf[i*2+j] = geom.const_lines()[i][j];
	  
	glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,scene_geom.vboLineElementArrayBufferID);
	glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER_ARB,
			bufsize*sizeof(int),
			element_buf.get(),
			GL_STATIC_DRAW_ARB);
	scene_geom.vboLineSize = bufsize;
	  
	bufsize = geom.const_triangles().size()*3;
	element_buf.reset(new int[bufsize]);
	for(unsigned int i = 0; i < geom.const_triangles().size(); i++)
	  for(unsigned int j = 0; j < 3; j++)
	    element_buf[i*3+j] = geom.const_triangles()[i][j];
	  
	glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,scene_geom.vboTriElementArrayBufferID);
	glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER_ARB,
			bufsize*sizeof(int),
			element_buf.get(),
			GL_STATIC_DRAW_ARB);
	scene_geom.vboTriSize = bufsize;
	  
	bufsize = geom.const_quads().size()*4;
	element_buf.reset(new int[bufsize]);
	for(unsigned int i = 0; i < geom.const_quads().size(); i++)
	  for(unsigned int j = 0; j < 4; j++)
	    element_buf[i*4+j] = geom.const_quads()[i][j];
	  
	glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,scene_geom.vboQuadElementArrayBufferID);
	glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER_ARB,
			bufsize*sizeof(int),
			element_buf.get(),
			GL_STATIC_DRAW_ARB);
	scene_geom.vboQuadSize = bufsize;
      }
    
    _usingVBO = true;
    _vboUpdated = true;
  }

  void VolumeViewer::setClipPlanes()
  {
    VolMagick::BoundingBox bbox = boundingBox();
    
    double plane0[] = { 0.0, 0.0, -1.0, bbox.maxz };
    glClipPlane(GL_CLIP_PLANE0, plane0);
    glEnable(GL_CLIP_PLANE0);

    double plane1[] = { 0.0, 0.0, 1.0, -bbox.minz };
    glClipPlane(GL_CLIP_PLANE1, plane1);
    glEnable(GL_CLIP_PLANE1);

    double plane2[] = { 0.0, -1.0, 0.0, bbox.maxy };
    glClipPlane(GL_CLIP_PLANE2, plane2);
    glEnable(GL_CLIP_PLANE2);

    double plane3[] = { 0.0, 1.0, 0.0, -bbox.miny };
    glClipPlane(GL_CLIP_PLANE3, plane3);
    glEnable(GL_CLIP_PLANE3);

    double plane4[] = { -1.0, 0.0, 0.0, bbox.maxx };
    glClipPlane(GL_CLIP_PLANE4, plane4);
    glEnable(GL_CLIP_PLANE4);

    double plane5[] = { 1.0, 0.0, 0.0, -bbox.minx };
    glClipPlane(GL_CLIP_PLANE5, plane5);
    glEnable(GL_CLIP_PLANE5);
  }

  void VolumeViewer::disableClipPlanes()
  {
    glDisable(GL_CLIP_PLANE0);
    glDisable(GL_CLIP_PLANE1);
    glDisable(GL_CLIP_PLANE2);
    glDisable(GL_CLIP_PLANE3);
    glDisable(GL_CLIP_PLANE4);
    glDisable(GL_CLIP_PLANE5);
  }

  void VolumeViewer::normalizeScene()
  {
    _doNormalizeScene = true;
  }

  // 11/04/2011 - transfix - renaming to doNormalizeScene for delayed scene normalization
  void VolumeViewer::doNormalizeScene()
  {
    static log4cplus::Logger logger = FUNCTION_LOGGER;

    // 2/22/2012 - cha - adding feature to disable scene normalization
    if (cvcapp.hasProperty(getObjectName("normalize_scene"))) {
      bool flag = cvcapp.properties<bool>(getObjectName("normalize_scene"));
      if( !flag ) { fprintf( stderr, "normalize scene disabled\n"); return; }
    }

    switch(volumeRenderingType())
      {
      case RGBA:
	setSceneCenter(qglviewer::Vec((_rgbaVolumes[0].XMax()-_rgbaVolumes[0].XMin())/2.0+_rgbaVolumes[0].XMin(),
				      (_rgbaVolumes[0].YMax()-_rgbaVolumes[0].YMin())/2.0+_rgbaVolumes[0].YMin(),
				      (_rgbaVolumes[0].ZMax()-_rgbaVolumes[0].ZMin())/2.0+_rgbaVolumes[0].ZMin()));
	setSceneRadius(std::max(_rgbaVolumes[0].XMax()-_rgbaVolumes[0].XMin(),
			   std::max(_rgbaVolumes[0].YMax()-_rgbaVolumes[0].YMin(),
			       _rgbaVolumes[0].ZMax()-_rgbaVolumes[0].ZMin())));
	break;
      case ColorMapped:
	setSceneCenter(qglviewer::Vec((_cmVolume.XMax()-_cmVolume.XMin())/2.0+_cmVolume.XMin(),
				      (_cmVolume.YMax()-_cmVolume.YMin())/2.0+_cmVolume.YMin(),
				      (_cmVolume.ZMax()-_cmVolume.ZMin())/2.0+_cmVolume.ZMin()));
	setSceneRadius(std::max(_cmVolume.XMax()-_cmVolume.XMin(),
			   std::max(_cmVolume.YMax()-_cmVolume.YMin(),
			       _cmVolume.ZMax()-_cmVolume.ZMin()))/2.0  );
	break;
      }

    if(copyCameraOnNormalize())
      copyCameraToPropertyMap();
    

    if(hasBeenInitialized() &&
       showEntireSceneOnNormalize()) 
      {
        LOG4CPLUS_TRACE(logger, getObjectName() << ": showing entire scene");
        // cvcapp.log(
        //   2,
        //   boost::str(
        //     boost::format("%s: %s: showing entire scene\n")
        //     % BOOST_CURRENT_FUNCTION
        //     % getObjectName()
        //   )
        // );
        showEntireScene();
      }
  }

  void VolumeViewer::drawHandle(const qglviewer::Vec& from, const qglviewer::Vec& to, 
			  float radius, int nbSubdivisions, 
			  int headName, int shaftName)
  {
    GLUquadric* quadric = gluNewQuadric();

    glPushMatrix();
    glTranslatef(from[0],from[1],from[2]);
    glMultMatrixd(qglviewer::Quaternion(qglviewer::Vec(0,0,1), to-from).matrix());

    double length = (to-from).norm();
    if (radius < 0.0) radius = 0.01 * length;
    const float head = 2.5*(radius / length) + 0.1;
    const float coneRadiusCoef = 4.0 - 5.0 * head;
    
    glPushName(shaftName);
    gluCylinder(quadric, radius, radius, length * (1.0 - head/coneRadiusCoef), nbSubdivisions, 1);
    glPopName();
    glTranslatef(0.0, 0.0, length * (1.0 - head));
    glPushName(headName);
    gluCylinder(quadric, coneRadiusCoef * radius, 0.0, head * length, nbSubdivisions, 1);
    glPopName();
    glTranslatef(0.0, 0.0, -length * (1.0 - head));
    
    glPopMatrix();

    gluDeleteQuadric(quadric);
  }

  void VolumeViewer::drawHandle(const qglviewer::Vec& from, const qglviewer::Vec& to, 
			  float r, float g, float b,
			  int headName, int shaftName,
			  float pointSize, float lineWidth)
  {
    float bgcolor[4];

    glPushAttrib(GL_ENABLE_BIT);
    glDisable(GL_LIGHTING);

    //draw shaft
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(lineWidth);
    glColor3f(r,g,b);
    glPushName(shaftName);
    glBegin(GL_LINES);
    glVertex3f(from[0],from[1],from[2]);
    glVertex3f(to[0],to[1],to[2]);
    glEnd();
    glPopName();
    glDisable(GL_LINE_SMOOTH);

    //draw head
    glPointSize(pointSize);
    glGetFloatv(GL_COLOR_CLEAR_VALUE, bgcolor);
    glColor3f(1.0-bgcolor[0],1.0-bgcolor[1],1.0-bgcolor[2]);
    glPushName(headName);
    glBegin(GL_POINTS);
    glVertex3f(to[0],to[1],to[2]);
    glEnd();
    glPopName();

    glPopAttrib();
  }

  void VolumeViewer::mousePressEvent(QMouseEvent *e)
  {
    select(e->pos());
    if(_selectedObj != -1) {
      _selectedPoint = e->pos();
    }
    else {
      QGLViewer::mousePressEvent(e);
    }
    // to test handling XmlRpc delay
    _mousePressed = true;
  }

  void VolumeViewer::mouseMoveEvent(QMouseEvent *e)
  {
    using namespace qglviewer;

    if(_selectedObj != -1)
      {
	double minx, miny, minz;
	double maxx, maxy, maxz;
	
	switch(volumeRenderingType())
	  {
	  case RGBA:
	    minx = _rgbaSubVolumeSelector.minx;
	    miny = _rgbaSubVolumeSelector.miny;
	    minz = _rgbaSubVolumeSelector.minz;
	    maxx = _rgbaSubVolumeSelector.maxx;
	    maxy = _rgbaSubVolumeSelector.maxy;
	    maxz = _rgbaSubVolumeSelector.maxz;
	    break;
	  default:;
	  case ColorMapped:
	    minx = _cmSubVolumeSelector.minx;
	    miny = _cmSubVolumeSelector.miny;
	    minz = _cmSubVolumeSelector.minz;
	    maxx = _cmSubVolumeSelector.maxx;
	    maxy = _cmSubVolumeSelector.maxy;
	    maxz = _cmSubVolumeSelector.maxz;
	    break;
	  }

	Vec center((maxx-minx)/2+minx,
		   (maxy-miny)/2+miny,
		   (maxz-minz)/2+minz);
	
	Vec screenCenter = camera()->projectedCoordinatesOf(center);
	Vec transPoint(e->pos().x(),e->pos().y(),screenCenter.z);
	Vec selectPoint(_selectedPoint.x(),_selectedPoint.y(),screenCenter.z);

	switch(_selectedObj)
	  {
	  case MaxXHead:
	  case MaxXShaft:
	    {
	      Vec maxXrayDir = Vec(maxx,center.y,center.z) - center;

	      Vec transPointWorld = camera()->unprojectedCoordinatesOf(transPoint);
	      Vec selectPointWorld = camera()->unprojectedCoordinatesOf(selectPoint);

	      Vec transPointWorldDir = transPointWorld - center;
	      Vec selectPointWorldDir = selectPointWorld - center;

	      transPointWorldDir.projectOnAxis(maxXrayDir);
	      selectPointWorldDir.projectOnAxis(maxXrayDir);
	      
	      transPointWorldDir -= selectPointWorldDir;

	      maxx += transPointWorldDir.x;

	      if(_selectedObj == MaxXShaft)
		minx += transPointWorldDir.x;
	      else
		minx -= transPointWorldDir.x;
	    }
	    break;
	  case MinXHead:
	  case MinXShaft:
	    {
	      Vec minXrayDir = Vec(minx,center.y,center.z) - center;

	      Vec transPointWorld = camera()->unprojectedCoordinatesOf(transPoint);
	      Vec selectPointWorld = camera()->unprojectedCoordinatesOf(selectPoint);

	      Vec transPointWorldDir = transPointWorld - center;
	      Vec selectPointWorldDir = selectPointWorld - center;

	      transPointWorldDir.projectOnAxis(minXrayDir);
	      selectPointWorldDir.projectOnAxis(minXrayDir);
	      
	      transPointWorldDir -= selectPointWorldDir;

	      minx += transPointWorldDir.x;

	      if(_selectedObj == MinXShaft)
		maxx += transPointWorldDir.x;
	      else
		maxx -= transPointWorldDir.x;
	    }
	    break;
	  case MaxYHead:
	  case MaxYShaft:
	    {
	      Vec maxYrayDir = Vec(center.x,maxy,center.z) - center;

	      Vec transPointWorld = camera()->unprojectedCoordinatesOf(transPoint);
	      Vec selectPointWorld = camera()->unprojectedCoordinatesOf(selectPoint);

	      Vec transPointWorldDir = transPointWorld - center;
	      Vec selectPointWorldDir = selectPointWorld - center;

	      transPointWorldDir.projectOnAxis(maxYrayDir);
	      selectPointWorldDir.projectOnAxis(maxYrayDir);
	      
	      transPointWorldDir -= selectPointWorldDir;

	      maxy += transPointWorldDir.y;

	      if(_selectedObj == MaxYShaft)
		miny += transPointWorldDir.y;
	      else
		miny -= transPointWorldDir.y;
	    }
	    break;
	  case MinYHead:
	  case MinYShaft:
	    {
	      Vec minYrayDir = Vec(center.x,miny,center.z) - center;

	      Vec transPointWorld = camera()->unprojectedCoordinatesOf(transPoint);
	      Vec selectPointWorld = camera()->unprojectedCoordinatesOf(selectPoint);

	      Vec transPointWorldDir = transPointWorld - center;
	      Vec selectPointWorldDir = selectPointWorld - center;

	      transPointWorldDir.projectOnAxis(minYrayDir);
	      selectPointWorldDir.projectOnAxis(minYrayDir);
	      
	      transPointWorldDir -= selectPointWorldDir;

	      miny += transPointWorldDir.y;

	      if(_selectedObj == MinYShaft)
		maxy += transPointWorldDir.y;
	      else
		maxy -= transPointWorldDir.y;
	    }
	    break;
	  case MaxZHead:
	  case MaxZShaft:
	    {
	      Vec maxZrayDir = Vec(center.x,center.y,maxz) - center;

	      Vec transPointWorld = camera()->unprojectedCoordinatesOf(transPoint);
	      Vec selectPointWorld = camera()->unprojectedCoordinatesOf(selectPoint);

	      Vec transPointWorldDir = transPointWorld - center;
	      Vec selectPointWorldDir = selectPointWorld - center;

	      transPointWorldDir.projectOnAxis(maxZrayDir);
	      selectPointWorldDir.projectOnAxis(maxZrayDir);
	      
	      transPointWorldDir -= selectPointWorldDir;

	      maxz += transPointWorldDir.z;

	      if(_selectedObj == MaxZShaft)
		minz += transPointWorldDir.z;
	      else
		minz -= transPointWorldDir.z;
	    }
	    break;
	  case MinZHead:
	  case MinZShaft:
	    {
	      Vec minZrayDir = Vec(center.x,center.y,minz) - center;

	      Vec transPointWorld = camera()->unprojectedCoordinatesOf(transPoint);
	      Vec selectPointWorld = camera()->unprojectedCoordinatesOf(selectPoint);

	      Vec transPointWorldDir = transPointWorld - center;
	      Vec selectPointWorldDir = selectPointWorld - center;

	      transPointWorldDir.projectOnAxis(minZrayDir);
	      selectPointWorldDir.projectOnAxis(minZrayDir);
	      
	      transPointWorldDir -= selectPointWorldDir;

	      minz += transPointWorldDir.z;

	      if(_selectedObj == MinZShaft)
		maxz += transPointWorldDir.z;
	      else
		maxz -= transPointWorldDir.z;
	    }
	    break;
	  }

	//if we ran into trouble with our above math, abort this bounding box change
#ifdef _WIN32
	if(_isnan(minx) || _isnan(miny) || _isnan(minz) ||
	   _isnan(maxx) || _isnan(maxy) || _isnan(maxz))
	  return;
#else
	if(std::isnan(minx) || std::isnan(miny) || std::isnan(minz) ||
	   std::isnan(maxx) || std::isnan(maxy) || std::isnan(maxz))
	  return;
#endif

	switch(volumeRenderingType())
	  {
	  case RGBA:
	    _rgbaSubVolumeSelector.minx = std::min(std::max(minx,_rgbaVolumes[0].XMin()),_rgbaVolumes[0].XMax());
	    _rgbaSubVolumeSelector.miny = std::min(std::max(miny,_rgbaVolumes[0].YMin()),_rgbaVolumes[0].YMax());
	    _rgbaSubVolumeSelector.minz = std::min(std::max(minz,_rgbaVolumes[0].ZMin()),_rgbaVolumes[0].ZMax());
	    _rgbaSubVolumeSelector.maxx = std::max(std::min(maxx,_rgbaVolumes[0].XMax()),_rgbaVolumes[0].XMin());
	    _rgbaSubVolumeSelector.maxy = std::max(std::min(maxy,_rgbaVolumes[0].YMax()),_rgbaVolumes[0].YMin());
	    _rgbaSubVolumeSelector.maxz = std::max(std::min(maxz,_rgbaVolumes[0].ZMax()),_rgbaVolumes[0].ZMin());
	    break;
	  default:;
	  case ColorMapped:
	    _cmSubVolumeSelector.minx = std::min(std::max(minx,_cmVolume.XMin()),_cmVolume.XMax());
	    _cmSubVolumeSelector.miny = std::min(std::max(miny,_cmVolume.YMin()),_cmVolume.YMax());
	    _cmSubVolumeSelector.minz = std::min(std::max(minz,_cmVolume.ZMin()),_cmVolume.ZMax());
	    _cmSubVolumeSelector.maxx = std::max(std::min(maxx,_cmVolume.XMax()),_cmVolume.XMin());
	    _cmSubVolumeSelector.maxy = std::max(std::min(maxy,_cmVolume.YMax()),_cmVolume.YMin());
	    _cmSubVolumeSelector.maxz = std::max(std::min(maxz,_cmVolume.ZMax()),_cmVolume.ZMin());

// 	    qDebug("(%f,%f,%f) (%f,%f,%f)",
// 	       _cmSubVolumeSelector.minx,_cmSubVolumeSelector.miny,_cmSubVolumeSelector.minz,
// 	       _cmSubVolumeSelector.maxx,_cmSubVolumeSelector.maxy,_cmSubVolumeSelector.maxz);
	    break;
	  }

	_selectedPoint = e->pos();

	if(hasBeenInitialized()) scheduleUpdateGL();
      }
    else
      {
        QGLViewer::mouseMoveEvent(e);

        //copyCameraToPropertyMap();

        if(_mousePressed)
          _mouseMoved = true;
      }
  }

  void VolumeViewer::wheelEvent(QWheelEvent* e)
  {
    if(_mousePressed)
      _mouseMoved = true;
    QGLViewer::wheelEvent(e);
  }

  std::vector<cvcraw_geometry::cvcgeom_t> VolumeViewer::getGeometriesFromDatamap() const
  {
    using namespace std;
    using namespace boost;
    using namespace boost::algorithm;
    
    vector<cvcraw_geometry::cvcgeom_t> geoms =
      cvcapp.propertyData<cvcraw_geometry::cvcgeom_t>(getObjectName("geometries"),true);

#ifdef USING_TILING
    vector<cvcraw_geometry::contours_t> contours_sets =
      cvcapp.propertyData<cvcraw_geometry::contours_t>(getObjectName("geometries"),true);
    for (const auto& contours : contours_sets) {
      geoms.push_back(contours.geom());
    }
#endif

    return geoms;
  }

  void VolumeViewer::copyBoxToDatamapAndEmitChange()
  {
    VolMagick::BoundingBox emitbox;
    switch(volumeRenderingType())
      {
      case RGBA:
        emitbox = _rgbaSubVolumeSelector;
        break;
      case ColorMapped:
        emitbox = _cmSubVolumeSelector;
        break;
      }
    
    copyBoundingBoxToDataMap();
    emit subVolumeSelectorChanged(emitbox);
  }

  void VolumeViewer::mouseReleaseEvent(QMouseEvent *e)
  {
    if(_selectedObj != -1)
      {
	_selectedObj = -1;
	
	//lets normalize the bounding boxes in case they're inverted (i.e. min > max)
	_rgbaSubVolumeSelector.normalize();
	_cmSubVolumeSelector.normalize();
	
        cvcapp.log(
          6,
          boost::str(
            boost::format("%s: (%f,%f,%f) (%f,%f,%f)\n")
            % BOOST_CURRENT_FUNCTION
            % _cmSubVolumeSelector.minx 
            % _cmSubVolumeSelector.miny 
            % _cmSubVolumeSelector.minz
            % _cmSubVolumeSelector.maxx 
            % _cmSubVolumeSelector.maxy 
            % _cmSubVolumeSelector.maxz
          )
        );

        copyBoxToDatamapAndEmitChange();
      }

    QGLViewer::mouseReleaseEvent(e);
    // to test handling XmlRpc delay
    _mousePressed = false;

    // for quaternion updates
    _mouseMoved = false;
  }

  void VolumeViewer::customEvent(QEvent *event)
  {
    CVCEvent *mwe = dynamic_cast<CVCEvent*>(event);
    if(!mwe) return;

    if(mwe->name == "handlePropertiesChanged")
      handlePropertiesChanged(boost::any_cast<std::string>(mwe->data));
    if(mwe->name == "handleDataChanged")
      handleDataChanged(boost::any_cast<std::string>(mwe->data));
    if(mwe->name == "handleThreadsChanged")
      handleThreadsChanged(boost::any_cast<std::string>(mwe->data));
    if(mwe->name == "handleStateChanged")
      handleStateChanged(boost::any_cast<std::string>(mwe->data));
    if(mwe->name == "handleGeomStateChanged")
      handleGeomStateChanged(boost::any_cast<std::string>(mwe->data));
    if(mwe->name == "updateGL" && !_glUpdated)
      {
        update();
        _glUpdated = true;
      }
  }

  void VolumeViewer::copyBoundingBoxToDataMap()
  {
    VolMagick::BoundingBox box;
    switch(volumeRenderingType())
      {
      case RGBA:
	box = _rgbaSubVolumeSelector;
	break;
      case ColorMapped:
	box = _cmSubVolumeSelector;
	break;
      }

    //copy bounding box to data map if property is set
    if(cvcapp.hasProperty(getObjectName("subvolume_box_data")))
      cvcapp.data(cvcapp.properties(getObjectName("subvolume_box_data")),
		  box);
  }

  float VolumeViewer::getWireframeLineWidth()
  {
    float lw = 1.2;
    if (cvcapp.hasProperty(getObjectName("geometry_line_width"))) {
      lw = cvcapp.properties<float>(getObjectName("geometry_line_width"));
    }
    else if (cvcapp.hasProperty("viewers.geometry_line_width")) {
      lw = cvcapp.properties<float>("viewers.geometry_line_width");
    }
    else {
      cvcapp.properties<float>(getObjectName("geometry_line_width"), lw);
    }
    return lw;
  }

  void VolumeViewer::propertiesChanged(const std::string& key)
  {
    QCoreApplication::postEvent(this,
      new CVCEvent("handlePropertiesChanged",key)
    );
  }

  void VolumeViewer::dataChanged(const std::string& key)
  {
    QCoreApplication::postEvent(this,
      new CVCEvent("handleDataChanged",key)
    );
  }

  void VolumeViewer::threadsChanged(const std::string& key)
  {
    QCoreApplication::postEvent(this,
      new CVCEvent("handleThreadsChanged",key)
    );
  }

  void VolumeViewer::stateChanged(const std::string& key)
  {
    QCoreApplication::postEvent(this,
      new CVCEvent("handleStateChanged",key)
    );
  }

  void VolumeViewer::geomStateChanged(const std::string& key)
  {
    QCoreApplication::postEvent(this,
      new CVCEvent("handleGeomStateChanged",key)
    );
  }

  // 12/02/2011 -- transfix -- setting default colortable if transfer function is 'none'
  void VolumeViewer::handlePropertiesChanged(const std::string& key)
  {
    using namespace std;
    using namespace boost;
    using namespace boost::algorithm;

    static log4cplus::Logger logger = log4cplus::Logger::getInstance("VolumeViewer.handlePropertiesChanged");

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

    try
      {
        vector<string> key_idents;
        split(key_idents, key, is_any_of("."));
        //check if the key was meant for this viewer
        if(key_idents.size() == 2 && key_idents[0] == getObjectName())
          {
            cvcapp.log(5,str(boost::format("%s :: viewer %s, prop %s\n")
                             % BOOST_CURRENT_FUNCTION
                             % key_idents[0]
                             % key_idents[1]));
	   

	    if(key_idents[1] == "take_snapshot") {
	      if(cvcapp.properties(key) == "true") {		
		
		string tmp = cvcapp.properties(key_idents[0]+".snapshot_filename");
		setSnapshotFileName(tmp.c_str());
		setSnapshotFormat("png");
		saveSnapshot(QString::fromStdString(tmp), true);
	      }
	    }



            if(key_idents[1] == "rendering_mode")
              {
                VolumeRenderingType t = ColorMapped;
                if(cvcapp.properties(key) == "colormapped")
                  t = ColorMapped;
                else if(cvcapp.properties(key) == "rgba" ||
                        cvcapp.properties(key) == "RGBA")
                  t = RGBA;
                else
                  cvcapp.log(1,str(boost::format("%s :: Error: invalid rendering_mode %s\n")
				   % BOOST_CURRENT_FUNCTION
                                   % cvcapp.properties(key)));
                volumeRenderingType(t);
              }
            else if(key_idents[1] == "shaded_rendering_enabled")
              {
                bool val = false;
                if(cvcapp.properties(key) == "true")
                  val = true;
                else if(cvcapp.properties(key) == "false")
                  val = false;
                else
                  cvcapp.log(1,str(boost::format("%s :: Error: invalid shaded_rendering_enabled\n")
				   % BOOST_CURRENT_FUNCTION));              

                setShadedRenderingMode(val);
              }
            if(key_idents[1] == "geometry_rendering_mode")
	      {

		//scene_geometry_t::render_mode_t rm;
		int rm = 0;

		if(cvcapp.properties(key) == "wireframe") {
		  rm=1;
		} else if(cvcapp.properties(key) == "filled_wireframe") {
		  rm=2;
		} else if(cvcapp.properties(key) == "flat") {
		  rm=3;
		} else if(cvcapp.properties(key) == "flat_wireframe") {
		  rm=4;
		} 
		// loop through the geometries and set rendering mode...
		for (const auto& mypair : _geometries) {
                  if(cvcapp.hasProperty(mypair.first+".render_mode"))
                    {
                      cvcapp.log(3,str(boost::format("%s :: %s has a render mode set, ignoring\n")
                                       % BOOST_CURRENT_FUNCTION
                                       % mypair.first));
                      continue;
                    }

		  if (mypair.second->render_mode == scene_geometry_t::TRIANGLES || 
		      mypair.second->render_mode == scene_geometry_t::TRIANGLE_WIREFRAME || 
		      mypair.second->render_mode == scene_geometry_t::TRIANGLE_FILLED_WIRE ||
		      mypair.second->render_mode == scene_geometry_t::TRIANGLES_FLAT || 
		      mypair.second->render_mode == scene_geometry_t::TRIANGLE_FLAT_FILLED_WIRE ) {
		    if (rm == 0) {
		      mypair.second->render_mode = scene_geometry_t::TRIANGLES;
		    } else if (rm == 1) {
		      mypair.second->render_mode = scene_geometry_t::TRIANGLE_WIREFRAME;
		    } else if (rm == 2) {
		      mypair.second->render_mode = scene_geometry_t::TRIANGLE_FILLED_WIRE;
		    } else if (rm == 3) {
		      // arand, 9-12-2011: adding "flat" rendering mode for triangles
		      mypair.second->render_mode = scene_geometry_t::TRIANGLES_FLAT;
		    } else if (rm == 4) {
		      mypair.second->render_mode = scene_geometry_t::TRIANGLE_FLAT_FILLED_WIRE;
		    }
		  } else if (mypair.second->render_mode == scene_geometry_t::QUADS || 
			     mypair.second->render_mode == scene_geometry_t::QUAD_WIREFRAME || 
			     mypair.second->render_mode == scene_geometry_t::QUAD_FILLED_WIRE) {
		    // arand TODO: implement "flat" rendering for quads

		    if (rm == 0) {
		      mypair.second->render_mode = scene_geometry_t::QUADS;
		    } else if (rm == 1) {
		      mypair.second->render_mode = scene_geometry_t::QUAD_WIREFRAME;
		    } else if (rm == 2) {
		      mypair.second->render_mode = scene_geometry_t::QUAD_FILLED_WIRE;
		    }
		  }
		  // WARNING: I am ignoring hex/tet meshes with this property...
		}
	      }
            else if(key_idents[1] == "draw_bounding_box")
              {
                bool val = false;
                if(cvcapp.properties(key) == "true")
                  val = true;
                else if(cvcapp.properties(key) == "false")
                  val = false;
                else
                  cvcapp.log(1,str(boost::format("%s :: Error: invalid draw_bounding_box\n")
				   % BOOST_CURRENT_FUNCTION));              

                drawBoundingBox(val);
              }
            else if(key_idents[1] == "draw_subvolume_selector")
              {
                bool val = false;
                if(cvcapp.properties(key) == "true")
                  val = true;
                else if(cvcapp.properties(key) == "false")
                  val = false;
                else
                  cvcapp.log(1,str(boost::format("%s :: Error: invalid draw_subvolume_selector\n")
				   % BOOST_CURRENT_FUNCTION));              

                drawSubVolumeSelector(val);
              }
            else if(key_idents[1] == "volume_rendering_quality")
              {
                double val = cvcapp.properties<double>(key);
                quality(val);
              }
            else if(key_idents[1] == "volume_rendering_near_plane")
              {
                double val = cvcapp.properties<double>(key);
                nearPlane(val);
              }
            else if(key_idents[1] == "draw_corner_axis")
              {
                bool val = false;
                if(cvcapp.properties(key) == "true")
                  val = true;
                else if(cvcapp.properties(key) == "false")
                  val = false;
                else
                  cvcapp.log(1,str(boost::format("%s :: Error: invalid draw_corner_axis\n")
				   % BOOST_CURRENT_FUNCTION));              

                drawCornerAxis(val);
              }
            else if(key_idents[1] == "projection_mode")
              {
                if(cvcapp.properties(key) == "orthographic")
		  camera()->setType(qglviewer::Camera::ORTHOGRAPHIC);
                else if(cvcapp.properties(key) == "perspective")
		  camera()->setType(qglviewer::Camera::PERSPECTIVE);
                else
                  cvcapp.log(1,str(boost::format("%s :: Error: invalid projection_mode\n")
				   % BOOST_CURRENT_FUNCTION));              
              }

            else if(key_idents[1] == "draw_geometry")
              {
                bool val = false;
                if(cvcapp.properties(key) == "true")
                  val = true;
                else if(cvcapp.properties(key) == "false")
                  val = false;
                else
                  cvcapp.log(1,str(boost::format("%s :: Error: invalid draw_geometry\n")
				   % BOOST_CURRENT_FUNCTION));              

                drawGeometry(val);
              }
            else if(key_idents[1] == "draw_volumes")
              {
                bool val = false;
                if(cvcapp.properties(key) == "true")
                  val = true;		       
                else if(cvcapp.properties(key) == "false")
                  val = false;
                else
                  cvcapp.log(1,str(boost::format("%s :: Error: invalid draw_volumes\n")
				   % BOOST_CURRENT_FUNCTION));              

                drawVolumes(val);
              }
            else if(key_idents[1] == "clip_geometry")
              {
                bool val = false;
                if(cvcapp.properties(key) == "true")
                  val = true;
                else if(cvcapp.properties(key) == "false")
                  val = false;
                else
                  cvcapp.log(1,str(boost::format("%s :: Error: invalid clip_geometry\n")
				   % BOOST_CURRENT_FUNCTION));              

                clipGeometryToBoundingBox(val);
              }
            else if(key_idents[1] == "draw_geometry_normals")
              {
                bool val = false;
                if(cvcapp.properties(key) == "true")
                  val = true;
                else if(cvcapp.properties(key) == "false")
                  val = false;
                else
                  cvcapp.log(1,str(boost::format("%s :: Error: invalid draw_geometry_normals\n")
				   % BOOST_CURRENT_FUNCTION));              

                drawGeometryNormals(val);
              }
            //---------------------------------------------------------
            // Multi-Tile Server properties
            else if(key_idents[1] == "syncCamera_with_multiTileServer")
              {
                bool val = false;
                if(cvcapp.properties(key) == "true")
                  val = true;
                else if(cvcapp.properties(key) == "false")
                  val = false;
                else
                  cvcapp.log(1,str(boost::format("%s :: Error: invalid syncCamera_with_multiTileServer\n")
				   % BOOST_CURRENT_FUNCTION));              

                syncCameraWithMultiTileServer(val);
                LOG4CPLUS_TRACE(logger, "syncCamera_with_multiTileServer changed");
              }
            else if(key_idents[1] == "syncTransferFunc_with_multiTileServer")
              {
                bool val = false;
                if(cvcapp.properties(key) == "true")
                  val = true;
                else if(cvcapp.properties(key) == "false")
                  val = false;
                else
                  cvcapp.log(1,str(boost::format("%s :: Error: invalid syncTransferFunc_with_multiTileServer\n")
				   % BOOST_CURRENT_FUNCTION));              

                syncTransferFuncWithMultiTileServer(val);
                LOG4CPLUS_TRACE(logger, "syncTransferFunc_with_multiTileServer changed");
              }
            else if(key_idents[1] == "syncShadedRender_with_multiTileServer")
              {
                bool val = false;
                if(cvcapp.properties(key) == "true")
                  val = true;
                else if(cvcapp.properties(key) == "false")
                  val = false;
                else
                  cvcapp.log(1,str(boost::format("%s :: Error: invalid syncShadedRender_with_multiTileServer\n")
				   % BOOST_CURRENT_FUNCTION));              

                syncShadedRenderWithMultiTileServer(val);
                LOG4CPLUS_TRACE(logger, "syncShadedRender_with_multiTileServer changed");
              }
            else if(key_idents[1] == "syncRenderMode_with_multiTileServer")
              {
                bool val = false;
                if(cvcapp.properties(key) == "true")
                  val = true;
                else if(cvcapp.properties(key) == "false")
                  val = false;
                else
                  cvcapp.log(1,str(boost::format("%s :: Error: invalid syncRenderMode_with_multiTileServer\n")
				   % BOOST_CURRENT_FUNCTION));              

                syncRenderModeWithMultiTileServer(val);
                LOG4CPLUS_TRACE(logger, "syncRenderMode_with_multiTileServer changed");
              }
            else if(key_idents[1] == "syncInitial_multiTileServer")
              {
                if( cvcapp.properties<int>(key) > 0 ) {
                   initializeMultiTileServer();
                   fprintf( stderr, "syncInitial_multiTileServer\n");
                }
              }
            else if(key_idents[1] == "interactiveMode_with_multiTileServer")
              {
                bool val = false;
                if(cvcapp.properties(key) == "true")
                  val = true;
                else if(cvcapp.properties(key) == "false")
                  val = false;
                else
                  cvcapp.log(1,str(boost::format("%s :: Error: invalid interactiveMode_with_multiTileServer\n")
				   % BOOST_CURRENT_FUNCTION));

                interactiveSyncWithMultiTileServer( val );
                LOG4CPLUS_TRACE(logger, "interactiveMode_with_multiTileServer changed");
              }
            else if(key_idents[1] == "syncMode_with_multiTileServer")
              {
                if( cvcapp.properties<int>(key) > 0 ) {
                   syncCurrentWithMultiTileServer();
                   fprintf( stderr, "syncCurrent_multiTileServer\n");
                }
              }
            //---------------------------------------------------------
            else if(key_idents[1] == "io_distance")
              {
                camera()->setIODistance(cvcapp.properties<float>(key));
              }
            else if(key_idents[1] == "physical_distance_to_screen")
              {
                camera()->setPhysicalDistanceToScreen(cvcapp.properties<float>(key));
              }
            else if(key_idents[1] == "physical_screen_width")
              {
                camera()->setPhysicalScreenWidth(cvcapp.properties<float>(key));
              }
            else if(key_idents[1] == "focus_distance")
              {
                camera()->setFocusDistance(cvcapp.properties<float>(key));
              }
            else if(key_idents[1] == "background_color")
              {
                string valstring =
                  cvcapp.properties(key);

		QColor color(QString::fromStdString(valstring));
                if(color.isValid())
                  {
                    LOG4CPLUS_TRACE(logger, "setting background color to " << color.name().toStdString());
                    // string msg = str(boost::format("%s :: setting background color to %s")
                    //                  % BOOST_CURRENT_FUNCTION
                    //                  % color.name().toStdString());
                    // cvcapp.log(3,msg+"\n");
                    if(hasBeenInitialized()) makeCurrent();
                    setBackgroundColor(color);
                  }
                else
                  {
                    vector<string> color_components =
                      cvcapp.listProperty(key);
                    QColor color;
                    qreal comp[4] = { 0.0, 0.0, 0.0, 1.0 };
                    vector<qreal> components;
                    for (const auto& val : color_components)
                      components.push_back(lexical_cast<qreal>(val));
                    int i=0;
                    for (const auto& val : components)
                      if(i < 4) comp[i++] = val;
                    color = QColor::fromRgbF(comp[0],comp[1],comp[2],comp[3]);
                  }
              }
            else if(key_idents[1] == "volumes")
              {
                vector<VolMagick::Volume> vols = cvcapp.propertyData<VolMagick::Volume>(key,true);
                if(!vols.empty())
                  {
                    colorMappedVolume(vols[0]);
                    rgbaVolumes(vols);
                  }
              }
            else if(key_idents[1] == "geometries")
              {
                vector<string> geo_names = cvcapp.listProperty(key,true);

                //Tell the viewer to update VBOs for any new geometry
                _vboUpdated = false;

                //Get the set of geometries to remove from the internal map
                set<string> geo_names_set;
                copy(geo_names.begin(), geo_names.end(), 
                     inserter(geo_names_set, geo_names_set.begin()));

                set<string> orig_geo_names_set;
                for (const auto& val : _geometries)
                  orig_geo_names_set.insert(val.first);
                  
                set<string> remove_names;
                for (const auto& name : orig_geo_names_set)
                  if(geo_names_set.find(name)==geo_names_set.end())
                    remove_names.insert(name);

                //Now remove the geometries from the internal map
                for (const auto& name : remove_names)
                  _geometries.erase(name);
               
                for (auto name : geo_names) {
                    //No need to add if already present
                    if(_geometries.find(name)!=_geometries.end())
                      continue;

		    // arand: fixed this, 4-12-2011
		    trim(name);

		    // pick the correct mode here...
		    /*
                      enum render_mode_t { POINTS, LINES,
                      TRIANGLES, TRIANGLE_WIREFRAME,
                      TRIANGLE_FILLED_WIRE, QUADS,
                      QUAD_WIREFRAME, QUAD_FILLED_WIRE,		       
                      TETRA, HEXA };
		    */

                  }
              }
            //check to see if the transfer function we are watching has changed
            else if(key_idents[1] == "transfer_function")
              {
                if(cvcapp.hasProperty(getObjectName("transfer_function")) &&
                   cvcapp.properties(getObjectName("transfer_function"))!="none")
                  {
#ifdef VOLUMEROVER2_FLOAT_COLORTABLE
                    shared_array<float> table =
                      cvcapp.data<shared_array<float> >(cvcapp.properties(getObjectName("transfer_function")));
#else
                    shared_array<unsigned char> table =
                      cvcapp.data<shared_array<unsigned char> >(cvcapp.properties(getObjectName("transfer_function")));
#endif
                    colorTable(table.get());
                  }
                else
                  {
                    setDefaultColorTable();
                  }
                
                //for multiTileServer
                if(interactiveSyncWithMultiTileServer() && syncTransferFuncWithMultiTileServer())
                   #ifdef USE_XmlRpc
                   needToUpdateMultiTileServer( -1, true );
                   #else
                   syncColorTable();
                   #endif
              }
            else if(key_idents[1] == "orientation")
              {
                if(!_mouseMoved)
                  {
                    string quaternion = cvcapp.properties(getObjectName("orientation"));
                    stringstream ss(quaternion);
                    double q0,q1,q2,q3;
                    ss >> q0 >> q1 >> q2 >> q3;
                    cvcapp.log(4,str(boost::format("%1% :: Q: %2% %3% %4% %5%\n")
                                     % BOOST_CURRENT_FUNCTION
                                     % q0 % q1 % q2 % q3));
                    camera()->setOrientation(qglviewer::Quaternion(q0,q1,q2,q3));
                    scheduleUpdateGL();
                  }
              }
            else if(key_idents[1] == "position")
              {
                if(!_mouseMoved)
                  {
                    string position = cvcapp.properties(getObjectName("position"));
                    vector<string> posstrs;
                    split(posstrs, position, is_any_of(","));
                    posstrs.resize(3); //make sure we have 3
                    double 
                      x = lexical_cast<double>(posstrs[0]),
                      y = lexical_cast<double>(posstrs[1]),
                      z = lexical_cast<double>(posstrs[2]);
                    cvcapp.log(4,str(boost::format("%1% :: pos: %2% %3% %4%\n")
                                     % BOOST_CURRENT_FUNCTION
                                     % x % y % z));
                    camera()->setPosition(qglviewer::Vec(x,y,z));
                    scheduleUpdateGL();
                  }
              }
            else if(key_idents[1] == "fov")
              {
                camera()->setFieldOfView(cvcapp.properties<float>(getObjectName("fov")));
              }
          }
        else
          {
            //Check if the property is for any of the geometries we care about
            std::vector<std::string> geo_names =
              cvcapp.listProperty(getObjectName("geometries"),true);
            if(key_idents.size()>1)
              for (const auto& name : geo_names)
                if(name == key_idents[0] && _geometries[name])
                  {
                    if(key_idents[1] == "render_mode")
                      {
                        scene_geometry_ptr scene_geom_ptr = _geometries[name];
                        scene_geom_ptr->render_mode = 
                          scene_geometry_t::render_mode_enum(cvcapp.properties(key));
                            
                        //generate the render_geometry for this geometry if it doesn't
                        //exist, and we're rendering TETRA or HEXA - Joe, 4-19-2011
                        if(scene_geom_ptr->render_geometry.empty() &&
                           (scene_geom_ptr->render_mode == scene_geometry_t::TETRA ||
                            scene_geom_ptr->render_mode == scene_geometry_t::HEXA))
                          {
                            cvcapp.log(3,str(boost::format("%s :: generating wire interior for geometry '%s'\n")
                                             % BOOST_CURRENT_FUNCTION
                                             % name));
                            scene_geom_ptr->render_geometry =
                              scene_geom_ptr->geometry.generate_wire_interior();
                          }
                            
                        scheduleUpdateGL();
                      }
                    else if(key_idents[1] == "visible")
                      {
                        scene_geometry_ptr scene_geom_ptr = _geometries[name];
                        scene_geom_ptr->visible = cvcapp.properties(key) != "false";
                        scheduleUpdateGL();
                      }
                    break;
                  }
          }
      }
    catch(std::exception& e)
      {
        string msg = str(boost::format("%s :: Error: %s")
                         % BOOST_CURRENT_FUNCTION
                         % e.what());
        cvcapp.log(1,msg+"\n");
        QMessageBox::critical(this,"Error",
                              QString("Error setting value for key %1 : %2")
                              .arg(QString::fromStdString(key))
                              .arg(QString::fromStdString(msg)));
      }
  }

  void VolumeViewer::handleDataChanged(const std::string& key)
  {
    using namespace std;
    using namespace boost;
    using namespace boost::algorithm;

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
	//check to see if the transfer function we are watching has changed
	if(cvcapp.properties(getObjectName("transfer_function")) == key)
	  {
	    cvcapp.log(4,str(boost::format("%s :: viewer %s, data %s\n")
			     % BOOST_CURRENT_FUNCTION
			     % getObjectName()
			     % key));

#ifdef VOLUMEROVER2_FLOAT_COLORTABLE
            shared_array<float> table =
              cvcapp.data<shared_array<float> >(key);
#else
            shared_array<unsigned char> table =
              cvcapp.data<shared_array<unsigned char> >(key);
#endif
	    colorTable(table.get());
            
            //for multiTileServer
            if(interactiveSyncWithMultiTileServer() && syncTransferFuncWithMultiTileServer())
               syncColorTable();
	  }

	//check if any of the volumes we care about changed
	{
	  bool volChanged = false;
	  vector<string> vol_names =
            cvcapp.listProperty(getObjectName("volumes"));
	  vector<VolMagick::Volume> vols;
	  for (const auto& name : vol_names) {
	      if(cvcapp.isData<VolMagick::Volume>(name))
		{
		  vols.push_back(cvcapp.data<VolMagick::Volume>(name));
		  if(name == key)
		    {
		      volChanged = true;
		      break;
		    }
		}
	    }
	  if(!vols.empty() && volChanged)
	    {
	      cvcapp.log(4,boost::str(boost::format("%s :: viewer %s, data %s\n")
                                      % BOOST_CURRENT_FUNCTION
                                      % getObjectName()
                                      % key));

	      colorMappedVolume(vols[0]);
	      rgbaVolumes(vols);

#if 0
              if(cvcapp.properties(getObjectName("rendering_mode"))=="colormapped")
                uploadColorMappedVolume();
              else if(cvcapp.properties(getObjectName("rendering_mode"))=="rgba" ||
                      cvcapp.properties(getObjectName("rendering_mode"))=="RGBA")
                uploadRGBAVolume();
#endif
            }
	}

	//check if any of the geometries we care about changed
	{
	  vector<string> geo_names =
            cvcapp.listProperty(getObjectName("geometries"));
	  for (auto name : geo_names) {
              trim(name);
              if(name == key)
                {
                  cvcapp.log(5,boost::str(boost::format("%s :: Checking name = %s and key = %s\n")
                                          % BOOST_CURRENT_FUNCTION
                                          % name 
                                          % key));
                  if(cvcapp.isData<cvcraw_geometry::cvcgeom_t>(key))
                    {
                      cvcapp.log(4,boost::str(boost::format("%s :: viewer %s, data %s\n")
                                              % BOOST_CURRENT_FUNCTION
                                              % getObjectName()
                                              % key));

		      using namespace std;		      
		      cvcraw_geometry::cvcgeom_t myGeom = cvcapp.data<cvcraw_geometry::cvcgeom_t>(key);
		      int np = myGeom.const_points().size();
		      int nl = myGeom.const_lines().size();
		      int nt = myGeom.const_triangles().size();
		      int nq = myGeom.const_quads().size();

		      if (nl == 0 && nt == 0 && nq == 0) {
			addGeometry(myGeom, key, scene_geometry_t::POINTS);
		      } else if (nt == 0 && nq == 0) {
			addGeometry(myGeom, key, scene_geometry_t::LINES);
		      } else if (nt == 0) {
			addGeometry(myGeom, key, scene_geometry_t::QUADS);
		      } else {		      
			addGeometry(myGeom, key);
		      }
		      
		      /*
                      addGeometry(
				  cvcapp.data<cvcraw_geometry::cvcgeom_t>(key), key
				  );
		      */
				  
				  
                      break;
                    }
#ifdef USING_TILING
                  else if(cvcapp.isData<cvcraw_geometry::contours_t>(key)) {
                    cvcapp.log(5,str(boost::format("%s :: Adding %s as a contour dataset")
                                     % BOOST_CURRENT_FUNCTION % name));
		    using namespace std;		      
		    cvcraw_geometry::cvcgeom_t myGeom = cvcapp.data<cvcraw_geometry::contours_t>(key).geom();
		    addGeometry(myGeom, key, scene_geometry_t::LINES);
		  }
#endif
                }
            }
	}
      }
    catch(std::exception& e)
      {
	string msg = str(boost::format("%s :: Error: %s")
			 % BOOST_CURRENT_FUNCTION
			 % e.what());
	cvcapp.log(1,msg+"\n");
	QMessageBox::critical(this,"Error",
			      QString("Error setting value for key %1 : %2")
			      .arg(QString::fromStdString(key))
			      .arg(QString::fromStdString(msg)));
      }
  }

  // Purpose: to provide thread feedback to the user in the volume viewers
  // 09/18/2011 - transfix - initial implementation.
  void VolumeViewer::handleThreadsChanged(const std::string& key)
  {
    using namespace std;
    using namespace boost;

    //ignore empty keys.  TODO: find out where they're coming from...
    if(key.empty()) return;

    cvcapp.log(5,str(boost::format("%s :: key == %s\n")
                     % BOOST_CURRENT_FUNCTION
                     % key));

    ThreadPtr tp = cvcapp.threads(key);
    std::stringstream ss;
    if(tp)
      ss << tp->get_id();
    else
      ss << "null";
    string tid_str = ss.str();

    double progress = cvcapp.threadProgress(key)*100.0;

    string message = str(boost::format("%1% :: %2% :: %3% :: %4%")
                         % key
                         % tid_str
                         % cvcapp.threadInfo(key)
                         % progress);

    //Not calling displayMessage() here directly because this function
    //gets called very rapidly sometimes and paint would slow everything down.
    _displayMessage = QString::fromStdString(message);
    _doDisplayMessage = true;
  }

  //03/14/2012 -- transfix -- creation
  void VolumeViewer::handleStateChanged(const std::string& childState)
  {
    using namespace std;
    using namespace boost;
    
    static log4cplus::Logger logger = log4cplus::Logger::getInstance("VolumeViewer.handleStateChanged");

    cvcapp.log(5,
               str(
                   boost::format("%s :: state: %s\n")
                   % BOOST_CURRENT_FUNCTION
                   % childState));

    std::string val = cvcstate(stateName(childState)).value();

    try
      {
        if(childState == "take_snapshot" &&
           val == "true")
          {
            string tmp = cvcstate(stateName("snapshot_filename")).value();
            setSnapshotFileName(tmp.c_str());
            setSnapshotFormat("png");
            saveSnapshot(QString::fromStdString(tmp), true);
            //TODO: set back to false?
          }

        else if(childState == "rendering_mode")
          {
            VolumeRenderingType t = ColorMapped;
            if(val == "colormapped")
              t = ColorMapped;
            else if(val == "rgba" ||
                    val == "RGBA")
              t = RGBA;
            else
              cvcapp.log(1,str(boost::format("%s :: Error: invalid rendering_mode %s\n")
                               % BOOST_CURRENT_FUNCTION
                               % val));
            volumeRenderingType(t);
          }

        else if(childState == "shaded_rendering_enabled")
          {
            bool flag = false;
            if(val == "true")
              flag = true;
            else if(val == "false")
              flag = false;
            else
              cvcapp.log(1,str(boost::format("%s :: Error: invalid shaded_rendering_enabled\n")
                               % BOOST_CURRENT_FUNCTION));              
            setShadedRenderingMode(flag);
          }

        else if(childState == "geometry_rendering_mode")
          {
            int rm = 0;

            if(val == "wireframe") {
              rm=1;
            } else if(val == "filled_wireframe") {
              rm=2;
            } else if(val == "flat") {
              rm=3;
            } else if(val == "flat_wireframe") {
              rm=4;
            } 
            // loop through the geometries and set rendering mode...
            for (const auto& mypair : _geometries) {
              if(cvcstate(mypair.first+".render_mode").initialized())
                {
                  cvcapp.log(3,str(boost::format("%s :: %s has a render mode set, ignoring\n")
                                   % BOOST_CURRENT_FUNCTION
                                   % mypair.first));
                  continue;
                }

              if (mypair.second->render_mode == scene_geometry_t::TRIANGLES || 
                  mypair.second->render_mode == scene_geometry_t::TRIANGLE_WIREFRAME || 
                  mypair.second->render_mode == scene_geometry_t::TRIANGLE_FILLED_WIRE ||
                  mypair.second->render_mode == scene_geometry_t::TRIANGLES_FLAT || 
                  mypair.second->render_mode == scene_geometry_t::TRIANGLE_FLAT_FILLED_WIRE ) {
                if (rm == 0) {
                  mypair.second->render_mode = scene_geometry_t::TRIANGLES;
                } else if (rm == 1) {
                  mypair.second->render_mode = scene_geometry_t::TRIANGLE_WIREFRAME;
                } else if (rm == 2) {
                  mypair.second->render_mode = scene_geometry_t::TRIANGLE_FILLED_WIRE;
                } else if (rm == 3) {
                  // arand, 9-12-2011: adding "flat" rendering mode for triangles
                  mypair.second->render_mode = scene_geometry_t::TRIANGLES_FLAT;
                } else if (rm == 4) {
                  mypair.second->render_mode = scene_geometry_t::TRIANGLE_FLAT_FILLED_WIRE;
                }
              } else if (mypair.second->render_mode == scene_geometry_t::QUADS || 
                         mypair.second->render_mode == scene_geometry_t::QUAD_WIREFRAME || 
                         mypair.second->render_mode == scene_geometry_t::QUAD_FILLED_WIRE) {
                // arand TODO: implement "flat" rendering for quads

                if (rm == 0) {
                  mypair.second->render_mode = scene_geometry_t::QUADS;
                } else if (rm == 1) {
                  mypair.second->render_mode = scene_geometry_t::QUAD_WIREFRAME;
                } else if (rm == 2) {
                  mypair.second->render_mode = scene_geometry_t::QUAD_FILLED_WIRE;
                }
              }
              // WARNING: I am ignoring hex/tet meshes with this property...
            }
          }

        else if(childState == "draw_bounding_box")
          {
            bool flag = false;
            if(val == "true")
              flag = true;
            else if(val == "false")
              flag = false;
            else
              cvcapp.log(1,str(boost::format("%s :: Error: invalid draw_bounding_box: %s\n")
                               % BOOST_CURRENT_FUNCTION
                               % flag));
            drawBoundingBox(flag);
          }

        else if(childState == "draw_subvolume_selector")
          {
            bool flag = false;
            if(val == "true")
              flag = true;
            else if(val == "false")
              flag = false;
            else
              cvcapp.log(1,str(boost::format("%s :: Error: invalid draw_subvolume_selector\n")
                               % BOOST_CURRENT_FUNCTION));
            drawSubVolumeSelector(flag); 
          }

        else if(childState == "volume_rendering_quality")
          {
            double val = cvcstate(stateName(childState)).value<double>();
            quality(val);
          }

        else if(childState == "volume_rendering_near_plane")
          {
            double val = cvcstate(stateName(childState)).value<double>();
            nearPlane(val);
          }

        else if(childState == "draw_corner_axis")
          {
            bool flag = false;
            if(val == "true")
              flag = true;
            else if(val == "false")
              flag = false;
            else
              cvcapp.log(1,str(boost::format("%s :: Error: invalid draw_corner_axis\n")
                               % BOOST_CURRENT_FUNCTION));              
            drawCornerAxis(flag);
          }

        else if(childState == "projection_mode")
          {
            if(val == "orthographic")
              camera()->setType(qglviewer::Camera::ORTHOGRAPHIC);
            else if(val == "perspective")
              camera()->setType(qglviewer::Camera::PERSPECTIVE);
            else
              cvcapp.log(1,str(boost::format("%s :: Error: invalid projection_mode\n")
                               % BOOST_CURRENT_FUNCTION));              
          }

        else if(childState == "draw_geometry")
          {
            bool flag = false;
            if(val == "true")
              flag = true;
            else if(val == "false")
              flag = false;
            else
              cvcapp.log(1,str(boost::format("%s :: Error: invalid draw_geometry\n")
                               % BOOST_CURRENT_FUNCTION));              
            drawGeometry(flag);
          }

        else if(childState == "draw_volumes")
          {
            bool flag = false;
            if(val == "true")
              flag = true;		       
            else if(val == "false")
              flag = false;
            else
              cvcapp.log(1,str(boost::format("%s :: Error: invalid draw_volumes\n")
                               % BOOST_CURRENT_FUNCTION));              
            drawVolumes(flag);
          }

        else if(childState == "clip_geometry")
          {
            bool flag = false;
            if(val == "true")
              flag = true;		       
            else if(val == "false")
              flag = false;
            else
              cvcapp.log(1,str(boost::format("%s :: Error: invalid clip_geometry\n")
                               % BOOST_CURRENT_FUNCTION));              
            clipGeometryToBoundingBox(flag);
          }

        else if(childState == "draw_geometry_normals")
          {
            bool flag = false;
            if(val == "true")
              flag = true;		       
            else if(val == "false")
              flag = false;
            else
              cvcapp.log(1,str(boost::format("%s :: Error: invalid draw_geometry_normals\n")
                               % BOOST_CURRENT_FUNCTION));              
            drawGeometryNormals(flag);
          }

        //---------------------------------------------------------
        // Multi-Tile Server properties
        else if(childState == "syncCamera_with_multiTileServer")
          {
            bool flag = false;
            if(val == "true")
              flag = true;
            else if(val == "false")
              flag = false;
            else
              cvcapp.log(1,str(boost::format("%s :: Error: invalid syncCamera_with_multiTileServer\n")
                               % BOOST_CURRENT_FUNCTION));              
            syncCameraWithMultiTileServer(flag);
            LOG4CPLUS_TRACE(logger, "syncCamera_with_multiTileServer changed");
          }
        else if(childState == "syncTransferFunc_with_multiTileServer")
          {
            bool flag = false;
            if(val == "true")
              flag = true;
            else if(val == "false")
              flag = false;
            else
              cvcapp.log(1,str(boost::format("%s :: Error: invalid syncTransferFunc_with_multiTileServer\n")
                               % BOOST_CURRENT_FUNCTION));              
            syncTransferFuncWithMultiTileServer(flag);
            LOG4CPLUS_TRACE(logger, "syncTransferFunc_with_multiTileServer changed");
          }
        else if(childState == "syncShadedRender_with_multiTileServer")
          {
            bool flag = false;
            if(val == "true")
              flag = true;
            else if(val == "false")
              flag = false;
            else
              cvcapp.log(1,str(boost::format("%s :: Error: invalid syncShadedRender_with_multiTileServer\n")
                               % BOOST_CURRENT_FUNCTION));              
            syncShadedRenderWithMultiTileServer(flag);
            LOG4CPLUS_TRACE(logger, "syncShadedRender_with_multiTileServer changed");
          }
        else if(childState == "syncRenderMode_with_multiTileServer")
          {
            bool flag = false;
            if(val == "true")
              flag = true;
            else if(val == "false")
              flag = false;
            else
              cvcapp.log(1,str(boost::format("%s :: Error: invalid syncRenderMode_with_multiTileServer\n")
                               % BOOST_CURRENT_FUNCTION));              
            syncRenderModeWithMultiTileServer(flag);
            LOG4CPLUS_TRACE(logger, "syncRenderMode_with_multiTileServer changed");
          }
        else if(childState == "syncInitial_multiTileServer")
          {
            if( cvcstate(stateName(childState)).value<int>() > 0 ) {
              initializeMultiTileServer();
              fprintf( stderr, "syncInitial_multiTileServer\n");
            }
          }
        else if(childState == "interactiveMode_with_multiTileServer")
          {
            bool flag = false;
            if(val == "true")
              flag = true;
            else if(val == "false")
              flag = false;
            else
              cvcapp.log(1,str(boost::format("%s :: Error: invalid interactiveMode_with_multiTileServer\n")
                               % BOOST_CURRENT_FUNCTION));
            interactiveSyncWithMultiTileServer( flag );
            LOG4CPLUS_TRACE(logger, "interactiveMode_with_multiTileServer changed");
          }
        else if(childState == "syncMode_with_multiTileServer")
          {
            if( cvcstate(stateName(childState)).value<int>() > 0 ) {
              syncCurrentWithMultiTileServer();
              fprintf( stderr, "syncCurrent_multiTileServer\n");
            }
          }
        //---------------------------------------------------------
        else if(childState == "io_distance")
          {
            camera()->setIODistance(cvcstate(stateName(childState)).value<float>());
          }
        else if(childState == "physical_distance_to_screen")
          {
            camera()->setPhysicalDistanceToScreen(cvcstate(stateName(childState)).value<float>());
          }
        else if(childState == "physical_screen_width")
          {
            camera()->setPhysicalScreenWidth(cvcstate(stateName(childState)).value<float>());
          }
        else if(childState == "focus_distance")
          {
            camera()->setFocusDistance(cvcstate(stateName(childState)).value<float>());
          }
        else if(childState == "background_color")
          {
            QColor color(QString::fromStdString(val));
            if(color.isValid())
              {
                LOG4CPLUS_TRACE(logger, "setting background color to " << color.name().toStdString());
                // string msg = str(boost::format("%s :: setting background color to %s")
                //                  % BOOST_CURRENT_FUNCTION
                //                  % color.name().toStdString());
                // cvcapp.log(3,msg+"\n");
                if(hasBeenInitialized()) makeCurrent();
                setBackgroundColor(color);
              }
            else
              {
                vector<string> color_components =
                  cvcstate(stateName(childState)).values();
                QColor color;
                qreal comp[4] = { 0.0, 0.0, 0.0, 1.0 };
                vector<qreal> components;
                for (const auto& val : color_components)
                  components.push_back(lexical_cast<qreal>(val));
                int i=0;
                for (const auto& val : components)
                  if(i < 4) comp[i++] = val;
                color = QColor::fromRgbF(comp[0],comp[1],comp[2],comp[3]);
              }
          }

        else if(childState == "volumes")
          {
            vector<VolMagick::Volume> vols = 
              cvcstate(stateName(childState)).valueData<VolMagick::Volume>();
            if(!vols.empty())
              {
                colorMappedVolume(vols[0]);
                rgbaVolumes(vols);
              }
          }

        else if(childState == "geometries")
          {
            vector<string> geo_names = cvcstate(stateName(childState)).values();

            //Tell the viewer to update VBOs for any new geometry
            _vboUpdated = false;

            //Get the set of geometries to remove from the internal map
            set<string> geo_names_set;
            copy(geo_names.begin(), geo_names.end(), 
                 inserter(geo_names_set, geo_names_set.begin()));

            set<string> orig_geo_names_set;
            for (const auto& val : _geometries)
              orig_geo_names_set.insert(val.first);
                  
            set<string> remove_names;
            for (const auto& name : orig_geo_names_set)
              if(geo_names_set.find(name)==geo_names_set.end())
                remove_names.insert(name);

            //Now remove the geometries from the internal map
            for (const auto& name : remove_names)
              _geometries.erase(name);
               
            for (auto name : geo_names) {
                //No need to add if already present
                if(_geometries.find(name)!=_geometries.end())
                  continue;

                // arand: fixed this, 4-12-2011
                trim(name);

                // pick the correct mode here...
                /*
                  enum render_mode_t { POINTS, LINES,
                  TRIANGLES, TRIANGLE_WIREFRAME,
                  TRIANGLE_FILLED_WIRE, QUADS,
                  QUAD_WIREFRAME, QUAD_FILLED_WIRE,		       
                  TETRA, HEXA };
                */

              }
          }

        else if(childState == "transfer_function")
          {
            if(val != "none")
              {
#ifdef VOLUMEROVER2_FLOAT_COLORTABLE
                shared_array<float> table =
                  cvcstate(val).data<shared_array<float> >();
#else
                shared_array<unsigned char> table =
                  cvcstate(val).data<shared_array<unsigned char> >();
#endif
                colorTable(table.get());
              }
            else
              {
                setDefaultColorTable();
              }
                
            //for multiTileServer
            if(interactiveSyncWithMultiTileServer() && syncTransferFuncWithMultiTileServer())
#ifdef USE_XmlRpc
              needToUpdateMultiTileServer( -1, true );
#else
            syncColorTable();
#endif
          }

        else if(childState == "orientation")
          {
            if(!_mouseMoved)
              {
                string quaternion = val;
                stringstream ss(quaternion);
                double q0,q1,q2,q3;
                ss >> q0 >> q1 >> q2 >> q3;
                cvcapp.log(4,str(boost::format("%1% :: Q: %2% %3% %4% %5%\n")
                                 % BOOST_CURRENT_FUNCTION
                                 % q0 % q1 % q2 % q3));
                camera()->setOrientation(qglviewer::Quaternion(q0,q1,q2,q3));
                scheduleUpdateGL();
              }
          }

        else if(childState == "position")
          {
            if(!_mouseMoved)
              {
                string position = val;
                vector<string> posstrs;
                split(posstrs, position, is_any_of(","));
                posstrs.resize(3); //make sure we have 3
                double 
                  x = lexical_cast<double>(posstrs[0]),
                  y = lexical_cast<double>(posstrs[1]),
                  z = lexical_cast<double>(posstrs[2]);
                cvcapp.log(4,str(boost::format("%1% :: pos: %2% %3% %4%\n")
                                 % BOOST_CURRENT_FUNCTION
                                 % x % y % z));
                camera()->setPosition(qglviewer::Vec(x,y,z));
                scheduleUpdateGL();
              }
          }

        else if(childState == "fov")
          {
            camera()->setFieldOfView(cvcstate(stateName(childState)).value<float>());
          }
      }
    catch(std::exception& e)
      {
        string msg = str(boost::format("%s :: Error: %s")
                         % BOOST_CURRENT_FUNCTION
                         % e.what());
        cvcapp.log(1,msg+"\n");
        QMessageBox::critical(this,"Error",
                              QString("Error setting value for key %1 : %2")
                              .arg(QString::fromStdString(stateName(childState)))
                              .arg(QString::fromStdString(msg)));
      }
  }

  void VolumeViewer::handleGeomStateChanged(const std::string& childState)
  {
    cvcapp.log(2,
               str(
                   boost::format("%s :: state: %s\n")
                   % BOOST_CURRENT_FUNCTION
                   % childState));
  }

  //03/14/2012 -- transfix -- using cvcstate.
  void VolumeViewer::copyCameraToPropertyMap()
  {
    //TODO: get rid of properties usage

    //set orientation property
    {
      std::stringstream ss;
      qglviewer::Quaternion Q = camera()->orientation();
      ss<< Q[0] << '\t' << Q[1] << '\t' << Q[2] << '\t' << Q[3];
      //ss << camera()->orientation();

      if(cvcapp.properties(getObjectName("orientation"))!=ss.str())
        cvcapp.properties(getObjectName("orientation"),ss.str());

      if(cvcstate(stateName("orientation")).value()!=ss.str())
        cvcstate(stateName("orientation")).value(ss.str());
    }
    
    //set camera position property
    {
      using namespace boost;
      using namespace std;

      //ss << camera()->position();

      string posstr = str(boost::format("%1%,%2%,%3%")
                          % lexical_cast<string>(camera()->position().x)
                          % lexical_cast<string>(camera()->position().y)
                          % lexical_cast<string>(camera()->position().z));

      cvcapp.log(5,
                 str(
                     boost::format("%s :: %s\n")
                     % BOOST_CURRENT_FUNCTION
                     % posstr));

      if(cvcapp.properties(getObjectName("position"))!=posstr)
        cvcapp.properties(getObjectName("position"),posstr);

      if(cvcstate(stateName("position")).value()!=posstr)
        cvcstate(stateName("position")).value(posstr);
    }
  }
}
