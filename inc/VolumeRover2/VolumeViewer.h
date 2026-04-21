/*
  Copyright 2011 The University of Texas at Austin
  
	Authors: Jose Rivera <transfix@ices.utexas.edu>
	Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of Volume Rover.

  Volume Rover is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  Volume Rover is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with iotree; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/* $Id$ */

#ifndef __VOLUMEVIEWER_H__
#define __VOLUMEVIEWER_H__

//*********************
// select to use XmlRpc, otherwise use socket
#define USE_XmlRpc
//*********************

// GLEW must be included before any OpenGL headers
#include <GL/glew.h>

#include <CVC/Namespace.h>
#include <CVC/App.h>
#include <CVC/State.h>

#include <cstring>
#include <string>
#include <assert.h>
#include <vector>
#include <boost/shared_array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>
#include <QGLViewer/qglviewer.h>
#include <QGLViewer/vec.h>
#include <VolMagick/VolMagick.h>
#include <VolumeRenderer/VolumeRenderer.h>

#include <cvcraw_geometry/cvcraw_geometry.h>
#include <cvcraw_geometry/cvcgeom.h>

#ifdef USE_XmlRpc
#include <xmlrpc/XmlRpc.h>
#else
#include <Socket/CVCSocket.h>
#endif

class QTimer;

namespace CVC_NAMESPACE
{
#ifdef USE_XmlRpc
  class XmlRpcThread;
#endif

  // 10/28/2011 - transfix - added copyCameraOnNormalize()
  // 12/02/2011 - transfix - added setDefaultColorTable()
  class VolumeViewer : public QGLViewer
  {
  Q_OBJECT

  public:

    enum VolumeRenderingType { ColorMapped, RGBA };

    explicit VolumeViewer(QWidget* parent=nullptr, 
		    Qt::WindowFlags flags=Qt::WindowFlags())
      : QGLViewer(parent, flags)
    { defaultConstructor(); }

    ~VolumeViewer();

    bool hasBeenInitialized() const { return _hasBeenInitialized; }

    void setDefaultColorTable();

    void colorMappedVolume(const VolMagick::Volume& v);
    const VolMagick::Volume& colorMappedVolume() const { return _cmVolume; }
    void rgbaVolumes(const std::vector<VolMagick::Volume>& v);
    const std::vector<VolMagick::Volume>& rgbaVolumes() const { return _rgbaVolumes; }

    //experimental accessors to references...
    VolMagick::Volume& colorMappedVolume()
      {
	_cmVolumeUploaded = false;
	return _cmVolume;
      }

    std::vector<VolMagick::Volume>& rgbaVolumes()
      {
	_rgbaVolumeUploaded = false;
	return _rgbaVolumes;
      }

    bool drawBoundingBox() const { return _drawBoundingBox; }
    VolumeRenderingType volumeRenderingType() const { return _volumeRenderingType; }
#ifdef VOLUMEROVER2_FLOAT_COLORTABLE
    const float* colorTable() const { return _colorTable; }
#else
    const unsigned char* colorTable() const { return _colorTable; }
#endif
    bool drawSubVolumeSelector() const { return _drawSubVolumeSelector;}

    bool drawable() const { return _drawable; }
    bool usingVBO() const { return _usingVBO; }

    VolMagick::BoundingBox cmSubVolumeSelector() const { return _cmSubVolumeSelector; }
    VolMagick::BoundingBox rgbaSubVolumeSelector() const { return _rgbaSubVolumeSelector; }
    VolMagick::BoundingBox subVolumeSelector() const
      {
	switch(volumeRenderingType())
	  {
	  default:
	  case ColorMapped: return cmSubVolumeSelector();
	  case RGBA: return rgbaSubVolumeSelector();
	  }
      }
    VolMagick::BoundingBox boundingBox() const
      {
	switch(volumeRenderingType())
	  {
	  default:
	  case ColorMapped: return _cmVolume.boundingBox();
	  case RGBA: return !_rgbaVolumes.empty() ? 
	      _rgbaVolumes[0].boundingBox() : VolMagick::BoundingBox();
	  }
      }

    bool drawCornerAxis() const { return _drawCornerAxis; }
    bool drawGeometry() const { return _drawGeometry; }
    bool drawVolumes() const { return _drawVolumes; }

    typedef cvcraw_geometry::geometry_t geometry_t;
    typedef cvcraw_geometry::cvcgeom_t cvcgeom_t;
    struct scene_geometry_t
    {
      enum render_mode_t { POINTS, LINES,
			   TRIANGLES, TRIANGLE_WIREFRAME,
			   TRIANGLE_FILLED_WIRE, QUADS,
			   QUAD_WIREFRAME, QUAD_FILLED_WIRE,		       
                           TETRA, HEXA,
			   TRIANGLES_FLAT, TRIANGLE_FLAT_FILLED_WIRE };

      // old
      //geometry_t geometry; 

      // new, arand 4-12-2011
      cvcgeom_t geometry;

      // Used for storing geometry that is to be rendered,
      // rather than the real geometry.  This is used to support
      // simple tetra and hexa rendering - Joe, 4-19-2011
      cvcgeom_t render_geometry;

      std::string name;
      render_mode_t render_mode;

      bool visible;

      //VBO stuff.. unused if VBO is not supported
      unsigned int vboArrayBufferID;
      unsigned int vboArrayOffsets[3];
      unsigned int vboVertSize;
      unsigned int vboLineElementArrayBufferID, vboLineSize;
      unsigned int vboTriElementArrayBufferID, vboTriSize;
      unsigned int vboQuadElementArrayBufferID, vboQuadSize;

      /*
      scene_geometry_t(const geometry_t& geom = geometry_t(),
		       const std::string& n = std::string(),
		       render_mode_t mode = TRIANGLES)
	: geometry(geom), name(n), render_mode(mode) { reset_vbo_info(); }
      */

      scene_geometry_t(const cvcgeom_t& geom = cvcgeom_t(),
		       const std::string& n = std::string(),
		       render_mode_t mode = TRIANGLES)
	: geometry(geom), name(n), render_mode(mode), visible(true) 
      { 
        reset_vbo_info(); 
      }

      
      scene_geometry_t(const scene_geometry_t& sg)
	: geometry(sg.geometry), name(sg.name), 
          render_mode(sg.render_mode), visible(sg.visible) 
      { 
        reset_vbo_info();
      }

      ~scene_geometry_t()
      {
      }

      scene_geometry_t& operator=(const scene_geometry_t& geom)
      {
	geometry = geom.geometry;
	name = geom.name;
	render_mode = geom.render_mode;
        visible = geom.visible;
	reset_vbo_info();
	return *this;
      }

      void reset_vbo_info()
      {
	vboArrayBufferID = 0;
	std::fill(vboArrayOffsets, vboArrayOffsets + 3, 0);
	vboVertSize = 0;
	vboLineElementArrayBufferID = 0;
	vboTriElementArrayBufferID = 0;
	vboQuadElementArrayBufferID = 0;
	vboLineSize = 0;
	vboTriSize = 0;
	vboQuadSize = 0;
      }
      
      //Use this to easily get the name of this scene geometry's state object.
      std::string stateName(const std::string& childState = std::string()) const { 
        std::string geo_root = "VolumeViewer::scene_geometry_t"+
          CVC::State::SEPARATOR+
          boost::lexical_cast<std::string>(this);
        return 
          !childState.empty() ? 
          geo_root + CVC::State::SEPARATOR + childState :
          geo_root;
      }

      static render_mode_t render_mode_enum(const std::string& str);
    };

    typedef boost::shared_ptr<scene_geometry_t> scene_geometry_ptr;
    typedef std::map<std::string, scene_geometry_ptr> scene_geometry_collection;


    scene_geometry_collection geometries() const;
    void geometries(const scene_geometry_collection& coll);
#if 0
    scene_geometry_collection& geometries() 
      { 
	_vboUpdated = false;
	return _geometries; 
      }
#endif

    void addGeometry(const cvcraw_geometry::cvcgeom_t& geom,
		     const std::string& name = std::string("geometry"),
		     scene_geometry_t::render_mode_t mode = scene_geometry_t::TRIANGLES,
		     bool do_updategl = true);

    bool clipGeometryToBoundingBox() const { return _clipGeometryToBoundingBox; }
    bool drawGeometryNormals() const { return _drawGeometryNormals; }
    
    //Use this function as a shortcut for producing a string to look up properties
    //with for this object on the property map.
    std::string getObjectName(const std::string& property = std::string()) const;

    std::vector<cvcraw_geometry::cvcgeom_t> getGeometriesFromDatamap() const;

    bool shadedRenderingEnabled() const { return _shadedRenderingEnabled; }
    //multi-tile server related
    bool syncCameraWithMultiTileServer() const { return _syncCameraWithMultiTileServer; }
    bool syncTransferFuncWithMultiTileServer() const { return _syncTransferFuncWithMultiTileServer; }
    bool syncShadedRenderWithMultiTileServer() const { return _syncShadedRenderWithMultiTileServer; }
    bool syncRenderModeWithMultiTileServer() const { return _syncRenderModeWithMultiTileServer; }
    bool interactiveSyncWithMultiTileServer() const { return _interactiveSyncWithMultiTileServer; }
    bool multiTileServerInitialized() const { return _multiTileServerIinitialized; }
    std::string multiTileServerVolumeFile() const { return _multiTileServer_volumeFile; }
    void multiTileServerNtiles(int *x, int *y) const { *x = _multiTileServerNtiles[0]; *y = _multiTileServerNtiles[1]; }
    bool needToUpdateMultiTileServer( int id ) const {
       if( _updateMultiTileServer ) return _updateMultiTileServer[ id ];
       else return false;
    }
    bool terminateMultiTileClient( int id ) const {
       if( _terminateMultiTileClient ) return _terminateMultiTileClient[ id ];
       else return false;
    }

    //Use this to easily get the name of this viewer's state object.
    std::string stateName(const std::string& childState = std::string()) const { 
      std::string viewer_root = "VolumeViewer"+
        CVC::State::SEPARATOR+
        boost::lexical_cast<std::string>(this);
      return 
        !childState.empty() ? 
        viewer_root + CVC::State::SEPARATOR + childState :
        viewer_root;
    }

  public slots:
    void drawBoundingBox(bool draw);
    void volumeRenderingType(VolumeRenderingType t);
#ifdef VOLUMEROVER2_FLOAT_COLORTABLE
    void colorTable(const float * table);
#else
    void colorTable(const unsigned char * table);
#endif
    void drawSubVolumeSelector(bool show);
    void quality(double q);
    void nearPlane(double n);
    void drawCornerAxis(bool draw);
    void drawGeometry(bool draw);
    void drawVolumes(bool draw);
    void clipGeometryToBoundingBox(bool clip);
    void drawGeometryNormals(bool draw);
    void setShadedRenderingMode(bool val);
    
    //multi-tile server related
    void syncCameraWithMultiTileServer(bool val) { _syncCameraWithMultiTileServer = val; }
    void syncTransferFuncWithMultiTileServer(bool val) { _syncTransferFuncWithMultiTileServer = val; }
    void syncShadedRenderWithMultiTileServer(bool val) { _syncShadedRenderWithMultiTileServer = val; }
    void syncRenderModeWithMultiTileServer(bool val) { _syncRenderModeWithMultiTileServer = val; }
    void interactiveSyncWithMultiTileServer(bool val) { _interactiveSyncWithMultiTileServer = val; }
    void multiTileServerInitialized(bool val) { _multiTileServerIinitialized = val; }
    void multiTileServerVolumeFile(std::string fname) { _multiTileServer_volumeFile = fname; }
    void multiTileServerNtiles(int x, int y) { _multiTileServerNtiles[0] = x; _multiTileServerNtiles[1] = y; }
    void needToUpdateMultiTileServer( int id, bool flag ) {
       if( _updateMultiTileServer ) {
         if( id < 0 )
           for(int i=0; i<_multiTileServer_nServer; i++) _updateMultiTileServer[ i ] = flag;
         else
           _updateMultiTileServer[ id ] = flag;
       }
    }
    void terminateMultiTileClient( int id, bool flag ) {
       if( _terminateMultiTileClient ) {
         if( id < 0 )
           for(int i = 0; i<_multiTileServer_nServer; i++ ) _terminateMultiTileClient[i] = flag;
         else
           _terminateMultiTileClient[ id ] = flag;
       }
    }

    void reset(); //resets the viewer to it's default state

    //Use this instead of updateGL() if you want to trigger a redraw
    //but do not need it immediately.  This is safe to call from other
    //threads as well.
    void scheduleUpdateGL();

    bool copyCameraOnNormalize() const { return _copyCameraOnNormalize; }
    void copyCameraOnNormalize(bool flag) { _copyCameraOnNormalize = flag; }
    void copyCameraToPropertyMap();
    bool normalizeOnVolume() const { return _normalizeOnVolume; }
    void normalizeOnVolume(bool flag) { _normalizeOnVolume = flag; }
    bool showEntireSceneOnNormalize() const { return _showEntireSceneOnNormalize; }
    void showEntireSceneOnNormalize(bool flag) { _showEntireSceneOnNormalize = flag; }

    bool glUpdated() const { return _glUpdated; }

  protected slots:
    void timeout();
    void normalizeScene(); //sets the flag, will call doNormalizeScene() on timeout().

  signals:
    void subVolumeSelectorChanged(const VolMagick::BoundingBox& subvolbox);
    void postInit(); //emitted after init() was called

  private:
    void defaultConstructor();
    
  protected:
    virtual void init();
    virtual void draw();
    virtual void postDraw();
    virtual void drawWithNames();
    virtual void endSelection(const QPoint&);
    virtual void postSelection(const QPoint& point);
    virtual QString helpString() const;

    //all drawing functions must be called AFTER init() has been called...
    void doDrawBoundingBox();
    void doDrawSubVolumeBoundingBox();
    void doDrawSubVolumeSelector(bool withNames = false); //must be called last because of the depth buffer clear
    void doDrawCornerAxis();
    void doDrawGeometry();
    void uploadColorMappedVolume();
    void uploadRGBAVolume();
    void uploadColorTable();
    void updateVBO();
    void setClipPlanes();
    void disableClipPlanes();

    void initializeMultiTileServer();
    void syncViewInformation();
    void syncColorTable();
    void syncShadedRenderFlag();
    void syncRenderModeFlag();
    void syncCurrentWithMultiTileServer();

    void doNormalizeScene(); //call to re-center scene around the volume... 
                             //only needed when volume changes, or rendering type changes

    void drawHandle(const qglviewer::Vec& from, const qglviewer::Vec& to, 
		    float radius = -1.0, int nbSubdivisions = 12, 
		    int headName = 0, int shaftName = 0);

    void drawHandle(const qglviewer::Vec& from, const qglviewer::Vec& to,
		    float r = 1.0f, float g = 1.0f, float b = 1.0f,
		    int headName = 0, int shaftName = 0,
		    float pointSize = 2.0f, float lineWidth = 2.0f);
    
    enum HandleName
      {
	MaxXHead = 1, MaxXShaft,
	MinXHead, MinXShaft,
	MaxYHead, MaxYShaft,
	MinYHead, MinYShaft,
	MaxZHead, MaxZShaft,
	MinZHead, MinZShaft,
      };

    void copyBoxToDatamapAndEmitChange();

    // Mouse events functions
    virtual void mousePressEvent(QMouseEvent *e);
    virtual void mouseMoveEvent(QMouseEvent *e);
    virtual void mouseReleaseEvent(QMouseEvent *e);
    virtual void wheelEvent(QWheelEvent* e);

    void customEvent(QEvent *event);

    void copyBoundingBoxToDataMap();

    float getWireframeLineWidth();

    void propertiesChanged(const std::string&);
    void handlePropertiesChanged(const std::string&);
    boost::signals2::connection _propertiesConnection;

    void dataChanged(const std::string&);
    void handleDataChanged(const std::string&);
    boost::signals2::connection _dataConnection;

    void threadsChanged(const std::string&);
    void handleThreadsChanged(const std::string&);
    boost::signals2::connection _threadsConnection;

    void stateChanged(const std::string&);
    void handleStateChanged(const std::string&);
    boost::signals2::connection _stateConnection;

    void geomStateChanged(const std::string&);
    void handleGeomStateChanged(const std::string&);
    boost::signals2::connection _geomStateConnection;

    bool _hasBeenInitialized; //this is set to true once init() has been called.
                              //Do not do any GL calls before this is true!!!!!!!
    
    VolumeRenderer _cmRenderer; //the color mapped volume renderer
    VolumeRenderer _rgbaRenderer; //the RGBA renderer
    VolMagick::Volume _cmVolume; //the color mapped volume (each voxel is unsigned char)
    std::vector<VolMagick::Volume> _rgbaVolumes; //the volumes that make up RGBA components
#ifdef VOLUMEROVER2_FLOAT_COLORTABLE
    float _colorTable[256*4]; //the color table
#else
    unsigned char _colorTable[256*4]; //the color table
#endif
    bool _cmVolumeUploaded;
    bool _rgbaVolumeUploaded;

    bool _shadedRenderingEnabled;

    bool _drawBoundingBox;
    VolumeRenderingType _volumeRenderingType;

    bool _drawable;

    bool _drawSubVolumeSelector;
    VolMagick::BoundingBox _cmSubVolumeSelector;
    VolMagick::BoundingBox _rgbaSubVolumeSelector;

    //the following are used to maintain the selection state
    int _selectedObj;
    QPoint _selectedPoint;

    bool _drawCornerAxis;

    bool _drawGeometry;
    scene_geometry_collection _geometries;

    bool _drawVolumes;

    //If this is true, then we have updated the on-card vertex buffer objects
    //with the contents of _geometries since it was last changed.
    bool _usingVBO;
    bool _vboUpdated; //set to false any time the _geometries map is changed
    //id's of allocated buffers that we should clean up when necessary
    std::vector<unsigned int> _allocatedBuffers;

    typedef std::map<std::string, std::vector<unsigned int> > AllocatedBuffersPerGeom;
    AllocatedBuffersPerGeom _allocatedBuffersPerGeom;
    std::vector<std::string> _allocatedBuffersToReInit;

    //if this is true, geometry is clipped by the bounding box of the current volume
    bool _clipGeometryToBoundingBox;

    bool _drawGeometryNormals;

    QTimer *_timer; //used to track object name changes and possibly other stuff
    QString _oldObjectName;

    //Used when scheduling updateGL calls for later via a custom event.
    //This is useful to compress many redundant updateGL calls into 1 call.
    bool _glUpdated;

    //for parallel multiTileServer
    bool _syncCameraWithMultiTileServer;
    bool _syncTransferFuncWithMultiTileServer;
    bool _syncShadedRenderWithMultiTileServer;
    bool _syncRenderModeWithMultiTileServer;
    bool _interactiveSyncWithMultiTileServer;
    bool _multiTileServerIinitialized;
    bool *_updateMultiTileServer;
    bool *_terminateMultiTileClient;

    std::string _multiTileServer_volumeFile;
    std::string *_multiTileServer_hostList;
    int _multiTileServer_nServer;
    int *_multiTileServer_portList;
    int _multiTileServerNtiles[2];

#ifdef USE_XmlRpc
    std::vector<XmlRpcThread*> _xmlrpcClients;
#else
    CVCSocketTCP *m_socket;
#endif

    bool _mousePressed;

    bool _doDisplayMessage;
    QString _displayMessage;

    bool _mouseMoved;

    bool _copyCameraOnNormalize;
    bool _normalizeOnVolume;
    bool _showEntireSceneOnNormalize;
    bool _doNormalizeScene;

  };

  typedef VolumeViewer::scene_geometry_ptr scene_geometry_ptr;
  typedef VolumeViewer::scene_geometry_t   scene_geometry_t;

#ifdef USE_XmlRpc
  class XmlRpcThread {
  public:
     XmlRpcThread(){ m_client = NULL; m_savedShadedRenderFlag = false; m_savedRenderModeFlag = 0; }
     ~XmlRpcThread() { delete m_client; }

     static void* clientThreadEntryPoint( void * pthis );
     void run(void);
     int start( VolumeViewer *_ptr, std::string _host, int _port, int _id );
     void end(void);

  private:
     VolumeViewer *m_VolumeViewerPtr;
     std::string m_hostname;
     int m_port;
     int m_id;
     XmlRpc::XmlRpcClient *m_client;
     pthread_t m_thread_id;
     bool m_savedShadedRenderFlag;
     int m_savedRenderModeFlag;
  };
#endif
}

#endif
