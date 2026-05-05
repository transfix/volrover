#include <cvc_compat.h>
#include <cvc_compat.h>
/*
  Copyright 2012 The University of Texas at Austin

        Authors: Jose Rivera <transfix@ices.utexas.edu>
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

#include <VolumeRover2/VolumeViewer2.h>
#include <CVCEvent.h>
#include <cvc_compat.h>
#include <VolMagick/Volume.h>
#include <cvcraw_geometry/cvcgeom.h>
#include <QGLViewer/frame.h>
#include <QApplication>
#include <boost/shared_array.hpp>
#include <boost/bind.hpp>
#include <boost/format.hpp>
#include <cmath>

#include <QPoint>
#include <QMouseEvent>

namespace
{
  static inline unsigned int upToPowerOfTwoLocal(unsigned int value)
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
}

namespace CVC_NAMESPACE
{
  VolumeViewer2::~VolumeViewer2()
  {
    _threadsConnection.disconnect();
    _stateConnection.disconnect();
  }

  bool VolumeViewer2::glUpdated() const
  {
    return bool(cvcstate("glUpdated").value<int>());
  }

  bool VolumeViewer2::hasBeenInitialized() const
  {
    return bool(cvcstate("hasBeenInitialized").value<int>());
  }

  void VolumeViewer2::defaultConstructor()
  {
    //NO GL CALLS ALLOWED HERE!!!! do them in init()
    //I mean it, you will have a nasty crash in the GL driver otherwise! This is because
    //this function is called without the guarantee that the GL context yet exists.

    //State objects meant to be changed by users are between =====

    _stateConnection.disconnect();

    //===============================
    //set this object's initial state
    cvcstate("rendering_mode")
      .value("colormapped")
      .comment("The volume viewer's volume rendering mode: 'colormapped' or 'rgba'");
    cvcstate("shaded_rendering_enabled")
      .value(int(false))
      .comment("Toggle shaded volume rendering if hardware support available.");
    cvcstate("draw_bounding_box")
      .value(int(true))
      .comment("Toggle the overall bounding box");
    cvcstate("draw_subvolume_selector")
      .value(int(true))
      .comment("Toggle the subvolume selector visibility for specifiying subvolumes.");
    cvcstate("volume_rendering_quality")
      .value(double(0.5))
      .comment("A double value from [0.0,1.0]. Lower number, the faster the rendering but lower the quality.");
    cvcstate("volume_rendering_near_plane")
      .value("0.0")
      .comment("A double value from [0.0,1.0]. If 0.0, the volume rendering near plane is at the view plane. "
               "1.0 means the entire volume is clipped out.");
    cvcstate("projection_mode")
      .value("perspective")
      .comment("The OpenGL projection mode to use: 'perspective' or 'orthographic.'");
    cvcstate("draw_corner_axis")
      .value(int(true))
      .comment("Toggle the corner axis.  Used as a visual hint of orientation in the world coordinate system.");
    cvcstate("draw_geometry")
      .value(int(true))
      .comment("Toggle drawing of geometry.");
    cvcstate("draw_volumes")
      .value(int(true))
      .comment("Toggle drawing of volumes.");
    cvcstate("clip_geometry")
      .value(int(true))
      .comment("Toggle whether the rendered geometry should be clipped by the overall volume bounding box.");
    cvcstate("background_color")
      .value("#000000")
      .comment("Set the background color.  Uses HTML color format.");
    cvcstate("fov")
      .value(boost::lexical_cast<std::string>(M_PI/4.0f))
      .comment("The volume viewer's field of view.");
    cvcstate("show_entire_scene_on_normalize")
      .value(int(true))
      .comment("Toggle whether to reset the camera such that the entire scene is visible whenever the scene"
               " is re-normalized because something was added or removed.");

    //stereo related properties
    cvcstate("io_distance")
      .value("0.062")
      .comment("Stereo: intra-ocular distance in meters");
    cvcstate("physical_distance_to_screen")
      .value("2.0")
      .comment("Stereo: physical distance to screen in meters");
    cvcstate("physical_screen_width")
      .value("1.8")
      .comment("Stereo: physical screen width in meters");
    cvcstate("focus_distance")
      .value("1000.0")
      .comment("Stereo: focus distance in meters");

    //===============================

#if 0 //the following should be handled by scene geometry objects
    cvcstate("draw_geometry_normals")
      .value(int(false));
    cvcstate("geometry_line_width").value("1.2");
#endif

    cvcstate("syncCamera_with_multiTileServer").value(int(false));
    cvcstate("syncTransferFunc_with_multiTileServer").value(int(false));
    cvcstate("syncShadedRender_with_multiTileServer").value(int(false));
    cvcstate("syncRenderMode_with_multiTileServer").value(int(false));
    cvcstate("interactiveMode_with_multiTileServer").value(int(false));
    cvcstate("syncMode_with_multiTileServer").value("0");
    cvcstate("syncInitial_multiTileServer").value("0");

    cvcstate("glUpdated")
      .value(int(false))
      .comment("Used when scheduling updateGL calls for later via a custom event. "
               "This is useful to compress many redundant updateGL calls into 1 call. "
               "When someone sets this to 0, an openGL redraw has been scheduled. "
               "When it is set back to non-zero, the redraw has finished.")
      .hidden(true);

    //this is set to true once init() has been called.
    //Do not do any GL calls on this widget's GL context before this is true!!!!!!!
    cvcstate("hasBeenInitialized")
      .value(int(false))
      .comment("This is set to non-zero once init() has been called. "
               "Do not do any GL calls on this widget's GL context before this is true!!!!!!!")
      .hidden(true);

    cvcstate("cmVolumeUploaded")
      .value(int(false))
      .comment("Flag set when the graphics card has the latest colormapped volume.")
      .hidden(true);

    cvcstate("rgbaVolumeUploaded")
      .value(int(false))
      .comment("Flag set when the graphics card has the latest rgba volume.")
      .hidden(true);

    cvcstate("drawable")
      .value(int(false))
      .comment("Flag set when initializing the colormapped or rgba volume renderers without error.")
      .hidden(true);

    cvcstate("usingVBO")
      .value(int(false))
      .comment("Flag set when this viewer is using vertex buffer objects for fast geometry rendering.")
      .hidden(true);
    cvcstate("vboUpdated")
      .value(int(false))
      .comment("Flag set when the VBOs for geometry have been updated with the latest data from system memory.")
      .hidden(true);

    cvcstate("colortableUploaded")
      .value(int(false))
      .comment("Flag set when the color table has been uploaded to the graphics card.")
      .hidden(true);

    cvcstate("sub_volume_selector")
      .value(cvcstate("scene.region_of_interest").fullName())
      .comment("This is set to the name of the scene object containing a bounding box to use as the sub volume region of interest selector.")
      .hidden(true);

    cvcstate("selected_object")
      .value(int(-1))
      .comment("The id of the selected object, -1 if none")
      .hidden(true);

    cvcstate("selected_point")
      .value("QPoint")
      .data(QPoint())
      .comment("The screen point where the last selected object was selected from.")
      .hidden(true);

    cvcstate("mouse_pressed")
      .value(int(false))
      .comment("Flag set when the mouse button is pressed for dragging.")
      .hidden(true);

    cvcstate("mouse_moved")
      .value(int(false))
      .comment("Flag set when the mouse has been clicked and dragged.")
      .hidden(true);

    //watch the thread map to print messages for the user
    _threadsConnection.disconnect();
    _threadsConnection = 
      cvcapp.threadsChanged.connect(
          MapChangeSignal::slot_type(
            &VolumeViewer2::threadsChanged, this, _1
          )
        );
  }

  void VolumeViewer2::customEvent(QEvent *event)
  {
    CVCEvent *mwe = dynamic_cast<CVCEvent*>(event);
    if(!mwe) return;

    try
      {
        if(mwe->name == "handleThreadsChanged")
          handleThreadsChanged(boost::any_cast<std::string>(mwe->data));
        else if(mwe->name == "updateGL" && !cvcstate("glUpdated").value<int>())
          {
            update();
            cvcstate("glUpdated").value(int(true));
          }
      }
    catch(Exception& e)
      {
        using namespace boost; using boost::format;
        cvcapp.log(2,str(boost::format("%s :: Exception: %s\n")
                         % BOOST_CURRENT_FUNCTION
                         % e.what()));
      }
  }

  void VolumeViewer2::mousePressEvent(QMouseEvent *e)
  {
    try
      {
        select(e->pos());
        if(cvcstate("selected_object").value<int>() != -1) {
          cvcstate("selected_point").data(e->pos());
        }
        else {
          QGLViewer::mousePressEvent(e);
        }
        // to test handling XmlRpc delay
        cvcstate("mouse_pressed").value(int(true));
      }
    catch(Exception& e)
      {
        using namespace boost; using boost::format;
        cvcapp.log(2,str(boost::format("%s :: Exception: %s\n")
                         % BOOST_CURRENT_FUNCTION
                         % e.what()));
      }
  }

  void VolumeViewer2::mouseMoveEvent(QMouseEvent *e)
  {
    using namespace qglviewer;

    try
      {
        int selectedObj = cvcstate("selected_object").value<int>();

        if(selectedObj != -1)
          {
            double minx, miny, minz;
            double maxx, maxy, maxz;
	
            BoundingBox bbox = cvcstate(cvcstate("sub_volume_selector").value()).data<BoundingBox>();
            minx = bbox.minx; miny = bbox.miny; minz = bbox.minz;
            maxx = bbox.maxx; maxy = bbox.maxy; maxz = bbox.maxz;

            Vec center((maxx-minx)/2+minx,
                       (maxy-miny)/2+miny,
                       (maxz-minz)/2+minz);
	
            QPoint selectedPoint = cvcstate("selected_point").data<QPoint>();
            Vec screenCenter = camera()->projectedCoordinatesOf(center);
            Vec transPoint(e->pos().x(),e->pos().y(),screenCenter.z);
            Vec selectPoint(selectedPoint.x(),selectedPoint.y(),screenCenter.z);

            switch(selectedObj)
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

                  if(selectedObj == MaxXShaft)
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

                  if(selectedObj == MinXShaft)
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

                  if(selectedObj == MaxYShaft)
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

                  if(selectedObj == MinYShaft)
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

                  if(selectedObj == MaxZShaft)
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

                  if(selectedObj == MinZShaft)
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

            BoundingBox global_bbox = globalBoundingBox();

            bbox.minx = std::min(std::max(minx,global_bbox.XMin()),global_bbox.XMax());
            bbox.miny = std::min(std::max(miny,global_bbox.YMin()),global_bbox.YMax());
            bbox.minz = std::min(std::max(minz,global_bbox.ZMin()),global_bbox.ZMax());
            bbox.maxx = std::max(std::min(maxx,global_bbox.XMax()),global_bbox.XMin());
            bbox.maxy = std::max(std::min(maxy,global_bbox.YMax()),global_bbox.YMin());
            bbox.maxz = std::max(std::min(maxz,global_bbox.ZMax()),global_bbox.ZMin());

            cvcstate(cvcstate("sub_volume_selector").value()).data(bbox);
            cvcstate("selected_point").data(e->pos());

            if(hasBeenInitialized()) scheduleUpdateGL();
          }
        else
          {
            QGLViewer::mouseMoveEvent(e);

            if(cvcstate("mouse_pressed").value<int>())
              cvcstate("mouse_moved").value(int(true));
          }
      }
    catch(Exception& e)
      {
        using namespace boost; using boost::format;
        cvcapp.log(2,str(boost::format("%s :: Exception: %s\n")
                         % BOOST_CURRENT_FUNCTION
                         % e.what()));
      }
  }

  void VolumeViewer2::mouseReleaseEvent(QMouseEvent *e)
  {
    try
      {
        if(cvcstate("selected_object").value<int>())
          {
            cvcstate("selected_object").value(int(-1));
	
            //lets normalize the bounding box in case it is inverted (i.e. min > max)
            BoundingBox bbox = cvcstate(cvcstate("sub_volume_selector").value()).data<BoundingBox>();
            bbox.normalize();
            cvcstate(cvcstate("sub_volume_selector").value()).data(bbox);
            
            cvcapp.log(
                       6,
                       boost::str(
                                  boost::format("%s: (%f,%f,%f) (%f,%f,%f)\n")
                                  % BOOST_CURRENT_FUNCTION
                                  % bbox.minx 
                                  % bbox.miny 
                                  % bbox.minz
                                  % bbox.maxx 
                                  % bbox.maxy 
                                  % bbox.maxz
                                  )
                       );
          }

        QGLViewer::mouseReleaseEvent(e);
        // to test handling XmlRpc delay
        cvcstate("mouse_pressed").value(int(false));

        // for quaternion updates
        cvcstate("mouse_moved").value(int(false));
      }
    catch(Exception& e)
      {
        using namespace boost; using boost::format;
        cvcapp.log(2,str(boost::format("%s :: Exception: %s\n")
                         % BOOST_CURRENT_FUNCTION
                         % e.what()));
      }
  }

  void VolumeViewer2::wheelEvent(QWheelEvent* e)
  {
    try
      {
        if(cvcstate("mouse_pressed").value<int>())
          cvcstate("mouse_moved").value(int(true));
        QGLViewer::wheelEvent(e);
      }
    catch(Exception& e)
      {
        using namespace boost; using boost::format;
        cvcapp.log(2,str(boost::format("%s :: Exception: %s\n")
                         % BOOST_CURRENT_FUNCTION
                         % e.what()));
      }
  }

  void VolumeViewer2::setDefaultColorTable()
  {
    boost::shared_array<unsigned char> table(new unsigned char[256*4]);

    //set the color table to greyscale for now
    for(int i=0; i<256; i++)
      {
	table[i*4+0] = i;
	table[i*4+1] = i;
	table[i*4+2] = i;
	table[i*4+3] = i;
      }

    cvcstate("colortable").data(table);
  }

  void VolumeViewer2::setDefaultScene()
  {
    cvcstate("scene").reset();
    cvcstate("scene.region_of_interest").data(BoundingBox(-0.25,-0.25,-0.25,0.25,0.25,0.25));

    cvcstate("scene.volume").data(VolMagick::Volume(VolMagick::Dimension(4,4,4),
                                                 UChar,
                                                 VolMagick::BoundingBox(-0.5,-0.5,-0.5,0.5,0.5,0.5)));
    cvcstate("scene.volume").value("colormapped"); //effects the way this volume is rendered. can be "colormapped" or "rgba"
    cvcstate("scene.volume.uploaded")
      .value(int(false))
      .comment("Flag set when the graphics card has the latest volume");
  }

  void VolumeViewer2::normalizeScene(bool forceShowAll)
  {
    BoundingBox global_bbox = globalBoundingBox();
    setSceneCenter(qglviewer::Vec((global_bbox.XMax()-global_bbox.XMin())/2.0+global_bbox.XMin(),
                                  (global_bbox.YMax()-global_bbox.YMin())/2.0+global_bbox.YMin(),
                                  (global_bbox.ZMax()-global_bbox.ZMin())/2.0+global_bbox.ZMin()));
    setSceneRadius(std::max(global_bbox.XMax()-global_bbox.XMin(),
                            std::max(global_bbox.YMax()-global_bbox.YMin(),
                                     global_bbox.ZMax()-global_bbox.ZMin())));
    if(forceShowAll ||
       cvcstate("show_entire_scene_on_normalize").value<int>())
      showEntireScene();
  }

  BoundingBox VolumeViewer2::globalBoundingBox() const
  {
    BoundingBox global_bbox;
    if(!cvcstate("global_bounding_box").isData<BoundingBox>())
      calculateGlobalBoundingBox();
    return cvcstate("global_bounding_box").data<BoundingBox>();
  }

  void VolumeViewer2::reset()
  {
    if(hasBeenInitialized()) makeCurrent();
    defaultConstructor();
    init();
    if(hasBeenInitialized()) scheduleUpdateGL();
  }

  void VolumeViewer2::scheduleUpdateGL()
  {
    cvcstate("glUpdated").value(int(false));
    if(hasBeenInitialized())
      QCoreApplication::postEvent(this,new CVCEvent("updateGL"));
  }

  void VolumeViewer2::resetROI()
  {
    BoundingBox global_bbox = globalBoundingBox();

    double minx, miny, minz;
    double maxx, maxy, maxz;
    minx = global_bbox.XMin();
    miny = global_bbox.YMin();
    minz = global_bbox.ZMin();
    maxx = global_bbox.XMax();
    maxy = global_bbox.YMax();
    maxz = global_bbox.ZMax();

    double centerx, centery, centerz;
    centerx = (maxx-minx)/2.0+minx;
    centery = (maxy-miny)/2.0+miny;
    centerz = (maxz-minz)/2.0+minz;

    BoundingBox default_roi(centerx - (centerx - minx)/2.0,
                            centery - (centery - miny)/2.0,
                            centerz - (centerz - minz)/2.0,
                            centerx + (maxx - centerx)/2.0,
                            centery + (maxy - centery)/2.0,
                            centerz + (maxz - centerz)/2.0);

    cvcstate("scene.region_of_interest").data(default_roi);
  }

  void VolumeViewer2::init()
  {
    setSelectRegionWidth(10);
    setSelectRegionHeight(10);
    setTextIsEnabled(true);

    setMouseTracking(true);

    if(!cvcstate("drawable").value<int>() &&
       _renderer.initRenderer())
      cvcstate("drawable").value(int(true));

    setDefaultColorTable();
    setDefaultScene();

    cvcstate("hasBeenInitialized").value(int(true));
    emit postInit();
  }

  void VolumeViewer2::draw()
  {
    copyCameraToState();

    //main scene rendering via renderStates
    cvcstate("scene").traverse(
      boost::bind(&VolumeViewer2::renderStates, 
        boost::ref(*this), 
        _1));
  }

  void VolumeViewer2::drawWithNames()
  {
    //main scene rendering via renderStates
    cvcstate("scene").traverse(
      boost::bind(&VolumeViewer2::renderStatesWithNames, 
        boost::ref(*this), 
        _1));
  }

  void VolumeViewer2::postDraw()
  {
    QGLViewer::postDraw();
    if(cvcstate("draw_corner_axis").value<int>()) 
      doDrawCornerAxis();
  }

  void VolumeViewer2::doDrawCornerAxis()
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
  }

  void VolumeViewer2::renderStates(const std::string& s)
  {
    cvcapp.log(4,
               str(
                   boost::format("%s :: state: %s\n")
                   % BOOST_CURRENT_FUNCTION
                   % s));

    if(cvcstate("draw_volumes").value<int>() &&             //single volume
       cvcstate(s).isData<VolMagick::Volume>())
      {
        cvcapp.log(4,
                   str(
                       boost::format("%s :: rendering volume %s\n")
                       % BOOST_CURRENT_FUNCTION
                       % s));

        if(cvcstate("draw_bounding_box").value<int>())
          doDrawBoundingBox(s);

        //TODO: volume renderer generates geometry
      }
    if(cvcstate("draw_volumes").value<int>() &&             //vector of volumes
       cvcstate(s).isData<std::vector<VolMagick::Volume> >())
      {
        cvcapp.log(4,
                   str(
                       boost::format("%s :: rendering volumes %s\n")
                       % BOOST_CURRENT_FUNCTION
                       % s));

        //get overall bounding box
        std::vector<VolMagick::Volume> vols = 
          cvcstate(s).data<std::vector<VolMagick::Volume> >();
        BoundingBox bbox = vols.empty() ? BoundingBox() : vols[0].boundingBox();
        for (auto& vol : vols)
          bbox += vol.boundingBox();

        if(cvcstate("draw_bounding_box").value<int>())
          doDrawBoundingBox(bbox);

        //TODO: volume renderer generates geometry
      }
    else if(cvcstate("draw_geometry").value<int>() &&
            cvcstate(s).isData<cvcraw_geometry::cvcgeom_t>())
      {
        cvcapp.log(4,
                   str(
                       boost::format("%s :: rendering geometry %s\n")
                       % BOOST_CURRENT_FUNCTION
                       % s));

        //TODO: replace the render pass with a call to this and a call to the volume
        //renderer to generate geometry to render
        _geometryRenderer.render(s); 
      }
    else if(cvcstate(s).isData<BoundingBox>())
      {
        cvcapp.log(4,
                   str(
                       boost::format("%s :: bounding box %s\n")
                       % BOOST_CURRENT_FUNCTION
                       % s));

        doDrawBoundingBox(s);

        //draw subvolume selector if flag set
        if(cvcstate("draw_subvolume_selector").value<int>() ||
           cvcstate(s)("draw_subvolume_selector").value<int>())
          doDrawSubVolumeSelector(s);
      }
    else if(cvcstate(s).value() == "glPushMatrix")
      {
        cvcapp.log(4,
                   str(
                       boost::format("%s :: glPushMatrix %s\n")
                       % BOOST_CURRENT_FUNCTION
                       % s));

        glPushMatrix();
      }
    else if(cvcstate(s).value() == "glPopMatrix")
      {
        cvcapp.log(4,
                   str(
                       boost::format("%s :: glPopMatrix %s\n")
                       % BOOST_CURRENT_FUNCTION
                       % s));

        glPopMatrix();
      }
    else if(cvcstate(s).isData<qglviewer::Frame>())
      {
        cvcapp.log(4,
                   str(
                       boost::format("%s :: Frame %s\n")
                       % BOOST_CURRENT_FUNCTION
                       % s));

        qglviewer::Frame fr = cvcstate(s).data<qglviewer::Frame>();
        glMultMatrixd(fr.matrix());
      }
  }

  void VolumeViewer2::renderStatesWithNames(const std::string& s)
  {
    cvcapp.log(4,
               str(
                   boost::format("%s :: state: %s\n")
                   % BOOST_CURRENT_FUNCTION
                   % s));

    if(cvcstate(s).isData<BoundingBox>())
      {
        cvcapp.log(4,
                   str(
                       boost::format("%s :: bounding box %s\n")
                       % BOOST_CURRENT_FUNCTION
                       % s));
        
        //draw subvolume selector if flag set
        if(cvcstate("draw_subvolume_selector").value<int>() ||
           cvcstate(s)("draw_subvolume_selector").value<int>())
          doDrawSubVolumeSelector(s,true);
      }
  }

  void VolumeViewer2::calculateGlobalBoundingBox() const
  {
    BoundingBox global_bbox;
    cvcstate("scene").traverse(
      boost::bind(&VolumeViewer2::combineBoundingBox,
                  boost::ref(global_bbox),
                  _1));
    if(global_bbox.isNull()) //disallow null bboxes
      global_bbox = BoundingBox(-0.5,-0.5,-0.5,0.5,0.5,0.5);
    cvcstate("global_bounding_box").data(global_bbox);
  }

  void VolumeViewer2::combineBoundingBox(BoundingBox& bbox, const std::string& s)
  {
    if(cvcstate(s).isData<VolMagick::Volume>())
      {
        VolMagick::Volume vol = cvcstate(s).data<VolMagick::Volume>();
        bbox += vol.boundingBox();
      }
    if(cvcstate(s).isData<std::vector<VolMagick::Volume> >())
      {
        //get overall bounding box
        std::vector<VolMagick::Volume> vols = 
          cvcstate(s).data<std::vector<VolMagick::Volume> >();
        BoundingBox all_bbox = vols.empty() ? BoundingBox() : vols[0].boundingBox();
        for (auto& vol : vols)
          all_bbox += vol.boundingBox();

        bbox += all_bbox;
      }
    else if(cvcstate(s).isData<cvcraw_geometry::cvcgeom_t>())
      {
        //TODO: get bounding box from geometry!
      }
    else if(cvcstate(s).isData<BoundingBox>())
      {
        bbox += cvcstate(s).data<BoundingBox>();
      }
  }

  void VolumeViewer2::endSelection(const QPoint&)
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

  void VolumeViewer2::postSelection(const QPoint& point)
  {
    cvcstate("selected_point").data(point);
    cvcstate("selected_object").value(int(selectedName()));
  }

  void VolumeViewer2::setClipPlanes()
  {
    BoundingBox bbox = globalBoundingBox();
    
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

  void VolumeViewer2::threadsChanged(const std::string& key)
  {
    QCoreApplication::postEvent(this,
      new CVCEvent("handleThreadsChanged",key)
    );
  }

  void VolumeViewer2::handleThreadsChanged(const std::string& key)
  {
    using namespace std;
    using namespace boost; using boost::format;

    cvcapp.log(5,
               str(
                   boost::format("%s :: state: %s\n")
                   % BOOST_CURRENT_FUNCTION
                   % key));
  }

  void VolumeViewer2::doDrawBoundingBox(const std::string& s)
  {
    BoundingBox bbox;

    if(cvcstate(s).isData<VolMagick::Volume>())
      bbox = cvcstate(s).data<VolMagick::Volume>().boundingBox();
    else
      bbox = cvcstate(s).data<BoundingBox>();

    doDrawBoundingBox(bbox);
  }

  void VolumeViewer2::doDrawBoundingBox(const BoundingBox& bbox)
  {
    double minx, miny, minz;
    double maxx, maxy, maxz;
    float bgcolor[4];

    minx = bbox.XMin();
    miny = bbox.YMin();
    minz = bbox.ZMin();
    maxx = bbox.XMax();
    maxy = bbox.YMax();
    maxz = bbox.ZMax();

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

  void VolumeViewer2::doDrawSubVolumeSelector(const std::string& s, bool withNames)
  {
    using namespace qglviewer;

    double minx, miny, minz;
    double maxx, maxy, maxz;
    
    BoundingBox bbox;

    if(cvcstate(s).isData<VolMagick::Volume>())
      bbox = cvcstate(s).data<VolMagick::Volume>().boundingBox();
    else
      bbox = cvcstate(s).data<BoundingBox>();

    minx = bbox.XMin();
    miny = bbox.YMin();
    minz = bbox.ZMin();
    maxx = bbox.XMax();
    maxy = bbox.YMax();
    maxz = bbox.ZMax();

    double centerx, centery, centerz;
    centerx = (maxx-minx)/2+minx;
    centery = (maxy-miny)/2+miny;
    centerz = (maxz-minz)/2+minz;

    /*** draw the selector handle arrows ***/
    float pointSize = withNames ? 6.0 : 6.0;
    float lineWidth = withNames ? 4.0 : 2.0;

    glClear(GL_DEPTH_BUFFER_BIT); //show arrows on top of everything

    //TODO: make the names per-object so we can have multiple selectable objects
    //in the scene...
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

  void VolumeViewer2::drawHandle(const qglviewer::Vec& from, const qglviewer::Vec& to, 
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

  void VolumeViewer2::drawHandle(const qglviewer::Vec& from, const qglviewer::Vec& to, 
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

  void VolumeViewer2::handleStateChanged(const std::string& childState)
  {
    using namespace std;
    using namespace boost; using boost::format;

    boost::mutex::scoped_lock lock(_handleStateChangedMutex);

    cvcapp.log(5,
               str(
                   boost::format("%s :: state: %s\n")
                   % BOOST_CURRENT_FUNCTION
                   % childState));

    bool updateGLNeeded = true;
    if(childState == "colortable")
      {
        cvcstate("colortableUploaded").value(int(false));
      }
    else if(childState == "scene")
      {
        cvcstate("vboUpdated").value(int(false));
        calculateGlobalBoundingBox();
      }
    else if(childState == "shaded_rendering_enabled")
      {
        if(hasBeenInitialized()) makeCurrent();
        if(cvcstate(childState).value<int>())
          _renderer.enableShadedRendering();
        else
          _renderer.disableShadedRendering();
      }
    else if(childState == "global_bounding_box")
      {
        normalizeScene();
      }
    else if(childState == "volume_rendering_quality")
      {
        if(hasBeenInitialized()) makeCurrent();
        _renderer.setQuality(cvcstate(childState).value<double>());
      }
    else if(childState == "volume_rendering_near_plane")
      {
        if(hasBeenInitialized()) makeCurrent();
        _renderer.setNearPlane(cvcstate(childState).value<double>());
      }
    else if(childState == "projection_mode")
      {
        if(cvcstate("projection_mode").value()=="perspective")
          camera()->setType(qglviewer::Camera::PERSPECTIVE);
        else if(cvcstate("projection_mode").value()=="orthographic")
          camera()->setType(qglviewer::Camera::ORTHOGRAPHIC);
        else
          cvcapp.log(1,
                     str(
                         boost::format("%s :: Unknown projection_mode\n")
                         % BOOST_CURRENT_FUNCTION
                         % childState));
      }
    else if(childState == "background_color")
      {
        setBackgroundColor(QColor(cvcstate(childState).value().c_str()));
      }
    else if(childState == "fov")
      {
        camera()->setFieldOfView(cvcstate(childState).value<double>());
      }
    else if(childState == "io_distance")
      {
        camera()->setIODistance(cvcstate(childState).value<double>());
      }
    else if(childState == "physical_distance_to_screen")
      {
        camera()->setPhysicalDistanceToScreen(cvcstate(childState).value<double>());
      }
    else if(childState == "physical_screen_width")
      {
        camera()->setPhysicalScreenWidth(cvcstate(childState).value<double>());
      }
    else if(childState == "focus_distance")
      {
        camera()->setFocusDistance(cvcstate(childState).value<double>());
      }

    if(updateGLNeeded) scheduleUpdateGL();
  }

  void VolumeViewer2::uploadColorTable()
  {
    if(hasBeenInitialized()) makeCurrent();
    if(cvcstate("colortableUploaded").value<int>()) return;
    boost::shared_array<unsigned char> table =
      cvcstate("colortable").data<boost::shared_array<unsigned char> >();
    _renderer.uploadColorMap(table.get());
    cvcstate("colortableUploaded").value(int(true));
  }

  void VolumeViewer2::uploadVolume(const std::string& s)
  {

  }

  void VolumeViewer2::updateVBO()
  {

  }

  void VolumeViewer2::copyCameraToState()
  {
    //set orientation property
    {
      std::stringstream ss;
      qglviewer::Quaternion Q = camera()->orientation();
      ss<< Q[0] << '\t' << Q[1] << '\t' << Q[2] << '\t' << Q[3];
      //ss << camera()->orientation();

      if(cvcstate("orientation").value()!=ss.str())
        cvcstate("orientation").value(ss.str());
    }
    
    //set camera position property
    {
      using namespace boost; using boost::format;
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

      if(cvcstate("position").value()!=posstr)
        cvcstate("position").value(posstr);
    }
  }
}
