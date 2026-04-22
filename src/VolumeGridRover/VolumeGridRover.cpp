/*
  Copyright 2005-2008 The University of Texas at Austin

        Authors: Joe Rivera <transfix@ices.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of VolumeGridRover.

  VolumeGridRover is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  VolumeGridRover is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

/* $Id: VolumeGridRover.cpp 5276 2012-03-15 18:06:20Z deukhyun $ */

// MUST include GLEW before any OpenGL headers
#include <GL/glew.h>

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <fstream>
#include <algorithm>
#include <utility>
#include <set>
#include <sstream>
#include <stdexcept>
//#include <boost/lambda/lambda.hpp>
//#include <boost/lambda/construct.hpp>
#include <iterator>
#include <boost/range.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/utility.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>

//used in solving for control points in b-spline fitting
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/scoped_array.hpp>
#include <boost/regex.hpp>

#include <boost/lexical_cast.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
//#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/exception.hpp>

#include <boost/format.hpp>

#include <boost/current_function.hpp>

#include <QApplication>
#include <QMouseEvent>
#include <QStyleFactory>

#include <qpainter.h>
#include <qcolor.h>
#include <QString>
#include <qstringlist.h>
#include <qlayout.h>
#include <qimage.h>
#include <qcolordialog.h>
#include <qpushbutton.h>
#include <qcombobox.h>
#include <qspinbox.h>
#include <qcheckbox.h>
#include <qvalidator.h>
#include <QFileDialog>
#include <QFileInfo>
#include <QMessageBox>
#include <qprogressdialog.h>
#include <QFile>
#include <QTextStream>
#include <qbuttongroup.h>
#include <qsettings.h>
#include <QDir>
#include <qlistview.h>

#include <QXmlStreamReader>
#include <QXmlStreamWriter>
#include <qdom.h>

#include <xmlrpc/XmlRpc.h>
#include <Segmentation/GenSeg/genseg.h>

#include <QGLViewer/manipulatedCameraFrame.h>
#include <QGLViewer/constraint.h>
#include <QSurfaceFormat>

#include <VolumeGridRover/VolumeGridRover.h>
#include <VolumeGridRover/bspline_opt.h>
#include <VolumeGridRover/sdf_opt.h>

#include <CVC/App.h>

#ifdef USING_VOLUMEGRIDROVER_MEDAX
#include <VolumeGridRover/medax.h>
#endif

#include <GL/glew.h>

// using GSL library for 2D contour interpolation
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

#ifdef USING_TILING
#include <ContourTiler/Contour.h>
#include <ContourTiler/reader_ser.h>
// #include <Tiling/tiling.h>
// #include <cvcraw_geometry/Geometry.h>
// #include <GeometryFileTypes/GeometryLoader.h>
// #include <Tiling/SeriesFileReader.h>
#endif

#include <CGAL/Simple_cartesian.h> //for calculating planes for taoju's contour data format

//for outputting contour edge normals for Lei Na's b-spline ==> a-spline conversion
#include <CGAL/Point_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Vector_2.h>

#ifdef VOLUMEGRIDROVER_ISOCONTOURING
#include <Contour/contour.h>
#include <Contour/datasetreg2.h>
#endif

#include <ui_VolumeGridRoverBase.h>

using namespace XmlRpc;
using namespace qglviewer;

#define TEXTURE_TILE_X  256
#define TEXTURE_TILE_Y  256

//Define the following to filter out contour points that greatly change the
//calculated bounding box when generating a new volume upon contour load.
//#define CONTOUR_BOUNDING_BOX_FILTER_HEURISTIC

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

static inline float clamp(float val)
{
  return 
    val > 1.0 ? 1.0 :
    val < 0.0 ? 0.0 :
    val;
}

PointClass::PointClass(const QColor &color, const QString &name)
  : m_PointClassColor(color), m_Name(name)
{
  m_PointList.clear();
}

PointClass::~PointClass()
{
  m_PointList.clear();
}

void PointClass::addPoint(unsigned int x, unsigned int y, unsigned int z)
{
  if(pointIsInClass(x,y,z)) return; /* avoid duplicates */
  m_PointList.append(new GridPoint(x,y,z));
}

bool PointClass::removePoint(unsigned int x, unsigned int y, unsigned int z)
{
  int size = m_PointList.size();
  for( int i = 0; i < size; i++ )
  {
    GridPoint *gp = m_PointList[ i ];
    if(gp->x == x && gp->y == y && gp->z == z)
      {
	m_PointList.removeOne( gp );
	return true;
      }
   }

  return false;
}

bool PointClass::pointIsInClass(unsigned int x, unsigned int y, unsigned int z)
{
  int size = m_PointList.size();
  for( int i = 0; i < size; i++ )
  {
    GridPoint *gp = m_PointList[ i ];
    if(gp->x == x && gp->y == y && gp->z == z) return true;
  }

  return false;
}

GridPoint *PointClass::getClosestPoint(unsigned int x, unsigned int y, unsigned int z, Plane onPlane)
{
  int pointListSize = m_PointList.size();

  if( pointListSize < 1 ) return NULL;

  GridPoint *closest = m_PointList[ 0 ];
  double dist = sqrt((double)(closest->x - x)*(closest->x - x)+(closest->y - y)*(closest->y - y)+(closest->z - z)*(closest->z - z));
  for( int i = 1; i < pointListSize; i++ )
  {
      GridPoint *gp = m_PointList[ i ];
      if(gp->x == x && gp->y == y && gp->z == z) return gp;
      double newdist = sqrt((double)(gp->x - x)*(gp->x - x)+(gp->y - y)*(gp->y - y)+(gp->z - z)*(gp->z - z));
      if(newdist > 0.0 && newdist < dist &&
	 (onPlane == XY ? gp->z == z : 
	  onPlane == XZ ? gp->y == y : 
	  onPlane == ZY ? gp->x == x : 
	  onPlane == ANY ? true : false))
	{
	  dist = newdist;
	  closest = gp;
	}
   }
  return closest;
}

SliceCanvas::SliceCanvas(SliceCanvas::SliceAxis a, QWidget * parent, const char * name)
  : QGLViewer(parent), m_Depth(0), m_PointClassList(NULL), 
    m_CurrentPointClass(NULL), m_SliceAxis(a), m_PointSize(1), m_Variable(0), m_Timestep(0), m_Drawable(true),
    m_SliceTiles(NULL), m_NumSliceTiles(0), m_MouseIsDown(false), m_SelectionMode(0), m_Drag(false),
    /*m_SliceRenderer(NULL),*/ m_RenderSDF(false), m_RenderControlPoints(true), m_SliceDirty(true),
    m_MouseZoomStarted(false), m_UpdateSliceOnRelease(false), m_CameraInitialized(false)
{
  // Request OpenGL compatibility profile for ARB fragment program support
  QSurfaceFormat format;
  format.setProfile(QSurfaceFormat::CompatibilityProfile);
  format.setVersion(3, 3);
  setFormat(format);
  
  /* initialize the color maps */
  m_Palette = m_ByteMap;
  memset(m_ByteMap,0,256*4);
  for(int i=0; i<256; i++)
    {
      m_GreyMap[i*4+0] = i;
      m_GreyMap[i*4+1] = i;
      m_GreyMap[i*4+2] = i;
      m_GreyMap[i*4+3] = 255;
    }
  
  setMouseTracking(true);
  setFocusPolicy(Qt::ClickFocus);
  
  m_TexCoordMinX = 0.0; m_TexCoordMaxX = 1.0; m_TexCoordMinY = 0.0; m_TexCoordMaxY = 1.0;
  m_VertCoordMinX = -0.5; m_VertCoordMaxX = 0.5; m_VertCoordMinY = -0.5; m_VertCoordMaxY = 0.5;

  gsl_set_error_handler_off();

  // Use middle button for camera zoom (Qt6 doesn't support button combinations)
  setMouseBinding(Qt::NoModifier, Qt::MiddleButton, CAMERA, ZOOM);
}

// TODO: Qt6 - QGLFormat constructor removed, use the standard constructor instead
#if 0
SliceCanvas::SliceCanvas(SliceCanvas::SliceAxis a, const QGLFormat& format, QWidget * parent, const char * name)
  : QGLViewer(parent), m_Depth(0), m_PointClassList(NULL), 
    m_CurrentPointClass(NULL), m_SliceAxis(a), m_PointSize(1), m_Variable(0), m_Timestep(0), m_Drawable(true),
    m_SliceTiles(NULL), m_NumSliceTiles(0), m_MouseIsDown(false), m_SelectionMode(0), m_Drag(false),
    /*m_SliceRenderer(NULL),*/ m_RenderSDF(false), m_RenderControlPoints(true), m_SliceDirty(true),
    m_MouseZoomStarted(false), m_UpdateSliceOnRelease(false)
{
  /* initialize the color maps */
  m_Palette = m_ByteMap;
  memset(m_ByteMap,0,256*4);
  for(int i=0; i<256; i++)
    {
      m_GreyMap[i*4+0] = i;
      m_GreyMap[i*4+1] = i;
      m_GreyMap[i*4+2] = i;
      m_GreyMap[i*4+3] = 255;
    }
  
  setMouseTracking(true);
  setFocusPolicy(Qt::ClickFocus);
  
  m_TexCoordMinX = 0.0; m_TexCoordMaxX = 1.0; m_TexCoordMinY = 0.0; m_TexCoordMaxY = 1.0;
  m_VertCoordMinX = -0.5; m_VertCoordMaxX = 0.5; m_VertCoordMinY = -0.5; m_VertCoordMaxY = 0.5;

  gsl_set_error_handler_off();

  // Use middle button for camera zoom (Qt6 doesn't support button combinations)
  setMouseBinding(Qt::NoModifier, Qt::MiddleButton, CAMERA, ZOOM);
}
#endif

SliceCanvas::~SliceCanvas()
{
  unsetVolume();
}

//This function sets the current shown slice according to the direction being viewed at (XY,XZ,ZY)
//Changes:
//2010.09.17 - Joe R. - moving update slice logic to its own function
void SliceCanvas::setDepth(int d)
{
  m_Depth = static_cast<unsigned int>(d);

  m_SelectedPoints.clear();
  m_SDF.reset();

  m_SliceDirty = true;

  update();
  emit depthChanged(static_cast<int>(m_Depth));
}

void SliceCanvas::setTransferFunction(unsigned char *func)
{
  if(!m_Drawable) return;
  memcpy(m_ByteMap,func,256*4);
  makeCurrent();
  uploadColorTable();
  update();
}

void SliceCanvas::addPoint(unsigned int x, unsigned int y, unsigned int z)
{
  if(m_CurrentPointClass)
    {
      m_CurrentPointClass->addPoint(x,y,z);
      emit pointAdded(m_Variable,m_Timestep,m_CurrentPointClassIndex,x,y,z);
    }
}

void SliceCanvas::removePoint(unsigned int x, unsigned int y, unsigned int z)
{
  if(m_CurrentPointClass)
    {
      m_CurrentPointClass->removePoint(x,y,z);
      emit pointRemoved(m_Variable,m_Timestep,m_CurrentPointClassIndex,x,y,z);
    }
}

void SliceCanvas::setCurrentClass(int index)
{
  if(m_PointClassList)
    {
      if( m_PointClassList[m_Variable][m_Timestep]->size() <= index ) {
	if (m_PointClassList[m_Variable][m_Timestep]->size() == 0) {
	  // when there are no point classes, do nothing and don't complain so there
	  // are not messages on startup
	} else {
         fprintf( stderr, "SliceCanvas::setCurrentClass(): point class list index is out of range\n");
	}
      }
      else {
         m_CurrentPointClass = m_PointClassList[m_Variable][m_Timestep]->at(index);
         m_CurrentPointClassIndex = index;
      }
    }	
}

void SliceCanvas::setVolume(const VolMagick::VolumeFileInfo& vfi)
{
  float aspect_x,aspect_y,slicex,slicey,distx,disty,dimx,dimy;
  unsigned int imgx,imgy,imgz;
  unsetVolume();
  if(!vfi.isSet()) return;
  m_VolumeFileInfo = vfi;
  imgx = upToPowerOfTwoLocal(m_VolumeFileInfo.XDim());
  imgy = upToPowerOfTwoLocal(m_VolumeFileInfo.YDim());
  imgz = upToPowerOfTwoLocal(m_VolumeFileInfo.ZDim());

  m_TexCoordMinX = 0.0; m_TexCoordMaxX = 1.0; m_TexCoordMinY = 0.0; m_TexCoordMaxY = 1.0;
  m_VertCoordMinX = -0.5; m_VertCoordMaxX = 0.5; m_VertCoordMinY = -0.5; m_VertCoordMaxY = 0.5;
 
  setDepth(0);
  updateSlice();
  resetView();


}

void SliceCanvas::unsetVolume()
{
  m_VolumeFileInfo = VolMagick::VolumeFileInfo();
  m_Slice.reset();
  if(m_SliceTiles)
    {
      makeCurrent();
      for(unsigned int i=0; i<m_NumSliceTiles; i++) glDeleteTextures(1,&(m_SliceTiles[i].texture));
      delete [] m_SliceTiles;
    }
  m_SliceTiles = NULL;
  m_NumSliceTiles = 0;
  m_PointClassList = NULL;
  m_CurrentPointClass = NULL;
  m_Depth = 0;
  m_Variable = 0;
  m_Timestep = 0;
  m_SDF.reset();
}

void SliceCanvas::resetView()
{
//  VolMagick::BoundingBox bbox;

  Vec minPt, maxPt;

  switch(m_SliceAxis)
    {
    case XY:
      {
	  
        minPt = 
          Vec(m_VolumeFileInfo.XMin(),
              m_VolumeFileInfo.YMin(),
              -0.5);
        maxPt =
          Vec(m_VolumeFileInfo.XMax(),
              m_VolumeFileInfo.YMax(),
              0.5);
        break;
      }
    case XZ:
      {
        minPt = 
          Vec(m_VolumeFileInfo.XMin(),
              m_VolumeFileInfo.ZMin(),
              -0.5);
        maxPt =
          Vec(m_VolumeFileInfo.XMax(),
              m_VolumeFileInfo.ZMax(),
              0.5);
        break;
      }
    case ZY:
      {
        minPt = 
          Vec(m_VolumeFileInfo.ZMin(),
              m_VolumeFileInfo.YMin(),
              -0.5);
        maxPt =
          Vec(m_VolumeFileInfo.ZMax(),
              m_VolumeFileInfo.YMax(),
              0.5);
        break;
      }
    }

  camera()->fitBoundingBox(minPt,maxPt);
  update();
}

void SliceCanvas::setCellMarkingMode(int m)
{
  m_CellMarkingMode = CellMarkingMode(m);
}

void SliceCanvas::setCurrentContour(const std::string& name)
{
  m_CurrentContour = m_Contours[m_Variable][m_Timestep][name];
}

/*
void SliceCanvas::setInterpolationType(int interp)
{
  if(m_CurrentContour != NULL) m_CurrentContour->interpolationType(SurfRecon::InterpolationType(interp));
  update();
}
*/

Vec SliceCanvas::getObjectSpaceCoords(const Vec& screenCoords) const
{
  Vec objCoords;
  Vec coords = camera()->unprojectedCoordinatesOf(screenCoords);
  
  //map the slice vertex coordinates to actual bounding box coordinates
  switch(m_SliceAxis)
    {
    case XY:
      {
	objCoords.x = ((coords.x-m_VertCoordMinX)/(m_VertCoordMaxX-m_VertCoordMinX))*
	  (m_VolumeSlice.XMax()-m_VolumeSlice.XMin())+m_VolumeSlice.XMin();
	objCoords.y = ((coords.y-m_VertCoordMinY)/(m_VertCoordMaxY-m_VertCoordMinY))*
	  (m_VolumeSlice.YMax()-m_VolumeSlice.YMin())+m_VolumeSlice.YMin();
	objCoords.z = m_VolumeSlice.ZDim()-1 == 0 ? m_VolumeSlice.ZMin() : 
	  (double(m_Depth)/(double(m_VolumeSlice.ZDim())-1.0))*
	  (m_VolumeSlice.ZMax()-m_VolumeSlice.ZMin())+m_VolumeSlice.ZMin();
      }
      break;
    case XZ:
      {
	objCoords.x = ((coords.x-m_VertCoordMinX)/(m_VertCoordMaxX-m_VertCoordMinX))*
	  (m_VolumeSlice.XMax()-m_VolumeSlice.XMin())+m_VolumeSlice.XMin();
	objCoords.y = m_VolumeSlice.YDim()-1 == 0 ? m_VolumeSlice.YMin() : 
	  (double(m_Depth)/(double(m_VolumeSlice.YDim())-1.0))*
	  (m_VolumeSlice.YMax()-m_VolumeSlice.YMin())+m_VolumeSlice.YMin();
	objCoords.z = ((coords.y-m_VertCoordMinY)/(m_VertCoordMaxY-m_VertCoordMinY))*
	  (m_VolumeSlice.ZMax()-m_VolumeSlice.ZMin())+m_VolumeSlice.ZMin();
      }
      break;
    case ZY:
      {
	objCoords.x = m_VolumeSlice.XDim()-1 == 0 ? m_VolumeSlice.XMin() : 
	  (double(m_Depth)/(double(m_VolumeSlice.XDim())-1.0))*
	  (m_VolumeSlice.XMax()-m_VolumeSlice.XMin())+m_VolumeSlice.XMin();
	objCoords.y = ((coords.y-m_VertCoordMinY)/(m_VertCoordMaxY-m_VertCoordMinY))*
	  (m_VolumeSlice.YMax()-m_VolumeSlice.YMin())+m_VolumeSlice.YMin();
	objCoords.z = ((coords.x-m_VertCoordMinX)/(m_VertCoordMaxX-m_VertCoordMinX))*
	  (m_VolumeSlice.ZMax()-m_VolumeSlice.ZMin())+m_VolumeSlice.ZMin();
      }
      break;
    default: break;
    }

  return objCoords;
}

void SliceCanvas::updateSlice()
{
  unsigned int i,j,imgx,imgy,imgz,size;
  float aspect_x,aspect_y,slicex,slicey,distx,disty,dimx,dimy;

  if(!m_VolumeFileInfo.isSet()) return;

  m_TexCoordMinX = 0.0; m_TexCoordMaxX = 1.0; m_TexCoordMinY = 0.0; m_TexCoordMaxY = 1.0;
  m_VertCoordMinX = -0.5; m_VertCoordMaxX = 0.5; m_VertCoordMinY = -0.5; m_VertCoordMaxY = 0.5;



  /*
    // 9-20-2011
    // arand: commenting this section: it is causing a compile issue...
  //Set the new HDF5 CVC file handler to use a higher dimension for slices.
  //This is kind of ugly, streamline this later!
  int maxdim = 256;
  if(!VolMagick::VolumeFile_IO::handlerMap()[".cvc"].empty())
    VolMagick::VolumeFile_IO::handlerMap()[".cvc"][0]->
      maxdim(VolMagick::Dimension(maxdim,maxdim,maxdim));
  */


  //Window edges in world coordinates
  Vec topRight = 
    camera()->unprojectedCoordinatesOf(Vec(width(),0.0,0.0));
  Vec bottomLeft = 
    camera()->unprojectedCoordinatesOf(Vec(0.0,height(),0.0));

  switch(m_SliceAxis)
    {
    case XY:
      {
	m_Depth = m_Depth <= m_VolumeFileInfo.ZDim()-1 ? m_Depth : m_VolumeFileInfo.ZDim()-1;
 
   // cout<<"Xmin, boleft, xmax, tor: " << m_VolumeFileInfo.XMin()<< " " <<  (double(bottomLeft.x)-m_VolumeFileInfo.XSpan()) <<"  "<< m_VolumeFileInfo.XMax() << " " << ( double(topRight.x)+m_VolumeFileInfo.XSpan()) << endl;

	VolMagick::BoundingBox bbox(/*std::max(*/m_VolumeFileInfo.XMin(),
                                            // double(bottomLeft.x)-m_VolumeFileInfo.XSpan()), //a little bit outside the screen so edges render correctly 
				    /*std::max(*/m_VolumeFileInfo.YMin(),
                                 //            double(bottomLeft.y)-m_VolumeFileInfo.YSpan()),
				    m_VolumeFileInfo.ZMin() + m_Depth*m_VolumeFileInfo.ZSpan(),
				    /*std::min(*/m_VolumeFileInfo.XMax(),
                                 //            double(topRight.x)+m_VolumeFileInfo.XSpan()),
				    /*std::min(*/m_VolumeFileInfo.YMax(),
                                  //           double(topRight.y)+m_VolumeFileInfo.YSpan()),
				    m_VolumeFileInfo.ZMin() + m_Depth*m_VolumeFileInfo.ZSpan());

	/* read a single slice into m_VolumeSlice */
	VolMagick::readVolumeFile(cvcapp, m_VolumeSlice,
				  m_VolumeFileInfo.filename(),
				  m_Variable,m_Timestep,
				  bbox);
				  
	/* duplicate the slice and map it to unsigned char range */
	VolMagick::Volume mappedSlice(m_VolumeSlice);
	if(mappedSlice.voxelType() != VolMagick::UChar)
	  {
	    /* get min and max from overall file, else it will be calculated from subvolume */
	    mappedSlice.min(m_VolumeFileInfo.min(m_Variable,m_Timestep));
	    mappedSlice.max(m_VolumeFileInfo.max(m_Variable,m_Timestep));
	    mappedSlice.map(0.0,255.0);
	    mappedSlice.voxelType(VolMagick::UChar);
	  }
      
        //initialize slice upload buffer
        imgx = upToPowerOfTwoLocal(m_VolumeSlice.XDim());
        imgy = upToPowerOfTwoLocal(m_VolumeSlice.YDim());
        imgz = upToPowerOfTwoLocal(m_VolumeSlice.ZDim());
        size = imgx*imgy;
        m_Slice.reset(new unsigned char[size]);

	/* now copy it to the uploadable slice buffer */
	for(i=0; i<m_VolumeSlice.XDim(); i++)
	  for(j=0; j<m_VolumeSlice.YDim(); j++)
	    m_Slice[i + imgx*j] = mappedSlice(i,j,0);

        // Set up quad verts for this slice
        m_VertCoordMinX = bbox.XMin();
        m_VertCoordMaxX = bbox.XMax();
        m_VertCoordMinY = bbox.YMin();
        m_VertCoordMaxY = bbox.YMax();
      }
      break;
    case XZ:
      {
	m_Depth = m_Depth <= m_VolumeFileInfo.YDim()-1 ? m_Depth : m_VolumeFileInfo.YDim()-1;

	VolMagick::BoundingBox bbox(/*std::max(*/m_VolumeFileInfo.XMin(),
                                            // double(bottomLeft.x)-m_VolumeFileInfo.XSpan()),
				    m_VolumeFileInfo.YMin() + m_Depth*m_VolumeFileInfo.YSpan(),
				    /*std::max(*/m_VolumeFileInfo.ZMin(),
                                 //            double(bottomLeft.y)-m_VolumeFileInfo.ZSpan()),
				    /*std::min(*/m_VolumeFileInfo.XMax(),
                                 //            double(topRight.x)+m_VolumeFileInfo.XSpan()),
				    m_VolumeFileInfo.YMin() + m_Depth*m_VolumeFileInfo.YSpan(),
				   /* std::min(*/m_VolumeFileInfo.ZMax()//,
//                                             double(topRight.y)+m_VolumeFileInfo.ZSpan())
								 );
      
	/* read a single slice into m_VolumeSlice */
	VolMagick::readVolumeFile(cvcapp, m_VolumeSlice,
				  m_VolumeFileInfo.filename(),
				  m_Variable,m_Timestep,
				  bbox);
				  
	/* duplicate the slice and map it to the unsigned char range */
	VolMagick::Volume mappedSlice(m_VolumeSlice);
	if(mappedSlice.voxelType() != VolMagick::UChar)
	  {
	    /* get min and max from overall file, else it will be calculated from subvolume */
	    mappedSlice.min(m_VolumeFileInfo.min(m_Variable,m_Timestep));
	    mappedSlice.max(m_VolumeFileInfo.max(m_Variable,m_Timestep));
	    mappedSlice.map(0.0,255.0);
	    mappedSlice.voxelType(VolMagick::UChar);
	  }
      
        //initialize slice upload buffer
        imgx = upToPowerOfTwoLocal(m_VolumeSlice.XDim());
        imgy = upToPowerOfTwoLocal(m_VolumeSlice.YDim());
        imgz = upToPowerOfTwoLocal(m_VolumeSlice.ZDim());
        size = imgx*imgz;
        m_Slice.reset(new unsigned char[size]);

	/* now copy it to the uploadable slice buffer */
	for(i=0; i<m_VolumeSlice.XDim(); i++)
	  for(j=0; j<m_VolumeSlice.ZDim(); j++)
	    m_Slice[i + imgx*j] = mappedSlice(i,0,j);

        // Set up quad verts for this slice
        m_VertCoordMinX = bbox.XMin();
        m_VertCoordMaxX = bbox.XMax();
        m_VertCoordMinY = bbox.ZMin();
        m_VertCoordMaxY = bbox.ZMax();
      }
      break;
    case ZY:
      {
	m_Depth = m_Depth <= m_VolumeFileInfo.XDim()-1 ? m_Depth : m_VolumeFileInfo.XDim()-1;

	VolMagick::BoundingBox bbox(m_VolumeFileInfo.XMin() + m_Depth*m_VolumeFileInfo.XSpan(),
				    /*std::max(*/m_VolumeFileInfo.YMin(),
                                 //            double(bottomLeft.x)-m_VolumeFileInfo.YSpan()),
				    /*std::max(*/m_VolumeFileInfo.ZMin(),
                                 //            double(bottomLeft.y)-m_VolumeFileInfo.ZSpan()),
				    m_VolumeFileInfo.XMin() + m_Depth*m_VolumeFileInfo.XSpan(),
				   /* std::min(*/m_VolumeFileInfo.YMax(),
                                 //            double(topRight.x)+m_VolumeFileInfo.YSpan()),
				   /* std::min(*/m_VolumeFileInfo.ZMax() //,
//                                             double(topRight.y)+m_VolumeFileInfo.ZSpan())
								 );
      
	/* read a single slice into m_VolumeSlice */
	VolMagick::readVolumeFile(cvcapp, m_VolumeSlice,
				  m_VolumeFileInfo.filename(),
				  m_Variable,m_Timestep,
				  bbox);
				  
	/* duplicate the slice and map it to the unsigned char range */
	VolMagick::Volume mappedSlice(m_VolumeSlice);
	if(mappedSlice.voxelType() != VolMagick::UChar)
	  {
	    /* get min and max from overall file, else it will be calculated from subvolume */
	    mappedSlice.min(m_VolumeFileInfo.min(m_Variable,m_Timestep));
	    mappedSlice.max(m_VolumeFileInfo.max(m_Variable,m_Timestep));
	    mappedSlice.map(0.0,255.0);
	    mappedSlice.voxelType(VolMagick::UChar);
	  }

        //initialize slice upload buffer
        imgx = upToPowerOfTwoLocal(m_VolumeSlice.XDim());
        imgy = upToPowerOfTwoLocal(m_VolumeSlice.YDim());
        imgz = upToPowerOfTwoLocal(m_VolumeSlice.ZDim());
        size = imgz*imgy;
        m_Slice.reset(new unsigned char[size]);
      
	/* now copy it to the uploadable slice buffer */
	for(i=0; i<m_VolumeSlice.ZDim(); i++)
	  for(j=0; j<m_VolumeSlice.YDim(); j++)
	    m_Slice[i + imgz*j] = mappedSlice(0,j,i);

        // Set up quad verts for this slice
        m_VertCoordMinX = bbox.ZMin();
        m_VertCoordMaxX = bbox.ZMax();
        m_VertCoordMinY = bbox.YMin();
        m_VertCoordMaxY = bbox.YMax();
      }
      break;
    }

  /* calculate the texture rectangle dimensions */
  aspect_x = m_SliceAxis == XY ? m_VolumeSlice.XDim() 
    : m_SliceAxis == XZ ? m_VolumeSlice.XDim() 
    : m_SliceAxis == ZY ? m_VolumeSlice.ZDim() : 0;
  aspect_y = m_SliceAxis == XY ? m_VolumeSlice.YDim() 
    : m_SliceAxis == XZ ? m_VolumeSlice.ZDim() 
    : m_SliceAxis == ZY ? m_VolumeSlice.YDim() : 0;
  dimx = m_SliceAxis == XY ? imgx
    : m_SliceAxis == XZ ? imgx
    : m_SliceAxis == ZY ? imgz : 0;
  dimy = m_SliceAxis == XY ? imgy 
    : m_SliceAxis == XZ ? imgz
    : m_SliceAxis == ZY ? imgy : 0;
  aspect_x /= dimx;
  aspect_y /= dimy;
  m_TexCoordMaxX *= aspect_x;
  m_TexCoordMaxY *= aspect_y;

  uploadSlice(m_Slice.get(), m_SliceAxis == XY ? imgx : m_SliceAxis == XZ ? imgx : m_SliceAxis == ZY ? imgz : 0,
	      m_SliceAxis == XY ? imgy : m_SliceAxis == XZ ? imgz : m_SliceAxis == ZY ? imgy : 0);
  
  // Update scene bounding box for proper camera framing
  float centerX = (m_VertCoordMinX + m_VertCoordMaxX) / 2.0;
  float centerY = (m_VertCoordMinY + m_VertCoordMaxY) / 2.0;
  float width = m_VertCoordMaxX - m_VertCoordMinX;
  float height = m_VertCoordMaxY - m_VertCoordMinY;
  float radius = std::sqrt(width*width + height*height) / 2.0;
  
  setSceneCenter(Vec(centerX, centerY, 0.0));
  setSceneRadius(radius);
  
  // Only reset camera view on first slice load, not when changing slices
  if(!m_CameraInitialized && m_VolumeFileInfo.isSet()) {
    camera()->showEntireScene();
    m_CameraInitialized = true;
  }
  
  m_SliceDirty = false;
}

void SliceCanvas::mouseMoveEvent(QMouseEvent *event)
{
  int x,y,realx,realy,realz;
  unsigned char *slice = m_Slice.get();
  unsigned int imgx, imgy;
  
  if(m_SelectionMode != 0)
    {
      m_SelectionRectangle.setX(event->x());
      m_SelectionRectangle.setY(event->y());
      update();
    }

  //If we started to zoom, tell the slice canvas to update the slice when we're done
  //zooming so we can get a new level of detail from VolMagick if possible.
  if(m_MouseZoomStarted)
    m_UpdateSliceOnRelease = true;

  if(m_Slice && m_VolumeFileInfo.isSet())
    {
      //update cell info
      {
	Vec newp = getObjectSpaceCoords(Vec(event->x(),event->y(),0.0));

	Vec coords = camera()->unprojectedCoordinatesOf(Vec(event->x(),event->y(),0.0));
	//printf("coords.x: %f, coords.y: %f, coords.z: %f\n",coords.x,coords.y,coords.z);
	switch(m_SliceAxis)
	  {
	  case XY: imgx = m_VolumeSlice.XDim(); imgy = m_VolumeSlice.YDim(); break;
	  case XZ: imgx = m_VolumeSlice.XDim(); imgy = m_VolumeSlice.ZDim(); break;
	  case ZY: imgx = m_VolumeSlice.ZDim(); imgy = m_VolumeSlice.YDim(); break;
	  default: return;
	  }
    
	if(m_VertCoordMinX <= coords.x && coords.x < m_VertCoordMaxX && /* dont count the furthest edge as part of the slice */
	   m_VertCoordMinY <= coords.y && coords.y < m_VertCoordMaxY)
	  {
	    x = int(((coords.x-m_VertCoordMinX)/(m_VertCoordMaxX-m_VertCoordMinX))*float(imgx));
	    y = int(((coords.y-m_VertCoordMinY)/(m_VertCoordMaxY-m_VertCoordMinY))*float(imgy));
	  }
	else
	  {
	    x = 0;
	    y = 0;
	  }
	//printf("x: %d, y: %d\n",x,y);
	//printf("coords.x: %f, coords.y: %f\n",coords.x,coords.y);
	slice+=x+y*upToPowerOfTwoLocal(imgx);
    
	switch(m_SliceAxis)
	  {
	  case XY: realx = x; realy = y; realz = m_Depth; break;
	  case XZ: realx = x; realy = m_Depth; realz = y; break;
	  case ZY: realx = m_Depth; realy = y; realz = x; break;
	  default: return;
	  }
    
	/* make sure realx,realy,realz are within the volume dimensions... sometimes the mouse location can be right on the
	   edge of the slice */
	if(realx >= int(m_VolumeSlice.XDim())) realx = int(m_VolumeSlice.XDim()-1);
	if(realy >= int(m_VolumeSlice.YDim())) realy = int(m_VolumeSlice.YDim()-1);
	if(realz >= int(m_VolumeSlice.ZDim())) realz = int(m_VolumeSlice.ZDim()-1);
    
	emit mouseOver(realx,realy,realz,
		       newp.x,newp.y,newp.z,
		       m_ByteMap[(*slice)*4+0],m_ByteMap[(*slice)*4+1],m_ByteMap[(*slice)*4+2],m_ByteMap[(*slice)*4+3],
		       m_SliceAxis == XY ? m_VolumeSlice(realx,realy,0) :
		       m_SliceAxis == XZ ? m_VolumeSlice(realx,0,realz) :
		       m_SliceAxis == ZY ? m_VolumeSlice(0,realy,realz) : 0.0);
      }

      //translate selected points if translation keys are pressed
      if(m_Drag)
	{
	  SurfRecon::PointPtr p;
	  Vec orig = getObjectSpaceCoords(Vec(m_LastPoint.x(),m_LastPoint.y(),0.0));
	  Vec newp = getObjectSpaceCoords(Vec(event->x(),event->y(),0.0));
	  Vec diff = newp - orig;

	  for(std::set<SurfRecon::WeakPointPtr>::iterator cur = m_SelectedPoints.begin();
	      cur != m_SelectedPoints.end();
	      cur++)
	    if(p = (*cur).lock())
	      *p = *p + SurfRecon::Vector(diff.x,diff.y,diff.z);
	  
	  m_LastPoint = event->pos();

	  update();
	}
    }

  QGLViewer::mouseMoveEvent(event);
}

void SliceCanvas::mouseDoubleClickEvent(QMouseEvent *event)
{
  if(m_Slice && m_VolumeFileInfo.isSet())
    {
      if(m_CurrentPointClass != NULL && m_CellMarkingMode == PointClassMarking)
	{
	  int x,y;
	  unsigned char *slice = m_Slice.get();
	  unsigned int imgx, imgy, realx, realy, realz;
	  Vec coords = camera()->unprojectedCoordinatesOf(Vec(event->x(),event->y(),0.0));
	  //printf("coords.x: %f, coords.y: %f, coords.z: %f\n",coords.x,coords.y,coords.z);
	  switch(m_SliceAxis)
	    {
	    case XY: imgx = m_VolumeSlice.XDim(); imgy = m_VolumeSlice.YDim(); break;
	    case XZ: imgx = m_VolumeSlice.XDim(); imgy = m_VolumeSlice.ZDim(); break;
	    case ZY: imgx = m_VolumeSlice.ZDim(); imgy = m_VolumeSlice.YDim(); break;
	    default: return;
	    }
	  
	  if(m_VertCoordMinX <= coords.x && coords.x < m_VertCoordMaxX && /* dont count the furthest edge as part of the slice */
	     m_VertCoordMinY <= coords.y && coords.y < m_VertCoordMaxY)
	    {
	      x = int(((coords.x-m_VertCoordMinX)/(m_VertCoordMaxX-m_VertCoordMinX))*float(imgx));
	      y = int(((coords.y-m_VertCoordMinY)/(m_VertCoordMaxY-m_VertCoordMinY))*float(imgy));
	      
	      cvcapp.log(5, boost::str(boost::format("x: %d, y: %d") % x % y));
	      cvcapp.log(5, boost::str(boost::format("coords.x: %f, coords.y: %f")%coords.x%coords.y));
	      slice+=x+y*upToPowerOfTwoLocal(imgx);
	      
	      switch(m_SliceAxis)
		{
		case XY: realx = x; realy = y; realz = m_Depth; break;
		case XZ: realx = x; realy = m_Depth; realz = y; break;
		case ZY: realx = m_Depth; realy = y; realz = x; break;
		default: return;
		}
	      
	      cvcapp.log(5, boost::str(boost::format("realx: %d, realy: %d, realz: %d\n")%realx%realy%realz));
	      
	      /* check to see if we should remove a point */
              if( !m_CurrentPointClass )
                fprintf( stderr, "current class point is NULL\n");

	      GridPoint *gp = m_CurrentPointClass->getClosestPoint(realx,realy,realz, 
								   m_SliceAxis == XY ? PointClass::XY :
								   m_SliceAxis == XZ ? PointClass::XZ :
								   m_SliceAxis == ZY ? PointClass::ZY : PointClass::ANY);
	      if((gp != NULL) && 
		 (sqrt((double)(gp->x - realx)*(gp->x - realx)+(gp->y - realy)*(gp->y - realy)+(gp->z - realz)*(gp->z - realz)) < m_PointSize) &&
		 (m_SliceAxis == XY ? gp->z == realz : m_SliceAxis == XZ ? gp->y == realy : m_SliceAxis == ZY ? gp->x == realx : false))
		removePoint(gp->x,gp->y,gp->z); /* make sure we only delete points that are visible (i.e. on the current slice) */
	      else
		addPoint(realx,realy,realz);
	      update();
	    }
	}
#if 0
      else if(m_CurrentContour != NULL && m_CellMarkingMode == ContourMarking)
	{
	  Vec coords = camera()->unprojectedCoordinatesOf(Vec(event->x(),event->y(),0.0));
	  //double realx, realy, realz;
	  Vec realCoords;

	  if(m_VertCoordMinX <= coords.x && coords.x < m_VertCoordMaxX && /* dont count the furthest edge as part of the slice */
	     m_VertCoordMinY <= coords.y && coords.y < m_VertCoordMaxY)
	    {
	      realCoords = getObjectSpaceCoords(Vec(event->x(),event->y(),0.0));

	      //cvcapp.log(5, "Adding curve node");
	      
	      SurfRecon::PointPtr p;

	      m_CurrentContour->add(p = SurfRecon::PointPtr(new SurfRecon::Point(realCoords.x,
										 realCoords.y,
										 realCoords.z)), 
				    m_Depth, 
				    SurfRecon::Orientation(m_SliceAxis));
	      
	      m_SelectedPoints.clear();
	      m_SelectedPoints.insert(SurfRecon::WeakPointPtr(p));

	      //cvcapp.log(5, "Current curve size: %d",int(m_CurrentContour->currentCurvePoints().size()));
	    }
	}
#endif
    }
  
  QGLViewer::mouseDoubleClickEvent(event);
}

void SliceCanvas::mousePressEvent(QMouseEvent *event)
{
  // Start selection. Mode is ADD with Shift key and TOGGLE with Alt key.
  m_SelectionRectangle = QRect(event->pos(), event->pos());
  if(event->button() == Qt::LeftButton)
    {
      Qt::KeyboardModifiers mods = event->modifiers();
      if(mods & Qt::ShiftModifier && (mods & (Qt::ControlModifier | Qt::AltModifier)))
      {
        // we dont want to select in this case.. instead we want to translate selected points
        m_SelectionMode = 0;
        m_Drag = true;
        m_LastPoint = event->pos();
      }
      else if(mods & Qt::ShiftModifier)
        m_SelectionMode = 2; //add to selection
      else if(mods & (Qt::ControlModifier | Qt::AltModifier))
        m_SelectionMode = 3; //remove from selection
      else
        m_SelectionMode = 1; //new selection
      return;
    }

  if(event->button() == Qt::MiddleButton ||
     event->button() == Qt::RightButton)
    m_MouseZoomStarted = true;
  
  QGLViewer::mousePressEvent(event);
}

void SliceCanvas::mouseReleaseEvent(QMouseEvent *event)
{
//  m_Drag = false;
//  m_MouseZoomStarted = false;

//  if(m_UpdateSliceOnRelease)
//    {
 //     m_SliceDirty = true;
 //     m_UpdateSliceOnRelease = false;
 //     update();
 //   }

#if 0
  if(m_RenderControlPoints)
    {
      m_SelectionMode = 0;
      return; //dont do anything if user cannot see control points
    }
#endif

  if(m_SelectionMode != 0)
    {
      // Backup the current selection
      std::set<SurfRecon::WeakPointPtr> previousSelection(m_SelectedPoints);

      m_SelectionRectangle = m_SelectionRectangle.normalized();
      // Define selection window dimensions
      setSelectRegionWidth(m_SelectionRectangle.width());
      setSelectRegionHeight(m_SelectionRectangle.height());
      // Compute rectangle center and perform selection
      select(m_SelectionRectangle.center());
      // Update display to show new selected objects

      m_SelectionMode = 0;

      // If no new points were selected, add a point
      if(m_SelectedPoints.empty())
	{
	  //m_SelectedPoints = previousSelection;

	  if(m_CurrentContour != NULL && m_CellMarkingMode == ContourMarking)
	    {
	      Vec coords = camera()->unprojectedCoordinatesOf(Vec(event->x(),event->y(),0.0));
	      //double realx, realy, realz;
	      Vec realCoords;

	      /* dont count the furthest edge as part of the slice */
	      if(m_VertCoordMinX <= coords.x && coords.x < m_VertCoordMaxX &&
		 m_VertCoordMinY <= coords.y && coords.y < m_VertCoordMaxY)
		{
		  realCoords = getObjectSpaceCoords(Vec(event->x(),event->y(),0.0));

		  //cvcapp.log(5, "Adding curve node");
	      
		  SurfRecon::PointPtr p;
		  SurfRecon::CurvePtr c;
		  bool prepend = false;

		  //check to see if we've previously selected a first or last point on a curve
		  if(previousSelection.size() == 1)
		    {
		      // Determine if the selected point is a first or last point of a curve.
		      // If it is, then add a new point to that curve
		      for(SurfRecon::ContourPtrMap::iterator i = m_Contours[m_Variable][m_Timestep].begin();
			  i != m_Contours[m_Variable][m_Timestep].end();
			  i++)
			{
			  if((*i).second == NULL) continue;
			  for(std::vector<SurfRecon::CurvePtr>::iterator cur = (*i).second->curves().begin();
			      cur != (*i).second->curves().end();
			      cur++)
			    {
			      if(!boost::get<2>(**cur).empty())
				{
				  if(previousSelection.find(*boost::get<2>(**cur).begin()) != previousSelection.end())
				    {
				      c = *cur;
				      prepend = true;
				    }
				  else if(previousSelection.find(*boost::get<2>(**cur).rbegin()) != previousSelection.end())
				    {
				      c = *cur;
				      prepend = false;
				    }
				}
			    }
			}

		      if(c != NULL) m_CurrentContour->currentCurve(c); //set the current curve to be the one we selected
		    }

		  //if we haven't previously selected a first or last point of a curve, start a new curve
		  if(c == NULL)
		    m_CurrentContour->add(SurfRecon::CurvePtr(new SurfRecon::Curve()));
		  
		  if(prepend)
		    m_CurrentContour->prepend(p = SurfRecon::PointPtr(new SurfRecon::Point(realCoords.x,
											   realCoords.y,
											   realCoords.z)), 
					      m_Depth, 
					      SurfRecon::Orientation(m_SliceAxis));
		  else
		    m_CurrentContour->add(p = SurfRecon::PointPtr(new SurfRecon::Point(realCoords.x,
										       realCoords.y,
										       realCoords.z)), 
					  m_Depth, 
					  SurfRecon::Orientation(m_SliceAxis));
	      
		  //set the current selection to the point we just added
		  m_SelectedPoints.insert(SurfRecon::WeakPointPtr(p));
		}
	    }
	}

      //If the front and back points of a curve have been selected, close the loop.
      //If a loop has already been closed, it cannot be added to.  Drawing the first and last points
      //in the same location has the effect that selecting 1 selects both, preventing the following code's execution.
      else if(m_SelectedPoints.size() == 1 && previousSelection.size() == 1)
	{
	  SurfRecon::CurvePtr c[2];
	  SurfRecon::PointPtr p;

	  // Determine if the selected point is a first or last point of a curve.
	  for(SurfRecon::ContourPtrMap::iterator i = m_Contours[m_Variable][m_Timestep].begin();
	      i != m_Contours[m_Variable][m_Timestep].end();
	      i++)
	    {
	      if((*i).second == NULL) continue;
	      for(std::vector<SurfRecon::CurvePtr>::iterator cur = (*i).second->curves().begin();
		  cur != (*i).second->curves().end();
		  cur++)
		{
		  if(!boost::get<2>(**cur).empty())
		    {
		      if((previousSelection.find(*boost::get<2>(**cur).begin()) != previousSelection.end()) ||
			 (previousSelection.find(*boost::get<2>(**cur).rbegin()) != previousSelection.end()))
			{
			  c[0] = *cur;
			}
		      if((m_SelectedPoints.find(*boost::get<2>(**cur).begin()) != m_SelectedPoints.end()) ||
			 (m_SelectedPoints.find(*boost::get<2>(**cur).rbegin()) != m_SelectedPoints.end()))
			{
			  c[1] = *cur;
			}
		    }
		}
	    }

	  if(c[0] != NULL &&
	     c[0] == c[1] &&
	     m_SelectedPoints.find(*previousSelection.begin()) == m_SelectedPoints.end())
	    {
	      m_SelectedPoints.clear();
	      //add a copy of the first point of the curve to the end of the curve
	      boost::get<2>(*c[0]).push_back(p = SurfRecon::PointPtr(new SurfRecon::Point(**boost::get<2>(*c[0]).begin())));
	      m_SelectedPoints.insert(SurfRecon::WeakPointPtr(p));
	    }
	}

      update();

      return;
    }

  QGLViewer::mouseReleaseEvent(event);
}

void SliceCanvas::keyPressEvent(QKeyEvent *event)
{
  if(m_VolumeFileInfo.isSet())
    {
      unsigned int dim = 
	m_SliceAxis == XY ? 
	m_VolumeFileInfo.ZDim() :
	m_SliceAxis == XZ ?
	m_VolumeFileInfo.YDim() :
	m_VolumeFileInfo.XDim();
      
      switch(event->key())
	{
	case Qt::Key_Left:
	  setDepth(getDepth()-1 > 0 ? getDepth()-1 : 0);
	  break;
	case Qt::Key_Right:
	  setDepth(getDepth()+1 < dim ? getDepth()+1 : getDepth());
	  break;
	case Qt::Key_PageDown:
	  setDepth(getDepth()-5 > 0 ? getDepth()-5 : 0);
	  break;
	case Qt::Key_PageUp:
	  setDepth(getDepth()+5 < dim ? getDepth()+5 : getDepth());
	  break;
	case Qt::Key_Delete:
	  deleteSelected();
	  update();
	  break;
	default:
	  QGLViewer::keyPressEvent(event);
	  break;
	}
    }
  else
    QGLViewer::keyPressEvent(event);
}

void SliceCanvas::setPointSize(int r)
{
  m_PointSize = r;
  update();
}

void SliceCanvas::setGreyScale(bool set)
{
  if(!m_Drawable) return;

  if(set)
    m_Palette = m_GreyMap;
  else
    m_Palette = m_ByteMap;
	
  makeCurrent();

  uploadColorTable();

  update();
}

void SliceCanvas::setRenderSDF(bool set)
{
  if(!m_Drawable) return;

  m_RenderSDF = set;

  update();
}

void SliceCanvas::setRenderControlPoints(bool set)
{
  if(!m_Drawable) return;
  m_RenderControlPoints = set;
  //m_SelectedPoints.clear();
  update();
}

void SliceCanvas::setCurrentVariable(int var)
{
  m_Variable = static_cast<unsigned int>(var);
  setDepth(m_Depth); /* this will get the new slice and updateGL */
  setCurrentClass(0);
  setCurrentContour("");
}

void SliceCanvas::setCurrentTimestep(int time)
{
  m_Timestep = static_cast<unsigned int>(time);
  setDepth(m_Depth); /* this will get the new slice and updateGL */
  setCurrentClass(0);
  setCurrentContour("");
}

void SliceCanvas::init()
{
  GLenum err = glewInit();
  if(GLEW_OK != err)
  {
     fprintf(stderr,"Error: %s\n", glewGetErrorString(err));
  }
  fprintf(stdout,"Status: Using GLEW %s\n",glewGetString(GLEW_VERSION));

  glDisable(GL_LIGHTING);

  /* initialize the slice renderer */
  fprintf(stderr, "SliceCanvas::init(): Trying ARB Fragment Program Renderer...\n");
  m_SliceRenderer.reset(new ARBFragmentProgramSliceRenderer(this));
  if(m_SliceRenderer->init())
    {
      fprintf(stderr, "SliceCanvas::init(): SUCCESS - using ARB Fragment Program Renderer\n");
      cvcapp.log(5, "SliceCanvas::init(): using ARB Fragment Program Renderer");
      goto initscene;
    }

  fprintf(stderr, "SliceCanvas::init(): Trying Paletted Slice Renderer...\n");
  m_SliceRenderer.reset(new PalettedSliceRenderer(this));
  if(m_SliceRenderer->init())
    {
      fprintf(stderr, "SliceCanvas::init(): SUCCESS - using Paletted Slice Renderer\n");
      cvcapp.log(5, "SliceCanvas::init(): using Paletted Slice Renderer");
      goto initscene;
    }
  
  fprintf(stderr, "SliceCanvas::init(): Trying SGI Color Table Slice Renderer...\n");
  m_SliceRenderer.reset(new SGIColorTableSliceRenderer(this));
  if(m_SliceRenderer->init())
    {
      fprintf(stderr, "SliceCanvas::init(): SUCCESS - using SGI Color Table Slice Renderer\n");
      cvcapp.log(5, "SliceCanvas::init(): using SGI Color Table Slice Renderer");
      goto initscene;
    }

  fprintf(stderr, "SliceCanvas::init(): ERROR - could not find suitable slice renderer!\n");
  cvcapp.log(5, "SliceCanvas::init(): could not find suitable slice renderer!");
  if(m_SliceAxis == XY) //only pop up a message box for 1 of the slice canvases
    QMessageBox::critical(this,"Error","Cannot render slices using hardware! No suitable renderer available.",QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
  m_SliceRenderer.reset();

 initscene:
  setSceneRadius(1.0);
  camera()->setPosition(Vec(0.0,0.0,0.5));
  camera()->lookAt(sceneCenter());
  camera()->setType(Camera::ORTHOGRAPHIC);
  camera()->showEntireScene();
  
  /* Forbid rotation */
  m_Constraint.reset(new qglviewer::WorldConstraint());
  m_Constraint->setRotationConstraintType(AxisPlaneConstraint::FORBIDDEN);
  camera()->frame()->setConstraint(m_Constraint.get());

  if(m_SliceAxis == XY)
    {
      /* print some GL info */
      cvcapp.log(5, boost::str(boost::format("gl: GL Vendor: %s")%glGetString(GL_VENDOR)));
      cvcapp.log(5, boost::str(boost::format("gl: GL Renderer: %s")%glGetString(GL_RENDERER)));
      cvcapp.log(5, boost::str(boost::format("gl: GL Version: %s")%glGetString(GL_VERSION)));
      
      // In OpenGL 3.0+, glGetString(GL_EXTENSIONS) may return NULL
      // Use glGetStringi instead for modern contexts
      const GLubyte* extensions = glGetString(GL_EXTENSIONS);
      if(extensions)
        cvcapp.log(5, boost::str(boost::format("gl: GL Extensions: %s")%extensions));
      else
        {
          GLint numExtensions = 0;
          glGetIntegerv(GL_NUM_EXTENSIONS, &numExtensions);
          std::string extList;
          for(GLint i = 0; i < numExtensions; i++)
            {
              if(i > 0) extList += " ";
              extList += reinterpret_cast<const char*>(glGetStringi(GL_EXTENSIONS, i));
            }
          cvcapp.log(5, boost::str(boost::format("gl: GL Extensions: %s")%extList));
        }
      
      /* initialize GL extensions */
      cvcapp.log(5, boost::str(boost::format("gl: Initializing GL_VERSION_1_2: %s")%(glewIsSupported("GL_VERSION_1_2") ? "OK" : "FAILED")));
      cvcapp.log(5, boost::str(boost::format("gl: Initializing GL_VERSION_1_3: %s")%(glewIsSupported("GL_VERSION_1_3") ? "OK" : "FAILED")));
      cvcapp.log(5, boost::str(boost::format("gl: Initializing GL_SGIS_texture_edge_clamp: %s")%(glewIsSupported("GL_SGIS_texture_edge_clamp") ? "OK" : "FAILED")));
      cvcapp.log(5, boost::str(boost::format("gl: Initializing GL_EXT_texture3D: %s")%(glewIsSupported("GL_EXT_texture3D") ? "OK" : "FAILED")));
      cvcapp.log(5, boost::str(boost::format("gl: Initializing GL_EXT_paletted_texture: %s")%(glewIsSupported("GL_EXT_paletted_texture") ? "OK" : "FAILED")));
      cvcapp.log(5, boost::str(boost::format("gl: Initializing GL_NV_vertex_program: %s")%(glewIsSupported("GL_NV_vertex_program") ? "OK" : "FAILED")));
      cvcapp.log(5, boost::str(boost::format("gl: Initializing GL_NV_fragment_program: %s")%(glewIsSupported("GL_NV_fragment_program") ? "OK" : "FAILED")));
      cvcapp.log(5, boost::str(boost::format("gl: Initializing GL_ARB_multitexture: %s")%(glewIsSupported("GL_ARB_multitexture") ? "OK" : "FAILED")));
      cvcapp.log(5, boost::str(boost::format("gl: Initializing GL_ARB_vertex_program: %s")%(glewIsSupported("GL_ARB_vertex_program") ? "OK" : "FAILED")));
      cvcapp.log(5, boost::str(boost::format("gl: Initializing GL_ARB_fragment_program: %s")%(glewIsSupported("GL_ARB_fragment_program") ? "OK" : "FAILED")));
      cvcapp.log(5, boost::str(boost::format("gl: Initializing GL_SGI_texture_color_table: %s")%(glewIsSupported("GL_SGI_texture_color_table") ? "OK" : "FAILED")));
      cvcapp.log(5, boost::str(boost::format("gl: Initializing GL_SGI_color_table: %s")%(glewIsSupported("GL_SGI_color_table") ? "OK" : "FAILED")));
      cvcapp.log(5, boost::str(boost::format("gl: Initializing GL_NV_register_combiners: %s")%(glewIsSupported("GL_NV_register_combiners") ? "OK" : "FAILED")));
      cvcapp.log(5, boost::str(boost::format("gl: Initializing GL_NV_texture_shader: %s")%(glewIsSupported("GL_NV_texture_shader") ? "OK" : "FAILED")));
      cvcapp.log(5, boost::str(boost::format("gl: Initializing GL_NV_texture_shader2: %s")%(glewIsSupported("GL_NV_texture_shader2") ? "OK" : "FAILED")));
      cvcapp.log(5, boost::str(boost::format("gl: Initializing GL_NV_texture_shader3: %s")%(glewIsSupported("GL_NV_texture_shader3") ? "OK" : "FAILED")));
      cvcapp.log(5, boost::str(boost::format("gl: Initializing GL_ARB_vertex_buffer_object: %s")%(glewIsSupported("GL_ARB_vertex_buffer_object") ? "OK" : "FAILED")));
    }
}

void SliceCanvas::draw()
{
  if(m_SliceDirty)
    updateSlice();
  
  if(m_SliceRenderer) m_SliceRenderer->draw();
}

void SliceCanvas::uploadSlice(unsigned char *slice, unsigned int dimx, unsigned int dimy)
{
  if(m_SliceRenderer) m_SliceRenderer->uploadSlice(slice,dimx,dimy);
}

void SliceCanvas::uploadColorTable()
{
  if(m_SliceRenderer) m_SliceRenderer->uploadColorTable();
}

void SliceCanvas::drawPointClasses()
{
  GLUquadric *q = gluNewQuadric();
  gluQuadricNormals(q,GLU_NONE);
  /* cycle through the point classes and draw all the points */
  if(m_PointClassList) {
    int sizeClass = m_PointClassList[m_Variable][m_Timestep]->size();
    for( int i = 0; i < sizeClass; i ++ )
    {
       PointClass *pc = m_PointClassList[m_Variable][m_Timestep]->at(i);
       int sizePoint = pc->getPointList().size();
       for( int j=0; j<sizePoint; j++)
       {
          GridPoint *gp = pc->getPointList()[j];
	  glColor3f(float(pc->getColor().red())/255.0,
		       float(pc->getColor().green())/255.0,
		       float(pc->getColor().blue())/255.0);
	  
	  switch(m_SliceAxis)
	    {
	    case XY:
	      if(gp->z == m_Depth)
		{
		  glPushMatrix();
		  /* Translate to the right position above the slice quad such that the disk draws over the center of it's voxel */
		  glTranslatef((float(gp->x)/float(m_VolumeFileInfo.XDim()))*(m_VertCoordMaxX-m_VertCoordMinX)+m_VertCoordMinX + (1.0f/float(m_VolumeFileInfo.XDim()))*(m_VertCoordMaxX-m_VertCoordMinX)*(m_VolumeFileInfo.XSpan()/2.0f),
				  (float(gp->y)/float(m_VolumeFileInfo.YDim()))*(m_VertCoordMaxY-m_VertCoordMinY)+m_VertCoordMinY + (1.0f/float(m_VolumeFileInfo.YDim()))*(m_VertCoordMaxY-m_VertCoordMinY)*(m_VolumeFileInfo.YSpan()/2.0f),
				  0.001f); /* position slightly in front of the slice quad */
		  glScalef((m_VertCoordMaxX-m_VertCoordMinX)/m_VolumeFileInfo.XDim(),
			      (m_VertCoordMaxY-m_VertCoordMinY)/m_VolumeFileInfo.YDim(),
			      1.0);
		  gluDisk(q,0.0,m_PointSize,8,1);
		  glPopMatrix();
		}
	      break;
	    case XZ:
	      if(gp->y == m_Depth)
		{
		  glPushMatrix();
		  /* Translate to the right position above the slice quad such that the disk draws over the center of it's voxel */
		  glTranslatef((float(gp->x)/float(m_VolumeFileInfo.XDim()))*(m_VertCoordMaxX-m_VertCoordMinX)+m_VertCoordMinX + (1.0f/float(m_VolumeFileInfo.XDim()))*(m_VertCoordMaxX-m_VertCoordMinX)*(m_VolumeFileInfo.XSpan()/2.0f),
				  (float(gp->z)/float(m_VolumeFileInfo.ZDim()))*(m_VertCoordMaxY-m_VertCoordMinY)+m_VertCoordMinY + (1.0f/float(m_VolumeFileInfo.ZDim()))*(m_VertCoordMaxY-m_VertCoordMinY)*(m_VolumeFileInfo.ZSpan()/2.0f),
				  0.001f); /* position slightly in front of the slice quad */
		  glScalef((m_VertCoordMaxX-m_VertCoordMinX)/m_VolumeFileInfo.XDim(),
			      (m_VertCoordMaxY-m_VertCoordMinY)/m_VolumeFileInfo.ZDim(),
			      1.0);
		  gluDisk(q,0.0,m_PointSize,8,1);
		  glPopMatrix();
		}
	      break;
	    case ZY:
	      if(gp->x == m_Depth)
		{
		  glPushMatrix();
		  /* Translate to the right position above the slice quad such that the disk draws over the center of it's voxel */
		  glTranslatef((float(gp->z)/float(m_VolumeFileInfo.ZDim()))*(m_VertCoordMaxX-m_VertCoordMinX)+m_VertCoordMinX + (1.0f/float(m_VolumeFileInfo.ZDim()))*(m_VertCoordMaxX-m_VertCoordMinX)*(m_VolumeFileInfo.ZSpan()/2.0f),
				  (float(gp->y)/float(m_VolumeFileInfo.YDim()))*(m_VertCoordMaxY-m_VertCoordMinY)+m_VertCoordMinY + (1.0f/float(m_VolumeFileInfo.YDim()))*(m_VertCoordMaxY-m_VertCoordMinY)*(m_VolumeFileInfo.YSpan()/2.0f),
				  0.001f); /* position slightly in front of the slice quad */
		  glScalef((m_VertCoordMaxX-m_VertCoordMinX)/m_VolumeFileInfo.ZDim(),
			      (m_VertCoordMaxY-m_VertCoordMinY)/m_VolumeFileInfo.YDim(),
			      1.0);
		  gluDisk(q,0.0,m_PointSize,8,1);
		  glPopMatrix();
		}
	      break;
	    } //end of switch
	} // end of for
     } // end of for
  } // end of if
  gluDeleteQuadric(q);
}

void SliceCanvas::drawContours(bool withNames)
{
  using namespace std;
  
  GLUquadric *q = gluNewQuadric();
  gluQuadricNormals(q,GLU_NONE);
  
  const gsl_interp_type *interp;//gsl_interp_cspline; /* cubic spline with natural boundary conditions */
  gsl_interp_accel *acc_x, *acc_y, *acc_z;
  gsl_spline *spline_x, *spline_y, *spline_z;
  std::vector<double> vec_x, vec_y, vec_z, vec_t; //used to build the input arrays for gsl_spline_init()
  double t = 0.0;
  double interval = 0.5;
  int contourName = 0;

  m_AllPoints.clear();

  if(!m_Contours.empty())
    {
      for(SurfRecon::ContourPtrMap::iterator i = m_Contours[m_Variable][m_Timestep].begin();
	  i != m_Contours[m_Variable][m_Timestep].end();
	  i++)
	{
	  if((*i).second == NULL) continue;
	  for(std::vector<SurfRecon::CurvePtr>::iterator cur = (*i).second->curves().begin();
	      cur != (*i).second->curves().end();
	      cur++)
	    {
	      //if the curve lies on this slice, draw it
	      if(boost::get<0>(**cur) == m_Depth &&
		 boost::get<1>(**cur) == SurfRecon::Orientation(m_SliceAxis))
		{
		  vec_x.clear(); vec_y.clear(); vec_z.clear(); vec_t.clear(); t = 0.0;
		  
		  /*
		  glColor3f(m_Contours[m_Variable][m_Timestep][i]->color().get<0>(),
			       m_Contours[m_Variable][m_Timestep][i]->color().get<1>(),
			       m_Contours[m_Variable][m_Timestep][i]->color().get<2>());
		  */

		  //render the control points
		  for(SurfRecon::PointPtrList::iterator pcur = boost::get<2>(**cur).begin();
		      pcur != boost::get<2>(**cur).end();
		      pcur++, contourName++)
		    {
		      SurfRecon::Point p;
		      qglviewer::Vec v, screenv, discv;
		      double radius;
		      
		      p = **pcur;
		      
		      //collect the point in the arrays
		      vec_x.push_back(p.x());
		      vec_y.push_back(p.y());
		      vec_z.push_back(p.z());
		      vec_t.push_back(t); t+=1.0;

		      //get the world coordinates for the control point
		      switch(m_SliceAxis)
			{
			case XY:
			  v = qglviewer::Vec(((p.x()-m_VolumeSlice.XMin())/(m_VolumeSlice.XMax()-m_VolumeSlice.XMin()))*(m_VertCoordMaxX-m_VertCoordMinX)+m_VertCoordMinX,
					     ((p.y()-m_VolumeSlice.YMin())/(m_VolumeSlice.YMax()-m_VolumeSlice.YMin()))*(m_VertCoordMaxY-m_VertCoordMinY)+m_VertCoordMinY,
					     0.01f);
			  break;
			case XZ:
			  v = qglviewer::Vec(((p.x()-m_VolumeSlice.XMin())/(m_VolumeSlice.XMax()-m_VolumeSlice.XMin()))*(m_VertCoordMaxX-m_VertCoordMinX)+m_VertCoordMinX,
					     ((p.z()-m_VolumeSlice.ZMin())/(m_VolumeSlice.ZMax()-m_VolumeSlice.ZMin()))*(m_VertCoordMaxY-m_VertCoordMinY)+m_VertCoordMinY,
					     0.01f);
			  break;
			case ZY:
			  v = qglviewer::Vec(((p.z()-m_VolumeSlice.ZMin())/(m_VolumeSlice.ZMax()-m_VolumeSlice.ZMin()))*(m_VertCoordMaxX-m_VertCoordMinX)+m_VertCoordMinX,
					     ((p.y()-m_VolumeSlice.YMin())/(m_VolumeSlice.YMax()-m_VolumeSlice.YMin()))*(m_VertCoordMaxY-m_VertCoordMinY)+m_VertCoordMinY,
					     0.01f);
			  break;
			}

		      //keep the disc radius constant relative to screen coordinates
		      screenv = camera()->projectedCoordinatesOf(v);
		      discv = camera()->unprojectedCoordinatesOf(qglviewer::Vec(screenv.x + m_PointSize+3,
										screenv.y,
										screenv.z));
		      radius = (v - discv).norm();

		      glPushMatrix();
		      if(withNames) 
			{
			  glPushName(contourName);
			  //we must update this list even when not drawing
			  m_AllPoints.push_back(SurfRecon::WeakPointPtr(*pcur));
			}
		      glTranslatef(v.x,v.y,v.z);


		      if(m_SelectedPoints.find(*pcur) != m_SelectedPoints.end())
			glColor3f(1.0-(*i).second->color().get<0>(),
				     1.0-(*i).second->color().get<1>(),
				     1.0-(*i).second->color().get<2>());
		      else
			glColor3f((*i).second->color().get<0>(),
				  (*i).second->color().get<1>(),
				  (*i).second->color().get<2>());

		      if(m_RenderControlPoints) gluDisk(q,0.0,radius,8,1);
		      if(withNames) 
			glPopName();
		      glPopMatrix();
		    }
	
		  //now render the curve (note that it cannot be selected...)
		  if(!withNames)
		    {
		      /* set interpolation type */
		      switch((*i).second->interpolationType())
			{
			case 0: interp = gsl_interp_linear; break;
			case 1: interp = gsl_interp_polynomial; break;
			case 2: interp = gsl_interp_cspline; break;
			case 3: interp = gsl_interp_cspline_periodic; break;
			case 4: interp = gsl_interp_akima; break;
			case 5: interp = gsl_interp_akima_periodic; break;
			  
			default: interp = gsl_interp_cspline; break;
			}

		      gsl_spline *test_spline = gsl_spline_alloc(interp, 1000); //this is hackish but how else can i find the min size?
		      if(vec_t.size() >= gsl_spline_min_size(test_spline))
			{
			  acc_x = gsl_interp_accel_alloc();
			  acc_y = gsl_interp_accel_alloc();
			  acc_z = gsl_interp_accel_alloc();
			  spline_x = gsl_spline_alloc(interp, vec_t.size());
			  spline_y = gsl_spline_alloc(interp, vec_t.size());
			  spline_z = gsl_spline_alloc(interp, vec_t.size());
			  gsl_spline_init(spline_x, &(vec_t[0]), &(vec_x[0]), vec_t.size());
			  gsl_spline_init(spline_y, &(vec_t[0]), &(vec_y[0]), vec_t.size());
			  gsl_spline_init(spline_z, &(vec_t[0]), &(vec_z[0]), vec_t.size());

			  SurfRecon::Point p;
		      
			  glColor4f((*i).second->color().get<0>(),
				       (*i).second->color().get<1>(),
				       (*i).second->color().get<2>(),0.2f);

			  glPushMatrix();
			  glBegin(GL_LINE_STRIP);
			  interval = 1.0/(1+(*i).second->numberOfSamples());
			  for(double time_i = 0.0; time_i<=t-1.0; time_i+=interval)
			    {
			      p = SurfRecon::Point(gsl_spline_eval(spline_x, time_i, acc_x),
						   gsl_spline_eval(spline_y, time_i, acc_y),
						   gsl_spline_eval(spline_z, time_i, acc_z));
			      switch(m_SliceAxis)
				{
				case XY:
				  {
				    glVertex3f(((p.x()-m_VolumeSlice.XMin())/(m_VolumeSlice.XMax()-m_VolumeSlice.XMin()))*(m_VertCoordMaxX-m_VertCoordMinX)+m_VertCoordMinX,
						  ((p.y()-m_VolumeSlice.YMin())/(m_VolumeSlice.YMax()-m_VolumeSlice.YMin()))*(m_VertCoordMaxY-m_VertCoordMinY)+m_VertCoordMinY,
						  0.01f);
				  }
				  break;
				case XZ:
				  {
				    glVertex3f(((p.x()-m_VolumeSlice.XMin())/(m_VolumeSlice.XMax()-m_VolumeSlice.XMin()))*(m_VertCoordMaxX-m_VertCoordMinX)+m_VertCoordMinX,
						  ((p.z()-m_VolumeSlice.ZMin())/(m_VolumeSlice.ZMax()-m_VolumeSlice.ZMin()))*(m_VertCoordMaxY-m_VertCoordMinY)+m_VertCoordMinY,
						  0.01f);
				  }
				  break;
				case ZY:
				  {
				    glVertex3f(((p.z()-m_VolumeSlice.ZMin())/(m_VolumeSlice.ZMax()-m_VolumeSlice.ZMin()))*(m_VertCoordMaxX-m_VertCoordMinX)+m_VertCoordMinX,
						  ((p.y()-m_VolumeSlice.YMin())/(m_VolumeSlice.YMax()-m_VolumeSlice.YMin()))*(m_VertCoordMaxY-m_VertCoordMinY)+m_VertCoordMinY,
						  0.01f);
				  }
				  break;
				}
			    }

			  //do the last point..
			  p = SurfRecon::Point(gsl_spline_eval(spline_x, t-1.0, acc_x),
					       gsl_spline_eval(spline_y, t-1.0, acc_y),
					       gsl_spline_eval(spline_z, t-1.0, acc_z));
			  switch(m_SliceAxis)
			    {
			    case XY:
			      {
				glVertex3f(((p.x()-m_VolumeSlice.XMin())/(m_VolumeSlice.XMax()-m_VolumeSlice.XMin()))*(m_VertCoordMaxX-m_VertCoordMinX)+m_VertCoordMinX,
					      ((p.y()-m_VolumeSlice.YMin())/(m_VolumeSlice.YMax()-m_VolumeSlice.YMin()))*(m_VertCoordMaxY-m_VertCoordMinY)+m_VertCoordMinY,
					      0.01f);
			      }
			      break;
			    case XZ:
			      {
				glVertex3f(((p.x()-m_VolumeSlice.XMin())/(m_VolumeSlice.XMax()-m_VolumeSlice.XMin()))*(m_VertCoordMaxX-m_VertCoordMinX)+m_VertCoordMinX,
					      ((p.z()-m_VolumeSlice.ZMin())/(m_VolumeSlice.ZMax()-m_VolumeSlice.ZMin()))*(m_VertCoordMaxY-m_VertCoordMinY)+m_VertCoordMinY,
					      0.01f);
			      }
			      break;
			    case ZY:
			      {
				glVertex3f(((p.z()-m_VolumeSlice.ZMin())/(m_VolumeSlice.ZMax()-m_VolumeSlice.ZMin()))*(m_VertCoordMaxX-m_VertCoordMinX)+m_VertCoordMinX,
					      ((p.y()-m_VolumeSlice.YMin())/(m_VolumeSlice.YMax()-m_VolumeSlice.YMin()))*(m_VertCoordMaxY-m_VertCoordMinY)+m_VertCoordMinY,
					      0.01f);
			      }
			      break;
			    }
			  
			  glEnd();

			  //now draw a blended polygon for those curves who's contours have been selected and are closed...
			  //using the stencil buffer to draw convex polygons
			  if(i->second->selected() && 
			     (vec_x.size()>3) &&
			     (vec_x.front()==vec_x.back()) && 
			     (vec_y.front()==vec_y.back()) &&
			     (vec_z.front()==vec_z.back()))
			    {
			      //get all the interpolated points
			      vector<SurfRecon::Point> interp_points;
			      vector<double> trans_points;
			      double trans_points_minx, trans_points_maxx,
				trans_points_miny, trans_points_maxy,
				trans_points_centerx, trans_points_centery;

			      if((*i).second->interpolationType()==0)
				{
				  for(vector<double>::iterator vec_x_itr = vec_x.begin(),
					vec_y_itr = vec_y.begin(),
					vec_z_itr = vec_z.begin();
				      vec_x_itr != vec_x.end();
				      vec_x_itr++,vec_y_itr++,vec_z_itr++)
				    interp_points.push_back(SurfRecon::Point(*vec_x_itr,
									     *vec_y_itr,
									     *vec_z_itr));
				}
			      else
				{
				  interval = 1.0/(1+(*i).second->numberOfSamples());
				  for(double time_i = 0.0; time_i<=t-1.0; time_i+=interval)
				    interp_points.push_back(SurfRecon::Point(gsl_spline_eval(spline_x, time_i, acc_x),
									     gsl_spline_eval(spline_y, time_i, acc_y),
									     gsl_spline_eval(spline_z, time_i, acc_z)));
				}
			      
			      //transform the points to the world coordinate system...
			      switch(m_SliceAxis)
				{
				case XY:     
				  for(vector<SurfRecon::Point>::iterator interp = interp_points.begin();
				      interp != interp_points.end();
				      interp++)
				    {
				      trans_points.push_back(((interp->x()-m_VolumeSlice.XMin())/(m_VolumeSlice.XMax()-m_VolumeSlice.XMin()))*(m_VertCoordMaxX-m_VertCoordMinX)+m_VertCoordMinX);
				      trans_points.push_back(((interp->y()-m_VolumeSlice.YMin())/(m_VolumeSlice.YMax()-m_VolumeSlice.YMin()))*(m_VertCoordMaxY-m_VertCoordMinY)+m_VertCoordMinY);
				      trans_points.push_back(0.01f);
				    }
				  break;
				case XZ:
				  for(vector<SurfRecon::Point>::iterator interp = interp_points.begin();
				      interp != interp_points.end();
				      interp++)
				    {
				      trans_points.push_back(((interp->x()-m_VolumeSlice.XMin())/(m_VolumeSlice.XMax()-m_VolumeSlice.XMin()))*(m_VertCoordMaxX-m_VertCoordMinX)+m_VertCoordMinX);
				      trans_points.push_back(((interp->z()-m_VolumeSlice.ZMin())/(m_VolumeSlice.ZMax()-m_VolumeSlice.ZMin()))*(m_VertCoordMaxY-m_VertCoordMinY)+m_VertCoordMinY);
				      trans_points.push_back(0.01f);
				    }
				  break;
				case ZY:
				  for(vector<SurfRecon::Point>::iterator interp = interp_points.begin();
				      interp != interp_points.end();
				      interp++)
				    {
				      trans_points.push_back(((interp->z()-m_VolumeSlice.ZMin())/(m_VolumeSlice.ZMax()-m_VolumeSlice.ZMin()))*(m_VertCoordMaxX-m_VertCoordMinX)+m_VertCoordMinX);
				      trans_points.push_back(((interp->y()-m_VolumeSlice.YMin())/(m_VolumeSlice.YMax()-m_VolumeSlice.YMin()))*(m_VertCoordMaxY-m_VertCoordMinY)+m_VertCoordMinY);
				      trans_points.push_back(0.01f);
				    }
				  break;
				}

			      //get the bounding box of this polygon
			      trans_points_minx = trans_points_maxx = trans_points[0];
			      trans_points_miny = trans_points_maxy = trans_points[1];
			      for(vector<double>::iterator trans_points_itr = trans_points.begin();
				  trans_points_itr != trans_points.end();
				  trans_points_itr += 3)
				{
				  if(*(trans_points_itr+0) < trans_points_minx) trans_points_minx = *(trans_points_itr+0);
				  if(*(trans_points_itr+0) > trans_points_maxx) trans_points_maxx = *(trans_points_itr+0);
				  if(*(trans_points_itr+1) < trans_points_miny) trans_points_miny = *(trans_points_itr+1);
				  if(*(trans_points_itr+1) > trans_points_maxy) trans_points_maxy = *(trans_points_itr+1);
				}

			      trans_points_centerx = trans_points_minx + (trans_points_maxx - trans_points_minx)/2;
			      trans_points_centery = trans_points_miny + (trans_points_maxy - trans_points_miny)/2;

			      glEnable(GL_STENCIL_TEST);
			      glClearStencil(0x0);
			      glClear(GL_STENCIL_BUFFER_BIT); //clear the stencil buffer for this polygon

			      glColorMask(GL_FALSE,GL_FALSE,GL_FALSE,GL_FALSE); //disable framebuffer writes
			      glDepthMask(GL_FALSE); //disable depth buffer writes
			      			      			      
			      //check this out - http://www.glprogramming.com/red/chapter14.html#name13
			      glStencilFunc(GL_ALWAYS,0x1,0x1);
			      glStencilOp(GL_INVERT,GL_INVERT,GL_INVERT);
			      
			      glBegin(GL_TRIANGLE_FAN);
			      glVertex3f(trans_points_centerx,
					    trans_points_centery,
					    0.01f);
			      for(vector<double>::iterator trans_points_itr = trans_points.begin();
				  trans_points_itr != (trans_points.end()); //the last point duplicates the first
				  trans_points_itr += 3)
				glVertex3f(*(trans_points_itr+0),
					      *(trans_points_itr+1),
					      *(trans_points_itr+2));
			      glEnd();

			      glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE); //re-enable framebuffer writes
			      glDepthMask(GL_TRUE); //re-enable depth buffer writes
			      
			      glEnable(GL_BLEND);
			      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			      //draw the polygon's bounding box... the magical stencil buffer will only show our
			      //concave polygon
			      glStencilFunc(GL_EQUAL,0x1,0x1);
			      glStencilOp(GL_KEEP,GL_KEEP,GL_KEEP);
			      glBegin(GL_QUADS);
			      glVertex3f(trans_points_minx, trans_points_miny, 0.01f);
			      glVertex3f(trans_points_maxx, trans_points_miny, 0.01f);
			      glVertex3f(trans_points_maxx, trans_points_maxy, 0.01f);
			      glVertex3f(trans_points_minx, trans_points_maxy, 0.01f);
			      glEnd();

			      glDisable(GL_BLEND);
			      glDisable(GL_STENCIL_TEST);
			    }
			  
			  glPopMatrix();

			  /* clean up our mess */
			  gsl_spline_free(spline_x);
			  gsl_spline_free(spline_y);
			  gsl_spline_free(spline_z);
			  gsl_interp_accel_free(acc_x);
			  gsl_interp_accel_free(acc_y);
			  gsl_interp_accel_free(acc_z);
			}
		      gsl_spline_free(test_spline);
		    }
		}
	    }
	}
    }

  gluDeleteQuadric(q);
}

void SliceCanvas::drawSelectionRectangle()
{
  startScreenCoordinatesSystem();
  
  glPushAttrib(GL_ALL_ATTRIB_BITS);

  glDisable(GL_LIGHTING);
  glEnable(GL_BLEND);
  glBlendFunc(GL_ONE, GL_ONE);

  glColor4f(0.0, 0.0, 0.3, 0.3);
  glBegin(GL_QUADS);
  glVertex2i(m_SelectionRectangle.left(),  m_SelectionRectangle.top());
  glVertex2i(m_SelectionRectangle.right(), m_SelectionRectangle.top());
  glVertex2i(m_SelectionRectangle.right(), m_SelectionRectangle.bottom());
  glVertex2i(m_SelectionRectangle.left(),  m_SelectionRectangle.bottom());
  glEnd();

  glLineWidth(2.0);
  glColor4f(0.4, 0.4, 0.5, 0.5);
  glBegin(GL_LINE_LOOP);
  glVertex2i(m_SelectionRectangle.left(),  m_SelectionRectangle.top());
  glVertex2i(m_SelectionRectangle.right(), m_SelectionRectangle.top());
  glVertex2i(m_SelectionRectangle.right(), m_SelectionRectangle.bottom());
  glVertex2i(m_SelectionRectangle.left(),  m_SelectionRectangle.bottom());
  glEnd();

  glDisable(GL_BLEND);
  glEnable(GL_LIGHTING);

  glPopAttrib();

  stopScreenCoordinatesSystem();
}

void SliceCanvas::postDraw()
{
  drawPointClasses();
  drawContours();
  if(m_SelectionMode != 0)
    drawSelectionRectangle();
  QGLViewer::postDraw();
}

void SliceCanvas::drawWithNames()
{
  drawContours(true);
}

void SliceCanvas::endSelection(const QPoint& p)
{
  Q_UNUSED(p);

  glFlush();

  // Get the number of objects that were seen through the pick matrix frustum. Reset GL_RENDER mode.
  GLint nbHits = glRenderMode(GL_RENDER);

  if(m_SelectionMode == 1) m_SelectedPoints.clear(); //new selection, so clear old selected points

  if(nbHits > 0)
    {
      // Interpret results : each object created 4 values in the selectBuffer().
      // (selectBuffer())[4*i+3] is the id pushed on the stack.
      for(int i=0; i<nbHits; ++i)
	switch (m_SelectionMode)
	  {
	  case 1: 
	  case 2: addIdToSelection((selectBuffer())[4*i+3]); break;
	  case 3: removeIdFromSelection((selectBuffer())[4*i+3]); break;
	  default : break;
	  }
    }

  m_SelectionMode = 0;
}

void SliceCanvas::addIdToSelection(int id)
{
  //cvcapp.log(5, "adding id: %d",id);

  m_SelectedPoints.insert(m_AllPoints[id]);
}

void SliceCanvas::removeIdFromSelection(int id)
{
  //cvcapp.log(5, "removing id: %d", id);

  m_SelectedPoints.erase(m_AllPoints[id]);
}

void SliceCanvas::deleteSelected()
{
  if(!m_Contours.empty() && !m_SelectedPoints.empty())
    for(SurfRecon::ContourPtrMap::iterator i = m_Contours[m_Variable][m_Timestep].begin();
	i != m_Contours[m_Variable][m_Timestep].end();
	i++)
      {
	if((*i).second == NULL) continue;
	for(std::vector<SurfRecon::CurvePtr>::iterator cur = (*i).second->curves().begin();
	    cur != (*i).second->curves().end();
	    cur++)
	  {
	    if(SurfRecon::getCurvePoints(**cur).size() == 0) continue;

	    for(SurfRecon::PointPtrList::iterator pcur = boost::get<2>(**cur).begin();
		pcur != boost::get<2>(**cur).end();
		pcur++)
	      if(m_SelectedPoints.find(SurfRecon::WeakPointPtr(*pcur)) != m_SelectedPoints.end())
		{
		  m_SelectedPoints.erase(SurfRecon::WeakPointPtr(*pcur));
		  boost::get<2>(**cur).erase(pcur); //erase the current point from the current curve
		}
	  }
      }
}

/*********** ARB Fragment Program Slice Renderer ************/
ARBFragmentProgramSliceRenderer::ARBFragmentProgramSliceRenderer(SliceCanvas *canvas)
  : SliceRenderer(canvas), m_FragmentProgram(0), m_PaletteTexture(0) {}
 
ARBFragmentProgramSliceRenderer::~ARBFragmentProgramSliceRenderer() {}

bool ARBFragmentProgramSliceRenderer::init()
{
  /* make sure the proper extensions are initialized */
  fprintf(stderr, "ARBFragmentProgramSliceRenderer::init(): Checking extensions - ARB_vertex_program=%d ARB_fragment_program=%d ARB_multitexture=%d\n",
    glewIsSupported("GL_ARB_vertex_program"),
    glewIsSupported("GL_ARB_fragment_program"), 
    glewIsSupported("GL_ARB_multitexture"));
  
  if(//glewIsSupported("GL_VERSION_1_3") &&
     glewIsSupported("GL_ARB_vertex_program") &&
     glewIsSupported("GL_ARB_fragment_program") &&
     glewIsSupported("GL_ARB_multitexture"))
    {
      /* Initialize the fragment program */
      const GLubyte program[] = 
        "!!ARBfp1.0\n"
        "PARAM c0 = {0.5, 1, 2.7182817, 0};\n"
        "TEMP R0;\n"
        "TEX R0.x, fragment.texcoord[0].xyzx, texture[0], 2D;\n"
        "TEX result.color, R0.x, texture[1], 1D;\n"
        "END\n";
      /* initialize the fragment program */
      glEnable(GL_FRAGMENT_PROGRAM_ARB);
      glGenProgramsARB(1,&(m_FragmentProgram));
      glBindProgramARB(GL_FRAGMENT_PROGRAM_ARB, m_FragmentProgram);
      glProgramStringARB(GL_FRAGMENT_PROGRAM_ARB, GL_PROGRAM_FORMAT_ASCII_ARB, strlen((const char *)program), program);
      glDisable(GL_FRAGMENT_PROGRAM_ARB);
    }
  else
    {
      cvcapp.log(5, "ARBFragmentProgramSliceRenderer::init(): Unable to initialize renderer");
      return false;
    }

  return true;
}

void ARBFragmentProgramSliceRenderer::draw()
{
  if(!m_SliceCanvas->m_Drawable || !m_SliceCanvas->m_VolumeFileInfo.isSet()) return;

  fprintf(stderr, "ARBFragmentProgramSliceRenderer::draw(): NumSliceTiles=%u, Drawable=%d, VolumeFileInfo.isSet()=%d\n",
    m_SliceCanvas->m_NumSliceTiles, m_SliceCanvas->m_Drawable, m_SliceCanvas->m_VolumeFileInfo.isSet());

  glColor3f(1.0,1.0,1.0);
  
  glPushAttrib(GL_ENABLE_BIT);
  
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  
  glEnable(GL_FRAGMENT_PROGRAM_ARB);
  glBindProgramARB(GL_FRAGMENT_PROGRAM_ARB, m_FragmentProgram);
  
  static int debugCount = 0;
  if(debugCount < 3) {  // Only print first 3 frames to avoid spam
    fprintf(stderr, "ARBFragmentProgramSliceRenderer::draw(): FragmentProgram=%u, PaletteTexture=%u\n",
      m_FragmentProgram, m_PaletteTexture);
    debugCount++;
  }
  
  for(unsigned int i=0; i<m_SliceCanvas->m_NumSliceTiles; i++)
    {
      /* bind the transfer function */
      glActiveTextureARB(GL_TEXTURE1_ARB);
      glEnable(GL_TEXTURE_1D);
      glBindTexture(GL_TEXTURE_1D, m_PaletteTexture);
	
      /* bind the data texture */
      glActiveTextureARB(GL_TEXTURE0_ARB);
      glEnable(GL_TEXTURE_2D);
      glBindTexture(GL_TEXTURE_2D, m_SliceCanvas->m_SliceTiles[i].texture);
      
      if(debugCount <= 3) {
        fprintf(stderr, "  Tile[%u]: texture=%u, tex_coords=[%.3f,%.3f %.3f,%.3f %.3f,%.3f %.3f,%.3f], vert_coords=[%.3f,%.3f %.3f,%.3f %.3f,%.3f %.3f,%.3f]\n",
          i, m_SliceCanvas->m_SliceTiles[i].texture,
          m_SliceCanvas->m_SliceTiles[i].tex_coords[0], m_SliceCanvas->m_SliceTiles[i].tex_coords[1],
          m_SliceCanvas->m_SliceTiles[i].tex_coords[2], m_SliceCanvas->m_SliceTiles[i].tex_coords[3],
          m_SliceCanvas->m_SliceTiles[i].tex_coords[4], m_SliceCanvas->m_SliceTiles[i].tex_coords[5],
          m_SliceCanvas->m_SliceTiles[i].tex_coords[6], m_SliceCanvas->m_SliceTiles[i].tex_coords[7],
          m_SliceCanvas->m_SliceTiles[i].vert_coords[0], m_SliceCanvas->m_SliceTiles[i].vert_coords[1],
          m_SliceCanvas->m_SliceTiles[i].vert_coords[2], m_SliceCanvas->m_SliceTiles[i].vert_coords[3],
          m_SliceCanvas->m_SliceTiles[i].vert_coords[4], m_SliceCanvas->m_SliceTiles[i].vert_coords[5],
          m_SliceCanvas->m_SliceTiles[i].vert_coords[6], m_SliceCanvas->m_SliceTiles[i].vert_coords[7]);
      }
	  
      glBegin(GL_QUADS);
      for(unsigned int j=0; j<4; j++)
	{
	  //glTexCoord2f(m_SliceTiles[i].tex_coords[j*2+0],m_SliceTiles[i].tex_coords[j*2+1]);
	  glMultiTexCoord2fARB(GL_TEXTURE0_ARB,m_SliceCanvas->m_SliceTiles[i].tex_coords[j*2+0],m_SliceCanvas->m_SliceTiles[i].tex_coords[j*2+1]);
	  glVertex2f(m_SliceCanvas->m_SliceTiles[i].vert_coords[j*2+0],m_SliceCanvas->m_SliceTiles[i].vert_coords[j*2+1]);
	}
      glEnd();
    }
	
  glBindProgramARB(GL_FRAGMENT_PROGRAM_ARB, 0);
  glDisable(GL_FRAGMENT_PROGRAM_ARB);
 
  /* now draw the point class points as discs */
  glDisable(GL_TEXTURE_1D);
  glDisable(GL_TEXTURE_2D);
  glDisable(GL_BLEND);
   
  glPopAttrib();
}

void ARBFragmentProgramSliceRenderer::uploadSlice(unsigned char *slice, unsigned int dimx, unsigned int dimy)
{
  unsigned int i;
  //  unsigned int num_tile_x, num_tile_y,a,b,j,index;
  //  unsigned char tex_buf[TEXTURE_TILE_X*TEXTURE_TILE_Y];
  //  GLuint *texture_ids;
  /* 
     float tex_square[8] = { 0.0, 1.0, 1.0, 1.0,        
     1.0, 0.0, 0.0, 0.0 };      
     float vert_square[8] = { -0.5, 0.5, 0.5, 0.5,      
     0.5, -0.5, -0.5, -0.5 };  
			   
  */
  fprintf(stderr, "ARBFragmentProgramSliceRenderer::uploadSlice(): dimx=%u, dimy=%u, slice=%p\n", dimx, dimy, slice);
  
  if(!m_SliceCanvas->m_Drawable) return;
  
  m_SliceCanvas->makeCurrent();
  

#ifndef MULTI_TILE
  if(m_SliceCanvas->m_SliceTiles)
    {
      for(i=0; i<m_SliceCanvas->m_NumSliceTiles; i++) glDeleteTextures(1,&(m_SliceCanvas->m_SliceTiles[i]).texture);
      delete [] m_SliceCanvas->m_SliceTiles;
    }
  m_SliceCanvas->m_NumSliceTiles = 1;
  m_SliceCanvas->m_SliceTiles = new SliceCanvas::SliceTile[m_SliceCanvas->m_NumSliceTiles];
  
  glGenTextures(1,&(m_SliceCanvas->m_SliceTiles[0].texture));
  
  m_SliceCanvas->m_SliceTiles[0].tex_coords[0] = m_SliceCanvas->m_TexCoordMinX; m_SliceCanvas->m_SliceTiles[0].tex_coords[1] = m_SliceCanvas->m_TexCoordMaxY;      /* top left */
  m_SliceCanvas->m_SliceTiles[0].tex_coords[2] = m_SliceCanvas->m_TexCoordMaxX; m_SliceCanvas->m_SliceTiles[0].tex_coords[3] = m_SliceCanvas->m_TexCoordMaxY;      /* top right */
  m_SliceCanvas->m_SliceTiles[0].tex_coords[4] = m_SliceCanvas->m_TexCoordMaxX; m_SliceCanvas->m_SliceTiles[0].tex_coords[5] = m_SliceCanvas->m_TexCoordMinY;      /* bottom right */
  m_SliceCanvas->m_SliceTiles[0].tex_coords[6] = m_SliceCanvas->m_TexCoordMinX; m_SliceCanvas->m_SliceTiles[0].tex_coords[7] = m_SliceCanvas->m_TexCoordMinY;      /* bottom left */
  
  m_SliceCanvas->m_SliceTiles[0].vert_coords[0] = m_SliceCanvas->m_VertCoordMinX; m_SliceCanvas->m_SliceTiles[0].vert_coords[1] = m_SliceCanvas->m_VertCoordMaxY;  /* top left */
  m_SliceCanvas->m_SliceTiles[0].vert_coords[2] = m_SliceCanvas->m_VertCoordMaxX; m_SliceCanvas->m_SliceTiles[0].vert_coords[3] = m_SliceCanvas->m_VertCoordMaxY;  /* top right */
  m_SliceCanvas->m_SliceTiles[0].vert_coords[4] = m_SliceCanvas->m_VertCoordMaxX; m_SliceCanvas->m_SliceTiles[0].vert_coords[5] = m_SliceCanvas->m_VertCoordMinY;  /* bottom right */
  m_SliceCanvas->m_SliceTiles[0].vert_coords[6] = m_SliceCanvas->m_VertCoordMinX; m_SliceCanvas->m_SliceTiles[0].vert_coords[7] = m_SliceCanvas->m_VertCoordMinY;  /* bottom left */
  
  glBindTexture(GL_TEXTURE_2D, m_SliceCanvas->m_SliceTiles[0].texture);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_R8, dimx, dimy, 0, GL_RED, GL_UNSIGNED_BYTE, slice);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  // Set swizzle to replicate red channel to all components for compatibility
  GLint swizzleMask[] = {GL_RED, GL_RED, GL_RED, GL_ONE};
  glTexParameteriv(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_RGBA, swizzleMask);
  
  uploadColorTable();
#else
  
  //dimx and dimy must be 2^n so that makes things easier!

  const unsigned int max_dimx = 256;
  const unsigned int max_dimy = 256;
  
  unsigned int num_x_splits = max_dimx > dimx ? 1 : dimx/max_dimx;
  unsigned int num_y_splits = max_dimy > dimy ? 1 : dimy/max_dimy;

  float vec_x_inc = (m_SliceCanvas->m_VertCoordMaxX - m_SliceCanvas->m_VertCoordMinX)/num_x_splits;
  float vec_y_inc = (m_SliceCanvas->m_VertCoordMaxY - m_SliceCanvas->m_VertCoordMinY)/num_y_splits;

  if(m_SliceCanvas->m_SliceTiles)
    {
      for(i=0; i<m_SliceCanvas->m_NumSliceTiles; i++) glDeleteTextures(1,&(m_SliceCanvas->m_SliceTiles[i]).texture);
      delete [] m_SliceCanvas->m_SliceTiles;
    }
  m_SliceCanvas->m_NumSliceTiles = num_x_splits*num_y_splits;
  m_SliceCanvas->m_SliceTiles = new SliceCanvas::SliceTile[m_SliceCanvas->m_NumSliceTiles];

  std::vector<GLuint> tex_ids(m_SliceCanvas->m_NumSliceTiles);
  glGenTextures(tex_ids.size(),&(tex_ids[0]));
  for(std::vector<GLuint>::iterator id = tex_ids.begin();
      id != tex_ids.end();
      id++)
    m_SliceCanvas->m_SliceTiles[std::distance(tex_ids.begin(),id)].texture = *id;
  
  for(i=0; i < m_SliceCanvas->m_NumSliceTiles)
    {
      m_SliceCanvas->m_SliceTiles[i].tex_coords[0] = m_SliceCanvas->m_TexCoordMinX; m_SliceCanvas->m_SliceTiles[i].tex_coords[1] = m_SliceCanvas->m_TexCoordMaxY;      /* top left */
      m_SliceCanvas->m_SliceTiles[i].tex_coords[2] = m_SliceCanvas->m_TexCoordMaxX; m_SliceCanvas->m_SliceTiles[i].tex_coords[3] = m_SliceCanvas->m_TexCoordMaxY;      /* top right */
      m_SliceCanvas->m_SliceTiles[i].tex_coords[4] = m_SliceCanvas->m_TexCoordMaxX; m_SliceCanvas->m_SliceTiles[i].tex_coords[5] = m_SliceCanvas->m_TexCoordMinY;      /* bottom right */
      m_SliceCanvas->m_SliceTiles[i].tex_coords[6] = m_SliceCanvas->m_TexCoordMinX; m_SliceCanvas->m_SliceTiles[i].tex_coords[7] = m_SliceCanvas->m_TexCoordMinY;      /* bottom left */

      for(unsigned int x_split = 0; x_split < num_x_splits; x_split++)
	for(unsigned int y_split = 0; y_split < num_y_splits; y_split++)
	  {
	    float minx = m_SliceCanvas->m_VertCoordMinX + x_split * vec_x_inc;
	    float maxx = m_SliceCanvas->m_VertCoordMinX + (x_split+1) * vec_x_inc;
	    float miny = m_SliceCanvas->m_VertCoordMinY + y_split * vec_y_inc;
	    float maxy = m_SliceCanvas->m_VertCoordMinY + (y_split+1) * vec_y_inc;

	    m_SliceCanvas->m_SliceTiles[0].vert_coords[0] = minx; m_SliceCanvas->m_SliceTiles[0].vert_coords[1] = maxy;  /* top left */
	    m_SliceCanvas->m_SliceTiles[0].vert_coords[2] = maxx; m_SliceCanvas->m_SliceTiles[0].vert_coords[3] = maxy;  /* top right */
	    m_SliceCanvas->m_SliceTiles[0].vert_coords[4] = maxx; m_SliceCanvas->m_SliceTiles[0].vert_coords[5] = miny;  /* bottom right */
	    m_SliceCanvas->m_SliceTiles[0].vert_coords[6] = minx; m_SliceCanvas->m_SliceTiles[0].vert_coords[7] = miny;  /* bottom left */
	  }

      glBindTexture(GL_TEXTURE_2D, m_SliceCanvas->m_SliceTiles[i].texture);

      //extract a tile from the slice to upload, then upload it!
      unsigned char sub_dimx = (x_split+1)*max_dimx - x_split*max_dimx;

      glTexImage2D(GL_TEXTURE_2D, 0, GL_R8, dimx, dimy, 0, GL_RED, GL_UNSIGNED_BYTE, slice);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      // Set swizzle to replicate red channel to all components for compatibility
      GLint swizzleMask[] = {GL_RED, GL_RED, GL_RED, GL_ONE};
      glTexParameteriv(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_RGBA, swizzleMask);
    }

#endif
}

void ARBFragmentProgramSliceRenderer::uploadColorTable()
{
  if(!m_SliceCanvas->m_Drawable) return;
  if(glIsTexture(m_PaletteTexture))
    glDeleteTextures(1,&(m_PaletteTexture));
  glGenTextures(1,&(m_PaletteTexture));
  glBindTexture(GL_TEXTURE_1D, m_PaletteTexture);
  glTexImage1D(GL_TEXTURE_1D, 0, 4, 256, 0, GL_RGBA, GL_UNSIGNED_BYTE, m_SliceCanvas->m_Palette);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
}

/*************** The Paletted renderer ***************/
PalettedSliceRenderer::PalettedSliceRenderer(SliceCanvas *canvas)
  : SliceRenderer(canvas) {}

PalettedSliceRenderer::~PalettedSliceRenderer() {}

bool PalettedSliceRenderer::init()
{
  /* make sure the proper extensions are initialized */
  if(glewIsSupported("GL_SGIS_texture_edge_clamp") &&
     glewIsSupported("GL_EXT_paletted_texture") &&
     glewIsSupported("GL_ARB_multitexture"))
    {
      return true;
    }
  else
    {
      //QMessageBox::critical(this,"Error","Cannot render slices using hardware!",QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
      cvcapp.log(5, "PalettedSliceRenderer::init(): Unable to initialize renderer");
      return false;
    }
}

void PalettedSliceRenderer::draw()
{
  if(!m_SliceCanvas->m_Drawable || !m_SliceCanvas->m_VolumeFileInfo.isSet()) return;

  glColor3f(1.0,1.0,1.0);
  
  glPushAttrib(GL_ENABLE_BIT);
  
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  
  glEnable(GL_COLOR_TABLE);
  glEnable(GL_TEXTURE_2D);
  
  for(unsigned int i=0; i<m_SliceCanvas->m_NumSliceTiles; i++)
    {
      glBindTexture(GL_TEXTURE_2D, m_SliceCanvas->m_SliceTiles[i].texture);
	  
      glBegin(GL_QUADS);
      for(unsigned int j=0; j<4; j++)
	{
	  glTexCoord2f(m_SliceCanvas->m_SliceTiles[i].tex_coords[j*2+0],m_SliceCanvas->m_SliceTiles[i].tex_coords[j*2+1]);
	  glVertex2f(m_SliceCanvas->m_SliceTiles[i].vert_coords[j*2+0],m_SliceCanvas->m_SliceTiles[i].vert_coords[j*2+1]);
	}
      glEnd();
    }
 
  /* now draw the point class points as discs */
  glDisable(GL_TEXTURE_2D);
  glDisable(GL_COLOR_TABLE);
  glDisable(GL_BLEND);
   
  glPopAttrib();
}

void PalettedSliceRenderer::uploadSlice(unsigned char *slice, unsigned int dimx, unsigned int dimy)
{
  //  unsigned int num_tile_x, num_tile_y,j,a,b,index;
  unsigned int i;
  //  unsigned char tex_buf[TEXTURE_TILE_X*TEXTURE_TILE_Y];
  //  GLuint *texture_ids;
  /*
    float tex_square[8] = { 0.0, 1.0, 1.0, 1.0,        
    1.0, 0.0, 0.0, 0.0 };      
    float vert_square[8] = { -0.5, 0.5, 0.5, 0.5,      
    0.5, -0.5, -0.5, -0.5 };  
			   
  */
  if(!m_SliceCanvas->m_Drawable) return;
  
  m_SliceCanvas->makeCurrent();
  
  m_SliceCanvas->m_NumSliceTiles = 1;
  if(m_SliceCanvas->m_SliceTiles)
    {
      for(i=0; i<m_SliceCanvas->m_NumSliceTiles; i++) glDeleteTextures(1,&(m_SliceCanvas->m_SliceTiles[i]).texture);
      delete [] m_SliceCanvas->m_SliceTiles;
    }
  m_SliceCanvas->m_SliceTiles = new SliceCanvas::SliceTile[m_SliceCanvas->m_NumSliceTiles];
  
  glGenTextures(1,&(m_SliceCanvas->m_SliceTiles[0].texture));
  
  m_SliceCanvas->m_SliceTiles[0].tex_coords[0] = m_SliceCanvas->m_TexCoordMinX; m_SliceCanvas->m_SliceTiles[0].tex_coords[1] = m_SliceCanvas->m_TexCoordMaxY;      /* top left */
  m_SliceCanvas->m_SliceTiles[0].tex_coords[2] = m_SliceCanvas->m_TexCoordMaxX; m_SliceCanvas->m_SliceTiles[0].tex_coords[3] = m_SliceCanvas->m_TexCoordMaxY;      /* top right */
  m_SliceCanvas->m_SliceTiles[0].tex_coords[4] = m_SliceCanvas->m_TexCoordMaxX; m_SliceCanvas->m_SliceTiles[0].tex_coords[5] = m_SliceCanvas->m_TexCoordMinY;      /* bottom right */
  m_SliceCanvas->m_SliceTiles[0].tex_coords[6] = m_SliceCanvas->m_TexCoordMinX; m_SliceCanvas->m_SliceTiles[0].tex_coords[7] = m_SliceCanvas->m_TexCoordMinY;      /* bottom left */
  
  m_SliceCanvas->m_SliceTiles[0].vert_coords[0] = m_SliceCanvas->m_VertCoordMinX; m_SliceCanvas->m_SliceTiles[0].vert_coords[1] = m_SliceCanvas->m_VertCoordMaxY;  /* top left */
  m_SliceCanvas->m_SliceTiles[0].vert_coords[2] = m_SliceCanvas->m_VertCoordMaxX; m_SliceCanvas->m_SliceTiles[0].vert_coords[3] = m_SliceCanvas->m_VertCoordMaxY;  /* top right */
  m_SliceCanvas->m_SliceTiles[0].vert_coords[4] = m_SliceCanvas->m_VertCoordMaxX; m_SliceCanvas->m_SliceTiles[0].vert_coords[5] = m_SliceCanvas->m_VertCoordMinY;  /* bottom right */
  m_SliceCanvas->m_SliceTiles[0].vert_coords[6] = m_SliceCanvas->m_VertCoordMinX; m_SliceCanvas->m_SliceTiles[0].vert_coords[7] = m_SliceCanvas->m_VertCoordMinY;  /* bottom left */
  
  glBindTexture(GL_TEXTURE_2D, m_SliceCanvas->m_SliceTiles[0].texture);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_COLOR_INDEX8_EXT, dimx, dimy, 0, GL_COLOR_INDEX, GL_UNSIGNED_BYTE, slice);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

  uploadColorTable();
}

void PalettedSliceRenderer::uploadColorTable()
{
  if(m_SliceCanvas->m_SliceTiles)
    {
      glBindTexture(GL_TEXTURE_2D, m_SliceCanvas->m_SliceTiles[0].texture);
      glColorTableEXT(GL_TEXTURE_2D, GL_RGBA8, 256, GL_RGBA, GL_UNSIGNED_BYTE, m_SliceCanvas->m_Palette);
    }
}

/************ The SGI color table renderer ************/
SGIColorTableSliceRenderer::SGIColorTableSliceRenderer(SliceCanvas *canvas)
  : SliceRenderer(canvas) {}

SGIColorTableSliceRenderer::~SGIColorTableSliceRenderer() {}

bool SGIColorTableSliceRenderer::init()
{
  /* make sure the proper extensions are initialized */
  if(glewIsSupported("GL_SGIS_texture_edge_clamp") &&
     glewIsSupported("GL_SGI_texture_color_table") &&
     glewIsSupported("GL_SGI_color_table"))
    {
      return true;
    }
  else
    {
      //QMessageBox::critical(this,"Error","Cannot render slices using hardware!",QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
      cvcapp.log(5, "SGIColorTableSliceRenderer::init(): Unable to initialize renderer");
      return false;
    }
}

void SGIColorTableSliceRenderer::draw()
{
  if(!m_SliceCanvas->m_Drawable || !m_SliceCanvas->m_VolumeFileInfo.isSet()) return;

  glColor3f(1.0,1.0,1.0);
  
  glPushAttrib(GL_ENABLE_BIT);
  
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  
  glEnable(GL_TEXTURE_COLOR_TABLE_SGI);
  glEnable(GL_TEXTURE_2D);
  
  for(unsigned int i=0; i<m_SliceCanvas->m_NumSliceTiles; i++)
    {
      glBindTexture(GL_TEXTURE_2D, m_SliceCanvas->m_SliceTiles[i].texture);
	  
      glBegin(GL_QUADS);
      for(unsigned int j=0; j<4; j++)
	{
	  glTexCoord2f(m_SliceCanvas->m_SliceTiles[i].tex_coords[j*2+0],m_SliceCanvas->m_SliceTiles[i].tex_coords[j*2+1]);
	  glVertex2f(m_SliceCanvas->m_SliceTiles[i].vert_coords[j*2+0],m_SliceCanvas->m_SliceTiles[i].vert_coords[j*2+1]);
	}
      glEnd();
    }
 
  /* now draw the point class points as discs */
  glDisable(GL_TEXTURE_2D);
  glDisable(GL_TEXTURE_COLOR_TABLE_SGI);
  glDisable(GL_BLEND);
   
  glPopAttrib();
}

void SGIColorTableSliceRenderer::uploadSlice(unsigned char *slice, unsigned int dimx, unsigned int dimy)
{
  //  unsigned int num_tile_x, num_tile_y,j,a,b,index;
  unsigned int i;
  //  unsigned char tex_buf[TEXTURE_TILE_X*TEXTURE_TILE_Y];
  //  GLuint *texture_ids;
  /*
    float tex_square[8] = { 0.0, 1.0, 1.0, 1.0,        
    1.0, 0.0, 0.0, 0.0 };      
    float vert_square[8] = { -0.5, 0.5, 0.5, 0.5,      
    0.5, -0.5, -0.5, -0.5 };  
			   
  */
  if(!m_SliceCanvas->m_Drawable) return;
  
  m_SliceCanvas->makeCurrent();
  
  m_SliceCanvas->m_NumSliceTiles = 1;
  if(m_SliceCanvas->m_SliceTiles)
    {
      for(i=0; i<m_SliceCanvas->m_NumSliceTiles; i++) glDeleteTextures(1,&(m_SliceCanvas->m_SliceTiles[i]).texture);
      delete [] m_SliceCanvas->m_SliceTiles;
    }
  m_SliceCanvas->m_SliceTiles = new SliceCanvas::SliceTile[m_SliceCanvas->m_NumSliceTiles];
  
  glGenTextures(1,&(m_SliceCanvas->m_SliceTiles[0].texture));
  
  m_SliceCanvas->m_SliceTiles[0].tex_coords[0] = m_SliceCanvas->m_TexCoordMinX; m_SliceCanvas->m_SliceTiles[0].tex_coords[1] = m_SliceCanvas->m_TexCoordMaxY;      /* top left */
  m_SliceCanvas->m_SliceTiles[0].tex_coords[2] = m_SliceCanvas->m_TexCoordMaxX; m_SliceCanvas->m_SliceTiles[0].tex_coords[3] = m_SliceCanvas->m_TexCoordMaxY;      /* top right */
  m_SliceCanvas->m_SliceTiles[0].tex_coords[4] = m_SliceCanvas->m_TexCoordMaxX; m_SliceCanvas->m_SliceTiles[0].tex_coords[5] = m_SliceCanvas->m_TexCoordMinY;      /* bottom right */
  m_SliceCanvas->m_SliceTiles[0].tex_coords[6] = m_SliceCanvas->m_TexCoordMinX; m_SliceCanvas->m_SliceTiles[0].tex_coords[7] = m_SliceCanvas->m_TexCoordMinY;      /* bottom left */
  
  m_SliceCanvas->m_SliceTiles[0].vert_coords[0] = m_SliceCanvas->m_VertCoordMinX; m_SliceCanvas->m_SliceTiles[0].vert_coords[1] = m_SliceCanvas->m_VertCoordMaxY;  /* top left */
  m_SliceCanvas->m_SliceTiles[0].vert_coords[2] = m_SliceCanvas->m_VertCoordMaxX; m_SliceCanvas->m_SliceTiles[0].vert_coords[3] = m_SliceCanvas->m_VertCoordMaxY;  /* top right */
  m_SliceCanvas->m_SliceTiles[0].vert_coords[4] = m_SliceCanvas->m_VertCoordMaxX; m_SliceCanvas->m_SliceTiles[0].vert_coords[5] = m_SliceCanvas->m_VertCoordMinY;  /* bottom right */
  m_SliceCanvas->m_SliceTiles[0].vert_coords[6] = m_SliceCanvas->m_VertCoordMinX; m_SliceCanvas->m_SliceTiles[0].vert_coords[7] = m_SliceCanvas->m_VertCoordMinY;  /* bottom left */
  
  glBindTexture(GL_TEXTURE_2D, m_SliceCanvas->m_SliceTiles[0].texture);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_R8, dimx, dimy, 0, GL_RED, GL_UNSIGNED_BYTE, slice);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  // Set swizzle to replicate red channel to all components for compatibility
  GLint swizzleMask[] = {GL_RED, GL_RED, GL_RED, GL_ONE};
  glTexParameteriv(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_RGBA, swizzleMask);

  uploadColorTable();
}

void SGIColorTableSliceRenderer::uploadColorTable()
{
  if(m_SliceCanvas->m_SliceTiles)
    {
      glBindTexture(GL_TEXTURE_2D, m_SliceCanvas->m_SliceTiles[0].texture);
      glColorTableSGI(GL_TEXTURE_COLOR_TABLE_SGI, GL_RGBA8, 256, GL_RGBA, GL_UNSIGNED_BYTE, m_SliceCanvas->m_Palette);
    }
}

/*********** The VolumeGridRover *****************/
VolumeGridRover::VolumeGridRover(QWidget* parent, Qt::WindowFlags fl)
  : QWidget(parent, fl), m_PointClassList(nullptr),
#ifdef USING_EM_CLUSTERING
    m_Clusters(256,50), m_HistogramCalculated(false),
#endif
    m_RemoteGenSegThread(this), m_LocalGenSegThread(this), m_TilingThread(this),
    m_DoIsocontouring(false)
{
  _ui = new Ui::VolumeGridRoverBase;
  _ui->setupUi( this );

  // TODO: Qt6 - QGLFormat removed, configure OpenGL context via QSurfaceFormat if needed
  // QGLFormat format;
  // format.setAlpha(true);
  // format.setStencil(true);
  
  m_XYSliceCanvas = new SliceCanvas(SliceCanvas::XY,_ui->m_XYSliceFrame);
  m_XYSliceCanvas->setObjectName("m_XYSliceCanvas");
  m_XZSliceCanvas = new SliceCanvas(SliceCanvas::XZ,_ui->m_XZSliceFrame);
  m_XZSliceCanvas->setObjectName("m_XZSliceCanvas");
  m_ZYSliceCanvas = new SliceCanvas(SliceCanvas::ZY,_ui->m_ZYSliceFrame);
  m_ZYSliceCanvas->setObjectName("m_ZYSliceCanvas");
	
  QGridLayout *XYSliceCanvasLayout = new QGridLayout(_ui->m_XYSliceFrame);
  XYSliceCanvasLayout->addWidget(m_XYSliceCanvas);
  
  QGridLayout *XZSliceCanvasLayout = new QGridLayout(_ui->m_XZSliceFrame);
  XZSliceCanvasLayout->addWidget(m_XZSliceCanvas);
  
  QGridLayout *ZYSliceCanvasLayout = new QGridLayout(_ui->m_ZYSliceFrame);
  ZYSliceCanvasLayout->addWidget(m_ZYSliceCanvas);

  _ui->m_X->setValidator(new QIntValidator(_ui->m_X));
  _ui->m_Y->setValidator(new QIntValidator(_ui->m_Y));
  _ui->m_Z->setValidator(new QIntValidator(_ui->m_Z));
  _ui->m_ObjX->setReadOnly(true);
  _ui->m_ObjY->setReadOnly(true);
  _ui->m_ObjZ->setReadOnly(true);

  _ui->m_XYDepthSlide->setMinimum(0);
  _ui->m_XYDepthSlide->setMaximum(0);
  _ui->m_XZDepthSlide->setMinimum(0);
  _ui->m_XZDepthSlide->setMaximum(0);
  _ui->m_ZYDepthSlide->setMinimum(0);
  _ui->m_ZYDepthSlide->setMaximum(0);
  
  setGridCellInfo(0,0,0,
		  0.0,0.0,0.0,
		  0,0,0,0,
		  0);
  setColorName();

 
  connect(m_XYSliceCanvas,SIGNAL(mouseOver(int,int,int,float,float,float,int,int,int,int,double)),
	  SLOT(setGridCellInfo(int,int,int,float,float,float,int,int,int,int,double)));
  connect(m_XZSliceCanvas,SIGNAL(mouseOver(int,int,int,float,float,float,int,int,int,int,double)),
	  SLOT(setGridCellInfo(int,int,int,float,float,float,int,int,int,int,double)));
  connect(m_ZYSliceCanvas,SIGNAL(mouseOver(int,int,int,float,float,float,int,int,int,int,double)),
	  SLOT(setGridCellInfo(int,int,int,float,float,float,int,int,int,int,double)));
  connect(m_XYSliceCanvas,SIGNAL(pointAdded(unsigned int,unsigned int,int,unsigned int,unsigned int,unsigned int)),
	  SIGNAL(pointAdded(unsigned int,unsigned int,int,unsigned int,unsigned int,unsigned int)));
  connect(m_XYSliceCanvas,SIGNAL(pointRemoved(unsigned int,unsigned int,int,unsigned int,unsigned int,unsigned int)),
	  SIGNAL(pointRemoved(unsigned int,unsigned int,int,unsigned int,unsigned int,unsigned int)));
  connect(m_XZSliceCanvas,SIGNAL(pointAdded(unsigned int,unsigned int,int,unsigned int,unsigned int,unsigned int)),
	  SIGNAL(pointAdded(unsigned int,unsigned int,int,unsigned int,unsigned int,unsigned int)));
  connect(m_XZSliceCanvas,SIGNAL(pointRemoved(unsigned int,unsigned int,int,unsigned int,unsigned int,unsigned int)),
	  SIGNAL(pointRemoved(unsigned int,unsigned int,int,unsigned int,unsigned int,unsigned int)));
  connect(m_ZYSliceCanvas,SIGNAL(pointAdded(unsigned int,unsigned int,int,unsigned int,unsigned int,unsigned int)),
	  SIGNAL(pointAdded(unsigned int,unsigned int,int,unsigned int,unsigned int,unsigned int)));
  connect(m_ZYSliceCanvas,SIGNAL(pointRemoved(unsigned int,unsigned int,int,unsigned int,unsigned int,unsigned int)),
	  SIGNAL(pointRemoved(unsigned int,unsigned int,int,unsigned int,unsigned int,unsigned int)));
	  

  connect(_ui->m_Variable,SIGNAL(activated(int)),SLOT(setCurrentVariable(int)));
  connect(_ui->m_Timestep,SIGNAL(valueChanged(int)),SLOT(setCurrentTimestep(int)));
  connect(_ui->m_Variable,SIGNAL(activated(int)),m_XYSliceCanvas,SLOT(setCurrentVariable(int)));
  connect(_ui->m_Variable,SIGNAL(activated(int)),m_XZSliceCanvas,SLOT(setCurrentVariable(int)));
  connect(_ui->m_Variable,SIGNAL(activated(int)),m_ZYSliceCanvas,SLOT(setCurrentVariable(int)));
  connect(_ui->m_Timestep,SIGNAL(valueChanged(int)),m_XYSliceCanvas,SLOT(setCurrentTimestep(int)));
  connect(_ui->m_Timestep,SIGNAL(valueChanged(int)),m_XZSliceCanvas,SLOT(setCurrentTimestep(int)));
  connect(_ui->m_Timestep,SIGNAL(valueChanged(int)),m_ZYSliceCanvas,SLOT(setCurrentTimestep(int)));
	  
  connect(_ui->m_XYDepthSlide,SIGNAL(sliderMoved(int)),m_XYSliceCanvas,SLOT(setDepth(int)));
  connect(_ui->m_XZDepthSlide,SIGNAL(sliderMoved(int)),m_XZSliceCanvas,SLOT(setDepth(int)));
  connect(_ui->m_ZYDepthSlide,SIGNAL(sliderMoved(int)),m_ZYSliceCanvas,SLOT(setDepth(int)));
  connect(_ui->m_XYDepthSlide,SIGNAL(valueChanged(int)),m_XYSliceCanvas,SLOT(setDepth(int)));
  connect(_ui->m_XZDepthSlide,SIGNAL(valueChanged(int)),m_XZSliceCanvas,SLOT(setDepth(int)));
  connect(_ui->m_ZYDepthSlide,SIGNAL(valueChanged(int)),m_ZYSliceCanvas,SLOT(setDepth(int)));
  connect(m_XYSliceCanvas,SIGNAL(depthChanged(int)),SLOT(setZ(int)));
  connect(m_XZSliceCanvas,SIGNAL(depthChanged(int)),SLOT(setY(int)));
  connect(m_ZYSliceCanvas,SIGNAL(depthChanged(int)),SLOT(setX(int)));
  connect(m_XYSliceCanvas,SIGNAL(depthChanged(int)),_ui->m_XYDepthSlide,SLOT(setValue(int)));
  connect(m_XZSliceCanvas,SIGNAL(depthChanged(int)),_ui->m_XZDepthSlide,SLOT(setValue(int)));
  connect(m_ZYSliceCanvas,SIGNAL(depthChanged(int)),_ui->m_ZYDepthSlide,SLOT(setValue(int)));
  connect(m_XYSliceCanvas,SIGNAL(depthChanged(int)),SLOT(ZDepthChangedSlot(int)));
  connect(m_XZSliceCanvas,SIGNAL(depthChanged(int)),SLOT(YDepthChangedSlot(int)));
  connect(m_ZYSliceCanvas,SIGNAL(depthChanged(int)),SLOT(XDepthChangedSlot(int)));
  connect(_ui->m_XYResetViewButton,SIGNAL(clicked()),m_XYSliceCanvas,SLOT(resetView()));
  connect(_ui->m_XZResetViewButton,SIGNAL(clicked()),m_XZSliceCanvas,SLOT(resetView()));
  connect(_ui->m_ZYResetViewButton,SIGNAL(clicked()),m_ZYSliceCanvas,SLOT(resetView()));
  
  //===============================================================================================
  // connect variable selection tab
  connect(_ui->m_Variable,SIGNAL(activated(int)),SLOT(setCurrentVariable(int)));
  connect(_ui->m_Timestep,SIGNAL(valueChanged(int)),SLOT(setCurrentTimestep(int)));
  connect(_ui->m_Variable,SIGNAL(activated(int)),m_XYSliceCanvas,SLOT(setCurrentVariable(int)));
  connect(_ui->m_Variable,SIGNAL(activated(int)),m_XZSliceCanvas,SLOT(setCurrentVariable(int)));
  connect(_ui->m_Variable,SIGNAL(activated(int)),m_ZYSliceCanvas,SLOT(setCurrentVariable(int)));
  connect(_ui->m_Timestep,SIGNAL(valueChanged(int)),m_XYSliceCanvas,SLOT(setCurrentTimestep(int)));
  connect(_ui->m_Timestep,SIGNAL(valueChanged(int)),m_XZSliceCanvas,SLOT(setCurrentTimestep(int)));
  connect(_ui->m_Timestep,SIGNAL(valueChanged(int)),m_ZYSliceCanvas,SLOT(setCurrentTimestep(int)));

  //-----------------------------------------------------------------------------------------------
  // connet Grid Cell Marking Tab
     // connect point class tab
     connect(_ui->m_PointClass,SIGNAL(activated(int)),m_XYSliceCanvas,SLOT(setCurrentClass(int)));
     connect(_ui->m_PointClass,SIGNAL(activated(int)),m_XZSliceCanvas,SLOT(setCurrentClass(int)));
     connect(_ui->m_PointClass,SIGNAL(activated(int)),m_ZYSliceCanvas,SLOT(setCurrentClass(int)));
     connect(_ui->m_PointClass,SIGNAL(activated(int)),SLOT(showColor(int)));

     _ui->m_PointClassColor->setStyle( QStyleFactory::create("plastique") );
     QPalette palette = _ui->m_PointClassColor->palette();
     palette.setColor( QPalette::Normal, QPalette::Button, QColor(13, 208, 200) );
     _ui->m_PointClassColor->setPalette( palette );
     connect((const class QObject*)(_ui->m_PointClassColor), SIGNAL(clicked()), this, SLOT(colorSlot()));

     connect(_ui->m_AddPointClass, SIGNAL( clicked() ), SLOT( addPointClassSlot() ) );
     connect(_ui->m_DeletePointClass, SIGNAL( clicked() ), SLOT( deletePointClassSlot() ) );
     connect(_ui->m_PointClassesSaveButton, SIGNAL( clicked() ), SLOT( savePointClassesSlot() ) );
     connect(_ui->m_PointClassesLoadButton, SIGNAL( clicked() ), SLOT( loadPointClassesSlot() ) );
 
     // connect contour tab
     connect(_ui->m_AddContour, SIGNAL( clicked() ), SLOT( addContourSlot() ) );
     connect(_ui->m_DeleteContour, SIGNAL( clicked() ), SLOT( deleteContourSlot() ) );
     _ui->m_ContourColor->setStyle( QStyleFactory::create("plastique") );
     palette = _ui->m_ContourColor->palette();
     palette.setColor( QPalette::Normal, QPalette::Button, QColor(128, 255, 128) );
     _ui->m_ContourColor->setPalette( palette );
     connect(_ui->m_ContourColor, SIGNAL( clicked() ), SLOT( contourColorSlot() ) );
     connect(_ui->m_InterpolationType, SIGNAL( activated(int) ), SLOT( setInterpolationTypeSlot(int) ) );
     _ui->m_InterpolationType->addItem( QString("Linear") );
     connect(_ui->m_InterpolationSampling, SIGNAL( valueChanged(int) ), SLOT( setInterpolationSamplingSlot(int) ) );
 
     connect(_ui->m_SaveContoursButton, SIGNAL( clicked() ), SLOT( saveContoursSlot() ) );
     // connect(_ui->m_LoadContourButton, SIGNAL( clicked() ), SLOT( loadContoursSlot() ) );
     // connect(_ui->m_tileLoadButton, SIGNAL( clicked() ), SLOT( loadContoursSlot() ) );

  connect(_ui->m_GreyScale,SIGNAL(toggled(bool)),m_XYSliceCanvas,SLOT(setGreyScale(bool)));
  connect(_ui->m_GreyScale,SIGNAL(toggled(bool)),m_XZSliceCanvas,SLOT(setGreyScale(bool)));
  connect(_ui->m_GreyScale,SIGNAL(toggled(bool)),m_ZYSliceCanvas,SLOT(setGreyScale(bool)));

  connect(_ui->m_RenderSDF,SIGNAL(toggled(bool)),m_XYSliceCanvas,SLOT(setRenderSDF(bool)));
  connect(_ui->m_RenderSDF,SIGNAL(toggled(bool)),m_XZSliceCanvas,SLOT(setRenderSDF(bool)));
  connect(_ui->m_RenderSDF,SIGNAL(toggled(bool)),m_ZYSliceCanvas,SLOT(setRenderSDF(bool)));

  connect(_ui->m_RenderControlPoints,SIGNAL(toggled(bool)),m_XYSliceCanvas,SLOT(setRenderControlPoints(bool)));
  connect(_ui->m_RenderControlPoints,SIGNAL(toggled(bool)),m_XZSliceCanvas,SLOT(setRenderControlPoints(bool)));
  connect(_ui->m_RenderControlPoints,SIGNAL(toggled(bool)),m_ZYSliceCanvas,SLOT(setRenderControlPoints(bool)));
  
  connect(_ui->m_PointSize,SIGNAL(sliderMoved(int)),SLOT(setPointSize(int)));

  _ui->m_BackgroundColor->setStyle( QStyleFactory::create("plastique") );
  palette = _ui->m_BackgroundColor->palette();
  palette.setColor( QPalette::Normal, QPalette::Button, QColor(64, 64, 64) );
  _ui->m_BackgroundColor->setPalette( palette );
  connect(_ui->m_BackgroundColor,SIGNAL(clicked()),this, SLOT(backgroundColorSlot()));

  //-----------------------------------------------------------------------------------------------
  // connet Segmentation Tab

  connect(_ui->m_LocalRun,SIGNAL(clicked()), this, SLOT(localSegmentationRunSlot()));
  connect(_ui->m_RemoteRun,SIGNAL(clicked()), this, SLOT(remoteSegmentationRunSlot()));

  //-----------------------------------------------------------------------------------------------

  connect(_ui->m_Objects,SIGNAL(currentItemChanged(QListWidgetItem*,QListWidgetItem*)),SLOT(currentObjectSelectionChanged()));
  
  _ui->m_ThresholdLow->setValidator(new QIntValidator(0,255, _ui->m_ThresholdLow));
  _ui->m_ThresholdHigh->setValidator(new QIntValidator(0,255, _ui->m_ThresholdHigh));
  _ui->m_Port->setValidator(new QIntValidator(0,65535,_ui->m_Port));

  sliceAxisChangedSlot(); /* call this to enable the correct xyz line edit widget */

  connect(this,SIGNAL(cellMarkingModeChanged(int)),SLOT(setCellMarkingMode(int)));
  connect(this,SIGNAL(cellMarkingModeChanged(int)),m_XYSliceCanvas,SLOT(setCellMarkingMode(int)));
  connect(this,SIGNAL(cellMarkingModeChanged(int)),m_XZSliceCanvas,SLOT(setCellMarkingMode(int)));
  connect(this,SIGNAL(cellMarkingModeChanged(int)),m_ZYSliceCanvas,SLOT(setCellMarkingMode(int)));

  // TODO: setInterpolationType slot is commented out, uncomment if needed
  // connect(_ui->m_InterpolationType,SIGNAL(activated(int)),m_XYSliceCanvas,SLOT(setInterpolationType(int)));
  // connect(_ui->m_InterpolationType,SIGNAL(activated(int)),m_XZSliceCanvas,SLOT(setInterpolationType(int)));
  // connect(_ui->m_InterpolationType,SIGNAL(activated(int)),m_ZYSliceCanvas,SLOT(setInterpolationType(int)));

  setCellMarkingMode(0); //start with point class marking
  m_XYSliceCanvas->setCellMarkingMode(0);
  m_XZSliceCanvas->setCellMarkingMode(0);
  m_ZYSliceCanvas->setCellMarkingMode(0);

  connect(_ui->m_Objects,SIGNAL(itemSelectionChanged()),SLOT(setSelectedContours()));
  connect(_ui->m_Isocontouring,SIGNAL(toggled(bool)),SLOT(doIsocontouring(bool)));

  // set varialbles
  m_hasVolume = false;

  // arand, 4-21-2011: get initial colortable info
  m_XYSliceCanvas->setGreyScale(false);
  m_ZYSliceCanvas->setGreyScale(false);
  m_XZSliceCanvas->setGreyScale(false);

}

VolumeGridRover::~VolumeGridRover()
{
  unsetVolume();
}

void VolumeGridRover::colorSlot()
{
  if( (m_PointClassList == NULL) || (_ui->m_PointClass->currentIndex() < 0) ) return;

  QPalette palette = _ui->m_PointClassColor->palette();
  QColor color = QColorDialog::getColor( palette.color(QPalette::Button) );
  if(color.isValid())
    {
      palette.setColor( QPalette::Normal, QPalette::Button, color );
      _ui->m_PointClassColor->setPalette( palette );
      if(m_PointClassList[_ui->m_Variable->currentIndex()][_ui->m_Timestep->value()]->at(_ui->m_PointClass->currentIndex()))
	m_PointClassList[_ui->m_Variable->currentIndex()][_ui->m_Timestep->value()]->at(_ui->m_PointClass->currentIndex())->setColor(color);
    }
}

void VolumeGridRover::addPointClassSlot()
{
  static int classNum = 0;
  if(m_PointClassList == NULL) return;

  m_PointClassList[_ui->m_Variable->currentIndex()][_ui->m_Timestep->value()]->append(new PointClass(_ui->m_PointClassColor->palette().color( QPalette::Button ),QString("Class %1").arg(classNum)));
  _ui->m_PointClass->insertItem(_ui->m_PointClass->count(), QString("Class %1").arg(classNum));
  _ui->m_PointClass->setCurrentIndex(_ui->m_PointClass->count()-1);

  m_XYSliceCanvas->setCurrentClass(_ui->m_PointClass->currentIndex());
  m_XZSliceCanvas->setCurrentClass(_ui->m_PointClass->currentIndex());
  m_ZYSliceCanvas->setCurrentClass(_ui->m_PointClass->currentIndex());
  showColor(_ui->m_PointClass->currentIndex() );

  classNum++;
}

void VolumeGridRover::deletePointClassSlot()
{
  if(m_PointClassList == NULL) return;

  m_PointClassList[_ui->m_Variable->currentIndex()][_ui->m_Timestep->value()]->removeAt(_ui->m_PointClass->currentIndex());
  _ui->m_PointClass->removeItem(_ui->m_PointClass->currentIndex());

  int currentItem = _ui->m_PointClass->currentIndex();
  if( currentItem >= 0 )
  {
    m_XYSliceCanvas->setCurrentClass( currentItem ); /* make sure other objects dont point to the removed point class */
    m_XZSliceCanvas->setCurrentClass( currentItem );
    m_ZYSliceCanvas->setCurrentClass( currentItem );
    showColor( currentItem);
  }

  m_XYSliceCanvas->update();
  m_XZSliceCanvas->update();
  m_ZYSliceCanvas->update();
}

void VolumeGridRover::xChangedSlot()
{
  m_ZYSliceCanvas->setDepth(_ui->m_X->text().toInt());
}

void VolumeGridRover::yChangedSlot()
{
  m_XZSliceCanvas->setDepth(_ui->m_Y->text().toInt());
}

void VolumeGridRover::zChangedSlot()
{
  m_XYSliceCanvas->setDepth(_ui->m_Z->text().toInt());
}

void VolumeGridRover::backgroundColorSlot()
{
  QPalette palette = _ui->m_BackgroundColor->palette();
  QColor color = QColorDialog::getColor( palette.color(QPalette::Button) );
  if(color.isValid())
    {
      palette.setColor( QPalette::Normal, QPalette::Button, color );
      _ui->m_BackgroundColor->setPalette( palette );

      m_XYSliceCanvas->makeCurrent();
      m_XYSliceCanvas->setBackgroundColor( color );
      m_XZSliceCanvas->makeCurrent();
      m_XZSliceCanvas->setBackgroundColor( color );
      m_ZYSliceCanvas->makeCurrent();
      m_ZYSliceCanvas->setBackgroundColor( color );
      m_XYSliceCanvas->update();
      m_XZSliceCanvas->update();
      m_ZYSliceCanvas->update();
    }
}

void VolumeGridRover::addContour(const SurfRecon::Contour& c)
{
  if(m_Contours.empty()) return;
  m_Contours[_ui->m_Variable->currentIndex()][_ui->m_Timestep->value()][c.name()] = SurfRecon::ContourPtr(new SurfRecon::Contour(c));
  _ui->m_Objects->insertItem( _ui->m_Objects->count(), new QListWidgetItem( QString(c.name().data()), _ui->m_Objects, 0));
  _ui->m_Objects->setCurrentItem(_ui->m_Objects->findItems( QString(c.name().data()), Qt::MatchExactly).first());

  updateSliceContours();
  showCurrentObject();
  update();
}

void VolumeGridRover::removeContour(const std::string& name)
{
  m_Contours[_ui->m_Variable->currentIndex()][_ui->m_Timestep->value()].erase(name);
  delete _ui->m_Objects->findItems( QString(name.data()) , Qt::MatchExactly ).first();

  updateSliceContours();
  showCurrentObject();
  update();
}

void VolumeGridRover::showCurrentObject()
{
  m_XYSliceCanvas->setCurrentContour(_ui->m_Objects->currentItem() ? _ui->m_Objects->currentItem()->text().toStdString().c_str() : "");
  m_XZSliceCanvas->setCurrentContour(_ui->m_Objects->currentItem() ? _ui->m_Objects->currentItem()->text().toStdString().c_str() : "");
  m_ZYSliceCanvas->setCurrentContour(_ui->m_Objects->currentItem() ? _ui->m_Objects->currentItem()->text().toStdString().c_str() : "");

  showContourColor(_ui->m_Objects->currentItem() ? _ui->m_Objects->currentItem()->text().toStdString().c_str() : "");
  showContourInterpolationType(_ui->m_Objects->currentItem() ? _ui->m_Objects->currentItem()->text().toStdString().c_str() : "");
  showContourInterpolationSampling(_ui->m_Objects->currentItem() ? _ui->m_Objects->currentItem()->text().toStdString().c_str() : "");
}

void VolumeGridRover::updateSliceContours()
{
  m_XYSliceCanvas->setContours(m_Contours);
  m_XZSliceCanvas->setContours(m_Contours);
  m_ZYSliceCanvas->setContours(m_Contours);
}

QString VolumeGridRover::cacheDir() const
{
  QSettings settings(QDir::homePath() + "/.CCV/settings.ini", QSettings::IniFormat);
  QString cacheDirString = settings.value("/Volume Rover/CacheDir", QString()).toString();
  if(cacheDirString.isEmpty()) cacheDirString = ".";
  cacheDirString += "/VolumeCache";
  return cacheDirString;
}

void VolumeGridRover::isocontourValues(const IsocontourValues& vals)
{
  m_IsocontourValues = vals;
  if(doIsocontouring())
    {
      clearIsocontours();
      generateIsocontours();
    }
}

void VolumeGridRover::doIsocontouring(bool flag)
{
  m_DoIsocontouring = flag;
  if(flag)
    generateIsocontours();
  else
    clearIsocontours();
}

void VolumeGridRover::clearIsocontours()
{
  //remove all contours with the name "isocontour_%f"
  int size = _ui->m_Objects->count();

  for( int i=0; i < size; i++ )
  {
      QListWidgetItem *item = _ui->m_Objects->item( i );
      QStringList split = item->text().split("_");
     
      if( split[ 0 ] == "isocontour" )
	{
	  m_Contours[_ui->m_Variable->currentIndex()][_ui->m_Timestep->value()].erase( item->text().toStdString().c_str());
	  delete item;
	}
    }
}

void VolumeGridRover::generateIsocontours()
{
#ifdef VOLUMEGRIDROVER_ISOCONTOURING
  if(m_IsocontourValues.empty())
    {
      cvcapp.log(5, boost::str(boost::format("%s: no isocontour values")% BOOST_CURRENT_FUNCTION));
      return;
    };

  std::string val_list;
  for(IsocontourValues::iterator i = m_IsocontourValues.begin();
      i != m_IsocontourValues.end();
      i++)
    val_list += boost::str(boost::format("%1% ") % i->value);
  boost::algorithm::trim(val_list);
  cvcapp.log(5, 
	     boost::str(boost::format("%s\n")%
    boost::str(
      boost::format("VolumeGridRover: isocontouring the following values :: %1%")
      % val_list
    ).c_str()
			));

  SliceCanvas *curCanvas = getCurrentSliceCanvas();
  VolMagick::Volume sliceData = curCanvas->sliceData();
  VolMagick::Volume normalizedSliceData;

  //need to orient data such that it is always XY.. convert coords back after
  //doing the isocontouring.
  switch(curCanvas->getSliceAxis())
    {
    case SliceCanvas::XY:
      {
	//nothing to do here except copy.
	normalizedSliceData = sliceData;
	break;
      }
    case SliceCanvas::XZ:
      {
	normalizedSliceData.dimension(VolMagick::Dimension(sliceData.XDim(),
							   sliceData.ZDim(),
							   1));
	normalizedSliceData.boundingBox(VolMagick::BoundingBox(sliceData.XMin(),
							       sliceData.ZMin(),
							       0.0,
							       sliceData.XMax(),
							       sliceData.ZMax(),
							       0.0));
	for(VolMagick::uint64 i = 0; i < sliceData.XDim(); i++)
	  for(VolMagick::uint64 j = 0; j < sliceData.ZDim(); j++)
	    normalizedSliceData(i,j,0, sliceData(i,0,j));
	break;
      }
    case SliceCanvas::ZY:
      {
	normalizedSliceData.dimension(VolMagick::Dimension(sliceData.ZDim(),
							   sliceData.YDim(),
							   1));
	normalizedSliceData.boundingBox(VolMagick::BoundingBox(sliceData.ZMin(),
							       sliceData.YMin(),
							       0.0,
							       sliceData.ZMax(),
							       sliceData.YMax(),
							       0.0));
	for(VolMagick::uint64 i = 0; i < sliceData.ZDim(); i++)
	  for(VolMagick::uint64 j = 0; j < sliceData.YDim(); j++)
	    normalizedSliceData(i,j,0, sliceData(0,j,i));
	break;
      }
    }

  ConDataset *the_data;
  int dim[2] = { normalizedSliceData.XDim(),normalizedSliceData.YDim() };

  //convert it to a supported libcontour type and load it into libisocontour
  switch(normalizedSliceData.voxelType())
    {
    case VolMagick::UInt:
    case VolMagick::Double:
    case VolMagick::UInt64:
    case VolMagick::Float:
      {
	normalizedSliceData.voxelType(VolMagick::Float);
	the_data = newDatasetReg(CONTOUR_FLOAT, CONTOUR_REG_2D, 1, 1, dim, *normalizedSliceData);
	break;
      }
    case VolMagick::UChar:
      {
	the_data = newDatasetReg(CONTOUR_UCHAR, CONTOUR_REG_2D, 1, 1, dim, *normalizedSliceData);
	break;
      }
    case VolMagick::UShort:
      {
	the_data = newDatasetReg(CONTOUR_USHORT, CONTOUR_REG_2D, 1, 1, dim, *normalizedSliceData);
	break;
      }
    }

  float span[2] = { normalizedSliceData.XSpan(), normalizedSliceData.YSpan() };
  float orig[2] = { normalizedSliceData.XMin(), normalizedSliceData.YMin() };
  ((Datareg2 *)the_data->data->getData(0))->setOrig(orig);
  ((Datareg2 *)the_data->data->getData(0))->setSpan(span);

  for(IsocontourValues::iterator i = m_IsocontourValues.begin();
      i != m_IsocontourValues.end();
      i++)
    {
      Contour2dData *isocontour = getContour2d(the_data,
					       0, 0,
      					       i->value);
      SurfRecon::Contour contour(SurfRecon::Color(255.0,255.0,255.0),
				 boost::str(boost::format("isocontour_%1%")%i->value),
				 SurfRecon::InterpolationType(_ui->m_InterpolationType->currentIndex()),
				 _ui->m_InterpolationSampling->value());
      for(int edge_idx = 0; edge_idx < isocontour->nedge; edge_idx++)
	{
	  SurfRecon::Curve c(curCanvas->getDepth(),
			     SurfRecon::Orientation(curCanvas->getSliceAxis()),
			     SurfRecon::PointPtrList(),
			     contour.name());
	  
	  float *pt0,*pt1;
	  u_int *p_idx = isocontour->edge[edge_idx];
	  pt0 = isocontour->vert[p_idx[0]];
	  pt1 = isocontour->vert[p_idx[1]];

	  switch(curCanvas->getSliceAxis())
	    {
	    case SliceCanvas::XY:
	      {
		double depth_coord = curCanvas->getDepth()*m_VolumeFileInfo.ZSpan()+m_VolumeFileInfo.ZMin();
		getCurvePoints(c).push_back(SurfRecon::PointPtr(new SurfRecon::Point(pt0[0],
										     pt0[1],
										     depth_coord)));
		getCurvePoints(c).push_back(SurfRecon::PointPtr(new SurfRecon::Point(pt1[0],
										     pt1[1],
										     depth_coord)));
	      }
	      break;
	    case SliceCanvas::XZ:
	      {
		double depth_coord = curCanvas->getDepth()*m_VolumeFileInfo.YSpan()+m_VolumeFileInfo.YMin();
		getCurvePoints(c).push_back(SurfRecon::PointPtr(new SurfRecon::Point(pt0[0],
										     depth_coord,
										     pt0[1])));
		getCurvePoints(c).push_back(SurfRecon::PointPtr(new SurfRecon::Point(pt1[0],
										     depth_coord,
										     pt1[1])));
	      }
	      break;
	    case SliceCanvas::ZY:
	      {
		double depth_coord = curCanvas->getDepth()*m_VolumeFileInfo.XSpan()+m_VolumeFileInfo.XMin();
		getCurvePoints(c).push_back(SurfRecon::PointPtr(new SurfRecon::Point(depth_coord,
										     pt0[1],
										     pt0[0])));
		getCurvePoints(c).push_back(SurfRecon::PointPtr(new SurfRecon::Point(depth_coord,
										     pt1[1],
										     pt1[0])));
	      }
	      break;
	    }

	  contour.add(SurfRecon::CurvePtr(new SurfRecon::Curve(c)));
	}

      addContour(contour);
      delete isocontour;
    }

  clearDataset(the_data);
  
#else
  cvcapp.log(5, "VolumeGridRover: Missing isocontouring support");
#endif
}

//Returns the current visible slice canvas
SliceCanvas *VolumeGridRover::getCurrentSliceCanvas()
{
  switch(_ui->m_SliceCanvasTab->currentIndex() )
    {
    default: return NULL;
    case 0: return m_XYSliceCanvas; break;
    case 1: return m_XZSliceCanvas; break;
    case 2: return m_ZYSliceCanvas; break;
    }
}

void VolumeGridRover::saveContoursSlot()
{
  //cvcapp.log(5, "VolumeGridRover::saveContoursSlot()");
  unsigned int var, time;

  QString filename = QFileDialog::getSaveFileName(this,
						  "Choose a filename to save under",
						  QString(),
						  "VolumeRover Contour files (*.cnt);;"
						  "ctr2suf Contour format (*.contour);;"
						  "Reconstruction point cloud (*.pcd);;"
						  "B-Spline contours (*.bspline);;"
						  "MAT-1.0.x Skeletonization Software contours (*.matc);;"
						  "All Files (*)");
  
  if(filename == QString())
    {
      cvcapp.log(5, "VolumeGridRover::saveContoursSlot(): save cancelled (or null filename)");
      return;
    }
  
  if(filename.endsWith(".pcd")) //samrat's point cloud dataset.. basically the raw points of contours without connectivity
    {
      if(!m_Contours.empty())
	{
	  unsigned int numSteps = 0;
	  
	  //calculate the number of steps
	  for(var = 0; var < m_VolumeFileInfo.numVariables(); var++)
	    for(time = 0; time < m_VolumeFileInfo.numTimesteps(); time++)
	      numSteps += m_Contours[var][time].size();

	  QProgressDialog progressDialog("Writing point cloud data for each objects set of contours...", 
					 "Cancel", 0, numSteps, this);
	  progressDialog.setValue(0);

	  for(var = 0; var < m_VolumeFileInfo.numVariables(); var++)
	    for(time = 0; time < m_VolumeFileInfo.numTimesteps(); time++)
	      {
		for(SurfRecon::ContourPtrMap::iterator i = m_Contours[var][time].begin();
		    i != m_Contours[var][time].end();
		    i++)
		  {
		    progressDialog.setValue(progressDialog.value()+1);
		    qApp->processEvents();
		    if(progressDialog.wasCanceled())
		      return;

		    if((*i).second == NULL) continue;

		    QString fileContents;
		    unsigned int pointCount = 0;
		    std::vector<double> vec_x, vec_y, vec_z, vec_t; //used to build the input arrays for gsl_spline_init()
		    double t, interval;
		    const gsl_interp_type *interp;
		    gsl_interp_accel *acc_x, *acc_y, *acc_z;
		    gsl_spline *spline_x, *spline_y, *spline_z;

		    switch((*i).second->interpolationType())
		      {
		      case 0: interp = gsl_interp_linear; break;
		      case 1: interp = gsl_interp_polynomial; break;
		      case 2: interp = gsl_interp_cspline; break;
		      case 3: interp = gsl_interp_cspline_periodic; break;
		      case 4: interp = gsl_interp_akima; break;
		      case 5: interp = gsl_interp_akima_periodic; break;
		      default: interp = gsl_interp_cspline; break;
		      }

		    for(std::vector<SurfRecon::CurvePtr>::iterator cur = (*i).second->curves().begin();
			cur != (*i).second->curves().end();
			cur++)
		      {
			if(SurfRecon::getCurvePoints(**cur).empty()) continue;
			
			t = 0.0;
			vec_x.clear(); vec_y.clear(); vec_z.clear(); vec_t.clear();
			

			for(SurfRecon::PointPtrList::iterator pcur = SurfRecon::getCurvePoints(**cur).begin();
			    pcur != SurfRecon::getCurvePoints(**cur).end();
			    pcur++)
			  {
			    vec_x.push_back((*pcur)->x());
			    vec_y.push_back((*pcur)->y());
			    vec_z.push_back((*pcur)->z());
			    vec_t.push_back(t); t+=1.0;

			    /*
			    fileContents += QString("%1 %2 %3\n")
			      .arg((*pcur)->x())
			      .arg((*pcur)->y())
			      .arg((*pcur)->z());
			    */
			  }
			t-=1.0;
			//this is hackish but how else can i find the min size?
			gsl_spline *test_spline = gsl_spline_alloc(interp, 1000);
			if(vec_t.size() >= gsl_spline_min_size(test_spline))
			  {
			    double time_i;
			    SurfRecon::PointPtrList tmplist;
			    //int count;
			      
			    acc_x = gsl_interp_accel_alloc();
			    acc_y = gsl_interp_accel_alloc();
			    acc_z = gsl_interp_accel_alloc();
			    spline_x = gsl_spline_alloc(interp, vec_t.size());
			    spline_y = gsl_spline_alloc(interp, vec_t.size());
			    spline_z = gsl_spline_alloc(interp, vec_t.size());
			    gsl_spline_init(spline_x, &(vec_t[0]), &(vec_x[0]), vec_t.size());
			    gsl_spline_init(spline_y, &(vec_t[0]), &(vec_y[0]), vec_t.size());
			    gsl_spline_init(spline_z, &(vec_t[0]), &(vec_z[0]), vec_t.size());
			    
			    //Just sample some number of points for non-linear interpolation for now.
			    interval = 1.0/(1+(*i).second->numberOfSamples());
			    //In the future, sample according to the spline's curvature between 2
			    // control points...

			    //this bothers me... please just calculate the number of points in the future!
			    //for(time_i = 0.0, count=0; time_i<=t; time_i+=interval, count++);
			    //outfile << count << std::endl;
			    
			    for(time_i = 0.0; time_i<=t; time_i+=interval)
			      {
				/*outfile << gsl_spline_eval(spline_x, time_i, acc_x) << " " 
					<< gsl_spline_eval(spline_y, time_i, acc_y) << " " 
					<< gsl_spline_eval(spline_z, time_i, acc_z) << std::endl;*/
				tmplist.push_back(SurfRecon::PointPtr(new SurfRecon::Point(gsl_spline_eval(spline_x, time_i, acc_x),
											   gsl_spline_eval(spline_y, time_i, acc_y),
											   gsl_spline_eval(spline_z, time_i, acc_z))));
			      }
			      
			    /* clean up our mess */
			    gsl_spline_free(spline_x);
			    gsl_spline_free(spline_y);
			    gsl_spline_free(spline_z);
			    gsl_interp_accel_free(acc_x);
			    gsl_interp_accel_free(acc_y);
			    gsl_interp_accel_free(acc_z);

			    SurfRecon::PointPtr pfront = tmplist.front();
			    SurfRecon::PointPtr pback = tmplist.back();

			    //remove the last point if the same as first
			    if((fabs(pfront->x() - pback->x()) <= 0.00001) &&
			       (fabs(pfront->y() - pback->y()) <= 0.00001) &&
			       (fabs(pfront->z() - pback->z()) <= 0.00001))
			      tmplist.pop_back();

			    //actually output the points
			    for(SurfRecon::PointPtrList::const_iterator pcur = tmplist.begin();
				pcur != tmplist.end();
				pcur++)
			      {
				//outfile << (*pcur)->x() << " " << (*pcur)->y() << " " << (*pcur)->z() << std::endl;
				fileContents += QString("%1 %2 %3\n")
				  .arg((*pcur)->x())
				  .arg((*pcur)->y())
				  .arg((*pcur)->z());
			      }
			    pointCount += tmplist.size();
			  }
			gsl_spline_free(test_spline);
		      }

		    QString new_filename(filename);
		    new_filename.replace(".pcd",QString("_%1.pcd").arg("%1f").arg( QString((*i).second->name().data() ) ) );
		    QFile f(new_filename);
		    if(!f.open(QIODevice::WriteOnly))
		      {
			QMessageBox::critical(this,"Error",QString("Could not open the file %1!").arg(new_filename),
					      QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
			return;
		      }

		    QTextStream stream(&f);
		    stream << QString("%1\n").arg(pointCount) + fileContents;
		    f.close();
		  }
	      }
	}
    }
  else if(filename.endsWith(".contour")) //taoju's contour data for his program ctr2suf...
    {
      //this format allows contours on arbitrary planes rather than parallel planes
      //however here in VolumeRover land we're only dealing with parallel planes...
      //TODO: resample using interpolation as we did above

      if(!m_Contours.empty())
	{
	  typedef CGAL::Simple_cartesian<double> Kernel;
	  typedef Kernel::Point_3 Point_3;
	  typedef Kernel::Plane_3 Plane_3;

	  typedef std::set<Point_3> Contour_vertex_set;
	  typedef std::vector<Point_3> Contour_vertex_vector;
	  typedef boost::tuple<Contour_vertex_set::iterator,Contour_vertex_set::iterator,int,int> Contour_edge;
	  typedef std::vector<Contour_edge> Contour_edges;
	  typedef boost::tuple<double,double,double,double> Contour_basic_plane; //use this to take advantage of operator< 
	  typedef boost::tuple<Contour_basic_plane, Contour_vertex_set, Contour_edges> Contour_plane;
	  typedef std::map<Contour_basic_plane, Contour_plane> Contour_planes;
	  
	  unsigned int numSteps = 0;
	  
	  //calculate the number of steps
	  for(var = 0; var < m_VolumeFileInfo.numVariables(); var++)
	    for(time = 0; time < m_VolumeFileInfo.numTimesteps(); time++)
	      numSteps += m_Contours[var][time].size();

	  QProgressDialog progressDialog("Writing contour files for each object...", "Cancel", 0, numSteps, this );
	  progressDialog.setValue(0);

	  for(var = 0; var < m_VolumeFileInfo.numVariables(); var++)
	    for(time = 0; time < m_VolumeFileInfo.numTimesteps(); time++)
	      {
		for(SurfRecon::ContourPtrMap::iterator i = m_Contours[var][time].begin();
		    i != m_Contours[var][time].end();
		    i++)
		  {
		    progressDialog.setValue(progressDialog.value()+1);
		    qApp->processEvents();
		    if(progressDialog.wasCanceled())
		      return;

		    if((*i).second == NULL) continue;

		    QString fileContents;
		    //int curve_count = 0;

		    Contour_planes planes;
		    for(std::vector<SurfRecon::CurvePtr>::iterator cur = (*i).second->curves().begin();
			cur != (*i).second->curves().end();
			cur++)
		      {
			if(SurfRecon::getCurvePoints(**cur).empty()) continue;
			
			if(SurfRecon::getCurvePoints(**cur).size() < 4) //the last point duplicates the first
			  {
			    cvcapp.log(5, "Warning: Skipping curve with less than 3 points.");
			    continue;
			  }

			for(SurfRecon::PointPtrList::iterator pcur = SurfRecon::getCurvePoints(**cur).begin();
			    pcur != SurfRecon::getCurvePoints(**cur).end();
			    pcur++)
			  {
			    Contour_planes::iterator cur_plane;

			    //check to see if the plane for this curve exists already,
			    //if it does, add points and edges to that plane, else create a new plane and add to that
			    cur_plane = 
			      planes.insert(Contour_planes::value_type(Contour_basic_plane(0.0,0.0,
											   1.0,(*pcur)->z()),
								       Contour_plane(Contour_basic_plane(0.0,0.0,
													 1.0,(*pcur)->z()))))
			      .first;

			    Contour_vertex_set &curve_points = cur_plane->second.get<1>();
			    Contour_edges &curve_edges = cur_plane->second.get<2>();

			    std::pair<Contour_vertex_set::iterator,bool> p0 = 
			      curve_points.insert(Point_3((*pcur)->x(),
							  (*pcur)->y(),
							  (*pcur)->z()));
			    
			    std::pair<Contour_vertex_set::iterator,bool> p1;
			    //keep inserting edges until pcur gets to the duplicate point at the end of the curve
			    if(boost::next(pcur) != SurfRecon::getCurvePoints(**cur).end())
			      {
				p1 = curve_points.insert(Point_3((*boost::next(pcur))->x(),
								 (*boost::next(pcur))->y(),
								 (*boost::next(pcur))->z()));

				curve_edges.push_back(Contour_edge(p0.first,p1.first,1,2)); //left mat 1, right mat 2
			      }
			  }
		      }

#if 0
		    //plane description
		    fileContents += QString("%1 %2 %3 %4\n")
		      .arg(0.0)
		      .arg(0.0)
		      .arg(1.0)
		      .arg(curve_points[0].z());

		    //number of vertices and edges (same number in VolGridRover's data structure)
		    fileContents += QString("%1 %2\n")
		      .arg(curve_points.size())
		      .arg(curve_points.size());

		    for(std::vector<Point_3>::iterator cur_point = curve_points.begin();
			cur_point != curve_points.end();
			cur_point++)
		      {
			fileContents += QString("%1 %2 %3\n")
			  .arg(cur_point->x())
			  .arg(cur_point->y())
			  .arg(cur_point->z());
		      }

		    for(std::vector<Point_3>::iterator cur_point = curve_points.begin();
			cur_point != curve_points.end();
			cur_point++)
		      {
			fileContents += QString("%1 %2 %3 %4\n")
			  .arg(std::distance(curve_points.begin(),cur_point))
			  .arg(cur_point+1 == curve_points.end() ? 0 : 
			       std::distance(curve_points.begin(),cur_point)+1)
			  .arg(1) //not sure what material refers to...
			  .arg(2);
		      }
#endif
		    
		    //number of planes being written to file
		    fileContents = QString("%1\n").arg(planes.size());
		    //generate file contents with our magical plane map
		    for(Contour_planes::iterator cur_plane = planes.begin();
			cur_plane != planes.end();
			cur_plane++)
		      {
			//plane description
			fileContents += QString("%1 %2 %3 %4\n")
			  .arg(cur_plane->first.get<0>())  //a
			  .arg(cur_plane->first.get<1>())  //b
			  .arg(cur_plane->first.get<2>())  //c
			  .arg(cur_plane->first.get<3>()); //d

			//number of vertices and edges
			fileContents += QString("%1 %2\n")
			  .arg(cur_plane->second.get<1>().size())
			  .arg(cur_plane->second.get<2>().size());

			//vertices...
			for(Contour_vertex_set::iterator cur_point = cur_plane->second.get<1>().begin();
			    cur_point != cur_plane->second.get<1>().end();
			    cur_point++)
			  {
			    fileContents += QString("%1 %2 %3\n")
			      .arg(cur_point->x())
			      .arg(cur_point->y())
			      .arg(cur_point->z());
			  }

			//lets convert the set to a vector for efficiency of vertex index calculation
			Contour_vertex_vector vertices(cur_plane->second.get<1>().begin(),
						       cur_plane->second.get<1>().end());
			//edges
			for(Contour_edges::iterator cur_edge = cur_plane->second.get<2>().begin();

			    cur_edge != cur_plane->second.get<2>().end();
			    cur_edge++)
			  {
			    fileContents += QString("%1 %2 %3 %4\n")
			      .arg(std::lower_bound(vertices.begin(),
						    vertices.end(),
						    *(cur_edge->get<0>())) - vertices.begin())
			      .arg(std::lower_bound(vertices.begin(),
						    vertices.end(),
						    *(cur_edge->get<1>())) - vertices.begin())
			      .arg(cur_edge->get<2>())
			      .arg(cur_edge->get<3>());
			  }
		      }

		    QString new_filename(filename);
		    new_filename.replace(".contour",QString("_%1.contour").arg("%lf").arg( QString( (*i).second->name().data()) ) );
		    QFile f(new_filename);
		    if(!f.open(QIODevice::WriteOnly))
		      {
			QMessageBox::critical(this,"Error",QString("Could not open the file %1!").arg(new_filename),
					      QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
			return;
		      }

		    QTextStream stream(&f);
		    stream << fileContents;
		    f.close();
		  }
	      }
	}
    }
  else if(filename.endsWith(".bspline"))
    {
      fprintf( stderr, "Need to be ported\n");
/*
      int degree = 3;
      bool only_output_current_slice = false;
      //output each contour's point and the normal at that point instead of knots and control points
      bool points_and_normals = false;
      bool fit_splines = false;

      if(!m_Contours.empty())
	{
	  //create a dialog to get some info for building the b-spline files
	  {
	    bspline_opt bsd(this,0,true);
	    	    
	    if(bsd.exec() == QDialog::Accepted)
	      {
		degree = bsd.degree->text().toInt();
		only_output_current_slice = bsd.curslice->isChecked();
		points_and_normals = bsd.spline_output_type->selectedId() == 1;

		if(bsd.control_type->selectedId() == 0)
		  {
		
		 //   QMessageBox::critical(this,"Error","Fitting splines to contour points not yet implemented!",
		 //			  QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
		 //   return;

		    fit_splines = true;

		    if(degree != 3)
		      {
			QMessageBox::warning(this,
					     "Warning",
					     "Fitting splines is currently limited to using a degree of 3! Forcing...",
					     QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
			degree = 3;
		      }
		  }
	      }
	    else return;
	  }

	  unsigned int numSteps = 0;
	  
	  //calculate the number of steps
	  for(var = 0; var < m_VolumeFileInfo.numVariables(); var++)
	    for(time = 0; time < m_VolumeFileInfo.numTimesteps(); time++)
	      numSteps += m_Contours[var][time].size();

	  QProgressDialog progressDialog("Writing b-spline data for each objects set of contours...", 
					 "Cancel", numSteps, this, "Contours", true);
	  progressDialog.setProgress(0);

	  for(var = 0; var < m_VolumeFileInfo.numVariables(); var++)
	    for(time = 0; time < m_VolumeFileInfo.numTimesteps(); time++)
	      {
		for(SurfRecon::ContourPtrMap::iterator i = m_Contours[var][time].begin();
		    i != m_Contours[var][time].end();
		    i++)
		  {
		    progressDialog.setProgress(progressDialog.progress()+1);
		    qApp->processEvents();
		    if(progressDialog.wasCanceled())
		      return;

		    if((*i).second == NULL) continue;

		    //QString fileContents;
		    std::vector<double> vec_x, vec_y, vec_z, vec_t; //used to build the input arrays for gsl_spline_init()
		    double t, interval;
		    const gsl_interp_type *interp;
		    gsl_interp_accel *acc_x, *acc_y, *acc_z;
		    gsl_spline *spline_x, *spline_y, *spline_z;

		    switch((*i).second->interpolationType())
		      {
		      case 0: interp = gsl_interp_linear; break;
		      case 1: interp = gsl_interp_polynomial; break;
		      case 2: interp = gsl_interp_cspline; break;
		      case 3: interp = gsl_interp_cspline_periodic; break;
		      case 4: interp = gsl_interp_akima; break;
		      case 5: interp = gsl_interp_akima_periodic; break;
		      default: interp = gsl_interp_cspline; break;
		      }

		    for(std::vector<SurfRecon::CurvePtr>::iterator cur = (*i).second->curves().begin();
			cur != (*i).second->curves().end();
			cur++)
		      {
			if(SurfRecon::getCurvePoints(**cur).empty()) continue;
			if(only_output_current_slice)
			  {
			    //get the current shown canvas
			    SliceCanvas *cur_canvas = NULL;
			    switch(m_SliceCanvasStack->id(m_SliceCanvasStack->visibleWidget()))
			      {
			      default:
			      case 0: cur_canvas = m_XYSliceCanvas; break;
			      case 1: cur_canvas = m_XZSliceCanvas; break;
			      case 2: cur_canvas = m_ZYSliceCanvas; break;
			      }
			    
			    //now check if this curve is currently drawn on that canvas
			    if(SurfRecon::getCurveSlice(**cur) != cur_canvas->getDepth() ||
			       SurfRecon::getCurveOrientation(**cur) != SurfRecon::Orientation(cur_canvas->getSliceAxis()))
			      continue;
			  }

			t = 0.0;
			vec_x.clear(); vec_y.clear(); vec_z.clear(); vec_t.clear();
			

			for(SurfRecon::PointPtrList::iterator pcur = SurfRecon::getCurvePoints(**cur).begin();
			    pcur != SurfRecon::getCurvePoints(**cur).end();
			    pcur++)
			  {
			    vec_x.push_back((*pcur)->x());
			    vec_y.push_back((*pcur)->y());
			    vec_z.push_back((*pcur)->z());
			    vec_t.push_back(t); t+=1.0;

			    
			 //   fileContents += QString("%1 %2 %3\n")
			 //     .arg((*pcur)->x())
			 //     .arg((*pcur)->y())
			 //     .arg((*pcur)->z());
			    
			  }
			t-=1.0;
			//this is hackish but how else can i find the min size?
			gsl_spline *test_spline = gsl_spline_alloc(interp, 1000);
			if(vec_t.size() >= gsl_spline_min_size(test_spline))
			  {
			    double time_i;
			    SurfRecon::PointPtrList tmplist;
			    //int count;
			      
			    acc_x = gsl_interp_accel_alloc();
			    acc_y = gsl_interp_accel_alloc();
			    acc_z = gsl_interp_accel_alloc();
			    spline_x = gsl_spline_alloc(interp, vec_t.size());
			    spline_y = gsl_spline_alloc(interp, vec_t.size());
			    spline_z = gsl_spline_alloc(interp, vec_t.size());
			    gsl_spline_init(spline_x, &(vec_t[0]), &(vec_x[0]), vec_t.size());
			    gsl_spline_init(spline_y, &(vec_t[0]), &(vec_y[0]), vec_t.size());
			    gsl_spline_init(spline_z, &(vec_t[0]), &(vec_z[0]), vec_t.size());
			    
			    //Just sample some number of points for non-linear interpolation for now.
			    interval = 1.0/(1+(*i).second->numberOfSamples());
			    //In the future, sample according to the spline's curvature between 2
			    // control points...

			    //this bothers me... please just calculate the number of points in the future!
			    //for(time_i = 0.0, count=0; time_i<=t; time_i+=interval, count++);
			    //outfile << count << std::endl;
			    
			    for(time_i = 0.0; time_i<=t; time_i+=interval)
			      {
				//outfile << gsl_spline_eval(spline_x, time_i, acc_x) << " " 
				//	<< gsl_spline_eval(spline_y, time_i, acc_y) << " " 
				//	<< gsl_spline_eval(spline_z, time_i, acc_z) << std::endl;
				tmplist.push_back(SurfRecon::PointPtr(new SurfRecon::Point(gsl_spline_eval(spline_x, time_i, acc_x),
											   gsl_spline_eval(spline_y, time_i, acc_y),
											   gsl_spline_eval(spline_z, time_i, acc_z))));
			      }
			      
			    // clean up our mess
			    gsl_spline_free(spline_x);
			    gsl_spline_free(spline_y);
			    gsl_spline_free(spline_z);
			    gsl_interp_accel_free(acc_x);
			    gsl_interp_accel_free(acc_y);
			    gsl_interp_accel_free(acc_z);

			    SurfRecon::PointPtr pfront = tmplist.front();
			    SurfRecon::PointPtr pback = tmplist.back();

			    //remove the last point if the same as first
			    if((fabs(pfront->x() - pback->x()) <= 0.00001) &&
			       (fabs(pfront->y() - pback->y()) <= 0.00001) &&
			       (fabs(pfront->z() - pback->z()) <= 0.00001))
			      tmplist.pop_back();

			    //finally output the b-spline specification using the points we've collected as control points
			    //http://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-curve-closed.html

			    if(tmplist.size() < degree)
			      {
				cvcapp.log(5, "Warning: dropping contour because it doesn't have enough points");
			      }
			    else
			      {
				typedef CGAL::Simple_cartesian<double> Kernel;
				typedef Kernel::Point_2 Point;
				typedef Kernel::Line_2 Line;
				typedef Kernel::Vector_2 Vector;

				//write one file for each set of b-splines.. may change later
				QString fileContents;

				if(points_and_normals)
				  {
				    fileContents += QString("1\n");
				    fileContents += QString("%1\n").arg(tmplist.size());

				    //output points and averaged normals at each point on the contour
				    //We must have at least 2 points for this loop!
				    for(SurfRecon::PointPtrList::iterator k = tmplist.begin();
					k != tmplist.end();
					k++)
				      {
					using namespace boost;
					Line ante, des;

					//border cases and regular case
					if(k == tmplist.begin())
					  {
					    //the anterior edge is the last edge because this is a loop...
					    ante = Line(Point((*prior(tmplist.end()))->x(),
							      (*prior(tmplist.end()))->y()),
							Point((*k)->x(),
							      (*k)->y()));
					    des = Line(Point((*k)->x(),
							     (*k)->y()),
						       Point((*next(k))->x(),
							     (*next(k))->y()));
					  }
					else if(k == prior(tmplist.end()))
					  {
					    ante = Line(Point((*prior(k))->x(),
							      (*prior(k))->y()),
							Point((*k)->x(),
							      (*k)->y()));
					    des = Line(Point((*k)->x(),
							     (*k)->y()),
						       Point((*tmplist.begin())->x(),
							     (*tmplist.begin())->y()));
					  }
					else
					  {
					    ante = Line(Point((*prior(k))->x(),
							      (*prior(k))->y()),
							Point((*k)->x(),
							      (*k)->y()));
					    des = Line(Point((*k)->x(),
							     (*k)->y()),
						       Point((*next(k))->x(),
							     (*next(k))->y()));
					  }

					//normal to the vectors of the lines, counter clockwise
#if 0
					Vector ante_n(ante.to_vector().y(),
						      -ante.to_vector().x());
					Vector des_n(des.to_vector().y(),
						      -des.to_vector().x());
					Vector avg((ante_n.x()+des_n.x())/2,
						   (ante_n.y()+des_n.y())/2);
#endif
					Vector ante_n(ante.to_vector().y(),
						      -ante.to_vector().x());
					Vector des_n(des.to_vector().y(),
						      -des.to_vector().x());
					//make sure these are normalized
					ante_n = ante_n / (sqrt(ante_n.squared_length()));
					des_n = des_n / (sqrt(des_n.squared_length()));
					Vector avg = (ante_n + des_n)/2;
					
					//finally write out the string
					fileContents += 
					  QString("%1 %2 %3 %4\n")
					  .arg((*k)->x())
					  .arg((*k)->y())
					  .arg(avg.x())
					  .arg(avg.y());
				      }
				  }
				else
				  {
				    if(fit_splines)
				      {
#if 0
					using namespace boost::numeric::ublas;

					//duplicate the last point to close the loop...
					tmplist.push_back(*tmplist.begin());

					//http://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/INT-APP/CURVE-INT-global.html

					int n = tmplist.size()-1;
					int m = n + degree + 1;

					//uniform parameters
					std::vector<double> params;
					for(int k = 0; k <= n; k++)
					  params.push_back(double(k)/double(n));

					//construct uniform knot sequence of m+1 knots (clamped)
					//TODO: Probably should let the user choose different ways of constructing knot vector
					//beyond simply uniform spacing...
					std::vector<double> knots;
					for(int k = 0; k <= degree; k++) knots.push_back(0.0);
					for(int k = 1; k <= n-degree; k++) knots.push_back(double(k)/(double(n-degree+1)));
					for(int k = m - degree; k <= m; k++) knots.push_back(1.0);
					
					//compute basis function coefficients
					matrix<double> N(n+1,n+1);
					for(int k = 0; k < n+1; k++)
					  for(int j = 0; j < n+1; j++)
					    N(j,k) = 0.0;
					
					for(int k = 0; k < n+1; k++)
					  {
					    for(int j = 0; j < n+1; j++)
					      {
						
					      }
					  }
#endif

#if 0
					//since currently fit_splines returns a spline for each edge of the contour
					//instead of a single spline for the whole contour, we must output here each
					//spline separately for Lei Na's program!
					std::vector<SurfRecon::B_Spline> splines = SurfRecon::fit_spline(tmplist);
					fileContents += QString("%1\n").arg(splines.size());
					for(std::vector<SurfRecon::B_Spline>::iterator k = splines.begin();
					    k != splines.end();
					    k++)
					  {
					    fileContents += QString("%1\n").arg(degree); //should be 3 for now!
					    fileContents += QString("%1\n").arg(k->get<0>().size()); //should be 8
					    fileContents += QString("%1\n").arg(k->get<1>().size()); //should be 4
					    
					    for(std::vector<double>::iterator l = k->get<0>().begin();
						l != k->get<0>().end();
						l++)
					      fileContents += QString("%1\n").arg(*l);

					    for(std::vector<SurfRecon::Point_2>::iterator q = k->get<1>().begin();
						q != k->get<1>().end();
						q++)
					      fileContents += QString("%1 %2 %3\n")
						.arg(q->x())
						.arg(q->y())
						.arg((*tmplist.begin())->z()); //use the z value from the original points
					  }
#endif

					//combine fit spline's output into a single spline
					std::vector<SurfRecon::B_Spline> splines = SurfRecon::fit_spline(tmplist);
					SurfRecon::B_Spline new_spline;
					int n = splines.size() - 1;
					for(int k = -3; k <= n+4; k++)
					  SurfRecon::getKnots(new_spline).push_back(double(k));
					for(std::vector<SurfRecon::B_Spline>::iterator k = splines.begin();
					    k != splines.end();
					    k++)
					  {
					    std::vector<SurfRecon::Point_2> &pts = SurfRecon::getPoints(new_spline);
					    assert(SurfRecon::getPoints(*k).size() == 4);
					    pts.insert(pts.end(),
						       SurfRecon::getPoints(*k).begin(),
						       SurfRecon::getPoints(*k).end());
					    if(boost::next(k) != splines.end()) pts.erase(pts.end()-3,pts.end());
					  }

					fileContents += QString("1\n"); //number of splines
					fileContents += QString("%1\n").arg(degree);
					fileContents += QString("%1\n").arg(SurfRecon::getKnots(new_spline).size());
					fileContents += QString("%1\n").arg(SurfRecon::getPoints(new_spline).size());
					for(std::vector<double>::const_iterator k = SurfRecon::getKnots(new_spline).begin();
					    k != SurfRecon::getKnots(new_spline).end();
					    k++)
					  fileContents += QString("%1\n").arg(*k);
					for(std::vector<SurfRecon::Point_2>::const_iterator k = SurfRecon::getPoints(new_spline).begin();
					    k != SurfRecon::getPoints(new_spline).end();
					    k++)
					  fileContents += QString("%1 %2 %3\n")
					    .arg(k->x())
					    .arg(k->y())
					    .arg((*tmplist.begin())->z()); //use the z value from the original points
				      }
				    else
				      {
					//Since this should be a closed curve, lets wrap the control points
					for(SurfRecon::PointPtrList::iterator k = tmplist.begin();
					    std::distance(tmplist.begin(),k) < degree;
					    k++)
					  tmplist.push_back(*k);

					int n = tmplist.size()-1;
					int m = n + degree + 1;

					//construct uniform knot sequence of m+1 knots
					//TODO: Probably should let the user choose different ways of constructing knot vector
					//beyond simply uniform spacing...
					std::vector<double> knots;
					for(int k = 0; k <= m; k++)
					  knots.push_back(double(k)/double(m));

					
					//Since this should be a closed curve, lets wrap the control points
				//	for(SurfRecon::PointPtrList::iterator k = tmplist.begin();
				//	std::distance(tmplist.begin(),k) < degree;
				//	k++)
				//	tmplist.push_back(*k);
					

					fileContents += QString("1\n"); //number of splines
					fileContents += QString("%1\n").arg(degree);
					fileContents += QString("%1\n").arg(knots.size());
					fileContents += QString("%1\n").arg(tmplist.size());
				
					for(std::vector<double>::const_iterator kcur = knots.begin();
					    kcur != knots.end();
					    kcur++)
					  fileContents += QString("%1\n").arg(*kcur);

					for(SurfRecon::PointPtrList::const_iterator pcur = tmplist.begin();
					    pcur != tmplist.end();
					    pcur++)
					  fileContents += QString("%1 %2 %3\n")
					    .arg((*pcur)->x())
					    .arg((*pcur)->y())
					    .arg((*pcur)->z());
				      }
				  }

				if(!fileContents.isEmpty())
				  {
				    QString new_filename(filename);
				    new_filename.replace(".bspline",
							 QString("_%1_%2.bspline")
							 .arg((*i).second->name())
							 .arg(std::distance((*i).second->curves().begin(),
									    cur)));
				    QFile f(new_filename);
				    if(!f.open(QIODevice::WriteOnly))
				      {
					QMessageBox::critical(this,"Error",QString("Could not open the file %1!").arg(new_filename),
							      QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
					return;
				      }
				    
				    QTextStream stream(&f);
				    stream << fileContents;
				    f.close();
				  }
			      }
			  }
			gsl_spline_free(test_spline);
		      }

		  
		//    QString new_filename(filename);
		//    new_filename.replace(".bspline",QString("_%1.bspline").arg((*i).second->name()));
		//    QFile f(new_filename);
		//    if(!f.open(QIODevice::WriteOnly))
		//     {
		//	QMessageBox::critical(this,"Error",QString("Could not open the file %1!").arg(new_filename),
		//			      QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
		//	return;
		//      }

		//    QTextStream stream(&f);
		//    stream << fileContents;
		//    f.close();
		    
		  }
	      }
	}
      else
	{
	  QMessageBox::critical(this,"Error",QString("No contours!"),
				QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
	}
*/
    }
  else if(filename.endsWith(".matc")) // MAT-1.0.x Skeletonization Software contour format
    {
      fprintf( stderr, "need to be ported\n");
/*
      bool only_output_current_slice = true;

      if(!m_Contours.empty())
	{
	  //TODO: create a dialog to get some info for building the MAT contour files.
	  {
	    
	  }

	  unsigned int numSteps = 0;
	  
	  //calculate the number of steps
	  for(var = 0; var < m_VolumeFileInfo.numVariables(); var++)
	    for(time = 0; time < m_VolumeFileInfo.numTimesteps(); time++)
	      numSteps += m_Contours[var][time].size();

	  QProgressDialog progressDialog("Writing contour data in MAT's format for each objects set of contours...", 
					 "Cancel", numSteps, this, "Contours", true);
	  progressDialog.setProgress(0);

	  QString fileContents;

	  for(var = 0; var < m_VolumeFileInfo.numVariables(); var++)
	    for(time = 0; time < m_VolumeFileInfo.numTimesteps(); time++)
	      {
		for(SurfRecon::ContourPtrMap::iterator i = m_Contours[var][time].begin();
		    i != m_Contours[var][time].end();
		    i++)
		  {
		    progressDialog.setProgress(progressDialog.progress()+1);
		    qApp->processEvents();
		    if(progressDialog.wasCanceled())
		      return;

		    if((*i).second == NULL) continue;

		    //QString fileContents;
		    std::vector<double> vec_x, vec_y, vec_z, vec_t; //used to build the input arrays for gsl_spline_init()
		    double t, interval;
		    const gsl_interp_type *interp;
		    gsl_interp_accel *acc_x, *acc_y, *acc_z;
		    gsl_spline *spline_x, *spline_y, *spline_z;

		    switch((*i).second->interpolationType())
		      {
		      case 0: interp = gsl_interp_linear; break;
		      case 1: interp = gsl_interp_polynomial; break;
		      case 2: interp = gsl_interp_cspline; break;
		      case 3: interp = gsl_interp_cspline_periodic; break;
		      case 4: interp = gsl_interp_akima; break;
		      case 5: interp = gsl_interp_akima_periodic; break;
		      default: interp = gsl_interp_cspline; break;
		      }

		    for(std::vector<SurfRecon::CurvePtr>::iterator cur = (*i).second->curves().begin();
			cur != (*i).second->curves().end();
			cur++)
		      {
			if(SurfRecon::getCurvePoints(**cur).empty()) continue;
			if(only_output_current_slice)
			  {
			    //get the current shown canvas
			    SliceCanvas *cur_canvas = NULL;
			    switch(m_SliceCanvasStack->id(m_SliceCanvasStack->visibleWidget()))
			      {
			      default:
			      case 0: cur_canvas = m_XYSliceCanvas; break;
			      case 1: cur_canvas = m_XZSliceCanvas; break;
			      case 2: cur_canvas = m_ZYSliceCanvas; break;
			      }
			    
			    //now check if this curve is currently drawn on that canvas
			    if(SurfRecon::getCurveSlice(**cur) != cur_canvas->getDepth() ||
			       SurfRecon::getCurveOrientation(**cur) != SurfRecon::Orientation(cur_canvas->getSliceAxis()))
			      continue;
			  }

			t = 0.0;
			vec_x.clear(); vec_y.clear(); vec_z.clear(); vec_t.clear();
			
			QString object_name = QString("%1___%2_%3")
			  .arg(getCurveName(**cur))
			  .arg(std::distance((*i).second->curves().begin(),cur))
			  .arg(getCurveSlice(**cur));

			for(SurfRecon::PointPtrList::iterator pcur = SurfRecon::getCurvePoints(**cur).begin();
			    pcur != SurfRecon::getCurvePoints(**cur).end();
			    pcur++)
			  {
			    vec_x.push_back((*pcur)->x());
			    vec_y.push_back((*pcur)->y());
			    vec_z.push_back((*pcur)->z());
			    vec_t.push_back(t); t+=1.0;

			    
			 //   fileContents += QString("%1 %2 %3\n")
			 //     .arg((*pcur)->x())
			 //     .arg((*pcur)->y())
			 //     .arg((*pcur)->z());
			    
			  }
			t-=1.0;
			//this is hackish but how else can i find the min size?
			gsl_spline *test_spline = gsl_spline_alloc(interp, 1000);
			if(vec_t.size() >= gsl_spline_min_size(test_spline))
			  {
			    double time_i;
			    SurfRecon::PointPtrVector tmplist;
			    //int count;
			      
			    acc_x = gsl_interp_accel_alloc();
			    acc_y = gsl_interp_accel_alloc();
			    acc_z = gsl_interp_accel_alloc();
			    spline_x = gsl_spline_alloc(interp, vec_t.size());
			    spline_y = gsl_spline_alloc(interp, vec_t.size());
			    spline_z = gsl_spline_alloc(interp, vec_t.size());
			    gsl_spline_init(spline_x, &(vec_t[0]), &(vec_x[0]), vec_t.size());
			    gsl_spline_init(spline_y, &(vec_t[0]), &(vec_y[0]), vec_t.size());
			    gsl_spline_init(spline_z, &(vec_t[0]), &(vec_z[0]), vec_t.size());
			    
			    //Just sample some number of points for non-linear interpolation for now.
			    interval = 1.0/(1+(*i).second->numberOfSamples());
			    //In the future, sample according to the spline's curvature between 2
			    // control points...

			    //this bothers me... please just calculate the number of points in the future!
			    //for(time_i = 0.0, count=0; time_i<=t; time_i+=interval, count++);
			    //outfile << count << std::endl;
			    
			    for(time_i = 0.0; time_i<=t; time_i+=interval)
			      {
				//outfile << gsl_spline_eval(spline_x, time_i, acc_x) << " " 
				//	<< gsl_spline_eval(spline_y, time_i, acc_y) << " " 
				//	<< gsl_spline_eval(spline_z, time_i, acc_z) << std::endl;
				tmplist.push_back(SurfRecon::PointPtr(new SurfRecon::Point(gsl_spline_eval(spline_x, time_i, acc_x),
											   gsl_spline_eval(spline_y, time_i, acc_y),
											   gsl_spline_eval(spline_z, time_i, acc_z))));
			      }
			      
			    // clean up our mess
			    gsl_spline_free(spline_x);
			    gsl_spline_free(spline_y);
			    gsl_spline_free(spline_z);
			    gsl_interp_accel_free(acc_x);
			    gsl_interp_accel_free(acc_y);
			    gsl_interp_accel_free(acc_z);

			    SurfRecon::PointPtr pfront = tmplist.front();
			    SurfRecon::PointPtr pback = tmplist.back();

			    //remove the last point if the same as first
			    if((fabs(pfront->x() - pback->x()) <= 0.00001) &&
			       (fabs(pfront->y() - pback->y()) <= 0.00001) &&
			       (fabs(pfront->z() - pback->z()) <= 0.00001))
			      tmplist.pop_back();

			    fileContents += QString("object %1 {\n"
						    "type = contour_unord; format = ascii;\n"
						    "contour_id = %2;\n"
						    "closed = true;\n"
						    "pointdim = %3;\n"
						    "x_extent = 0; y_extent = 0;\n"
						    "beginloadseq (Sy_GenLoader)\n"
						    "#2718sync314beg#\n")
			      .arg(object_name)
			      .arg(std::distance((*i).second->curves().begin(),cur))
			      .arg(SurfRecon::getCurvePoints(**cur).size());
			    //fileContents += "\n";
			    
			    //get the current shown canvas
			    SliceCanvas *cur_canvas = NULL;
			    switch(m_SliceCanvasStack->id(m_SliceCanvasStack->visibleWidget()))
			      {
			      default:
			      case 0: cur_canvas = m_XYSliceCanvas; break;
			      case 1: cur_canvas = m_XZSliceCanvas; break;
			      case 2: cur_canvas = m_ZYSliceCanvas; break;
			      }
			    
			    //get the current shown slice's resolution
			    unsigned int xres, yres;
			    switch(cur_canvas->getSliceAxis())
			      {
			      case SliceCanvas::XY:
				xres = m_VolumeFileInfo.XDim();
				yres = m_VolumeFileInfo.YDim();
				break;
			      case SliceCanvas::XZ:
				xres = m_VolumeFileInfo.XDim();
				yres = m_VolumeFileInfo.ZDim();
			      case SliceCanvas::ZY:
				xres = m_VolumeFileInfo.ZDim();
				yres = m_VolumeFileInfo.YDim();
			      }

			    for(SurfRecon::PointPtrVector::iterator k = tmplist.begin();
				k != tmplist.end();
				k++)
			      {
				if((k-tmplist.begin())%10 == 0) fileContents += "\n";

				//TODO: let the user select what resolution they want to map to
				//For now, just use the slice resolution for the current shown slice

				switch(cur_canvas->getSliceAxis())
				  {
				  case SliceCanvas::XY:
				    fileContents += QString("%1 %2 ")
				      .arg(int((((*k)->x()-m_VolumeFileInfo.XMin())/
						(m_VolumeFileInfo.XMax()-m_VolumeFileInfo.XMin()))*double(xres)))
				      .arg(int((((*k)->y()-m_VolumeFileInfo.YMin())/
						(m_VolumeFileInfo.YMax()-m_VolumeFileInfo.YMin()))*double(yres)));
				    break;
				  case SliceCanvas::XZ:
				    fileContents += QString("%1 %2 ")
				      .arg(int((((*k)->x()-m_VolumeFileInfo.XMin())/
						(m_VolumeFileInfo.XMax()-m_VolumeFileInfo.XMin()))*double(xres)))
				      .arg(int((((*k)->z()-m_VolumeFileInfo.ZMin())/
						(m_VolumeFileInfo.ZMax()-m_VolumeFileInfo.ZMin()))*double(yres)));
				    break;
				  case SliceCanvas::ZY:
				    fileContents += QString("%1 %2 ")
				      .arg(int((((*k)->z()-m_VolumeFileInfo.ZMin())/
						(m_VolumeFileInfo.ZMax()-m_VolumeFileInfo.ZMin()))*double(xres)))
				      .arg(int((((*k)->y()-m_VolumeFileInfo.YMin())/
						(m_VolumeFileInfo.YMax()-m_VolumeFileInfo.YMin()))*double(yres)));
				    break;
				  }
			      }
			    fileContents += "\n\n";
			    fileContents += QString("#2718sync314end#\n"
						    "endloadseq;\n"
						    "}\n\n");
			    
			    
			    //  if(!fileContents.isEmpty())
			    //  {
			    //  QString new_filename(filename);
			    //  new_filename.replace(".bspline",
			    //  QString("_%1_%2.bspline")
			    //  .arg((*i).second->name())
			    //  .arg(std::distance((*i).second->curves().begin(),
			    //  cur)));
			    //  QFile f(new_filename);
			    //  if(!f.open(QIODevice::WriteOnly))
			    //  {
			    //  QMessageBox::critical(this,"Error",QString("Could not open the file %1!").arg(new_filename),
			    //  QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
			    //  return;
			    //  }
				
			    //  QTextStream stream(&f);
			    //  stream << fileContents;
			    //  f.close();
			    //  }
			   
			  }
			gsl_spline_free(test_spline);
		      }

		   
		  //  QString new_filename(filename);
		  //  new_filename.replace(".bspline",QString("_%1.bspline").arg((*i).second->name()));
		  //  QFile f(new_filename);
		  //  if(!f.open(QIODevice::WriteOnly))
		  //   {
		  //	QMessageBox::critical(this,"Error",QString("Could not open the file %1!").arg(new_filename),
		  //			      QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
		  //	return;
		  //    }

		  // QTextStream stream(&f);
		  //  stream << fileContents;
		  //  f.close();
		    
		  }
	      }

	  QFile f(filename);
	  if(!f.open(QIODevice::WriteOnly))
	    {
	      QMessageBox::critical(this,"Error",QString("Could not open the file %1!").arg(filename),
				    QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
	      return;
	    }
	  
	  QTextStream stream(&f);
	  stream << fileContents;
	  f.close();
	}
      else
	{
	  QMessageBox::critical(this,"Error",QString("No contours!"),
				QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
	}
*/
    }
  //if(filename.endsWith(".cnt"))
  else //default with cnt format
    {
      fprintf( stderr, "need to be ported\n");
/*
      QDomDocument doc("contourdoc");
      QDomElement root = doc.createElement("contours");
      doc.appendChild(root);
      
      if(!m_Contours.empty())
	{
	  unsigned int numSteps=0;
	  
	  //calculate the number of steps
	  for(var = 0; var < m_VolumeFileInfo.numVariables(); var++)
	    for(time = 0; time < m_VolumeFileInfo.numTimesteps(); time++)
	      numSteps += m_Contours[var][time].size();
	  
	  QProgressDialog progressDialog("Creating a DOM of contours for XML output...", "Cancel", numSteps, this, "Contours", true);
	  progressDialog.setProgress(0);
	  
	  for(var = 0; var < m_VolumeFileInfo.numVariables(); var++)
	    for(time = 0; time < m_VolumeFileInfo.numTimesteps(); time++)
	      {
		for(SurfRecon::ContourPtrMap::iterator i = m_Contours[var][time].begin();
		    i != m_Contours[var][time].end();
		    i++)
		  {
		    progressDialog.setProgress(progressDialog.progress()+1);
		    qApp->processEvents();
		    if(progressDialog.wasCanceled())
		      return;
		    
		    if((*i).second == NULL) continue;
		    
		    QDomElement ctag;
		    
		    ctag = doc.createElement("contour");
		    ctag.setAttribute("name", (*i).second->name());
		    ctag.setAttribute("color", QColor(int(clamp((*i).second->color().get<0>())*255.0),
						      int(clamp((*i).second->color().get<1>())*255.0),
						      int(clamp((*i).second->color().get<2>())*255.0)).name());
		    ctag.setAttribute("variable", QString("%1").arg(var));
		    ctag.setAttribute("timestep", QString("%1").arg(time));
		    ctag.setAttribute("interpolationType", 
				      SurfRecon::InterpolationTypeStrings[(*i).second->interpolationType()]);
		    ctag.setAttribute("numSamples", QString("%1").arg((*i).second->numberOfSamples()));
		    root.appendChild(ctag);
		    
		    for(std::vector<SurfRecon::CurvePtr>::iterator cur = (*i).second->curves().begin();
			cur != (*i).second->curves().end();
			cur++)
		      {
			QDomElement curveTag;
			
			if(SurfRecon::getCurvePoints(**cur).empty()) continue;
			
			curveTag = doc.createElement("curve");
			curveTag.setAttribute("slice", QString("%1").arg(SurfRecon::getCurveSlice(**cur)));
			curveTag.setAttribute("orientation", SurfRecon::getCurveOrientationString(**cur));
			curveTag.setAttribute("name", SurfRecon::getCurveName(**cur));
			ctag.appendChild(curveTag);
			
			for(SurfRecon::PointPtrList::iterator pcur = SurfRecon::getCurvePoints(**cur).begin();
			    pcur != SurfRecon::getCurvePoints(**cur).end();
			    pcur++)
			  {
			    QDomElement ptag;
			    
			    //if(*pcur != NULL)
			    {
			      ptag = doc.createElement("point");
			      ptag.setAttribute("x", QString("%1").arg((*pcur)->x()));
			      ptag.setAttribute("y", QString("%1").arg((*pcur)->y()));
			      ptag.setAttribute("z", QString("%1").arg((*pcur)->z()));
			      curveTag.appendChild(ptag);
			    }
			  }
		      }
		  }
	      }
	  
	  progressDialog.setProgress(numSteps);
	}
      
      QFile f(filename);
      if(!f.open(QIODevice::WriteOnly))
	{
	  QMessageBox::critical(this,"Error",QString("Could not open the file %1!").arg(filename),
				QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
	  return;
	}
      
      QTextStream stream(&f);
      stream << doc.toString();
      f.close();
*/
    }
}

void VolumeGridRover::loadContoursSlot()
{
  printf(" in load Contours slot\n");

  // std::list<SurfRecon::ContourPtr> contours;

  QStringList filenames = QFileDialog::getOpenFileNames(this,
							"Choose a file to open",
							QString(),
							"Reconstruct series (*.ser);;"
							"Contour files (*.cnt);;"
							"ContourTiler config (*.config);;"
							"All Files (*)");

  if(filenames.size() == 0)
    {
      cvcapp.log(5, "VolumeGridRover::loadContoursSlot(): load cancelled (or null filename)");
      return;
    }

  if(m_Contours.empty()) m_Contours.resize(1);
  if(m_Contours[0].empty()) m_Contours[0].resize(1);

  updateSliceContours(); //make sure the slices have a valid m_Contours

  //m_Contours[m_Variable->currentIndex()][m_Timestep->value()].clear();
  //m_Contour->clear();

  for(QStringList::iterator it = filenames.begin();
      it != filenames.end();
      it++)
    {
      QString &filename = *it;
      if(filename.endsWith(".ser"))
	{
#ifdef USING_TILING

	 // //  * Series XML file format:
	 // //  *
	 // //  * In filename <name>.ser:
	 // //  *
	 // //  *  <Series> - global series attributes, similar to <Section>, has index attrib for child 2D contours to have a Z val
	 // //  *    <Contour/> - contour description, point attrib has point coords.  Each point is 2D.
	 // //  *           name - name of object that this contour belongs to
	 // //  *           closed - whether contour is closed or open (bool)
	 // //  *           border - border color I think
	 // //  *           fill - fill color I think (ignoring for now since we use same color for fill and border
	 // //  *           mode - drawing mode? ignore for now.
	 // //  *    <ZContour/> - contour description with 3D points, same attribs as Contour
	 // //  *  </Series>
	 // //  *
	 // //  * In section file <name>.<section index>
	 // //  *
	 // //  *  <Section> - set of contours for each slice. index attribute
	 // //  *    <Transform> - transform applied to children elements
	 // //  *      <Image/>  - An image for this section (optional)
	 // //  *      <Contour/> - same as above
	 // //  *    </Transform>
	 // //  *  </Section>
	 // //  *
	   

	 //  contours = SeriesFileReader::readSeries(filename.toStdString().c_str(),m_VolumeFileInfo);//1.0//,1000.0);
	 //  //cvcapp.log(5, "contours.size(): %d", contours.size());

	 //  for(std::list<SurfRecon::ContourPtr>::iterator i = contours.begin();
	 //      i != contours.end();
	 //      i++)
	 //    {
	 //      m_Contours[m_Variable->currentIndex()][m_Timestep->value()][(*i)->name()] = *i;
	 //      //m_Contour->insertItem((*i)->name());
	 //      QListViewItem *tmp = m_Objects->findItem((*i)->name().c_str(),0);
	 //      if(tmp) delete tmp;
	 //      m_Objects->insertItem(new QListViewItem(m_Objects,(*i)->name().c_str()));
	 //    }
	  list<Contour_handle> contours;
	  list<Contour_exception> exceptions;
	  vector<string> empty_components;
	  const int start = 59;
	  const int end = 160;
	  const double smoothing_factor = -1;
	  list<string> components;
	  list<string> components_skip;
	  string fn = filename.toStdString();
	  fn = fn.substr(0, fn.length() - 4);
	  read_contours_ser(fn, back_inserter(contours), back_inserter(exceptions),
			    start, end, smoothing_factor,
			    components.begin(), components.end(),
			    components_skip.begin(), components_skip.end());

	  for (list<Contour_exception>::const_iterator it = exceptions.begin(); it != exceptions.end(); ++it) {
	    std::cout << "Error in reading contours: " << it->what();
	    // LOG4CPLUS_WARN(logger, "Error in reading contours: " << it->what());
	  }
	  set<string> component_names;
	  for (list<Contour_handle>::iterator it = contours.begin(); it != contours.end(); ++it) {
	    const Contour_handle contour = *it;
	    string component = contour->info().object_name();
	    component_names.insert(component);
	  }

	  for (set<string>::const_iterator it = component_names.begin(); it != component_names.end(); ++it) {
	    string name = *it;
	    // m_Contours[m_Variable->currentIndex()][m_Timestep->value()][name] = *i;
	    // QListWidgetItem *tmp = _ui->m_Objects->findItem(name.c_str(),0);
	    // if(tmp) delete tmp;
	    // _ui->m_Objects->insertItem(new QListWidgetItem(QString(name.c_str()), _ui->m_Objects, 0));
	    // _ui->m_Objects->addItem(QString(name.c_str()));
	    // _ui->m_tileObjectsList->addItem(QString(name.c_str()));
	  }
#else
	  QMessageBox::critical(this,"Error","Build with Tiling support to load Reconstruct Series files",
				QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
	  //return;
	  continue;
#endif
	}
      else if(filename.endsWith(".config"))
	{
	  using namespace std;
	  using namespace boost;
	  using namespace boost::algorithm;
     
	  //ContourTiler's pts/config format
	  try
	    {
	      int linenum = 0;
	      string config_file_contents;
	      ifstream config_file(filename.toStdString().c_str());
	      if(!config_file)
		throw runtime_error(str(format("Could not open file %1%") % filename.toStdString().c_str()));
	  
	      while(!config_file.eof())
		{
		  string line;
		  getline(config_file,line);
		  linenum++;
		  if(!config_file && !config_file.eof())
		    throw runtime_error(str(format("Could not read file %1% at line %2%")
					    % filename.toStdString().c_str()
					    % linenum));
		  //remove comments
		  line.replace(std::find(line.begin(), line.end(), '#'),
			       line.end(),
			       string(""));
		  //remove extra whitespace
		  trim(line);
		  if(line.empty()) continue;
		  config_file_contents += line + "\n";
		}

	      //cout << "Config file: " << endl;
	      //cout << config_file_contents << endl;

	      //search and collect iterators for the start of each section
	      regex expression("^(.*):(.*)$");
	      cmatch what;
	      stringstream config_file_stream(config_file_contents);
	      set<string::iterator> section_tokens;
	      while(!config_file_stream.eof())
		{
		  string line;
		  getline(config_file_stream, line);
		  if(regex_match(line.c_str(),what,expression))
		    section_tokens.insert(config_file_contents.begin() + 
					  config_file_contents.find(what[1]));
		}

	      //insert the end iterator
	      section_tokens.insert(config_file_contents.end());

	      if(section_tokens.size() <= 1)
		throw runtime_error(str(format("Invalid config file %1%") % filename.toStdString().c_str()));
	  
	      //iterate through each section and collect the input
	      bool seen_prefix = false, seen_suffix = false, seen_range = false;
	      string prefix, suffix;
	      int low_range=0, high_range=0;
	      for(set<string::iterator>::iterator i = section_tokens.begin();
		  i != prior(section_tokens.end());
		  i++)
		{
		  //find where the colon is and move 1 past it
		  string::iterator colon = std::find(*i,*std::next(i),':');
		  assert(colon != *std::next(i)); //this must be the case because of the regex above
		  string param(*i,colon);
		  string arg(std::next(colon),*std::next(i));
		  trim(arg);
		  //cout << "lolz: '" << arg << "'" << endl;

		  //we only care about the following
		  if(param == "PREFIX")
		    {
		      prefix = arg;
		      seen_prefix = true;
		    }
		  else if(param == "SUFFIX")
		    {
		      suffix = arg;
		      seen_suffix = true;
		    }
		  else if(param == "SLICE_RANGE")
		    {
		      vector<string> split_arg;
		      split(split_arg,arg,is_any_of(" "));
		      if(split_arg.size() != 2)
			throw runtime_error(str(format("Invalid config file %1% - bad slice range")
						% filename.toStdString().c_str()));
		      try
			{
			  low_range = lexical_cast<int>(split_arg[0]);
			  high_range = lexical_cast<int>(split_arg[1]);
			}
		      catch(bad_lexical_cast& e)
			{
			  throw runtime_error(str(format("Invalid config file %1% - %2%")
						  % filename.toStdString().c_str()
						  % e.what()));
			}
		      seen_range = true;
		    }
		}

	      if(!seen_prefix || !seen_suffix || !seen_range)
		throw runtime_error(str(format("Invalid config file %1% - must have at least prefix, suffix, and range")
					% filename.toStdString().c_str()));
	  
	      cout << "Prefix: " << prefix << endl;
	      cout << "Suffix: " << suffix << endl;
	      cout << "Slice range: " << low_range << " " << high_range << endl;

	      int var = _ui->m_Variable->currentIndex();
	      int time = _ui->m_Timestep->value();

	      m_Contours[var][time][prefix] =
		SurfRecon::ContourPtr(new SurfRecon::Contour(SurfRecon::Color(0.5,0.5,0.5),
							     prefix));
      
	      //iterate through each of the pts files, adding them to the contour data structure
	      for(int i = low_range; i <= high_range; i++)
		{
		  string pts_filename(str(format("%1%%2%%3%") % 
					  prefix % 
					  i %
					  suffix));

		  if(QFileInfo(pts_filename.c_str()).isRelative())
		    {
		      pts_filename = string((QFileInfo(filename)
					     .absolutePath() + "/" + pts_filename.c_str()).toStdString().c_str());
		      //cout << "abs filename: " << pts_filename << endl;
		    }

		  ifstream ptsfile(pts_filename.c_str());
		  if(!ptsfile)
		    throw runtime_error(str(format("Could not open file %1%") % pts_filename));
	      
		  linenum = 0;
		  try
		    {
		      while(!ptsfile.eof() && !ptsfile.fail())
			{
			  string line;
			  vector<string> line_split;
			  int numpts;
			  getline(ptsfile,line);
			  linenum++;
			  trim(line);
			  if(line.empty()) //if there is a blank line, skip it
			    continue;
			  if(!ptsfile)
			    throw runtime_error(str(format("Error reading number of points in file %1% line %2%")
						    % pts_filename
						    % linenum));
			  split(line_split,line,is_any_of(" "));
			  if(line_split.size() != 1)
			    throw runtime_error(str(format("Invalid number of points in file %1% line %2%")
						    % pts_filename
						    % linenum));
			  try
			    {
			      numpts = lexical_cast<int>(line_split[0]);
			    }
			  catch(bad_lexical_cast& e)
			    {
			      throw runtime_error(str(format("Error converting string in file %1% line %2% to integer (%3%)")
						      % pts_filename
						      % linenum
						      % line.size()));
			    }
		  
			  SurfRecon::CurvePtr curveptr(new SurfRecon::Curve(i,
									    SurfRecon::XY,
									    SurfRecon::PointPtrList(),
									    prefix));
			  m_Contours[var][time][prefix]->add(curveptr);
			  for(int j = 0; j < numpts; j++)
			    {
			      double x,y,z;
			      getline(ptsfile,line);
			      linenum++;
			      trim(line);
			      if(!ptsfile)
				throw runtime_error(str(format("Error reading point %1% in file %2% line %3%")
							% (j+1)
							% pts_filename
							% linenum));
			      if(line.empty()) //if there is a blank line, skip it
				{
				  j--; //don't count this line as a point
				  continue;
				}
			      split(line_split,line,is_any_of(" "));
			      if(line_split.size() != 3)
				throw runtime_error(str(format("Invalid point in file %1% line %2%")
							% pts_filename
							% linenum));
			      x = lexical_cast<double>(line_split[0]);
			      y = lexical_cast<double>(line_split[1]);
			      z = lexical_cast<double>(line_split[2]);
			      SurfRecon::getCurvePoints(*curveptr)
				.push_back(SurfRecon::PointPtr(new SurfRecon::Point(x,y,z)));
			    }

			  //duplicate the first point and put it in the end to make this a closed loop
			  if(SurfRecon::getCurvePoints(*curveptr).size() >= 3 &&
			     **SurfRecon::getCurvePoints(*curveptr).begin() != 
			     **prior(SurfRecon::getCurvePoints(*curveptr).end()))
			    SurfRecon::getCurvePoints(*curveptr)
			      .push_back(SurfRecon::PointPtr(new SurfRecon::Point(**SurfRecon::getCurvePoints(*curveptr).begin())));
			}
		    }
		  catch(bad_lexical_cast& e)
		    {
		      throw runtime_error(str(format("Error converting string to number in file %1% line %2% (%3%)")
					      % pts_filename
					      % linenum
					      % e.what()));
		    }
		}
	
	      if( _ui->m_Objects->findItems( QString(prefix.data()) ,Qt::MatchExactly).size() == 0)
		_ui->m_Objects->insertItem(_ui->m_Objects->count(), new QListWidgetItem( QString( prefix.data() ), _ui->m_Objects, 0 ));

	      _ui->m_Objects->setCurrentItem( _ui->m_Objects->findItems( QString(prefix.data()), Qt::MatchExactly ).first() );
	    }
	  catch(std::exception& e)
	    {
	      QMessageBox::critical(this,"Error",
				    QString(e.what()),
				    QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
	      //return;
	      continue;
	    }
	}
      else
	{
          printf("not implemented yet\n");
/*
	  //default with cnt file
	  QDomDocument doc("contourdoc");
	  QFile file(filename);
	  if(!file.open(QIODevice::ReadOnly))
	    {
	      QMessageBox::critical(this,"Error","Unable to open contour data file!",
				    QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
	      //return;
	      continue;
	    }
	  if(!doc.setContent(&file))
	    {
	      QMessageBox::critical(this,"Error","Unable to open contour data file!",
				    QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
	      file.close();
	      //return;
	      continue;
	    }
	  file.close();

	  QDomElement contours = doc.documentElement();
	  //if(!contours.isNull()) cvcapp.log(5, "contours.tagName() == %s",contours.tagName().toStdString().c_str());
	  unsigned int numCurves = contours.elementsByTagName("curve").count();

	  QProgressDialog progressDialog( QString("Building contours from XML input..."), QString("Cancel"), 0, numCurves, this);
	  progressDialog.setValue(0);

	  if(!contours.isNull() && contours.tagName() == "contours")
	    {
	      QDomNode contourNode = contours.firstChild();
	      while(!contourNode.isNull())
		{
		  QDomElement contour = contourNode.toElement();
		  if(!contour.isNull() && contour.tagName() == "contour")
		    {
		      //cvcapp.log(5, "contour.tagName() == %s", contour.tagName().toStdString().c_str());
		  
		      unsigned int time = contour.attribute("timestep").toInt();
		      unsigned int var = contour.attribute("variable").toInt();
		      QString name = contour.attribute("name");
		  
		      if(name.isEmpty() || 
			 ((time >= m_VolumeFileInfo.numTimesteps() ||
			   var >= m_VolumeFileInfo.numVariables()) && m_VolumeFileInfo.isSet()))
			{
			  contourNode = contourNode.nextSibling();
			  continue;
			}
		      else if(!m_VolumeFileInfo.isSet()) //if we don't have a volume loaded, use 0
			{
			  var = time = 0;
			}
		  
		      QString colorstr = contour.attribute("color");
		      QColor color(colorstr.isEmpty() ? "#FFFFFF" : colorstr);
		      QString sampstring = contour.attribute("numSamples");
		      int samples = sampstring.isEmpty() ? 5 : sampstring.toInt();
		  
		      m_Contours[var][time][name.toStdString().c_str()] =
			SurfRecon::ContourPtr(new SurfRecon::Contour(SurfRecon::Color(color.red()/255.0,
										      color.green()/255.0,
										      color.blue()/255.0),
								     std::string(name.toStdString().c_str()),
								     SurfRecon::getInterpolationTypeFromString( std::string(contour.attribute("interpolationType").toStdString().c_str()) ),
								     samples));
		  
		      if( _ui->m_Objects->findItems( name,Qt::MatchExactly).size() == 0 )
 			 _ui->m_Objects->insertItem( _ui->m_Objects->count(), new QListWidgetItem(name, _ui->m_Objects, 0) );
		      _ui->m_Objects->setCurrentItem( _ui->m_Objects->findItems( name, 0 ).first() );
		  
		      QDomNode curveNode = contour.firstChild();
		      while(!curveNode.isNull())
			{
			  QDomElement curve = curveNode.toElement();
			  if(!curve.isNull() && curve.tagName() == "curve")
			    {
			      progressDialog.setValue(progressDialog.value()+1);
			      qApp->processEvents();
			      if(progressDialog.wasCanceled())
				goto end; //we should still use what we already processed

			      QString slicestr = curve.attribute("slice");
			      QString orientationstr = curve.attribute("orientation");
			      QString curvenamestr = curve.attribute("name");

			      //cvcapp.log(5, "curve.tagName()=%s, slice=%s, orientation=%s, curvenamestr=%s",curve.tagName().toStdString().c_str(),
			      //	 slicestr.toStdString().c_str(),orientationstr.toStdString().c_str(),curvenamestr.toStdString().c_str());

			      if(slicestr.isEmpty() || orientationstr.isEmpty())
				{
				  curveNode = curveNode.nextSibling();
				  continue;
				}
			  
			      int slice = slicestr.toInt();
			      SurfRecon::Orientation orientation = 
				orientationstr == "XY" ? SurfRecon::XY :
				orientationstr == "XZ" ? SurfRecon::XZ :
				orientationstr == "ZY" ? SurfRecon::ZY : SurfRecon::XY;

			      SurfRecon::CurvePtr curveptr(new SurfRecon::Curve(slice,
										orientation,
										SurfRecon::PointPtrList(),
										curvenamestr));
			  
			      m_Contours[var][time][name.toStdString().c_str()]->add(curveptr);
			  
			      QDomNode pointNode = curve.firstChild();
			      while(!pointNode.isNull())
				{
				  QDomElement point = pointNode.toElement();
				  if(!point.isNull() && point.tagName() == "point")
				    {
				      QString xstr = point.attribute("x");
				      QString ystr = point.attribute("y");
				      QString zstr = point.attribute("z");
				  
				      if(xstr.isEmpty() || ystr.isEmpty() || zstr.isEmpty())
					{
					  pointNode = pointNode.nextSibling();
					  continue;
					}
				  
				      double x = xstr.toDouble();
				      double y = ystr.toDouble();
				      double z = zstr.toDouble();
				  
				      SurfRecon::getCurvePoints(*curveptr).push_back(SurfRecon::PointPtr(new SurfRecon::Point(x,y,z)));
				    }
				  pointNode = pointNode.nextSibling();
				}
			    }
			  curveNode = curveNode.nextSibling();
			}
		    }
		  contourNode = contourNode.nextSibling();
		}
	    }

	  progressDialog.setValue(numCurves);
*/
	}
    }

  //if a volume isn't loaded, lets generate one
  if(!m_VolumeFileInfo.isSet())
    {
      using namespace boost;
      using namespace std;
      //calculate the needed bounding box and dimension to fit all contours
      VolMagick::BoundingBox globalbox(-0.5,-0.5,-0.5,0.5,0.5,0.5);
      VolMagick::Dimension globaldim(4,4,4);
      unsigned int first_slice = 0;
      bool first_slice_initialized = false;
      for(SurfRecon::ContourPtrMap::iterator i = m_Contours[0][0].begin(); //if no vol, inserted at var = 0 time = 0
	  i != m_Contours[0][0].end();
	  i++)
	{
	  if((*i).second == NULL) continue;
	  cerr << "Curve: " << i->first << endl;
	  for(vector<SurfRecon::CurvePtr>::iterator cur = (*i).second->curves().begin();
	      cur != (*i).second->curves().end();
	      cur++)
	    {
	      if(!first_slice_initialized)
		{
		  first_slice = get<0>(**cur);
		  first_slice_initialized = true;
		}
	      
	      first_slice = min(first_slice,get<0>(**cur));
	      
	      //get the max dimension represented in the contour info via slice numbers
	      globaldim.zdim = max(static_cast<unsigned int>(globaldim.zdim),
				   get<0>(**cur));
	      //calc bounding box for contour points
	      if(get<2>(**cur).size() < 3) continue; //don't count non polygonal contours
	      for(SurfRecon::PointPtrList::iterator pcur = get<2>(**cur).begin();
		  pcur != get<2>(**cur).end();
		  pcur++)
		{
		  SurfRecon::Point p = **pcur;

#ifdef CONTOUR_BOUNDING_BOX_FILTER_HEURISTIC
		  const double max_distance = 10000.0; //distance of centerpoint of bounding box to current contour point
		  SurfRecon::Point center_pt(globalbox.minx + (globalbox.maxx-globalbox.minx)/2.0,
					     globalbox.miny + (globalbox.maxy-globalbox.miny)/2.0,
					     globalbox.minz + (globalbox.maxz-globalbox.minz)/2.0);
		  if(CGAL::to_double(CGAL::squared_distance(p,
							    center_pt)) > 
		     max_distance*max_distance)
		    continue;
#endif

		  globalbox.minx = min(globalbox.minx,p.x());
		  globalbox.miny = min(globalbox.miny,p.y());
		  globalbox.minz = min(globalbox.minz,p.z());
		  globalbox.maxx = max(globalbox.maxx,p.x());
		  globalbox.maxy = max(globalbox.maxy,p.y());
		  globalbox.maxz = max(globalbox.maxz,p.z());
		}
	    }
	}
	  
      cerr << str(format("Contours: %1%\nDimension(%2%,%3%,%4%)\n"
			 "BoundingBox(%5%,%6%,%7%,%8%,%9%,%10%)\n")
		  % m_Contours[0][0].size()
		  % globaldim.xdim
		  % globaldim.ydim
		  % globaldim.zdim
		  % globalbox.minx % globalbox.miny % globalbox.minz
		  % globalbox.maxx % globalbox.maxy % globalbox.maxz);

      if(globalbox.isNull())
	{
	  cvcapp.log(5, "Null bounding box for loaded contours..something is wrong!");
	  QMessageBox::critical(this,"Error opening the generated volume",
				"Null bounding box for loaded contours..something is wrong!");
	  m_Contours.clear();
	  _ui->m_Objects->clear();
	  return;
	}

      //now actually generate the volume file
      QString newfilename;
      try
	{
	  QSettings settings(QDir::homePath() + "/.CCV/settings.ini", QSettings::IniFormat);
	  VolMagick::Volume empty_vol;
	  empty_vol.boundingBox(globalbox);
	  empty_vol.dimension(globaldim);
	  QString cacheDirString = settings.value("/Volume Rover/CacheDir", QString()).toString();
	  if(cacheDirString.isEmpty()) cacheDirString = ".";
	  QDir tmpdir(cacheDirString + "/VolumeCache/tmp");
	  if(!tmpdir.exists()) tmpdir.mkpath(tmpdir.absolutePath());
	  QFileInfo tmpfile(tmpdir.absolutePath(),"tmp.rawiv");
	  cvcapp.log(5, boost::str(boost::format("Creating volume %s")%tmpfile.absoluteFilePath().toStdString().c_str()));
	  newfilename = QDir::toNativeSeparators(tmpfile.absoluteFilePath());
	  cvcapp.log(5, boost::str(boost::format("newfilename = %s")%newfilename.toStdString().c_str()));
	  VolMagick::createVolumeFile(cvcapp, newfilename.toStdString().c_str(),
				      empty_vol.boundingBox(),
				      empty_vol.dimension(),
				      vector<VolMagick::VoxelType>(1, empty_vol.voxelType()));
	  VolMagick::writeVolumeFile(cvcapp, empty_vol,newfilename.toStdString().c_str());
	}
      catch(const VolMagick::Exception& e)
	{
	  QMessageBox::critical(this,"Error opening the generated volume",e.what());
	  m_Contours.clear();
	  _ui->m_Objects->clear();
	  return;
	}

      //backup m_Contours because after emitting volumeGenerated it will be obliterated
      SurfRecon::ContourPtrArray contours = m_Contours;

      //do the reload! (hopefully)
      emit volumeGenerated(newfilename);

      m_Contours = contours;

      //rebuild m_Objects based on m_Contours contents
      for(SurfRecon::ContourPtrMap::iterator i = m_Contours[0][0].begin(); //if no vol, inserted at var = 0 time = 0
	  i != m_Contours[0][0].end();
	  i++)
	_ui->m_Objects->insertItem(_ui->m_Objects->count(), new QListWidgetItem(QString(i->first.c_str()), _ui->m_Objects, 0));

      //hop to a slice with contour data so the user isn't confused!
      m_XYSliceCanvas->setDepth(first_slice);
      m_XZSliceCanvas->setDepth(first_slice);
      m_ZYSliceCanvas->setDepth(first_slice);
    }

 end:;

  updateSliceContours();
  showCurrentObject();
  update();
}

void VolumeGridRover::sdfOptionsSlot()
{
/*
  if(!m_VolumeFileInfo.isSet())
    {
      QMessageBox::critical(this,"Error","Load a volume first.",
			    QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
      return;
    }

  sdf_opt sdfd(this,0,true);

  if(sdfd.exec() == QDialog::Accepted)
    {
      std::vector<SDF2D::Polygon_2> polygons;
      SDF2D::Image sdf;
      SDF2D::Dimension dim;
      SDF2D::BoundingBox bbox;

      //get the current shown canvas
      SliceCanvas *cur_canvas = NULL;
      switch(_ui->m_SliceCanvasTab->currentIndex())
	{
	default:
	case 0:
	  {
	    cur_canvas = m_XYSliceCanvas;
	    bbox = SDF2D::BoundingBox(m_VolumeFileInfo.XMin(),m_VolumeFileInfo.YMin(),
				      m_VolumeFileInfo.XMax(),m_VolumeFileInfo.YMax());
	    break;
	  }
	case 1:
	  {
	    cur_canvas = m_XZSliceCanvas;
	    bbox = SDF2D::BoundingBox(m_VolumeFileInfo.XMin(),m_VolumeFileInfo.ZMin(),
				      m_VolumeFileInfo.XMax(),m_VolumeFileInfo.ZMax());
	    break;
	  }
	case 2:
	  {
	    cur_canvas = m_ZYSliceCanvas;
	    bbox = SDF2D::BoundingBox(m_VolumeFileInfo.ZMin(),m_VolumeFileInfo.YMin(),
				      m_VolumeFileInfo.ZMax(),m_VolumeFileInfo.YMax());
	    break;
	  }
	}
            
      if(sdfd.x_sample_res->text().toInt() < 0 ||
	 sdfd.y_sample_res->text().toInt() < 0)
	{
	  QMessageBox::critical(this,"Error","Invalid x/y resolution!",
				QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
	  return;
	}

      if(sdfd.sign_bitmap_only->isChecked())
	{
	  QMessageBox::warning(this,"Warning","Sign bitmap only not yet implemented...",
			       QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
	}

      dim = SDF2D::Dimension(sdfd.x_sample_res->text().toInt(),
			     sdfd.y_sample_res->text().toInt());
      sdf.resize(boost::extents
		 [sdfd.x_sample_res->text().toInt()]
		 [sdfd.y_sample_res->text().toInt()]);

      //finally collect all selected contours and convert them to CGAL polygons, then run SDF2D
      QListViewItemIterator it(m_Objects, QListViewItemIterator::Selected);
      while(it.current()!=NULL)
	{
	  SurfRecon::CurvePtrVector 
	    curves(m_Contours[m_Variable->currentIndex()][m_Timestep->value()][it.current()->text(0).toStdString().c_str()]->curves());
	  
	  const gsl_interp_type *interp;
	  
	  switch(m_Contours
		 [m_Variable->currentIndex()]
		 [m_Timestep->value()]
		 [it.current()->text(0).toStdString().c_str()]->interpolationType())
	    {
	    case 0: interp = gsl_interp_linear; break;
	    case 1: interp = gsl_interp_polynomial; break;
	    case 2: interp = gsl_interp_cspline; break;
	    case 3: interp = gsl_interp_cspline_periodic; break;
	    case 4: interp = gsl_interp_akima; break;
	    case 5: interp = gsl_interp_akima_periodic; break;
	    default: interp = gsl_interp_cspline; break;
	    }

	  for(SurfRecon::CurvePtrVector::iterator j = curves.begin();
	      j != curves.end();
	      j++)
	    {
	      std::vector<double> vec_x, vec_y, vec_z, vec_t; //used to build the input arrays for gsl_spline_init()
	      double t=0.0, interval;
	      gsl_interp_accel *acc_x, *acc_y, *acc_z;
	      gsl_spline *spline_x, *spline_y, *spline_z;

	      //now check if this curve is currently drawn on that canvas
	      if(SurfRecon::getCurveSlice(**j) != cur_canvas->getDepth() ||
		 SurfRecon::getCurveOrientation(**j) != SurfRecon::Orientation(cur_canvas->getSliceAxis()))
		continue;
	      
	      if(SurfRecon::getCurvePoints(**j).size() < 3)
		{
		  cvcapp.log(5, "Warning: skipping curve with less than 3 points...");
		  continue;
		}	    
	      
	      for(SurfRecon::PointPtrList::iterator pcur = SurfRecon::getCurvePoints(**j).begin();
		  pcur != SurfRecon::getCurvePoints(**j).end();
		  pcur++)
		{
		  vec_x.push_back((*pcur)->x());
		  vec_y.push_back((*pcur)->y());
		  vec_z.push_back((*pcur)->z());
		  vec_t.push_back(t); t+=1.0;
		}
	      t-=1.0;

	      SurfRecon::PointPtrVector tmplist;
	      //this is hackish but how else can i find the min size?
	      gsl_spline *test_spline = gsl_spline_alloc(interp, 1000);
	      if(vec_t.size() >= gsl_spline_min_size(test_spline))
		{
		  double time_i;
		  //int count;
		  
		  acc_x = gsl_interp_accel_alloc();
		  acc_y = gsl_interp_accel_alloc();
		  acc_z = gsl_interp_accel_alloc();
		  spline_x = gsl_spline_alloc(interp, vec_t.size());
		  spline_y = gsl_spline_alloc(interp, vec_t.size());
		  spline_z = gsl_spline_alloc(interp, vec_t.size());
		  gsl_spline_init(spline_x, &(vec_t[0]), &(vec_x[0]), vec_t.size());
		  gsl_spline_init(spline_y, &(vec_t[0]), &(vec_y[0]), vec_t.size());
		  gsl_spline_init(spline_z, &(vec_t[0]), &(vec_z[0]), vec_t.size());
		  
		  //Just sample some number of points for non-linear interpolation for now.
		  interval = 1.0/(1+m_Contours
				  [m_Variable->currentIndex()]
				  [m_Timestep->value()]
				  [it.current()->text(0).toStdString().c_str()]->numberOfSamples());
		  //In the future, sample according to the spline's curvature between 2
		  // control points...
		  
		  //this bothers me... please just calculate the number of points in the future!
		  //for(time_i = 0.0, count=0; time_i<=t; time_i+=interval, count++);
		  //outfile << count << std::endl;
		  
		  for(time_i = 0.0; time_i<=t; time_i+=interval)
		    {
		      //outfile << gsl_spline_eval(spline_x, time_i, acc_x) << " " 
			//<< gsl_spline_eval(spline_y, time_i, acc_y) << " " 
			//<< gsl_spline_eval(spline_z, time_i, acc_z) << std::endl;
		      tmplist.push_back(SurfRecon::PointPtr(new SurfRecon::Point(gsl_spline_eval(spline_x, time_i, acc_x),
										 gsl_spline_eval(spline_y, time_i, acc_y),
										 gsl_spline_eval(spline_z, time_i, acc_z))));
		    }
		  
		  // clean up our mess
		  gsl_spline_free(spline_x);
		  gsl_spline_free(spline_y);
		  gsl_spline_free(spline_z);
		  gsl_interp_accel_free(acc_x);
		  gsl_interp_accel_free(acc_y);
		  gsl_interp_accel_free(acc_z);
		  
		  SurfRecon::PointPtr pfront = tmplist.front();
		  SurfRecon::PointPtr pback = tmplist.back();
		  
		  //remove the last point if the same as first
		  if((fabs(pfront->x() - pback->x()) <= 0.00001) &&
		     (fabs(pfront->y() - pback->y()) <= 0.00001) &&
		     (fabs(pfront->z() - pback->z()) <= 0.00001))
		    tmplist.pop_back();
		}
	      gsl_spline_free(test_spline);

	      std::list<SDF2D::Point_2> points;
	      switch(SurfRecon::getCurveOrientation(**j))
		{
		case SurfRecon::XY:
		  for(SurfRecon::PointPtrVector::iterator k = tmplist.begin();
		      k != tmplist.end();
		      k++)
		    points.push_back(SDF2D::Point_2((*k)->x(),(*k)->y()));
		  break;
		case SurfRecon::XZ:
		  for(SurfRecon::PointPtrVector::iterator k = tmplist.begin();
		      k != tmplist.end();
		      k++)
		    points.push_back(SDF2D::Point_2((*k)->x(),(*k)->z()));
		  break;
		case SurfRecon::ZY:
		  for(SurfRecon::PointPtrVector::iterator k = tmplist.begin();
		      k != tmplist.end();
		      k++)
		    points.push_back(SDF2D::Point_2((*k)->z(),(*k)->y()));
		  break;
		}
	      
	      //points.pop_back(); //the last point is a duplicate of the first
	      
	      SDF2D::Polygon_2 pgn(points.begin(),points.end());
	      if(!pgn.is_simple())
		{
		  cvcapp.log(5, "Warning: self intersecting polygon detected! skipping...");
		  ++it;
		  continue;
		}
	      polygons.push_back(pgn);
	    }
	  
	  ++it;
	}
      
      try
	{
	  sdf = SDF2D::signedDistanceFunction(polygons,dim,bbox,
					      SDF2D::SignMethod(sdfd.sign_method->currentIndex()),
					      SDF2D::DistanceMethod(sdfd.dist_method->currentIndex()));
	
	  //output a B&W image showing inside or out
	  for(SDF2D::ImageIndex i = 0; i < sdf.shape()[0]; i++)
	    for(SDF2D::ImageIndex j = 0; j < sdf.shape()[1]; j++)
	      {
		if(sdf[i][j] > 0) sdf[i][j] = 1.0;
		else sdf[i][j] = 0.0;
	      }
	  
	  // build image for writing ...
	  QImage distanceImage(sdf.shape()[0], sdf.shape()[1], 32);
	  cvcapp.log(5, "num colors: %d",distanceImage.numColors());
	  cvcapp.log(5, "w,h: %d,%d",int(sdf.shape()[0]),int(sdf.shape()[1]));
	  for(SDF2D::ImageIndex i = 0; i < sdf.shape()[0]; i++)
	    for(SDF2D::ImageIndex j = 0; j < sdf.shape()[1]; j++)
	      {
		QRgb color = QColor(sdf[i][j]*255.0,sdf[i][j]*255.0,sdf[i][j]*255.0).rgb();		
		int fixedY = sdf.shape()[1]-j-1;
		distanceImage.setPixel(i,fixedY,color);
	      }
	  
	  // save PPM image to filename TEMP_PPM_PATH
	  //distanceImage.save("./tmp2.png", "PNG");
	  emit showImage(distanceImage);
	}
      catch(std::exception& e)
	{
	  std::cerr << e.what() << std::endl;
	}
    }
*/
}

//calculate medial axis for selected closed contours, then output open contours for medial axis
//so they can be visualized
void VolumeGridRover::medialAxisSlot()
{
/*
#ifdef USING_VOLUMEGRIDROVER_MEDAX
  if(!m_VolumeFileInfo.isSet())
    {
      QMessageBox::critical(this,"Error","Load a volume first.",
			    QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
      return;
    }

  //get the current shown canvas
  SliceCanvas *cur_canvas = NULL;
  switch(_ui->m_SliceCanvasTab->currentIndex())
    {
    default:
    case 0:
      {
	cur_canvas = m_XYSliceCanvas;
	break;
      }
    case 1:
      {
	cur_canvas = m_XZSliceCanvas;
	break;
      }
    case 2:
      {
	cur_canvas = m_ZYSliceCanvas;
	break;
      }
    }

  //finally collect all selected contours and convert them to CGAL polygons, then run medax
  std::vector<medax::Polygon_2> polygons;
  QListViewItemIterator it(m_Objects, QListViewItemIterator::Selected);
  while(it.current()!=NULL)
    {
      SurfRecon::CurvePtrVector 
	curves(m_Contours[m_Variable->currentIndex()][m_Timestep->value()][it.current()->text(0).toStdString().c_str()]->curves());
	  
      const gsl_interp_type *interp;
	  
      switch(m_Contours
	     [m_Variable->currentIndex()]
	     [m_Timestep->value()]
	     [it.current()->text(0).toStdString().c_str()]->interpolationType())
	{
	case 0: interp = gsl_interp_linear; break;
	case 1: interp = gsl_interp_polynomial; break;
	case 2: interp = gsl_interp_cspline; break;
	case 3: interp = gsl_interp_cspline_periodic; break;
	case 4: interp = gsl_interp_akima; break;
	case 5: interp = gsl_interp_akima_periodic; break;
	default: interp = gsl_interp_cspline; break;
	}

      for(SurfRecon::CurvePtrVector::iterator j = curves.begin();
	  j != curves.end();
	  j++)
	{
	  std::vector<double> vec_x, vec_y, vec_z, vec_t; //used to build the input arrays for gsl_spline_init()
	  double t=0.0, interval;
	  gsl_interp_accel *acc_x, *acc_y, *acc_z;
	  gsl_spline *spline_x, *spline_y, *spline_z;

	  //now check if this curve is currently drawn on that canvas
	  if(SurfRecon::getCurveSlice(**j) != cur_canvas->getDepth() ||
	     SurfRecon::getCurveOrientation(**j) != SurfRecon::Orientation(cur_canvas->getSliceAxis()))
	    continue;
	      
	  if(SurfRecon::getCurvePoints(**j).size() < 3)
	    {
	      cvcapp.log(5, "Warning: skipping curve with less than 3 points...");
	      continue;
	    }	    
	      
	  for(SurfRecon::PointPtrList::iterator pcur = SurfRecon::getCurvePoints(**j).begin();
	      pcur != SurfRecon::getCurvePoints(**j).end();
	      pcur++)
	    {
	      vec_x.push_back((*pcur)->x());
	      vec_y.push_back((*pcur)->y());
	      vec_z.push_back((*pcur)->z());
	      vec_t.push_back(t); t+=1.0;
	    }
	  t-=1.0;

	  SurfRecon::PointPtrVector tmplist;

	  if(interp == gsl_interp_linear) //if linear, no need to actually interpolate, just fill tmplist with the actual points
	    {
	      for(SurfRecon::PointPtrList::iterator pcur = SurfRecon::getCurvePoints(**j).begin();
		  pcur != SurfRecon::getCurvePoints(**j).end();
		  pcur++)
		tmplist.push_back(*pcur);
	    }
	  else
	    {
	      //this is hackish but how else can i find the min size?
	      gsl_spline *test_spline = gsl_spline_alloc(interp, 1000);
	      if(vec_t.size() >= gsl_spline_min_size(test_spline))
		{
		  double time_i;
		  //int count;
		  
		  acc_x = gsl_interp_accel_alloc();
		  acc_y = gsl_interp_accel_alloc();
		  acc_z = gsl_interp_accel_alloc();
		  spline_x = gsl_spline_alloc(interp, vec_t.size());
		  spline_y = gsl_spline_alloc(interp, vec_t.size());
		  spline_z = gsl_spline_alloc(interp, vec_t.size());
		  gsl_spline_init(spline_x, &(vec_t[0]), &(vec_x[0]), vec_t.size());
		  gsl_spline_init(spline_y, &(vec_t[0]), &(vec_y[0]), vec_t.size());
		  gsl_spline_init(spline_z, &(vec_t[0]), &(vec_z[0]), vec_t.size());
		  
		  //Just sample some number of points for non-linear interpolation for now.
		  interval = 1.0/(1+m_Contours
				  [m_Variable->currentIndex()]
				  [m_Timestep->value()]
				  [it.current()->text(0).toStdString().c_str()]->numberOfSamples());
		  //In the future, sample according to the spline's curvature between 2
		  // control points...
		  
		  //this bothers me... please just calculate the number of points in the future!
		  //for(time_i = 0.0, count=0; time_i<=t; time_i+=interval, count++);
		  //outfile << count << std::endl;
		  
		  for(time_i = 0.0; time_i<=t; time_i+=interval)
		    {
		      //outfile << gsl_spline_eval(spline_x, time_i, acc_x) << " " 
			//<< gsl_spline_eval(spline_y, time_i, acc_y) << " " 
			//<< gsl_spline_eval(spline_z, time_i, acc_z) << std::endl;
		      tmplist.push_back(SurfRecon::PointPtr(new SurfRecon::Point(gsl_spline_eval(spline_x, time_i, acc_x),
										 gsl_spline_eval(spline_y, time_i, acc_y),
										 gsl_spline_eval(spline_z, time_i, acc_z))));
		    }
		  
		  // clean up our mess
		  gsl_spline_free(spline_x);
		  gsl_spline_free(spline_y);
		  gsl_spline_free(spline_z);
		  gsl_interp_accel_free(acc_x);
		  gsl_interp_accel_free(acc_y);
		  gsl_interp_accel_free(acc_z);
		  
		  SurfRecon::PointPtr pfront = tmplist.front();
		  SurfRecon::PointPtr pback = tmplist.back();
		  
		  //remove the last point if the same as first
		  if((fabs(pfront->x() - pback->x()) <= 0.00001) &&
		     (fabs(pfront->y() - pback->y()) <= 0.00001) &&
		     (fabs(pfront->z() - pback->z()) <= 0.00001))
		    tmplist.pop_back();
		}
	      gsl_spline_free(test_spline);
	    }

	  std::list<medax::Poly_point_2> points;
	  switch(SurfRecon::getCurveOrientation(**j))
	    {
	    case SurfRecon::XY:
	      for(SurfRecon::PointPtrVector::iterator k = tmplist.begin();
		  k != tmplist.end();
		  k++)
		points.push_back(medax::Poly_point_2((*k)->x(),(*k)->y()));
	      break;
	    case SurfRecon::XZ:
	      for(SurfRecon::PointPtrVector::iterator k = tmplist.begin();
		  k != tmplist.end();
		  k++)
		points.push_back(medax::Poly_point_2((*k)->x(),(*k)->z()));
	      break;
	    case SurfRecon::ZY:
	      for(SurfRecon::PointPtrVector::iterator k = tmplist.begin();
		  k != tmplist.end();
		  k++)
		points.push_back(medax::Poly_point_2((*k)->z(),(*k)->y()));
	      break;
	    }
	      
	  //points.pop_back(); //the last point is a duplicate of the first
	      
	  medax::Polygon_2 pgn(points.begin(),points.end());
#if 0
	  if(!pgn.is_simple())
	    {
	      cvcapp.log(5, "Warning: self intersecting polygon detected! skipping...");
	      ++it;
	      continue;
	    }
#endif
	  polygons.push_back(pgn);
	}
	  
      ++it;
    }

  //compute the medial axes
  medax::Edges result = medax::compute(polygons);
  //medax::Edges result = medax::computeVoronoi(polygons);
  
  SurfRecon::Contour medax_contour(SurfRecon::Color(m_ContourColor->paletteBackgroundColor().red()/255.0,
						    m_ContourColor->paletteBackgroundColor().green()/255.0,
						    m_ContourColor->paletteBackgroundColor().blue()/255.0),
				   "medax",
				   SurfRecon::InterpolationType(m_InterpolationType->currentIndex()),
				   m_InterpolationSampling->value());
  
  //1 curve per halfedge for now... TODO: combine half edges into curves to reduce overhead
  for(medax::Edges::iterator i = result.begin();
      i != result.end();
      i++)
    {
      SurfRecon::Curve c(cur_canvas->getDepth(),
			 SurfRecon::Orientation(cur_canvas->getSliceAxis()),
			 SurfRecon::PointPtrList(),
			 "medax");
      
      //if(result.get<1>().find(*i) == result.get<1>().end()) continue; //this edge not found in set
      //if(!i->source->keep()) continue;

      //we only want bounded edges
      //if(!i->is_segment()) continue;

      switch(cur_canvas->getSliceAxis())
	{
	case SliceCanvas::XY:
	  {
	    double depth_coord = cur_canvas->getDepth()*m_VolumeFileInfo.ZSpan()+m_VolumeFileInfo.ZMin();
	    getCurvePoints(c).push_back(SurfRecon::PointPtr(new SurfRecon::Point(i->get<0>().x(),
										 i->get<0>().y(),
										 depth_coord)));
	    getCurvePoints(c).push_back(SurfRecon::PointPtr(new SurfRecon::Point(i->get<1>().x(),
										 i->get<1>().y(),
										 depth_coord)));
	  }
	  break;
	case SliceCanvas::XZ:
	  {
	    double depth_coord = cur_canvas->getDepth()*m_VolumeFileInfo.YSpan()+m_VolumeFileInfo.YMin();
	    getCurvePoints(c).push_back(SurfRecon::PointPtr(new SurfRecon::Point(i->get<0>().x(),
										 depth_coord,
										 i->get<0>().y())));
	    getCurvePoints(c).push_back(SurfRecon::PointPtr(new SurfRecon::Point(i->get<1>().x(),
										 depth_coord,
										 i->get<1>().y())));
	  }
	  break;
	case SliceCanvas::ZY:
	  {
	    double depth_coord = cur_canvas->getDepth()*m_VolumeFileInfo.XSpan()+m_VolumeFileInfo.XMin();
	    getCurvePoints(c).push_back(SurfRecon::PointPtr(new SurfRecon::Point(depth_coord,
										 i->get<0>().y(),
										 i->get<0>().x())));
	    getCurvePoints(c).push_back(SurfRecon::PointPtr(new SurfRecon::Point(depth_coord,
										 i->get<1>().y(),
										 i->get<1>().x())));
	  }
	  break;
	}

      medax_contour.add(SurfRecon::CurvePtr(new SurfRecon::Curve(c)));
    }

  addContour(medax_contour);
  
  //convert axes into VolumeGridRover contours so they can be visualized
  //cvcapp.log(5, "num edges in medax: %d",result.get<1>().size());
#else
  QMessageBox::critical(this,"Error","Medial axis calculation disabled.",
			QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
#endif
*/
}

void VolumeGridRover::curateContoursSlot()
{
  cvcapp.log(5, "VolumeGridRover::curateContoursSlot()");
}

inline const QList<PointClass*>* VolumeGridRover::getPointClassList() const { return m_PointClassList[_ui->m_Variable->currentIndex()][_ui->m_Timestep->value()]; }

void VolumeGridRover::addContourSlot()
{
//  cvcapp.log(5, "VolumeGridRover::addContoursSlot()");

  static int classNum = 0;

  if(m_Contours.empty()) {
    fprintf( stderr, "Contour list is empty\n");
    return;
  }

  QString name;
  name = QString("Contour %1").arg(classNum++);

  //cvcapp.log(5, "VolumeGridRover::addContourSlot()");
  QColor color = _ui->m_ContourColor->palette().color( QPalette::Button );

  m_Contours[_ui->m_Variable->currentIndex()][_ui->m_Timestep->value()][name.toStdString().c_str()] =
    SurfRecon::ContourPtr(new SurfRecon::Contour(SurfRecon::Color( color.red()/255.0,
								   color.green()/255.0,
								   color.blue()/255.0),
						 std::string(name.toStdString().c_str()),
						 SurfRecon::InterpolationType(_ui->m_InterpolationType->currentIndex()),
						 _ui->m_InterpolationSampling->value()));
  //m_Contour->insertItem(name);
  _ui->m_Objects->insertItem( _ui->m_Objects->count(), new QListWidgetItem(name, _ui->m_Objects, 0));

  m_XYSliceCanvas->setContours(m_Contours);
  m_XZSliceCanvas->setContours(m_Contours);
  m_ZYSliceCanvas->setContours(m_Contours);

  //m_Contour->setCurrentItem(m_Contour->count()-1);
  _ui->m_Objects->setCurrentItem( _ui->m_Objects->findItems(name,Qt::MatchExactly).first() );

  m_XYSliceCanvas->setCurrentContour(name.toStdString().c_str());
  m_XZSliceCanvas->setCurrentContour(name.toStdString().c_str());
  m_ZYSliceCanvas->setCurrentContour(name.toStdString().c_str());
  showContourColor(name.toStdString().c_str());
  showContourInterpolationType(name.toStdString().c_str());
  showContourInterpolationSampling(name.toStdString().c_str());

  m_XYSliceCanvas->update();
  m_XZSliceCanvas->update();
  m_ZYSliceCanvas->update();

}

void VolumeGridRover::deleteContourSlot()
{
//  cvcapp.log(5, "VolumeGridRover::deleteContourSlot()");

  int size = _ui->m_Objects->count();
  if( size < 1 ) return;

  for( int i=0; i<size; i++ )
  {
      QListWidgetItem *item = _ui->m_Objects->item( i );
      m_Contours[_ui->m_Variable->currentIndex()][_ui->m_Timestep->value()].erase( item->text().toStdString().c_str());
      delete item;
  }
    
  m_XYSliceCanvas->setContours(m_Contours);
  m_XZSliceCanvas->setContours(m_Contours);
  m_ZYSliceCanvas->setContours(m_Contours);

  m_XYSliceCanvas->setCurrentContour(_ui->m_Objects->currentItem() ? _ui->m_Objects->currentItem()->text().toStdString().c_str() : "");
  m_XZSliceCanvas->setCurrentContour(_ui->m_Objects->currentItem() ? _ui->m_Objects->currentItem()->text().toStdString().c_str() : "");
  m_ZYSliceCanvas->setCurrentContour(_ui->m_Objects->currentItem() ? _ui->m_Objects->currentItem()->text().toStdString().c_str() : "");
  showContourColor(_ui->m_Objects->currentItem() ? _ui->m_Objects->currentItem()->text().toStdString().c_str() : "");
  showContourInterpolationType(_ui->m_Objects->currentItem() ? _ui->m_Objects->currentItem()->text().toStdString().c_str() : "");
  showContourInterpolationSampling(_ui->m_Objects->currentItem() ? _ui->m_Objects->currentItem()->text().toStdString().c_str() : "");

  m_XYSliceCanvas->update();
  m_XZSliceCanvas->update();
  m_ZYSliceCanvas->update();

}

void VolumeGridRover::contourColorSlot()
{
//  cvcapp.log(5, "VolumeGridRover::contourColorSlot()");

  QPalette palette = _ui->m_ContourColor->palette();
  QColor color = QColorDialog::getColor( palette.color(QPalette::Button) );
  if(color.isValid())
    {
      palette.setColor( QPalette::Normal, QPalette::Button, color );
      _ui->m_ContourColor->setPalette( palette );
     
      int size = _ui->m_Objects->count();
      for( int i=0; i<size; i++ )
      {
          QListWidgetItem *item = _ui->m_Objects->item( i );
 
	  if(m_Contours[_ui->m_Variable->currentIndex()][_ui->m_Timestep->value()][item->text().toStdString().c_str()] != NULL)
	    m_Contours[_ui->m_Variable->currentIndex()][_ui->m_Timestep->value()][item->text().toStdString().c_str()]->
	      color(SurfRecon::Color(color.red()/255.0,
				     color.green()/255.0,
				     color.blue()/255.0));
      }
    }

  m_XYSliceCanvas->update();
  m_XZSliceCanvas->update();
  m_ZYSliceCanvas->update();

}

void VolumeGridRover::setInterpolationTypeSlot(int t)
{
//  cvcapp.log(5, "VolumeGridRover::setInterpolationTypeSlot()");

  int size = _ui->m_Objects->count();
  for( int i=0; i < size; i++ )
  {
      QListWidgetItem *item = _ui->m_Objects->item( i );
      if(m_Contours[_ui->m_Variable->currentIndex()][_ui->m_Timestep->value()][item->text().toStdString().c_str()] != NULL)
	m_Contours[_ui->m_Variable->currentIndex()][_ui->m_Timestep->value()][item->text().toStdString().c_str()]->
	  interpolationType(SurfRecon::InterpolationType(t));
  }
   
  m_XYSliceCanvas->update();
  m_XZSliceCanvas->update();
  m_ZYSliceCanvas->update();
}

void VolumeGridRover::setInterpolationSamplingSlot(int s)
{
//  cvcapp.log(5, "VolumeGridRover::setInterpolationSamplingSlot()");

  int size = _ui->m_Objects->count();
  for( int i=0; i<size; i++ )
  {
    QListWidgetItem *item = _ui->m_Objects->item( i );
    if(m_Contours[_ui->m_Variable->currentIndex()][_ui->m_Timestep->value()][item->text().toStdString().c_str()] != NULL)
	m_Contours[_ui->m_Variable->currentIndex()][_ui->m_Timestep->value()][item->text().toStdString().c_str()]->
	  numberOfSamples(s);
  }

  m_XYSliceCanvas->update();
  m_XZSliceCanvas->update();
  m_ZYSliceCanvas->update();

}

void VolumeGridRover::tilingRunSlot()
{
  cvcapp.log(5, "tilingRunSlot(): start");
/*
  if(!m_VolumeFileInfo.isSet())
    {
      QMessageBox::critical(this,"Error","Load a volume first.",
			    QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
      return;
    }
  
  if(m_TilingThread.running())
    QMessageBox::critical(this,"Error","Tiling thread already running!",
			  QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
  else
    {
      std::list<std::string> names;

      //tile selected objects...
      QListViewItemIterator it(m_Objects, QListViewItemIterator::Selected);
      while(it.current()!=NULL)
	{
	  names.push_back(it.current()->text(0).toStdString().c_str());
	  ++it;
	}

      m_TilingThread.setSelectedNames(names);
      m_TilingThread.start();
    }

  cvcapp.log(5, "tilingRunSlot(): finish");
*/
}

void VolumeGridRover::getTilingOutputDirectorySlot()
{
  cvcapp.log(5, "VolumeGridRover::getTilingOutputDirectorySlot()");
/*
  m_TilingOutputDirectory->setText(QFileDialog::getExistingDirectory(QString(),this,"get existing directory","Choose a directory",TRUE));
*/
}

void VolumeGridRover::handleTilingOutputDestinationSelectionSlot(int selection)
{
  cvcapp.log(5, "VolumeGridRover::handleTilingOutputDestinationSelectionSlot()");
/*
  switch(selection)
    {
    case 0:
      m_TilingOutputFilenameGroup->setEnabled(false);
      break;
    case 1:
      m_TilingOutputFilenameGroup->setEnabled(true);
      break;
    default: break;
    }
*/
}

//TODO: finish this!!  It's supposed to calc SDF then magically curate using SDF info
//Right now it's got a bunch of garbage that I've now moved to the new sdfOptionsSlot()
void VolumeGridRover::sdfCurationSlot()
{
  cvcapp.log(5, "VolumeGridRover::sdfCurationSlot()");
/*
  std::vector<SDF2D::Polygon_2> polygons;
  SDF2D::Image sdf;
  SDF2D::Dimension dim;
  SDF2D::BoundingBox bbox;

  cvcapp.log(5, "VolumeGridRover::sdfCurationSlot(): starting SDF curation of selected 2D contours");

  //get the current shown canvas
  SliceCanvas *cur_canvas = NULL;
  switch(m_SliceCanvasStack->id(m_SliceCanvasStack->visibleWidget()))
    {
    default:
    case 0:
      {
	cur_canvas = m_XYSliceCanvas;
	dim = SDF2D::Dimension(m_VolumeFileInfo.XDim(),m_VolumeFileInfo.YDim());
	bbox = SDF2D::BoundingBox(m_VolumeFileInfo.XMin(),m_VolumeFileInfo.YMin(),
				  m_VolumeFileInfo.XMax(),m_VolumeFileInfo.YMax());
	sdf.resize(boost::extents[m_VolumeFileInfo.XDim()][m_VolumeFileInfo.YDim()]);
	break;
      }
    case 1:
      {
	cur_canvas = m_XZSliceCanvas;
	dim = SDF2D::Dimension(m_VolumeFileInfo.XDim(),m_VolumeFileInfo.ZDim());
	bbox = SDF2D::BoundingBox(m_VolumeFileInfo.XMin(),m_VolumeFileInfo.ZMin(),
				  m_VolumeFileInfo.XMax(),m_VolumeFileInfo.ZMax());
	sdf.resize(boost::extents[m_VolumeFileInfo.XDim()][m_VolumeFileInfo.ZDim()]);
	break;
      }
    case 2:
      {
	cur_canvas = m_ZYSliceCanvas;
	dim = SDF2D::Dimension(m_VolumeFileInfo.ZDim(),m_VolumeFileInfo.YDim());
	bbox = SDF2D::BoundingBox(m_VolumeFileInfo.ZMin(),m_VolumeFileInfo.YMin(),
				  m_VolumeFileInfo.ZMax(),m_VolumeFileInfo.YMax());
	sdf.resize(boost::extents[m_VolumeFileInfo.ZDim()][m_VolumeFileInfo.YDim()]);
	break;
      }
    }

  QListViewItemIterator it(m_Objects, QListViewItemIterator::Selected);
  while(it.current()!=NULL)
    {
#if 0
      //for every curve belonging to this contour, convert that curve into a CGAL polygon and run 2D SDF on that polygon
      SurfRecon::CurvePtrVector 
	curves(m_Contours[m_Variable->currentIndex()][m_Timestep->value()][it.current()->text(0).toStdString().c_str()]->curves());
      for(SurfRecon::CurvePtrVector::iterator j = curves.begin();
	  j != curves.end();
	  j++)
	{
	  std::list<SDF2D::Point_2> points;
	  SDF2D::Image sdf;
	  SDF2D::Dimension dim;
	  SDF2D::BoundingBox bbox;

	  switch(SurfRecon::getCurveOrientation(**j))
	    {
	    case SurfRecon::XY:
	      for(SurfRecon::PointPtrList::iterator k = SurfRecon::getCurvePoints(**j).begin();
		  k != SurfRecon::getCurvePoints(**j).end();
		  k++)
		points.push_back(SDF2D::Point_2((*k)->x(),(*k)->y()));
	      dim = SDF2D::Dimension(m_VolumeFileInfo.XDim(),m_VolumeFileInfo.YDim());
	      bbox = SDF2D::BoundingBox(m_VolumeFileInfo.XMin(),m_VolumeFileInfo.YMin(),
					m_VolumeFileInfo.XMax(),m_VolumeFileInfo.YMax());
	      sdf.resize(boost::extents[m_VolumeFileInfo.XDim()][m_VolumeFileInfo.YDim()]);
	      break;
	    case SurfRecon::XZ:
	      for(SurfRecon::PointPtrList::iterator k = SurfRecon::getCurvePoints(**j).begin();
		  k != SurfRecon::getCurvePoints(**j).end();
		  k++)
		points.push_back(SDF2D::Point_2((*k)->x(),(*k)->z()));
	      dim = SDF2D::Dimension(m_VolumeFileInfo.XDim(),m_VolumeFileInfo.ZDim());
	      bbox = SDF2D::BoundingBox(m_VolumeFileInfo.XMin(),m_VolumeFileInfo.ZMin(),
					m_VolumeFileInfo.XMax(),m_VolumeFileInfo.ZMax());
	      sdf.resize(boost::extents[m_VolumeFileInfo.XDim()][m_VolumeFileInfo.ZDim()]);
	      break;
	    case SurfRecon::ZY:
	      for(SurfRecon::PointPtrList::iterator k = SurfRecon::getCurvePoints(**j).begin();
		  k != SurfRecon::getCurvePoints(**j).end();
		  k++)
		points.push_back(SDF2D::Point_2((*k)->z(),(*k)->y()));
	      dim = SDF2D::Dimension(m_VolumeFileInfo.ZDim(),m_VolumeFileInfo.YDim());
	      bbox = SDF2D::BoundingBox(m_VolumeFileInfo.ZMin(),m_VolumeFileInfo.YMin(),
					m_VolumeFileInfo.ZMax(),m_VolumeFileInfo.YMax());
	      sdf.resize(boost::extents[m_VolumeFileInfo.ZDim()][m_VolumeFileInfo.YDim()]);
	      break;
	    }

	  SDF2D::Polygon_2 pgn(points.begin(),points.end());

	  try
	    {
	      sdf = SDF2D::signedDistanceFunction(pgn,dim,bbox);
	    }
	  catch(std::exception &e)
	    {
	      std::cerr << e.what() << std::endl;
	    }
	}
#endif

      //for every curve belonging to this contour that's on the current drawn slice, 
      //convert that curve into a CGAL polygon and run 2D SDF on all polygons
      SurfRecon::CurvePtrVector 
	curves(m_Contours[m_Variable->currentIndex()][m_Timestep->value()][it.current()->text(0).toStdString().c_str()]->curves());

      for(SurfRecon::CurvePtrVector::iterator j = curves.begin();
	  j != curves.end();
	  j++)
	{
	  std::list<SDF2D::Point_2> points;
	  
	  //now check if this curve is currently drawn on that canvas
	  if(SurfRecon::getCurveSlice(**j) != cur_canvas->getDepth() ||
	     SurfRecon::getCurveOrientation(**j) != SurfRecon::Orientation(cur_canvas->getSliceAxis()))
	    continue;

	  if(SurfRecon::getCurvePoints(**j).size() < 3)
	    {
	      cvcapp.log(5, "Warning: skipping curve with less than 3 points...");
	      continue;
	    }	    

	  switch(SurfRecon::getCurveOrientation(**j))
	    {
	    case SurfRecon::XY:
	      for(SurfRecon::PointPtrList::iterator k = SurfRecon::getCurvePoints(**j).begin();
		  k != SurfRecon::getCurvePoints(**j).end();
		  k++)
		points.push_back(SDF2D::Point_2((*k)->x(),(*k)->y()));
	      break;
	    case SurfRecon::XZ:
	      for(SurfRecon::PointPtrList::iterator k = SurfRecon::getCurvePoints(**j).begin();
		  k != SurfRecon::getCurvePoints(**j).end();
		  k++)
		points.push_back(SDF2D::Point_2((*k)->x(),(*k)->z()));
	      break;
	    case SurfRecon::ZY:
	      for(SurfRecon::PointPtrList::iterator k = SurfRecon::getCurvePoints(**j).begin();
		  k != SurfRecon::getCurvePoints(**j).end();
		  k++)
		points.push_back(SDF2D::Point_2((*k)->z(),(*k)->y()));
	      break;
	    }

	  points.pop_back(); //the last point is a duplicate of the first

	  SDF2D::Polygon_2 pgn(points.begin(),points.end());
	  if(!pgn.is_simple())
	    {
	      cvcapp.log(5, "Warning: self intersecting polygon detected! skipping...");
	      ++it;
	      continue;
	    }
	  polygons.push_back(pgn);
	}

      ++it;
    }

  //cvcapp.log(5, "polygons size: %d",polygons.size());

  try
    {
      sdf = SDF2D::signedDistanceFunction(polygons,dim,bbox,SDF2D::ANGLE_SUM);

      //output a B&W image showing inside or out
      for(SDF2D::ImageIndex i = 0; i < sdf.shape()[0]; i++)
	for(SDF2D::ImageIndex j = 0; j < sdf.shape()[1]; j++)
	  {
	    if(sdf[i][j] > 0) sdf[i][j] = 1.0;
	    else sdf[i][j] = 0.0;
	  }

      // build image for writing ...
      QImage distanceImage(sdf.shape()[0], sdf.shape()[1], 32);
      cvcapp.log(5, "num colors: %d",distanceImage.numColors());
      cvcapp.log(5, "w,h: %d,%d",int(sdf.shape()[0]),int(sdf.shape()[1]));
      for(SDF2D::ImageIndex i = 0; i < sdf.shape()[0]; i++)
	for(SDF2D::ImageIndex j = 0; j < sdf.shape()[1]; j++)
	  {
	    QRgb color = QColor(sdf[i][j]*255.0,sdf[i][j]*255.0,sdf[i][j]*255.0).rgb();

	    int fixedY = sdf.shape()[1]-j-1;
	    distanceImage.setPixel(i,fixedY,color);
	  }
      
      // save PPM image to filename TEMP_PPM_PATH
      distanceImage.save("./tmp2.png", "PNG");
    }
  catch(std::exception &e)
    {
      std::cerr << e.what() << std::endl;
    }
  
  cvcapp.log(5, "VolumeGridRover::sdfCurationSlot(): finished SDF curation");
*/
}

bool VolumeGridRover::setVolume(const VolMagick::VolumeFileInfo& vfi)
{
  unsigned int i,j;
  //unsigned int k,v,t;
  unsigned int xdepth, ydepth, zdepth; /* used to keep the current depth value */

  zdepth = _ui->m_XYDepthSlide->value();
  ydepth = _ui->m_XZDepthSlide->value();
  xdepth = _ui->m_ZYDepthSlide->value();

  unsetVolume();

  if(!vfi.isSet()) return false;
  m_VolumeFileInfo = vfi;


  if( m_PointClassList == NULL )
  {
    /* loading successful, now lets initialize the widget with new info */
    m_PointClassList = new QList<PointClass*>**[m_VolumeFileInfo.numVariables()];
    for(i=0; i<m_VolumeFileInfo.numVariables(); i++)
      {
        m_PointClassList[i] = new QList<PointClass*>*[m_VolumeFileInfo.numTimesteps()];
        for(j=0; j<m_VolumeFileInfo.numTimesteps(); j++)
          {
	    m_PointClassList[i][j] = new QList<PointClass*>();
//Q3err	    m_PointClassList[i][j]->setAutoDelete(true);
	  }
      }
  }
//  m_Contours.resize(m_VolumeFileInfo.numVariables());
//  for(i=0; i<m_VolumeFileInfo.numVariables(); i++)
//    m_Contours[i].resize(m_VolumeFileInfo.numTimesteps());

  _ui->m_Variable->clear();
  for(i=0; i<m_VolumeFileInfo.numVariables(); i++)
     _ui->m_Variable->insertItem(i, QString(m_VolumeFileInfo.name(i).data()));

  _ui->m_Variable->setCurrentIndex(0);
  _ui->m_Timestep->setMinimum(0);
  _ui->m_Timestep->setMaximum(m_VolumeFileInfo.numTimesteps()-1);
  _ui->m_Timestep->setValue(0);


  setMinValue(m_VolumeFileInfo.min(_ui->m_Variable->currentIndex(),_ui->m_Timestep->value()));
  setMaxValue(m_VolumeFileInfo.max(_ui->m_Variable->currentIndex(),_ui->m_Timestep->value()));

  _ui->m_XYDepthSlide->setMaximum(m_VolumeFileInfo.ZDim()-1);
  m_XYSliceCanvas->setVolume(m_VolumeFileInfo);
  m_XYSliceCanvas->setPointClassList(m_PointClassList);
//  m_XYSliceCanvas->setContours(m_Contours);

  _ui->m_XZDepthSlide->setMaximum(m_VolumeFileInfo.YDim()-1);
  m_XZSliceCanvas->setVolume(m_VolumeFileInfo);
  m_XZSliceCanvas->setPointClassList(m_PointClassList);
//  m_XZSliceCanvas->setContours(m_Contours);
 
  _ui->m_ZYDepthSlide->setMaximum(m_VolumeFileInfo.XDim()-1);
  m_ZYSliceCanvas->setVolume(m_VolumeFileInfo);
  m_ZYSliceCanvas->setPointClassList(m_PointClassList);
//  m_ZYSliceCanvas->setContours(m_Contours);
  setCurrentData(0,0);

  /* set the depth back to original values */
  m_XYSliceCanvas->setDepth(zdepth <= m_VolumeFileInfo.ZDim()-1 ? zdepth : m_VolumeFileInfo.ZDim()-1);
  m_XZSliceCanvas->setDepth(ydepth <= m_VolumeFileInfo.YDim()-1 ? ydepth : m_VolumeFileInfo.YDim()-1);
  m_ZYSliceCanvas->setDepth(xdepth <= m_VolumeFileInfo.XDim()-1 ? xdepth : m_VolumeFileInfo.XDim()-1);

  sliceAxisChangedSlot(); /* make sure the right xyz widget is set to read/write */

  _ui->m_LocalOutputFile->setText(m_VolumeFileInfo.filename().data());

#ifdef USING_EM_CLUSTERING
  m_HistogramCalculated = false;
#endif

  m_hasVolume = true;

  return true;
}

void VolumeGridRover::unsetVolume()
{
  m_XYSliceCanvas->unsetVolume();
  m_XZSliceCanvas->unsetVolume();
  m_ZYSliceCanvas->unsetVolume();
  _ui->m_Variable->clear();
  _ui->m_Timestep->setMinimum(0);
  _ui->m_Timestep->setMaximum(0);
  _ui->m_Timestep->setValue(0);
  _ui->m_PointClass->clear(); /* clear the combo box */
  /* reset the depth sliders */
  _ui->m_XYDepthSlide->setValue(0);
  _ui->m_XZDepthSlide->setValue(0);
  _ui->m_ZYDepthSlide->setValue(0);
  _ui->m_XYDepthSlide->setMinimum(0);
  _ui->m_XYDepthSlide->setMaximum(0);
  _ui->m_XZDepthSlide->setMinimum(0);
  _ui->m_XZDepthSlide->setMaximum(0);
  _ui->m_ZYDepthSlide->setMinimum(0);
  _ui->m_ZYDepthSlide->setMaximum(0);
  _ui->m_X->setReadOnly(true);
  _ui->m_Y->setReadOnly(true);
  _ui->m_Z->setReadOnly(true);
  
  if(m_PointClassList)
    {
      for(unsigned int i=0; i<m_VolumeFileInfo.numVariables(); i++)
	{
	  for(unsigned int j=0; j<m_VolumeFileInfo.numTimesteps(); j++)
	    delete m_PointClassList[i][j];
	  delete [] m_PointClassList[i];
	}
      delete [] m_PointClassList;
    }
  m_PointClassList = NULL;

  m_VolumeFileInfo = VolMagick::VolumeFileInfo();

  //m_Contour->clear(); //clear the combo box
//  m_Objects->clear();
//  m_Contours.clear();
}

bool VolumeGridRover::hasVolume(void)
{
   return m_hasVolume;
}

void VolumeGridRover::getLocalOutputFileSlot()
{
  _ui->m_LocalOutputFile->setText(QFileDialog::getSaveFileName(this,
							  "Choose a filename to save under",
							  QString(m_VolumeFileInfo.filename().data()),
							  "RawIV Volumes (*.rawiv)"));
}

void VolumeGridRover::getRemoteFileSlot()
{
  QMessageBox::information(this,"Notice",
			   "File selection acts only upon the local filesystem. "
			   "Ensure that a filename you select is correct for the remote server.",
			   QMessageBox::Ok);
  _ui->m_RemoteFile->setText(QFileDialog::getOpenFileName(this,
						     "Choose a file to load remotely",
						     QString(),
						     "RawIV Volumes (*.rawiv)"));
}

void VolumeGridRover::localSegmentationRunSlot()
{
  if(!m_VolumeFileInfo.isSet())
    {
      QMessageBox::critical(this,"Error","Load a volume first.",
			    QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
      return;
    }
  
  if(m_LocalGenSegThread.isRunning())
    QMessageBox::critical(this,"Error","Local segmentation thread already running!",
			  QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
  else
    m_LocalGenSegThread.start();
}

void VolumeGridRover::remoteSegmentationRunSlot()
{
  if(!m_VolumeFileInfo.isSet())
    {
      QMessageBox::critical(this,"Error","Load a volume first.",
			    QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
      return;
    }
  
  if(m_RemoteGenSegThread.isRunning())
    QMessageBox::critical(this,"Error","Remote segmentation thread already running!",
			  QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
  else
    m_RemoteGenSegThread.start();
}

void VolumeGridRover::sliceAxisChangedSlot()
{
  switch(_ui->m_SliceCanvasTab->currentIndex())
    {
    default:
    case 0: // XY
      _ui->m_X->setReadOnly(true);
      _ui->m_Y->setReadOnly(true);
      _ui->m_Z->setReadOnly(false);
      break;
    case 1: // XZ
      _ui->m_X->setReadOnly(true);
      _ui->m_Y->setReadOnly(false);
      _ui->m_Z->setReadOnly(true);
      break;
    case 2: // ZY
      _ui->m_X->setReadOnly(false);
      _ui->m_Y->setReadOnly(true);
      _ui->m_Z->setReadOnly(true);
      break;
    }
}

void VolumeGridRover::EMClusteringRunSlot()
{
#ifdef USING_EM_CLUSTERING
  /*
    The EM clustering code currently only works with the first variable and timestep.
  */

  /* re-calculate the histogram */
  if(m_HistogramCalculated == false)
    {
      unsigned int v,t,i,j,k;
      /*for(v=0; v<m_MappedVolumeFile->numVariables(); v++)
	for(t=0; t<m_MappedVolumeFile->numTimesteps(); t++)*/
      printf("VolumeGridRover::setVolume(): Calculating the histogram for var:0, time:0\n");
      QProgressDialog progress("Calculating the histogram for var:0, time:0","Abort",m_VolumeFileInfo.ZDim(),this,"progress",true);
      progress.show();
      v=t=0;
      for(k=0; k<m_VolumeFileInfo.ZDim(); k++)
	{
	  VolMagick::Volume tmpvol;

	  /* read a slice */
	  VolMagick::readVolumeFile(cvcapp, tmpvol,
				    m_VolumeFileInfo.filename(),
				    v,t,
				    0,0,k,
				    VolMagick::Dimension(m_VolumeFileInfo.XDim(),
							 m_VolumeFileInfo.YDim(),1));
     
	  /* If the volume isn't unsigned char, lets convert it because sangmin's library likes uchar volumes */
	  if(tmpvol.voxelType() != VolMagick::UChar)
	    {
	      tmpvol.map(0.0,255.0);
	      tmpvol.voxelType(VolMagick::UChar);
	    }
     
	  for(j=0; j<m_VolumeFileInfo.YDim(); j++)
	    for(i=0; i<m_VolumeFileInfo.XDim(); i++)
	      m_Histogram[(unsigned char)tmpvol(i,j,0)]++; 
	  progress.setProgress(k);
	  qApp->processEvents();
	  if(progress.wasCanceled())
	    return;
	  fprintf(stderr,"%5.2f%%\r",(k/(m_VolumeFileInfo.ZDim()-1))*100.0);
	}
      printf("\n");
      progress.setProgress(m_VolumeFileInfo.ZDim());

      for(i=0; i<256; i++)
	{
	  printf("%d: %d\n",i,m_Histogram[i]);
	}

      /*
	Setup Sangmin's EM clustering stuff.
      */
      m_Clusters.setHistogram(m_Histogram, 1);
      m_Clusters.setData(NULL, 0, 255);

      m_HistogramCalculated = true;
    }

  unsigned int i,j,ithCluster;
  float *Material_prob;
  unsigned char *RangeMin,*RangeMax;
  QValueList<unsigned char> InitialValues;
  VolMagick::Volume vol;

  RangeMin = new unsigned char[m_PointClassList[0][0]->count()];
  RangeMax = new unsigned char[m_PointClassList[0][0]->count()];
  
  /* load the volume for random access */
  VolMagick::readVolumeFile(cvcapp, vol,m_VolumeFileInfo.filename());
  if(vol.voxelType() != VolMagick::UChar)
    {
      vol.map(0.0,255.0);
      vol.voxelType(VolMagick::UChar);
    }
			    
  m_Clusters.InitializeEM(m_PointClassList[0][0]->count());
  for(i=0; i<m_PointClassList[0][0]->count(); i++)
    {
      double avg;
      /* get the average voxel value for point class 'i' */
      if(m_PointClassList[0][0]->at(i)->getPointList().count() == 0) continue;
      for(j=0,avg=0; j<m_PointClassList[0][0]->at(i)->getPointList().count(); j++)
	{
	  GridPoint *gp = m_PointClassList[0][0]->at(i)->getPointList().at(j);
	  avg += vol(gp->x,gp->y,gp->z);
	  //avg += (double)(255.0*((m_MappedVolumeFile->get(0,0)->get(gp->x,gp->y,gp->z) - m_MappedVolumeFile->get(0,0)->m_Min)/
	  //			 (m_MappedVolumeFile->get(0,0)->m_Max - m_MappedVolumeFile->get(0,0)->m_Min)));
	}
      avg /= m_PointClassList[0][0]->at(i)->getPointList().count();
      //m_Clusters.setMeanVariance(i, avg, (double)1200);
      InitialValues.append((unsigned char)avg);
    }
  qHeapSort(InitialValues); /* setMeanVariance expects initial values to be in acending order */
  for(i=0; i<m_PointClassList[0][0]->count(); i++)
    m_Clusters.setMeanVariance(i,InitialValues[i], (double)1200);

  /*
    m_Clusters.setMeanVariance(0,(double)3,(double)1200);
    m_Clusters.setMeanVariance(1,(double)61,(double)1200);
    m_Clusters.setMeanVariance(2,(double)125,(double)1200);
    m_Clusters.setMeanVariance(3,(double)200,(double)1200);
  */

  m_Clusters.iterate();
  Material_prob = m_Clusters.getProbability();

  for(ithCluster=0; ithCluster<m_PointClassList[0][0]->count(); ithCluster++) 
    {
      for(j=0,i=0; j<256; i++)
	{
	  // Find RangeMin
	  for (i=j; i<256; i++)
	    {
	      if (m_Histogram[i]>0)
		{
		  if (Material_prob[i*m_PointClassList[0][0]->count() + ithCluster]>=0.1) { RangeMin[ithCluster] = (unsigned char)i; break; }
		}
	    }
	  i++;
	  // Find RangeMax
	  for(; i<256; i++)
	    {
	      if (Material_prob[i*m_PointClassList[0][0]->count() + ithCluster]<0.1) { RangeMax[ithCluster] = (unsigned char)(i-1); break; }
	    }
	  j=i+1;
	}
    }
  // Considering the final data value
  RangeMax[m_PointClassList[0][0]->count()-1] = 255;

  for(j=0; j<256; j++)
    {
      printf("%d: ",j);
      for(i=0; i<m_PointClassList[0][0]->count(); i++)
    	printf("%f ",Material_prob[j*m_PointClassList[0][0]->count() + i]);
      printf("\n");
    }

  for (ithCluster=0; ithCluster<m_PointClassList[0][0]->count(); ithCluster++)
    {
      printf ("Ith Material = %d ", ithCluster);
      printf ("Min & Max = %3d %3d ", (int)RangeMin[ithCluster], (int)RangeMax[ithCluster]);
      printf ("\n"); fflush(stdout);
    }

  delete RangeMin;
  delete RangeMax;
#else
  QMessageBox::information(this,"Notice","EM clustering is unavailable in this build of VolumeGridRover.",QMessageBox::Ok);
#endif
}

void VolumeGridRover::savePointClassesSlot()
{
  if(!m_VolumeFileInfo.isSet())
    {
      QMessageBox::critical(this,"Error","Load a volume file first!",
			    QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
      return;
    }

  QDomDocument doc("pointclassdoc");
  QDomElement root = doc.createElement("pointclassdoc");
  doc.appendChild(root);

  for(unsigned int i=0; i<m_VolumeFileInfo.numVariables(); i++)
    {
      for(unsigned int j=0; j<m_VolumeFileInfo.numTimesteps(); j++)
	{
          int size = m_PointClassList[i][j]->size();
          for( int k=0; k<size; k++ )
	    {
	      PointClass *pc = m_PointClassList[i][j]->at(k);
	      QDomElement pctag;
	     	
		  pctag = doc.createElement("pointclass");
	      pctag.setAttribute("name", pc->getName());
	      pctag.setAttribute("color", pc->getColor().name());
	      pctag.setAttribute("variable", QString("%1").arg(i));
	      pctag.setAttribute("timestep", QString("%1").arg(j));
	      root.appendChild(pctag);
	
              int sizePoint = pc->getPointList().size();
              for( int l = 0; l<sizePoint; l++ )
                {
	          GridPoint *gp = pc->getPointList().at(l);
		  QDomElement ptag;
		  QDomText t;
		  
		  ptag = doc.createElement("point");
		  pctag.appendChild(ptag);
		  
		  t = doc.createTextNode(QString("%1 %2 %3").arg(gp->x).arg(gp->y).arg(gp->z));
		  ptag.appendChild(t);
		}
	    }
	}
    }

  
  QString filename = QFileDialog::getSaveFileName(this,
						  "Choose a filename to save under",
						  QString(),
						  "Point Class Files (*.vgr);;All Files (*)");
  if(filename.isEmpty())
    {
      cvcapp.log(5, "VolumeGridRover::savePointClassesSlot(): save cancelled (or null filename)");
      return;
    }
  
  QFile f(filename);
  if(!f.open(QIODevice::WriteOnly))
    {
      QMessageBox::critical(this,"Error",QString("Could not open the file %1!").arg(filename),
			    QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
      return;
    }
  
  QTextStream stream(&f);
  stream << doc.toString();
  f.close();
  
  //printf("%s\n",doc.toString().toStdString().c_str());
}

void VolumeGridRover::loadPointClassesSlot()
{
  if(!m_VolumeFileInfo.isSet())
    {
      QMessageBox::critical(this,"Error","Load a volume file first!",
			    QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
      return;
    }
  
  QString filename = QFileDialog::getOpenFileName(this,
						  "Choose a file to open",
						  QString(),
						  "Point Class Files (*.vgr);;Point List Files (*.pts);;All Files (*)");
  if(filename.isEmpty())
    {
      QMessageBox::critical(this,"Error","Filename must be specified!",
			    QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
      return;
    }
  
  /* reset the point class list */
  for(unsigned int i=0; i<m_VolumeFileInfo.numVariables(); i++)
    for(unsigned int j=0; j<m_VolumeFileInfo.numTimesteps(); j++)
      m_PointClassList[i][j]->clear();

  _ui->m_PointClass->clear(); /* clear the combo box */
  
  if(filename.endsWith(".pts"))
    {
      srand(time(NULL));

      QFile f(filename);
      if(!f.open(QIODevice::ReadOnly))
	{
	  QMessageBox::critical(this,"Error",QString("Could not open the file %1!").arg(filename),
				QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
	  return;
	}
      
      /* 
       *.pts coordinate values appear to be in object space rather than in voxel space, so we must convert 
       them to voxel space so they can correctly mark corresponding voxels.
      */
      int numpts, i, j;
      float objpoint[3];
      int point[3], curclass=0;
      PointClass *currentPointClass;
      QString classname;
      QColor classcolor;
      QTextStream stream(&f);
      while(!stream.atEnd())
	{
	  //invent some class name and color since such values are not defined in *.pts files
	  classname = QString("Class %1").arg(curclass++);
	  classcolor = QColor(rand()%255,rand()%255,rand()%255);

	  //now insert the point class according to the current variable and timestep
	  m_PointClassList[_ui->m_Variable->currentIndex()][_ui->m_Timestep->value()]->append(currentPointClass = new PointClass(classcolor,classname));

	  stream >> numpts;
	  for(i=0; i<numpts; i++)
	    {
	      //read the points
	      for(j=0; j<3; j++)
		stream >> objpoint[j];

	      //convert to voxel indices
	      point[0] = int(((objpoint[0] - m_VolumeFileInfo.XMin())/(m_VolumeFileInfo.XMax() - m_VolumeFileInfo.XMin()))*
			     (m_VolumeFileInfo.XDim()-1));
	      point[1] = int(((objpoint[1] - m_VolumeFileInfo.YMin())/(m_VolumeFileInfo.YMax() - m_VolumeFileInfo.YMin()))*
			     (m_VolumeFileInfo.YDim()-1));
	      point[2] = int(((objpoint[2] - m_VolumeFileInfo.ZMin())/(m_VolumeFileInfo.ZMax() - m_VolumeFileInfo.ZMin()))*
			     (m_VolumeFileInfo.ZDim()-1));
	      
	      //make sure that point is a valid voxel index (i.e. within bounds)
	      point[0] = point[0] < 0 ? 0 : point[0] >= int(m_VolumeFileInfo.XDim()) ? m_VolumeFileInfo.XDim()-1 : point[0];
	      point[1] = point[1] < 0 ? 0 : point[1] >= int(m_VolumeFileInfo.YDim()) ? m_VolumeFileInfo.YDim()-1 : point[1];
	      point[2] = point[2] < 0 ? 0 : point[2] >= int(m_VolumeFileInfo.ZDim()) ? m_VolumeFileInfo.ZDim()-1 : point[2];

	      currentPointClass->addPoint(point[0],point[1],point[2]);
	    }
	}
    }
  else // interpret the file as *.vgr by default
    {
      // TODO: Qt6 migration - QXmlSimpleReader/QXmlDefaultHandler removed
      // Need to rewrite using QXmlStreamReader (pull-based API)
      QMessageBox::warning(this, "Not Implemented", 
                          "VGR file loading requires Qt6 XML migration.\n"
                          "Please use .pts files for now.",
                          QMessageBox::Ok);
      /*
      QFile f(filename);
      QXmlInputSource qxis(f);
      QXmlSimpleReader xsr;
      PointClassFileContentHandler pcfch(this);
      xsr.setContentHandler(&pcfch);
      if(!xsr.parse(&qxis,false))
	cvcapp.log(5, boost::str(boost::format("VolumeGridRover::loadPointClassesSlot(): Parse error: %s")%pcfch.errorString().toStdString().c_str()));
      */
    }
  
  m_XYSliceCanvas->update();
  m_XZSliceCanvas->update();
  m_ZYSliceCanvas->update();

  setCurrentData(_ui->m_Variable->currentIndex(),_ui->m_Timestep->value()); //reset the point class combo box and load the new values for it
}

void VolumeGridRover::setCurrentVariable(int variable)
{
  setCurrentData(variable,_ui->m_Timestep->value());
}

void VolumeGridRover::setCurrentTimestep(int timestep)
{
  setCurrentData(_ui->m_Variable->currentIndex(),timestep);
}

void VolumeGridRover::setCurrentData(int variable, int timestep)
{
  unsigned int i;

  setMinValue(m_VolumeFileInfo.min(variable,timestep));
  setMaxValue(m_VolumeFileInfo.max(variable,timestep));

  _ui->m_PointClass->clear();

  int size = m_PointClassList[variable][timestep]->size();
  for(i=0; i<size; i++)
    _ui->m_PointClass->insertItem(i, QString(m_PointClassList[variable][timestep]->at(i)->getName()));

  _ui->m_PointClass->setCurrentIndex(0);
/*
  //m_Contour->clear();
  m_Objects->clear();
  for(SurfRecon::ContourPtrMap::iterator i = m_Contours[m_Variable->currentIndex()][m_Timestep->value()].begin();
      i != m_Contours[m_Variable->currentIndex()][m_Timestep->value()].end();
      i++)
    {
      if(i->second.get() == NULL) continue;
      //m_Contour->insertItem(m_Contours[variable][timestep][i]->name());
      m_Objects->insertItem(new QListViewItem(m_Objects,(*i).second->name()));
      //m_Contour->setCurrentItem(0);
    }
*/

  //update the viewer objects
  m_XYSliceCanvas->setCurrentClass(0);
  m_XZSliceCanvas->setCurrentClass(0);
  m_ZYSliceCanvas->setCurrentClass(0);
  showColor(0);
/*
  m_XYSliceCanvas->setCurrentContour(m_Objects->currentItem() ? m_Objects->currentItem()->text(0).toStdString().c_str() : "");
  m_XZSliceCanvas->setCurrentContour(m_Objects->currentItem() ? m_Objects->currentItem()->text(0).toStdString().c_str() : "");
  m_ZYSliceCanvas->setCurrentContour(m_Objects->currentItem() ? m_Objects->currentItem()->text(0).toStdString().c_str() : "");
  showContourColor(m_Objects->currentItem() ? m_Objects->currentItem()->text(0).toStdString().c_str() : "");
  showContourInterpolationType(m_Objects->currentItem() ? m_Objects->currentItem()->text(0).toStdString().c_str() : "");
  showContourInterpolationSampling(m_Objects->currentItem() ? m_Objects->currentItem()->text(0).toStdString().c_str() : "");
*/
}

void VolumeGridRover::setX(int x)
{
  _ui->m_X->setText(QString("%1").arg(x));
}

void VolumeGridRover::setY(int y)
{
  _ui->m_Y->setText(QString("%1").arg(y));
}

void VolumeGridRover::setZ(int z)
{
  _ui->m_Z->setText(QString("%1").arg(z));
}

void VolumeGridRover::setXYZ(int x, int y, int z)
{
  setX(x);
  setY(y);
  setZ(z);
}

void VolumeGridRover::setObjX(float objx)
{
  _ui->m_ObjX->setText(QString("%1").arg(objx));
}

void VolumeGridRover::setObjY(float objy)
{
  _ui->m_ObjY->setText(QString("%1").arg(objy));
}

void VolumeGridRover::setObjZ(float objz)
{
  _ui->m_ObjZ->setText(QString("%1").arg(objz));
}

void VolumeGridRover::setObjXYZ(float objx, float objy, float objz)
{
  setObjX(objx);
  setObjY(objy);
  setObjZ(objz);
}

void VolumeGridRover::setR(int r)
{
  _ui->m_R->setText(QString("%1").arg(r));
  setColorName();
}

void VolumeGridRover::setG(int g)
{
  _ui->m_G->setText(QString("%1").arg(g));
  setColorName();
}

void VolumeGridRover::setB(int b)
{
  _ui->m_B->setText(QString("%1").arg(b));
  setColorName();
}

void VolumeGridRover::setA(int a)
{
  _ui->m_A->setText(QString("%1").arg(a));
}

void VolumeGridRover::setRGBA(int r, int g, int b, int a)
{
  setR(r); setG(g); setB(b); setA(a);
}

void VolumeGridRover::setValue(double value)
{
  _ui->m_Value->setText(QString("%1").arg(value));
  if(m_VolumeFileInfo.isSet() && m_VolumeFileInfo.voxelTypes(_ui->m_Variable->currentIndex()) != VolMagick::UChar)
    {
      _ui->m_MappedValue->setText(QString("%1")
			     .arg((unsigned char)(255.0*
						  ((value - m_VolumeFileInfo.min(_ui->m_Variable->currentIndex(),_ui->m_Timestep->value()))/
						   (m_VolumeFileInfo.max(_ui->m_Variable->currentIndex(),_ui->m_Timestep->value()) - 
						    m_VolumeFileInfo.min(_ui->m_Variable->currentIndex(),_ui->m_Timestep->value()))))));
    }
  else
    _ui->m_MappedValue->setText(QString("%1").arg(value));
}

void VolumeGridRover::setMinValue(double value)
{
  _ui->m_MinimumValue->setText(QString("%1").arg(value));
}

void VolumeGridRover::setMaxValue(double value)
{
  _ui->m_MaximumValue->setText(QString("%1").arg(value));
}

void VolumeGridRover::setGridCellInfo(int x, int y, int z,
				      float objx, float objy, float objz,
				      int r, int g, int b, int a, double value)
{
  setXYZ(x,y,z); setObjXYZ(objx,objy,objz); setRGBA(r,g,b,a); setValue(value);
}

void VolumeGridRover::showColor(int i)
{
  // if there are no point classes, do nothig and return without a warning message
  if (m_PointClassList[_ui->m_Variable->currentIndex()][_ui->m_Timestep->value()]->size() == 0) {
    return;
  }


  if( m_PointClassList[_ui->m_Variable->currentIndex()][_ui->m_Timestep->value()]->size() <= i ) {
     fprintf( stderr, "VolumeGridRover::showColor(int i): i is larger than point class list size\n");
     return;
  }

  if(m_PointClassList[_ui->m_Variable->currentIndex()][_ui->m_Timestep->value()]->at(i))
    _ui->m_PointClassColor->setPalette( QPalette(m_PointClassList[_ui->m_Variable->currentIndex()][_ui->m_Timestep->value()]->at(i)->getColor()) );
}

void VolumeGridRover::showContourColor(const QString& name)
{
/*
  if(m_Contours[m_Variable->currentIndex()][m_Timestep->value()][name.toStdString().c_str()] != NULL)
    m_ContourColor->setPaletteBackgroundColor(QColor(m_Contours[m_Variable->currentIndex()][m_Timestep->value()][name.toStdString().c_str()]->color().get<0>()*255.0,
						     m_Contours[m_Variable->currentIndex()][m_Timestep->value()][name.toStdString().c_str()]->color().get<1>()*255.0,
						     m_Contours[m_Variable->currentIndex()][m_Timestep->value()][name.toStdString().c_str()]->color().get<2>()*255.0));
*/
}

void VolumeGridRover::showContourInterpolationType(const QString& name)
{
/*
  if(m_Contours[m_Variable->currentIndex()][m_Timestep->value()][name.toStdString().c_str()] != NULL)
    m_InterpolationType->setCurrentItem(int(m_Contours[m_Variable->currentIndex()][m_Timestep->value()][name.toStdString().c_str()]->interpolationType()));
*/
}

void VolumeGridRover::showContourInterpolationSampling(const QString& name)
{
/*
   if(m_Contours[m_Variable->currentIndex()][m_Timestep->value()][name.toStdString().c_str()] != NULL)
     m_InterpolationSampling->setValue(m_Contours[m_Variable->currentIndex()][m_Timestep->value()][name.toStdString().c_str()]->numberOfSamples());
*/
}

void VolumeGridRover::currentObjectSelectionChanged()
{
/*
  m_XYSliceCanvas->setCurrentContour(lvi ? lvi->text(0).toStdString().c_str() : "");
  m_XZSliceCanvas->setCurrentContour(lvi ? lvi->text(0).toStdString().c_str() : "");
  m_ZYSliceCanvas->setCurrentContour(lvi ? lvi->text(0).toStdString().c_str() : "");
  showContourColor(lvi ? lvi->text(0).toStdString().c_str() : "");
  showContourInterpolationType(lvi ? lvi->text(0).toStdString().c_str() : "");
  showContourInterpolationSampling(lvi ? lvi->text(0).toStdString().c_str() : "");
*/
}

void VolumeGridRover::cellMarkingModeTabChangedSlot(QWidget *w)
{
  emit cellMarkingModeChanged(_ui->m_GridCellMarkingToolSelection->currentIndex());
}

void VolumeGridRover::setCellMarkingMode(int m)
{
  m_CellMarkingMode = CellMarkingMode(m);
  m_XYSliceCanvas->setCellMarkingMode(m);
  m_XZSliceCanvas->setCellMarkingMode(m);
  m_ZYSliceCanvas->setCellMarkingMode(m);
  //cvcapp.log(5, "VolumeGridRover::setCellMarkingMode(): mode changed to %d",m);
}

void VolumeGridRover::hideEvent(QHideEvent *e)
{
  //cvcapp.log(5, "VolumeGridRover::hideEvent(): %p\n",e);
  emit showToggle(false);
}

void VolumeGridRover::showEvent(QShowEvent *e)
{
  //cvcapp.log(5, "VolumeGridRover::showEvent(): %p\n",e);
  emit showToggle(true);
}

void VolumeGridRover::setColorName()
{
  _ui->m_ColorName->setText(QColor(_ui->m_R->text().toInt(),_ui->m_G->text().toInt(),_ui->m_B->text().toInt()).name());
}

void VolumeGridRover::setPointSize(int r)
{
  m_XYSliceCanvas->setPointSize(r);
  m_XZSliceCanvas->setPointSize(r);
  m_ZYSliceCanvas->setPointSize(r);
}

void VolumeGridRover::ZDepthChangedSlot(int d)
{
  emit depthChanged(m_XYSliceCanvas, d);
}

void VolumeGridRover::YDepthChangedSlot(int d)
{
  emit depthChanged(m_XZSliceCanvas, d);
}

void VolumeGridRover::XDepthChangedSlot(int d)
{
  emit depthChanged(m_ZYSliceCanvas, d);
}

void VolumeGridRover::setSelectedContours()
{
/*
  QListViewItemIterator it(m_Objects, QListViewItemIterator::Selected);
  QListViewItem *tmp;
  while((tmp=it.current())!=NULL)
    {
      ++it;
      m_Contours[m_Variable->currentIndex()][m_Timestep->value()][tmp->text(0).toStdString().c_str()]->selected(true);
    }

  QListViewItemIterator it2(m_Objects, QListViewItemIterator::Unselected);
  while((tmp=it2.current())!=NULL)
    {
      ++it2;
      m_Contours[m_Variable->currentIndex()][m_Timestep->value()][tmp->text(0).toStdString().c_str()]->selected(false);
    }

  m_XYSliceCanvas->update();
  m_XZSliceCanvas->update();
  m_ZYSliceCanvas->update();
*/
}

/***** segmentation stuff (Similar to the segmentation code in VolumeProperties.cpp) ******/

/* Need these to signal the GUI thread the segmentation result so it can popup a message to the user */
class SegmentationFailedEvent : public QEvent
{
public:
#ifdef WIN32
  SegmentationFailedEvent(const QString &m) : QEvent(QEvent::Type(1100)), msg(m) {}
#else
  SegmentationFailedEvent(const QString &m) : QEvent(QEvent::Type(QEvent::User+100)), msg(m) {}
#endif
  QString message() const { return msg; }
private:
  QString msg;
};

class SegmentationFinishedEvent : public QEvent
{
public:
#ifdef WIN32
  SegmentationFinishedEvent(const QString &m,
			    const QString& newvol = QString())
    : QEvent(QEvent::Type(1101)), _msg(m), _newvol(newvol) {}
#else
  SegmentationFinishedEvent(const QString &m,
			    const QString& newvol = QString())
    : QEvent(QEvent::Type(QEvent::User+101)), _msg(m), _newvol(newvol) {}
#endif
  QString message() const { return _msg; }
  QString outputVolumeFilename() const { return _newvol; }
private:
  QString _msg;
  QString _newvol;
};

class TilingFailedEvent : public QEvent
{
public:
#ifdef WIN32
  TilingFailedEvent(const QString &m) : QEvent(QEvent::Type(1102)), msg(m) {}
#else
  TilingFailedEvent(const QString &m) : QEvent(QEvent::Type(QEvent::User+102)), msg(m) {}
#endif
  QString message() const { return msg; }
private:
  QString msg;
};

class TilingFinishedEvent : public QEvent
{
public:
#ifdef WIN32
  TilingFinishedEvent(const QString &m) : QEvent(QEvent::Type(1103)), msg(m) {}
#else
  TilingFinishedEvent(const QString &m) : QEvent(QEvent::Type(QEvent::User+103)), msg(m) {}
#endif
  QString message() const { return msg; }
private:
  QString msg;
};

VolumeGridRover::RemoteGenSegThread::RemoteGenSegThread(VolumeGridRover *vgr, unsigned int stackSize)
  : QThread(vgr), m_VolumeGridRover(vgr)
{}

void VolumeGridRover::RemoteGenSegThread::run()
{
  unsigned int i,j,index=0;
  // bool ok;
 
  if(m_VolumeGridRover->_ui->m_Hostname->text().isEmpty() || m_VolumeGridRover->_ui->m_Port->text().isEmpty() ||
     m_VolumeGridRover->_ui->m_ThresholdLow->text().isEmpty() || m_VolumeGridRover->_ui->m_ThresholdHigh->text().isEmpty())
    {
      QApplication::postEvent(m_VolumeGridRover,new SegmentationFailedEvent("Values for Hostname, Port, Threshold low and Threshold high must be set!"));
      return;
    }

  if(m_VolumeGridRover->_ui->m_RemoteFile->text().isEmpty())
    {
      QApplication::postEvent(m_VolumeGridRover,new SegmentationFailedEvent("Please specify a remote filename."));
      return;
    }
 
  XmlRpc::setVerbosity(5);
  XmlRpcClient c(m_VolumeGridRover->_ui->m_Hostname->text().toStdString().c_str(), m_VolumeGridRover->_ui->m_Port->text().toInt());
  XmlRpcValue args, result;
 
  cvcapp.log(5, boost::str(boost::format("VolumeGridRover::remoteSegmentationRunSlot(): segmenting '%s', loading '%s' remotely.")%m_VolumeGridRover->m_VolumeFileInfo.filename().c_str()%m_VolumeGridRover->_ui->m_RemoteFile->text().toStdString().c_str()));
  /*
    GenSegmentation( filename, threshold low, threshold high, number of seed point classes
    num of points in class1, class1 seed points x y z ... ,
    num of points in class2, class2 seed points x y z ... )
  */
  args[0] = m_VolumeGridRover->_ui->m_RemoteFile->text().toStdString().c_str();
  args[1] = m_VolumeGridRover->_ui->m_ThresholdLow->text().toInt();
  args[2] = m_VolumeGridRover->_ui->m_ThresholdHigh->text().toInt();
  args[3] = int(m_VolumeGridRover->m_PointClassList[m_VolumeGridRover->_ui->m_Variable->currentIndex()][m_VolumeGridRover->_ui->m_Timestep->value()]->count());
  args[4] = m_VolumeGridRover->_ui->m_RemoteFile->text().toStdString().c_str();
  index = 5;
  for(i=0; i<m_VolumeGridRover->m_PointClassList[m_VolumeGridRover->_ui->m_Variable->currentIndex()][m_VolumeGridRover->_ui->m_Timestep->value()]->count(); i++)
    {
#ifdef DEBUG
      cvcapp.log(5, boost::str(boost::format("%s")%m_VolumeGridRover->m_PointClassList[m_VolumeGridRover->_ui->m_Variable->currentIndex()][m_VolumeGridRover->_ui->m_Timestep->value()]->at(i)->getName().toStdString().c_str()));
#endif
      QList<GridPoint*> points(m_VolumeGridRover->m_PointClassList[m_VolumeGridRover->_ui->m_Variable->currentIndex()][m_VolumeGridRover->_ui->m_Timestep->value()]->at(i)->getPointList());
#ifdef DEBUG
      cvcapp.log(5, boost::str(boost::format("Num points: %d")%points.count());
#endif
      args[index] = int(points.count());
      index++;
      for(j=0; j<points.count(); j++)
	{
	  args[index+0] = int(points.at(j)->x);
	  args[index+1] = int(points.at(j)->y);
	  args[index+2] = int(points.at(j)->z);
#ifdef DEBUG
	  printf("(%d, %d, %d) ",points.at(j)->x,points.at(j)->y,points.at(j)->z);
#endif
	  if(j % 10 == 0) printf("\n");
	  index+=3;
	}
#ifdef DEBUG
      printf("\n");
#endif
    }
 
  if(!c.execute("GenSegmentation", args, result))
    QApplication::postEvent(m_VolumeGridRover,new SegmentationFailedEvent("Error calling GenSegmentation on remote server!"));
  else
    QApplication::postEvent(m_VolumeGridRover,new SegmentationFinishedEvent("Remote general segmentation finished."));
}

VolumeGridRover::LocalGenSegThread::LocalGenSegThread(VolumeGridRover *vgr, unsigned int stackSize)
  : QThread(vgr), m_VolumeGridRover(vgr)
{}

void VolumeGridRover::LocalGenSegThread::run()
{
  unsigned int i,j,index=0;
  //  bool ok;
  
  if(m_VolumeGridRover->_ui->m_ThresholdLow->text().isEmpty() || m_VolumeGridRover->_ui->m_ThresholdHigh->text().isEmpty())
    {
      QApplication::postEvent(m_VolumeGridRover,new SegmentationFailedEvent("Values for Threshold low and Threshold high must be set!"));
      return;
    }
  if(m_VolumeGridRover->_ui->m_LocalOutputFile->text().isEmpty())
    {
      QApplication::postEvent(m_VolumeGridRover,new SegmentationFailedEvent("Please specify an output filename."));
      return;
    }
  
  XmlRpc::setVerbosity(5);
  XmlRpcValue args, result;
  
  cvcapp.log(5, boost::str(boost::format("VolumeGridRover::localSegmentationRunSlot(): segmenting '%s', loading '%s' locally.")%m_VolumeGridRover->m_VolumeFileInfo.filename().c_str()%m_VolumeGridRover->_ui->m_RemoteFile->text().toStdString().c_str()));
  /*
    GenSegmentation( filename, threshold low, threshold high, number of seed point classes
    num of points in class1, class1 seed points x y z ... ,
    num of points in class2, class2 seed points x y z ... )
  */
  args[0] = m_VolumeGridRover->_ui->m_LocalOutputFile->text().toStdString().c_str();
  args[1] = m_VolumeGridRover->_ui->m_ThresholdLow->text().toInt();
  args[2] = m_VolumeGridRover->_ui->m_ThresholdHigh->text().toInt();
  args[3] = int(m_VolumeGridRover->m_PointClassList[m_VolumeGridRover->_ui->m_Variable->currentIndex()][m_VolumeGridRover->_ui->m_Timestep->value()]->count());
  //args[4] = std::string(QString(m_VolumeGridRover->cacheDir() + "/tmp/tmp.rawiv").toStdString().c_str());
  args[4] = std::string(QString("SegTmp.rawiv").toStdString().c_str());
  index = 5;
  for(i=0; i<m_VolumeGridRover->m_PointClassList[m_VolumeGridRover->_ui->m_Variable->currentIndex()][m_VolumeGridRover->_ui->m_Timestep->value()]->count(); i++)
    {
#ifdef DEBUG
      cvcapp.log(5, boost::str(boost::format("%s")%m_VolumeGridRover->m_PointClassList[m_VolumeGridRover->_ui->m_Variable->currentIndex()][m_VolumeGridRover->_ui->m_Timestep->value()]->at(i)->getName().toStdString().c_str()));
#endif
      QList<GridPoint*> points(m_VolumeGridRover->m_PointClassList[m_VolumeGridRover->_ui->m_Variable->currentIndex()][m_VolumeGridRover->_ui->m_Timestep->value()]->at(i)->getPointList());
#ifdef DEBUG
      cvcapp.log(5, boost::str(boost::format("Num points: %d")%points.count()));
#endif
      args[index] = int(points.count());
      index++;
      for(j=0; j<points.count(); j++)
	{
	  args[index+0] = int(points.at(j)->x);
	  args[index+1] = int(points.at(j)->y);
	  args[index+2] = int(points.at(j)->z);
#ifdef DEBUG
	  printf("(%d, %d, %d) ",points.at(j)->x,points.at(j)->y,points.at(j)->z);
#endif
	  if(j % 10 == 0) printf("\n");
	  index+=3;
	}
#ifdef DEBUG
      printf("\n");
#endif
    }
 
  if(!generalSegmentation(args, result))
    QApplication::postEvent(m_VolumeGridRover,new SegmentationFailedEvent("Error running local general segmentation!"));
  else
    {
      using namespace std;
      using namespace boost;
      using namespace boost::algorithm;

      //successfully ran general segmentation.. now convert output to a rawv volume and load that and remove the temporary subunit files
      vector<VolMagick::Volume> volumes;
      QDir d(".");
      regex expression("(.*_subunit_\\d*\\.rawiv)");
      cmatch what;
      int point_class_idx = 0;
      for(int i = 0; i < d.count(); i++)
	{
	  if(regex_match(d[i].toStdString().c_str(),what,expression))
	    {
	      //cvcapp.log(5, "Found: %s",d.absFilePath(d[i]).toStdString().c_str());
              std::cout << d.absoluteFilePath(d[i]).toStdString() << std::endl;

	      VolMagick::Volume vol;
	      VolMagick::readVolumeFile(cvcapp, vol,d.absoluteFilePath(d[i]).toStdString().c_str());
	      //set each volume's description to whatever the point class was named in the VGR UI
	      vol.desc(m_VolumeGridRover->_ui->m_PointClass->itemText(point_class_idx).toStdString().c_str());
	      volumes.push_back(vol);
	      boost::filesystem::remove(d.absoluteFilePath(d[i]).toStdString().c_str());
	      point_class_idx++;
	    }
	}

      //QString newvol_filename(m_VolumeGridRover->cacheDir() + "/tmp/tmp.rawv");
      QString newvol_filename("SegTmp.rawv");
      VolMagick::writeVolumeFile(cvcapp, volumes,string(newvol_filename.toStdString().c_str()));

      QApplication::postEvent(m_VolumeGridRover,
			      new SegmentationFinishedEvent("Local general segmentation finished.",newvol_filename));
      std::cout << "Local general segmentation finished: " << newvol_filename.toStdString() << std::endl;
    }
}


VolumeGridRover::TilingThread::TilingThread(VolumeGridRover *vgr, unsigned int stackSize)
  : QThread(vgr), m_VolumeGridRover(vgr)
{}

#ifdef USING_TILING
// typedef boost::tuple<std::string, boost::shared_ptr<Geometry> > GeometryPkg; //couple a string with a geometry object
#endif

void VolumeGridRover::TilingThread::run()
{
#ifdef USING_TILING
  // QSettings settings;
  // settings.insertSearchPath(QSettings::Windows, "/CCV");
  
  // bool result;
  // QString dirString = settings.readEntry("/Volume Rover/CacheDir", QString(), &result);
  // if(!result) dirString = ".";
  // QDir dir(dirString);
  // dir.cd("VolumeCache"); //write intermediate files to VolumeCache dir

  // QValueList<GeometryPkg> *geometries;
  // for(std::list<std::string>::iterator i = m_SelectedNames.begin();
  //     i != m_SelectedNames.end();
  //     i++)
  //   {
  //     printf("Tiling %s ...\n",i->c_str());
  //     geometries = new QValueList<GeometryPkg>();
// #ifdef USING_TILING
//       geometries->push_back(GeometryPkg(*i,
// 					Tiling::surfaceFromContour(m_VolumeGridRover->
// 								   m_Contours[m_VolumeGridRover->m_Variable->currentIndex()][m_VolumeGridRover->m_Timestep->value()][*i],
// 								   m_VolumeGridRover->m_VolumeFileInfo,
// 								   dir.absPath())));
// #endif
//       QApplication::postEvent(m_VolumeGridRover,
// 			      new QCustomEvent(QEvent::Type(QEvent::User+104),geometries));
//     }

//   QApplication::postEvent(m_VolumeGridRover,
// 			  new TilingFinishedEvent("Tiling complete"));
// #else
  QApplication::postEvent(m_VolumeGridRover,new TilingFailedEvent("Tiling not enabled."));
#endif
}


void VolumeGridRover::customEvent(QEvent *e)
{
  if(e->type() == QEvent::User+100) /* SegmentationFailedEvent */
    QMessageBox::critical(this,"Error",static_cast<SegmentationFailedEvent*>(e)->message(),QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
  else if(e->type() == QEvent::User+101) /* SegmentationFinishedEvent */
    {
      SegmentationFinishedEvent *event = static_cast<SegmentationFinishedEvent*>(e);
      if(event->outputVolumeFilename() != QString())
	emit volumeGenerated(event->outputVolumeFilename());
      QMessageBox::information(this,"Notice",event->message(),QMessageBox::Ok);
    }
#ifdef USING_TILING
  // else if(e->type() == QEvent::User+102) /* TilingFailedEvent */
  //   QMessageBox::critical(this,"Error",static_cast<TilingFailedEvent*>(e)->message(),QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
  // else if(e->type() == QEvent::User+103) /* TilingFinishedEvent */
  //   QMessageBox::information(this,"Notice",static_cast<TilingFinishedEvent*>(e)->message(),QMessageBox::Ok);
  // else if(e->type() == QEvent::User+104) /* geometry data from tiler */
  //   {
  //     GeometryLoader loader;

  //     cvcapp.log(5, "Notice: received tiled geometry");
  //     QValueList<GeometryPkg>::iterator it;
  //     QValueList<GeometryPkg> *geometries = static_cast<QValueList<GeometryPkg> *>(e->data());
  //     switch(m_TilingOutputOptions->selectedId())
  // 	{
  // 	case 0: /* keep in-core */
  // 	default:
  // 	  for(it = geometries->begin(); it != geometries->end(); ++it)
  // 	    emit tilingComplete((*it).get<1>());
  // 	  break;
  // 	case 1: /* write to file */
  // 	  {
  // 	    QDir outputDir(m_TilingOutputDirectory->text());
	    
  // 	    for(it = geometries->begin(); it != geometries->end(); ++it)
  // 	      {
  // 		//if for some reason the object's name is a whole path, then reduce it down to just a filename
  // 		std::string filename(QFileInfo(QString((*it).get<0>().c_str())).fileName().toStdString().c_str());
  // 		if(!loader.saveFile(outputDir.absFilePath(filename + ".rawc").toStdString().c_str(),
  // 				    "Rawc files (*.rawc)",
  // 				    (*it).get<1>().get()))
  // 		  cvcapp.log(5, "Error: Could not save %s!",(filename + ".rawc").c_str());
  // 		else
  // 		  cvcapp.log(5, "Notice: Saved %s",(filename + ".rawc").c_str());
  // 	      }
  // 	  }
  // 	  break;
  // 	}
  //     delete geometries;
  //   }
#endif
}

