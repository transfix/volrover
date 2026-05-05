/*
  Copyright 2012 The University of Texas at Austin
  
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

/* $Id: VolumeViewer2.h 5883 2012-07-20 19:52:38Z transfix $ */

#ifndef __VOLUMEVIEWER2_H__
#define __VOLUMEVIEWER2_H__

#include <cvc/state_object.h>
#include <cvc_compat.h>
#include <cvc_compat.h>
#include <cvc_compat.h>
#include <VolumeRenderer/VolumeRenderer.h>
#include <QGLViewer/qglviewer.h>
#include <string>

#include <GeometryRenderer/GeometryRenderer.h>

namespace CVC_NAMESPACE
{
  // ------------------
  // CVC::VolumeViewer2
  // ------------------
  // Purpose: 
  //  Updated version of VolumeViewer that uses cvcstate.
  // ---- Change History ----
  // 03/30/2012 -- Joe R. -- Started.
  // 06/15/2012 -- Joe R. -- Using state_object<>
  class VolumeViewer2 : public QGLViewer, public state_object<VolumeViewer2>
  {
  Q_OBJECT

  public:

    enum VolumeRenderingType { ColorMapped, RGBA };

    explicit VolumeViewer2(QWidget* parent=nullptr, 
                          Qt::WindowFlags flags=Qt::WindowFlags())
      : QGLViewer(parent, flags),
        state_object<VolumeViewer2>(volrover_app_instance())
    { defaultConstructor(); }

    ~VolumeViewer2();

    bool glUpdated() const;
    bool hasBeenInitialized() const;

    //void setDefaultVolumes();
    virtual void setDefaultColorTable();
    virtual void setDefaultScene();
    virtual void normalizeScene(bool forceShowAll = false);

    //get/set vector of volumes
    //get/set color table
    //get/set subvolume selector
    //get/set vector of geometries

    //main volume rendering bounding box

    //volume rendering type

    //has been initialized flag
    //using VBO flag
    //draw bounding box flag
    //draw subvolume selector flag
    //draw corner axis flag
    //draw geometry flag
    //draw volumes flag
    //clip geometry to bounding box flag
    //shaded rendering flag
    
    //draw geometry normals flag (consider doing this per-geometry)

    //render quality (double)
    //near plane (double)

    BoundingBox globalBoundingBox() const;

  public slots:
    void reset(); //resets the viewer to its default state

    //Use this instead of updateGL() if you want to trigger a redraw
    //but do not need it immediately.  This is safe to call from other
    //threads as well.
    void scheduleUpdateGL();

    void resetROI(); //resets the ROI to the default based on the global bbox

  signals:
    void postInit(); //emitted after init() was called

  protected:
    virtual void init();
    virtual void draw();
    virtual void drawWithNames();
    virtual void postDraw();
    virtual void doDrawCornerAxis();
    virtual void renderStates(const std::string& s);
    virtual void renderStatesWithNames(const std::string& s);
    virtual void calculateGlobalBoundingBox() const;
    virtual void endSelection(const QPoint&);
    virtual void postSelection(const QPoint& point);
    virtual void setClipPlanes();

    static void combineBoundingBox(BoundingBox& bbox, const std::string& s);
    
    enum HandleName
      {
	MaxXHead = 1, MaxXShaft,
	MinXHead, MinXShaft,
	MaxYHead, MaxYShaft,
	MinYHead, MinYShaft,
	MaxZHead, MaxZShaft,
	MinZHead, MinZShaft,
      };
    void doDrawBoundingBox(const std::string& s);
    void doDrawBoundingBox(const BoundingBox& bbox);
    void doDrawSubVolumeSelector(const std::string& s, bool withNames = false);
    void drawHandle(const qglviewer::Vec& from, const qglviewer::Vec& to, 
		    float radius = -1.0, int nbSubdivisions = 12, 
		    int headName = 0, int shaftName = 0);
    void drawHandle(const qglviewer::Vec& from, const qglviewer::Vec& to,
		    float r = 1.0f, float g = 1.0f, float b = 1.0f,
		    int headName = 0, int shaftName = 0,
		    float pointSize = 2.0f, float lineWidth = 2.0f);

    virtual void customEvent(QEvent *event);

    // Mouse events functions
    virtual void mousePressEvent(QMouseEvent *e);
    virtual void mouseMoveEvent(QMouseEvent *e);
    virtual void mouseReleaseEvent(QMouseEvent *e);
    virtual void wheelEvent(QWheelEvent* e);

    //responding to state changes
    virtual void handleStateChanged(const std::string&);

    //responding to thread changes
    void threadsChanged(const std::string&);
    void handleThreadsChanged(const std::string&);
    boost::signals2::connection _threadsConnection;

    virtual void uploadColorTable();
    virtual void uploadVolume(const std::string& s);
    virtual void updateVBO();

    void copyCameraToState();

    VolumeRenderer _renderer; //the volume renderer
    GeometryRenderer _geometryRenderer;

    boost::mutex _handleStateChangedMutex;

  private:
    void defaultConstructor();
  };
}

#endif
