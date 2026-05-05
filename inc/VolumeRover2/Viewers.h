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

#ifndef __CVC_VIEWERS_H__
#define __CVC_VIEWERS_H__

#include <cvc_compat.h>
#include <VolumeRover2/VolumeViewer.h>
#include <ColorTable2/ColorTable.h>

#include <QString>
#include <QWidget>

#include <cvcraw_geometry/cvcgeom.h>
#include <VolMagick/VolMagick.h>

#include <string>
#include <map>

class QSlider;
class QTimer;
namespace Ui
{
  class VolumeViewerPage;
  class VolumeViewerPageManipulators;
}

#ifdef USING_VOLUMEGRIDROVER
class VolumeGridRover;
#endif

namespace CVC_NAMESPACE
{

  class Viewers : public QWidget
  {
    Q_OBJECT

   public:
    Viewers(QWidget *parent = 0,
                  Qt::WindowFlags flags = Qt::WindowFlags());
    virtual ~Viewers();

    VolumeViewer* thumbnailViewer()
      { return _thumbnailViewer; }
    VolumeViewer* subvolumeViewer()
      { return _subvolumeViewer; }

    CVCColorTable::ColorTable* colorTable()
      { return _colorTable; }

#ifdef USING_VOLUMEGRIDROVER
    void setVolumeGridRoverPtr( VolumeGridRover *ptr );
#endif

    //used to easily append a property sub-key to this object's name
    std::string getObjectName(const std::string& property = std::string()) const;

  public slots:
    virtual void markThumbnailDirty(bool flag = true);
    virtual void markSubVolumeDirty(bool flag = true);

    virtual void loadThumbnail();
    virtual void loadSubVolume();

  protected slots:
    virtual void setDefaultSubVolumeViewerState();
    virtual void setDefaultThumbnailViewerState();
    virtual void setThumbnailQuality(int);
    virtual void setSubVolumeQuality(int);
    virtual void setThumbnailNearPlane(int);
    virtual void setSubVolumeNearPlane(int);
    virtual void setDefaultScene();
    virtual void setDefaultOptions();
    virtual void timeout();
    virtual void updateColorTable();
    virtual void saveImage(int);

  protected:
    void ensureVolumeAvailability();
    void ensureSubVolumeAvailability();
    void syncViewers(const std::string& key);

    void customEvent(QEvent *event);
    void propertiesChanged(const std::string&);
    void handlePropertiesChanged(const std::string&);
    boost::signals2::connection _propertiesConnection;
    void dataChanged(const std::string&);
    void handleDataChanged(const std::string&);
    boost::signals2::connection _dataConnection;

    bool _thumbnailVolumeDirty;
    bool _subVolumeDirty;

    VolumeViewer *_thumbnailViewer;
    VolumeViewer *_subvolumeViewer;

    std::map<std::string,VolumeViewer*> _viewerMap;

    QSlider *_subvolumeRenderQualitySlider;
    QSlider *_subvolumeNearClipPlaneSlider;
    QSlider *_thumbnailRenderQualitySlider;
    QSlider *_thumbnailNearClipPlaneSlider;

    CVCColorTable::ColorTable *_colorTable;
    CVCColorTable::ColorTable::isocontour_nodes _oldNodes;

    QTimer *_timer; //used to check if we need to update volumes

    Ui::VolumeViewerPage *_ui;
    Ui::VolumeViewerPageManipulators *_uiManip;

    QString _oldObjectName;

#ifdef USING_VOLUMEGRIDROVER
    VolumeGridRover *_volumeGridRoverPtr;
#endif

    //flag for first time initialization of this object on the property map
    bool _defaultSceneSet;

    //flags to help sync the viewers once they have been initialized
    bool _thumbnailPostInitFinished;
    bool _subvolumePostInitFinished;
  };
}

#endif
