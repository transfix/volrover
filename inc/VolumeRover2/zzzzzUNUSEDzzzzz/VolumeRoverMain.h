/*
  Copyright 2008 The University of Texas at Austin
  
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

/* $Id: VolumeRoverMain.h 4087 2011-05-06 14:54:32Z arand $ */

#ifndef __VOLUMEROVERMAIN_H__
#define __VOLUMEROVERMAIN_H__

#include <qglobal.h>

#if QT_VERSION < 0x040000
#include <qmainwindow.h>
#else
#include <QMainWindow>
#endif

#include <map>
#include <string>
#include <VolMagick/VolMagick.h>

#include <boost/shared_ptr.hpp>
#include <boost/any.hpp>

class QTabWidget;
class VolumeInterface;
class VolumeViewerPage;

#if QT_VERSION < 0x040000
class QWidgetStack;
class QListView;
class QListViewItem;
#else
class QStackedWidget;
class QTreeWidget;
class QTreeWidgetItem;
class QWidget;
#endif

//VolumeGridRover is currently not supported in Qt4
#if QT_VERSION >= 0x040000
#undef VOLUMEGRIDROVER
#endif

#ifdef VOLUMEGRIDROVER
class VolumeGridRover;
#endif

//class UI forward declaration
#if QT_VERSION < 0x040000
class VolumeRoverMainBase;
#else
namespace Ui
{
  class VolumeRoverMain;
}
#endif

// ---------------
// VolumeRoverMain
// ---------------
// Purpose: 
//   Main window of VolumeRover.  Contains the list of objects
//   For rendering/manipulation.  Use VolumeRoverMain::instance() to
//   retrieve the singleton instance of the main window from anywhere
//   in the application.
// ---- Change History ----
// ??/??/2008 -- Joe R. -- Initial implementation.
class VolumeRoverMain : public QMainWindow
{
  Q_OBJECT

 public:
  typedef boost::shared_ptr<VolumeRoverMain> VolumeRoverMainPtr;
  typedef std::map<std::string,boost::any> ObjectMap;

  static VolumeRoverMain& instance();

  virtual ~VolumeRoverMain();

  // ***** main API

  // Add objects to the object map.  If volrover cannot handle
  // a particular object type, it should complain...
  virtual void addObject(const std::string& name,
                         const boost::any& obj);

  // Remove the object from the map.
  virtual void removeObject(const std::string& name);

  // Removes all objects.
  virtual void clearObjects();

  // Retrieve an object from the map;
  virtual boost::any getObject(const std::string& name);

  virtual bool hasObject(const std::string& name) const;

  // Will attempt to load a file and stick it in the object map.
  // If you specify something for objname, it will use that for the object
  // name in the map.  Else the object name will be the filename.
  // Throws an exception upon failure.
  virtual bool addFile(const std::string& filename,
                       const std::string& objname = std::string());

  // 'Selects' the object.  A selected object will have it's interface
  // shown.  If it is a volume, it will be shown in the volume viewer.
  virtual boost::any selectObject(const std::string& name);

  // Returns the selected object if any.
  virtual boost::any getSelectedObject();

 protected:
  VolumeRoverMain( QWidget* parent = 0,
#if QT_VERSION < 0x040000
                   const char* name = 0, WFlags f = WType_TopLevel
#else
                   Qt::WFlags flags=0
#endif
                   );

  void addVolumeToStack(const VolMagick::VolumeFileInfo &vfi, bool force = false);

 public slots:
  
#ifdef VOLUMEGRIDROVER
  void toggleVolumeGridRoverSlot(bool show);
#endif

 protected slots:
  void newVolumeSlot();
  void openVolumeSlot();
  void closeVolumeSlot();

  void bilateralFilterSlot();

  void raiseSelected();
  void loadShownVolume(int);

#ifdef VOLUMEGRIDROVER
  void syncIsocontourValuesWithVolumeGridRover();
#endif

 protected:
  QTabWidget *_mainTabs;

#if QT_VERSION < 0x040000
  QWidgetStack *_volumeStack;
  QListView *_volumeStackList;
#else
  QStackedWidget *_volumeStack;
  QTreeWidget *_volumeStackList;
#endif

  VolumeViewerPage *_volumeViewerPage;

  std::map<void*,VolumeInterface*> _itemToInterface;
  std::map<VolumeInterface*,void*> _interfaceToItem;

#ifdef VOLUMEGRIDROVER
  VolumeGridRover *_volumeGridRover;
#endif

#if QT_VERSION < 0x040000
  VolumeRoverMainBase *_ui;
#else
  Ui::VolumeRoverMain *_ui;
  //Placeholder frame for the main window central widget
  QWidget *_centralWidget;
#endif

  static VolumeRoverMainPtr _instance;
  ObjectMap _objects; //stores all loaded objects
  ObjectMap::iterator _selectedObject;

 private:
  static VolumeRoverMainPtr instancePtr();
  
  //disallow copy construction
  VolumeRoverMain(const VolumeRoverMain&);
};

#endif
