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

#ifndef __CVCMAINWINDOW_H__
#define __CVCMAINWINDOW_H__

#include <QMainWindow>
#include <QEvent>
#include <QThread>

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/any.hpp>

#include <cvc_compat.h>
#include <cvc_compat.h>
#include <CVCEvent.h>

#include <xmlrpc/XmlRpc.h>

#include <string>
#include <map>
#include <typeinfo>

class QWidget;
class QTreeWidget;
class QTreeWidgetItem;
class QStackedWidget;
class QMenuBar;
class QTabWidget;
class QDialog;
#ifdef USING_VOLUMEGRIDROVER
class VolumeGridRover;
#endif

#ifdef USING_SECONDARYSTRUCTURES
class Skel;
class Histogram;
#endif

#if defined USING_MSLEVELSET && defined USING_MPSEGMENTATION
class MPLevelSetParams;
#endif 

namespace Ui
{
  class CVCMainWidget;
  class SegmentVirusMapDialog;
  class SecondaryStructureDialog;
  class generateRawVDialog;
  class generateMeshDialog;
  class HistogramDialog;
  class MSLevelSetDialog;
  class MPLevelSetDialog;
}

namespace CVC_NAMESPACE
{
    class windowbuf;
    class CVCMainWindow;

    class RemoteSegThread : public QThread 
      {
	public:
		RemoteSegThread(Ui::SegmentVirusMapDialog *dialog, CVCMainWindow *cvcmw, unsigned int stackSize = 0);
		~RemoteSegThread() {}
		virtual void run();
	private:
		CVCMainWindow *m_CVCMainWindow;
		XmlRpc::XmlRpcValue m_Params;
		int type;
		XmlRpc::XmlRpcClient m_XmlRpcClient;
      };
	
    class LocalSegThread : public QThread 
      {
        public:
		LocalSegThread(std::string filename, Ui::SegmentVirusMapDialog *dialog, CVCMainWindow *cvcmw, unsigned int stackSize = 0);
		~LocalSegThread() {}

		virtual void run();
        private:
		CVCMainWindow *m_CVCMainWindow;
		XmlRpc::XmlRpcValue m_Params;
		int type;
      };

    //12/16/2011 -- transfix -- added savePropertyMapSlot()
  class CVCMainWindow : public QMainWindow
  {
    Q_OBJECT
  public:
    typedef boost::shared_ptr<CVCMainWindow> CVCMainWindowPtr;
    typedef CVCEvent CVCMainWindowEvent;

    static CVCMainWindow& instance();
    static void terminate();

    template<class T>
      void insertDataWidget(QWidget* w)
      {
        boost::mutex::scoped_lock lock(_dataWidgetMapMutex);
        _dataWidgetMap[typeid(T).name()] = w;
      }

    virtual ~CVCMainWindow();

    //Use this to get this window's menu bar.  On Windows/X11 it
    //just calls QMainWindow::menuBar() but on MacOSX it returns the
    //parentless default menu bar below.
    QMenuBar *menu();

    //Use this to add your own tabs
    QTabWidget *tabWidget();

    //Get VolumeGridRover pointer to copy to Viewers
    #ifdef USING_VOLUMEGRIDROVER
    VolumeGridRover *volumeGridRover();
    #endif

    windowbuf* createWindowStreamBuffer();

    RemoteSegThread *m_RemoteSegThread;
    LocalSegThread *m_LocalSegThread;

  protected slots:
    void propertyMapItemChangedSlot(QTreeWidgetItem *item, int column);
    void addPropertySlot();
    void deleteSelectedPropertiesSlot();
    void dataMapItemClickedSlot(QTreeWidgetItem *item, int);


    void unimplementedSlot();

    // file menu
    void openFileSlot();
    void saveFileSlot();
    void closeDataSlot();
    void savePropertyMapSlot();
    void loadPropertyMapSlot();
    void saveImageSlot();
    void mainOptionsSlot();

    // view menu
    void selectCurrentVolumeSlot();
    void setSliceRenderingSlot();
#ifdef USING_VOLUMEGRIDROVER
    void volumeGridRoverSlot();
#endif
    void setBackgroundColorSlot();
    void setThumbnailBackgroundColorSlot();
    void geometryViewOptionsSlot();
    void toggleWireCubeSlot();
    void toggleClipBBoxSlot();

    // tools menu
    void anisotropicDiffusionSlot();
    void bilateralFilterSlot();
    void contrastEnhancementSlot();
    void gdtvFilterSlot();
    void virusSegmentationSlot();
    void secondaryStructureSlot();


#ifdef USING_MMHLS
	void generateRawVSlot();
	void generateMeshSlot();
	void openManifestFileSlot();
#endif

#ifdef USING_MSLEVELSET
	void MSLevelSetSlot();
#ifdef USING_HOSEGMENTATION
	void HOSegmentationSlot();
#ifdef USING_MPSEGMENTATION
	void MPSegmentationSlot();
	void on_UserSegReadButton_clickedSlot();
#endif
#endif
#endif

    void curationSlot();
    void hlsSurfaceSlot();
    void pocketTunnelSlot();
    void skeletonizationSlot();
    void tightCoconeSlot();
	void superSecondaryStructuresSlot();
	void contourTilerSlot();
    void LBIE_Slot();
    void LBIE_qualityImprovement_Slot();
#ifdef USING_RECONSTRUCTION
    void reconstructionSlot();
#endif    
    void raytraceSlot();
    void multiTileServerSlot();
 
    // help menu
    void helpSlot();
    void aboutVolRoverSlot();

#ifdef USING_SECONDARYSTRUCTURES
    void uploadgeometrySlot();
    void showHistogramDialogSlot();
#endif

  protected:
    CVCMainWindow(QWidget *parent = nullptr, Qt::WindowFlags flags=Qt::WindowFlags());

    void customEvent(QEvent *event);
    void initializePropertyMapWidget();
    void initializeDataMapWidget();
    void initializeThreadMapWidget();
    void initializeVariables();
#ifdef USING_SECONDARYSTRUCTURES
    Skel* m_Skeleton;
    Ui::HistogramDialog *SecondaryStructureHistogramUi;
    Histogram *m_AlphaHistogram;
    Histogram *m_BetaHistogram;
    QDialog *histogram_dialog;
    bool hasSecondaryStructure;
#endif
    void setupGridRover();

    //boost::signals2 slots for monitoring global data management changes
    void propertyMapChanged(const std::string&);
    void dataMapChanged(const std::string&);
    void threadMapChanged(const std::string&);

    boost::scoped_ptr<Ui::CVCMainWidget> _ui;

    //Placeholder frame for the main window central widget
    QWidget *_centralWidget;

    //Data map tab widgets
    QTreeWidget *_dataMap;
    QStackedWidget *_dataWidgetStack;

    //Mapping of data type name to a widget that can manipulate that
    //kind of data.
    std::map<std::string,QWidget*> _dataWidgetMap;
    boost::mutex _dataWidgetMapMutex;

    //Use a parentless menu bar on mac osx so all windows share it
#ifdef Q_WS_MAC
    QMenuBar *_defaultMenuBar;
#endif

    // keep track of the checkable menu items here...
    QAction * _viewWireCubeMenuAction;
    QAction * _clipBBoxMenuAction;

    static CVCMainWindowPtr instancePtr();
    static CVCMainWindowPtr _instance;
    static boost::mutex     _instanceMutex;

#ifdef USING_VOLUMEGRIDROVER
    VolumeGridRover *m_VolumeGridRover;
#endif

   Ui::SecondaryStructureDialog *SecondaryStructureUi;

private:

#ifdef USING_MMHLS
   Ui::generateRawVDialog *generateRawVUi;
   Ui::generateMeshDialog *generateMeshUi;
   QString m_manifestFile;
#endif

#ifdef USING_MSLEVELSET
	Ui::MSLevelSetDialog *MSLevelSetDialogUi;
	bool m_DoInitCUDA;
#ifdef  USING_MPSEGMENTATION
	Ui::MPLevelSetDialog *MPLevelSetDialogUi;
    MPLevelSetParams *MPLSParams;
#endif
#endif

  private:
    //disallow copy construction
    CVCMainWindow(const CVCMainWindow&);
  };
}

#endif
