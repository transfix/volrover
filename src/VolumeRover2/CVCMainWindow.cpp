/*
  Copyright 2008-2011 The University of Texas at Austin

        Authors: Joe Rivera <transfix@ices.utexas.edu>
                 Alex Rand <arand@ices.utexas.edu>
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

#include <log4cplus_compat.h>

#include <VolumeRover2/CVCMainWindow.h>
#include <VolumeRover2/DataWidget.h>
#include <VolumeRover2/SaveImageInfo.h>
#include <VolumeRover2/windowbuf.h>
#include <cvc_compat.h>

#include <ColorTable2/ColorTable.h>

#include <boost/current_function.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <QFrame>
#include <QList>
#include <QString>
#include <QStringList>
#include <QTreeWidget>
#include <QMessageBox>
#include <QMenu>
#include <QMenuBar>
#include <QFileDialog>
#include <QTabWidget>
#include <QSettings>
#include <QColorDialog>
#include <QGridLayout>
#include <QStackedWidget>
#include <QSplitter>

#include "ui_CVCMainWidget.h"
#include "ui_AddPropertyDialog.h"
#include "ui_SaveFileDialog.h"
#include "ui_CloseDataDialog.h"
#include "ui_MainOptions.h"
#include "ui_SaveImageDialog.h"
#include "ui_ViewGeometryOptionsDialog.h"
#include "ui_SliceRenderingDialog.h"
#include "ui_SelectCurrentVolumeDialog.h"
#include "ui_UnknownData.h"
//#include "ui_AnisotropicDiffusionDialog.h"
//#include "ui_BilateralFilterDialog.h"
//#include "ui_ContrastEnhancementDialog.h"
//#include "ui_GDTVFilterDialog.h"
#include "ui_SegmentVirusMapDialog.h"
#include "ui_SecondaryStructureDialog.h"

#include "VolumeRover2/AnisotropicDiffusionDialog.h"
#include "VolumeRover2/BilateralFilterDialog.h"
#include "VolumeRover2/ContourTilerDialog.h"
#include "VolumeRover2/ContrastEnhancementDialog.h"
#include "VolumeRover2/GDTVFilterDialog.h"
#include "VolumeRover2/MultiTileServerDialog.h"

//#include <Filters/OOCBilateralFilter.h>

#include <VolMagick/VolMagick.h>
#include <cvcraw_geometry/cvcgeom.h>
#include <cvcraw_geometry/io.h>
#include <cvcraw_geometry/cvcraw_geometry.h>
#include <cvcraw_geometry/Geometry.h>
#ifdef USING_TILING
#include <cvcraw_geometry/contours.h>
#endif

#include <xmlrpc/XmlRpc.h>

#ifdef USING_VOLUMEGRIDROVER
#include <VolumeGridRover/VolumeGridRover.h>
#endif

#ifdef USING_SEGMENTATION
#include <Segmentation/SegCapsid/segcapsid.h>
#include <Segmentation/SegMonomer/segmonomer.h>
#include <Segmentation/SegSubunit/segsubunit.h>
#include <Segmentation/SecStruct/secstruct.h>
#endif

#ifdef USING_TIGHT_COCONE
#include <VolumeRover2/TightCoconeDialog.h>
//#include "ui_TightCoconeDialog.h"
//#include <TightCocone/tight_cocone.h>
//#include <VolumeRover/TightCoconeDialog.h>
#endif

#ifdef USING_SUPERSECONDARY_STRUCTURES
#include <VolumeRover2/SuperSecondaryStructuresDialog.h>
#endif

#ifdef USING_CURATION
#include <VolumeRover2/CurationDialog.h>
#endif

//LBIE stuff
#include "ui_LBIE_dialog.h"
//#include "ui_LBIE_qualityImprovement.h"
#include <VolumeRover2/LBIEQualityImprovementDialog.h>
#include <LBIE/LBIE_Mesher.h>
//#include <LBIE/quality_improve.h>
#include <LBIE/octree.h>


#ifdef USING_SWEETMESH
#include <SweetMesh/triangle.h>
#include <SweetMesh/tetrahedra.h>
#include <SweetMesh/hexmesh.h>
#include <SweetMesh/volRoverDisplay.h>
#include <SweetMesh/vertex.h>
#include <SweetMesh/meshIO.h>
#endif

#ifdef USING_SWEETLBIE
#include "ui_LBIE_dialog.h"
#include <SweetLBIE/octree.h>
#endif

#ifdef USING_HLEVELSET
#include "ui_HLSSurfaceDialog.h"
#include <HLevelSet/HLevelSet.h>
#endif

#ifdef USING_POCKET_TUNNEL
#include <VolumeRover2/PocketTunnelDialog.h>
#endif

#ifdef USING_SKELETONIZATION
#include <VolumeRover2/SkeletonizationDialog.h>
//#include "ui_SkeletonizationDialog.h"
//#include <Skeletonization/Skeletonization.h>
//#include <VolumeRover/SkeletonizationDialog.h>
#endif

#ifdef USING_SECONDARYSTRUCTURES
#include "ui_HistogramDialog.h"
#include <SecondaryStructures/skel.h>
#include <Histogram/histogram.h>
using namespace SecondaryStructures;
#endif

#ifdef USING_RECONSTRUCTION
#include "VolumeRover2/ReconstructionDialog.h"
#endif


#ifdef USING_MMHLS
#include <MMHLS/generateRawV.h>
#include <MMHLS/generateMesh.h>
#include "ui_generateRawvDialog.h"
#include "ui_generateMeshDialog.h"
#endif


#ifdef USING_MSLEVELSET
#include "ui_mslevelsetDialog.h"
#include <MSLevelSet/levelset3D.h>
#ifdef USING_HOSEGMENTATION
#include <HigherOrderSegmentation/higherorderseg.h>
#ifdef USING_MPSEGMENTATION
#include "ui_mplevelsetDialog.h"
#include <MultiphaseSegmentation/multiphaseseg.h>
#endif
#endif
#endif

#ifdef USING_VOLUMEGRIDROVER
#include <VolumeGridRover/VolumeGridRover.h>
#endif

#include <VolMagick/VolumeFile_IO.h>

#include <list>
#include <sstream>

#ifndef VOLUMEROVER_VERSION_STRING
#define VOLUMEROVER_VERSION_STRING "2.0.0"
#endif

// arand: implemented, 5-2-2011
//   should we move this elsewhere?
class OpenFileThread
{
public:
  OpenFileThread(QWidget * mainWindow,QStringList filenames)
    : _mainWindow(mainWindow), _filenames(filenames) {}
  
  void operator()()
  {
    // this is a counter for saving full path of files to test muti-tile server
    static int readFileCnt = 0;

    CVC::ThreadFeedback feedback(BOOST_CURRENT_FUNCTION);
     
    QStringList failedToLoad;
    size_t fileIdx = 0;
    for(QStringList::iterator it = _filenames.begin();
        it != _filenames.end();
	++it)
      {
        cvcapp.threadProgress(double(fileIdx++)/_filenames.size());

	if((*it).isEmpty()) continue;
	
	if(!cvcapp.readData((*it).toStdString()))
	  failedToLoad.push_back(*it);
        // this is code to save&get full path of volume files for sync with multiTileServer
        else {
          char str[1024];
          sprintf( str, "file_%d_fullPath", readFileCnt++);
          cvcapp.properties( str, (*it).toStdString());
          sprintf( str, "%d", readFileCnt );
          cvcapp.properties("number_of_file_read", str);
        }
      }

    cvcapp.threadProgress(1.0);
    if(!failedToLoad.isEmpty())
      QCoreApplication::postEvent(_mainWindow,
	    new CVC_NAMESPACE::CVCMainWindow::CVCMainWindowEvent("openFileFailed")
          );
 
#ifdef USING_VOLUMEGRIDROVER
    QCoreApplication::postEvent(_mainWindow,
	new CVC_NAMESPACE::CVCMainWindow::CVCMainWindowEvent("setupGridRover")
    );
#endif       

  }
  
private:
  QWidget * _mainWindow;
  QStringList _filenames;
};

class OpenImageThread
{
public:
  OpenImageThread(QWidget * mainWindow,QStringList filenames)
    : _mainWindow(mainWindow), _filenames(filenames) {}
  
  void operator()()
  {
    // this is a counter for saving full path of files to test muti-tile server
    static int readFileCnt = 0;

    CVC::ThreadFeedback feedback(BOOST_CURRENT_FUNCTION);
     
    QStringList failedToLoad;
    size_t fileIdx = 0;
    for(QStringList::iterator it = _filenames.begin();
        it != _filenames.end();
	++it)
      {
        cvcapp.threadProgress(double(fileIdx++)/_filenames.size());

	if((*it).isEmpty()) continue;
	
	if(!cvcapp.readData((*it).toStdString()))
	  failedToLoad.push_back(*it);
        // this is code to save&get full path of volume files for sync with multiTileServer
        else {
          char str[1024];
          sprintf( str, "file_%d_fullPath", readFileCnt++);
          cvcapp.properties( str, (*it).toStdString());
          sprintf( str, "%d", readFileCnt );
          cvcapp.properties("number_of_file_read", str);
        }
      }

    cvcapp.threadProgress(1.0);
    if(!failedToLoad.isEmpty())
      QCoreApplication::postEvent(_mainWindow,
	    new CVC_NAMESPACE::CVCMainWindow::CVCMainWindowEvent("openFileFailed")
          );
 
#ifdef USING_VOLUMEGRIDROVER
    QCoreApplication::postEvent(_mainWindow,
	new CVC_NAMESPACE::CVCMainWindow::CVCMainWindowEvent("setupGridRover")
    );
#endif       

  }
  
private:
  QWidget * _mainWindow;
  QStringList _filenames;
};

namespace CVC_NAMESPACE
{
  /* Need these to signal the GUI thread the segmentation result so it can popup a message to the user */
  class SegmentationFailedEvent : public QEvent
  {
    public:
      SegmentationFailedEvent(const QString &m) : QEvent(static_cast<QEvent::Type>(QEvent::User+100)), msg(m) {}
      QString message() const { return msg; }
    private:
      QString msg;
  };

  class SegmentationFinishedEvent : public QEvent
  {
    public:
      SegmentationFinishedEvent(const QString &m) : QEvent(static_cast<QEvent::Type>(QEvent::User+101)), msg(m) {}
      QString message() const { return msg; }
    private:
      QString msg;
  };


  CVCMainWindow::CVCMainWindowPtr CVCMainWindow::_instance;
  boost::mutex CVCMainWindow::_instanceMutex;
  CVCMainWindow::CVCMainWindowPtr CVCMainWindow::instancePtr()
  {
    boost::mutex::scoped_lock lock(_instanceMutex);
    if(!_instance)
      {
        _instance.reset(new CVCMainWindow);

        //Set up slot connections in this factory function because
        //objects cannot be set up to track themselves in their
        //constructor.
        using namespace boost::placeholders;
        cvcapp.propertiesChanged.connect(
          MapChangeSignal::slot_type(
            &CVCMainWindow::propertyMapChanged, _instance.get(), _1
          ).track(_instance)
        );
        
        cvcapp.dataChanged.connect(
          MapChangeSignal::slot_type(
            &CVCMainWindow::dataMapChanged, _instance.get(), _1
          ).track(_instance)
        );
        
        cvcapp.threadsChanged.connect(
          MapChangeSignal::slot_type(
            &CVCMainWindow::threadMapChanged, _instance.get(), _1
          ).track(_instance)
        );
      }
        
    return _instance;
  }

  CVCMainWindow& CVCMainWindow::instance()
  {
    return *instancePtr();
  }

  void CVCMainWindow::terminate()
  {
    boost::mutex::scoped_lock lock(_instanceMutex);
    _instance.reset();
  }

  CVCMainWindow::CVCMainWindow(QWidget *parent, Qt::WindowFlags flags)
    : QMainWindow(parent,flags),
      _centralWidget(NULL),
      _dataMap(NULL),
      _dataWidgetStack(NULL)
  {
    _centralWidget = new QWidget;
    setCentralWidget(_centralWidget);
    _ui.reset(new Ui::CVCMainWidget);
    _ui->setupUi(_centralWidget);

#ifdef Q_WS_MAC
    _defaultMenuBar = new QMenuBar(0);
#endif

    QMenu *fileMenu = menu()->addMenu(tr("&File"));
    fileMenu->addAction(tr("&Open"), this, SLOT(openFileSlot()));
    fileMenu->addAction(tr("&Save"), this, SLOT(saveFileSlot()));
    fileMenu->addAction(tr("&Close"), this, SLOT(closeDataSlot()));
    fileMenu->addSeparator();
    fileMenu->addAction(tr("Save &Property Map"), this, SLOT(savePropertyMapSlot()));
    fileMenu->addAction(tr("Load &Property Map"), this, SLOT(loadPropertyMapSlot()));
    fileMenu->addSeparator();
    fileMenu->addAction(tr("S&ave Image"), this, SLOT(saveImageSlot()));
    fileMenu->addSeparator();
    fileMenu->addAction(tr("O&ptions"), this, SLOT(mainOptionsSlot()));
    fileMenu->addSeparator();
    fileMenu->addAction(tr("&Quit"), this, SLOT(close()));

    QMenu *viewMenu = menu()->addMenu(tr("&View"));
    viewMenu->addAction(tr("&Select Current Volume"), this, SLOT(selectCurrentVolumeSlot()));
    viewMenu->addSeparator();
    // arand: moved grid rover to the view menu
	viewMenu->addAction(tr("Set Slice &Rendering"),this, SLOT(setSliceRenderingSlot())); 
#ifdef USING_VOLUMEGRIDROVER
    viewMenu->addAction(tr("&Show Volume Grid Rover"), this, SLOT(volumeGridRoverSlot()));
    viewMenu->addSeparator();
    m_VolumeGridRover = new VolumeGridRover();
#endif

    viewMenu->addAction(tr("&Set Background Color"), this, SLOT(setBackgroundColorSlot()));
    viewMenu->addAction(tr("&Set Thumbnail Background Color"), this, SLOT(setThumbnailBackgroundColorSlot()));
    viewMenu->addSeparator();
    viewMenu->addAction(tr("&Set Geometry View Options"), this, SLOT(geometryViewOptionsSlot()));
    _viewWireCubeMenuAction = viewMenu->addAction(tr("&View Wire Cube"), this, SLOT(toggleWireCubeSlot()));
    _clipBBoxMenuAction = viewMenu->addAction(tr("&Clip Geometry to Bounding Box"), this, SLOT(toggleClipBBoxSlot()));

    _viewWireCubeMenuAction->setCheckable(true);
    _viewWireCubeMenuAction->setChecked(true);

    _clipBBoxMenuAction->setCheckable(true);
    _clipBBoxMenuAction->setChecked(true);



    QMenu *toolsMenu = menu()->addMenu(tr("&Tools"));
    // volumetric filtering functions
    toolsMenu->addAction(tr("&Anisotropic Diffusion"), this, SLOT(anisotropicDiffusionSlot()));
    toolsMenu->addAction(tr("&Bilateral Filtering"), this, SLOT(bilateralFilterSlot()));
    toolsMenu->addAction(tr("&Contrast Enhancement"), this, SLOT(contrastEnhancementSlot()));
    toolsMenu->addAction(tr("&GDTV Filtering"), this, SLOT(gdtvFilterSlot()));

    toolsMenu->addSeparator();
    toolsMenu->addAction(tr("&Segment Virus Map"), this, SLOT(virusSegmentationSlot()));
#ifdef USING_SECONDARYSTRUCTURES
    toolsMenu->addAction(tr("&Secondary Structure Elucidation"), this, SLOT(secondaryStructureSlot()));
    toolsMenu->addAction(tr("  > &Show Histogram Dialog"), this, SLOT(showHistogramDialogSlot()));
#endif
    toolsMenu->addSeparator();
    toolsMenu->addAction(tr("&Surface Curation"), this, SLOT(curationSlot()));
    toolsMenu->addAction(tr("&HLS Surface"), this, SLOT(hlsSurfaceSlot()));
    toolsMenu->addAction(tr("&Pocket/Tunnel Detection"), this, SLOT(pocketTunnelSlot()));
    toolsMenu->addAction(tr("&Skeletonization"), this, SLOT(skeletonizationSlot()));
    toolsMenu->addAction(tr("&Tight Cocone"), this, SLOT(tightCoconeSlot()));
    toolsMenu->addSeparator();
    toolsMenu->addAction(tr("&Contour Tiler"), this, SLOT(contourTilerSlot()));
    toolsMenu->addAction(tr("&LBIE Mesh Generation"), this, SLOT(LBIE_Slot()));
    toolsMenu->addAction(tr("&LBIE Quality Improvement"), this, SLOT(LBIE_qualityImprovement_Slot()));
#ifdef USING_SUPERSECONDARY_STRUCTURES
	toolsMenu->addAction(tr("&SuperSecondary Structures"), this, SLOT(superSecondaryStructuresSlot()));
#endif

#ifdef USING_MMHLS
	QMenu* mmhlsMenu = toolsMenu->addMenu(tr("&MMHLS"));
	mmhlsMenu->addAction(tr("&GenerateRawV"), this, SLOT(generateRawVSlot()));
	mmhlsMenu->addAction(tr("Generate&Mesh"), this, SLOT(generateMeshSlot()));
#endif
#ifdef USING_MSLEVELSET
	toolsMenu->addAction(tr("MS&LevelSet Segmentation"), this, SLOT(MSLevelSetSlot()));
#ifdef USING_HOSEGMENTATION
	toolsMenu->addAction(tr("H&OLevelSet Segmentation"), this, SLOT(HOSegmentationSlot()));
#ifdef USING_MPSEGMENTATION
	toolsMenu->addAction(tr("Multi&Phase Segmentation"), this, SLOT(MPSegmentationSlot()));
#endif
#endif
#endif

#ifdef USING_RECONSTRUCTION
    toolsMenu->addSeparator();
    toolsMenu->addAction(tr("&ET Reconstruction"), this, SLOT(reconstructionSlot()));
#endif
    toolsMenu->addSeparator();
    toolsMenu->addAction(tr("&Ray Traced Image"), this, SLOT(raytraceSlot()));

    toolsMenu->addAction(tr("&Sync Multi-Tile Server"), this, SLOT(multiTileServerSlot()));

    QMenu *helpMenu = menu()->addMenu(tr("&Help"));
    helpMenu->addAction(tr("&Help"), this, SLOT(helpSlot()));
    helpMenu->addAction(tr("&About"), this, SLOT(aboutVolRoverSlot()));

    //create the data map tab and associated widgets
    QWidget *dataTab = new QWidget;
    QGridLayout *dataTabLayout = new QGridLayout(dataTab);
    QSplitter *dataTabSplitter = new QSplitter(dataTab);
    dataTabLayout->addWidget(dataTabSplitter);
    _dataMap = new QTreeWidget;
    QTreeWidgetItem *dataitem = _dataMap->headerItem();
    dataitem->setText(1, QApplication::translate("CVCMainWidget", "Type"));
    dataitem->setText(0, QApplication::translate("CVCMainWidget", "Object Key"));
    dataTabSplitter->addWidget(_dataMap);
    _dataWidgetStack = new QStackedWidget;
    dataTabSplitter->addWidget(_dataWidgetStack);
    _ui->_tabs->insertTab(1,dataTab,"Data");
    _ui->_tabs->
      setTabText(_ui->_tabs->indexOf(dataTab),
                 QApplication::translate("CVCMainWidget", 
                                         "Data"));

    connect(_ui->_propertyMap,
            SIGNAL(itemChanged(QTreeWidgetItem*,int)),
            SLOT(propertyMapItemChangedSlot(QTreeWidgetItem*,int)));
    connect(_dataMap,
            SIGNAL(itemClicked(QTreeWidgetItem*,int)),
            SLOT(dataMapItemClickedSlot(QTreeWidgetItem*,int)));
    connect(_ui->_addProperty,
            SIGNAL(clicked()),
            SLOT(addPropertySlot()));
    connect(_ui->_deleteSelectedProperties,
            SIGNAL(clicked()),
            SLOT(deleteSelectedPropertiesSlot()));
    
    initializePropertyMapWidget();
    initializeDataMapWidget();
    initializeThreadMapWidget();
    initializeVariables();

    setWindowTitle("VolumeRover - " VOLUMEROVER_VERSION_STRING);
  }

  CVCMainWindow::~CVCMainWindow()
  {
#ifdef Q_WS_MAC
    delete _defaultMenuBar;
#endif
    if(m_RemoteSegThread) delete m_RemoteSegThread;
    if(m_LocalSegThread) delete m_LocalSegThread;

#ifdef USING_MMHLS
	if(generateRawVUi) delete generateRawVUi;
	if(generateMeshUi) delete generateMeshUi;
#endif


#ifdef USING_MSLEVELSET
	if(MSLevelSetDialogUi) delete MSLevelSetDialogUi;
#ifdef USING_MPSEGMENTATION
	if(MPLevelSetDialogUi) delete MPLevelSetDialogUi;
	if(MPLSParams) delete MPLSParams;
#endif
#endif

  }

  QMenuBar *CVCMainWindow::menu()
  {
#ifdef Q_WS_MAC
    return _defaultMenuBar;
#else
    return QMainWindow::menuBar();
#endif
  }

  void CVCMainWindow::initializeVariables()
  {
     m_RemoteSegThread = NULL;
     m_LocalSegThread = NULL;

#ifdef USING_SECONDARYSTRUCTURES
     m_Skeleton = NULL;
     SecondaryStructureHistogramUi = NULL;
     histogram_dialog = NULL;
     hasSecondaryStructure = false;
#endif

#ifdef USING_MMHLS
	generateRawVUi = NULL;
	generateMeshUi = NULL;
#endif

#ifdef USING_MSLEVELSET
	MSLevelSetDialogUi  = NULL;
    m_DoInitCUDA = true;
#ifdef USING_MPSEGMENTATION
	MPLevelSetDialogUi = NULL;
	MPLSParams = NULL;
#endif
#endif
  }

  QTabWidget *CVCMainWindow::tabWidget()
  {
    return _ui->_tabs;
  }

#ifdef USING_VOLUMEGRIDROVER
  VolumeGridRover *CVCMainWindow::volumeGridRover()
  {
    return m_VolumeGridRover;
  }
#endif

  // windowbuf* CVCMainWindow::createWindowStreamBuffer()
  // {
  //   return new CVC::windowbuf(_ui->_logTextEdit);
  // }

  void CVCMainWindow::propertyMapItemChangedSlot(QTreeWidgetItem *item, int column)
  {
    switch(column)
      {
      case 0: //key column
        {
          QMessageBox::warning(this,"Warning",
                               "Cannot change key directly. Create a new property instead.");
          //Restore the widget as it was before edit.  Can't do it here directly because we'll segfault,
          //so do it after this slot returns.
          QCoreApplication::postEvent(this,
             new CVCMainWindowEvent("initializePropertyMapWidget")
          );
          break;
        }
      case 1: //value column
        {
          cvcapp.properties(item->text(0).toStdString(),
                            item->text(1).toStdString());
          break;
        }
      }
  }

  void CVCMainWindow::addPropertySlot()
  {
    Ui::AddPropertyDialog ui;
    QDialog dialog;
    ui.setupUi(&dialog);
    if(dialog.exec() == QDialog::Accepted)
      {
        std::string key = ui._keyInput->text().toStdString();
        std::string val = ui._valueInput->text().toStdString();
        
        if(!cvcapp.properties(key).empty())
          {
            if(QMessageBox::Ok != 
               QMessageBox::warning(this,"Warning",
                                    QString("Key %1 exists, overwrite?").arg(QString::fromStdString(key)),
                                    QMessageBox::Ok | QMessageBox::Cancel,
                                    QMessageBox::Cancel))
              return;
          }

        if(val.empty())
          {
            QMessageBox::critical(this,"Error","No value set");
            return;
          }

        cvcapp.properties(key,val);
      }
  }

  void CVCMainWindow::deleteSelectedPropertiesSlot()
  {
    QList<QTreeWidgetItem*> items = _ui->_propertyMap->selectedItems();
    PropertyMap map = cvcapp.properties();
    for (const auto& item : items)
      map[item->text(0).toStdString()] = "";
    cvcapp.properties(map);
  }

  void CVCMainWindow::dataMapItemClickedSlot(QTreeWidgetItem *item, int)
  {
    boost::mutex::scoped_lock lock(_dataWidgetMapMutex);

    //Get the raw type name from the datamap for this item and use it to lookup
    //a registered widget for that type.
    std::string typestr = cvcapp.data(item->text(0).toStdString()).type().name();

    //If a widget isn't found that can handle the data type, use the 'UnknownData' widget
    if(_dataWidgetMap.find(typestr)==_dataWidgetMap.end())
      {
        QWidget* newwidget = new QWidget;
        _dataWidgetMap[typestr] = newwidget;
        Ui::UnknownData ud;
        ud.setupUi(newwidget);
      }

    //If the data type's widget isn't in the stack, add it
    if(_dataWidgetStack->indexOf(_dataWidgetMap[typestr]) == -1)
      _dataWidgetStack->insertWidget(0,_dataWidgetMap[typestr]);

    //now raise the widget so it is visible
    _dataWidgetStack->setCurrentWidget(_dataWidgetMap[typestr]);

    //finally, initialize the widget with the data in the map
    DataWidget *dw = 
      dynamic_cast<DataWidget*>(_dataWidgetMap[typestr]);
    if(dw) dw->initialize(item->text(0).toStdString());
  }

  void CVCMainWindow::openFileSlot()
  {
    static log4cplus::Logger logger = FUNCTION_LOGGER;

    std::vector<std::string> extensions = VolMagick::VolumeFile_IO::getExtensions();
    QString list("volume (");
    for (const auto& ext : extensions) {
        list += QString("*%1 ").arg(QString::fromStdString(ext));
      }
    list += ");;";

    QSettings settings;
    QString selectedFilter = settings.value("openFile/filter", "").toString();
    QString dir = settings.value("openFile/dir", "").toString();
    
    QStringList filenames = QFileDialog::getOpenFileNames(this,
                                                          "Open File",
                                                          dir,
                                                          list +
                                                          "geometry (*.raw *.rawn *.rawc *.rawnc *.off);;"
                                                          "Reconstruct series (*.ser);;"
							  "Images (*.png *.jpg *.jpeg *.bmp *.xbm *.xpm *.pnm *.mng *.gif);;"
						          "2D MRC (*.mrc *.map);;"
                                                          "All Files (*)",
                                                          &selectedFilter);

    settings.setValue("openFile/filter", selectedFilter);
    if (!filenames.empty()) {
      QFileInfo fileInfo(*filenames.begin());
      dir = fileInfo.dir().absolutePath();
      settings.setValue("openFile/dir", dir);
      LOG4CPLUS_TRACE(logger, "dir = " << dir.toStdString());
    }
   
    cvcapp.startThread("open_file_thread", OpenFileThread(this,filenames));
   
  }

  // arand, 4-21-2011: implemented
  // TODO: currently this only handles geometries... in the future handle volumes also
  void CVCMainWindow::saveFileSlot() {

    Ui::SaveFileDialog ui;
    QDialog dialog;
    ui.setupUi(&dialog);


    bool found = false;
    DataMap map = cvcapp.data();
    for (const auto& val : map) { 
      std::string myname(val.first);
      std::string mytype(val.second.type().name());
      // only deal with geometry files for now...
      if (cvcapp.isData<cvcraw_geometry::cvcgeom_t>(myname)) {
	ui.GeometryList->addItem(QString::fromStdString(myname));
	found = true;      
      }
      
      if (cvcapp.isData<VolMagick::Volume>(myname)) {
	ui.VolumeList->addItem(QString::fromStdString(myname));
	found = true;      
      }
      
      if (cvcapp.isData<LBIE::Mesher>(myname)) {
	ui.meshList->addItem(QString::fromStdString(myname));
	found = true;
      }

    }

    if (!found) {
      QMessageBox::information(this, tr("VolRover"),
			   tr("No geometries loaded."), QMessageBox::Ok);
      return;
    }

    if(dialog.exec() == QDialog::Accepted) {

      if (ui.tabWidget->currentIndex() == 0) {
	std::string geomSelected = ui.GeometryList->currentText().toStdString();
	
	cvcraw_geometry::cvcgeom_t geom = boost::any_cast<cvcraw_geometry::cvcgeom_t>(map[geomSelected]);
	
	string outputFile = ui.outputFileEdit->displayText().toStdString();
	if (outputFile.size() == 0) {
	  QMessageBox::information(this, tr("VolRover"),
				   tr("No output filename given."), QMessageBox::Ok);
	  return;
	}
	
	int type = ui.fileTypeComboBox->currentIndex();
	// 0: raw; 1: rawn; 2: rawc; 3: rawnc 
	
	bool problem = false;

	if (type == 3) {
	  // rawnc file type...
	  if (geom.const_normals().size() == 0) {
	    // no normals, revert to rawc
	    type = 2; 
	    problem = true;
	  } else if (geom.const_colors().size() == 0) {
	    // no colors, revert to rawn
	    type = 1;
	    problem = true;
	  }
	}
	if (type == 2 && geom.const_colors().size() == 0) {
	  // no colors, revert to raw
	  type = 0;
	  problem = true;
	}
	if (type == 1 && geom.const_normals().size() == 0) {
	  // no normals, revert to raw
	  type = 0;
	  problem = true;
	}

	string fileExtension;
	if (type == 0) {
	  fileExtension = ".raw";
	} else if (type == 1) {
	  fileExtension = ".rawn";
	} else if (type == 2) {
	  fileExtension = ".rawc";
	} else if (type == 3) {
	  fileExtension = ".rawnc";
	}

	if (problem) {
	  string message = "Requested file type does not match available data.  Despite the filename, your actual file will be of type " + fileExtension;
	  QMessageBox::information(this, tr("VolRover"),
				   tr(message.c_str()), QMessageBox::Ok);
	}
        
	// TODO: eliminate uneccessary conversion...
	cvcraw_geometry::write(cvcraw_geometry::geometry_t(geom), outputFile+fileExtension);
      }
      if (ui.tabWidget->currentIndex() == 1) {
	std::string volSelected = ui.VolumeList->currentText().toStdString();
	
	VolMagick::Volume vol = boost::any_cast<VolMagick::Volume>(map[volSelected]);

	string outputFile = ui.volOutputFileEdit->displayText().toStdString();
	
	std::string fileExtension = ui.volFileTypeComboBox->currentText().toStdString();
	createVolumeFile(cvcapp, vol,outputFile + "." + fileExtension);
      }
      if (ui.tabWidget->currentIndex() == 2) {
	std::string meshSelected = ui.meshList->currentText().toStdString();
	string outputFile = ui.meshOutputFileEdit->displayText().toStdString() + ".raw";
	LBIE::Mesher mesher = boost::any_cast<LBIE::Mesher>(map[meshSelected]);
	switch(mesher.mesh().mesh_type){
	  case LBIE::geoframe::SINGLE:{
	    mesher.mesh().saveTriangle(outputFile.c_str());
	  }break;
	  case LBIE::geoframe::TETRA:{
	    outputFile += "ts";
	    mesher.mesh().saveTetra(outputFile.c_str());
	  }break;
	  case LBIE::geoframe::QUAD:{
	    mesher.mesh().saveQuad(outputFile.c_str());
	  }break;
	  case LBIE::geoframe::HEXA:{
	    outputFile += "hs";
	    mesher.mesh().saveHexa(outputFile.c_str());
	  }break;
	  case LBIE::geoframe::DOUBLE:{
	    mesher.mesh().saveTriangle(outputFile.c_str());
	  }break;
	  case LBIE::geoframe::TETRA2:{
	    outputFile += "ts";
	    mesher.mesh().saveTetra(outputFile.c_str());
	  }break;
	  default: break;
	}
      }

    } 

    
  }

  // arand, 4-21-2011: implemented
  void CVCMainWindow::closeDataSlot() {

    Ui::CloseDataDialog ui;
    QDialog dialog;
    ui.setupUi(&dialog);


    bool found = false;
    DataMap map = cvcapp.data();
    for (const auto& val : map) { 
      std::string myname(val.first);

  if (cvcapp.isData<VolMagick::Volume>(myname) || 
      cvcapp.isData<VolMagick::VolumeFileInfo>(myname) || 
#ifdef USING_TILING
      cvcapp.isData<cvcraw_geometry::contours_t>(myname) ||
#endif
      cvcapp.isData<cvcraw_geometry::cvcgeom_t>(myname)) {
	if (myname.compare("thumbnail_isocontour") != 0 && 
	    myname.compare("zoomed_isocontour") != 0 &&
	    myname.compare("thumbnail_volume") != 0  && 
	    myname.compare("zoomed_volume") != 0 ) {
	  ui.DataList->addItem(QString::fromStdString(myname));
	  found = true;
	}
      }
    }

    if (!found) {
      QMessageBox::information(this, tr("VolRover"),
			   tr("No data loaded."), QMessageBox::Ok);
      return;
    }

    if(dialog.exec() == QDialog::Accepted) { 

      std::string dataSelected = ui.DataList->currentText().toStdString();
      std::string typeSelected = map[dataSelected].type().name();

      if (cvcapp.isData<cvcraw_geometry::cvcgeom_t>(dataSelected)
#ifdef USING_TILING
	  || cvcapp.isData<cvcraw_geometry::contours_t>(dataSelected)
#endif
	  ) {
	cout << "Close Geom." << endl;	

	// arand: 5-3-2011, simple functions below added eliminating manual code
	cvcapp.listPropertyRemove("thumbnail.geometries", dataSelected);
	cvcapp.listPropertyRemove("zoomed.geometries", dataSelected);

	/*
	vector<string> thumbGeoNames = cvcapp.listProperty("thumbnail.geometries",true);
	string thumbGeoNames1;
	for (int i=0; i<thumbGeoNames.size(); i++) {
	  if (dataSelected.compare(thumbGeoNames[i]) != 0) {
	      if (thumbGeoNames1.size() > 0) {
		thumbGeoNames1 = thumbGeoNames1 + ",";
	      }
	      thumbGeoNames1 = thumbGeoNames1 + thumbGeoNames[i];
	  }
	}
	cvcapp.properties("thumbnail.geometries", thumbGeoNames1);

	vector<string> zoomedGeoNames = cvcapp.listProperty("zoomed.geometries",true);
	string zoomedGeoNames1;
	for (int i=0; i<zoomedGeoNames.size(); i++) {
	  if (dataSelected.compare(zoomedGeoNames[i]) != 0) {
	      if (zoomedGeoNames1.size() > 0) {
		zoomedGeoNames1 = zoomedGeoNames1 + ",";
	      }
	      zoomedGeoNames1 = zoomedGeoNames1 + zoomedGeoNames[i];
	  }
	}
	cvcapp.properties("zoomed.geometries", zoomedGeoNames1);
	*/

	// close the data file
	boost::any tmp;
	cvcapp.data(dataSelected,tmp);

      } else if (cvcapp.isData<VolMagick::VolumeFileInfo>(dataSelected)) {
	if (dataSelected.compare(cvcapp.properties("viewers.vfi")) == 0) {
	  string tmp = "none";
	  cvcapp.properties("viewers.vfi",tmp);
	}	
	
	// close the data file
	boost::any tmp;
	cvcapp.data(dataSelected,tmp);

     
      } else if (cvcapp.isData<VolMagick::Volume>(dataSelected)) {
	// close the data file, this _should_ free the memory
	boost::any tmp;
	cvcapp.data(dataSelected,tmp);

      }  else {
	QMessageBox::information(this, tr("VolRover"),
				 tr("Error: unknown data type."), QMessageBox::Ok);
      }
    }    
  }

  // 12-16-2011 -- transfix -- initial implementation
  void CVCMainWindow::savePropertyMapSlot()
  {
    QString filename = QFileDialog::getSaveFileName(this,
                                                    tr("Save Property Map"),
                                                    QString(),
                                                    "VolumeRover Property Map (*.vpm *.xml *.info);;"
                                                    "All Files (*)");
    if(filename.isEmpty()) return;
    try
      {
        cvcapp.writePropertyMap(filename.toStdString());
      }
    catch(std::exception& e)
      {
        QMessageBox::critical(this,tr("Error"),
                              tr("Could not save property map file ")+
                              QString::fromStdString(e.what()),
                              QMessageBox::Ok);
      }
  }

  // 12-16-2011 -- transfix -- initial implementation
  void CVCMainWindow::loadPropertyMapSlot()
  {
    QStringList filenames = QFileDialog::getOpenFileNames(this,
                                                          tr("Load Property Map"),
                                                          QString(),
                                                          "VolumeRover Property Map (*.vpm *.xml *.info);;"
                                                          "All Files (*)");
    if(filenames.isEmpty()) return;
    try
      {
        for (auto& str : filenames)
          cvcapp.readPropertyMap(str.toStdString());

	// load datafiles

	QStringList filenames;	
	if (cvcapp.hasProperty("file_0_fullPath")) {
	  QString fn((cvcapp.properties("file_0_fullPath")).c_str());
	  filenames.append(fn);
	}
	if (cvcapp.hasProperty("file_1_fullPath")) {
	  QString fn((cvcapp.properties("file_1_fullPath")).c_str());
	  filenames.append(fn);
	}
	if (cvcapp.hasProperty("file_2_fullPath")) {
	  QString fn((cvcapp.properties("file_2_fullPath")).c_str());
	  filenames.append(fn);
	}
	if (cvcapp.hasProperty("file_3_fullPath")) {
	  QString fn((cvcapp.properties("file_3_fullPath")).c_str());
	  filenames.append(fn);
	}
	if (cvcapp.hasProperty("file_4_fullPath")) {
	  QString fn((cvcapp.properties("file_4_fullPath")).c_str());
	  filenames.append(fn);
	}
	
	cvcapp.startThread("open_file_thread", OpenFileThread(this,filenames));

	// load vinay file
	if (cvcapp.hasProperty("transfer_function_fullPath")) {
	  // do nothing, the Viewers widget automatically detects this.
	}


      }
    catch(std::exception& e)
      {
        QMessageBox::critical(this,tr("Error"),
                              tr("Could not load property map file ")+
                              QString::fromStdString(e.what()),
                              QMessageBox::Ok);
      }
  }

  // arand, 6-8-2011: started
  void CVCMainWindow::saveImageSlot() {
    Ui::SaveImageDialog ui;
    QDialog dialog;
    ui.setupUi(&dialog);
    
    if(dialog.exec() == QDialog::Accepted) { 
      
      int panel = ui.PanelList->currentIndex();

      // arand: I am not totally sure if this is the "right" 
      //        way to pass info to the viewers...
      QCoreApplication::postEvent(tabWidget()->widget(0),
	    new CVCEvent("saveImage",panel));

    }

  }

  
  void CVCMainWindow::mainOptionsSlot() {
    Ui::MainOptionsDialog ui;
    QDialog dialog;
    ui.setupUi(&dialog);

    if (cvcapp.properties("viewers.colortable_interactive_updates").compare("true") == 0) {
      ui.UpdateComboBox->setCurrentIndex(1);
    }

    if (cvcapp.properties("thumbnail.shaded_rendering_enabled").compare("true") == 0) {
      ui.VolRendComboBox->setCurrentIndex(1);
    }

    if (cvcapp.properties("viewers.rendering_mode").compare("RGBA") == 0 ||
        cvcapp.properties("viewers.rendering_mode").compare("rgba") == 0) {
      ui.RendStyleComboBox->setCurrentIndex(1);
    }

    // also VolRendComboBox : Shaded or Unshaded
    //      {Zoomed,Thumb}RendStyleComboBox : Color Mapped or RGBA Combined

    if(dialog.exec() == QDialog::Accepted) { 

      string updateMode = ui.UpdateComboBox->currentText().toStdString();
      if (updateMode.compare("Interactive") == 0)
	cvcapp.properties("viewers.colortable_interactive_updates", "true");
      else 
	cvcapp.properties("viewers.colortable_interactive_updates", "false");

      string volRendMode = ui.VolRendComboBox->currentText().toStdString();
      if (volRendMode.compare("Shaded") == 0) {
	cvcapp.properties("thumbnail.shaded_rendering_enabled", "true");
	cvcapp.properties("zoomed.shaded_rendering_enabled", "true");
      }
      else {
	cvcapp.properties("thumbnail.shaded_rendering_enabled", "false");
	cvcapp.properties("zoomed.shaded_rendering_enabled", "false");
      }

      string rendMode = ui.RendStyleComboBox->currentText().toStdString();
      if (rendMode.compare("RGBA Combined") == 0)
	cvcapp.properties("viewers.rendering_mode", "rgba");
      else 
	cvcapp.properties("viewers.rendering_mode", "colormapped");
    }
  }


  // arand: functionality added, 5-2-2011
  //  TODO: fix the volume accessed by the grid rover!
  void CVCMainWindow::selectCurrentVolumeSlot() {
    Ui::SelectCurrentVolumeDialog ui;
    QDialog dialog;
    ui.setupUi(&dialog);

    std::vector<std::string> keys;    
    std::vector<std::string> vfiKeys = cvcapp.data<VolMagick::VolumeFileInfo>();
    std::vector<std::string> volKeys = cvcapp.data<VolMagick::Volume>();
    keys.insert(keys.end(),
		vfiKeys.begin(),vfiKeys.end());
    keys.insert(keys.end(),
		volKeys.begin(),volKeys.end());
    
    if(keys.empty())
      {
	QMessageBox::information(this, tr("VolumeRover"),
                                 tr("No volumes loaded."), QMessageBox::Ok);
	return;
      }
    
    for (const auto& key : keys) {
      if (key.compare("zoomed_volume") != 0  && 
	  key.compare("thumbnail_volume") != 0  )
      ui.VolumeList->addItem(QString::fromStdString(key));
    }

    ui.VolumeList->addItem(QString::fromStdString("none"));

    if(dialog.exec() == QDialog::Accepted) {  

      std::string newvolume = ui.VolumeList->currentText().toStdString();

      cvcapp.properties("viewers.vfi",newvolume);

    }    
  }
  
  
  void CVCMainWindow::setBackgroundColorSlot() {
    QColorDialog qd;
    // set current color to current background color

    QColor newColor = qd.getColor();
    // set the background color
    cvcapp.properties("zoomed.background_color",newColor.name().toStdString());

  }

  void CVCMainWindow::setThumbnailBackgroundColorSlot() {
    QColorDialog qd;
    // set current color to current background color

    QColor newColor = qd.getColor();
    // set the background color
    cvcapp.properties("thumbnail.background_color",newColor.name().toStdString());
  }

 void CVCMainWindow::setSliceRenderingSlot(){
//           QMessageBox::warning(this,"Warning",
  //                             "Not implemented right now.");
  	Ui::SliceRenderingDialog ui;
	QDialog dialog;
	ui.setupUi(&dialog);
	  if(dialog.exec() == QDialog::Accepted) {
QMessageBox::warning(this,"Warning", "here it is");

	  }

 }


  // arand, 4-21-2011: implmentation finished
  // arand, 4-25-2011: added isocontouring method control
  void CVCMainWindow::geometryViewOptionsSlot() {

    Ui::GeometryViewOptionsDialog ui;
    QDialog dialog;
    ui.setupUi(&dialog);

    if (cvcapp.properties("thumbnail.geometry_rendering_mode").compare("wireframe") == 0) {
      ui.ThumbnailViewCombo->setCurrentIndex(1);
    } else if (cvcapp.properties("thumbnail.geometry_rendering_mode").compare("filled_wireframe") == 0) {
      ui.ThumbnailViewCombo->setCurrentIndex(2);
    } else if (cvcapp.properties("thumbnail.geometry_rendering_mode").compare("flat") == 0) {
      ui.ThumbnailViewCombo->setCurrentIndex(3);
    } else if (cvcapp.properties("thumbnail.geometry_rendering_mode").compare("flat_wireframe") == 0) {
      ui.ThumbnailViewCombo->setCurrentIndex(4);
    }

    if (cvcapp.properties("zoomed.geometry_rendering_mode").compare("wireframe") == 0) {
      ui.ZoomedViewCombo->setCurrentIndex(1);
    } else if (cvcapp.properties("zoomed.geometry_rendering_mode").compare("filled_wireframe") == 0) {
      ui.ZoomedViewCombo->setCurrentIndex(2);
    } else if (cvcapp.properties("zoomed.geometry_rendering_mode").compare("flat") == 0) {
      ui.ZoomedViewCombo->setCurrentIndex(3);
    } else if (cvcapp.properties("zoomed.geometry_rendering_mode").compare("flat_wireframe") == 0) {
      ui.ZoomedViewCombo->setCurrentIndex(4);
    }
    
    // look for the isocontour in the geometry view list...
    ui.ThumbnailIsoCheckBox->setChecked(false);
    ui.ZoomedIsoCheckBox->setChecked(false);

    {      
      std::vector<std::string> thumbGeos = cvcapp.listProperty("thumbnail.geometries",true);;
      for (int i=0; i<thumbGeos.size(); i++) {
	if (thumbGeos[i].compare("thumbnail_isocontour") == 0) {
	  ui.ThumbnailIsoCheckBox->setChecked(true);
	}
      }
      
      std::vector<std::string> zoomedGeos = cvcapp.listProperty("zoomed.geometries",true);
      for (int i=0; i<zoomedGeos.size(); i++) {
	if (zoomedGeos[i].compare("zoomed_isocontour") == 0) {
	  ui.ZoomedIsoCheckBox->setChecked(true);
	}
      }

      if (cvcapp.properties("thumbnail_isocontouring.method").compare("FastContouring")==0) {
	ui.ThumbnailMethodComboBox->setCurrentIndex(1);
      }

      if (cvcapp.properties("zoomed_isocontouring.method").compare("FastContouring")==0) {
	ui.ZoomedMethodComboBox->setCurrentIndex(1);
      }
    }

    if(dialog.exec() == QDialog::Accepted) { 
      
      string thumbMode = ui.ThumbnailViewCombo->currentText().toStdString();

      if (thumbMode.compare("Solid") == 0)
	cvcapp.properties("thumbnail.geometry_rendering_mode", "solid");
      else if (thumbMode.compare("Wireframe") == 0)
	cvcapp.properties("thumbnail.geometry_rendering_mode", "wireframe");
      else if (thumbMode.compare("Filled Wireframe") == 0)
	cvcapp.properties("thumbnail.geometry_rendering_mode", "filled_wireframe");
      else if (thumbMode.compare("Flat") == 0)
	cvcapp.properties("thumbnail.geometry_rendering_mode", "flat");
      else if (thumbMode.compare("Flat Filled Wireframe") == 0)
	cvcapp.properties("thumbnail.geometry_rendering_mode", "flat_wireframe");

      string zoomedMode = ui.ZoomedViewCombo->currentText().toStdString();
      if (zoomedMode.compare("Solid") == 0)
	cvcapp.properties("zoomed.geometry_rendering_mode", "solid");
      else if (zoomedMode.compare("Wireframe") == 0)
	cvcapp.properties("zoomed.geometry_rendering_mode", "wireframe");
      else if (zoomedMode.compare("Filled Wireframe") == 0)
	cvcapp.properties("zoomed.geometry_rendering_mode", "filled_wireframe");
      else if (zoomedMode.compare("Flat") == 0)
	cvcapp.properties("zoomed.geometry_rendering_mode", "flat");
      else if (zoomedMode.compare("Flat Filled Wireframe") == 0)
	cvcapp.properties("zoomed.geometry_rendering_mode", "flat_wireframe");


      cvcapp.properties("thumbnail_isocontouring.method",ui.ThumbnailMethodComboBox->currentText().toStdString());
      cvcapp.properties("zoomed_isocontouring.method",ui.ZoomedMethodComboBox->currentText().toStdString());

    }

    // arand, 5-3-2011: simplified list property interface
    if (ui.ThumbnailIsoCheckBox->isChecked() ) {
      cvcapp.listPropertyAppend("thumbnail.geometries", "thumbnail_isocontour");
    } else {
      cvcapp.listPropertyRemove("thumbnail.geometries", "thumbnail_isocontour");
    }

    if (ui.ZoomedIsoCheckBox->isChecked() ) {
      cvcapp.listPropertyAppend("zoomed.geometries", "zoomed_isocontour");
    } else {
      cvcapp.listPropertyRemove("zoomed.geometries", "zoomed_isocontour");
    }

  }
  


  // arand, implemented 4-12-2011
  void CVCMainWindow::toggleWireCubeSlot() {
    //_viewWireCubeMenuAction->setChecked(true);
    if (!_viewWireCubeMenuAction->isChecked() ) {
      cvcapp.properties("zoomed.draw_subvolume_selector", "false");
      cvcapp.properties("zoomed.draw_bounding_box", "false");

      cvcapp.properties("thumbnail.draw_subvolume_selector", "false");
      cvcapp.properties("thumbnail.draw_bounding_box", "false");

      _viewWireCubeMenuAction->setChecked(false);
    } else {

      cvcapp.properties("zoomed.draw_subvolume_selector", "false");
      cvcapp.properties("zoomed.draw_bounding_box", "true");

      cvcapp.properties("thumbnail.draw_subvolume_selector", "true");
      cvcapp.properties("thumbnail.draw_bounding_box", "true");


      _viewWireCubeMenuAction->setChecked(true);
    }
  }

  // arand, implemented 4-12-2011
  void CVCMainWindow::toggleClipBBoxSlot() {
    //_clipBBoxMenuAction->setChecked(true);
    if (!_clipBBoxMenuAction->isChecked()) {
      cvcapp.properties("zoomed.clip_geometry", "false");
      cvcapp.properties("thumbnail.clip_geometry", "false");

      _clipBBoxMenuAction->setChecked(false);
    } else {
      cvcapp.properties("zoomed.clip_geometry", "true");
      cvcapp.properties("thumbnail.clip_geometry", "true");

      _clipBBoxMenuAction->setChecked(true);
    }
  }
  


  void CVCMainWindow::customEvent(QEvent *event)
  {
    CVCMainWindowEvent *mwe = dynamic_cast<CVCMainWindowEvent*>(event);
    if(!mwe) return;
    
    if(mwe->name == "initializePropertyMapWidget")
      initializePropertyMapWidget();
    else if(mwe->name == "initializeDataMapWidget")
      initializeDataMapWidget();
    else if(mwe->name == "initializeThreadMapWidget")
      initializeThreadMapWidget();
    else if(mwe->name == "setupGridRover")
      setupGridRover();
    else if(mwe->name == "openFileFailed")
      QMessageBox::critical(this,"Error Loading Files",
    			    QString("Error loading some files."));       // TODO: list files
      /*
      QMessageBox::critical(this,"Error Loading Files",
    			    QString("Error loading the following files: %1")
			    .arg(failedToLoad.join(" :: ")));       
      */
    else if(mwe->name == "logEntry") {
      try {
        string msg = boost::any_cast<string>(mwe->data);
        _ui->_logTextEdit->insertPlainText(QString::fromStdString(msg));
      }
      catch(const boost::bad_any_cast &) {
        throw std::logic_error("logEntry event data must be a std::string");
      }
    }
  }

  void CVCMainWindow::initializePropertyMapWidget()
  {
    _ui->_propertyMap->setSortingEnabled(false);
    _ui->_propertyMap->clear();

    //re-initialize the TreeWidget with the newest contents of the propertymap
    PropertyMap map = cvcapp.properties();
    QList<QTreeWidgetItem*> items;
    for (const auto& val : map) {    
      QStringList slval;
      slval.append(QString::fromStdString(val.first));
      slval.append(QString::fromStdString(val.second));
      
      QTreeWidgetItem *twi = new QTreeWidgetItem((QTreeWidget*)0, slval);
      twi->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable | Qt::ItemIsEditable);
      items.append(twi);
    }
    
    _ui->_propertyMap->insertTopLevelItems(0,items);
    _ui->_propertyMap->setSortingEnabled(true);
    _ui->_propertyMap->sortByColumn(0,Qt::AscendingOrder);
  }

  void CVCMainWindow::initializeDataMapWidget()
  {
    _dataMap->setSortingEnabled(false);
    _dataMap->clear();

    //re-initialize the TreeWidget with the newest contents of the datamap
    DataMap map = cvcapp.data();
    QList<QTreeWidgetItem*> items;
    for (auto& val : map) {    
      QStringList slval;
      slval.append(QString::fromStdString(val.first));
      slval.append(QString(QString::fromStdString(cvcapp.dataTypeName(val.first))));
      
      QTreeWidgetItem *twi = new QTreeWidgetItem((QTreeWidget*)0, slval);
      twi->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
      items.append(twi);
    }

    _dataMap->insertTopLevelItems(0,items);
    _dataMap->setSortingEnabled(true);
    _dataMap->sortByColumn(0,Qt::AscendingOrder);
  }

  // 07/29/2011 - transfix - adding fourth column for thread info string
  void CVCMainWindow::initializeThreadMapWidget()
  {
    _ui->_threadMap->setSortingEnabled(false);
    _ui->_threadMap->clear();

    //re-initialize the TreeWidget with the newest contents of the thread map
    ThreadMap map = cvcapp.threads();
    QList<QTreeWidgetItem*> items;
    for (const auto& val : map) {    
      QStringList slval;

      //First column
      slval.append(QString::fromStdString(val.first));

      //Second column
      ThreadPtr tp = cvcapp.threads(val.first);
      std::stringstream ss;
      if(tp)
        ss << tp->get_id();
      else
        ss << "null";
      slval.append(QString::fromStdString(ss.str()));

      //Third column
      double progress = cvcapp.threadProgress(val.first);
      slval.append(QString("%1%").arg(progress*100.0));

      //Fourth column
      QString info = QString::fromStdString(cvcapp.threadInfo(val.first));
      slval.append(info);
      
      QTreeWidgetItem *twi = new QTreeWidgetItem((QTreeWidget*)0, slval);
      twi->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
      items.append(twi);
    }    

    _ui->_threadMap->insertTopLevelItems(0,items);
    _ui->_threadMap->setSortingEnabled(true);
    _ui->_threadMap->sortByColumn(0,Qt::AscendingOrder);
  }


  // future: fix this so it behaves nicely when several volumes are loaded...
  void CVCMainWindow::setupGridRover() {
    // set volume for VolumeGridRover
#ifdef USING_VOLUMEGRIDROVER
     bool hasSource = false;
     std::string volSelected;

     CVC_NAMESPACE::DataMap map = cvcapp.data();
     for (const auto& val : map) { 
       std::string myname(val.first);
       // only deal with files for now...
       if (cvcapp.isData<VolMagick::VolumeFileInfo>(myname)) {
           volSelected = myname;
           hasSource = true;
	   break;
       }
     }
     
     // make sure data is loaded
     if ( hasSource ) {
	VolMagick::VolumeFileInfo vfi = boost::any_cast<VolMagick::VolumeFileInfo>(map[volSelected]); 
	volumeGridRover()->setVolume( vfi );
     }
     
#endif    
  }


  void CVCMainWindow::propertyMapChanged(const std::string& key)
  {
    std::stringstream ss;
    ss << BOOST_CURRENT_FUNCTION << ": " << key << " changed" << std::endl;
    cvcapp.log(5,ss.str());


    // arand, 2/17/2012
    //   This is kind of a hack to close the program via the property map.
    //   This gives a cleaner way for the "batch mode" to exit the program.
    if (key == "exitVolumeRover") {
      // this seems more reliable than just calling close()
      //  but I am not enough of a qt expert to know what is going on...
      QCoreApplication::quit();
      return;
    }

    //This boost::signals2 slot might get called from another thread so it is important
    //that we do a thread safe posting of an event instead of calling the initialization
    //function directly.
    QCoreApplication::postEvent(this,
      new CVCMainWindowEvent("initializePropertyMapWidget")
    );
  }

  void CVCMainWindow::dataMapChanged(const std::string& key)
  {
    std::stringstream ss;
    ss << BOOST_CURRENT_FUNCTION << ": " << key << std::endl;
    cvcapp.log(5,ss.str());

    QCoreApplication::postEvent(this,
      new CVCMainWindowEvent("initializeDataMapWidget")
    );
  }

  void CVCMainWindow::threadMapChanged(const std::string& key)
  {
    std::stringstream ss;
    ss << BOOST_CURRENT_FUNCTION << ": " << key << std::endl;
    cvcapp.log(5,ss.str());

    QCoreApplication::postEvent(this,
      new CVCMainWindowEvent("initializeThreadMapWidget")
    );    
  }

  // anisotropicDiffusion slot
  // arand: implemented 4-1-2011
  // transfix: refactored, added support for Volume 4-24-2011
  // arand: moved to AnisotropicDiffusionDialog.cpp, 4-28-2011
  //        I think this is a good way to (1) keep CVCMainWindow.cpp smaller
  //          and (2) implement slots/signals that are local to the dialog.
  //        If there is a better way, let me know before I refactor all the
  //          functionality.  Also, maybe there is a "good" way to integrate
  //          threading into the framework.
  void CVCMainWindow::anisotropicDiffusionSlot() {
    AnisotropicDiffusionDialog dialog;
    if(dialog.exec() == QDialog::Accepted) {
      dialog.RunAnisotropicDiffusion();
    }
  }

  // arand: contrastEnhancementSlot was rewritten on 3-31-2011
  // arand: improved implementation in separate cpp file, 5-2-2011
  void CVCMainWindow::contrastEnhancementSlot() {
    ContrastEnhancementDialog dialog;
    if(dialog.exec() == QDialog::Accepted) {
      dialog.RunContrastEnhancement();
    }
  }

  // arand, 3-31-3011: bilateralFilterSlot was rewritten
  // arand,  5-2-2011: external class written to improve interface 
  //                   and use threads
  void CVCMainWindow::bilateralFilterSlot() {
    BilateralFilterDialog dialog;
    if(dialog.exec() == QDialog::Accepted) {
      dialog.RunBilateralFilter();
    }    
  }

  // arand, function converted to Qt4, 4-1-2011
  // arand, updated to new style, 5-2-2011
  void CVCMainWindow::gdtvFilterSlot() {
    GDTVFilterDialog dialog;
    if(dialog.exec() == QDialog::Accepted) {
      dialog.RunGDTVFilter();
    }
  }

  // cha, function converted to Qt4, 4-1-11
  void CVCMainWindow::virusSegmentationSlot()
  {
#ifdef USING_SEGMENTATION
     Ui::SegmentVirusMapDialog ui;
     QDialog dialog;
     ui.setupUi(&dialog);

     bool hasSource = false;
     DataMap map = cvcapp.data();
     for (const auto& val : map) { 
       //std::cout << val.first << " " << val.second.type().name() << std::endl;
       std::string myname(val.first);
       std::string mytype(val.second.type().name());
       // only deal with files for now...
       if (cvcapp.isData<VolMagick::VolumeFileInfo>(myname)) {           
	   ui.m_VolumeList->addItem(QString::fromStdString(myname));
           hasSource = true;
       }
     }

	// make sure data is loaded
	if ( hasSource ) {
                
		//set up dialog with saved settings
		if (dialog.exec() == QDialog::Accepted) {
		
                        std::string volSelected = ui.m_VolumeList->currentText().toStdString();
                        VolMagick::VolumeFileInfo vfi = boost::any_cast<VolMagick::VolumeFileInfo>(map[volSelected]);

			if( ui.m_LocalSegmentationButton->isChecked() )
			{
				if( m_LocalSegThread )
				{
                                        if( m_LocalSegThread->isRunning() )
						QMessageBox::critical(this,"Error","Local Segmentation thread already running!",QMessageBox::Ok,Qt::NoButton,Qt::NoButton);
					else
					{
						delete m_LocalSegThread;
						m_LocalSegThread = new LocalSegThread(vfi.filename(),&ui,this);
						m_LocalSegThread->start();
					}
				}
				else
				{
					m_LocalSegThread = new LocalSegThread(vfi.filename(),&ui,this);
					m_LocalSegThread->start();
				}
			}
                        else
			{
				if(m_RemoteSegThread && m_RemoteSegThread->isRunning())
					QMessageBox::critical(this,"Error","Remote Segmentation thread already running!",QMessageBox::Ok,Qt::NoButton,Qt::NoButton);
				else
				{
					if(m_RemoteSegThread) delete m_RemoteSegThread;
					m_RemoteSegThread = new RemoteSegThread(&ui,this);
					m_RemoteSegThread->start();
				}
			}
		}	
	}
	else {
		QMessageBox::warning(this, tr("Segment Virus Map"),
			tr("A volume must be loaded for this feature to work."), 0,0);
	}
#else 
	QMessageBox::information(this, tr("Virus Segmentation"),
				 tr("Virus segmentation has been disabled in this VolumeRover build."), QMessageBox::Ok);

#endif
   }

#ifdef USING_SECONDARYSTRUCTURES
  void CVCMainWindow::uploadgeometrySlot()
  {
    int alphaValue = SecondaryStructureHistogramUi->m_SecondaryStructureAlphaSpinBox->value();
    int betaValue = SecondaryStructureHistogramUi->m_SecondaryStructureBetaSpinBox->value();

    float alphaMinWidth = m_AlphaHistogram->getMinWidth();
    float alphaMaxWidth = m_AlphaHistogram->getMaxWidth();
    float betaMinWidth =  m_BetaHistogram->getMinWidth();
    float betaMaxWidth =  m_BetaHistogram->getMaxWidth();
    bool alphahistogramChanged = SecondaryStructureHistogramUi->m_RenderAlphaCheckBox->isChecked();
    bool betahistogramChanged = SecondaryStructureHistogramUi->m_RenderBetaCheckBox->isChecked();

    cout<<"min, max, min, max" << alphaMinWidth <<" " << alphaMaxWidth <<" " << betaMinWidth <<" " << betaMaxWidth<<" " << alphahistogramChanged << "  " << betahistogramChanged<< endl;

    m_Skeleton->buildAllGeometry(alphaValue, betaValue, alphaMinWidth, alphaMaxWidth, betaMinWidth, betaMaxWidth, alphahistogramChanged, betahistogramChanged);
    
  	if(m_Skeleton->helixGeom !=0)
	{
	  string ptgeom = "helix_Geometry";
	  cvcapp.data(ptgeom, *(m_Skeleton->helixGeom));
      cvcapp.properties("thumbnail.geometries", cvcapp.properties("thumbnail.geometries")+","+ptgeom);
      cvcapp.properties("zoomed.geometries", cvcapp.properties("zoomed.geometries")+","+ptgeom);
	}

  	if(m_Skeleton->sheetGeom !=0)
	{
	  string ptgeom = "sheet_Geometry";
	  cvcapp.data(ptgeom, *(m_Skeleton->sheetGeom));
      cvcapp.properties("thumbnail.geometries", cvcapp.properties("thumbnail.geometries")+","+ptgeom);
      cvcapp.properties("zoomed.geometries", cvcapp.properties("zoomed.geometries")+","+ptgeom);
	}

  	if(m_Skeleton->skelGeom !=0)
	{
	  string ptgeom = "skeleteon_Geometry";
	  cvcapp.data(ptgeom, *(m_Skeleton->skelGeom));
      cvcapp.properties("thumbnail.geometries", cvcapp.properties("thumbnail.geometries")+","+ptgeom);
      cvcapp.properties("zoomed.geometries", cvcapp.properties("zoomed.geometries")+","+ptgeom);
	}

  	if(m_Skeleton->curveGeom !=0)
	{
	  string ptgeom = "curve_Geometry";
	  cvcapp.data(ptgeom, *(m_Skeleton->curveGeom));
      cvcapp.properties("thumbnail.geometries", cvcapp.properties("thumbnail.geometries")+","+ptgeom);
      cvcapp.properties("zoomed.geometries", cvcapp.properties("zoomed.geometries")+","+ptgeom);
	}
  }

void CVCMainWindow::showHistogramDialogSlot() {
     if( hasSecondaryStructure )
        histogram_dialog->show();
     else
        QMessageBox::information(this, tr("Secondary Structure Histogram"),
			   tr("You need to run ""Secondary Structure Elucidation"" first"), QMessageBox::Ok);

  }
#endif


  void CVCMainWindow::secondaryStructureSlot() {
  //  Ui::SecondaryStructureDialog ui;
    SecondaryStructureUi = new Ui::SecondaryStructureDialog();
    QDialog dialog;
    SecondaryStructureUi->setupUi(&dialog);
   
    SecondaryStructureUi->HelixWidthEdit->insert("2.0");
    SecondaryStructureUi->MinHelixWidthRatioEdit->insert("0.001");
    SecondaryStructureUi->MaxHelixWidthRatioEdit->insert("8.0");
    SecondaryStructureUi->MinHelixLengthEdit->insert("1.8");

    SecondaryStructureUi->SheetWidthEdit->insert("2.6");
    SecondaryStructureUi->MinSheetWidthRatioEdit->insert("0.01");
    SecondaryStructureUi->MaxSheetWidthRatioEdit->insert("8.0");
    SecondaryStructureUi->SheetExtendedEdit->insert("1.5");

    SecondaryStructureUi->ThresholdEdit->insert("125");


    DataMap map = volrover_app_instance().data();

    bool Isgeom = 0;
    for (const auto& val : map) { 
    //  std::cout << val.first << " " << val.second.type().name() << std::endl;
      std::string myname(val.first);
      std::string mytype(val.second.type().name());
      // only deal with files for now...
      if (mytype.compare("N9VolMagick14VolumeFileInfoE") == 0) {
		  SecondaryStructureUi->VolumeList->addItem(QString::fromStdString(myname));
		  SecondaryStructureUi->InputList->addItem(QString::fromStdString(myname));
      }
	  if(mytype.compare("N15cvcraw_geometry9cvcgeom_tE") == 0)
	  {
	  	  SecondaryStructureUi->InputList->addItem(QString::fromStdString(myname));
		  Isgeom = 1;
	  }

    }
 
         
     if(dialog.exec() == QDialog::Accepted) {
      // do stuff...
      if ( SecondaryStructureUi->tabWidget->currentIndex() == 0) {
#ifdef USING_SEGMENTATION
	// run Zeyun's method here ui.IterationsEdit->displayText().toInt();

	XmlRpc::XmlRpcValue m_Params;
	std::string volSelected =  SecondaryStructureUi->VolumeList->currentText().toStdString();
	VolMagick::VolumeFileInfo vfi = boost::any_cast<VolMagick::VolumeFileInfo>(map[volSelected]);

	m_Params[0] = vfi.filename(); 
	m_Params[1] =  SecondaryStructureUi->HelixWidthEdit->text().toDouble();
	m_Params[2] =  SecondaryStructureUi->MinHelixWidthRatioEdit->text().toDouble();
	m_Params[3] =  SecondaryStructureUi->MaxHelixWidthRatioEdit->text().toDouble();
	m_Params[4] =  SecondaryStructureUi->MinHelixLengthEdit->text().toDouble();
	m_Params[5] =  SecondaryStructureUi->SheetWidthEdit->text().toDouble();
	m_Params[6] =  SecondaryStructureUi->MinSheetWidthRatioEdit->text().toDouble();
	m_Params[7] =  SecondaryStructureUi->MaxSheetWidthRatioEdit->text().toDouble();
	m_Params[8] =  SecondaryStructureUi->SheetExtendedEdit->text().toDouble();
	if ( SecondaryStructureUi->ThresholdCheckbox->isChecked())
	  m_Params[9] =  SecondaryStructureUi->ThresholdEdit->text().toDouble();
	else m_Params[9] = double( -1.0);
	XmlRpc::XmlRpcValue result;
	if(secondaryStructureDetection(m_Params,result))
	{
		 std::string basename =  QFileInfo(QString::fromStdString(vfi.filename())).path().toStdString() +"/"+
		 	QFileInfo(QString::fromStdString(vfi.filename())).baseName().toStdString() ;
		
		 std::string helixC = basename +"_helix_c.rawc";
		 std::string helixS = basename  + "_helix_s.rawc";
		 std::string sheet =  basename  + "_sheet.rawc";
		 std::string skeleton = basename  + "_skeleton.rawc";
      
	  QStringList failedToLoad;
	  if(!cvcapp.readData(helixC))
          failedToLoad.push_back(QString::fromStdString(helixC));

      if(!cvcapp.readData(helixS))
          failedToLoad.push_back(QString::fromStdString(helixS));
         
	  if(!cvcapp.readData(sheet))
          failedToLoad.push_back(QString::fromStdString(sheet));
	  
	  if(!cvcapp.readData(skeleton))
          failedToLoad.push_back(QString::fromStdString(skeleton));


 
      if(!failedToLoad.isEmpty())
         QMessageBox::critical(this,"Error Loading Files",
                            QString("Error loading the following files: %1")
                            .arg(failedToLoad.join(" :: ")));

	}
	else
	  cout << "Failure" << endl;

	// TODO: run above as a thread
	//       get data from memory, not a file...

#else
  QMessageBox::information(this, tr("Secondary Structure Elucidation"),
			   tr("Volume based secondary structure elucidation has been disabled in this VolumeRover build.  You must enable the segmentation library in the build process for this method."), QMessageBox::Ok);
#endif
       }
      
      else if ( SecondaryStructureUi->tabWidget->currentIndex() == 1) {

#ifdef USING_SECONDARYSTRUCTURES

	CVCGEOM_NAMESPACE::cvcgeom_t geom;
      // handle geometry input
	  if(Isgeom == 1){
	 	 std::string geomSelected =  SecondaryStructureUi->InputList->currentText().toStdString();
	     geom= boost::any_cast<CVCGEOM_NAMESPACE::cvcgeom_t>(map[geomSelected]);
	  }
	  else  //handle volume input 
	  { 
	  	std::string volSelected =  SecondaryStructureUi->InputList->currentText().toStdString();
	  	VolMagick::VolumeFileInfo vfi = boost::any_cast<VolMagick::VolumeFileInfo>(map[volSelected]);
	 } 

    //Computation

    //vector<cvcraw_geometry::cvcgeom_t> geoms;
	m_Skeleton = new Skel();
	m_Skeleton->compute_secondary_structures(&geom);

 	HistogramData alphaData = m_Skeleton->get_alpha_histogram_data();
	HistogramData betaData = m_Skeleton->get_beta_histogram_data();

        if( SecondaryStructureHistogramUi == NULL )
          SecondaryStructureHistogramUi = new Ui::HistogramDialog();
        if( histogram_dialog == NULL )
          histogram_dialog = new QDialog();

        SecondaryStructureHistogramUi->setupUi(histogram_dialog);
        
        QSizePolicy sizePolicy(static_cast<QSizePolicy::Policy>(7), static_cast<QSizePolicy::Policy>(7));
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth( SecondaryStructureHistogramUi->m_AlphaHistogram->sizePolicy().hasHeightForWidth());

        m_AlphaHistogram = new Histogram(SecondaryStructureHistogramUi->m_AlphaHistogram);
        m_BetaHistogram = new Histogram(SecondaryStructureHistogramUi->m_BetaHistogram);

        sizePolicy.setHeightForWidth( m_AlphaHistogram->sizePolicy().hasHeightForWidth());
        m_AlphaHistogram->setSizePolicy( sizePolicy );
        m_AlphaHistogram->resize( SecondaryStructureHistogramUi->m_AlphaHistogram->frameSize());

        sizePolicy.setHeightForWidth( m_BetaHistogram->sizePolicy().hasHeightForWidth());
        m_BetaHistogram->setSizePolicy( sizePolicy );
        m_BetaHistogram->resize( SecondaryStructureHistogramUi->m_BetaHistogram->frameSize());

        SecondaryStructureHistogramUi->m_SecondaryStructureAlphaSpinBox->setValue(m_Skeleton->getDefaultAlphaCount());
        SecondaryStructureHistogramUi->m_SecondaryStructureBetaSpinBox->setValue(m_Skeleton->getDefaultBetaCount());

        m_AlphaHistogram->setData(alphaData);
        m_BetaHistogram->setData(betaData);

    	connect(SecondaryStructureHistogramUi->m_PushButtonUpdate, SIGNAL(clicked()), this, SLOT(uploadgeometrySlot()));
	    uploadgeometrySlot(); 	
	
        histogram_dialog->show();
        hasSecondaryStructure = true;

#else
  QMessageBox::information(this, tr("Secondary Structure Elucidation"),
			   tr("Surface based secondary structure elucidation has been disabled in this VolumeRover build."), QMessageBox::Ok);

#endif

       }
    }
  }

#ifdef USING_RECONSTRUCTION
  void CVCMainWindow::reconstructionSlot() {
    ReconstructionDialog dialog;
    if(dialog.exec() == QDialog::Accepted) {
      dialog.RunReconstruction();
    }
  }
#endif


  // FUTURE: set up a dialog and make this actually work
  // 10/09/2011 -- transfix -- Exception handling
  void CVCMainWindow::raytraceSlot() {
    try
      {
        ofstream fout("test.cnf");

        // combined rendering mode (?)
        fout << 8 << endl;
   
        // get the file info
        string myvol = cvcapp.properties("viewers.vfi");

        VolMagick::VolumeFileInfo vfi = cvcapp.data<VolMagick::VolumeFileInfo>(myvol);
        fout << vfi.XSpan() << " " << vfi.YSpan() << " " << vfi.ZSpan() << endl;

        fout << vfi.XDim() << " " << vfi.YDim() << " " << vfi.ZDim() << endl;

        // arand: fix this if you want a subvolume
        fout << 0 << " " << 0 << " " << 0 << " "
             << vfi.XDim() << " " << vfi.YDim() << " " << vfi.ZDim() << endl;

        //fout << vfi.XMin() << " " << vfi.YMin() << " " << vfi.ZMin() << " " 
        //	 << vfi.XMax() << " " << vfi.YMax() << " " << vfi.ZMax() << endl;
    
        fout << endl;

        // arand: need viewing parameters here...


        // no materials
        fout << "0" << endl << endl;

        // no isocontours
        fout << "0" << endl << endl;
    
        boost::shared_array<unsigned char> trans =
          cvcapp.data<boost::shared_array<unsigned char> >("viewers_transfer_function");

        fout << 256 << endl;
        for (int i=0; i<256; i++) {
          fout << i << " " << (int)trans[4*i] << " " << (int)trans[4*i+1] << " " 
               << (int)trans[4*i+2] << " " << (int)trans[4*i+3] << endl; 
        }

        // add info for lights...
    
        // no cutting planes
        fout << "0" << endl << endl;

        // rendering modes 
        fout << "1 0 1 3" << endl << endl;

        // misc parameters 
        fout << "1500" << endl // step size (?) higher is slower
             << "0 20 255" << endl // gradient ramp (?)
             << "255 255 255" << endl;// last is background color
      }
    catch(std::exception& e)
      {
        std::string msg = boost::str(boost::format("%s :: ERROR :: %s\n")
                                     % BOOST_CURRENT_FUNCTION
                                     % e.what());
        cvcapp.log(2,msg);
        QMessageBox::critical(this,"Exception",QString::fromStdString(msg));
      }
  }

  void CVCMainWindow::multiTileServerSlot()
  {
     static MultiTileServerDialog *dialog = new MultiTileServerDialog();
     dialog->updateList();
     dialog->show();
  }


  void CVCMainWindow::hlsSurfaceSlot() {
#ifdef USING_HLEVELSET
    Ui::HLSSurfaceDialog ui;
    QDialog dialog;
    ui.setupUi(&dialog);

    ui.DimensionEdit->insert("128");
    ui.IterationsEdit->insert("1");
    ui.EdgeLengthEdit->insert("1.0");
    ui.MaxDimEdit->insert("100");

   DataMap map = cvcapp.data();
    for (const auto& val : map) { 
      //std::cout << val.first << " " << val.second.type().name() << std::endl;
      std::string myname(val.first);
      std::string mytype(val.second.type().name());
      // only deal with files for now...
      if (cvcapp.isData<cvcraw_geometry::cvcgeom_t>(myname)) {
	ui.GeometryList->addItem(QString::fromStdString(myname));
      }
    }

    if(dialog.exec() == QDialog::Accepted) {

      unsigned int dim[3];
      dim[0] = dim[1] = dim[2] = ui.DimensionEdit->displayText().toInt();
      double it = ui.IterationsEdit->displayText().toInt();
      double edge = ui.EdgeLengthEdit->displayText().toDouble();
      double maxdim = ui.MaxDimEdit->displayText().toInt();
      
      
      // find the volume file to work with
      std::string geomSelected = ui.GeometryList->currentText().toStdString();
      // get the selected geometry... 
      cvcraw_geometry::cvcgeom_t geom = boost::any_cast<cvcraw_geometry::cvcgeom_t>(map[geomSelected]);

      // TODO: actually call the function.
      boost::shared_ptr<Geometry> geometry (new Geometry(Geometry::conv(geom)));
      std::vector<float> vertex_Positions(geometry->m_NumPoints*3);
      for(unsigned int i = 0; i < geometry->m_NumPoints; i++) {
	vertex_Positions[i*3+0] = geometry->m_Points[i*3+0];
	vertex_Positions[i*3+1] = geometry->m_Points[i*3+1];
	vertex_Positions[i*3+2] = geometry->m_Points[i*3+2];
      }

      clock_t t;
      double time;
      t=clock();
      time=(double)t/(double)CLOCKS_PER_SEC;

      HLevelSetNS::HLevelSet* hLevelSet = new HLevelSetNS::HLevelSet();
      boost::tuple<bool,VolMagick::Volume> result = 
	hLevelSet->getHigherOrderLevelSetSurface_Xu_Li(vertex_Positions,dim,edge,it, maxdim);
      delete hLevelSet;

      clock_t t_end;
      double time_end;
      t_end=clock();
      time_end=(double)t_end/(double)CLOCKS_PER_SEC;
      time=time_end-time;
      printf("\n Highlevelset Time = %f \n",time);
      
      //  printf("OK2"); 
      try {
	
	/*
	VolMagick::Volume result_vol(result.get<1>());
	QDir tmpdir(m_CacheDir.absPath() + "/tmp");
	if(!tmpdir.exists()) tmpdir.mkdir(tmpdir.absPath());
	QFileInfo tmpfile(tmpdir,"tmp.rawiv");
	std::stringstream ss;
	ss << "Creating volume " << tmpfile.absFilePath().toStdString().c_str() <<std::endl;
	cvcapp.log(5,ss.str());
	QString newfilename = QDir::convertSeparators(tmpfile.absFilePath());
	VolMagick::createVolumeFile(cvcapp, newfilename.toStdString(),
				    result_vol.boundingBox(),
				    result_vol.dimension(),
				    std::vector<VolMagick::VoxelType>(1, result_vol.voxelType()));
	VolMagick::writeVolumeFile(cvcapp, result_vol,newfilename.toStd.String());
	openFile(newfilename);
	*/
	// TODO: update code above...

      }
      catch(const VolMagick::Exception& e) {
	QMessageBox::critical(this,"Error loading HLS volume",e.what());
	return;
      }


    }


#else
  QMessageBox::information(this, tr("Higher-Order Level Set Surface"),
			   tr("The HLS surface feature has been disabled in this VolumeRover build."), QMessageBox::Ok);

#endif
  }


  // arand: refactor and threading, 5-4-2011
  void CVCMainWindow::curationSlot() {
#ifdef USING_CURATION
    CurationDialog dialog;
    if(dialog.exec() == QDialog::Accepted) {
      dialog.RunCuration();
    }
#else
    QMessageBox::information(this, tr("Curation"),
			     tr("Curation feature has been disabled in this VolumeRover build."), QMessageBox::Ok);
    
#endif

  }

  // arand: refactor and threading, 5-3-2011
  void CVCMainWindow::pocketTunnelSlot() {
#ifdef USING_POCKET_TUNNEL
    PocketTunnelDialog dialog;
    if(dialog.exec() == QDialog::Accepted) {
      dialog.RunPocketTunnel();
    }
#else
  QMessageBox::information(this, tr("Pocket/Tunnel Detection"),
			   tr("Pocket/Tunnel detection feature has been disabled in this VolumeRover build."), QMessageBox::Ok);
#endif
  }

  void CVCMainWindow::tightCoconeSlot() {
#ifdef USING_TIGHT_COCONE
    TightCoconeDialog dialog;
    if(dialog.exec() == QDialog::Accepted) {
      dialog.RunTightCocone();
    }
#else
  QMessageBox::information(this, tr("Tight Cocone"),
			   tr("Tight cocone has been disabled in this VolumeRover build."), QMessageBox::Ok);
#endif

  }

  void CVCMainWindow::superSecondaryStructuresSlot(){
#ifdef USING_SUPERSECONDARY_STRUCTURES
	SuperSecondaryStructuresDialog dialog;
	if(dialog.exec() == QDialog::Accepted)
	{
		dialog.RunSuperSecondaryStructures();
	}
#else
	QMessageBox::information(this, tr("Super Secondary Structures"),
				tr("Super Secondary Structures identification has been disabled in this VolumeRover build."), QMessageBox::Ok);
#endif

  }

  void CVCMainWindow::skeletonizationSlot() {
#ifdef USING_SKELETONIZATION
    SkeletonizationDialog dialog;
    if(dialog.exec() == QDialog::Accepted) {
      dialog.RunSkeletonization();
    }
#else
  QMessageBox::information(this, tr("Skeletonization"),
			   tr("Skeletonization has been disabled in this VolumeRover build."), QMessageBox::Ok);
#endif
  }

  void CVCMainWindow::contourTilerSlot(){
#ifdef USING_TILING
    ContourTilerDialog dialog;
    if(dialog.exec() == QDialog::Accepted) {
      dialog.RunContourTiler();
    }    
#else
  QMessageBox::information(this, tr("ContourTiler"),
			   tr("Contour tiling has been disabled in this VolumeRover build."), QMessageBox::Ok);
#endif
  }

  void CVCMainWindow::LBIE_Slot(){
    Ui::LBIE_dialog ui;
    QDialog dialog;
    ui.setupUi(&dialog);
    QString outputMessage;

    float outer_isoval, inner_isoval;

    // set up default values...
    ui.m_ErrorTolerance->setText("1.2501");
    ui.m_InnerErrorTolerance->setText("0.0001");
    ui.m_Iterations->setText("0");
    DataMap map = cvcapp.data();

    bool found = false;
    for (const auto& val : map) {
      std::string myname(val.first);
      std::string mytype(val.second.type().name());
      if(cvcapp.isData<VolMagick::VolumeFileInfo>(myname)) {
	ui.volumeList->addItem(QString::fromStdString(myname));
	found = true;
      }
    }
    
    if (!found) {
      QMessageBox::information(this, "LBIE Meshing", "No Volume Loaded", QMessageBox::Ok);
      return;
    }

    // find the volume files to work with
    std::string volSelected = ui.volumeList->currentText().toStdString();

    // read in the data if necessary
    VolMagick::VolumeFileInfo vfi = boost::any_cast<VolMagick::VolumeFileInfo>(map[volSelected]);

    CVCColorTable::ColorTable::color_table_info cti = boost::any_cast<CVCColorTable::ColorTable::color_table_info>(map["viewers_color_table_info"]);
    if (cti.isocontourNodes().size() <1) {
      ui.m_OuterIsoValue->setText("0.0");
      ui.m_InnerIsoValue->setText("0.0");
    }else{
      std::vector<double> values;
      for (const auto& node : cti.isocontourNodes()) {
	//cout << "Iso " << node.position << endl;
	values.push_back(vfi.min() + node.position*(vfi.max()-vfi.min()));
      }
      QString outerIsovalStr, innerIsovalStr;
      outerIsovalStr.setNum(*std::min_element(values.begin(),values.end()));
      innerIsovalStr.setNum(*std::max_element(values.begin(),values.end()));
      
      ui.m_OuterIsoValue->setText(outerIsovalStr);
      ui.m_InnerIsoValue->setText(innerIsovalStr);
    }
    
    if(dialog.exec() != QDialog::Accepted) { return; }
    
    outer_isoval = ui.m_OuterIsoValue->text().toDouble();
    inner_isoval = ui.m_InnerIsoValue->text().toDouble();
    
    string resultName = ui.OutputEdit->displayText().toStdString();
    if (resultName.empty()) {
      resultName = volSelected + "_LBIE";
    }
    string meshName = resultName + "_Mesh";
   
#ifdef USING_SWEETMESH
    CVCGEOM_NAMESPACE::cvcgeom_t geometry;
    sweetMesh::hexMesh hMesh;
    sweetMesh::runLBIE(vfi, outer_isoval, inner_isoval, ui.m_ErrorTolerance->text().toDouble(), ui.m_InnerErrorTolerance->text().toDouble(), LBIE::Mesher::MeshType(ui.m_MeshType->currentIndex()), LBIE::Mesher::NormalType(ui.m_NormalType->currentIndex()), ui.m_Iterations->text().toUInt(), outputMessage, geometry, hMesh);
    cvcapp.data(resultName, geometry);
    cvcapp.listPropertyAppend("thumbnail.geometries", resultName);
    cvcapp.listPropertyAppend("zoomed.geometries", resultName);
    if(LBIE::Mesher::MeshType(ui.m_MeshType->currentIndex() == LBIE::geoframe::HEXA)){
      cvcapp.data(meshName, hMesh);
      cvcapp.listPropertyAppend("thumbnail.geometries", meshName);
      cvcapp.listPropertyAppend("zoomed.geometries", meshName);
    }
#endif
    
//     VolMagick::Volume vol;
//     readVolumeFile(cvcapp, vol,vfi.filename());
//     LBIE::Mesher mesher(outer_isoval, inner_isoval, ui.m_ErrorTolerance->text().toDouble(), ui.m_InnerErrorTolerance->text().toDouble(), LBIE::Mesher::MeshType(ui.m_MeshType->currentItem()), LBIE::Mesher::GEO_FLOW, LBIE::Mesher::NormalType(ui.m_NormalType->currentItem()), LBIE::Mesher::DUALLIB, false);
//     mesher.extractMesh(vol);
//     mesher.qualityImprove(ui.m_Iterations->text().toUInt());
// 
//     switch(mesher.mesh().mesh_type){
//       case LBIE::geoframe::DOUBLE:{
// 	for(unsigned int i=0; i<mesher.mesh().numtris; i++){
// 	  CVCGEOM_NAMESPACE::cvcgeom_t::triangle_t newTri;
// 	  newTri[0] = mesher.mesh().triangles[i][0];
// 	  newTri[1] = mesher.mesh().triangles[i][1];
// 	  newTri[2] = mesher.mesh().triangles[i][2];
// 	  geometry.triangles().push_back(newTri);
// 	}
// 	outputMessage = "Total Vertices: " + QString::number(mesher.mesh().numverts) + "\nTotal Triangles: " + QString::number(mesher.mesh().numtris);
// 	#ifdef USING_SWEETMESH
// 	std::list<double> aspectRatios;
// 	for(unsigned int i=0; i<mesher.mesh().numtris; i++){
// 	  sweetMesh::sweetMeshVertex v0, v1, v2;
// 	  v0.set(geometry.points()[geometry.triangles()[i][0]][0], geometry.points()[geometry.triangles()[i][0]][1], geometry.points()[geometry.triangles()[i][0]][2]);
// 	  v1.set(geometry.points()[geometry.triangles()[i][1]][0], geometry.points()[geometry.triangles()[i][1]][1], geometry.points()[geometry.triangles()[i][1]][2]);
// 	  v2.set(geometry.points()[geometry.triangles()[i][2]][0], geometry.points()[geometry.triangles()[i][2]][1], geometry.points()[geometry.triangles()[i][2]][2]);
// 	  sweetMesh::triangle tri(v0, v1, v2);
// 	  aspectRatios.push_back(tri.aspectRatio());
// 	}
// 	double lowest = 2000000.0;
// 	double highest = -20.0;
// 	double average = 0.0;
// 	unsigned int bin1, bin2, bin3, bin4, bin5, bin6, bin7, bin8, bin9, bin10, bin11, bin12;
// 	bin1 = bin2 = bin3 = bin4 = bin5 = bin6 = bin7 = bin8 = bin9 = bin10 = bin11 = bin12 = 0;
// 	aspectRatios.sort();
// 	for(std::list<double>::iterator itr=aspectRatios.begin(); itr!=aspectRatios.end(); itr++){
// 	  average += *itr;
// 	  lowest = std::min(lowest, *itr);
// 	  highest = std::max(highest, *itr);
// 	  if(*itr < 1.5){ bin1++; }
// 	  else{
// 	    if(*itr < 2.0){ bin2++; }
// 	    else{
// 	      if(*itr < 2.5){ bin3++; }
// 	      else{
// 		if(*itr < 3.0){ bin4++; }
// 		else{
// 		  if(*itr < 4.0){ bin5++; }
// 		  else{
// 		    if(*itr < 6.0){ bin6++; }
// 		    else{
// 		      if(*itr < 10.0){ bin7++; }
// 		      else{
// 			if(*itr < 15.0){ bin8++; }
// 			else{
// 			  if(*itr < 25){ bin9++; }
// 			  else{
// 			    if(*itr < 50){ bin10++; }
// 			    else{
// 			      if(*itr < 100){ bin11++; }
// 			      else{
// 				bin12++;
// 			      }
// 			    }
// 			  }
// 			}
// 		      }
// 		    }
// 		  }
// 		}
// 	      }
// 	    }
// 	  }
// 	}
// 	outputMessage = outputMessage + "\nAverage Aspect Ratio: " + QString::number(average/(double)aspectRatios.size()) + "\nMinimal Aspect Ratio: " + QString::number(lowest) + "\nMaximal Aspect Ratio: " + QString::number(highest) + "\nAspect Ratio Histogram:"
// 	+ "\n      < 1.5  :  " + QString::number(bin1) + "\t|    6 - 10   :  " + QString::number(bin7)
// 	+ "\n  1.5 - 2.0  :  " + QString::number(bin2) + "\t|   10 - 15   :  " + QString::number(bin8)
// 	+ "\n  2.0 - 2.5  :  " + QString::number(bin3) + "\t|   15 - 25   :  " + QString::number(bin9)
// 	+ "\n  2.5 - 3.0  :  " + QString::number(bin4) + "\t|   25 - 50   :  " + QString::number(bin10)
// 	+ "\n  3.0 - 4.0  :  " + QString::number(bin5) + "\t|   50 - 100  :  " + QString::number(bin11)
// 	+ "\n  4.0 - 6.0  :  " + QString::number(bin6) + "\t|  100 -      :  " + QString::number(bin12);
// #endif
//       cvcapp.data(resultName, geometry);
//       cvcapp.listPropertyAppend("thumbnail.geometries", resultName);
//       cvcapp.listPropertyAppend("zoomed.geometries", resultName);
//       }break;
//       
//       case LBIE::geoframe::TETRA2:{
// 	for(unsigned int i=0; i<mesher.mesh().numtris/4; i++){
// 	  CVCGEOM_NAMESPACE::cvcgeom_t::line_t newLine;
// 	  newLine[0] = mesher.mesh().triangles[4*i][0];
// 	  newLine[1] = mesher.mesh().triangles[4*i][1];
// 	  geometry.lines().push_back(newLine);
// 	  newLine[0] = mesher.mesh().triangles[4*i][1];
// 	  newLine[1] = mesher.mesh().triangles[4*i][2];
// 	  geometry.lines().push_back(newLine);
// 	  newLine[0] = mesher.mesh().triangles[4*i][2];
// 	  newLine[1] = mesher.mesh().triangles[4*i][0];
// 	  geometry.lines().push_back(newLine);
// 	  newLine[0] = mesher.mesh().triangles[4*i+1][2];
// 	  newLine[1] = mesher.mesh().triangles[4*i][0];
// 	  geometry.lines().push_back(newLine);
// 	  newLine[0] = mesher.mesh().triangles[4*i+1][2];
// 	  newLine[1] = mesher.mesh().triangles[4*i][1];
// 	  geometry.lines().push_back(newLine);
// 	  newLine[0] = mesher.mesh().triangles[4*i+1][2];
// 	  newLine[1] = mesher.mesh().triangles[4*i][2];
// 	  geometry.lines().push_back(newLine);
// 	}
// 	outputMessage = "Total Vertices: " + QString::number(mesher.mesh().numverts) + "\nTotal Tetrahedra: " + QString::number(mesher.mesh().numtris/4);
// 	#ifdef USING_SWEETMESH
// 	std::list<double> aspectRatios;
// 	sweetMesh::sweetMeshVertex v0, v1, v2, v3;
// 	double lowest = 20000000.0;
// 	double highest = -20.0;
// 	double average = 0.0;
// 	double ratio = 0;
// 	unsigned int bin1, bin2, bin3, bin4, bin5, bin6, bin7, bin8, bin9, bin10, bin11, bin12;
// 	bin1 = bin2 = bin3 = bin4 = bin5 = bin6 = bin7 = bin8 = bin9 = bin10 = bin11 = bin12 = 0;
// 	
// 	for(unsigned int i=0; i<mesher.mesh().numtris/4; i++){
// 	  v0.set(geometry.points()[mesher.mesh().triangles[4*i][0]][0], geometry.points()[mesher.mesh().triangles[4*i][0]][1], geometry.points()[mesher.mesh().triangles[4*i][0]][2]);
// 	  v1.set(geometry.points()[mesher.mesh().triangles[4*i][1]][0], geometry.points()[mesher.mesh().triangles[4*i][1]][1], geometry.points()[mesher.mesh().triangles[4*i][1]][2]);
// 	  v2.set(geometry.points()[mesher.mesh().triangles[4*i][2]][0], geometry.points()[mesher.mesh().triangles[4*i][2]][1], geometry.points()[mesher.mesh().triangles[4*i][2]][2]);
// 	  v3.set(geometry.points()[mesher.mesh().triangles[4*i+1][2]][0], geometry.points()[mesher.mesh().triangles[4*i+1][2]][1], geometry.points()[mesher.mesh().triangles[4*i+1][2]][2]);
// 	   
// 	  sweetMesh::tetrahedron tet(v0, v1, v2, v3);
// 	  ratio = tet.aspectRatio();
// 	  average += ratio;
// 	  lowest = std::min(lowest, ratio);
// 	  highest = std::max(highest, ratio);
// 	  aspectRatios.push_back(tet.aspectRatio());
// 	}
// 	average /= (double)aspectRatios.size();
// 	
// 	aspectRatios.sort();
// 	for(std::list<double>::iterator itr=aspectRatios.begin(); itr!=aspectRatios.end(); itr++){
// 	  if(*itr < 1.5){ bin1++; }
// 	  else{
// 	    if(*itr < 2.0){ bin2++; }
// 	    else{
// 	      if(*itr < 2.5){ bin3++; }
// 	      else{
// 		if(*itr < 3.0){ bin4++; }
// 		else{
// 		  if(*itr < 4.0){ bin5++; }
// 		  else{
// 		    if(*itr < 6.0){ bin6++; }
// 		    else{
// 		      if(*itr < 10.0){ bin7++; }
// 		      else{
// 			if(*itr < 15.0){ bin8++; }
// 			else{
// 			  if(*itr < 25){ bin9++; }
// 			  else{
// 			    if(*itr < 50){ bin10++; }
// 			    else{
// 			      if(*itr < 100){ bin11++; }
// 			      else{
// 				bin12++;
// 			      }
// 			    }
// 			  }
// 			}
// 		      }
// 		    }
// 		  }
// 		}
// 	      }
// 	    }
// 	  }
// 	}
// 	outputMessage = outputMessage + "\nAverage Aspect Ratio: " + QString::number(average) + "\nMinimal Aspect Ratio: " + QString::number(lowest) + "\nMaximal Aspect Ratio: " + QString::number(highest) + "\nAspect Ratio Histogram:"
// 	+ "\n      < 1.5  :  " + QString::number(bin1) + "\t|    6 - 10   :  " + QString::number(bin7)
// 	+ "\n  1.5 - 2.0  :  " + QString::number(bin2) + "\t|   10 - 15   :  " + QString::number(bin8)
// 	+ "\n  2.0 - 2.5  :  " + QString::number(bin3) + "\t|   15 - 25   :  " + QString::number(bin9)
// 	+ "\n  2.5 - 3.0  :  " + QString::number(bin4) + "\t|   25 - 50   :  " + QString::number(bin10)
// 	+ "\n  3.0 - 4.0  :  " + QString::number(bin5) + "\t|   50 - 100  :  " + QString::number(bin11)
// 	+ "\n  4.0 - 6.0  :  " + QString::number(bin6) + "\t|  100 -      :  " + QString::number(bin12);
// #endif
//       cvcapp.data(resultName, geometry);
//       cvcapp.listPropertyAppend("thumbnail.geometries", resultName);
//       cvcapp.listPropertyAppend("zoomed.geometries", resultName);
//       cvcapp.data(meshName, mesher);
//       }break;
//       default: break;
//     }
    QMessageBox::information(this, "Finished LBIE Meshing", outputMessage, QMessageBox::Ok);
  }
  


  // TODO: there is a bug in this function... it seems to crash sometimes
  //       I suspect this occurs when normal vectors are not available
  //  5-26-2011, above issue has been fixed?
  void CVCMainWindow::LBIE_qualityImprovement_Slot(){
    LBIEQualityImprovementDialog dialog;
    if(dialog.exec() == QDialog::Accepted) {
      dialog.RunLBIEQualityImprovement();
    }


    //arand, 5-26-2011, moved to LBIEQualityImprovementDialog.cpp
    //                  delete commented code from this file soon
    /*
    Ui::LBIEQualityImprovementDialogBase ui;
    QDialog dialog;
    ui.setupUi(&dialog);
   
    ui.m_Iterations->clear();
    ui.m_Iterations->insert("2");
    
    DataMap map = cvcapp.data();
    bool found = false;
    for (const auto& val : map) { 
      //std::cout << val.first << " " << val.second.type().name() << std::endl;
      std::string myname(val.first);
      std::string mytype(val.second.type().name());
      // only deal with files for now...
      if (cvcapp.isData<cvcraw_geometry::cvcgeom_t>(myname)) {
	ui.GeometryList->addItem(QString::fromStdString(myname));
	found = true;
      }
    }
    if (!found) {
      QMessageBox::information(this, tr("LBIE quality improvement."),
			   tr("No geometries loaded."), QMessageBox::Ok);
      return;
    }

    if(dialog.exec() != QDialog::Accepted) { return; }

    std::string geomSelected = ui.GeometryList->currentText().toStdString();
    cvcraw_geometry::cvcgeom_t geom = boost::any_cast<cvcraw_geometry::cvcgeom_t>(map[geomSelected]);
   
    string resultName = ui.OutputEdit->displayText().toStdString();
    if (resultName.empty()) {
      resultName = geomSelected + "_LBIE_imp";
    }
    int iterations = ui.m_Iterations->text().toUInt();

    // calculate normals if we don't have them...
    if(geom.normals().size() != geom.points().size()) {
      geom.calculate_surf_normals();
    }

    //convert to LBIE::geoframe
    //      FUTURE: rewrite LBIE and avoid conversions
    //      FUTURE: support different types of improvement
    LBIE::Octree oc;
    Geometry geo1 = Geometry::conv(geom);
    oc.setMeshType(geo1.m_GeoFrame->mesh_type);


    for(unsigned int i = 0; i < iterations; i++) {
      oc.quality_improve(*geo1.m_GeoFrame,LBIE::Mesher::GEO_FLOW);
    }   

    // TODO: fix the colors here so that the mesh preserves the original color (done)
    //       include the normals if available
    //       combine this conversion with the one in LBIE_slot()
    CVCGEOM_NAMESPACE::cvcgeom_t geometry;
    CVCGEOM_NAMESPACE::cvcgeom_t::color_t meshColor;
    meshColor[0] = 0.0; meshColor[1] = 1.0; meshColor[2] = 0.001;
    for(unsigned int i=0; i<geo1.m_GeoFrame->numverts; i++){
      CVCGEOM_NAMESPACE::cvcgeom_t::point_t newVertex;
      newVertex[0] = geo1.m_GeoFrame->verts[i][0];
      newVertex[1] = geo1.m_GeoFrame->verts[i][1];
      newVertex[2] = geo1.m_GeoFrame->verts[i][2];
      geometry.points().push_back(newVertex);

      // IMPORTANT: this only works if the number of vertices doesn't change
      //            and they are not reordered.
      //            If there is vertex removal or other complex operations
      //            the color information is going to be all messed up.
      if (geom.const_colors().size()>0) {
	geometry.colors().push_back(geom.const_colors()[i]);
      }
    }

    switch(geo1.m_GeoFrame->mesh_type){
      case LBIE::geoframe::SINGLE:{
	for(unsigned int i=0; i<geo1.m_GeoFrame->numtris; i++){
	  CVCGEOM_NAMESPACE::cvcgeom_t::triangle_t newTri;
	  newTri[0] = geo1.m_GeoFrame->triangles[i][2];
	  newTri[1] = geo1.m_GeoFrame->triangles[i][1];
	  newTri[2] = geo1.m_GeoFrame->triangles[i][0];
	  geometry.triangles().push_back(newTri);
	}	
      }break;
      
      case LBIE::geoframe::TETRA:{
	for(unsigned int i=0; i<geo1.m_GeoFrame->numtris/4; i++){
	  CVCGEOM_NAMESPACE::cvcgeom_t::line_t newLine;
	  newLine[0] = geo1.m_GeoFrame->triangles[4*i][0];
	  newLine[1] = geo1.m_GeoFrame->triangles[4*i][1];
	  geometry.lines().push_back(newLine);
	  newLine[0] = geo1.m_GeoFrame->triangles[4*i][1];
	  newLine[1] = geo1.m_GeoFrame->triangles[4*i][2];
	  geometry.lines().push_back(newLine);
	  newLine[0] = geo1.m_GeoFrame->triangles[4*i][2];
	  newLine[1] = geo1.m_GeoFrame->triangles[4*i][0];
	  geometry.lines().push_back(newLine);
	  newLine[0] = geo1.m_GeoFrame->triangles[4*i+1][2];
	  newLine[1] = geo1.m_GeoFrame->triangles[4*i][0];
	  geometry.lines().push_back(newLine);
	  newLine[0] = geo1.m_GeoFrame->triangles[4*i+1][2];
	  newLine[1] = geo1.m_GeoFrame->triangles[4*i][1];
	  geometry.lines().push_back(newLine);
	  newLine[0] = geo1.m_GeoFrame->triangles[4*i+1][2];
	  newLine[1] = geo1.m_GeoFrame->triangles[4*i][2];
	  geometry.lines().push_back(newLine);
	}
      }break;
      
      case LBIE::geoframe::QUAD:{
	for(unsigned int i=0; i<geo1.m_GeoFrame->numquads; i++){
	  CVCGEOM_NAMESPACE::cvcgeom_t::quad_t newQuad;
	  newQuad[0] = geo1.m_GeoFrame->quads[i][0];
	  newQuad[1] = geo1.m_GeoFrame->quads[i][1];
	  newQuad[2] = geo1.m_GeoFrame->quads[i][2];
	  newQuad[3] = geo1.m_GeoFrame->quads[i][3];
	  geometry.quads().push_back(newQuad);
	}
      }break;
      
      case LBIE::geoframe::HEXA:{
	for(unsigned int i=0; i<geo1.m_GeoFrame->numquads; i++){
	  CVCGEOM_NAMESPACE::cvcgeom_t::line_t newLine;
	  newLine[0] = geo1.m_GeoFrame->quads[i][0];
	  newLine[1] = geo1.m_GeoFrame->quads[i][1];
	  geometry.lines().push_back(newLine);
	  newLine[0] = geo1.m_GeoFrame->quads[i][2];
	  geometry.lines().push_back(newLine);
	  newLine[1] = geo1.m_GeoFrame->quads[i][3];
	  geometry.lines().push_back(newLine);
	  newLine[0] = geo1.m_GeoFrame->quads[i][0];
	  geometry.lines().push_back(newLine);
	}
      }break;
      
      case LBIE::geoframe::DOUBLE:{
	for(unsigned int i=0; i<geo1.m_GeoFrame->numtris; i++){
	  CVCGEOM_NAMESPACE::cvcgeom_t::triangle_t newTri;
	  newTri[0] = geo1.m_GeoFrame->triangles[i][0];
	  newTri[1] = geo1.m_GeoFrame->triangles[i][1];
	  newTri[2] = geo1.m_GeoFrame->triangles[i][2];
	  geometry.triangles().push_back(newTri);
	}
      }break;
      
      case LBIE::geoframe::TETRA2:{
	for(unsigned int i=0; i<geo1.m_GeoFrame->numtris/4; i++){
	  CVCGEOM_NAMESPACE::cvcgeom_t::line_t newLine;
	  newLine[0] = geo1.m_GeoFrame->triangles[4*i][0];
	  newLine[1] = geo1.m_GeoFrame->triangles[4*i][1];
	  geometry.lines().push_back(newLine);
	  newLine[0] = geo1.m_GeoFrame->triangles[4*i][1];
	  newLine[1] = geo1.m_GeoFrame->triangles[4*i][2];
	  geometry.lines().push_back(newLine);
	  newLine[0] = geo1.m_GeoFrame->triangles[4*i][2];
	  newLine[1] = geo1.m_GeoFrame->triangles[4*i][0];
	  geometry.lines().push_back(newLine);
	  newLine[0] = geo1.m_GeoFrame->triangles[4*i+1][2];
	  newLine[1] = geo1.m_GeoFrame->triangles[4*i][0];
	  geometry.lines().push_back(newLine);
	  newLine[0] = geo1.m_GeoFrame->triangles[4*i+1][2];
	  newLine[1] = geo1.m_GeoFrame->triangles[4*i][1];
	  geometry.lines().push_back(newLine);
	  newLine[0] = geo1.m_GeoFrame->triangles[4*i+1][2];
	  newLine[1] = geo1.m_GeoFrame->triangles[4*i][2];
	  geometry.lines().push_back(newLine);
	}
      }break;
      default: break;
    }

    cvcapp.data(resultName, geometry);


    cvcapp.listPropertyRemove("thumbnail.geometries", geomSelected);
    cvcapp.listPropertyRemove("zoomed.geometries", geomSelected);

    cvcapp.listPropertyAppend("thumbnail.geometries", resultName);
    cvcapp.listPropertyAppend("zoomed.geometries", resultName);

    QMessageBox::information(this, "Notice", "Finished LBIE Quality Improvement", QMessageBox::Ok);

    */

  }

  void CVCMainWindow::helpSlot() {
    QMessageBox::warning(this,"Warning",
			 "This feature has not been updated.");    
  }

  void CVCMainWindow::aboutVolRoverSlot() {
    QMessageBox::information(this, tr("VolRover 2.0a"),
			     tr("VolRover was developed at the Computational Visualization Center, University of Texas at Austin.  For more information, see our website http://cvcweb.ices.utexas.edu or contact Professor Chandrajit Bajaj, bajaj@cs.utexas.edu."), QMessageBox::Ok);   
  }

#ifdef USING_VOLUMEGRIDROVER
  void CVCMainWindow::volumeGridRoverSlot(){

     // make sure data is loaded
     if ( m_VolumeGridRover->hasVolume() )
	m_VolumeGridRover->show();
     else {
        QMessageBox::warning(this, tr("Segment Virus Map"),
	    tr("A volume must be loaded for this feature to work."), 0,0);
     }  
  }
#endif



#ifdef USING_MMHLS
void CVCMainWindow::generateRawVSlot() {
    if( !generateRawVUi )
  	   generateRawVUi = new Ui::generateRawVDialog();
	QDialog dialog;
	generateRawVUi->setupUi(&dialog);

	DataMap map = volrover_app_instance().data();


    for (const auto& val : map) { 
  //    std::cout << val.first << " " << val.second.type().name() << std::endl;
      std::string myname(val.first);
      std::string mytype(val.second.type().name());
	  //only deal with rawiv file now.
 	   if (mytype.compare("N9VolMagick14VolumeFileInfoE") == 0) {
	  	  generateRawVUi->m_VolumeList->addItem(QString::fromStdString(myname));
		  generateRawVUi->m_VolumeListForMesh->addItem(QString::fromStdString(myname));
		}
    }
	if(generateRawVUi->m_VolumeList->currentText().isEmpty())
	{
        QMessageBox::warning(this, tr("MMHLS"),
	    tr("A volume must be loaded for this feature to work."), 0,0);
		return;
	}

    if(dialog.exec() == QDialog::Accepted){
		if(generateRawVUi->tabWidget->currentIndex() == 0)
		{

			std::string volSelected =  generateRawVUi->m_VolumeList->currentText().toStdString();
			VolMagick::VolumeFileInfo vif = boost::any_cast<VolMagick::VolumeFileInfo>(map[volSelected]);

			int dimension = generateRawVUi->m_Dimension->text().toInt();
			float edgeLength = generateRawVUi->m_Edgelength->text().toFloat();
			string prefix = generateRawVUi->m_Prefix->text().toStdString();
			generateRawVFromVolume(vif.filename(), dimension, edgeLength, prefix);

		}
		else if(generateRawVUi->tabWidget->currentIndex() == 1)
		{
			std::string volSelected =  generateRawVUi->m_VolumeList->currentText().toStdString();
			VolMagick::VolumeFileInfo vif = boost::any_cast<VolMagick::VolumeFileInfo>(map[volSelected]);
			
			string meshPrefix = generateRawVUi->m_MeshPref->text().toStdString();
			int meshStart= generateRawVUi->m_MeshStart->text().toInt();
			int meshEnd = generateRawVUi->m_MeshEnd->text().toInt();
			int dimension = generateRawVUi->m_MeshDimension->text().toInt();
			float edgeLength = generateRawVUi->m_MeshEdgelength->text().toFloat();
			string prefix = generateRawVUi->m_MeshOutPrefix->text().toStdString();
			
			generateRawVFromMesh(vif.filename(), meshStart, meshEnd, meshPrefix,  dimension, edgeLength, prefix);

		}
	}
	

 }

void CVCMainWindow::openManifestFileSlot()
{
	m_manifestFile = QFileDialog::getOpenFileName(this, tr("Open File"), "./", tr("Manifest File (*.manifest)"));
}


void CVCMainWindow::generateMeshSlot(){
	if(!generateMeshUi)
	  	generateMeshUi = new Ui::generateMeshDialog();
	QDialog dialog;
	generateMeshUi-> setupUi(&dialog);
    
	connect((QPushButton*)(generateMeshUi->m_ManifestFile), SIGNAL(clicked()), this, SLOT(openManifestFileSlot()));

	if(dialog.exec() == QDialog::Accepted){
		 float isoratio = generateMeshUi->m_IsoRatio->text().toFloat();
		 float tolerance = generateMeshUi->m_Tolerance->text().toFloat();
		 float volthresh = generateMeshUi->m_VolThresh->text().toFloat();
		 int meshStart = generateMeshUi->m_MeshStart->text().toInt();
		 int meshEnd = generateMeshUi->m_MeshEnd->text().toInt();
		 string outpref = generateMeshUi->m_MeshOutPrefix->text().toStdString();		
		 string originalFile=m_manifestFile.toStdString();
		 generateMesh(originalFile, isoratio, tolerance, volthresh, meshStart, meshEnd, outpref);
		
		 QFileInfo fi = m_manifestFile;
		 string name = fi.fileName().toStdString();
		 originalFile.erase(originalFile.length()-name.length(), name.length());
		 for(int i=meshStart; i<=meshEnd; i++)
		   {
		     stringstream outGeom;
		     outGeom<<originalFile<<outpref<<"-"<<i<<".rawc";
		     if(cvcapp.readData(outGeom.str()))
		       {
		     cvcapp.properties("thumbnail.geometries",cvcapp.properties("thumbnail.geometries")+","+outGeom.str());
		     cvcapp.properties("zoomed.geometries",cvcapp.properties("zoomed.geometries")+","+outGeom.str());
		       }
		   }

	}

 }

#endif


#ifdef USING_MSLEVELSET
void CVCMainWindow::MSLevelSetSlot()
{
  using namespace MumfordShahLevelSet;
  if(!MSLevelSetDialogUi)
	  MSLevelSetDialogUi = new Ui::MSLevelSetDialog();
  QDialog dialog;
  MSLevelSetDialogUi->setupUi(&dialog);

  string VolFileName;

  DataMap map = volrover_app_instance().data();
  VolMagick::VolumeFileInfo vif;

  for (const auto& val : map) {
  	std::string myname(val.first);
	std::string mytype(val.second.type().name());
	if(mytype.compare("N9VolMagick14VolumeFileInfoE") == 0) {
	  vif = boost::any_cast<VolMagick::VolumeFileInfo> (map[myname]); 
	  VolFileName = myname;
	}
  }

  
  if(VolFileName.empty())
  {
    QMessageBox::warning(this, tr("MS Level Set"),
			 tr("A volume must be loaded for this feature to work."), 0,0);
    return;
  } 

  /*if(m_RGBARendering)
    {
      QMessageBox::warning(this, tr("MS Level Set"),
			   tr("This feature is not available for RGBA volumes."), 0,0);
      return;
    } */

  VolMagick::Volume vol;
  float *data;

  MSLevelSet mslsSolver(m_DoInitCUDA);
  if(m_DoInitCUDA) m_DoInitCUDA = false;  

  MSLevelSetParams MSLSParams;
  MSLSParams.lambda1 = 1.0; MSLSParams.lambda2 = 1.0;
  MSLSParams.mu = 32.5125; MSLSParams.nu = 0.0;
  MSLSParams.deltaT = 0.01; MSLSParams.epsilon = 1.0;
  MSLSParams.nIter = 15; MSLSParams.medianIter = 10;
  MSLSParams.DTWidth = 10; MSLSParams.medianTolerance = 0.0005f;
  MSLSParams.subvolDim = 128;
  MSLSParams.PDEBlockDim = 4; MSLSParams.avgBlockDim = 4;
  MSLSParams.SDTBlockDim = 4; MSLSParams.medianBlockDim = 4;
  MSLSParams.superEllipsoidPower = 8.0;
  MSLSParams.BBoxOffset = 5;
  MSLSParams.init_interface_method = 0;
  MSLSParams.volval_min = vif.min();
  MSLSParams.volval_max = vif.max();

  

  if(dialog.exec() != QDialog::Accepted)
   return; 
  /*
    Check if we are doing in place filtering of the actual volume data file
    instead of simply the current subvolume buffer
    TODO: need to implement out-of-core filtering using VolMagick
 */
 
  
  
      MSLSParams.lambda1 =  MSLevelSetDialogUi->m_Lambda1Edit->text().toDouble();
	  MSLSParams.lambda2 =  MSLevelSetDialogUi->m_Lambda2Edit->text().toDouble();
	  MSLSParams.mu =  MSLevelSetDialogUi->m_MuEdit->text().toDouble();
	  MSLSParams.nu =  MSLevelSetDialogUi->m_NuEdit->text().toDouble();
	  MSLSParams.epsilon =  MSLevelSetDialogUi->m_EpsilonEdit->text().toDouble();
	  MSLSParams.deltaT =  MSLevelSetDialogUi->m_DeltaTEdit->text().toDouble();
	  MSLSParams.medianTolerance =  MSLevelSetDialogUi->m_MedianTolEdit->text().toDouble();
	  MSLSParams.nIter =  MSLevelSetDialogUi->m_MaxSolverIterEdit->text().toInt();
	  MSLSParams.DTWidth =  MSLevelSetDialogUi->m_DTWidthEdit->text().toDouble();
	  MSLSParams.medianIter =  MSLevelSetDialogUi->m_MaxMedianIterEdit->text().toInt();
	  MSLSParams.subvolDim =  MSLevelSetDialogUi->m_SubvolDimEdit->text().toInt();
	  MSLSParams.SDTBlockDim =  MSLevelSetDialogUi->m_BlockDimEdit->text().toInt();
	  MSLSParams.superEllipsoidPower =  MSLevelSetDialogUi->m_EllipsoidPowerEdit->text().toDouble();
	  MSLSParams.init_interface_method =  MSLevelSetDialogUi->m_DTInitComboBox->currentIndex();
	  int blockDim =  MSLevelSetDialogUi->m_BlockDimEdit->text().toInt();

	  switch(MSLevelSetDialogUi->m_BlockDimComboBox->currentIndex())
	  {
	  	case 0:
		MSLSParams.SDTBlockDim = blockDim;
		break;
	  	case 1:
		MSLSParams.avgBlockDim = blockDim;
		break;
	  	case 2:
		MSLSParams.medianBlockDim = blockDim;
		break;
	  	case 3:
		MSLSParams.PDEBlockDim = blockDim;
		break;
	  }

	  MSLSParams.BBoxOffset =  MSLevelSetDialogUi->m_BBoxOffsetEdit->text().toInt();

	 if(!MSLevelSetDialogUi->m_Preview->isChecked())
     {
      if(QMessageBox::warning(this,
			      "MSLevelSet",
			      "Are you sure you want to do this?\n"
			      "This operation will change the current loaded volume file\n"
			      "and cannot be un-done!",
			      QMessageBox::Cancel | QMessageBox::Default | QMessageBox::Escape,
			      QMessageBox::Ok) != QMessageBox::Ok)
	    return;
    
      try
  	 {

	  VolMagick::readVolumeFile(cvcapp, vol, vif.filename());
	  vol.voxelType(VolMagick::Float); //forces float voxel type
	  data = reinterpret_cast<float*>(*vol);
	  if(mslsSolver.runSolver(data, vol.XDim(), vol.YDim(), vol.ZDim(), &MSLSParams))
	    {
		  QFileInfo fi = QString::fromStdString(vif.filename());
		  string name = fi.fileName().toStdString();
		  string originalName = vif.filename();
		  originalName.erase(originalName.length()-name.length(), name.length());
		  QDir tmpdir(QString::fromStdString(originalName + "tmp"));

	      if(!tmpdir.exists()) tmpdir.mkdir(tmpdir.absPath());
	      QFileInfo tmpfile(tmpdir,QFileInfo("MSLevelSet_outputPhi.rawiv").fileName());
	      qDebug("Writting volume %s to %s (%s)", "MSLevelSet_outputPhi.rawiv", tmpfile.absFilePath().toStdString().c_str(), QString::fromStdString(originalName).toStdString().c_str());
	      QString newfilename = QDir::convertSeparators(tmpfile.absFilePath());
	      if(tmpfile.exists()) QFile::remove(newfilename);	    
	        VolMagick::createVolumeFile(cvcapp, vol,
					 newfilename.toStdString());
			cvcapp.readData(newfilename.toStdString());
	    }
	  else
	    QMessageBox::warning(this,
				 "MS Level Set",
				 "Solver failed to execute!", 0, 0);
	}
      catch(std::exception& e)
	{
	  QMessageBox::critical(this,"Error: ",e.what());
	  return;
	}
	

      return;
    }

    bool hasVolume = false;
    VolMagick::Volume vol1;

	std::string volselect("thumbnail_volume"); //zoomed_volume");
	if(cvcapp.isData<VolMagick::Volume>(volselect))
    {
        vol1 = cvcapp.data<VolMagick::Volume>(volselect);
		hasVolume = true;
    }

  if(!hasVolume)
  {
  	QMessageBox::information(this, tr("MSLevelSet"), 
		tr("No subvolume in zoomed in window!"), QMessageBox::Ok);
  	return;
   }

  unsigned int xdim, ydim, zdim, len;
  
  xdim = vol1.XDim();
  ydim = vol1.YDim();
  zdim = vol1.ZDim();
  len = xdim*ydim*zdim;
  vol1.voxelType(VolMagick::Float);
  float *org_dat = reinterpret_cast<float*>(*vol1);

  VolMagick::Volume volout(VolMagick::Dimension(xdim,ydim,zdim),VolMagick::Float);

  memcpy(*volout,org_dat,len*sizeof(float));

  data = reinterpret_cast<float*>(*volout);

  if(!mslsSolver.runSolver(data, volout.XDim(), volout.YDim(), volout.ZDim(), &MSLSParams))
  {
    QMessageBox::critical(this,"Error","MSLevelSet Solver failed to execute!");
    return;
  }  
 
  cvcapp.data("MSLevelSet_seg.rawiv", volout);

}


#ifdef USING_HOSEGMENTATION
void CVCMainWindow::HOSegmentationSlot()
{
	using namespace HigherOrderSegmentation;
	if(!MSLevelSetDialogUi)
	  	MSLevelSetDialogUi = new Ui::MSLevelSetDialog();
  	QDialog dialog;
  	MSLevelSetDialogUi->setupUi(&dialog);
  	m_DoInitCUDA = true;
  	string VolFileName;

  	DataMap map = volrover_app_instance().data();
	VolMagick::VolumeFileInfo vif;

  	for (const auto& val : map) {
  		std::string myname(val.first);
		std::string mytype(val.second.type().name());
		if(mytype.compare("N9VolMagick14VolumeFileInfoE") == 0) {
	    	vif = boost::any_cast<VolMagick::VolumeFileInfo> (map[myname]); 
	  	    VolFileName = myname;
	    }
  	}

  
  if(VolFileName.empty())
  {
    QMessageBox::warning(this, tr("Higher Order Level Set Segmentation"),
			 tr("A volume must be loaded for this feature to work."), 0,0);
    return;
  } 
     VolMagick::Volume vol;
  float *data;


  HOSegmentation hosegSolver(m_DoInitCUDA);
  if(m_DoInitCUDA) m_DoInitCUDA = false;  

  MSLevelSetParams MSLSParams;
  MSLSParams.lambda1 = 1.0; MSLSParams.lambda2 = 1.0;
  MSLSParams.mu = 32.5125; MSLSParams.nu = 0.0;
  MSLSParams.deltaT = 0.01; MSLSParams.epsilon = 1.0;
  MSLSParams.nIter = 15; MSLSParams.medianIter = 10;
  MSLSParams.DTWidth = 10; MSLSParams.medianTolerance = 0.0005f;
  MSLSParams.subvolDim = 128;
  MSLSParams.PDEBlockDim = 4; MSLSParams.avgBlockDim = 4;
  MSLSParams.SDTBlockDim = 4; MSLSParams.medianBlockDim = 4;
  MSLSParams.superEllipsoidPower = 8.0;
  MSLSParams.BBoxOffset = 5;
  MSLSParams.init_interface_method = 0;
  MSLSParams.volval_min = vif.min();
  MSLSParams.volval_max = vif.max();

  if(dialog.exec() != QDialog::Accepted)
     return; 
  MSLSParams.lambda1 =  MSLevelSetDialogUi->m_Lambda1Edit->text().toDouble();
  MSLSParams.lambda2 =  MSLevelSetDialogUi->m_Lambda2Edit->text().toDouble();
  MSLSParams.mu =  MSLevelSetDialogUi->m_MuEdit->text().toDouble();
  MSLSParams.nu =  MSLevelSetDialogUi->m_NuEdit->text().toDouble();
  MSLSParams.epsilon =  MSLevelSetDialogUi->m_EpsilonEdit->text().toDouble();
  MSLSParams.deltaT =  MSLevelSetDialogUi->m_DeltaTEdit->text().toDouble();
  MSLSParams.medianTolerance =  MSLevelSetDialogUi->m_MedianTolEdit->text().toDouble();
  MSLSParams.nIter =  MSLevelSetDialogUi->m_MaxSolverIterEdit->text().toInt();
  MSLSParams.DTWidth =  MSLevelSetDialogUi->m_DTWidthEdit->text().toDouble();
  MSLSParams.medianIter =  MSLevelSetDialogUi->m_MaxMedianIterEdit->text().toInt();
  MSLSParams.subvolDim =  MSLevelSetDialogUi->m_SubvolDimEdit->text().toInt();
  MSLSParams.SDTBlockDim =  MSLevelSetDialogUi->m_BlockDimEdit->text().toInt();
  MSLSParams.superEllipsoidPower =  MSLevelSetDialogUi->m_EllipsoidPowerEdit->text().toDouble();
  MSLSParams.init_interface_method =  MSLevelSetDialogUi->m_DTInitComboBox->currentIndex();
  int blockDim =  MSLevelSetDialogUi->m_BlockDimEdit->text().toInt();

  switch(MSLevelSetDialogUi->m_BlockDimComboBox->currentIndex())
  {
  	case 0:
	MSLSParams.SDTBlockDim = blockDim;
	break;
  	case 1:
	MSLSParams.avgBlockDim = blockDim;
	break;
  	case 2:
	MSLSParams.medianBlockDim = blockDim;
	break;
  	case 3:
	MSLSParams.PDEBlockDim = blockDim;
	break;
  }

  MSLSParams.BBoxOffset =  MSLevelSetDialogUi->m_BBoxOffsetEdit->text().toInt();

  if(!MSLevelSetDialogUi->m_Preview->isChecked())
    {
      if(QMessageBox::warning(this,
			      "Higher Order LevelSet Segmenation",
			      "Are you sure you want to do this?\n"
			      "This operation will change the current loaded volume file\n"
			      "and cannot be un-done!",
			      QMessageBox::Cancel | QMessageBox::Default | QMessageBox::Escape,
			      QMessageBox::Ok) != QMessageBox::Ok)
	    return;
    
 		try
  	 	{
	  	   VolMagick::readVolumeFile(cvcapp, vol, vif.filename());
		   vol.voxelType(VolMagick::Float); //forces float voxel type
	 	   data = reinterpret_cast<float*>(*vol);
		   if(hosegSolver.runSolver(data, vol.XDim(), vol.YDim(), vol.ZDim(), &MSLSParams))
	    	{
			  QFileInfo fi = QString::fromStdString(vif.filename());
			  string name = fi.fileName().toStdString();
			  string originalName = vif.filename();
			  originalName.erase(originalName.length()-name.length(), name.length());
			  QDir tmpdir(QString::fromStdString(originalName + "tmp"));

		      if(!tmpdir.exists()) tmpdir.mkdir(tmpdir.absPath());
	    	  QFileInfo tmpfile(tmpdir,QFileInfo("HOLevelSet_outputPhi.rawiv").fileName());
		      qDebug("Writting volume %s to %s (%s)", "HOLevelSet_outputPhi.rawiv", tmpfile.absFilePath().toStdString().c_str(), QString::fromStdString(originalName).toStdString().c_str());
		      QString newfilename = QDir::convertSeparators(tmpfile.absFilePath());
	    	  if(tmpfile.exists()) QFile::remove(newfilename);	    
	        	VolMagick::createVolumeFile(cvcapp, vol,
					 newfilename.toStdString());
				cvcapp.readData(newfilename.toStdString());
	    	}
	  	else
	    	QMessageBox::warning(this,
					 "Higher Order Level Set Segmentation",
					 "Solver failed to execute!", 0, 0);
	}
      catch(std::exception& e)
	{
	  QMessageBox::critical(this,"Error: ",e.what());
	  return;
	}
	

      return;
    }

    bool hasVolume = false;
    VolMagick::Volume vol1;

	std::string volselect("thumbnail_volume"); //zoomed_volume");
	if(cvcapp.isData<VolMagick::Volume>(volselect))
    {
        vol1 = cvcapp.data<VolMagick::Volume>(volselect);
		hasVolume = true;
    }

 	 if(!hasVolume)
     {
  		QMessageBox::information(this, tr("HOLevelSet Segmentation"), 
		 tr("No subvolume in zoomed in window!"), QMessageBox::Ok);
	  	return;
     }

  	unsigned int xdim, ydim, zdim, len;
  
	xdim = vol1.XDim();
	ydim = vol1.YDim();
    zdim = vol1.ZDim();
    len = xdim*ydim*zdim;
    vol1.voxelType(VolMagick::Float);
    float *org_dat = reinterpret_cast<float*>(*vol1);

    VolMagick::Volume volout(VolMagick::Dimension(xdim,ydim,zdim),VolMagick::Float);

    memcpy(*volout,org_dat,len*sizeof(float));

    data = reinterpret_cast<float*>(*volout);

    if(!hosegSolver.runSolver(data, volout.XDim(), volout.YDim(), volout.ZDim(), &MSLSParams))
    {
      QMessageBox::critical(this,"Error","Higher Order Segmentation Solver failed to execute!");
      return;
    }  
 
      cvcapp.data("HOLevelSet_seg.rawiv", volout); 
}


#ifdef USING_MPSEGMENTATION

void CVCMainWindow::on_UserSegReadButton_clickedSlot()
{
    QString s =  QFileDialog::getOpenFileName(this, tr("open File"), "./", tr("All (*.*)"));
	if(s == QString())
		strcpy(MPLSParams->userSegFileName, "");
	else
	{
		strncpy(MPLSParams->userSegFileName, s.toStdString().c_str(), s.length());
		MPLSParams->userSegFileName[s.length()]='\0';
	}
	MPLevelSetDialogUi->m_userSegFileEdit->setText(QString("%1").arg(MPLSParams->userSegFileName));

	if(strlen(MPLSParams->userSegFileName)!=0){
	    FILE *fid = fopen(MPLSParams->userSegFileName, "rb");
	    if(!fid) {
    	  debugPrint("Cannot open user segmentation file");
      	strcpy(MPLSParams->userSegFileName, "");
      	return;
    }
    unsigned int m_nClasses;
    fread(&m_nClasses, sizeof(unsigned int), 1, fid);//Read 32 bit integer
    MPLSParams->nImplicits = (unsigned int)ceil(log2(float(m_nClasses))); 
    fclose(fid);
    MPLevelSetDialogUi->m_nImplicitsSpin->setValue(MPLSParams->nImplicits);
  }
}

void CVCMainWindow::MPSegmentationSlot()
{
	using namespace MultiphaseSegmentation;
	if(!MPLevelSetDialogUi)
		MPLevelSetDialogUi = new Ui::MPLevelSetDialog();
	QDialog dialog;
	MPLevelSetDialogUi->setupUi(&dialog);
  	m_DoInitCUDA = true;
  	string VolFileName;

	connect((QPushButton*)(MPLevelSetDialogUi->m_UserSegReadButton), SIGNAL(clicked()), this, SLOT(on_UserSegReadButton_clickedSlot()));

  	DataMap map = volrover_app_instance().data();
	VolMagick::VolumeFileInfo vif;

  	for (const auto& val : map) {
  		std::string myname(val.first);
		std::string mytype(val.second.type().name());
		if(mytype.compare("N9VolMagick14VolumeFileInfoE") == 0) {
	    	vif = boost::any_cast<VolMagick::VolumeFileInfo> (map[myname]); 
	  	    VolFileName = myname;
	    }
  	}

  if(VolFileName.empty())
  {
    QMessageBox::warning(this, tr("Multiple Phase Segmentation"),
			 tr("A volume must be loaded for this feature to work."), 0,0);
    return;
  } 
  VolMagick::Volume vol;
  float *data;
  float **l_PHI;

  MPSegmentation mpsegSolver(m_DoInitCUDA);
  if(m_DoInitCUDA) m_DoInitCUDA = false;


  //Disable preview mode. It is difficult to get the vector stack in the preview mode.
  MPLevelSetDialogUi->m_Preview->setEnabled(FALSE);
  MPLevelSetDialogUi->m_Preview->setChecked(FALSE);

  if(!MPLSParams)
  	MPLSParams = new MPLevelSetParams();

  MPLSParams->lambda1 = 1.0; MPLSParams->lambda2 = 1.0;
  MPLSParams->mu = 32.5125; MPLSParams->nu = 0.0;
  MPLSParams->deltaT = 0.01; MPLSParams->epsilon = 1.0;
  MPLSParams->nIter = 15; MPLSParams->medianIter = 10;
  MPLSParams->DTWidth = 10; MPLSParams->medianTolerance = 0.0005f;
  MPLSParams->subvolDim = 128;// 256; //128;
  MPLSParams->PDEBlockDim = 4; MPLSParams->avgBlockDim = 4;
  MPLSParams->SDTBlockDim = 4; MPLSParams->medianBlockDim = 4;
  MPLSParams->superEllipsoidPower = 8.0;
  MPLSParams->BBoxOffset = 5;
  MPLSParams->init_interface_method = 0;
  MPLSParams->nImplicits = 3;
  MPLevelSetDialogUi->m_nImplicitsSpin->setValue(3);
  strcpy(MPLSParams->userSegFileName, "");
  MPLSParams->volval_min = vif.min();
  MPLSParams->volval_max = vif.max();
  MPLSParams->subvolDimSDT = 128;
  MPLSParams->subvolDimAvg = 128;
  MPLSParams->subvolDimPDE = 128;

  if(dialog.exec() != QDialog::Accepted)
     return; 
  
  MPLSParams->lambda1 =  MPLevelSetDialogUi->m_Lambda1Edit->text().toDouble();
  MPLSParams->lambda2 =  MPLevelSetDialogUi->m_Lambda2Edit->text().toDouble();
  MPLSParams->mu =  MPLevelSetDialogUi->m_MuEdit->text().toDouble();
  MPLSParams->nu =  MPLevelSetDialogUi->m_NuEdit->text().toDouble();
  MPLSParams->epsilon =  MPLevelSetDialogUi->m_EpsilonEdit->text().toDouble();
  MPLSParams->deltaT =  MPLevelSetDialogUi->m_DeltaTEdit->text().toDouble();
  MPLSParams->medianTolerance =  MPLevelSetDialogUi->m_MedianTolEdit->text().toDouble();
  MPLSParams->nIter =  MPLevelSetDialogUi->m_MaxSolverIterEdit->text().toInt();
  MPLSParams->DTWidth =  MPLevelSetDialogUi->m_DTWidthEdit->text().toDouble();
  MPLSParams->medianIter =  MPLevelSetDialogUi->m_MaxMedianIterEdit->text().toInt();
  MPLSParams->superEllipsoidPower =  MPLevelSetDialogUi->m_EllipsoidPowerEdit->text().toDouble();
  MPLSParams->init_interface_method =  MPLevelSetDialogUi->m_DTInitComboBox->currentIndex();
  MPLSParams->BBoxOffset =  MPLevelSetDialogUi->m_BBoxOffsetEdit->text().toInt();
  MPLSParams->nImplicits =  MPLevelSetDialogUi->m_nImplicitsSpin->value();
  int blockDim =  MPLevelSetDialogUi->m_BlockDimEdit->text().toInt();

  switch(MPLevelSetDialogUi->m_BlockDimComboBox->currentIndex())
  {
  	case 0:
	MPLSParams->SDTBlockDim = blockDim;
	break;
  	case 1:
	MPLSParams->avgBlockDim = blockDim;
	break;
  	case 2:
	MPLSParams->medianBlockDim = blockDim;
	break;
  	case 3:
	MPLSParams->PDEBlockDim = blockDim;
	break;
  }

  int subvolDim = MPLevelSetDialogUi->m_SubvolDimEdit->text().toInt();
  switch(MPLevelSetDialogUi->m_SubvolDimComboBox->currentIndex())
  {
  	case 0:
	MPLSParams->subvolDimSDT = subvolDim;
	break;
  	case 1:
	MPLSParams->subvolDimAvg = subvolDim;
	break;
  	case 2:
	MPLSParams->subvolDimPDE = subvolDim;
	break;
  }
  

  if(!MPLevelSetDialogUi->m_Preview->isChecked())
  {
  	 if(QMessageBox::warning(this,
			      "Multi-phase Segmentation",
			      "Are you sure you want to do this?\n"
			      "This operation will change the current loaded volume file\n"
			      "and cannot be un-done!",
			      QMessageBox::Cancel | QMessageBox::Default | QMessageBox::Escape,
			      QMessageBox::Ok) != QMessageBox::Ok)
	return;

    try 
    {
	  VolMagick::readVolumeFile(cvcapp, vol, vif.filename());
	  vol.voxelType(VolMagick::Float); //forces float voxel type
	  data = reinterpret_cast<float*>(*vol);
	  int nvols = MPLSParams->nImplicits;
	  l_PHI = new float*[nvols];
	  for(int i=0; i< nvols; i++)
	     l_PHI[i] = new float[vol.XDim()*vol.YDim()*vol.ZDim()];
	  if(mpsegSolver.runSolver(data, l_PHI, vol.XDim(), vol.YDim(), vol.ZDim(), MPLSParams))
	  {
		  QFileInfo fi = QString::fromStdString(vif.filename());
		  string name = fi.fileName().toStdString();
		  string originalName = vif.filename();
		  originalName.erase(originalName.length()-name.length(), name.length());
		  QDir tmpdir(QString::fromStdString(originalName + "tmp"));

	      if(!tmpdir.exists()) tmpdir.mkdir(tmpdir.absPath());

		  VolMagick::Dimension d;
	      d.xdim = vol.XDim();
	      d.ydim = vol.YDim();
	      d.zdim = vol.ZDim();
	 	
		  //Dump all the implicit functions
	      QFileInfo tmpfile(tmpdir,QFileInfo("MPLevelSet_outputPhi.rawv").fileName());
	      qDebug("Writting volume %s to %s (%s)", "MPLevelSet_outputPhi.rawv", tmpfile.absFilePath().toStdString().c_str(), QString::fromStdString(originalName).toStdString().c_str());
	      QString newfilename = QDir::convertSeparators(tmpfile.absFilePath());
	      if(tmpfile.exists()) QFile::remove(newfilename);	    
		  std::vector<VolMagick::Volume> implicit_vols(nvols);
		  char _str[256];

		for(int m = 0; m<nvols; m++)
		  {
			  sprintf(_str, "Implicit %d",m);
	        implicit_vols[m].desc(_str);
	        implicit_vols[m].voxelType(VolMagick::Float);
	        implicit_vols[m].dimension(d);  
	        implicit_vols[m].boundingBox(vol.boundingBox());
 	 	   }
		  
		   for(unsigned int m = 0; m < nvols; m++)
		   	 for(unsigned int k = 0; k < vol.ZDim(); k++)
		   	   for(unsigned int j = 0; j < vol.YDim(); j++)
		   	     for(unsigned int i = 0; i < vol.XDim(); i++)
				    implicit_vols[m](i,j,k, l_PHI[m][i+ vol.XDim()*(j+vol.YDim()*k)]);

	        VolMagick::writeVolumeFile(cvcapp, implicit_vols, newfilename.toStdString());

			implicit_vols.clear();
			for(int j=0; j<nvols; j++)
			  delete l_PHI[j];
			delete []l_PHI;

			cvcapp.readData(newfilename.toStdString());
	   }
	   else
	   QMessageBox::warning(this,
				 "Multi-phase Segmentation",
				 "Solver failed to execute!", 0, 0);
	}
    catch(std::exception& e)
    {
        QMessageBox::critical(this,"Error: ",e.what());
        return;
    }

    return;
  }

   bool hasVolume = false;
   VolMagick::Volume vol1;

	std::string volselect("thumbnail_volume"); //zoomed_volume");
	if(cvcapp.isData<VolMagick::Volume>(volselect))
    {
        vol1 = cvcapp.data<VolMagick::Volume>(volselect);
		hasVolume = true;
    }

 	 if(!hasVolume)
     {
  		QMessageBox::information(this, tr("MP Segmentation"), 
		 tr("No subvolume in zoomed in window!"), QMessageBox::Ok);
	  	return;
     }

  	unsigned int xdim, ydim, zdim, len;
  

	xdim = vol1.XDim();
	ydim = vol1.YDim();
    zdim = vol1.ZDim();
    len = xdim*ydim*zdim;
    vol1.voxelType(VolMagick::Float);
    float *org_dat = reinterpret_cast<float*>(*vol1);

	int nvols = MPLSParams->nImplicits;
	l_PHI = new float*[nvols];


    VolMagick::Volume volout(VolMagick::Dimension(xdim,ydim,zdim),VolMagick::Float);

    memcpy(*volout,org_dat,len*sizeof(float));

    data = reinterpret_cast<float*>(*volout);

    for(int i=0; i<nvols; i++)
	   l_PHI[i] = new float[xdim*ydim*zdim];
	
    if(!mpsegSolver.runSolver(data,l_PHI, volout.XDim(), volout.YDim(), volout.ZDim(), MPLSParams))
    {
      QMessageBox::critical(this,"Error","Multiple Phase Segmentation Solver failed to execute!");
      return;
    }  
   //TODO: Do something with the returned buffer l_PHI. This is currently unused.
    for(int j=0; j<nvols; j++)
       delete l_PHI[j];
    delete []l_PHI;

//      cvcapp.data("MP_seg.rawv", volout); 


}
#endif
#endif
#endif







void CVCMainWindow::unimplementedSlot() {
    QMessageBox::warning(this,"Warning",
			 "This feature has not been updated.");    
  }

   LocalSegThread::LocalSegThread(std::string filename, Ui::SegmentVirusMapDialog *dialog, CVCMainWindow *nvmw, unsigned int stackSize)
    : QThread(nvmw), m_CVCMainWindow(nvmw)
   {
	 type = dialog->m_TabSegmentationType->currentIndex();

	 m_Params[0] = filename;
	 
	 switch(dialog->m_TabSegmentationType->currentIndex())
	 {
		 case 0: // Capsid
			 
			 switch(dialog->m_TabCapsidLayerType->currentIndex())
			 {
				 case 0:
					 m_Params[1] = int(0);
					 m_Params[2] = dialog->m_TLowEditType0->text().toDouble();
					 m_Params[3] = dialog->m_X0EditType0->text().toInt();
					 m_Params[4] = dialog->m_Y0EditType0->text().toInt();
					 m_Params[5] = dialog->m_Z0EditType0->text().toInt();
					 // necessary to make VC++6.0 happy
					 m_Params[6] = XmlRpc::XmlRpcValue(dialog->m_RunDiffusionType0->isChecked());
					 break;
				 case 1:
					 m_Params[1] = int(1);
					 m_Params[2] = dialog->m_TLowEditType1->text().toDouble();
					 m_Params[3] = dialog->m_X0EditType1->text().toInt();
					 m_Params[4] = dialog->m_Y0EditType1->text().toInt();
					 m_Params[5] = dialog->m_Z0EditType1->text().toInt();
					 m_Params[6] = dialog->m_X1EditType1->text().toInt();
					 m_Params[7] = dialog->m_Y1EditType1->text().toInt();
					 m_Params[8] = dialog->m_Z1EditType1->text().toInt();
					 m_Params[9] = XmlRpc::XmlRpcValue(dialog->m_RunDiffusionType1->isChecked());
					 break;
				 case 2:
					 m_Params[1] = int(2);
					 m_Params[2] = dialog->m_TLowEditType2->text().toDouble();
					 m_Params[3] = dialog->m_SmallRadiusEditType2->text().toDouble();
					 m_Params[4] = dialog->m_LargeRadiusEditType2->text().toDouble();
					 break;
				 case 3:
					 m_Params[1] = int(3);
					 m_Params[2] = dialog->m_TLowEditType3->text().toDouble();
					 m_Params[3] = dialog->m_3FoldEditType3->text().toInt();
					 m_Params[4] = dialog->m_5FoldEditType3->text().toInt();
					 m_Params[5] = dialog->m_6FoldEditType3->text().toInt();
					 m_Params[6] = dialog->m_SmallRadiusEditType3->text().toDouble();
					 m_Params[7] = dialog->m_LargeRadiusEditType3->text().toDouble();
					 break;
			 }
			 
			 break;
		 case 1: // Monomer
			 
			 m_Params[1] = dialog->m_FoldNumEdit->text().toInt();
			 
			 break;
		 case 2: // Subunit
			 
			 m_Params[1] = dialog->m_HNumEdit->text().toInt();
			 m_Params[2] = dialog->m_KNumEdit->text().toInt();
			 m_Params[3] = dialog->m_3FoldEdit->text().toInt();
			 m_Params[4] = dialog->m_5FoldEdit->text().toInt();
			 m_Params[5] = dialog->m_6FoldEdit->text().toInt();
			 m_Params[6] = dialog->m_InitRadiusEdit->text().toInt();
			 
			 break;
	 }
   }

   void LocalSegThread::run()
   {

#ifdef USING_SEGMENTATION
	using namespace XmlRpc; /* we're only going to use XmlRpcValue so we can make the code more uniform... */
	printf("NewVolumeMainWindow::LocalSegThread::run(): Local segmentation thread started.\n");
	XmlRpc::XmlRpcValue result;
	
	switch(type)
	{
		case 0: /* Capsid */
			
			if(!virusSegCapsid(m_Params,result))
				QApplication::postEvent(m_CVCMainWindow,new SegmentationFailedEvent("Error segmenting virus capsid!"));
			else
				QApplication::postEvent(m_CVCMainWindow,new SegmentationFinishedEvent("Local segmentation complete."));
			
			break;
		case 1: /* Monomer */
			
			if(!virusSegMonomer(m_Params,result))
				QApplication::postEvent(m_CVCMainWindow,new SegmentationFailedEvent("Error segmenting virus monomer!"));
			else
				QApplication::postEvent(m_CVCMainWindow,new SegmentationFinishedEvent("Local segmentation complete."));
			
			break;
		case 2: /* Subunit */
		
			if(!virusSegSubunit(m_Params,result))
				QApplication::postEvent(m_CVCMainWindow,new SegmentationFailedEvent("Error segmenting virus subunit!"));
			else
				QApplication::postEvent(m_CVCMainWindow,new SegmentationFinishedEvent("Local segmentation complete."));

			break;
	}

	printf("NewVolumeMainWindow::LocalSegThread::run(): Local segmentation thread finished.\n");
#endif

   }

   RemoteSegThread::RemoteSegThread(Ui::SegmentVirusMapDialog *dialog, CVCMainWindow *nvmw, unsigned int stackSize)
    : QThread(nvmw), m_CVCMainWindow(nvmw),m_XmlRpcClient(dialog->m_RemoteSegmentationHostname->text().toUtf8().constData(), dialog->m_RemoteSegmentationPort->text().toInt())
   {

#ifdef USING_SEGMENTATION
	 type = dialog->m_TabSegmentationType->currentIndex();
	 
	 m_Params[0] = std::string(dialog->m_RemoteSegmentationFilename->text().toStdString().c_str());
	 
	 switch(dialog->m_TabSegmentationType->currentIndex())
	 {
		 case 0: // Capsid
			 
			 switch(dialog->m_TabCapsidLayerType->currentIndex())
			 {
				 case 0:
					 m_Params[1] = int(0);
					 m_Params[2] = dialog->m_TLowEditType0->text().toDouble();
					 m_Params[3] = dialog->m_X0EditType0->text().toInt();
					 m_Params[4] = dialog->m_Y0EditType0->text().toInt();
					 m_Params[5] = dialog->m_Z0EditType0->text().toInt();
					 m_Params[6] = XmlRpc::XmlRpcValue(dialog->m_RunDiffusionType0->isChecked());
					 break;
				 case 1:
					 m_Params[1] = int(1);
					 m_Params[2] = dialog->m_TLowEditType1->text().toDouble();
					 m_Params[3] = dialog->m_X0EditType1->text().toInt();
					 m_Params[4] = dialog->m_Y0EditType1->text().toInt();
					 m_Params[5] = dialog->m_Z0EditType1->text().toInt();
					 m_Params[6] = dialog->m_X1EditType1->text().toInt();
					 m_Params[7] = dialog->m_Y1EditType1->text().toInt();
					 m_Params[8] = dialog->m_Z1EditType1->text().toInt();
					 m_Params[9] = XmlRpc::XmlRpcValue(dialog->m_RunDiffusionType1->isChecked());
					 break;
				 case 2:
					 m_Params[1] = int(2);
					 m_Params[2] = dialog->m_TLowEditType2->text().toDouble();
					 m_Params[3] = dialog->m_SmallRadiusEditType2->text().toDouble();
					 m_Params[4] = dialog->m_LargeRadiusEditType2->text().toDouble();
					 break;
				 case 3:
					 m_Params[1] = int(3);
					 m_Params[2] = dialog->m_TLowEditType3->text().toDouble();
					 m_Params[3] = dialog->m_3FoldEditType3->text().toInt();
					 m_Params[4] = dialog->m_5FoldEditType3->text().toInt();
					 m_Params[5] = dialog->m_6FoldEditType3->text().toInt();
					 m_Params[6] = dialog->m_SmallRadiusEditType3->text().toDouble();
					 m_Params[7] = dialog->m_LargeRadiusEditType3->text().toDouble();
					 break;
			 }
			 
			 break;
		 case 1: // Monomer
			 
			 m_Params[1] = dialog->m_FoldNumEdit->text().toInt();
			 
			 break;
		 case 2: // Subunit
			 
			 m_Params[1] = dialog->m_HNumEdit->text().toInt();
			 m_Params[2] = dialog->m_KNumEdit->text().toInt();
			 m_Params[3] = dialog->m_3FoldEdit->text().toInt();
			 m_Params[4] = dialog->m_5FoldEdit->text().toInt();
			 m_Params[5] = dialog->m_6FoldEdit->text().toInt();
			 m_Params[6] = dialog->m_InitRadiusEdit->text().toInt();
			 
			 break;
	 }
#endif
   }

   void RemoteSegThread::run()
   {

#ifdef USING_SEGMENTATION
	using namespace XmlRpc; /* we're only going to use XmlRpcValue so we can make the code more uniform... */
	printf("NewVolumeMainWindow::RemoteSegThread::run(): Remote segmentation thread started.\n");
	XmlRpc::XmlRpcValue result;
	
	switch(type)
	{
		case 0: /* Capsid */
			
			if(!m_XmlRpcClient.execute("SegmentCapsid",m_Params,result))
				QApplication::postEvent(m_CVCMainWindow,new SegmentationFailedEvent("Error segmenting virus capsid!"));
			else
				QApplication::postEvent(m_CVCMainWindow,new SegmentationFinishedEvent("Remote segmentation complete."));
			
			break;
		case 1: /* Monomer */
			
			if(!m_XmlRpcClient.execute("SegmentMonomer",m_Params,result))
				QApplication::postEvent(m_CVCMainWindow,new SegmentationFailedEvent("Error segmenting virus monomer!"));
			else
				QApplication::postEvent(m_CVCMainWindow,new SegmentationFinishedEvent("Remote segmentation complete."));
			
			break;
		case 2: /* Subunit */
		
			if(!m_XmlRpcClient.execute("SegmentSubunit",m_Params,result))
				QApplication::postEvent(m_CVCMainWindow,new SegmentationFailedEvent("Error segmenting virus subunit!"));
			else
				QApplication::postEvent(m_CVCMainWindow,new SegmentationFinishedEvent("Remote segmentation complete."));

			break;
	}

	printf("NewVolumeMainWindow::RemoteSegThread::run(): Remote segmentation thread finished.\n");

#endif
   }

}


