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

#include <QString>
#include <QFileDialog>
#include <QLayout>
#include <qtoolbar.h>
#include <ColorTable/ColorTable.h>
//#include <VolumeGridRover/MappedRawIVFile.h>
//#include <VolumeGridRover/MappedRawVFile.h>
#include <VolMagick/VolMagick.h>

#include <VolumeGridRover/VolumeGridRover.h>
#include <VolumeGridRover/VolumeGridRoverMainWindow.h>
#include <cvc_compat.h>

static inline unsigned char mapToChar(double val)
{
  int inval = int(val*255);
  inval = ( inval<255 ? inval : 255 );
  inval = ( inval>0 ? inval : 0);
  return (unsigned char) inval;
}


VolumeGridRoverMainWindow::VolumeGridRoverMainWindow(QWidget* parent, Qt::WindowFlags fl)
  : QWidget(parent,fl)
{
  m_VolumeGridRoverLayout = new QGridLayout( this );
  m_VolumeGridRover = new VolumeGridRover( this);
  m_VolumeGridRoverLayout->addWidget(m_VolumeGridRover);

  m_ColorToolbar = new QToolBar("Transfer Function",this);
  m_ColorToolbar->setEnabled(true);
//Q3Err:not support in Qt4
//  m_ColorToolbar->setResizeEnabled(true);
//  m_ColorToolbar->setMovingEnabled(true);
//  m_ColorToolbar->setHorizontallyStretchable(true);
//  m_ColorToolbar->setOpaqueMoving(false);
  m_ColorTable = new CVC::ColorTable(m_ColorToolbar,"m_ColorTable");
  m_ColorTable->setSizePolicy(QSizePolicy(QSizePolicy::Preferred,QSizePolicy::Preferred,0,0,m_ColorTable->sizePolicy().hasHeightForWidth()));
  m_ColorTable->setSpectrumFunctions(NULL,NULL,NULL,NULL,NULL);
  connect(m_ColorTable,SIGNAL(functionExploring()),SLOT(functionChangedSlot()));
  connect(m_ColorTable,SIGNAL(functionExploring()),SLOT(functionChangedSlot()));
  connect(m_ColorTable,SIGNAL(functionChanged()),SLOT(functionChangedSlot()));
  connect(m_ColorTable,SIGNAL(everythingChanged()),SLOT(functionChangedSlot()));
}

VolumeGridRoverMainWindow::~VolumeGridRoverMainWindow(){}

void VolumeGridRoverMainWindow::fileOpen()
{
  QString filename = QFileDialog::getOpenFileName(QString::null,
					   "Volume Files (*.rawiv *.rawv *.cvc *.mrc)",
					   this,
					   "open file dialog",
					   "Choose a volume file" );
  if(filename == QString::null) return;
  
  functionChangedSlot(); /* set the current transfer function so the slice appears correct */
 
#if 0 
  if(filename.endsWith(".rawiv",false))
    {
      m_MappedVolumeFile = new MappedRawIVFile(filename.ascii(),true,true);
      if(!m_MappedVolumeFile->isValid())
	cvcapp.log(5, boost::str(boost::format("VolumeFile::VolumeFile(): Could not load '%s'")%filename.ascii()));
    }
  else if(filename.endsWith(".rawv",false))
    {
      m_MappedVolumeFile = new MappedRawVFile(filename.ascii(),true,true);
      if(!m_MappedVolumeFile->isValid())
	cvcapp.log(5, boost::str(boost::format("VolumeFile::VolumeFile(): Could not load '%s'")%filename.ascii()));
    }
  else /* try to figure out the volume type */
    {
      m_MappedVolumeFile = new MappedRawIVFile(filename.ascii(),true,false);
      if(m_MappedVolumeFile->isValid()) return;
      delete m_MappedVolumeFile;

      m_MappedVolumeFile = new MappedRawVFile(filename.ascii(),true,false);
      if(m_MappedVolumeFile->isValid()) return;
      delete m_MappedVolumeFile;

      m_MappedVolumeFile = NULL;
      cvcapp.log(5, "VolumeFile::VolumeFile(): m_MappedVolumeFile == NULL");
    }

  
  if(m_VolumeGridRover->setVolume(m_MappedVolumeFile))
    {
      functionChangedSlot();
      setCaption(filename + " - Volume Grid Rover");
    }
  else
    setCaption("Volume Grid Rover");
#endif

  VolMagick::VolumeFileInfo vfi(cvcapp, filename.ascii());
  m_VolumeFileInfo = vfi;
  m_VolumeGridRover->setVolume(vfi);
  setCaption(filename + " - Volume Grid Rover");
}

void VolumeGridRoverMainWindow::functionChangedSlot()
{
  double map[256*4];
  unsigned char byte_map[256*4];
  unsigned int c;
  m_ColorTable->GetTransferFunction(map, 256);
  for (c=0; c<256; c++)
    {
      byte_map[c*4+0] = mapToChar(map[c*4+0]);
      byte_map[c*4+1] = mapToChar(map[c*4+1]);
      byte_map[c*4+2] = mapToChar(map[c*4+2]);
      byte_map[c*4+3] = mapToChar(map[c*4+3]);
    }
  m_VolumeGridRover->setTransferFunction(byte_map);
}
