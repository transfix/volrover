/*
  Copyright 2008-2011 The University of Texas at Austin

        Authors: Joe Rivera <transfix@ices.utexas.edu>
                 Deukhyun Cha <deukhyun@ices.utexas.edu>
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

#ifndef __MULTITILESERVERDIALOG_CPP__
#define __MULTITILESERVERDIALOG_CPP__

#include <cvc_compat.h>
#include <VolumeRover2/MultiTileServerDialog.h>

#include <VolMagick/VolMagick.h>

#include <QFileDialog>
#include <QMessageBox>

#include "ui_MultiTileServerDialog.h"

#include <iostream> 
using namespace std;
using namespace CVC_NAMESPACE;

MultiTileServerDialog::MultiTileServerDialog(QWidget *parent,Qt::WindowFlags flags) 
  : QDialog(parent, flags) {

  _ui = new Ui::MultiTileServerDialog;
  _ui->setupUi(this);
  
  // connect slots and signals?
  connect(_ui->m_UpdateList,
          SIGNAL(clicked()),
          SLOT(updateList()));

  connect(_ui->m_ServerConfigFileOpen,
	  SIGNAL(clicked()),
	  SLOT(OutputFileSlot()));

  connect(_ui->m_SyncCamera,
          SIGNAL(clicked()),
          SLOT(syncCameraSlot()));
  connect(_ui->m_SyncTransferFunc,
          SIGNAL(clicked()),
          SLOT(syncTransFuncSlot()));
  connect(_ui->m_SyncShadedRender,
          SIGNAL(clicked()),
          SLOT(syncShadedRendSlot()));
  connect(_ui->m_SyncRenderMode,
          SIGNAL(clicked()),
          SLOT(syncRendModeSlot()));
 
  connect(_ui->m_Initialize,
          SIGNAL(clicked()),
          SLOT(initializeSlot()));
  connect(_ui->m_InteractiveMode,
          SIGNAL(clicked()),
          SLOT(interactiveModeSlot()));
  connect(_ui->m_SyncCurrent,
          SIGNAL(clicked()),
          SLOT(syncCurrentSlot()));

  _ui->m_SyncCamera->setChecked( false );
  _ui->m_SyncTransferFunc->setChecked( false );
  _ui->m_SyncShadedRender->setChecked( false );
  _ui->m_SyncRenderMode->setChecked( false );

  m_hasMultiTileServerConfigFile = false;
  m_firstCall = true;
  m_initialized = false;
  m_hasInputVolumeFile = false;

}

MultiTileServerDialog::~MultiTileServerDialog() {
  delete _ui;
}

void MultiTileServerDialog::updateList() {
  int nFiles;
  std::string nFileString = cvcapp.properties("number_of_file_read");
  if( nFileString.size() )
  {
     sscanf( nFileString.data(), "%d", &nFiles);
     fprintf( stderr, "nfiles: %d\n", nFiles );
     for( int i = 0; i < nFiles; i++ )
     {
        _ui->m_volumeFileList->clear();
        char str[1024];
        sprintf( str, "file_%d_fullPath", i );
        std::string filepath = cvcapp.properties( str );
        fprintf( stderr, "path: %s\n", filepath.data() );
        _ui->m_volumeFileList->addItem( QString::fromStdString( filepath ) );
     }
  }

  if( m_hasInputVolumeFile )
    _ui->m_InputVolumeFile->setText( QString( m_inputVolumeFile.data() ) );

  if( m_hasMultiTileServerConfigFile )
    _ui->m_ServerConfigFile->setText( QString( m_MultiTileServerConfigFile.data() ) );
}

void MultiTileServerDialog::OutputFileSlot() {
  _ui->m_ServerConfigFile->setText(QFileDialog::getOpenFileName(this,
						    "Multi-tile Server Config File"));
}

void MultiTileServerDialog::syncCameraSlot() {
   #ifdef SYNC_THUMBNAIL_WITH_MULTITILESERVER
   if( _ui->m_SyncCamera->isChecked() ) cvcapp.properties("thumbnail.syncCamera_with_multiTileServer","true");
   else cvcapp.properties("thumbnail.syncCamera_with_multiTileServer","false");
   #else
   if( _ui->m_SyncCamera->isChecked() ) cvcapp.properties("zoomed.syncCamera_with_multiTileServer","true");
   else cvcapp.properties("zoomed.syncCamera_with_multiTileServer","false");
   #endif
}

void MultiTileServerDialog::syncTransFuncSlot() {
   #ifdef SYNC_THUMBNAIL_WITH_MULTITILESERVER
   if( _ui->m_SyncTransferFunc->isChecked() ) cvcapp.properties("thumbnail.syncTransferFunc_with_multiTileServer", "true");
   else cvcapp.properties("thumbnail.syncTransferFunc_with_multiTileServer", "false");
   #else
   if( _ui->m_SyncTransferFunc->isChecked() ) cvcapp.properties("zoomed.syncTransferFunc_with_multiTileServer", "true");
   else cvcapp.properties("zoomed.syncTransferFunc_with_multiTileServer", "false");
   #endif
}

void MultiTileServerDialog::syncShadedRendSlot() {
   #ifdef SYNC_THUMBNAIL_WITH_MULTITILESERVER
   if( _ui->m_SyncShadedRender->isChecked() ) cvcapp.properties("thumbnail.syncShadedRender_with_multiTileServer", "true");
   else cvcapp.properties("thumbnail.syncShadedRender_with_multiTileServer", "false");
   #else
   if( _ui->m_SyncShadedRender->isChecked() ) cvcapp.properties("zoomed.syncShadedRender_with_multiTileServer", "true");
   else cvcapp.properties("zoomed.syncShadedRender_with_multiTileServer", "false");
   #endif
}

void MultiTileServerDialog::syncRendModeSlot() {
   #ifdef SYNC_THUMBNAIL_WITH_MULTITILESERVER
   if( _ui->m_SyncRenderMode->isChecked() ) cvcapp.properties("thumbnail.syncRenderMode_with_multiTileServer", "true");
   else cvcapp.properties("thumbnail.syncRenderMode_with_multiTileServer", "false");
   #else
   if( _ui->m_SyncRenderMode->isChecked() ) cvcapp.properties("zoomed.syncRenderMode_with_multiTileServer", "true");
   else cvcapp.properties("zoomed.syncRenderMode_with_multiTileServer", "false");
   #endif
}

void MultiTileServerDialog::initializeSlot() {
   static int calls = 0;
   // get initial paramters from UI
   setProperties();
   calls = (++calls) % 10 + 1;

   #ifdef SYNC_THUMBNAIL_WITH_MULTITILESERVER
   cvcapp.properties<int>("thumbnail.syncInitial_multiTileServer", calls);
   #else
   cvcapp.properties<int>("zoomed.syncInitial_multiTileServer", calls);
   #endif
   m_initialized = true;
}

void MultiTileServerDialog::interactiveModeSlot() {
   if( !m_initialized ) {
      if( _ui->m_InteractiveMode->isChecked() )
         _ui->m_InteractiveMode->setChecked( false );

      fprintf( stderr, "initialize first\n");
      return;
   }

   if( _ui->m_InteractiveMode->isChecked() ) {
      _ui->m_SyncCurrent->setEnabled( false );
      #ifdef SYNC_THUMBNAIL_WITH_MULTITILESERVER
      cvcapp.properties("thumbnail.interactiveMode_with_multiTileServer", "true");
      #else
      cvcapp.properties("zoomed.interactiveMode_with_multiTileServer", "true");
      #endif
   }
   else {
      _ui->m_SyncCurrent->setEnabled( true );
      #ifdef SYNC_THUMBNAIL_WITH_MULTITILESERVER
      cvcapp.properties("thumbnail.interactiveMode_with_multiTileServer", "false");
      #else
      cvcapp.properties("zoomed.interactiveMode_with_multiTileServer", "false");
      #endif
   }
}

void MultiTileServerDialog::syncCurrentSlot() {
   static int calls = 0;

   if( !m_initialized ) {
      fprintf( stderr, "initialize first\n");
      return;
   }

   calls = (++calls) % 10 + 1;
   #ifdef SYNC_THUMBNAIL_WITH_MULTITILESERVER
   cvcapp.properties<int>("thumbnail.syncMode_with_multiTileServer", calls);
   #else
   cvcapp.properties<int>("zoomed.syncMode_with_multiTileServer", calls);
   #endif
}

void MultiTileServerDialog::setProperties()
{
  m_MultiTileServerConfigFile = _ui->m_ServerConfigFile->text().toStdString();
  m_hasMultiTileServerConfigFile = true;
  if( _ui->m_UseInputVolumeFile->isChecked() ) {
     m_VolumeFile = m_inputVolumeFile = _ui->m_InputVolumeFile->text().toStdString();
     m_hasInputVolumeFile = true;
  }
  else
     m_VolumeFile = _ui->m_volumeFileList->currentItem()->text().toStdString();

  FILE *fp = fopen( m_MultiTileServerConfigFile.data(), "r");
  if( !fp ) {
     fprintf( stderr, "multi tile server config file open fail [%s]\n", m_MultiTileServerConfigFile.data() );
     #ifdef SYNC_THUMBNAIL_WITH_MULTITILESERVER
     cvcapp.properties("thumbnail.sync_with_multiTileServer", "false");
     #else
     cvcapp.properties("zoomed.sync_with_multiTileServer", "false");
     #endif
     return;
  }

  int nTilesX = 1, nTilesY = 1;
  int nServer = 0;
  char line[1025];
  char *ret = fgets( line, 1024, fp );

  sscanf(line, "%d %d", &nTilesX, &nTilesY );
  nServer = nTilesX * nTilesY;
  char **host = new char*[nServer];
  int *port   = new int[nServer];
  for( int i = 0; i < nServer;  i++ )
  {
     host[ i ] = new char[1024];
     ret = fgets( line, 1024, fp );
     sscanf( line, "%s %d", host[i], &port[i] );
     fprintf( stderr, "host: %s, port: %d\n", host[i], port[i]);
  }

  if( m_firstCall )
  {
     PropertyMap properties;
     properties["MultiTileServer.volFile"] = m_VolumeFile.data();
     #ifdef SYNC_THUMBNAIL_WITH_MULTITILESERVER
     properties["thumbnail.syncVolume_with_multiTileServer"] = m_VolumeFile.data();
     #else
     properties["zoomed.syncVolume_with_multiTileServer"] = m_VolumeFile.data();
     #endif
     char tmpc[2048], tmpc1[2048];
     sprintf( tmpc, "%d", nTilesX );
     properties["MultiTileServer.nTilesX"] = tmpc;
     sprintf( tmpc, "%d", nTilesY );
     properties["MultiTileServer.nTilesY"] = tmpc;
     sprintf( tmpc, "%d", nServer );
     properties["MultiTileServer.nServer"] = tmpc;
     for( int i = 0; i < nServer; i++ )
     {
        sprintf( tmpc, "MultiTileServer.host_%d", i );
        properties[tmpc] = host[i];
        sprintf( tmpc, "MultiTileServer.port_%d", i );
        sprintf( tmpc1, "%d", port[i] );
        properties[tmpc] = tmpc1;
     }
     cvcapp.addProperties( properties );
     m_firstCall = false;
  }
  else
  {
     cvcapp.properties("MultiTileServer.volFile", m_VolumeFile.data());
     #ifdef SYNC_THUMBNAIL_WITH_MULTITILESERVER
     cvcapp.properties("thumbnail.syncVolume_with_multiTileServer", m_VolumeFile.data());
     #else
     cvcapp.properties("zoomed.syncVolume_with_multiTileServer", m_VolumeFile.data());
     #endif
     char tmpc[2048], tmpc1[2048];
     sprintf( tmpc, "%d", nServer );
     cvcapp.properties("MultiTileServer.nServer", tmpc);
     for( int i = 0; i < nServer; i++ )
     {
        sprintf( tmpc, "MultiTileServer.host_%d", i );
        cvcapp.properties(tmpc , host[i]);
        sprintf( tmpc, "MultiTileServer.port_%d", i );
        sprintf( tmpc1, "%d", port[i] );
        cvcapp.properties(tmpc, tmpc1);
     }
  }
}

#endif
