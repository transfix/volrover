/*
  Copyright 2011 The University of Texas at Austin
  
	Authors: Deukhyun Cha <deukhyun@ices.utexas.edu>
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

#ifndef __MULTITILESERVERDIALOG_H__
#define __MULTITILESERVERDIALOG_H__

//#include <qglobal.h>

#include <QDialog>
#include <QString>

#define SYNC_THUMBNAIL_WITH_MULTITILESERVER

namespace Ui
{
  class MultiTileServerDialog;
}

#include "ui_MultiTileServerDialog.h"

class MultiTileServerDialog : public QDialog
{
  Q_OBJECT

 public:
  MultiTileServerDialog(QWidget *parent=nullptr,Qt::WindowFlags flags=Qt::WindowFlags());

  virtual ~MultiTileServerDialog();

  void setProperties();
 
  public slots:

  void updateList();
  void OutputFileSlot();
  void syncCameraSlot();
  void syncTransFuncSlot();
  void syncShadedRendSlot();
  void syncRendModeSlot();
  void initializeSlot();
  void interactiveModeSlot();
  void syncCurrentSlot();

  protected:

  Ui::MultiTileServerDialog *_ui;

  bool m_hasMultiTileServerConfigFile;
  bool m_firstCall;
  std::string m_MultiTileServerConfigFile;
  std::string m_VolumeFile;
  bool m_initialized;
  
  bool m_hasInputVolumeFile;
  std::string m_inputVolumeFile;
};

#endif
