/*
  Copyright 2011 The University of Texas at Austin
  
	Authors: Alex Rand <arand@ices.utexas.edu>
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

#ifndef __RECONSTRUCTIONDIALOG_H__
#define __RECONSTRUCTIONDIALOG_H__

//#include <qglobal.h>

#include <QDialog>
#include <QString>

namespace Ui
{
  class ReconstructionDialog;
}

#include "ui_ReconstructionDialog.h"

class ReconstructionDialog : public QDialog
{
  Q_OBJECT

 public:
  ReconstructionDialog(QWidget *parent=nullptr,Qt::WindowFlags flags=Qt::WindowFlags());

  virtual ~ReconstructionDialog();

  void RunReconstruction();

  
  public slots:

  void InitFileSlot();
  void DataFileSlot();

  protected:

  Ui::ReconstructionDialog *_ui;

};

#endif
