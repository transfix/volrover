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

#ifdef USING_TILING

#ifndef __CONTOURTILERDIALOG_H__
#define __CONTOURTILERDIALOG_H__

//#include <qglobal.h>

#include <QDialog>
#include <QString>

namespace Ui
{
  class ContourTilerDialog;
}

#include "ui_ContourTilerDialog.h"

class ContourTilerDialog : public QDialog
{
  Q_OBJECT

 public:
  ContourTilerDialog(QWidget *parent=nullptr,Qt::WindowFlags flags=Qt::WindowFlags());

  virtual ~ContourTilerDialog();

  void RunContourTiler();

  
  public slots:
  void geometryChangedSlot(int index);
  void outputDirSlot();

  protected:
  Ui::ContourTilerDialog *_ui;

};

#endif

#endif
