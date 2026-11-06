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

#ifndef __LBIEQUALITYIMPROVEMENTDIALOG_H__
#define __LBIEQUALITYIMPROVEMENTDIALOG_H__

//#include <qglobal.h>

#include <QDialog>
#include <QString>

namespace Ui
{
  class LBIEQualityImprovementDialogBase;
}

#include "ui_LBIE_qualityImprovement.h"

class LBIEQualityImprovementDialog : public QDialog
{
  Q_OBJECT

 public:
  LBIEQualityImprovementDialog(QWidget *parent=nullptr,Qt::WindowFlags flags=Qt::WindowFlags());

  virtual ~LBIEQualityImprovementDialog();

  void RunLBIEQualityImprovement();

  
  public slots:

  //void OutputFileSlot();

  protected:

  Ui::LBIEQualityImprovementDialogBase *_ui;

};

#endif
