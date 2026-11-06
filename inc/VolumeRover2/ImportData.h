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

/* $Id: ImportData.h 2159 2010-07-01 22:37:42Z transfix $ */

#ifndef __IMPORTDATA_H__
#define __IMPORTDATA_H__

#include <QDialog>

namespace Ui
{
  class ImportData;
}

class ImportData : public QDialog
{
  Q_OBJECT

 public:
  ImportData(QWidget* parent = 0,
#if QT_VERSION < 0x040000
                  const char* name = 0, WFlags f = WType_TopLevel
#else
                  Qt::WindowFlags flags={}
#endif
                  );
  virtual ~ImportData();

 protected slots:
  virtual void importFileSlot();
  virtual void okSlot();

 protected:
#if QT_VERSION < 0x040000
  ImportDataBase *_ui;
#else
  Ui::ImportData *_ui;
#endif
};

#endif
