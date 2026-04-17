/*
  Copyright 2008-2011 The University of Texas at Austin

        Authors: Joe Rivera <transfix@ices.utexas.edu>
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

/* $Id: VolumeInterface.cpp 3602 2011-02-12 00:02:44Z transfix $ */

#ifdef USING_TILING

#include <qglobal.h>

#include <QString>
#include <QLabel>
#include <QLineEdit>
#include <QTreeWidget>
#include <QPushButton>
#include <QMessageBox>
#include <QButtonGroup>
#include <QComboBox>
#include <QCheckBox>
#include <QFileInfo>
#include "ui_ContoursInterface.h"

#include <VolumeRover2/ContoursInterface.h>
#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>

#include <CVC/App.h>
#include <cvcraw_geometry/cvcgeom.h>
#include <cvcraw_geometry/io.h>
#include <cvcraw_geometry/cvcraw_geometry.h>

#include <boost/filesystem.hpp>

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iterator>
#include <iostream>

ContoursInterface::ContoursInterface(const cvcraw_geometry::contours_t & geom,
                  QWidget* parent,Qt::WindowFlags flags)  
  : DataWidget(parent,flags),_ui(NULL)
{
  _ui = new Ui::ContoursInterface;
  _ui->setupUi(this);

  setInterfaceInfo(geom);

  connect(_ui->_z_spacing, SIGNAL(editingFinished()), SLOT(changeZ()));
}

ContoursInterface::~ContoursInterface() {}

void ContoursInterface::setInterfaceInfo(const cvcraw_geometry::contours_t &contours)
{
  _ui->_z_spacing->setText(QString("%1").arg(contours.z_scale()));
  _mapName = contours.name();
}

void ContoursInterface::changeZ()
{
  log4cplus::Logger logger = log4cplus::Logger::getInstance("VolumeRover2.ContoursInterface.changeZ");

  double zspacing = _ui->_z_spacing->displayText().toDouble();
  LOG4CPLUS_TRACE(logger, "zSpacing = " << zspacing);
  // string dataset = key.substr(0, key.rfind('.'));
  // LOG4CPLUS_TRACE(logger, "dataset = " << zspacing);
  cvcraw_geometry::contours_t c = cvcapp.data<cvcraw_geometry::contours_t>(_mapName);
  c.set_z_scale(zspacing);
  cvcapp.data(_mapName, c);

  // LOG4CPLUS_TRACE(logger, "here");
  // double z = _ui->_z_spacing->displayText().toDouble();
  // LOG4CPLUS_TRACE(logger, "got z: " << z);
  // _contours->set_z_spacing(z);
  // LOG4CPLUS_TRACE(logger, "set z spacing");
}

#endif
