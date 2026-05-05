/*
  Copyright 2011 The University of Texas at Austin

        Authors: Alex Rand <arand@ices.utexas.edu>
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

#include <cvc_compat.h>

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
#include <QGridLayout>
#include <QSplitter>
#include "ui_VolumeMemoryInterface.h"

#include <VolumeRover2/VolumeMemoryInterface.h>
#include <VolumeRover2/VolumeViewer.h>

#include <VolMagick/VolMagick.h>

#include <boost/filesystem.hpp>

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iterator>
#include <iostream>

VolumeMemoryInterface::VolumeMemoryInterface(const VolMagick::Volume &vol,
                  QWidget* parent,Qt::WindowFlags flags)  
  : DataWidget(parent,flags),_ui(NULL)
{
  _ui = new Ui::VolumeMemoryInterface;
  _ui->setupUi(this);

  QGridLayout *volumeInterfaceViewerLayout = new QGridLayout(_ui->_volumeInterfaceViewer);
  _volumeInterfaceViewer = new CVC_NAMESPACE::VolumeViewer(_ui->_volumeInterfaceViewer);
  volumeInterfaceViewerLayout->addWidget(_volumeInterfaceViewer);

  setInterfaceInfo(vol);
}

VolumeMemoryInterface::~VolumeMemoryInterface() {}

void VolumeMemoryInterface::initialize(const std::string& datakey)
{
  DataWidget::initialize(datakey);
  cvcapp.properties("volumeInterfaceViewer.volumes",datakey);
}

void VolumeMemoryInterface::initialize(const boost::any& datum)
{
  setInterfaceInfo(boost::any_cast<VolMagick::Volume>(datum));
}

void VolumeMemoryInterface::setInterfaceInfo(const VolMagick::Volume &vol)
{
  _ui->_dimensionX->setText(QString("%1").arg(vol.XDim()));
  _ui->_dimensionY->setText(QString("%1").arg(vol.YDim()));
  _ui->_dimensionZ->setText(QString("%1").arg(vol.ZDim()));
  _ui->_boundingBoxMinX->setText(QString("%1").arg(vol.XMin()));
  _ui->_boundingBoxMinY->setText(QString("%1").arg(vol.YMin()));
  _ui->_boundingBoxMinZ->setText(QString("%1").arg(vol.ZMin()));
  _ui->_boundingBoxMaxX->setText(QString("%1").arg(vol.XMax()));
  _ui->_boundingBoxMaxY->setText(QString("%1").arg(vol.YMax()));
  _ui->_boundingBoxMaxZ->setText(QString("%1").arg(vol.ZMax()));
  _ui->_spanX->setText(QString("%1").arg(vol.XSpan()));
  _ui->_spanY->setText(QString("%1").arg(vol.YSpan()));
  _ui->_spanZ->setText(QString("%1").arg(vol.ZSpan()));
 
  _ui->_minValue->setText(QString("%1").arg(vol.min()));
  _ui->_maxValue->setText(QString("%1").arg(vol.max()));

  setViewerState();
}

void VolumeMemoryInterface::setViewerState()
{
  //so it responds to volumeInterfaceViewer.* property changes
  _volumeInterfaceViewer->setObjectName("volumeInterfaceViewer");

  CVC_NAMESPACE::PropertyMap properties;
  properties["volumeInterfaceViewer.rendering_mode"] = "colormapped";
  properties["volumeInterfaceViewer.shaded_rendering_enabled"] = "false";
  properties["volumeInterfaceViewer.draw_bounding_box"] = "true";
  properties["volumeInterfaceViewer.draw_subvolume_selector"] = "false";
  // arand, 6-14-2011: changed default below to 0.5
  properties["volumeInterfaceViewer.volume_rendering_quality"] = "0.5"; //[0.0,1.0]
  properties["volumeInterfaceViewer.volume_rendering_near_plane"] = "0.0";
  properties["volumeInterfaceViewer.projection_mode"] = "perspective";
  properties["volumeInterfaceViewer.draw_corner_axis"] = "true";
  properties["volumeInterfaceViewer.draw_geometry"] = "true";
  properties["volumeInterfaceViewer.draw_volumes"] = "true";
  properties["volumeInterfaceViewer.clip_geometry"] = "true";
  properties["volumeInterfaceViewer.geometries"] = "";
  properties["volumeInterfaceViewer.geometry_line_width"] = "1.2";
  properties["volumeInterfaceViewer.background_color"] = "#000000";

  //stereo related properties
  properties["volumeInterfaceViewer.io_distance"] = "0.062";
  properties["volumeInterfaceViewer.physical_distance_to_screen"] = "2.0";
  properties["volumeInterfaceViewer.physical_screen_width"] = "1.8";
  properties["volumeInterfaceViewer.focus_distance"] = "1000.0";

  //viewers_transfer_function should be a boost::shared_array<float>
  //on the data map.
  if(cvcapp.hasProperty("viewers.transfer_function"))
    properties["volumeInterfaceViewer.transfer_function"] = 
      cvcapp.properties("viewers.transfer_function");
  else
    properties["volumeInterfaceViewer.transfer_function"] = "none";

  cvcapp.addProperties(properties);
}

