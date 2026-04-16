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
#include "ui_GeometryInterface.h"

#include <VolumeRover2/GeometryInterface.h>
#include <VolumeRover2/VolumeViewer.h>

#include <cvcraw_geometry/cvcgeom.h>
#include <cvcraw_geometry/io.h>
#include <cvcraw_geometry/cvcraw_geometry.h>

#include <log4cplus/logger.h>
#include <CVC/log4cplus_compat.h>
#include <log4cplus/loggingmacros.h>

#include <boost/filesystem.hpp>

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iterator>
#include <iostream>

// 
// 3-31-2011: arand, initial implementation
//
GeometryInterface::GeometryInterface(const cvcraw_geometry::cvcgeom_t & geom,
                  QWidget* parent,Qt::WindowFlags flags)  
  : DataWidget(parent,flags),_ui(NULL)
{
  _ui = new Ui::GeometryInterface;
  _ui->setupUi(this);

  QGridLayout *geometryInterfaceViewerLayout = new QGridLayout(_ui->_geometryInterfaceViewer);
  _geometryInterfaceViewer = new CVC_NAMESPACE::VolumeViewer(_ui->_geometryInterfaceViewer);
  geometryInterfaceViewerLayout->addWidget(_geometryInterfaceViewer);

  setInterfaceInfo(geom);
}

GeometryInterface::~GeometryInterface() {}

void GeometryInterface::initialize(const std::string& datakey)
{
  DataWidget::initialize(datakey);
  cvcapp.properties("geometryInterfaceViewer.geometries",datakey);
}

void GeometryInterface::initialize(const boost::any& datum)
{
  setInterfaceInfo(boost::any_cast<cvcraw_geometry::cvcgeom_t>(datum));
}

void GeometryInterface::setInterfaceInfo(const cvcraw_geometry::cvcgeom_t &geom)
{
  _ui->_numVerts->setText(QString("%1").arg(geom.num_points()));
  _ui->_numLines->setText(QString("%1").arg(geom.num_lines()));
  _ui->_numTriangles->setText(QString("%1").arg(geom.num_triangles()));
  _ui->_numQuads->setText(QString("%1").arg(geom.num_quads()));
 
  CVCGEOM_NAMESPACE::cvcgeom_t::point_t tmpMin = geom.min_point();
  _ui->_minX->setText(QString("%1").arg(tmpMin[0]));
  _ui->_minY->setText(QString("%1").arg(tmpMin[1]));
  _ui->_minZ->setText(QString("%1").arg(tmpMin[2]));

  CVCGEOM_NAMESPACE::cvcgeom_t::point_t tmpMax = geom.max_point();
  _ui->_maxX->setText(QString("%1").arg(tmpMax[0]));
  _ui->_maxY->setText(QString("%1").arg(tmpMax[1]));
  _ui->_maxZ->setText(QString("%1").arg(tmpMax[2]));

  // TODO: deal with the normals/colors checkboxes...
  // arand: done, 4-21-2011
  _ui->_colors->setChecked(false);
  _ui->_normals->setChecked(false);
  
  if (geom.const_normals().size() > 0)
      _ui->_normals->setChecked(true);
  if (geom.const_colors().size() > 0)
      _ui->_colors->setChecked(true);

  setViewerState();
}

void GeometryInterface::setViewerState()
{
  using namespace std;
  using namespace boost;

  static log4cplus::Logger logger = FUNCTION_LOGGER;

  //so it responds to volumeInterfaceViewer.* property changes
  _geometryInterfaceViewer->setObjectName("geometryInterfaceViewer");

  CVC_NAMESPACE::PropertyMap properties;
  properties["geometryInterfaceViewer.rendering_mode"] = "colormapped";
  properties["geometryInterfaceViewer.shaded_rendering_enabled"] = "false";
  properties["geometryInterfaceViewer.draw_bounding_box"] = "true";
  properties["geometryInterfaceViewer.draw_subvolume_selector"] = "false";
  // arand, 6-14-2011: changed default below to 0.5
  properties["geometryInterfaceViewer.volume_rendering_quality"] = "0.5"; //[0.0,1.0]
  properties["geometryInterfaceViewer.volume_rendering_near_plane"] = "0.0";
  properties["geometryInterfaceViewer.projection_mode"] = "perspective";
  properties["geometryInterfaceViewer.draw_corner_axis"] = "true";
  properties["geometryInterfaceViewer.draw_geometry"] = "true";
  properties["geometryInterfaceViewer.draw_volumes"] = "true";
  properties["geometryInterfaceViewer.clip_geometry"] = "true";
  properties["geometryInterfaceViewer.volumes"] = "";
  properties["geometryInterfaceViewer.geometry_line_width"] = "1.2";
  properties["geometryInterfaceViewer.background_color"] = "#000000";

  //stereo related properties
  properties["geometryInterfaceViewer.io_distance"] = "0.062";
  properties["geometryInterfaceViewer.physical_distance_to_screen"] = "2.0";
  properties["geometryInterfaceViewer.physical_screen_width"] = "1.8";
  properties["geometryInterfaceViewer.focus_distance"] = "1000.0";

  //viewers_transfer_function should be a boost::shared_array<float>
  //on the data map.
  properties["geometryInterfaceViewer.transfer_function"] = "none";

  cvcapp.addProperties(properties);

  //set up the volume bounding box to encompass the geometry
  vector<cvcraw_geometry::cvcgeom_t> geoms =
    _geometryInterfaceViewer->getGeometriesFromDatamap();
  LOG4CPLUS_TRACE(logger, "geoms.size() == " << geoms.size());
  // cvcapp.log(3,str(format("%s :: geoms.size() == %s\n")
  //                  % BOOST_CURRENT_FUNCTION
  //                  % geoms.size()));

  if(geoms.empty())
    {
      //just use the default bounding box of a volume
      _geometryInterfaceViewer->colorMappedVolume(VolMagick::Volume());
      _geometryInterfaceViewer->rgbaVolumes(vector<VolMagick::Volume>());
    }
  else
    {
      VolMagick::BoundingBox bbox;
      cvcraw_geometry::cvcgeom_t initial_geom;
      
      //find a non empty geom to supply an initial bbox
      for (const auto& geom : geoms)
        if(!geom.empty())
          {
            initial_geom = geom;
            break;
          }
      
      cvcraw_geometry::cvcgeom_t::point_t minpt;
      cvcraw_geometry::cvcgeom_t::point_t maxpt;
      if(initial_geom.empty())
        {
          //no non empty geometry so lets just use the default box
          _geometryInterfaceViewer->colorMappedVolume(VolMagick::Volume());
          _geometryInterfaceViewer->rgbaVolumes(vector<VolMagick::Volume>());
        }
      else
        {
          minpt = initial_geom.min_point();
          maxpt = initial_geom.max_point();
          // arand: slightly enlarging
          double eps0 = (maxpt[0]-minpt[0])/20.0;
          double eps1 = (maxpt[1]-minpt[1])/20.0;
          double eps2 = (maxpt[2]-minpt[2])/20.0;
          
          bbox = VolMagick::BoundingBox(minpt[0]-eps0,minpt[1]-eps1,minpt[2]-eps2,
                                        maxpt[0]+eps0,maxpt[1]+eps1,maxpt[2]+eps2);
        }
      
      //build a bounding box that encompasses all bounding boxes of the geometries
      for (const auto& geom : geoms) {
          if(geom.empty()) continue;
          minpt = geom.min_point();
          maxpt = geom.max_point();
          VolMagick::BoundingBox geobox(minpt[0],minpt[1],minpt[2],
                                        maxpt[0],maxpt[1],maxpt[2]);
          
          cvcapp.log(5,str(boost::format("geobox: (%f,%f,%f) (%f,%f,%f)")
                           % geobox.minx % geobox.miny % geobox.minz
                           % geobox.maxx % geobox.maxy % geobox.maxz));
          cvcapp.log(5,str(boost::format("bbox: (%f,%f,%f) (%f,%f,%f)")
                           % bbox.minx % bbox.miny % bbox.minz
                           % bbox.maxx % bbox.maxy % bbox.maxz));
          
          bbox += geobox;
        }
      
      VolMagick::Volume newvol;
      newvol.boundingBox(bbox);
      _geometryInterfaceViewer->colorMappedVolume(newvol);
      _geometryInterfaceViewer->rgbaVolumes(vector<VolMagick::Volume>(4,newvol));
    }
}
