/*
  Copyright 2008-2011 The University of Texas at Austin
  
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

/* $Id: VolumeInterface.h 3602 2011-02-12 00:02:44Z transfix $ */

#ifndef __GEOMETRYINTERFACE_H__
#define __GEOMETRYINTERFACE_H__

#include <CVC/Namespace.h>
#include <VolumeRover2/DataWidget.h>
#include <cvcraw_geometry/cvcgeom.h>
#include <cvcraw_geometry/io.h>
#include <cvcraw_geometry/cvcraw_geometry.h>

namespace Ui
{
  class GeometryInterface;
}

namespace CVC_NAMESPACE
{
  class VolumeViewer;
}

// 12/02/2011 -- transfix -- added a viewer widget for visual data inspection
class GeometryInterface : public CVC_NAMESPACE::DataWidget
{
  Q_OBJECT

 public:
  GeometryInterface(const cvcraw_geometry::cvcgeom_t & geom= cvcraw_geometry::cvcgeom_t(),
                  QWidget* parent = nullptr, 
                  Qt::WindowFlags flags=Qt::WindowFlags());
  virtual ~GeometryInterface();

  virtual void initialize(const std::string& datakey);
  virtual void initialize(const boost::any& datum);

  void setInterfaceInfo(const cvcraw_geometry::cvcgeom_t &geom);

 protected slots:
  void setViewerState();

 signals:

 protected:

  Ui::GeometryInterface *_ui;
  CVC_NAMESPACE::VolumeViewer *_geometryInterfaceViewer;
};

#endif
