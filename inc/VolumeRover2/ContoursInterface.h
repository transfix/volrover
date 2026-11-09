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

/* $Id: VolumeInterface.h 3602 2011-02-12 00:02:44Z transfix $ */

#ifdef USING_TILING

#ifndef __CONTOURSINTERFACE_H__
#define __CONTOURSINTERFACE_H__


#include <VolumeRover2/DataWidget.h>
#include <cvcraw_geometry/contours.h>
#include <cvcraw_geometry/io.h>
#include <cvcraw_geometry/cvcraw_geometry.h>


namespace Ui
{
  class ContoursInterface;
}

class ContoursInterface : public CVC_NAMESPACE::DataWidget
{
  Q_OBJECT

 public:
  ContoursInterface(const cvcraw_geometry::contours_t & geom= cvcraw_geometry::contours_t(),
                  QWidget* parent = 0, 
                  Qt::WindowFlags flags=Qt::WindowFlags());
  virtual ~ContoursInterface();

  virtual void initialize(const boost::any& datum)
  {
    setInterfaceInfo(boost::any_cast<cvcraw_geometry::contours_t>(datum));
  }

  void setInterfaceInfo(const cvcraw_geometry::contours_t &geom);

 signals:

public slots:
  void changeZ();

 protected:

  Ui::ContoursInterface *_ui;

private:
  // cvcraw_geometry::contours_t* _contours;
  std::string _mapName;
};

#endif

#endif
