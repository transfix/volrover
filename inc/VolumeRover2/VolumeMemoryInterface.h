/*
  Copyright 2008 The University of Texas at Austin
  
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

/* $Id: VolumeInterface.h 3602 2011-02-12 00:02:44Z transfix $ */

#ifndef __VOLUMEMEMORYINTERFACE_H__
#define __VOLUMEMEMORYINTERFACE_H__

#include <CVC/Namespace.h>


#include <VolumeRover2/DataWidget.h>
#include <VolMagick/VolMagick.h>

namespace Ui
{
  class VolumeMemoryInterface;
}

namespace CVC_NAMESPACE
{
  class VolumeViewer;
}

// 12/02/2011 -- transfix -- added a viewer widget for visual data inspection
class VolumeMemoryInterface : public CVC_NAMESPACE::DataWidget
{
  Q_OBJECT

 public:
  VolumeMemoryInterface(const VolMagick::Volume & vol= VolMagick::Volume(),
                        QWidget* parent = 0, 
                        Qt::WindowFlags flags=Qt::WindowFlags());
  virtual ~VolumeMemoryInterface();

  virtual void initialize(const std::string& datakey);
  virtual void initialize(const boost::any& datum);

  void setInterfaceInfo(const VolMagick::Volume &vol);

 protected slots:
  void setViewerState();

 signals:

 protected:

  Ui::VolumeMemoryInterface *_ui;
  CVC_NAMESPACE::VolumeViewer *_volumeInterfaceViewer;
};

#endif
