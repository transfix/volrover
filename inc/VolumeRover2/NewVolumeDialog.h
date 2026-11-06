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

/* $Id: NewVolumeDialog.h 2122 2010-06-28 02:06:38Z transfix $ */

#ifndef __NEWVOLUMEDIALOG_H__
#define __NEWVOLUMEDIALOG_H__

#include <QDialog>

namespace Ui
{
  class NewVolumeDialog;
}

#include <VolMagick/BoundingBox.h>
#include <VolMagick/Dimension.h>

#include <string>

class NewVolumeDialog : public QDialog
{
  Q_OBJECT

 public:
  NewVolumeDialog(QWidget* parent = 0,
#if QT_VERSION < 0x040000
                  const char* name = 0, WFlags f = WType_TopLevel
#else
                  Qt::WindowFlags flags={}
#endif
                  );
  virtual ~NewVolumeDialog();

  bool createNewVolume() const;

  VolMagick::Dimension dimension() const;
  VolMagick::BoundingBox boundingBox() const;
  VolMagick::VoxelType variableType() const;
  std::string variableName() const;
  std::string filename() const;
  std::string volumeCopyFilename() const;
  bool extractSubVolume() const;
  
  enum ExtractSubVolumeMethod { INDICES, BOUNDING_BOX };
  ExtractSubVolumeMethod extractSubVolumeMethod() const;

  //INDICES
  VolMagick::IndexBoundingBox extractIndexSubVolume() const;

  //BOUNDING_BOX
  VolMagick::BoundingBox extractSubVolumeBoundingBox() const;
  VolMagick::Dimension extractSubVolumeDimension() const;
  
 protected slots:
  void okSlot();
  void fileSlot();
  void volumeCopyFilenameSlot();

  void acquireVolumeInfo(bool);

 protected:
#if QT_VERSION < 0x040000
  NewVolumeDialogBase *_ui;
#else
  Ui::NewVolumeDialog *_ui;
#endif
};

#endif
