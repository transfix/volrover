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

/* $Id: RemapVoxels.h 2159 2010-07-01 22:37:42Z transfix $ */

#ifndef __REMAPVOXELS_H__
#define __REMAPVOXELS_H__

#include <QDialog>

namespace Ui
{
  class RemapVoxels;
}

class RemapVoxels : public QDialog
{
  Q_OBJECT

 public:
  RemapVoxels(QWidget* parent = nullptr,
                  Qt::WindowFlags flags=Qt::WindowFlags());
  virtual ~RemapVoxels();

  double minValue() const;
  void minValue(double);
  double maxValue() const;
  void maxValue(double);

 protected slots:
  void okSlot();

 protected:
#if QT_VERSION < 0x040000
  RemapVoxelsBase *_ui;
#else
  Ui::RemapVoxels *_ui;
#endif
};

#endif
