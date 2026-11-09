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

/* $Id: DimensionModify.h 2159 2010-07-01 22:37:42Z transfix $ */

#ifndef __DIMENSIONMODIFY_H__
#define __DIMENSIONMODIFY_H__

#include <QDialog>
#include <VolMagick/Dimension.h>

namespace Ui
{
  class DimensionModify;
}

class DimensionModify : public QDialog
{
  Q_OBJECT

 public:
  DimensionModify(QWidget* parent = nullptr,
                  Qt::WindowFlags flags=Qt::WindowFlags());
  virtual ~DimensionModify();

  VolMagick::Dimension dimension() const;
  void dimension(const VolMagick::Dimension &dim);

 protected slots:
  virtual void okSlot();

 protected:
#if QT_VERSION < 0x040000
  DimensionModifyBase *_ui;
#else
  Ui::DimensionModify *_ui;
#endif
};

#endif
