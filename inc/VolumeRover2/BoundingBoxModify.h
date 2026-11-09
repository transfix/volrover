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

/* $Id: BoundingBoxModify.h 2015 2010-06-21 16:28:56Z transfix $ */

#ifndef __BOUNDINGBOXMODIFY_H__
#define __BOUNDINGBOXMODIFY_H__

#include <QDialog>
#include <VolMagick/BoundingBox.h>

namespace Ui
{
  class BoundingBoxModify;
}

class BoundingBoxModify : public QDialog
{
  Q_OBJECT

 public:
  BoundingBoxModify(QWidget* parent = nullptr, 
                    Qt::WindowFlags flags=Qt::WindowFlags());
  virtual ~BoundingBoxModify();

  VolMagick::BoundingBox boundingBox() const;
  void boundingBox(const VolMagick::BoundingBox& bbox);

  double centerPointX() const;
  double centerPointY() const;
  double centerPointZ() const;

  void centerPoint(double,double,double);

  bool usingCenterPoint() const;

 protected slots:
  void okSlot();

 protected:

#if QT_VERSION < 0x040000
  BoundingBoxModifyBase *_ui;
#else
  Ui::BoundingBoxModify *_ui;
#endif
};

#endif
