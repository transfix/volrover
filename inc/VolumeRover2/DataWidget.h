/*
  Copyright 2011 The University of Texas at Austin
  
	Authors: Jose Rivera <transfix@ices.utexas.edu>
	Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of Volume Rover.

  Volume Rover is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  Volume Rover is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with iotree; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/* $Id: DataWidget.h 4910 2011-12-03 02:51:47Z transfix $ */


#ifndef __DATAWIDGET_H__
#define __DATAWIDGET_H__

#include <cvc_compat.h>
#include <QFrame>
#include <boost/any.hpp>

namespace CVC_NAMESPACE
{
  // 12/02/2011 -- transfix -- added initialize with a string argument
  class DataWidget : public QFrame
  {
  public:
    DataWidget(QWidget *parent = nullptr, Qt::WindowFlags flags=Qt::WindowFlags()) : QFrame(parent,flags) {}
    virtual ~DataWidget() {}

    virtual void initialize(const std::string& datakey)
    {
      initialize(cvcapp.data(datakey));
    }
    virtual void initialize(const boost::any& datum) = 0;
  };
}

#endif
