/*
  Copyright 2011 The University of Texas at Austin

        Authors: Joe Rivera <transfix@ices.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of libCVC.

  libCVC is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  libCVC is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

/* $Id$ */

//NOTE: this header requires Qt for any translation unit including it

#include <cvc_compat.h>

#include <QEvent>

#include <string>
#include <boost/any.hpp>

namespace CVC_NAMESPACE
{
  struct CVCEvent : public QEvent
  {
    public:
    CVCEvent(const std::string n = std::string(),
             const boost::any d = boost::any())
      : QEvent(QEvent::User), name(n), data(d) {}
    
    std::string name;
    boost::any  data;
  };
}

