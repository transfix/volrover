/*
  Copyright 2012 The University of Texas at Austin

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

#include <VolumeRover2/Application.h>
#include <cvc_compat.h>
#include <boost/format.hpp>

#include <exception>

namespace CVC_NAMESPACE
{
  bool Application::notify(QObject *receiver, QEvent *e)
  {
    try
      {
        return QApplication::notify(receiver,e);
      }
    catch(std::exception& e)
      {
        using namespace boost;
        cvcapp.log(2,str(format("%s :: Exception: %s\n")
                         % BOOST_CURRENT_FUNCTION
                         % e.what()));
        return false;
      }
  }
}
