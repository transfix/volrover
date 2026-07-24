/*
  Copyright 2012 The University of Texas at Austin

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

/* $Id: BoundingBoxClip.h 5692 2012-06-01 22:39:25Z transfix $ */

#include <cvc_compat.h>
#include <cvc_compat.h>

namespace CVC_NAMESPACE
{
  //functor for clipping via a bounding box
  class BoundingBoxClip
  {
  public:
    BoundingBoxClip(const CVC::BoundingBox bbox)
      : _bbox(bbox) {}

    void operator()();

    CVC::BoundingBox boundingBox() const { return _bbox; }

  private:
    CVC::BoundingBox _bbox;
  };
}
