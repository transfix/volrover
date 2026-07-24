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

/* $Id: SceneLight.h 5746 2012-06-17 23:21:45Z transfix $ */

#include <cvc_compat.h>
#include <GL/glew.h>

namespace CVC_NAMESPACE
{
#if 0
  //functor for clipping via a bounding box
  class SceneLight
  {
  public:
    SceneLight()
      : position{0.0,0.0,0.0,1.0},
        diffuseColor{0.90f, 0.90f, 0.90f, 1.0f},
        specularColor{0.60f, 0.60f, 0.60f, 1.0f},
        ambientColor{0.0f,0.0f,0.0f,1.0f}
      {}

      void operator()();

    GLfloat position[4];
    GLfloat diffuseColor[4];
    GLfloat specularColor[4];
    GLfloat ambientColor[4];
  };
#endif
}
