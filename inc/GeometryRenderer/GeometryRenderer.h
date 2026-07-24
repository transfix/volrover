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

/* $Id: GeometryRenderer.h 5692 2012-06-01 22:39:25Z transfix $ */

#include <cvc/state_object.h>
#include <cvc_compat.h>
#include <cvc_compat.h>
#include <boost/function.hpp>

#include <cvcraw_geometry/cvcgeom.h>

namespace qglviewer
{
  class Camera;
}

namespace CVC_NAMESPACE
{
  class GeometryRenderer : public state_object<GeometryRenderer>
  {
  public:
    typedef boost::function<void ()> Callback;

    GeometryRenderer(const qglviewer::Camera* vc = NULL)
      : state_object<GeometryRenderer>(volrover_app_instance()),
        _viewerCamera(vc)
    { 
      defaultConstructor(); 
    }
    
    //Render's all the geometry found starting at sceneRoot
    virtual void render(const std::string& sceneRoot);

    const qglviewer::Camera *camera() const { return _viewerCamera; }
    void camera(const qglviewer::Camera* cam) { _viewerCamera = cam; }

  protected:
    void defaultConstructor();

    virtual void renderState(const std::string& root, const std::string& s);
    virtual void handleStateChanged(const std::string& childState);

    void doDrawBoundingBox(const std::string& s);
    void doDrawBoundingBox(const BoundingBox& bbox);

    const qglviewer::Camera *_viewerCamera;
  };
}
