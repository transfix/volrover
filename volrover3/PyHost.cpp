/*
  Copyright 2026 CyberPC Angel, LLC

  This file is part of VolumeRover.

  VolumeRover is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  VolumeRover is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.
*/

#include <volrover3/PyHost.h>

#include <cvc/gl/SceneGraph.h>

namespace volrover3 {

PyHost::PyHost(std::shared_ptr<cvc::app> app, std::shared_ptr<cvc::gl::SceneGraph> scene)
    : m_app(std::move(app)), m_scene(std::move(scene)) {}

PyHost::~PyHost() = default;

void PyHost::request_render() {
  if (!m_scene)
    return;
  // Marshal onto the scene's owner (UI) thread; safe from any thread.
  auto sg = m_scene;
  sg->postEvent([sg]() { sg->update(); });
}

} // namespace volrover3
