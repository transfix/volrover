#include <volrover3/PyHost.h>

#include <cvc/gl/SceneGraph.h>

namespace volrover3 {

PyHost::PyHost(std::shared_ptr<cvc::app> app, std::shared_ptr<SceneGraph> scene)
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
