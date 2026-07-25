#ifndef VOLROVER3_PYHOST_H
#define VOLROVER3_PYHOST_H

// --------------------------------------------------------------------
// volrover3::PyHost
// --------------------------------------------------------------------
// The host-control facade the embedded Python interpreter drives (SWIG module
// `vrhost`, see bindings/vrhost.i). Two jobs:
//
//   1. DELIVER the live cvc::app to pycvc scripts. pycvc (libcvc #136) has no
//      app singleton — every op takes an explicit std::shared_ptr<cvc::app>.
//      PyHost::app() hands scripts the ONE app volrover3 owns, so writes land
//      on the RUNNING state tree the UI widgets and SceneNodes already watch,
//      not a fresh disconnected pycvc.make_app() tree.
//
//   2. Control host state / scene / (later) the Qt UI through a typed surface.
//
// PyHost is NOT a singleton. It is constructed with injected handles (the
// owned app + scene graph) and bound into the interpreter as the `vrhost.host`
// proxy at boot -- the injected-handle pattern (contrast verlihub's
// cpiPython::me global + vh.myid reflection). See docs/EMBEDDED_PYTHON.md.
// --------------------------------------------------------------------

#include <cvc/core/app.h>

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

// cvcGL's scene graph (global namespace), the one pycvc_gl wraps and that
// carries the #136 injected-app ctor SceneGraph(cvc::app&, prefix). volrover3
// is converging onto cvcGL (dropping its private SceneGraph fork), so the C++
// app and Python scripts share one scene type over one app/state tree.
class SceneGraph;

namespace volrover3 {

class PyHost {
public:
  // Injected at construction by EmbeddedInterpreter (owned by MainWindow) --
  // the same app + scene graph the running application uses. `app` is the one
  // volrover3-owned shared_ptr<cvc::app> (no singleton); `scene` may be null in
  // headless/unit-test contexts.
  PyHost(std::shared_ptr<cvc::app> app, std::shared_ptr<SceneGraph> scene);
  ~PyHost();

  PyHost(const PyHost &) = delete;
  PyHost &operator=(const PyHost &) = delete;

  // -- app delivery (the crux) --------------------------------------------
  // THE live volrover3 app, as the same shared_ptr<cvc::app> proxy type pycvc
  // consumes (vrhost.i %import "pycvc.i"). Round-trips into pycvc.state_set /
  // volume / sdf with zero conversion:
  //     app = vrhost.host.app(); pycvc.state_set(app, "volrover3.camera.fov", "42")
  std::shared_ptr<cvc::app> app() const { return m_app; }

  // The state subtree prefix volrover3 owns, so scripts need not hardcode it.
  std::string state_prefix() const { return "volrover3"; }

  // NOTE: there are deliberately NO typed state/scene helpers here. State is
  // already SWIG-wrapped by pycvc (state_set/state_get/…) and the scene graph
  // by pycvc_gl — both operate on the app() this delivers. A script writes
  // `pycvc.state_set(vrhost.host.app(), "volrover3.camera.fov", "42")`, not a
  // bespoke PyHost method. PyHost only exposes what those bindings CANNOT get
  // for themselves: the live app handle, the Qt main window, and host-loop
  // actions.

  // -- scene / host-loop ---------------------------------------------------
  // The host's live SceneGraph (so a script can bridge it into pycvc_gl rather
  // than spin up a parallel scene).
  std::shared_ptr<SceneGraph> scene() const { return m_scene; }
  // Enqueue a render on the UI thread (via SceneGraph::postEvent, drained by
  // VTKRenderWidget's QTimer). A host-loop action pycvc_gl does not cover;
  // safe to call from Python worker threads.
  void request_render();

  // -- Qt bridge (Phase 6; see docs/EMBEDDED_PYTHON.md §7) -----------------
  // The live QMainWindow* as an integer, for PySide6/Shiboken:
  //     mw = shiboken6.wrapInstance(vrhost.host.main_window_ptr(),
  //                                 QtWidgets.QMainWindow)
  // Scripts then manipulate / replace the running UI with PySide6 widgets in
  // the one C++ QApplication. Returns 0 until a main window is bound.
  std::uintptr_t main_window_ptr() const { return m_mainWindow; }
  void set_main_window_ptr(std::uintptr_t mw) { m_mainWindow = mw; }

private:
  std::shared_ptr<cvc::app> m_app;
  std::shared_ptr<SceneGraph> m_scene;
  std::uintptr_t m_mainWindow = 0;
};

} // namespace volrover3

#endif // VOLROVER3_PYHOST_H
