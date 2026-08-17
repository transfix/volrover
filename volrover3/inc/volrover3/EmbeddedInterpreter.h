#ifndef VOLROVER3_EMBEDDEDINTERPRETER_H
#define VOLROVER3_EMBEDDEDINTERPRETER_H

// --------------------------------------------------------------------
// volrover3::EmbeddedInterpreter
// --------------------------------------------------------------------
// Owns the process's single embedded CPython interpreter (see
// docs/EMBEDDED_PYTHON.md). Boots CPython once, parks the GIL, and finalizes
// once in the destructor — never per script. One instance, owned by MainWindow,
// constructed right after the app + scene graph exist. Single-interpreter model
// (multi is unsafe for the C-extension packages pycvc/pycvc_gl/VTK).
//
// Modeled conceptually on verlihub's plugins/python lifecycle, ported to the
// Python 3 C-API: bring the interpreter up once, bracket every C++->Python
// crossing with the GIL, isolate script exceptions.
// --------------------------------------------------------------------

#include <cvc/core/app.h>

#include <cstdint>
#include <functional>
#include <memory>
#include <string>

class SceneGraph; // cvcGL (global namespace)

namespace volrover3 {

class PyHost;

// Single vs sub-interpreter policy (docs/EMBEDDED_PYTHON.md §12.2). The
// once-per-process CPython boot is identical in both modes; only the JOB
// LAUNCHER (Phase 4) branches on this. Single = the full-capability default
// (pycvc/pycvc_gl/VTK/numpy work); Multi = per-job Py_NewInterpreter isolation,
// pure-Python only (host/app bindings are unavailable + import-gated).
enum class InterpreterMode { Single, Multi };

class EmbeddedInterpreter {
public:
  // Startup config, sourced from volrover3::Settings and applied at boot
  // (restart-only — the interpreter topology can't be safely rebuilt mid-run).
  struct Config {
    InterpreterMode mode = InterpreterMode::Single;
    std::string python_home;          // "" -> $VOLROVER3_PYTHON_HOME / CPython default
    std::string module_path;          // dir with _vrhost.so + vrhost.py; "" -> $VOLROVER3_PYMODULE_PATH
    bool gate_multi_imports = true;   // deny pycvc/vtk/numpy in Multi sub-interpreters
  };

  // Boots the interpreter on the calling (UI) thread. `app` + `scene` are the
  // host's LIVE handles — used to build the PyHost that delivers the app to
  // scripts. `scene` may be null (headless/tests). `config` selects the mode +
  // python home; the 2-arg overload uses the defaults (single, no explicit home).
  EmbeddedInterpreter(std::shared_ptr<cvc::app> app, std::shared_ptr<SceneGraph> scene);
  EmbeddedInterpreter(std::shared_ptr<cvc::app> app, std::shared_ptr<SceneGraph> scene,
                      Config config);
  ~EmbeddedInterpreter();

  EmbeddedInterpreter(const EmbeddedInterpreter &) = delete;
  EmbeddedInterpreter &operator=(const EmbeddedInterpreter &) = delete;

  // True once CPython is initialized and ready.
  bool ok() const { return m_initialized; }

  // The active interpreter mode (drives the job launcher).
  InterpreterMode mode() const { return m_config.mode; }
  const Config &config() const { return m_config; }

  // Run a snippet in the interpreter's __main__ dict. Acquires the GIL around
  // the call; on a Python error prints the traceback and returns false (never
  // throws) so one bad script can't take down the host.
  bool run_string(const std::string &source);

  // REPL-style run for the console (docs/EMBEDDED_PYTHON.md §12.4): like
  // run_string, but under the GIL it swaps sys.stdout/sys.stderr to io.StringIO,
  // runs with Py_single_input so a bare expression echoes its repr (a real
  // REPL), captures both streams into `out`/`err`, and ALWAYS restores the real
  // streams (even on error). Returns false on a Python error (traceback in
  // `err`). run_string is left untouched.
  bool run_string_capture(const std::string &source, std::string &out, std::string &err);

  // The injected host facade (surfaced to scripts as `vrhost.host` once the
  // vrhost module is registered — Phase 1).
  std::shared_ptr<PyHost> host() const { return m_host; }

  // True once `import vrhost` succeeded AND the live PyHost was bound into it
  // (so `vrhost.host` / `vrhost.app()` deliver the running app). False in Multi
  // mode, when the vrhost module dir is unknown, or if the import failed.
  bool host_bound() const { return m_hostBound; }

  // Push the live QMainWindow's address to the vrhost shim so
  // `vrhost.main_window()` (shiboken6.wrapInstance) resolves it. Call after the
  // window exists; safe to call before boot (updates the C++ host only).
  void set_main_window_ptr(std::uintptr_t ptr);

  // Expose AppState::setWorldBounds to scripts as vrhost.set_world_bounds().
  //
  // `volrover3.world_bounds` documents itself as "computed from graphics bounds"
  // but is read-only in the state tree and is only ever written by the C++
  // generators — so a script that builds a scene through pycvc_gl leaves it at
  // the default and neither the grid nor the camera ever frames what it added
  // (pycvc.state_set on it raises cvc::read_only_error). Handing the setter to
  // Python closes that gap without making the node writable, and keeps the
  // script in charge of WHEN bounds are final — recomputing them on every
  // addGraphics would rebuild the grid once per node.
  //
  // `setter` is invoked on the caller's thread with the GIL held; MainWindow
  // supplies one that marshals onto the GUI thread. Pass nullptr to unbind.
  void set_world_bounds_hook(std::function<void(double, double, double, double, double, double)> setter);

private:
  // After boot (GIL held): put the vrhost module dir on sys.path, `import vrhost`
  // (which auto-imports pycvc, registering the shared cvc::app SWIG type), then
  // hand it the live PyHost via a PyCapsule so `vrhost.host` is the running app.
  // Best-effort: logs + returns false on failure, never throws. GIL must be held.
  bool bind_host();

  Config m_config;
  std::shared_ptr<PyHost> m_host;
  void *m_mainState = nullptr; // PyThreadState* parked after init
  bool m_initialized = false;
  bool m_booted = false;    // this instance called Py_Initialize (so it finalizes)
  bool m_hostBound = false; // vrhost imported + live host bound
  // Owns the world-bounds setter for as long as the injected Python callable can
  // be invoked; the callable holds a capsule pointing at it.
  std::shared_ptr<std::function<void(double, double, double, double, double, double)>>
      m_worldBoundsHook;
};

} // namespace volrover3

#endif // VOLROVER3_EMBEDDEDINTERPRETER_H
