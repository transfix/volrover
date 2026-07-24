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

  // The injected host facade (surfaced to scripts as `vrhost.host` once the
  // vrhost module is registered — Phase 1).
  std::shared_ptr<PyHost> host() const { return m_host; }

private:
  Config m_config;
  std::shared_ptr<PyHost> m_host;
  void *m_mainState = nullptr; // PyThreadState* parked after init
  bool m_initialized = false;
  bool m_booted = false; // this instance called Py_Initialize (so it finalizes)
};

} // namespace volrover3

#endif // VOLROVER3_EMBEDDEDINTERPRETER_H
