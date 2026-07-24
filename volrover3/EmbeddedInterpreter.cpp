// Python.h must be included before any system headers.
// clang-format off
#include <Python.h>
// clang-format on

#include <volrover3/EmbeddedInterpreter.h>
#include <volrover3/PyHost.h>

#include <cstdlib>
#include <iostream>

namespace volrover3 {

EmbeddedInterpreter::EmbeddedInterpreter(std::shared_ptr<cvc::app> app,
                                         std::shared_ptr<SceneGraph> scene)
    : EmbeddedInterpreter(std::move(app), std::move(scene), Config{}) {}

EmbeddedInterpreter::EmbeddedInterpreter(std::shared_ptr<cvc::app> app,
                                         std::shared_ptr<SceneGraph> scene, Config config)
    : m_config(std::move(config)),
      m_host(std::make_shared<PyHost>(std::move(app), std::move(scene))) {
  if (Py_IsInitialized()) {
    // Single-interpreter model: CPython is already up (another
    // EmbeddedInterpreter, or an embedding host). Reuse it; do not finalize it.
    m_initialized = true;
    m_booted = false;
    // The GIL is parked by whoever booted; take it to (re)bind our host.
    PyGILState_STATE gil = PyGILState_Ensure();
    m_hostBound = bind_host();
    PyGILState_Release(gil);
    return;
  }

  PyConfig pyconfig;
  PyConfig_InitPythonConfig(&pyconfig);

  // Point the interpreter at the hermetic cvcpkg Python home so it finds its
  // stdlib regardless of a possibly-non-relocatable compiled-in prefix. Prefer
  // the configured home (from Settings), else $VOLROVER3_PYTHON_HOME, else
  // CPython's own defaults.
  std::string home = m_config.python_home;
  if (home.empty()) {
    if (const char *env = std::getenv("VOLROVER3_PYTHON_HOME"))
      home = env;
  }
  if (!home.empty())
    PyConfig_SetBytesString(&pyconfig, &pyconfig.home, home.c_str());

  PyStatus status = Py_InitializeFromConfig(&pyconfig);
  PyConfig_Clear(&pyconfig);
  if (PyStatus_Exception(status)) {
    std::cerr << "volrover3: embedded CPython init failed"
              << (status.err_msg ? std::string(": ") + status.err_msg : "") << std::endl;
    m_initialized = false;
    return;
  }
  m_initialized = true;
  m_booted = true;

  // Deliver the live cvc::app to scripts: `import vrhost` + capsule-bind the host
  // while the GIL is still held from init (before we park it below).
  m_hostBound = bind_host();

  // Park the GIL so UI/worker threads acquire it per-crossing via
  // PyGILState_Ensure (the Py3 analogue of verlihub's PyThreadState_Swap(NULL) +
  // release-lock). Every C++->Python entry brackets Ensure/Release.
  m_mainState = PyEval_SaveThread();
}

// GIL must be held. Best-effort: never throws; logs + returns false on failure.
bool EmbeddedInterpreter::bind_host() {
  // Resolve the dir holding _vrhost.so + vrhost.py (Config, else env). pycvc
  // itself rides on the python-home site-packages, so only this dir is added.
  std::string dir = m_config.module_path;
  if (dir.empty()) {
    if (const char *env = std::getenv("VOLROVER3_PYMODULE_PATH"))
      dir = env;
  }
  if (!dir.empty()) {
    if (PyObject *sysPath = PySys_GetObject("path")) { // borrowed
      PyObject *p = PyUnicode_FromString(dir.c_str());
      if (p) {
        PyList_Insert(sysPath, 0, p); // prepend so our module wins
        Py_DECREF(p);
      }
    }
  }

  // `import vrhost` also imports pycvc (SWIG emits it), registering the shared
  // std::shared_ptr<cvc::app> SWIG type so PyHost::app() round-trips into pycvc.
  PyObject *vr = PyImport_ImportModule("vrhost"); // new
  if (!vr) {
    PyErr_Print(); // report + clear
    std::cerr << "volrover3: `import vrhost` failed — app scripting unavailable "
                 "(set VOLROVER3_PYMODULE_PATH to the dir containing _vrhost.so)\n";
    return false;
  }

  bool bound = false;
  // Hand the live PyHost* across as a name-matched PyCapsule (see vrhost.i).
  PyObject *cap = PyCapsule_New(m_host.get(), "volrover3.PyHost", nullptr); // new
  if (cap) {
    if (PyObject *r = PyObject_CallMethod(vr, "_bind_host_capsule", "O", cap)) { // new
      Py_DECREF(r);
      // `host` was captured as None at import (before this bind) — refresh it so
      // `vrhost.host` is the live host (vrhost.app() is already dynamic).
      if (PyObject *live = PyObject_CallMethod(vr, "_current_host", nullptr)) { // new
        PyObject_SetAttrString(vr, "host", live);
        Py_DECREF(live);
        bound = true;
      }
    } else {
      PyErr_Print();
    }
    Py_DECREF(cap);
  }
  Py_DECREF(vr);
  return bound;
}

EmbeddedInterpreter::~EmbeddedInterpreter() {
  if (!m_booted)
    return; // we did not boot CPython, so we must not finalize it
  if (m_mainState) {
    PyEval_RestoreThread(static_cast<PyThreadState *>(m_mainState));
    m_mainState = nullptr;
  }
  Py_FinalizeEx();
}

bool EmbeddedInterpreter::run_string(const std::string &source) {
  if (!m_initialized)
    return false;
  PyGILState_STATE gil = PyGILState_Ensure();
  int rc = PyRun_SimpleString(source.c_str());
  if (rc != 0 && PyErr_Occurred())
    PyErr_Print(); // report + clear; one bad script must not kill the host
  PyGILState_Release(gil);
  return rc == 0;
}

bool EmbeddedInterpreter::run_string_capture(const std::string &source, std::string &out,
                                             std::string &err) {
  out.clear();
  err.clear();
  if (!m_initialized)
    return false;

  PyGILState_STATE gil = PyGILState_Ensure();
  bool ok = false;

  PyObject *mainMod = PyImport_AddModule("__main__"); // borrowed
  PyObject *sysMod = PyImport_ImportModule("sys");    // new
  PyObject *ioMod = PyImport_ImportModule("io");      // new

  if (mainMod && sysMod && ioMod) {
    PyObject *mainDict = PyModule_GetDict(mainMod);              // borrowed
    PyObject *oldOut = PyObject_GetAttrString(sysMod, "stdout"); // new (or null)
    PyObject *oldErr = PyObject_GetAttrString(sysMod, "stderr"); // new (or null)
    PyObject *capOut = PyObject_CallMethod(ioMod, "StringIO", nullptr); // new
    PyObject *capErr = PyObject_CallMethod(ioMod, "StringIO", nullptr); // new

    if (capOut && capErr) {
      PyObject_SetAttrString(sysMod, "stdout", capOut);
      PyObject_SetAttrString(sysMod, "stderr", capErr);

      // Py_single_input: a bare expression echoes its repr via displayhook ->
      // our captured stdout, exactly like an interactive REPL.
      PyObject *res = PyRun_String(source.c_str(), Py_single_input, mainDict, mainDict);
      if (res) {
        Py_DECREF(res);
        ok = true;
      } else {
        PyErr_Print(); // -> capErr, and clears the error
      }

      auto drain = [](PyObject *sio) -> std::string {
        std::string s;
        if (PyObject *v = PyObject_CallMethod(sio, "getvalue", nullptr)) {
          if (const char *c = PyUnicode_AsUTF8(v))
            s = c;
          Py_DECREF(v);
        }
        return s;
      };
      out = drain(capOut);
      err = drain(capErr);

      // ALWAYS restore the real streams (even on error).
      if (oldOut)
        PyObject_SetAttrString(sysMod, "stdout", oldOut);
      if (oldErr)
        PyObject_SetAttrString(sysMod, "stderr", oldErr);
    }
    Py_XDECREF(oldOut);
    Py_XDECREF(oldErr);
    Py_XDECREF(capOut);
    Py_XDECREF(capErr);
  }
  Py_XDECREF(sysMod);
  Py_XDECREF(ioMod);

  PyGILState_Release(gil);
  return ok;
}

} // namespace volrover3
