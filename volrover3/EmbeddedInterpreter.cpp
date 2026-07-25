// Python.h must be included before any system headers.
// clang-format off
#include <Python.h>
// clang-format on

#include <volrover3/EmbeddedInterpreter.h>
#include <volrover3/PyHost.h>

#include <cvc/gl/SceneGraph.h> // complete type so the scene-capsule dtor can run ~SceneGraph

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

// PyCapsule("cvc.app") destructor: free the heap shared_ptr copy bind_host made.
static void app_capsule_dtor(PyObject *cap) {
  delete static_cast<std::shared_ptr<cvc::app> *>(PyCapsule_GetPointer(cap, "cvc.app"));
}

// PyCapsule("cvc.scenegraph") destructor: free the heap shared_ptr copy. The host
// keeps its own ref, so this drop is not the last one during normal operation.
static void scene_capsule_dtor(PyObject *cap) {
  delete static_cast<std::shared_ptr<SceneGraph> *>(PyCapsule_GetPointer(cap, "cvc.scenegraph"));
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

  // Import the PURE-PYTHON vrhost shim (it `import pycvc`s). No SWIG module.
  PyObject *vr = PyImport_ImportModule("vrhost"); // new
  if (!vr) {
    PyErr_Print(); // report + clear
    std::cerr << "volrover3: `import vrhost` failed — app scripting unavailable "
                 "(set VOLROVER3_PYMODULE_PATH to the dir containing vrhost.py)\n";
    return false;
  }

  bool bound = false;
  // Deliver the LIVE app as a PyCapsule("cvc.app") holding a heap shared_ptr COPY
  // (shares ownership; the capsule's destructor frees the copy). vrhost.app()
  // feeds it to pycvc.app_from_capsule, which wraps it into pycvc's OWN app type —
  // so volrover3 needs no SWIG and no pycvc.i (see docs/EMBEDDED_PYTHON.md §5).
  auto *sp = new std::shared_ptr<cvc::app>(m_host->app());
  PyObject *cap = PyCapsule_New(sp, "cvc.app", &app_capsule_dtor); // new
  if (cap) {
    if (PyObject_SetAttrString(vr, "_app_capsule", cap) == 0) {
      // Seed the window pointer (usually 0 at boot; set_main_window_ptr updates it
      // once MainWindow exists).
      if (PyObject *mw = PyLong_FromUnsignedLongLong(m_host->main_window_ptr())) {
        PyObject_SetAttrString(vr, "_main_window_ptr", mw);
        Py_DECREF(mw);
      }
      bound = true;
    } else {
      PyErr_Print();
    }
    Py_DECREF(cap); // vrhost._app_capsule holds the ref; dtor frees sp on teardown
  } else {
    delete sp; // capsule creation failed — free the copy ourselves
  }

  // Deliver the LIVE SceneGraph as PyCapsule("cvc.scenegraph") so vrhost.scene()
  // can adopt it into a pycvc_gl.Scene — add_geometry/add_volume then mutate the
  // RUNNING scene and appear in the live window (VTKRenderWidget's timer drains
  // the queued mutations). Same raw-shared_ptr-copy pattern as the app capsule; a
  // null scene (headless/tests) simply leaves vrhost._scene_capsule = None.
  if (auto scene = m_host->scene()) {
    auto *ssp = new std::shared_ptr<SceneGraph>(scene);
    PyObject *scap = PyCapsule_New(ssp, "cvc.scenegraph", &scene_capsule_dtor); // new
    if (scap) {
      if (PyObject_SetAttrString(vr, "_scene_capsule", scap) != 0)
        PyErr_Print();
      Py_DECREF(scap); // vrhost._scene_capsule holds the ref
    } else {
      delete ssp;
    }
  }

  Py_DECREF(vr);
  return bound;
}

// Push the live QMainWindow address to the (already-imported) vrhost shim so
// vrhost.main_window() resolves it. Called by MainWindow after boot (GIL parked).
void EmbeddedInterpreter::set_main_window_ptr(std::uintptr_t ptr) {
  if (m_host)
    m_host->set_main_window_ptr(ptr);
  if (!m_hostBound)
    return;
  PyGILState_STATE gil = PyGILState_Ensure();
  if (PyObject *vr = PyImport_ImportModule("vrhost")) { // cached import
    if (PyObject *mw = PyLong_FromUnsignedLongLong(ptr)) {
      PyObject_SetAttrString(vr, "_main_window_ptr", mw);
      Py_DECREF(mw);
    }
    Py_DECREF(vr);
  } else {
    PyErr_Clear();
  }
  PyGILState_Release(gil);
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
