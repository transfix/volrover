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

  // Park the GIL so UI/worker threads acquire it per-crossing via
  // PyGILState_Ensure (the Py3 analogue of verlihub's PyThreadState_Swap(NULL) +
  // release-lock). Every C++->Python entry brackets Ensure/Release.
  m_mainState = PyEval_SaveThread();
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

} // namespace volrover3
