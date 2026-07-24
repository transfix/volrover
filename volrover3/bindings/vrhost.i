// ---------------------------------------------------------------------------
// vrhost.i -- the volrover3 host-control SWIG module.
//
// `import vrhost` gives an embedded-interpreter script the injected `vrhost.host`
// proxy (a live volrover3::PyHost bound to the running MainWindow) and, through
// it, the live cvc::app. See docs/EMBEDDED_PYTHON.md.
//
// CRITICAL: %import "pycvc.i". pycvc (libcvc #136) registers
// %shared_ptr(cvc::app); SWIG keys swig_type_info by mangled type name in a
// per-process runtime table, so importing pycvc.i makes the
// std::shared_ptr<cvc::app> that PyHost::app() returns the IDENTICAL proxy type
// pycvc consumes -- the handle round-trips into pycvc.state_set / volume / sdf
// with no conversion. WITHOUT this %import, vrhost would mint a distinct,
// incompatible app type and pycvc.state_set(vrhost.host.app()) would be
// rejected (a silent, confusing failure). The two modules must also be built
// against the SAME libcvc (cvc::cvc) and SWIG version.
//
// Build note: SWIG needs pycvc.i on its -I path at generation time (the cvcpkg
// pycvc package must ship its .i, or volrover3 references the libcvc source).
// vrhost is built as a bundled `_vrhost.so` + `vrhost.py` shipped inside the
// volrover3 install; EmbeddedInterpreter puts that dir on sys.path and imports
// it at boot (after pycvc), then hands it the live host via a PyCapsule (below).
// It is a first-class, always-present environment feature -- not a user-loaded
// plugin -- even though it rides on a private .so the way pycvc does.
// ---------------------------------------------------------------------------
%module(directors="1") vrhost

// CRITICAL import ORDER: pycvc must initialize BEFORE the _vrhost C-extension.
// SWIG links cross-module type tables in import order, so _pycvc has to register
// std::shared_ptr<cvc::app> FIRST -- otherwise _vrhost (which SWIG imports at the
// top of vrhost.py, before the `import pycvc` its %import emits lower down) mints a
// private, unshared copy of that type and pycvc.state_set(vrhost.host.app())
// TypeErrors. %pythonbegin runs at the very top of vrhost.py, ahead of _vrhost.
%pythonbegin %{
import pycvc as _pycvc_preload  # noqa: F401  (force _pycvc to register its types first)
%}

%{
// cvc::exception must be complete here: pycvc's %exception typemap (inherited via
// the %import below) wraps every vrhost function in `catch (const cvc::exception&)`,
// but %import does NOT copy pycvc's own header block -- so we include it ourselves.
#include <cvc/core/exception.h>
#include <volrover3/PyHost.h>
using volrover3::PyHost;
%}

// %shared_ptr's TYPEMAPS are per-module. %import shares pycvc's cvc::app type
// REGISTRATION (same mangled name -> unified in the per-process SWIG runtime), but
// we must ALSO re-declare %shared_ptr(cvc::app) HERE so PyHost::app()'s
// std::shared_ptr<cvc::app> return actually uses the shared_ptr typemap. Without
// this, vrhost mints a distinct shared_ptr<cvc::app> proxy ("no destructor found")
// that pycvc.state_set rejects with a TypeError. std_shared_ptr.i must precede it.
%include <std_shared_ptr.i>
%import "pycvc.i"
%shared_ptr(cvc::app)

// pycvc's %exception typemap is INHERITED via %import, but the SWIG_exception
// macro it expands to is only emitted into a module that %includes exception.i
// DIRECTLY -- %import brings the handler, not the supporting runtime macro. So
// re-declare a self-contained handler for vrhost's own wrappers. SWIG_exception_fail
// is always emitted; cvc::exception is complete via the #include in the header block.
%exception {
  try {
    $action
  } catch (const cvc::exception &e) {
    SWIG_exception_fail(SWIG_RuntimeError, e.what());
  } catch (const std::exception &e) {
    SWIG_exception_fail(SWIG_RuntimeError, e.what());
  }
}

%include <std_string.i>
%include <stdint.i>
// std::vector<std::string> (StringVector) is already instantiated by pycvc.i and
// inherited via %import -- do not re-instantiate it here (SWIG Warning 404).

// PyHost is delivered as an injected proxy (vrhost.host); scripts never build
// one. Hide the ctor; expose the control surface.
%nodefaultctor volrover3::PyHost;
%nodefaultdtor volrover3::PyHost;

%include <volrover3/PyHost.h>

%inline %{
namespace volrover3 {
// Set by EmbeddedInterpreter at boot to the live host, then surfaced to Python
// as the module attribute `vrhost.host` (bound via SWIG_NewPointerObj). This
// accessor is the fallback the .py shim reads; it is NOT a singleton owner --
// the interpreter owns the PyHost and injects it.
static PyHost *g_current_host = nullptr;
PyHost *_current_host() { return g_current_host; }
void _set_current_host(PyHost *h) { g_current_host = h; }
} // namespace volrover3

// Called by the embedding host (EmbeddedInterpreter) at boot: it passes the live
// PyHost* wrapped in a PyCapsule, so no wrapper-internal SWIG type is needed from
// the host TU (vrhost may be a separate _vrhost.so). PyObject* params/returns
// pass through SWIG untouched.
static PyObject *_bind_host_capsule(PyObject *cap) {
  void *p = PyCapsule_GetPointer(cap, "volrover3.PyHost");
  if (!p)
    return NULL;
  volrover3::g_current_host = static_cast<volrover3::PyHost *>(p);
  Py_RETURN_NONE;
}
%}

// Convenience Python surface: `vrhost.host` (the injected live host) and
// `vrhost.app()` (== host.app()), so scripts read naturally:
//     import vrhost, pycvc
//     app = vrhost.app()
//     pycvc.state_set(app, "volrover3.camera.fov", "42")
%pythoncode %{
host = _current_host()

def app():
    """The live volrover3 cvc::app (shared_ptr<cvc::app>) for pycvc scripts."""
    return _current_host().app()

def main_window():
    """The running QMainWindow as a PySide6 object (Phase 6; needs pyside6/shiboken6).

    Bridges the SWIG world to Shiboken by adopting the existing C++ QWidget* via
    its address -- the same raw-pointer-as-int pattern as the pycvc<->VTK bridge.
    """
    import shiboken6
    from PySide6 import QtWidgets
    return shiboken6.wrapInstance(_current_host().main_window_ptr(),
                                  QtWidgets.QMainWindow)
%}
