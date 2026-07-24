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
// The generated vrhostPYTHON_wrap.cxx is compiled INTO volrover3_lib and
// registered via PyImport_AppendInittab("vrhost", PyInit_vrhost) -- no .so to
// find on disk. This is a baked-in first-class module, not a plugin.
// ---------------------------------------------------------------------------
%module(directors="1") vrhost

%{
#include <volrover3/PyHost.h>
using volrover3::PyHost;
%}

// Reuse pycvc's %shared_ptr(cvc::app) registration + exception translation so
// the app handle is one shared type across both modules.
%import "pycvc.i"

%include <std_string.i>
%include <std_vector.i>
%include <stdint.i>
%template(StringVector) std::vector<std::string>;

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
