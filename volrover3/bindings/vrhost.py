# vrhost — volrover3 host-control shim (PURE PYTHON; no SWIG, no pycvc.i).
#
# The embedded C++ host (EmbeddedInterpreter) injects the running app as a
# PyCapsule and the QMainWindow's address, then pycvc/shiboken6 turn those into
# usable objects. There is NO SWIG module here and NO %import of pycvc.i:
#
#   * app delivery — the host hands over `_app_capsule` (a PyCapsule named
#     "cvc.app" holding the C++ shared_ptr<cvc::app>). `pycvc.app_from_capsule`
#     wraps it into pycvc's OWN app proxy, so the handle is type-compatible with
#     pycvc.state_set / volume / geometry BY CONSTRUCTION — no cross-module SWIG
#     type sharing, no SWIG-runtime-version coupling.
#   * Qt bridge — the host sets `_main_window_ptr` (the QMainWindow* as an int);
#     shiboken6.wrapInstance adopts it (an int is all that crosses).
#
# See docs/EMBEDDED_PYTHON.md §5 (app delivery) and §7 (the Qt bridge).

import pycvc

# ── injected by EmbeddedInterpreter at boot / when the main window is set ──
_app_capsule = None    # PyCapsule("cvc.app") -> shared_ptr<cvc::app>*
_main_window_ptr = 0   # QMainWindow* as an int


def app():
    """The live volrover3 cvc::app (a pycvc app handle) for pycvc scripts.

        import vrhost, pycvc
        pycvc.state_set(vrhost.app(), "volrover3.camera.fov", "42")
    """
    if _app_capsule is None:
        raise RuntimeError(
            "vrhost: no host app bound — are you running inside the volrover3 "
            "embedded interpreter?")
    return pycvc.app_from_capsule(_app_capsule)


def main_window_ptr():
    """The running QMainWindow's address (int), or 0 if the window isn't up yet."""
    return _main_window_ptr


def main_window():
    """The running QMainWindow as a PySide6 object (needs pyside6/shiboken6).

    Adopts the existing C++ QWidget* by its address — the same raw-pointer-as-int
    pattern the pycvc<->VTK bridge uses. Requires the hermetic cvcpkg
    pyside6/shiboken6 linking the SAME Qt the app does (one libQt6Core in-process).
    """
    if not _main_window_ptr:
        raise RuntimeError("vrhost: main window not set yet")
    import shiboken6
    from PySide6 import QtWidgets
    return shiboken6.wrapInstance(_main_window_ptr, QtWidgets.QMainWindow)
