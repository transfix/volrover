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
_scene_capsule = None  # PyCapsule("cvc.scenegraph") -> shared_ptr<SceneGraph>*
_main_window_ptr = 0   # QMainWindow* as an int
_set_world_bounds = None  # PyCFunction -> AppState::setWorldBounds (see below)


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


def scene():
    """The LIVE volrover3 3D scene as a ``pycvc_gl.Scene`` bound to the running window.

    ``add_geometry`` / ``add_volume`` on the returned Scene mutate the RUNNING
    scene graph, so they appear in the live volrover3 viewport (the render widget's
    timer drains the queued mutations) — not a separate offscreen scene:

        import vrhost
        s = vrhost.scene()
        s.add_geometry("mesh", g)   # shows up in the live window

    The host hands the scene across as ``_scene_capsule`` (a PyCapsule named
    "cvc.scenegraph"); ``pycvc_gl.scene_from_capsule`` wraps it + the app capsule
    into a Scene that ADOPTS the live SceneGraph — no SWIG type sharing, no
    parallel scene. Requires ``pycvc_gl`` and single-interpreter mode (multi-mode
    denies pycvc_gl imports).
    """
    if _scene_capsule is None or _app_capsule is None:
        raise RuntimeError(
            "vrhost: no host scene bound — are you running inside the volrover3 "
            "embedded interpreter in single-interpreter mode?")
    import pycvc_gl
    return pycvc_gl.scene_from_capsule(_app_capsule, _scene_capsule)


def set_world_bounds(minx, miny, minz, maxx, maxy, maxz):
    """Set volrover3's world bounds — the grid and camera framing follow them.

    ``volrover3.world_bounds`` is read-only in the state tree (``pycvc.state_set``
    on it raises ``cvc::read_only_error``) and the C++ generators are the only
    things that ever wrote it, so a scene built from a script left the grid and
    camera framing a default 1-unit box. This is ``AppState::setWorldBounds``,
    which is what those generators call::

        import vrhost
        s = vrhost.scene()
        ...                                  # add all your graphics first
        vrhost.set_world_bounds(*s.compute_graphics_bounds())

    Call it ONCE, after the scene is built: each call rebuilds the world grid and
    re-centres the camera's orbit, so calling it per-node is needless work. It
    does not move the camera itself — the user stays in control of the view.
    """
    if _set_world_bounds is None:
        raise RuntimeError(
            "vrhost: no world-bounds setter bound — are you running inside the "
            "volrover3 embedded interpreter?")
    _set_world_bounds(float(minx), float(miny), float(minz),
                      float(maxx), float(maxy), float(maxz))


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
