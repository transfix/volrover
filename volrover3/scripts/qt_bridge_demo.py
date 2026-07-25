# qt_bridge_demo.py — drive the live volrover3 QMainWindow from embedded Python.
#
# The embedded interpreter hands Python the running app's QMainWindow* (a raw
# address) via `vrhost.main_window()`, which adopts it with
# shiboken6.wrapInstance(addr, QtWidgets.QMainWindow). Because the hermetic cvcpkg
# pyside6/shiboken6 link the SAME cvcpkg Qt the C++ app does (one libQt6Core in
# the process — verified by QtWidgets.QApplication.instance() being non-None
# without Python ever constructing one), the returned object IS the live window:
# reading its widgets and mutating its menus changes the real UI.
#
# Run it from the volrover3 Python console (PyConsoleDock REPL):
#     import qt_bridge_demo; qt_bridge_demo.install()
# See docs/EMBEDDED_PYTHON.md §7 (the Qt bridge).

import vrhost
from PySide6 import QtCore, QtWidgets


def read_widgets(window=None):
    """READ side: report the live window's title, menus and dock widgets."""
    mw = window or vrhost.main_window()
    menus = [m.title().replace("&", "") for m in mw.menuBar().findChildren(QtWidgets.QMenu)]
    docks = [d.objectName() or d.windowTitle() for d in mw.findChildren(QtWidgets.QDockWidget)]
    info = {"title": mw.windowTitle(), "menus": menus, "docks": docks}
    print("volrover3 window:", info)
    return info


def install(window=None):
    """MUTATE side: add a top-level 'Python' menu whose action fires a QMessageBox.

    Returns the created QAction. Idempotent-ish: re-running adds another menu.
    """
    mw = window or vrhost.main_window()

    menu = mw.menuBar().addMenu("&Python")
    menu.setObjectName("PythonMenu")
    act = menu.addAction("Say &Hello from Python")
    act.setObjectName("PythonHelloAction")

    def _hello():
        box = QtWidgets.QMessageBox(mw)
        box.setObjectName("PythonHelloBox")
        box.setWindowTitle("volrover3")
        box.setIcon(QtWidgets.QMessageBox.Information)
        box.setText("Hello from embedded Python!")
        box.setInformativeText(
            "This menu item and dialog were created by a Python script driving "
            "the live C++ QMainWindow through pyside6 + shiboken6.")
        # Headless/CI-safe: self-dismiss on the next event-loop turn so exec()
        # never blocks a test. A real user sees the dialog and clicks OK.
        QtCore.QTimer.singleShot(0, box.accept)
        box.exec()
        vrhost._bridge_hello_fired = True  # observable side effect for the test

    act.triggered.connect(_hello)

    # Keep Python refs alive on the module so the menu/action/closure aren't GC'd
    # (Qt holds the C++ side; these keep the Python wrappers + connection alive).
    vrhost._bridge_menu = menu
    vrhost._bridge_action = act
    vrhost._bridge_hello = _hello
    vrhost._bridge_hello_fired = False
    print("installed 'Python ▸ Say Hello from Python' menu item on the live window")
    return act


if __name__ == "__main__":
    read_widgets()
    install()
