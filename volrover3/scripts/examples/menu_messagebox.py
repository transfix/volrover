# menu_messagebox.py — the "hello world" of the volrover3 Qt-from-Python bridge.
#
# Adds a "Demo" menu with one item to the running application's menu bar; clicking
# it pops up a message box. Every line here drives the app's REAL C++ QMainWindow.
#
# HOW TO RUN (inside a running volrover3): open the "Python Console" dock, switch
# to the REPL tab, and paste:
#
#     exec(open("scripts/examples/menu_messagebox.py").read())
#
# HOW IT WORKS: the embedded interpreter's `vrhost` module hands Python the live
# QMainWindow as a raw C++ pointer, which pyside6/shiboken6 adopt with
# wrapInstance. Because the hermetic cvcpkg pyside6/shiboken6 link the SAME Qt the
# C++ app does (one libQt6Core in the process), the wrapper IS the live window —
# so `.menuBar().addMenu(...)` mutates the actual UI. See docs/EMBEDDED_PYTHON.md §7.

import vrhost
from PySide6 import QtWidgets

# 1. Adopt the app's live main window (vrhost.main_window() == shiboken6.wrapInstance
#    of the host's QMainWindow*). Grab it once and reuse it.
window = vrhost.main_window()

# 2. Add a "Demo" menu and a "Say Hello" item to the real menu bar.
menu = window.menuBar().addMenu("&Demo")
action = menu.addAction("Say &Hello")


# 3. When the item is triggered, show a message box parented to the app window.
def on_hello():
    QtWidgets.QMessageBox.information(
        window,
        "volrover3",
        "Hello from embedded Python!\n\n"
        "This menu item and dialog were created by a Python script "
        "driving the live C++ QMainWindow.",
    )


action.triggered.connect(on_hello)

# Keep the Python objects (menu, action, and the connected closure) alive beyond
# this script's scope — Qt owns the C++ side, but the signal→slot connection needs
# the Python callable to stay referenced.
vrhost._demo_menu = menu
vrhost._demo_action = action
vrhost._demo_on_hello = on_hello

print("Added the 'Demo ▸ Say Hello' menu item — click it to see the message box.")
