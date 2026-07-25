// Phase 6: the Qt-from-Python bridge. Proves embedded Python can READ the live
// app's QMainWindow widgets and MUTATE its UI — add a menu item that, when
// triggered, fires a QMessageBox — all through pyside6 + shiboken6 wrapping the
// real C++ QWidget* the host hands over (vrhost.main_window()).
//
// The load-bearing precondition is a SINGLE cvcpkg-Qt in the process: the
// hermetic pyside6/shiboken6 link the same libQt6Core the C++ app does (no
// bundled Qt), so `QtWidgets.QApplication.instance()` is non-None from Python
// without Python ever constructing one, and shiboken6.wrapInstance(addr, ...)
// yields the live window rather than a foreign-Qt proxy. See EMBEDDED_PYTHON.md §7.
//
// One interpreter for the suite (CPython can't re-Py_Initialize post-finalize).
// PySide6 + the vrhost module dir come from the environment set by CTest.

#include <cvc/core/app.h>
#include <gtest/gtest.h>
#include <volrover3/EmbeddedInterpreter.h>
#include <volrover3/PyHost.h>

#include <QAction>
#include <QApplication>
#include <QMainWindow>
#include <QMenu>
#include <QMenuBar>

#include <cstdint>
#include <cstdlib>
#include <memory>
#include <string>

class VrHostQtBridgeTest : public ::testing::Test {
protected:
  static std::shared_ptr<cvc::app> app;
  static std::unique_ptr<volrover3::EmbeddedInterpreter> interp;
  static void SetUpTestSuite() {
    app = std::make_shared<cvc::app>();
    interp = std::make_unique<volrover3::EmbeddedInterpreter>(app, nullptr);
    ASSERT_TRUE(interp->ok());
    ASSERT_TRUE(interp->host_bound()) << "vrhost not bound — is VOLROVER3_PYMODULE_PATH set?";
  }
  static void TearDownTestSuite() {
    interp.reset();
    app.reset();
  }
};
std::shared_ptr<cvc::app> VrHostQtBridgeTest::app;
std::unique_ptr<volrover3::EmbeddedInterpreter> VrHostQtBridgeTest::interp;

TEST_F(VrHostQtBridgeTest, ReadWidgetsAndAddMenuFiringMessageBox) {
  // A live QMainWindow standing in for the app's window, with a menu bar shaped
  // like volrover3's. Hand its raw address to the host, exactly as MainWindow
  // does, so vrhost.main_window() adopts THIS window via shiboken6.wrapInstance.
  QMainWindow win;
  win.setWindowTitle("volrover3");
  win.menuBar()->addMenu("&File");
  win.menuBar()->addMenu("&View");
  win.menuBar()->addMenu("&Tools");
  interp->set_main_window_ptr(reinterpret_cast<std::uintptr_t>(&win));

  // Single-Qt litmus + PySide6 present, and wrap the live window ONCE into `mw`
  // (repeated vrhost.main_window() would mint independent shiboken wrappers whose
  // GC invalidates each other's child objects). Everything below reuses `mw`.
  const char *scripts = std::getenv("VOLROVER3_SCRIPTS_DIR");
  ASSERT_NE(scripts, nullptr) << "VOLROVER3_SCRIPTS_DIR not set";
  ASSERT_TRUE(interp->run_string(
      std::string("import sys; sys.path.insert(0, r'") + scripts + "')\n" +
      "import vrhost\n"
      "from PySide6 import QtWidgets\n"
      "assert QtWidgets.QApplication.instance() is not None, 'dual-Qt: no shared qApp'\n"
      "mw = vrhost.main_window()\n"))
      << "PySide6 import / single-Qt litmus failed (see stderr)";

  // READ: Python reads the live window's title + menu titles (single expression
  // -> run_string_capture echoes the repr).
  std::string out, err;
  ASSERT_TRUE(interp->run_string_capture(
      "(mw.windowTitle(), [m.title().replace('&','') for m in "
      "  mw.menuBar().findChildren(QtWidgets.QMenu)])",
      out, err))
      << err;
  EXPECT_NE(out.find("volrover3"), std::string::npos) << out;
  EXPECT_NE(out.find("File"), std::string::npos) << out;
  EXPECT_NE(out.find("Tools"), std::string::npos) << out;

  // MUTATE: run the shipped demo script against the same wrapped window.
  ASSERT_TRUE(interp->run_string(
      "import qt_bridge_demo\n"
      "qt_bridge_demo.read_widgets(mw)\n"
      "qt_bridge_demo.install(mw)\n"))
      << "qt_bridge_demo failed (see stderr)";

  // Assert from C++ that the mutation reached the REAL QMainWindow.
  QMenu *pymenu = win.menuBar()->findChild<QMenu *>("PythonMenu");
  ASSERT_NE(pymenu, nullptr) << "Python did not add the menu to the live C++ window";
  EXPECT_EQ(pymenu->title(), QString("&Python"));
  QAction *act = win.findChild<QAction *>("PythonHelloAction");
  ASSERT_NE(act, nullptr) << "the Python action is not on the live window";

  // TRIGGER from C++: fires the Python slot -> QMessageBox (self-dismissing).
  act->trigger();
  QCoreApplication::processEvents();

  // The handler set vrhost._bridge_hello_fired — proves trigger -> Python -> box.
  ASSERT_TRUE(interp->run_string_capture("getattr(vrhost, '_bridge_hello_fired', False)", out, err))
      << err;
  EXPECT_NE(out.find("True"), std::string::npos)
      << "the menu action's handler (QMessageBox) did not fire; out=[" << out << "]";
}

int main(int argc, char **argv) {
  QApplication qapp(argc, argv); // widgets need a QApplication (offscreen via env)
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
