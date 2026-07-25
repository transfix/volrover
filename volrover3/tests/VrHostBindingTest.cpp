// Phase 1 capstone: the vrhost host-control module delivers the LIVE cvc::app to
// embedded pycvc scripts. This proves the whole chain end-to-end:
//
//   C++ owns `app` --> EmbeddedInterpreter boots CPython, imports the pure-Python
//   vrhost shim, and injects the app as a PyCapsule("cvc.app") --> a script calls
//   `pycvc.state_set(vrhost.app(), path, value)` (vrhost.app() feeds the capsule
//   to pycvc.app_from_capsule, which wraps it into pycvc's OWN type) --> we read
//   that value back out of the SAME C++ cvc::state tree.
//
// No SWIG in volrover3 and no pycvc.i: pycvc wraps the app with its own type, so
// the handle is compatible with state_set by construction (no cross-module type
// sharing, no SWIG-runtime-version coupling).
//
// One interpreter is booted for the whole suite (SetUpTestSuite) — CPython does
// not cleanly re-Py_Initialize after Py_Finalize within a process. The module
// dir + python home come from the environment set by CTest (see CMakeLists.txt).

#include <cvc/core/app.h>
#include <cvc/core/state.h>
#include <gtest/gtest.h>
#include <volrover3/EmbeddedInterpreter.h>

#include <memory>
#include <string>

class VrHostBindingTest : public ::testing::Test {
protected:
  static std::shared_ptr<cvc::app> app;
  static std::unique_ptr<volrover3::EmbeddedInterpreter> interp;
  static void SetUpTestSuite() {
    app = std::make_shared<cvc::app>();
    interp = std::make_unique<volrover3::EmbeddedInterpreter>(app, nullptr);
    ASSERT_TRUE(interp->ok());
  }
  static void TearDownTestSuite() {
    interp.reset();
    app.reset();
  }
};
std::shared_ptr<cvc::app> VrHostBindingTest::app;
std::unique_ptr<volrover3::EmbeddedInterpreter> VrHostBindingTest::interp;

// The host bound at boot, and `vrhost.host` / `vrhost.app()` are the live app.
// Exercises the console's REPL capture path — one statement each (run_string_capture
// uses Py_single_input, so it echoes a bare expression's repr like a real REPL).
TEST_F(VrHostBindingTest, HostBoundAtBoot) {
  ASSERT_TRUE(interp->host_bound()) << "vrhost did not import/bind — is "
                                       "VOLROVER3_PYMODULE_PATH set to the pymod dir?";
  ASSERT_TRUE(interp->run_string("import vrhost, pycvc"));
  std::string out, err;
  ASSERT_TRUE(interp->run_string_capture(
      "vrhost.app() is not None", out, err))
      << err;
  EXPECT_NE(out.find("True"), std::string::npos) << "out=[" << out << "] err=[" << err << "]";
}

// The capstone: a script writes state THROUGH vrhost.app(); C++ reads it
// back from the same app's cvc::state tree (run_string = Py_file_input, multi-line).
TEST_F(VrHostBindingTest, ScriptWritesLiveAppState) {
  ASSERT_TRUE(interp->run_string(
      "import vrhost, pycvc\n"
      "pycvc.state_set(vrhost.app(), 'volrover3.vrhost_demo', '99')\n"))
      << "the script raised — see the traceback on stderr";

  cvc::state *node = cvc::state::instance(*app).findDescendant("volrover3.vrhost_demo");
  ASSERT_NE(node, nullptr) << "script's state write did not reach the C++ app tree";
  EXPECT_EQ(node->value(), "99");
}

// The host handle really is the host's OWN app, not a process-global: a value set
// on vrhost.app() must NOT be visible on a fresh, independent make_app(). The
// boolean result is stashed back on the host app so C++ can read it (make_app()'s
// tree is not reachable from C++ here).
TEST_F(VrHostBindingTest, HostAppIsDistinctFromMakeApp) {
  ASSERT_TRUE(interp->run_string(
      "import vrhost, pycvc\n"
      "pycvc.state_set(vrhost.app(), 'volrover3.iso_probe', 'host')\n"
      "other = pycvc.make_app()\n"
      "pycvc.state_set(vrhost.app(), 'volrover3.iso_result',\n"
      "                str(pycvc.state_has(other, 'volrover3.iso_probe')))\n"));

  cvc::state *node = cvc::state::instance(*app).findDescendant("volrover3.iso_result");
  ASSERT_NE(node, nullptr);
  EXPECT_EQ(node->value(), "False") << "the host app was not distinct from make_app()";
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
