// Phase 0c/1 smoke test: the embedded CPython interpreter boots inside
// volrover3_lib, runs a script under the GIL, isolates a failing script, and
// hands back a PyHost that delivers the SAME live cvc::app it was built with.
// (The full vrhost.host.app() -> pycvc round-trip is Phase 1, gated on the
// post-#136 pycvc being importable.)

#include <cvc/core/app.h>
#include <gtest/gtest.h>
#include <volrover3/EmbeddedInterpreter.h>
#include <volrover3/PyHost.h>

#include <memory>

TEST(EmbeddedInterpreterTest, BootsAndRunsAScript) {
  auto app = std::make_shared<cvc::app>();
  volrover3::EmbeddedInterpreter interp(app, /*scene=*/nullptr);
  ASSERT_TRUE(interp.ok()) << "CPython failed to initialize (VOLROVER3_PYTHON_HOME?)";

  // Config default mode is Single (the full-capability mode).
  EXPECT_EQ(interp.mode(), volrover3::InterpreterMode::Single);

  // A well-formed script runs and returns true.
  EXPECT_TRUE(interp.run_string("x = 1 + 1\nassert x == 2\n"));

  // A raising script is caught + reported, returns false — never throws/crashes.
  EXPECT_FALSE(interp.run_string("raise RuntimeError('boom')\n"));

  // The host facade delivers the SAME app object we injected (identity, not a
  // copy) — the whole point of the injected-app design.
  ASSERT_TRUE(interp.host() != nullptr);
  EXPECT_EQ(interp.host()->app().get(), app.get());
}

TEST(EmbeddedInterpreterTest, ReplCaptureEchoesAndIsolatesErrors) {
  auto app = std::make_shared<cvc::app>();
  volrover3::EmbeddedInterpreter interp(app, /*scene=*/nullptr);
  ASSERT_TRUE(interp.ok());

  std::string out, err;
  // print() output is captured.
  EXPECT_TRUE(interp.run_string_capture("print('hi there')", out, err));
  EXPECT_NE(out.find("hi there"), std::string::npos);
  EXPECT_TRUE(err.empty());

  // A bare expression echoes its repr (REPL behavior, Py_single_input).
  out.clear();
  err.clear();
  EXPECT_TRUE(interp.run_string_capture("6 * 7", out, err));
  EXPECT_NE(out.find("42"), std::string::npos);

  // A raising snippet: false, traceback captured in err, streams restored.
  out.clear();
  err.clear();
  EXPECT_FALSE(interp.run_string_capture("raise ValueError('boom')", out, err));
  EXPECT_NE(err.find("ValueError"), std::string::npos);
  EXPECT_NE(err.find("boom"), std::string::npos);

  // Streams were restored: a subsequent normal capture still works.
  out.clear();
  err.clear();
  EXPECT_TRUE(interp.run_string_capture("print('again')", out, err));
  EXPECT_NE(out.find("again"), std::string::npos);
}
