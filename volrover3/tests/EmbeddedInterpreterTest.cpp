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

  // A well-formed script runs and returns true.
  EXPECT_TRUE(interp.run_string("x = 1 + 1\nassert x == 2\n"));

  // A raising script is caught + reported, returns false — never throws/crashes.
  EXPECT_FALSE(interp.run_string("raise RuntimeError('boom')\n"));

  // The host facade delivers the SAME app object we injected (identity, not a
  // copy) — the whole point of the injected-app design.
  ASSERT_TRUE(interp.host() != nullptr);
  EXPECT_EQ(interp.host()->app().get(), app.get());
}
