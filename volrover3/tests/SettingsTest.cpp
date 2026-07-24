// Phase 1: volrover3::Settings is backed by the cvc::state tree (per-instance
// section under the global root) and persists to ~/.volrover/settings.yaml +
// data.db. HOME is redirected to a temp dir by the test's ctest ENVIRONMENT.

#include <cvc/core/app.h>
#include <cvc/core/state.h>
#include <gtest/gtest.h>
#include <volrover3/Settings.h>

#include <QCoreApplication>

#include <cstdlib>
#include <filesystem>
#include <memory>

TEST(SettingsTest, StateBackedAndPersisted) {
  // First session: defaults, state-backing, then a change we persist.
  {
    auto app = std::make_shared<cvc::app>();
    volrover3::Settings s(app, "volrover3");
    s.load();

    // Defaults.
    EXPECT_EQ(s.mode(), volrover3::InterpreterMode::Single);
    EXPECT_EQ(s.tickMs(), 100);
    EXPECT_EQ(s.historySize(), 500);

    // State-backing: the value lives in cvc::state under the instance section,
    // so a Python script hitting volrover3.settings.* sees the same value.
    EXPECT_EQ(cvc::state::instance(*app)("volrover3")("settings")("interpreter.mode").value(),
              std::string("single"));

    // Change (persisted to yaml) + a KV blob (persisted to data.db).
    s.setMode(volrover3::InterpreterMode::Multi);
    s.kvPut("greeting", "hello");
    EXPECT_EQ(cvc::state::instance(*app)("volrover3")("settings")("interpreter.mode").value(),
              std::string("multi"));
  }

  // Second session (fresh app): loads from disk -> the change persisted.
  {
    auto app2 = std::make_shared<cvc::app>();
    volrover3::Settings s2(app2, "volrover3");
    s2.load();
    EXPECT_EQ(s2.mode(), volrover3::InterpreterMode::Multi); // from settings.yaml
    EXPECT_EQ(s2.kvGet("greeting"), std::string("hello"));   // from data.db
  }
}

TEST(SettingsTest, InstancesGetDistinctSections) {
  auto app = std::make_shared<cvc::app>();
  volrover3::Settings a(app, "volrover3");
  volrover3::Settings b(app, "volrover3.1");
  a.load();
  b.load();
  // Distinct sections of the one global state root — no collision.
  EXPECT_EQ(a.instancePrefix(), std::string("volrover3"));
  EXPECT_EQ(b.instancePrefix(), std::string("volrover3.1"));
  EXPECT_NE(&cvc::state::instance(*app)("volrover3")("settings"),
            &cvc::state::instance(*app)("volrover3.1")("settings"));
}

// QSqlDatabase (used by Settings::load for data.db) requires a QCoreApplication.
int main(int argc, char **argv) {
  QCoreApplication qapp(argc, argv);
  // Start from a clean ~/.volrover so a persisted value from a prior run (HOME is
  // the ctest-redirected temp dir) can't leak in and fail the defaults check.
  if (const char *home = std::getenv("HOME")) {
    std::error_code ec;
    std::filesystem::remove_all(std::filesystem::path(home) / ".volrover", ec);
  }
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
