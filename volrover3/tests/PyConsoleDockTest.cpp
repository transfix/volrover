// Phase 5: the PyConsoleDock REPL evaluates + displays, and its Jobs tab
// reflects the scheduler. Runs under QT_QPA_PLATFORM=offscreen with one
// interpreter for the suite (CPython won't cleanly re-init after finalize).

#include <cvc/core/app.h>
#include <gtest/gtest.h>
#include <volrover3/EmbeddedInterpreter.h>
#include <volrover3/JobScheduler.h>
#include <volrover3/PyConsoleDock.h>

#include <QApplication>

#include <memory>

class PyConsoleDockTest : public ::testing::Test {
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
std::shared_ptr<cvc::app> PyConsoleDockTest::app;
std::unique_ptr<volrover3::EmbeddedInterpreter> PyConsoleDockTest::interp;

TEST_F(PyConsoleDockTest, ReplEvaluatesAndDisplays) {
  volrover3::JobScheduler sched(interp.get(), 50);
  volrover3::PyConsoleDock dock(interp.get(), &sched);

  dock.evaluate("print('console works')");
  EXPECT_NE(dock.outputText().indexOf(">>> print('console works')"), -1); // echo
  EXPECT_NE(dock.outputText().indexOf("console works"), -1);              // captured stdout

  dock.evaluate("2 + 3");
  EXPECT_NE(dock.outputText().indexOf("5"), -1); // bare expression echoes its repr
}

TEST_F(PyConsoleDockTest, JobsTabReflectsScheduler) {
  volrover3::JobScheduler sched(interp.get(), 50);
  volrover3::PyConsoleDock dock(interp.get(), &sched);

  EXPECT_EQ(dock.jobRowCount(), 0);
  const int id = sched.submit("job1", "def step(dt):\n pass\n");
  ASSERT_GE(id, 0);
  dock.refreshJobs();
  EXPECT_EQ(dock.jobRowCount(), 1);
}

int main(int argc, char **argv) {
  QApplication qapp(argc, argv); // widgets need QApplication (offscreen via env)
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
