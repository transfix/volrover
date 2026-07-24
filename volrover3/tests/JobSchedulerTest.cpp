// Phase 4: JobScheduler cooperative tick over a per-job namespace registry.
// One interpreter is booted for the whole suite (SetUpTestSuite) — CPython does
// not cleanly re-Py_Initialize after Py_Finalize within a process.

#include <cvc/core/app.h>
#include <gtest/gtest.h>
#include <volrover3/EmbeddedInterpreter.h>
#include <volrover3/JobScheduler.h>

#include <QCoreApplication>

#include <chrono>
#include <memory>
#include <thread>

using volrover3::JobStatus;

class JobSchedulerTest : public ::testing::Test {
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
std::shared_ptr<cvc::app> JobSchedulerTest::app;
std::unique_ptr<volrover3::EmbeddedInterpreter> JobSchedulerTest::interp;

TEST_F(JobSchedulerTest, CooperativeTickAndControl) {
  volrover3::JobScheduler sched(interp.get(), 50);

  const int id = sched.submit("counter", "n = 0\ndef step(dt):\n    global n\n    n += 1\n");
  ASSERT_GE(id, 0);
  EXPECT_EQ(sched.size(), 1u);

  sched.tick();
  sched.tick();
  sched.tick();
  auto jobs = sched.listJobs();
  ASSERT_EQ(jobs.size(), 1u);
  EXPECT_EQ(jobs[0].steps, 3u);
  EXPECT_EQ(jobs[0].status, JobStatus::Running);

  // pause -> not stepped; resume -> stepped again.
  EXPECT_TRUE(sched.pause(id));
  sched.tick();
  EXPECT_EQ(sched.listJobs()[0].steps, 3u);
  EXPECT_EQ(sched.listJobs()[0].status, JobStatus::Paused);
  EXPECT_TRUE(sched.resume(id));
  sched.tick();
  EXPECT_EQ(sched.listJobs()[0].steps, 4u);

  // kill -> removed from the registry.
  EXPECT_TRUE(sched.kill(id));
  EXPECT_EQ(sched.size(), 0u);
}

TEST_F(JobSchedulerTest, RejectsNoStepAndIsolatesRaises) {
  volrover3::JobScheduler sched(interp.get(), 50);

  EXPECT_EQ(sched.submit("nostep", "x = 1\n"), -1);       // no step()
  EXPECT_EQ(sched.submit("bad", "def step(dt)\n"), -1);   // syntax error
  EXPECT_EQ(sched.size(), 0u);

  const int good = sched.submit("good", "n=0\ndef step(dt):\n global n\n n+=1\n");
  const int bad = sched.submit("raiser", "def step(dt):\n raise ValueError('boom')\n");
  ASSERT_GE(good, 0);
  ASSERT_GE(bad, 0);

  sched.tick(); // one bad job must not stop the loop
  for (const auto &j : sched.listJobs()) {
    if (j.name == "good") {
      EXPECT_EQ(j.steps, 1u);
      EXPECT_EQ(j.status, JobStatus::Running);
    }
    if (j.name == "raiser") {
      EXPECT_EQ(j.status, JobStatus::Error);
      EXPECT_NE(j.lastError.find("boom"), std::string::npos);
    }
  }
}

TEST_F(JobSchedulerTest, WorkerJobStepsOffThreadAndStopsCleanly) {
  volrover3::JobScheduler sched(interp.get(), 10); // 10 ms worker cadence

  const int id =
      sched.submit("worker", "n=0\ndef step(dt):\n global n\n n+=1\n", /*onWorker=*/true);
  ASSERT_GE(id, 0);

  // The worker self-drives on its OWN thread — no tick() calls here. Give it
  // time to step (it releases the GIL between steps).
  std::this_thread::sleep_for(std::chrono::milliseconds(200));

  auto jobs = sched.listJobs();
  ASSERT_EQ(jobs.size(), 1u);
  EXPECT_TRUE(jobs[0].onWorker);
  EXPECT_GT(jobs[0].steps, 0u); // ran off the UI thread without any scheduler tick

  // Clean kill: stop flag -> the worker exits between steps -> joined -> removed.
  EXPECT_TRUE(sched.kill(id));
  EXPECT_EQ(sched.size(), 0u);
}

int main(int argc, char **argv) {
  QCoreApplication qapp(argc, argv); // JobScheduler is a QObject / owns a QTimer
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
