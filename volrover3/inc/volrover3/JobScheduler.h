#ifndef VOLROVER3_JOBSCHEDULER_H
#define VOLROVER3_JOBSCHEDULER_H

// --------------------------------------------------------------------
// volrover3::JobScheduler — the embedded-interpreter job scheduler.
// --------------------------------------------------------------------
// UI-thread cooperative scheduler (docs/EMBEDDED_PYTHON.md §12.1), modeled on
// verlihub's OnTimer/CallAll tick: a QTimer drives tick(), which fans out over
// an ordered registry calling each job's cached step(dt) in order — no
// preemption, per-job exception isolation, one GIL acquisition per tick.
// Single-interpreter mode: jobs share the process interpreter, each isolated in
// its own module namespace (vr_job_<id>). A job "subscribes" by defining a
// top-level `step(dt)` (verlihub's OnTimer convention).
//
// Cooperative interrupt/stop live here; the hard-kill sacrificial-worker-thread
// path (for force-terminating a runaway C-extension loop) is layered on top.
// --------------------------------------------------------------------

#include <volrover3/EmbeddedInterpreter.h> // InterpreterMode

#include <QObject>

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

class QTimer;

namespace volrover3 {

// Mirrors libcvc state_exec's process_status vocabulary so the jobs tab renders
// uniformly across Python + DSL jobs.
enum class JobStatus { Ready, Running, Paused, Terminated, Killed, Error };

struct JobInfo {
  int id = -1;
  std::string name;
  JobStatus status = JobStatus::Ready;
  bool onWorker = false;
  std::uint64_t steps = 0;
  double elapsed = 0.0; // seconds since submit
  std::string lastError;
};

class JobScheduler : public QObject {
  Q_OBJECT
public:
  // `mode` selects the job launcher: Single = jobs share the process
  // interpreter (each in its own module namespace); Multi = each job gets its
  // own Py_NewInterpreter sub-interpreter with an enforced import gate
  // (pycvc/pycvc_gl/vrhost/vtk/numpy denied — they corrupt across
  // sub-interpreters; docs/EMBEDDED_PYTHON.md §12.2).
  JobScheduler(EmbeddedInterpreter *interp, int tickMs,
               InterpreterMode mode = InterpreterMode::Single, QObject *parent = nullptr);
  ~JobScheduler() override;

  InterpreterMode mode() const { return m_mode; }

  // Submit a Python job. `source` is exec'd ONCE into a fresh module namespace
  // to define a `step(dt)` callable (verlihub's OnTimer). The scheduler then
  // ticks step(dt) cooperatively. If `onWorker`, the job runs on a sacrificial
  // worker thread so it can be hard-killed (kill() force-terminates it).
  // Returns the job id (>=0), or -1 if the source failed to load / has no step.
  int submit(const std::string &name, const std::string &source, bool onWorker = false);

  std::vector<JobInfo> listJobs() const;
  bool pause(int id);
  bool resume(int id);
  // Cooperative interrupt: raise KeyboardInterrupt in the job's thread at the
  // next bytecode boundary (PyThreadState_SetAsyncExc). Best-effort — cannot
  // break a tight C-extension loop mid-call.
  bool interrupt(int id);
  // Stop + remove the job. Cooperative for tick jobs (drop from the registry +
  // clear the namespace); for worker jobs, requests cancel then, if it will not
  // yield, DETACHES the sacrificial thread and marks the job Killed (unclean).
  bool kill(int id);

  void start(); // start the QTimer
  void stop();  // stop the QTimer
  int tickMs() const { return m_tickMs; }
  std::size_t size() const { return m_jobs.size(); }

public slots:
  // One cooperative pass over the registry. The QTimer target; also callable
  // directly (tests, or a manual pump).
  void tick();

private:
  struct Job; // defined in the .cpp (holds PyObject* + optional worker thread)
  Job *find(int id);
  const Job *find(int id) const;
  // Run one step(dt) under the GIL (released before return). Static so they can
  // touch the private Job; workerLoop is the sacrificial-thread body.
  static bool runStepUnderGil(Job *j, double dt, std::string &err, unsigned long &tid);
  static void workerLoop(Job *j);
  // Free a job's Python resources: single -> decref; multi -> Py_EndInterpreter
  // its sub-interpreter. Its worker (if any) must already be stopped/joined.
  static void teardownJob(Job *j);

  EmbeddedInterpreter *m_interp;
  int m_tickMs;
  InterpreterMode m_mode;
  int m_nextId = 0;
  QTimer *m_timer = nullptr;
  std::int64_t m_lastTickNs = 0;
  std::vector<std::unique_ptr<Job>> m_jobs;
  // Hung worker jobs that could not be joined: detached + kept alive here so the
  // abandoned thread never dereferences freed Job memory (see EMBEDDED_PYTHON.md
  // §12.6 hard-kill finding). Never freed, never joined.
  std::vector<std::unique_ptr<Job>> m_zombies;
};

} // namespace volrover3

#endif // VOLROVER3_JOBSCHEDULER_H
