// Python.h must be included before any system/Qt headers.
// clang-format off
#include <Python.h>
// clang-format on

#include <volrover3/EmbeddedInterpreter.h>
#include <volrover3/JobScheduler.h>

#include <QTimer>

#include <atomic>
#include <chrono>
#include <mutex>
#include <string>
#include <thread>

namespace volrover3 {

using namespace std::chrono_literals;

struct JobScheduler::Job {
  int id = -1;
  std::string name;
  bool onWorker = false;
  PyObject *module = nullptr; // owned; the job's isolated namespace (vr_job_<id>)
  PyObject *stepFn = nullptr; // owned; the cached step(dt) callable

  mutable std::mutex mtx;     // guards the fields below (cross-thread for workers)
  JobStatus status = JobStatus::Ready;
  std::uint64_t steps = 0;
  std::string lastError;
  unsigned long threadId = 0; // Python thread id running this job (for SetAsyncExc)
  std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

  // worker path
  std::thread thread;
  std::atomic<bool> stopFlag{false};
  std::atomic<bool> finished{false};
  int tickMs = 100;
};

namespace {
// GIL must be held. Consumes + formats the current Python error, or "".
std::string fetchError() {
  if (!PyErr_Occurred())
    return {};
  PyObject *t = nullptr, *v = nullptr, *tb = nullptr;
  PyErr_Fetch(&t, &v, &tb);
  PyErr_NormalizeException(&t, &v, &tb);
  std::string msg;
  if (v) {
    if (PyObject *s = PyObject_Str(v)) {
      if (const char *c = PyUnicode_AsUTF8(s))
        msg = c;
      Py_DECREF(s);
    }
  }
  Py_XDECREF(t);
  Py_XDECREF(v);
  Py_XDECREF(tb);
  return msg;
}

std::int64_t nowNs() {
  return std::chrono::duration_cast<std::chrono::nanoseconds>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}
} // namespace

// Runs the job's step(dt) once (GIL acquired here). Returns true on success;
// on error fills `err`. Records the running thread's Python id. GIL is released
// before returning, so state updates by the caller take no GIL.
bool JobScheduler::runStepUnderGil(Job *j, double dt, std::string &err, unsigned long &tid) {
  PyGILState_STATE gil = PyGILState_Ensure();
  tid = PyThread_get_thread_ident();
  bool ok = false;
  if (PyObject *res = PyObject_CallFunction(j->stepFn, "d", dt)) {
    Py_DECREF(res);
    ok = true;
  } else {
    err = fetchError();
  }
  PyGILState_Release(gil);
  return ok;
}

// The sacrificial worker loop: step(dt) on its OWN thread, releasing the GIL
// between steps so a yielding slow job never freezes the app and can be stopped
// cleanly. dt is real elapsed time.
void JobScheduler::workerLoop(Job *j) {
  auto last = std::chrono::steady_clock::now();
  while (!j->stopFlag.load()) {
    bool paused;
    {
      std::lock_guard<std::mutex> lk(j->mtx);
      paused = (j->status == JobStatus::Paused);
    }
    if (!paused) {
      auto now = std::chrono::steady_clock::now();
      double dt = std::chrono::duration<double>(now - last).count();
      last = now;
      std::string err;
      unsigned long tid = 0;
      bool ok = runStepUnderGil(j, dt, err, tid); // GIL released inside
      std::lock_guard<std::mutex> lk(j->mtx);
      j->threadId = tid;
      if (ok) {
        ++j->steps;
        j->status = JobStatus::Running;
      } else {
        j->lastError = err;
        j->status = JobStatus::Error;
        break; // a raising worker stops itself
      }
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(j->tickMs));
  }
  j->finished.store(true);
}

JobScheduler::JobScheduler(EmbeddedInterpreter *interp, int tickMs, QObject *parent)
    : QObject(parent), m_interp(interp), m_tickMs(tickMs) {
  m_timer = new QTimer(this);
  m_timer->setInterval(tickMs);
  connect(m_timer, &QTimer::timeout, this, &JobScheduler::tick);
  m_lastTickNs = nowNs();
}

JobScheduler::~JobScheduler() {
  // Stop + join every worker (best-effort) before freeing anything.
  for (auto &j : m_jobs) {
    if (j->thread.joinable()) {
      j->stopFlag.store(true);
      j->thread.join();
    }
  }
  if (m_interp && m_interp->ok()) {
    PyGILState_STATE gil = PyGILState_Ensure();
    for (auto &j : m_jobs) {
      Py_XDECREF(j->stepFn);
      Py_XDECREF(j->module);
    }
    PyGILState_Release(gil);
  }
  m_jobs.clear();
  // m_zombies are intentionally leaked (their detached threads may still run).
  for (auto &z : m_zombies)
    (void)z.release();
}

int JobScheduler::submit(const std::string &name, const std::string &source, bool onWorker) {
  if (!m_interp || !m_interp->ok())
    return -1;
  const int id = m_nextId++;

  PyGILState_STATE gil = PyGILState_Ensure();
  const std::string modName = "vr_job_" + std::to_string(id);
  PyObject *mod = PyModule_New(modName.c_str()); // new
  PyObject *step = nullptr;
  bool loaded = false;
  if (mod) {
    PyObject *dict = PyModule_GetDict(mod); // borrowed
    if (PyDict_GetItemString(dict, "__builtins__") == nullptr)
      PyDict_SetItemString(dict, "__builtins__", PyEval_GetBuiltins());
    if (PyObject *res = PyRun_String(source.c_str(), Py_file_input, dict, dict)) {
      Py_DECREF(res);
      step = PyObject_GetAttrString(mod, "step"); // new or null
      loaded = step && PyCallable_Check(step);
    }
  }
  if (!loaded) {
    (void)fetchError();
    Py_XDECREF(step);
    Py_XDECREF(mod);
    PyGILState_Release(gil);
    return -1;
  }
  PyGILState_Release(gil);

  auto job = std::make_unique<Job>();
  job->id = id;
  job->name = name;
  job->onWorker = onWorker;
  job->module = mod;
  job->stepFn = step;
  job->status = JobStatus::Ready;
  job->tickMs = m_tickMs;
  Job *raw = job.get();
  m_jobs.push_back(std::move(job));
  if (onWorker)
    raw->thread = std::thread(&JobScheduler::workerLoop, raw); // self-driving
  return id;
}

void JobScheduler::tick() {
  if (!m_interp || !m_interp->ok() || m_jobs.empty())
    return;
  const std::int64_t t = nowNs();
  double dt = (t - m_lastTickNs) / 1e9;
  if (dt < 0)
    dt = 0;
  m_lastTickNs = t;

  for (auto &j : m_jobs) {
    if (j->onWorker || !j->stepFn) // worker jobs self-drive
      continue;
    {
      std::lock_guard<std::mutex> lk(j->mtx);
      if (j->status != JobStatus::Ready && j->status != JobStatus::Running)
        continue;
    }
    std::string err;
    unsigned long tid = 0;
    const bool ok = runStepUnderGil(j.get(), dt, err, tid);
    std::lock_guard<std::mutex> lk(j->mtx);
    j->threadId = tid;
    if (ok) {
      ++j->steps;
      j->status = JobStatus::Running;
    } else {
      j->lastError = err; // isolate: one bad job can't kill the loop
      j->status = JobStatus::Error;
    }
  }
}

std::vector<JobInfo> JobScheduler::listJobs() const {
  std::vector<JobInfo> v;
  v.reserve(m_jobs.size() + m_zombies.size());
  const auto now = std::chrono::steady_clock::now();
  auto snap = [&](const std::unique_ptr<Job> &j) {
    std::lock_guard<std::mutex> lk(j->mtx);
    JobInfo i;
    i.id = j->id;
    i.name = j->name;
    i.status = j->status;
    i.onWorker = j->onWorker;
    i.steps = j->steps;
    i.elapsed = std::chrono::duration<double>(now - j->start).count();
    i.lastError = j->lastError;
    v.push_back(std::move(i));
  };
  for (const auto &j : m_jobs)
    snap(j);
  for (const auto &z : m_zombies)
    snap(z);
  return v;
}

JobScheduler::Job *JobScheduler::find(int id) {
  for (auto &j : m_jobs)
    if (j->id == id)
      return j.get();
  return nullptr;
}
const JobScheduler::Job *JobScheduler::find(int id) const {
  for (const auto &j : m_jobs)
    if (j->id == id)
      return j.get();
  return nullptr;
}

bool JobScheduler::pause(int id) {
  Job *j = find(id);
  if (!j)
    return false;
  std::lock_guard<std::mutex> lk(j->mtx);
  if (j->status == JobStatus::Running || j->status == JobStatus::Ready)
    j->status = JobStatus::Paused;
  return true;
}

bool JobScheduler::resume(int id) {
  Job *j = find(id);
  if (!j)
    return false;
  std::lock_guard<std::mutex> lk(j->mtx);
  if (j->status == JobStatus::Paused)
    j->status = JobStatus::Ready;
  return true;
}

bool JobScheduler::interrupt(int id) {
  Job *j = find(id);
  if (!j)
    return false;
  unsigned long tid;
  {
    std::lock_guard<std::mutex> lk(j->mtx);
    tid = j->threadId;
  }
  if (tid == 0)
    return false;
  PyGILState_STATE gil = PyGILState_Ensure();
  const int n = PyThreadState_SetAsyncExc(tid, PyExc_KeyboardInterrupt);
  PyGILState_Release(gil);
  return n == 1;
}

bool JobScheduler::kill(int id) {
  auto it = m_jobs.begin();
  for (; it != m_jobs.end(); ++it)
    if ((*it)->id == id)
      break;
  if (it == m_jobs.end())
    return false;
  Job *j = it->get();

  if (!j->onWorker) {
    // Cooperative: just drop it + clear the namespace.
    PyGILState_STATE gil = PyGILState_Ensure();
    Py_XDECREF(j->stepFn);
    Py_XDECREF(j->module);
    PyGILState_Release(gil);
    m_jobs.erase(it);
    return true;
  }

  // Worker: request stop, nudge with KeyboardInterrupt, wait for a clean exit.
  j->stopFlag.store(true);
  interrupt(id); // best-effort SetAsyncExc into a Python loop
  const auto deadline = std::chrono::steady_clock::now() + 1500ms;
  while (!j->finished.load() && std::chrono::steady_clock::now() < deadline)
    std::this_thread::sleep_for(10ms);

  if (j->finished.load()) {
    j->thread.join();
    PyGILState_STATE gil = PyGILState_Ensure();
    Py_XDECREF(j->stepFn);
    Py_XDECREF(j->module);
    PyGILState_Release(gil);
    m_jobs.erase(it);
  } else {
    // Genuinely hung (a tight C loop holding the GIL). Cannot force-kill safely:
    // detach + quarantine the Job so the abandoned thread never dereferences
    // freed memory. Leaked by design; the PyObjects are NOT freed.
    {
      std::lock_guard<std::mutex> lk(j->mtx);
      j->status = JobStatus::Killed;
    }
    j->thread.detach();
    m_zombies.push_back(std::move(*it));
    m_jobs.erase(it);
  }
  return true;
}

void JobScheduler::start() {
  if (m_timer)
    m_timer->start();
}
void JobScheduler::stop() {
  if (m_timer)
    m_timer->stop();
}

} // namespace volrover3
