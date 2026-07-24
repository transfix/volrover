// Python.h must be included before any system/Qt headers.
// clang-format off
#include <Python.h>
// clang-format on

#include <volrover3/EmbeddedInterpreter.h>
#include <volrover3/JobScheduler.h>

#include <QTimer>

#include <chrono>
#include <string>

namespace volrover3 {

struct JobScheduler::Job {
  int id = -1;
  std::string name;
  JobStatus status = JobStatus::Ready;
  bool onWorker = false;
  PyObject *module = nullptr;  // owned; the job's isolated namespace (vr_job_<id>)
  PyObject *stepFn = nullptr;  // owned; the cached step(dt) callable
  std::uint64_t steps = 0;
  std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
  std::string lastError;
  unsigned long threadId = 0; // Python thread id running this job (for SetAsyncExc)
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

JobScheduler::JobScheduler(EmbeddedInterpreter *interp, int tickMs, QObject *parent)
    : QObject(parent), m_interp(interp), m_tickMs(tickMs) {
  m_timer = new QTimer(this);
  m_timer->setInterval(tickMs);
  connect(m_timer, &QTimer::timeout, this, &JobScheduler::tick);
  m_lastTickNs = nowNs();
}

JobScheduler::~JobScheduler() {
  if (m_interp && m_interp->ok() && !m_jobs.empty()) {
    PyGILState_STATE gil = PyGILState_Ensure();
    for (auto &j : m_jobs) {
      Py_XDECREF(j->stepFn);
      Py_XDECREF(j->module);
    }
    PyGILState_Release(gil);
  }
  m_jobs.clear();
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
  if (!loaded)
    (void)fetchError(); // clear any load error
  if (!loaded) {
    Py_XDECREF(step);
    Py_XDECREF(mod);
    PyGILState_Release(gil);
    return -1;
  }
  PyGILState_Release(gil);

  auto job = std::make_unique<Job>();
  job->id = id;
  job->name = name;
  job->onWorker = onWorker; // wired by the worker path (Phase 4b); ticked cooperatively for now
  job->module = mod;
  job->stepFn = step;
  job->status = JobStatus::Ready;
  m_jobs.push_back(std::move(job));
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

  PyGILState_STATE gil = PyGILState_Ensure();
  const unsigned long tid = PyThread_get_thread_ident();
  for (auto &j : m_jobs) {
    if (!j->stepFn)
      continue;
    if (j->status != JobStatus::Ready && j->status != JobStatus::Running)
      continue; // skip paused/terminated/killed/error
    j->threadId = tid;
    if (PyObject *res = PyObject_CallFunction(j->stepFn, "d", dt)) {
      Py_DECREF(res);
      ++j->steps;
      j->status = JobStatus::Running;
    } else {
      j->lastError = fetchError(); // isolate: one bad job can't kill the loop
      j->status = JobStatus::Error;
    }
  }
  PyGILState_Release(gil);
}

std::vector<JobInfo> JobScheduler::listJobs() const {
  std::vector<JobInfo> v;
  v.reserve(m_jobs.size());
  const auto now = std::chrono::steady_clock::now();
  for (const auto &j : m_jobs) {
    JobInfo i;
    i.id = j->id;
    i.name = j->name;
    i.status = j->status;
    i.onWorker = j->onWorker;
    i.steps = j->steps;
    i.elapsed = std::chrono::duration<double>(now - j->start).count();
    i.lastError = j->lastError;
    v.push_back(std::move(i));
  }
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
  if (j->status == JobStatus::Running || j->status == JobStatus::Ready)
    j->status = JobStatus::Paused;
  return true;
}

bool JobScheduler::resume(int id) {
  Job *j = find(id);
  if (!j)
    return false;
  if (j->status == JobStatus::Paused)
    j->status = JobStatus::Ready;
  return true;
}

bool JobScheduler::interrupt(int id) {
  Job *j = find(id);
  if (!j || j->threadId == 0)
    return false;
  PyGILState_STATE gil = PyGILState_Ensure();
  const int n = PyThreadState_SetAsyncExc(j->threadId, PyExc_KeyboardInterrupt);
  PyGILState_Release(gil);
  return n == 1;
}

bool JobScheduler::kill(int id) {
  for (auto it = m_jobs.begin(); it != m_jobs.end(); ++it) {
    if ((*it)->id != id)
      continue;
    Job *j = it->get();
    PyGILState_STATE gil = PyGILState_Ensure();
    Py_XDECREF(j->stepFn);
    Py_XDECREF(j->module);
    PyGILState_Release(gil);
    m_jobs.erase(it);
    return true;
  }
  return false;
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
