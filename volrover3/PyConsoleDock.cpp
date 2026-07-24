#include <volrover3/PyConsoleDock.h>

#include <volrover3/EmbeddedInterpreter.h>
#include <volrover3/JobScheduler.h>

#include <QHBoxLayout>
#include <QHeaderView>
#include <QKeyEvent>
#include <QLineEdit>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QTabWidget>
#include <QTableWidget>
#include <QTableWidgetItem>
#include <QTimer>
#include <QVBoxLayout>
#include <QWidget>

namespace volrover3 {

namespace {
const char *statusName(JobStatus s) {
  switch (s) {
  case JobStatus::Ready: return "ready";
  case JobStatus::Running: return "running";
  case JobStatus::Paused: return "paused";
  case JobStatus::Terminated: return "terminated";
  case JobStatus::Killed: return "killed";
  case JobStatus::Error: return "error";
  }
  return "?";
}
} // namespace

PyConsoleDock::PyConsoleDock(EmbeddedInterpreter *interp, JobScheduler *sched, QWidget *parent)
    : QDockWidget(tr("Python Console"), parent), m_interp(interp), m_sched(sched) {
  setObjectName("PyConsoleDock");
  auto *tabs = new QTabWidget(this);

  // ── REPL tab ──
  auto *replTab = new QWidget(tabs);
  auto *replLayout = new QVBoxLayout(replTab);
  m_output = new QPlainTextEdit(replTab);
  m_output->setReadOnly(true);
  {
    QFont mono("monospace");
    mono.setStyleHint(QFont::Monospace);
    m_output->setFont(mono);
  }
  m_input = new QLineEdit(replTab);
  m_input->setPlaceholderText(tr(">>> run Python (app scripting via vrhost/pycvc)"));
  m_input->installEventFilter(this); // Up/Down history
  replLayout->addWidget(m_output, 1);
  replLayout->addWidget(m_input);
  connect(m_input, &QLineEdit::returnPressed, this, &PyConsoleDock::onReplReturn);
  tabs->addTab(replTab, tr("REPL"));

  // ── Jobs tab ──
  auto *jobsTab = new QWidget(tabs);
  auto *jobsLayout = new QVBoxLayout(jobsTab);
  m_jobs = new QTableWidget(0, 5, jobsTab);
  m_jobs->setHorizontalHeaderLabels({tr("Id"), tr("Name"), tr("Status"), tr("Steps"), tr("Elapsed")});
  m_jobs->horizontalHeader()->setStretchLastSection(true);
  m_jobs->setSelectionBehavior(QAbstractItemView::SelectRows);
  m_jobs->setEditTriggers(QAbstractItemView::NoEditTriggers);
  auto *btnRow = new QHBoxLayout;
  auto *interruptBtn = new QPushButton(tr("Interrupt"), jobsTab);
  auto *stopBtn = new QPushButton(tr("Stop"), jobsTab);
  btnRow->addWidget(interruptBtn);
  btnRow->addWidget(stopBtn);
  btnRow->addStretch(1);
  jobsLayout->addWidget(m_jobs, 1);
  jobsLayout->addLayout(btnRow);
  connect(interruptBtn, &QPushButton::clicked, this, &PyConsoleDock::onInterruptClicked);
  connect(stopBtn, &QPushButton::clicked, this, &PyConsoleDock::onStopClicked);
  tabs->addTab(jobsTab, tr("Jobs"));

  setWidget(tabs);

  // Poll the jobs table ~1 Hz.
  m_refresh = new QTimer(this);
  m_refresh->setInterval(1000);
  connect(m_refresh, &QTimer::timeout, this, &PyConsoleDock::refreshJobs);
  m_refresh->start();
  refreshJobs();

  if (!m_interp || !m_interp->ok())
    m_output->appendPlainText(tr("[interpreter unavailable — scripting disabled]"));
  else if (m_sched && m_sched->mode() == InterpreterMode::Multi)
    m_output->appendPlainText(
        tr("[multi-interpreter mode: jobs are isolated pure-Python; app scripting "
           "via pycvc/vrhost is disabled — switch to single mode in ~/.volrover/settings.yaml]"));
}

PyConsoleDock::~PyConsoleDock() = default;

void PyConsoleDock::evaluate(const QString &source) {
  const QString src = source.trimmed();
  if (src.isEmpty())
    return;
  m_output->appendPlainText(">>> " + src);
  if (!m_interp || !m_interp->ok()) {
    m_output->appendPlainText("[interpreter unavailable]");
    return;
  }
  std::string out, err;
  m_interp->run_string_capture(src.toStdString(), out, err);
  if (!out.empty()) {
    QString s = QString::fromStdString(out);
    if (s.endsWith('\n'))
      s.chop(1);
    m_output->appendPlainText(s);
  }
  if (!err.empty()) {
    QString s = QString::fromStdString(err);
    if (s.endsWith('\n'))
      s.chop(1);
    m_output->appendPlainText(s);
  }
  m_history << src;
  m_historyPos = m_history.size();
}

QString PyConsoleDock::outputText() const { return m_output ? m_output->toPlainText() : QString(); }

int PyConsoleDock::jobRowCount() const { return m_jobs ? m_jobs->rowCount() : 0; }

void PyConsoleDock::onReplReturn() {
  const QString src = m_input->text();
  m_input->clear();
  evaluate(src);
}

void PyConsoleDock::refreshJobs() {
  if (!m_sched)
    return;
  const auto jobs = m_sched->listJobs();
  m_jobs->setRowCount(static_cast<int>(jobs.size()));
  for (int r = 0; r < static_cast<int>(jobs.size()); ++r) {
    const auto &j = jobs[r];
    m_jobs->setItem(r, 0, new QTableWidgetItem(QString::number(j.id)));
    m_jobs->setItem(r, 1, new QTableWidgetItem(QString::fromStdString(j.name)));
    m_jobs->setItem(r, 2, new QTableWidgetItem(QString::fromLatin1(statusName(j.status))));
    m_jobs->setItem(r, 3, new QTableWidgetItem(QString::number(qulonglong(j.steps))));
    m_jobs->setItem(r, 4, new QTableWidgetItem(QString::number(j.elapsed, 'f', 1) + "s"));
    // stash the id on the row for control actions
    m_jobs->item(r, 0)->setData(Qt::UserRole, j.id);
  }
}

int PyConsoleDock::selectedJobId() const {
  const auto sel = m_jobs->selectionModel()->selectedRows();
  if (sel.isEmpty())
    return -1;
  QTableWidgetItem *it = m_jobs->item(sel.first().row(), 0);
  return it ? it->data(Qt::UserRole).toInt() : -1;
}

void PyConsoleDock::onInterruptClicked() {
  const int id = selectedJobId();
  if (id >= 0 && m_sched)
    m_sched->interrupt(id);
}

void PyConsoleDock::onStopClicked() {
  const int id = selectedJobId();
  if (id >= 0 && m_sched) {
    m_sched->kill(id);
    refreshJobs();
  }
}

// Up/Down input history on the REPL line edit.
bool PyConsoleDock::eventFilter(QObject *obj, QEvent *ev) {
  if (obj == m_input && ev->type() == QEvent::KeyPress) {
    auto *ke = static_cast<QKeyEvent *>(ev);
    if (ke->key() == Qt::Key_Up && m_historyPos > 0) {
      m_input->setText(m_history.value(--m_historyPos));
      return true;
    }
    if (ke->key() == Qt::Key_Down) {
      if (m_historyPos < m_history.size() - 1)
        m_input->setText(m_history.value(++m_historyPos));
      else {
        m_historyPos = m_history.size();
        m_input->clear();
      }
      return true;
    }
  }
  return QDockWidget::eventFilter(obj, ev);
}

} // namespace volrover3
