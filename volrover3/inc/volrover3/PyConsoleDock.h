#ifndef VOLROVER3_PYCONSOLEDOCK_H
#define VOLROVER3_PYCONSOLEDOCK_H

// --------------------------------------------------------------------
// volrover3::PyConsoleDock — the embedded-Python console dock.
// --------------------------------------------------------------------
// A QDockWidget with two tabs (docs/EMBEDDED_PYTHON.md §12.4), cloned
// structurally from StateDashboardWidget:
//   - REPL tab: an output pane + an input line, driving
//     EmbeddedInterpreter::run_string_capture (captured stdout/stderr, bare
//     expressions echo their repr). Up/Down input history.
//   - Jobs tab: a table of the JobScheduler's jobs (id/name/status/steps),
//     polled ~1 Hz, with Interrupt / Stop acting on the selected row.
// Runs REPL input on the UI thread (the run_string_capture contract).
// --------------------------------------------------------------------

#include <QDockWidget>
#include <QString>
#include <QStringList>

class QPlainTextEdit;
class QLineEdit;
class QTableWidget;
class QTimer;
class QCheckBox;

namespace volrover3 {

class EmbeddedInterpreter;
class JobScheduler;

class PyConsoleDock : public QDockWidget {
  Q_OBJECT
public:
  PyConsoleDock(EmbeddedInterpreter *interp, JobScheduler *sched, QWidget *parent = nullptr);
  ~PyConsoleDock() override;

  // Run one REPL input: echo it, run under capture, append stdout (normal) +
  // stderr/traceback (error). Records history. Also the input line target.
  // Public for tests / programmatic use.
  void evaluate(const QString &source);
  QString outputText() const; // REPL output pane contents (tests)
  int jobRowCount() const;    // rows currently shown in the jobs table (tests)

  // Where "Load Script…" starts browsing (Settings::scriptsDir). MainWindow sets
  // this; falls back to $HOME when empty.
  void setScriptsDir(const QString &dir);

  // Load `path`, submit it to the scheduler as a job named after the file, and
  // report the outcome in the REPL pane. `onWorker` runs it on a hard-killable
  // sacrificial thread. Returns the new job id, or -1 on read/load failure.
  // Public so tests can drive it without a file dialog.
  int runScriptFile(const QString &path, bool onWorker = false);

public slots:
  void refreshJobs(); // repopulate the jobs table from the scheduler

protected:
  bool eventFilter(QObject *obj, QEvent *ev) override; // Up/Down REPL history

private slots:
  void onReplReturn();
  void onInterruptClicked();
  void onStopClicked();
  void onLoadScriptClicked(); // pick a .py and submit it as a job

private:
  int selectedJobId() const;

  EmbeddedInterpreter *m_interp;
  JobScheduler *m_sched;
  QPlainTextEdit *m_output = nullptr;
  QLineEdit *m_input = nullptr;
  QTableWidget *m_jobs = nullptr;
  QCheckBox *m_workerCheck = nullptr;
  QTimer *m_refresh = nullptr;
  QStringList m_history;
  int m_historyPos = 0;
  QString m_scriptsDir;
};

} // namespace volrover3

#endif // VOLROVER3_PYCONSOLEDOCK_H
