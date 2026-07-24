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
#include <QStringList>

class QPlainTextEdit;
class QLineEdit;
class QTableWidget;
class QTimer;

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

public slots:
  void refreshJobs(); // repopulate the jobs table from the scheduler

protected:
  bool eventFilter(QObject *obj, QEvent *ev) override; // Up/Down REPL history

private slots:
  void onReplReturn();
  void onInterruptClicked();
  void onStopClicked();

private:
  int selectedJobId() const;

  EmbeddedInterpreter *m_interp;
  JobScheduler *m_sched;
  QPlainTextEdit *m_output = nullptr;
  QLineEdit *m_input = nullptr;
  QTableWidget *m_jobs = nullptr;
  QTimer *m_refresh = nullptr;
  QStringList m_history;
  int m_historyPos = 0;
};

} // namespace volrover3

#endif // VOLROVER3_PYCONSOLEDOCK_H
