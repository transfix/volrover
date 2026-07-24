#ifndef VOLROVER3_SETTINGS_H
#define VOLROVER3_SETTINGS_H

// --------------------------------------------------------------------
// volrover3::Settings — per-instance settings, BACKED BY the cvc::state tree.
// --------------------------------------------------------------------
// Each volrover3 instance owns a dedicated section of the global state root
// (default prefix "volrover3", the same one the SceneGraph + AppState run
// under). Its settings live under <prefix>.settings.* in cvc::state, alongside
// the instance's camera/graphics state — so a Python script
// (`pycvc.state_set(app,"volrover3.settings.interpreter.mode","multi")`) and the
// C++ settings UI read/write the SAME live, reactive values.
//
// Persistence: the <prefix>.settings.* subtree round-trips to
// ~/.volrover/settings.yaml (yaml-cpp, human-editable); ~/.volrover/data.db
// (Qt6Sql/WAL) holds bulkier/arbitrary persisted state (console history, job
// audit, KV). A plain object (not a singleton); construct BEFORE the
// EmbeddedInterpreter so mode()/pythonHome() feed its Config.
// See docs/EMBEDDED_PYTHON.md §12.3.
// --------------------------------------------------------------------

#include <cvc/core/app.h>

#include <QSqlDatabase>
#include <memory>
#include <string>

namespace cvc {
class state;
}

namespace volrover3 {

enum class InterpreterMode { Single, Multi };

class Settings {
public:
  // `app` = the instance's owned app; `instancePrefix` = its state-root section.
  explicit Settings(std::shared_ptr<cvc::app> app, const std::string &instancePrefix = "volrover3");
  ~Settings();

  Settings(const Settings &) = delete;
  Settings &operator=(const Settings &) = delete;

  // Create ~/.volrover if missing, seed defaults, parse settings.yaml into the
  // state subtree, open data.db (WAL, CREATE TABLE IF NOT EXISTS). Idempotent.
  void load();
  // Emit the <prefix>.settings.* subtree back to ~/.volrover/settings.yaml.
  void save() const;

  // -- typed accessors over <prefix>.settings.* in cvc::state --------------
  InterpreterMode mode() const;
  void setMode(InterpreterMode m); // persisted; applies on next restart
  // "" -> caller falls back to $VOLROVER3_PYTHON_HOME / CPython default.
  std::string pythonHome() const;
  std::string scriptsDir() const; // default ~/.volrover/scripts
  int tickMs() const;             // scheduler tick, default 100
  int historySize() const;        // REPL history cap, default 500

  const std::string &instancePrefix() const { return m_prefix; }
  const std::string &homeDir() const { return m_home; } // ~/.volrover

  // -- arbitrary persisted state in data.db --------------------------------
  QSqlDatabase &db() { return m_db; }
  std::string kvGet(const std::string &key, const std::string &dflt = "") const;
  void kvPut(const std::string &key, const std::string &value);

private:
  cvc::state &node(const std::string &path) const; // <prefix>.settings.<path>
  std::string get(const std::string &path, const std::string &dflt) const;
  void seedDefaults();

  std::shared_ptr<cvc::app> m_app;
  std::string m_prefix;
  std::string m_home; // ~/.volrover
  QSqlDatabase m_db;
};

} // namespace volrover3

#endif // VOLROVER3_SETTINGS_H
