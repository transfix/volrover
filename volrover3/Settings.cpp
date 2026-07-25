#include <volrover3/Settings.h>

#include <cvc/core/state.h>

#include <QDir>
#include <QSqlError>
#include <QSqlQuery>
#include <QVariant>

#include <yaml-cpp/yaml.h>

#include <fstream>

namespace volrover3 {

namespace {
// The fixed settings schema: dotted state path under <prefix>.settings, the two
// yaml keys it maps to, and its default. Round-tripped explicitly (a small,
// human-editable file — no generic tree walk needed).
struct Key {
  const char *path;   // state path under <prefix>.settings
  const char *y0;     // yaml top-level key
  const char *y1;     // yaml nested key
  const char *dflt;   // default value (string)
};
const Key kKeys[] = {
    {"interpreter.mode", "interpreter", "mode", "single"},
    {"interpreter.python_home", "interpreter", "python_home", ""},
    {"interpreter.gate_multi_imports", "interpreter", "gate_multi_imports", "true"},
    {"interpreter.scripts_dir", "interpreter", "scripts_dir", ""},
    {"console.history_size", "console", "history_size", "500"},
    {"scheduler.tick_ms", "scheduler", "tick_ms", "100"},
};
} // namespace

Settings::Settings(std::shared_ptr<cvc::app> app, const std::string &instancePrefix)
    : m_app(std::move(app)), m_prefix(instancePrefix) {
  m_home = (QDir::homePath() + "/.volrover").toStdString();
}

Settings::~Settings() {
  if (m_db.isOpen())
    m_db.close();
}

cvc::state &Settings::node(const std::string &path) const {
  return cvc::state::instance(*m_app)(m_prefix)("settings")(path);
}

std::string Settings::get(const std::string &path, const std::string &dflt) const {
  std::string v = node(path).value();
  return v.empty() ? dflt : v;
}

void Settings::seedDefaults() {
  for (const auto &k : kKeys) {
    if (node(k.path).value().empty())
      node(k.path).value(std::string(k.dflt));
  }
  // scripts_dir default is derived from the home dir.
  if (node("interpreter.scripts_dir").value().empty())
    node("interpreter.scripts_dir").value(m_home + "/scripts");
}

void Settings::load() {
  QDir().mkpath(QString::fromStdString(m_home));
  QDir().mkpath(QString::fromStdString(m_home) + "/scripts");

  const std::string yamlPath = m_home + "/settings.yaml";
  if (QFile::exists(QString::fromStdString(yamlPath))) {
    try {
      YAML::Node root = YAML::LoadFile(yamlPath);
      for (const auto &k : kKeys) {
        const YAML::Node n = root[k.y0][k.y1];
        if (n && n.IsScalar())
          node(k.path).value(n.as<std::string>());
      }
    } catch (const std::exception &) {
      // A malformed file falls back to defaults rather than aborting startup.
    }
  }
  seedDefaults();
  save(); // normalize the file (create it if it was missing / add new keys)

  // data.db — arbitrary persisted state.
  const QString dbPath = QString::fromStdString(m_home + "/data.db");
  m_db = QSqlDatabase::addDatabase("QSQLITE", QString::fromStdString("volrover3:" + m_prefix));
  m_db.setDatabaseName(dbPath);
  if (m_db.open()) {
    QSqlQuery q(m_db);
    q.exec("PRAGMA journal_mode=WAL");
    q.exec("CREATE TABLE IF NOT EXISTS kv (key TEXT PRIMARY KEY, value BLOB, updated_at INTEGER)");
    q.exec("CREATE TABLE IF NOT EXISTS console_history (id INTEGER PRIMARY KEY AUTOINCREMENT, ts INTEGER, source TEXT)");
    q.exec("CREATE TABLE IF NOT EXISTS job_runs (id INTEGER PRIMARY KEY AUTOINCREMENT, name TEXT, mode TEXT, started INTEGER, ended INTEGER, status TEXT, error TEXT)");
  }
}

void Settings::save() const {
  YAML::Node root;
  root["version"] = 1;
  for (const auto &k : kKeys)
    root[k.y0][k.y1] = get(k.path, k.dflt);
  std::ofstream out(m_home + "/settings.yaml");
  if (out)
    out << root << "\n";
}

InterpreterMode Settings::mode() const {
  return get("interpreter.mode", "single") == "multi" ? InterpreterMode::Multi : InterpreterMode::Single;
}

void Settings::setMode(InterpreterMode m) {
  node("interpreter.mode").value(std::string(m == InterpreterMode::Multi ? "multi" : "single"));
  save();
}

std::string Settings::pythonHome() const { return get("interpreter.python_home", ""); }

std::string Settings::scriptsDir() const { return get("interpreter.scripts_dir", m_home + "/scripts"); }

int Settings::tickMs() const {
  try {
    return std::stoi(get("scheduler.tick_ms", "100"));
  } catch (...) {
    return 100;
  }
}

int Settings::historySize() const {
  try {
    return std::stoi(get("console.history_size", "500"));
  } catch (...) {
    return 500;
  }
}

std::string Settings::kvGet(const std::string &key, const std::string &dflt) const {
  if (!m_db.isOpen())
    return dflt;
  QSqlQuery q(m_db);
  q.prepare("SELECT value FROM kv WHERE key = ?");
  q.addBindValue(QString::fromStdString(key));
  if (q.exec() && q.next())
    return q.value(0).toString().toStdString();
  return dflt;
}

void Settings::kvPut(const std::string &key, const std::string &value) {
  if (!m_db.isOpen())
    return;
  QSqlQuery q(m_db);
  q.prepare("INSERT INTO kv(key, value, updated_at) VALUES(?, ?, strftime('%s','now')) "
            "ON CONFLICT(key) DO UPDATE SET value=excluded.value, updated_at=excluded.updated_at");
  q.addBindValue(QString::fromStdString(key));
  q.addBindValue(QString::fromStdString(value));
  q.exec();
}

} // namespace volrover3
