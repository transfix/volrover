#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QDockWidget>
#include <QLabel>
#include <QMainWindow>
#include <QProgressBar>
#include <QToolBar>
#include <boost/signals2.hpp>
#include <cvc/core/app.h>
#include <memory>
#include <vector>

class AppState;
class VTKRenderWidget;
class TransferFunctionWidget;
class SceneGraph;
class ThreadMonitorWidget;
class StateTreeWidget;
class StateDashboardWidget;
class GridOptionsDialog;
class SDFDialog;
class IsosurfaceDialog;
class GeometryDialog;
class VolumeDialog;
class ViewerOptionsDialog;
class CameraSettingsDialog;
namespace volrover3 {
class EmbeddedInterpreter;
class Settings;
class JobScheduler;
class PyConsoleDock;
}

class MainWindow : public QMainWindow {
  Q_OBJECT

public:
  explicit MainWindow(std::shared_ptr<cvc::app> app, QWidget *parent = nullptr);
  ~MainWindow();

  cvc::app &app() const { return *m_app; }
  AppState &appState() const { return *m_appState; }

  // Startup automation (--exec-script / --screenshot; see main.cpp). Run a Python
  // file in the embedded interpreter, and render the current scene to a PNG.
  void execStartupScript(const QString &path);
  bool saveScreenshot(const QString &path);

private slots:
  void openFile();
  void toggleGrid();
  void toggleAxis();
  void editBoundingBox();
  void editCameraSettings();
  void showGridOptions();
  void showViewerOptions();
  void showThreadMonitor();
  void showStateTree();
  void showStateDashboard();
  void showSDF();
  void showIsosurface();
  void showGeometry();
  void showVolume();
  void aboutVolRover();
  void updateThreadStatus();
  void resetCamera();
  void generateStanfordBunny();
  void generateSphere();
  void generateCube();
  void generateTorus();
  void generateCone();

protected:
  void closeEvent(QCloseEvent *event) override;

private:
  void createMenus();
  void createToolBar();
  void createDockWidgets();
  void setupConnections();
  void initializeCameraFromState();
  void setupStatusBar();

  // Owned app/state, declared first so they outlive the Qt child-widget
  // members below (which hold cvc::app&/AppState& references into them).
  std::shared_ptr<cvc::app> m_app;
  std::unique_ptr<AppState> m_appState;
  // Per-instance settings, backed by this instance's cvc::state section and
  // loaded from ~/.volrover before the interpreter (feeds its python home/mode).
  std::unique_ptr<volrover3::Settings> m_settings;

  VTKRenderWidget *m_renderWidget;
  TransferFunctionWidget *m_transferFunctionWidget;
  std::shared_ptr<SceneGraph> m_sceneGraph;
  // Embedded Python interpreter — built right after m_sceneGraph so it captures
  // the live app + scene; owns CPython's lifecycle + the injected PyHost.
  std::unique_ptr<volrover3::EmbeddedInterpreter> m_interp;
  // The job scheduler (owns the QTimer tick + job registry) and the console dock
  // that drives it + the REPL. m_consoleDock is Qt-owned (added via addDockWidget).
  std::unique_ptr<volrover3::JobScheduler> m_scheduler;
  volrover3::PyConsoleDock *m_consoleDock = nullptr;
  ThreadMonitorWidget *m_threadMonitor;
  StateTreeWidget *m_stateTreeWidget;
  StateDashboardWidget *m_stateDashboardWidget;
  GridOptionsDialog *m_gridOptionsDialog;
  SDFDialog *m_sdfDialog;
  IsosurfaceDialog *m_isosurfaceDialog;
  GeometryDialog *m_geometryDialog;
  VolumeDialog *m_volumeDialog;
  ViewerOptionsDialog *m_viewerOptionsDialog;
  CameraSettingsDialog *m_cameraDialog;

  // Toolbar
  QToolBar *m_mainToolBar;

  // Status bar widgets for thread monitoring
  QLabel *m_threadNameLabel;
  QLabel *m_threadInfoLabel;
  QProgressBar *m_threadProgressBar;

  std::vector<boost::signals2::connection> m_connections;

  bool m_gridVisible;
  bool m_axisVisible;
};

#endif // MAINWINDOW_H
