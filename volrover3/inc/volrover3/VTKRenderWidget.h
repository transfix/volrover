#ifndef VTKRENDERWIDGET_H
#define VTKRENDERWIDGET_H

#include <QTimer>
#include <QWidget>
#include <boost/signals2/connection.hpp>
#include <memory>
#include <vtkSmartPointer.h>

// Qt5/Qt6 compatibility
#if QT_VERSION >= QT_VERSION_CHECK(6, 0, 0)
#include <QVTKOpenGLNativeWidget.h>
#define QVTK_WIDGET_BASE QVTKOpenGLNativeWidget
#else
#include <QVTKOpenGLWidget.h>
#define QVTK_WIDGET_BASE QVTKOpenGLWidget
#endif

namespace cvc {
class app;
}
class AppState;
class vtkRenderer;
class vtkRenderWindow;
class vtkGenericOpenGLRenderWindow;
class vtkCornerAnnotation;
class SceneGraph;
class CameraController;
class InputState;
namespace volrover3 {
class JobScheduler;
}

class VTKRenderWidget : public QVTK_WIDGET_BASE {
  Q_OBJECT

public:
  explicit VTKRenderWidget(cvc::app &app, AppState &appState, QWidget *parent = nullptr);
  ~VTKRenderWidget();

  void setSceneGraph(std::shared_ptr<SceneGraph> sceneGraph);
  void resetCamera();
  void render(); // Force an immediate render

  // Render the current scene and write it to a PNG (for --screenshot / headless
  // demo capture). Pumps queued scene events + frames the camera first. Returns
  // false if there is no render window.
  bool saveScreenshot(const QString &path);

  // Continuous max-framerate mode: drive `scheduler`'s tick() off this widget's
  // render clock (so job step(dt) advances by real elapsed dt each frame) and
  // render every frame — decouples motion smoothness from the coarse job tick.
  // Pass the scheduler to enable; the caller should stop() the scheduler's own
  // timer first so it isn't double-driven. Pass nullptr to disable.
  void setContinuousMode(volrover3::JobScheduler *scheduler);

  // FPS display control
  void setShowFPS(bool show);
  bool showFPS() const { return m_showFPS; }

  CameraController *getCameraController() { return m_cameraController.get(); }
  InputState *getInputState() { return m_inputState.get(); }

protected:
  void keyPressEvent(QKeyEvent *event) override;
  void keyReleaseEvent(QKeyEvent *event) override;
  void mousePressEvent(QMouseEvent *event) override;
  void mouseReleaseEvent(QMouseEvent *event) override;
  void mouseMoveEvent(QMouseEvent *event) override;
  void wheelEvent(QWheelEvent *event) override;
  // Held keys/buttons must be dropped on focus loss, or a key held while
  // alt-tabbing stays "down" in the state tree forever.
  void focusOutEvent(QFocusEvent *event) override;
  void leaveEvent(QEvent *event) override;

private slots:
  void processSceneGraphEvents();
  void updateFPSDisplay();
  // Re-read viewer.max_fps and restart the frame/event timer at the new rate.
  // Invoked on the GUI thread (queued) since the state change may fire from a
  // job/script thread while QTimer must be driven from the widget's thread.
  void applyRenderRate();

private:
  void initializeVTK();
  void updateCamera();

  cvc::app &m_app;
  AppState &m_appState;
  vtkSmartPointer<vtkGenericOpenGLRenderWindow> m_renderWindow;
  vtkSmartPointer<vtkRenderer> m_renderer;
  std::shared_ptr<SceneGraph> m_sceneGraph;
  std::unique_ptr<CameraController> m_cameraController;
  // Mirrors mouse/keyboard into volrover3.input.* so scripts and nodes can
  // react to input through the state tree instead of a wrapped Qt API.
  std::unique_ptr<InputState> m_inputState;
  QTimer m_eventTimer; // Timer for processing SceneGraph events
  // Continuous-mode driving (see setContinuousMode): when set, m_eventTimer ticks
  // this scheduler and renders every frame.
  volrover3::JobScheduler *m_scheduler = nullptr;
  bool m_continuous = false;
  // Live "volrover3.viewer.max_fps" observer -> restarts m_eventTimer (GUI thread).
  boost::signals2::scoped_connection m_maxFpsConn;

  // FPS display
  vtkSmartPointer<vtkCornerAnnotation> m_fpsAnnotation;
  QTimer m_fpsTimer; // Timer for updating FPS display
  bool m_showFPS;

  QPoint m_lastMousePos;
};

#endif // VTKRENDERWIDGET_H
