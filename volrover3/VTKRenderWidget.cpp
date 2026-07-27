#include <QKeyEvent>
#include <QMouseEvent>
#include <QWheelEvent>
#include <volrover3/AppState.h>
#include <volrover3/CameraController.h>
#include <volrover3/InputState.h>
#include <volrover3/JobScheduler.h>
#include <cvc/gl/SceneGraph.h>
#include <volrover3/VTKRenderWidget.h>
#include <vtkCamera.h>
#include <vtkCornerAnnotation.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkNew.h>
#include <vtkPNGWriter.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkTextProperty.h>
#include <vtkWindowToImageFilter.h>

VTKRenderWidget::VTKRenderWidget(cvc::app &app, AppState &appState, QWidget *parent)
    : QVTK_WIDGET_BASE(parent), m_app(app), m_appState(appState),
      m_renderWindow(vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New()),
      m_renderer(vtkSmartPointer<vtkRenderer>::New()),
      m_cameraController(std::make_unique<CameraController>(app)),
      m_inputState(std::make_unique<InputState>(app)),
      m_fpsAnnotation(vtkSmartPointer<vtkCornerAnnotation>::New()), m_showFPS(false) {
  initializeVTK();

  // Set up timer to process SceneGraph events on main thread. Its interval is the
  // render/event refresh cap from state ("volrover3.viewer.max_fps"); applyRenderRate
  // reads it and (re)starts the timer. Observe the state key so the rate is tunable
  // live from the state tree / Python console.
  connect(&m_eventTimer, &QTimer::timeout, this, &VTKRenderWidget::processSceneGraphEvents);
  applyRenderRate(); // start the timer at the state-configured rate
  m_maxFpsConn = m_appState.onMaxFPSChanged([this]() {
    // May fire from a job/script thread; hop to the GUI thread to touch the QTimer.
    QMetaObject::invokeMethod(this, "applyRenderRate", Qt::QueuedConnection);
  });

  // Set up timer to update FPS display (every 500ms)
  connect(&m_fpsTimer, &QTimer::timeout, this, &VTKRenderWidget::updateFPSDisplay);

  // Load FPS display setting from state
  m_showFPS = m_appState.showFPS();
  if (m_showFPS) {
    m_fpsAnnotation->SetVisibility(true);
    m_fpsTimer.start(500);
  }
}

VTKRenderWidget::~VTKRenderWidget() {}

void VTKRenderWidget::initializeVTK() {
  // Set up render window
  setRenderWindow(m_renderWindow);
  m_renderWindow->AddRenderer(m_renderer);

  // Set background color (dark gray)
  m_renderer->SetBackground(0.2, 0.2, 0.2);

  // Set up camera
  vtkCamera *camera = m_renderer->GetActiveCamera();
  camera->SetPosition(0, 0, 10);
  camera->SetFocalPoint(0, 0, 0);
  camera->SetViewUp(0, 1, 0);

  // Initialize camera controller
  m_cameraController->setCamera(camera);

  // Set up FPS annotation in top-left corner
  m_fpsAnnotation->SetText(2, "FPS: --");                      // Position 2 = upper left
  m_fpsAnnotation->GetTextProperty()->SetColor(1.0, 1.0, 0.0); // Yellow
  m_fpsAnnotation->GetTextProperty()->SetFontSize(14);
  m_fpsAnnotation->SetVisibility(false); // Hidden by default
  m_renderer->AddViewProp(m_fpsAnnotation);

  // Enable focus for keyboard input
  setFocusPolicy(Qt::StrongFocus);
  // Publish hover position too, not just drags: without tracking, Qt only
  // delivers mouseMoveEvent while a button is held, so volrover3.input.mouse.x/y
  // would be stale whenever the user is not dragging.
  setMouseTracking(true);
}

void VTKRenderWidget::setSceneGraph(std::shared_ptr<SceneGraph> sceneGraph) {
  m_sceneGraph = sceneGraph;
  if (m_sceneGraph) {
    m_sceneGraph->setRenderer(m_renderer);
  }
}

void VTKRenderWidget::keyPressEvent(QKeyEvent *event) {
  if (m_inputState && !event->isAutoRepeat()) {
    // Skip auto-repeat: X11 synthesises press/release pairs for a held key, so
    // publishing them would make key.held flicker for a key that never moved.
    m_inputState->handleKeyPress(event->key(), static_cast<int>(event->modifiers()),
                                 event->text().toStdString());
  }
  if (m_cameraController) {
    m_cameraController->handleKeyPress(event->key());
    updateCamera();
    render(); // render() re-fits the clip range (see render())
  }
  QVTK_WIDGET_BASE::keyPressEvent(event);
}

void VTKRenderWidget::keyReleaseEvent(QKeyEvent *event) {
  if (m_inputState && !event->isAutoRepeat()) {
    m_inputState->handleKeyRelease(event->key(), static_cast<int>(event->modifiers()));
  }
  if (m_cameraController) {
    m_cameraController->handleKeyRelease(event->key());
  }
  QVTK_WIDGET_BASE::keyReleaseEvent(event);
}

void VTKRenderWidget::mousePressEvent(QMouseEvent *event) {
  m_lastMousePos = event->pos();
  if (m_inputState) {
    m_inputState->handleMousePress(static_cast<int>(event->button()), event->pos().x(),
                                   event->pos().y(), static_cast<int>(event->buttons()),
                                   static_cast<int>(event->modifiers()));
  }
  if (m_cameraController) {
    m_cameraController->handleMousePress(event->button());
  }
  // Don't pass middle mouse to VTK to avoid conflicts with our
  // CameraController
  if (event->button() != Qt::MiddleButton) {
    QVTK_WIDGET_BASE::mousePressEvent(event);
  }
}

void VTKRenderWidget::mouseReleaseEvent(QMouseEvent *event) {
  if (m_inputState) {
    m_inputState->handleMouseRelease(static_cast<int>(event->button()), event->pos().x(),
                                     event->pos().y(), static_cast<int>(event->buttons()),
                                     static_cast<int>(event->modifiers()));
  }
  if (m_cameraController) {
    m_cameraController->handleMouseRelease(event->button());
  }
  // Don't pass middle mouse to VTK to avoid conflicts with our
  // CameraController
  if (event->button() != Qt::MiddleButton) {
    QVTK_WIDGET_BASE::mouseReleaseEvent(event);
  }
}

void VTKRenderWidget::mouseMoveEvent(QMouseEvent *event) {
  QPoint delta = event->pos() - m_lastMousePos;
  if (m_inputState) {
    m_inputState->handleMouseMove(event->pos().x(), event->pos().y(), delta.x(), delta.y(),
                                  static_cast<int>(event->buttons()),
                                  static_cast<int>(event->modifiers()));
  }
  if (m_cameraController) {
    m_cameraController->handleMouseMove(delta.x(), delta.y());
    updateCamera();
    render(); // render() re-fits the clip range (see render())
  }
  m_lastMousePos = event->pos();
  // Don't pass middle mouse moves to VTK when middle button is pressed
  bool isMiddlePressed = (event->buttons() & Qt::MiddleButton);
  if (!isMiddlePressed) {
    QVTK_WIDGET_BASE::mouseMoveEvent(event);
  }
}

void VTKRenderWidget::wheelEvent(QWheelEvent *event) {
  if (m_inputState) {
    m_inputState->handleWheel(event->angleDelta().x(), event->angleDelta().y(),
                              static_cast<int>(event->position().x()),
                              static_cast<int>(event->position().y()),
                              static_cast<int>(event->modifiers()));
  }
  if (m_cameraController) {
    m_cameraController->handleMouseWheel(event->angleDelta().y());
    updateCamera();
    render(); // render() re-fits the clip range (see render())
  }
  QVTK_WIDGET_BASE::wheelEvent(event);
}

void VTKRenderWidget::focusOutEvent(QFocusEvent *event) {
  if (m_inputState) {
    m_inputState->clearHeld();
  }
  QVTK_WIDGET_BASE::focusOutEvent(event);
}

void VTKRenderWidget::leaveEvent(QEvent *event) {
  if (m_inputState) {
    m_inputState->getState("mouse.inside").value(false);
  }
  QVTK_WIDGET_BASE::leaveEvent(event);
}

void VTKRenderWidget::updateCamera() {
  if (m_cameraController) {
    m_cameraController->update();
  }
}

void VTKRenderWidget::resetCamera() {
  if (m_renderer) {
    m_renderer->ResetCamera();
    renderWindow()->Render();
  }
}

void VTKRenderWidget::render() {
  if (m_renderWindow) {
    // Re-fit the near/far clipping planes to the scene every frame. The camera is
    // driven directly (SetPosition/SetFocalPoint via CameraController / the
    // volrover3.camera state path), and VTK does NOT auto-adjust the clip range on
    // a manual camera move — so a large scene (e.g. the ~3 km Austin bundle) gets
    // sliced by a stale, too-tight far plane. ResetCameraClippingRange() recomputes
    // near/far from the current actor bounds + camera pose; it's cheap for our
    // handful of props and is exactly what VTK's own interactors do each render.
    if (m_renderer) {
      m_renderer->ResetCameraClippingRange();
    }
    m_renderWindow->Render();
  }
}

bool VTKRenderWidget::saveScreenshot(const QString &path) {
  if (!m_sceneGraph) {
    return false;
  }
  // Render the live scene to a dedicated OFFSCREEN window (reliable even under
  // QT_QPA_PLATFORM=offscreen, where the QVTK widget's own window is 0x0). We
  // temporarily hand the scene graph this renderer, capture, then restore the
  // widget's renderer so interactive rendering keeps working.
  vtkNew<vtkRenderer> renderer;
  vtkNew<vtkRenderWindow> window;
  window->SetOffScreenRendering(1);
  window->AddRenderer(renderer);
  window->SetSize(1024, 768);
  m_sceneGraph->setRenderer(renderer); // attach node actors to the offscreen renderer
  m_sceneGraph->processEvents();
  renderer->ResetCamera();
  // Tilt from the default top-down view to an oblique aerial angle so 3D
  // structure (terrain relief, draped tracks, markers) reads in the capture.
  if (vtkCamera *cam = renderer->GetActiveCamera()) {
    cam->Elevation(-60.0);
    cam->Azimuth(30.0);
    cam->OrthogonalizeViewUp();
    renderer->ResetCamera(); // refit at the oblique angle
  }
  window->Render();

  vtkNew<vtkWindowToImageFilter> w2i;
  w2i->SetInput(window);
  w2i->Update();
  vtkNew<vtkPNGWriter> writer;
  writer->SetFileName(path.toUtf8().constData());
  writer->SetInputConnection(w2i->GetOutputPort());
  writer->Write();

  // Restore the live view's renderer.
  if (m_renderer) {
    m_sceneGraph->setRenderer(m_renderer);
    m_sceneGraph->processEvents();
  }
  return true;
}

void VTKRenderWidget::processSceneGraphEvents() {
  if (!m_sceneGraph)
    return;
  if (m_continuous && m_scheduler) {
    // Render-synced pump: advance the jobs by REAL elapsed dt once per frame
    // (JobScheduler::tick() derives dt from a steady clock), drain scene events,
    // then render EVERY frame at the widget's cadence — smooth motion decoupled
    // from the coarse self-tick.
    m_scheduler->tick();
    m_sceneGraph->processEvents();
    m_sceneGraph->checkAndResetRenderNeeded(); // consume the flag; we render anyway
    render();
    return;
  }
  m_sceneGraph->processEvents();
  // Trigger render if any events modified the scene
  if (m_sceneGraph->checkAndResetRenderNeeded()) {
    render();
  }
}

void VTKRenderWidget::setContinuousMode(volrover3::JobScheduler *scheduler) {
  m_scheduler = scheduler;
  m_continuous = (scheduler != nullptr);
}

void VTKRenderWidget::applyRenderRate() {
  // Clamp to a sane range: <1 FPS would stall event pumping; >240 wastes CPU and
  // out-runs any display. Interval is milliseconds per frame.
  double fps = m_appState.maxFPS();
  if (!(fps >= 1.0))
    fps = 1.0; // also catches NaN
  if (fps > 240.0)
    fps = 240.0;
  int intervalMs = static_cast<int>(1000.0 / fps + 0.5);
  if (intervalMs < 1)
    intervalMs = 1;
  m_eventTimer.start(intervalMs); // QTimer::start restarts with the new interval
}

void VTKRenderWidget::setShowFPS(bool show) {
  m_showFPS = show;

  if (m_fpsAnnotation) {
    m_fpsAnnotation->SetVisibility(show ? 1 : 0);
  }

  if (show) {
    // Start timer to update FPS display
    m_fpsTimer.start(500); // Update every 500ms
    updateFPSDisplay();
  } else {
    m_fpsTimer.stop();
  }

  // Trigger a render to show/hide the annotation
  if (m_renderWindow) {
    m_renderWindow->Render();
  }
}

void VTKRenderWidget::updateFPSDisplay() {
  if (!m_showFPS || !m_renderer || !m_fpsAnnotation) {
    return;
  }

  // Get the last render time from the renderer
  double renderTime = m_renderer->GetLastRenderTimeInSeconds();

  QString fpsText;
  if (renderTime > 0.0) {
    double fps = 1.0 / renderTime;
    fpsText = QString("FPS: %1").arg(fps, 0, 'f', 1);
  } else {
    fpsText = "FPS: --";
  }

  m_fpsAnnotation->SetText(2, fpsText.toStdString().c_str());

  // Trigger render to update the display
  if (m_renderWindow) {
    m_renderWindow->Render();
  }
}
