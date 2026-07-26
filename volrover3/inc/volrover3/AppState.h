#ifndef APPSTATE_H
#define APPSTATE_H

#include <boost/signals2/connection.hpp>
#include <cvc/core/app.h>
#include <cvc/core/state.h>
#include <cvc/geometry/geometry.h>
#include <cvc/volume/bounding_box.h>
#include <cvc/volume/volume.h>
#include <memory>

// Application state manager using cvc::state for reactive updates
class AppState {
public:
  // Create an instance bound to an explicitly-owned cvc::app. The state prefix
  // roots this viewer's state under app's state tree (custom prefixes allow
  // multiple viewers or testing).
  explicit AppState(cvc::app &app, const std::string &statePrefix = "volrover3");

  // Get the state prefix for this instance
  std::string getStatePrefix() const { return m_statePrefix; }

  // State accessors with change notification
  cvc::bounding_box worldBounds();
  void setWorldBounds(const cvc::bounding_box &bounds);

  // Camera control mode (0 = orbit, 1 = fly)
  int cameraMode();
  void setCameraMode(int mode);

  // Camera settings
  double cameraSpeed();
  void setCameraSpeed(double speed);

  double cameraSensitivity();
  void setCameraSensitivity(double sensitivity);

  bool cameraInvertMouse();
  void setCameraInvertMouse(bool invert);

  // Camera key bindings
  int cameraKeyForward();
  void setCameraKeyForward(int key);

  int cameraKeyBackward();
  void setCameraKeyBackward(int key);

  int cameraKeyLeft();
  void setCameraKeyLeft(int key);

  int cameraKeyRight();
  void setCameraKeyRight(int key);

  int cameraKeyUp();
  void setCameraKeyUp(int key);

  int cameraKeyDown();
  void setCameraKeyDown(int key);

  // Camera position and orientation
  void getCameraPosition(double &x, double &y, double &z);
  void setCameraPosition(double x, double y, double z);

  void getCameraViewDirection(double &x, double &y, double &z);
  void setCameraViewDirection(double x, double y, double z);

  void getCameraUpVector(double &x, double &y, double &z);
  void setCameraUpVector(double x, double y, double z);

  double cameraFieldOfView();
  void setCameraFieldOfView(double fov);

  // Viewer options
  bool showFPS();
  void setShowFPS(bool show);

  // Render refresh rate cap (frames/sec). Drives the render widget's frame/event
  // timer — in continuous mode this is the target animation frame rate; otherwise
  // it is how often queued scene events are pumped. Tunable live via the state
  // key "volrover3.viewer.max_fps". Clamped to a sane range by the widget.
  double maxFPS();
  void setMaxFPS(double fps);

  // Register callbacks for state changes
  // Returns a connection object that can be used to disconnect the callback
  boost::signals2::connection onWorldBoundsChanged(const boost::function<void()> &callback);
  boost::signals2::connection onCameraModeChanged(const boost::function<void()> &callback);
  boost::signals2::connection onCameraChanged(const boost::function<void()> &callback);
  boost::signals2::connection onMaxFPSChanged(const boost::function<void()> &callback);

  // State tree access for debugging/inspection
  cvc::state &getRootState();

public:
  ~AppState() = default;

private:
  AppState(const AppState &) = delete;
  AppState &operator=(const AppState &) = delete;

  cvc::state &getState(const std::string &path);
  void initializeDefaults();

  cvc::app &m_app;
  std::string m_statePrefix;
};

#endif // APPSTATE_H
