#ifndef INPUTSTATE_H
#define INPUTSTATE_H

#include <set>
#include <string>

namespace cvc {
class app;
class state;
} // namespace cvc

class QKeyEvent;
class QMouseEvent;
class QWheelEvent;

/// Publishes mouse and keyboard state into the cvc state tree.
///
/// Rationale: the state tree is already the app's cross-cutting bus — the
/// camera syncs through it and Python reads/writes it via `pycvc.state_set` —
/// so putting input there lets scripts and nodes react to the mouse and
/// keyboard without anyone wrapping Qt's event API. A demo can poll
/// `volrover3.input.mouse.x`, or connect to a node's `valueChanged`, and never
/// touch C++.
///
/// Lives under `volrover3.input.*`, deliberately NOT under
/// `volrover3.camera.input.*` where CameraController already publishes three
/// button flags. Anything written beneath `volrover3.camera` synchronously
/// re-enters `CameraController::handleStateChanged`, which walks a linear
/// if/else chain; at 60+ mouse-move events per second with several writes each
/// that is hundreds of no-op dispatches per second. A sibling subtree costs
/// nothing.
class InputState {
public:
  InputState(cvc::app &app, const std::string &statePrefix = "volrover3.input");

  /// Seed every key with a default so the subtree is fully discoverable
  /// (in the state-tree browser, and to Python) before any input arrives.
  void initializeDefaults();

  void handleMouseMove(int x, int y, int dx, int dy, int buttons, int modifiers);
  void handleMousePress(int button, int x, int y, int buttons, int modifiers);
  void handleMouseRelease(int button, int x, int y, int buttons, int modifiers);
  void handleWheel(int deltaX, int deltaY, int x, int y, int modifiers);
  void handleKeyPress(int key, int modifiers, const std::string &text);
  void handleKeyRelease(int key, int modifiers);

  /// Clear transient/held state — call when the widget loses focus, otherwise a
  /// key held while alt-tabbing stays "down" forever.
  void clearHeld();

  cvc::state &getState(const std::string &path);

private:
  void setButtons(int buttons);
  void setModifiers(int modifiers);
  void publishKeysHeld();

  cvc::app &m_app;
  std::string m_statePrefix;
  std::set<int> m_keysHeld;
  /// Wheel is a discrete event; the tree holds an accumulator plus a
  /// monotonically increasing counter so a poller can detect "something
  /// happened" without racing on the delta itself.
  int m_wheelAccumY = 0;
  int m_wheelAccumX = 0;
  unsigned long long m_eventSeq = 0;
};

#endif // INPUTSTATE_H
