#include <QtCore/Qt>
#include <cvc/core/state.h>
#include <sstream>
#include <volrover3/InputState.h>

InputState::InputState(cvc::app &app, const std::string &statePrefix)
    : m_app(app), m_statePrefix(statePrefix) {
  initializeDefaults();
}

cvc::state &InputState::getState(const std::string &path) {
  return cvc::state::instance(m_app)(m_statePrefix)(path);
}

void InputState::initializeDefaults() {
  getState("mouse.x").value(0);
  getState("mouse.y").value(0);
  getState("mouse.dx").value(0);
  getState("mouse.dy").value(0);
  getState("mouse.left").value(false);
  getState("mouse.middle").value(false);
  getState("mouse.right").value(false);
  getState("mouse.buttons").value(0);
  getState("mouse.buttons").comment("Qt::MouseButtons bitmask of the buttons currently held");
  getState("mouse.inside").value(false);

  getState("wheel.dy").value(0);
  getState("wheel.dx").value(0);
  getState("wheel.accum_y").value(0);
  getState("wheel.accum_x").value(0);
  getState("wheel.accum_y").comment("Running total of vertical wheel delta since startup");

  getState("key.last_pressed").value(0);
  getState("key.last_released").value(0);
  getState("key.last_text").value(std::string());
  getState("key.held").value(std::string());
  getState("key.held").comment("Comma-separated Qt::Key codes currently held");
  getState("key.held_count").value(0);

  getState("modifiers").value(0);
  getState("modifiers").comment("Qt::KeyboardModifiers bitmask");
  getState("modifiers.shift").value(false);
  getState("modifiers.ctrl").value(false);
  getState("modifiers.alt").value(false);

  getState("event_seq").value(0);
  getState("event_seq").comment(
      "Increments on every input event; poll this to detect activity cheaply");
}

void InputState::setButtons(int buttons) {
  getState("mouse.buttons").value(buttons);
  getState("mouse.left").value((buttons & Qt::LeftButton) != 0);
  getState("mouse.middle").value((buttons & Qt::MiddleButton) != 0);
  getState("mouse.right").value((buttons & Qt::RightButton) != 0);
}

void InputState::setModifiers(int modifiers) {
  getState("modifiers").value(modifiers);
  getState("modifiers.shift").value((modifiers & Qt::ShiftModifier) != 0);
  getState("modifiers.ctrl").value((modifiers & Qt::ControlModifier) != 0);
  getState("modifiers.alt").value((modifiers & Qt::AltModifier) != 0);
}

void InputState::publishKeysHeld() {
  std::ostringstream oss;
  bool first = true;
  for (int key : m_keysHeld) {
    if (!first) {
      oss << ",";
    }
    oss << key;
    first = false;
  }
  getState("key.held").value(oss.str());
  getState("key.held_count").value(static_cast<int>(m_keysHeld.size()));
}

void InputState::handleMouseMove(int x, int y, int dx, int dy, int buttons, int modifiers) {
  getState("mouse.x").value(x);
  getState("mouse.y").value(y);
  getState("mouse.dx").value(dx);
  getState("mouse.dy").value(dy);
  getState("mouse.inside").value(true);
  setButtons(buttons);
  setModifiers(modifiers);
  getState("event_seq").value(static_cast<int>(++m_eventSeq));
}

void InputState::handleMousePress(int button, int x, int y, int buttons, int modifiers) {
  (void)button; // the full bitmask already carries which button went down
  getState("mouse.x").value(x);
  getState("mouse.y").value(y);
  setButtons(buttons);
  setModifiers(modifiers);
  getState("event_seq").value(static_cast<int>(++m_eventSeq));
}

void InputState::handleMouseRelease(int button, int x, int y, int buttons, int modifiers) {
  (void)button;
  getState("mouse.x").value(x);
  getState("mouse.y").value(y);
  setButtons(buttons);
  setModifiers(modifiers);
  getState("event_seq").value(static_cast<int>(++m_eventSeq));
}

void InputState::handleWheel(int deltaX, int deltaY, int x, int y, int modifiers) {
  m_wheelAccumX += deltaX;
  m_wheelAccumY += deltaY;
  getState("wheel.dx").value(deltaX);
  getState("wheel.dy").value(deltaY);
  getState("wheel.accum_x").value(m_wheelAccumX);
  getState("wheel.accum_y").value(m_wheelAccumY);
  getState("mouse.x").value(x);
  getState("mouse.y").value(y);
  setModifiers(modifiers);
  getState("event_seq").value(static_cast<int>(++m_eventSeq));
}

void InputState::handleKeyPress(int key, int modifiers, const std::string &text) {
  m_keysHeld.insert(key);
  getState("key.last_pressed").value(key);
  getState("key.last_text").value(text);
  publishKeysHeld();
  setModifiers(modifiers);
  getState("event_seq").value(static_cast<int>(++m_eventSeq));
}

void InputState::handleKeyRelease(int key, int modifiers) {
  m_keysHeld.erase(key);
  getState("key.last_released").value(key);
  publishKeysHeld();
  setModifiers(modifiers);
  getState("event_seq").value(static_cast<int>(++m_eventSeq));
}

void InputState::clearHeld() {
  m_keysHeld.clear();
  publishKeysHeld();
  setButtons(0);
  setModifiers(0);
  getState("mouse.inside").value(false);
}
