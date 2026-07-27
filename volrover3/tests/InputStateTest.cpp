#include <QtCore/Qt>
#include <cvc/core/app.h>
#include <cvc/core/state.h>
#include <gtest/gtest.h>
#include <sstream>
#include <volrover3/InputState.h>

class InputStateTest : public ::testing::Test {
protected:
  cvc::app ctx;

  void SetUp() override {
    // Unique subtree per test so state does not leak between cases.
    std::stringstream ss;
    ss << "volrover3.input.test" << testCounter++;
    prefix = ss.str();
    input = new InputState(ctx, prefix);
  }

  void TearDown() override { delete input; }

  cvc::state &node(const std::string &path) {
    return cvc::state::instance(ctx)(prefix)(path);
  }

  static int testCounter;
  std::string prefix;
  InputState *input = nullptr;
};

int InputStateTest::testCounter = 0;

TEST_F(InputStateTest, DefaultsAreSeededSoTheSubtreeIsDiscoverable) {
  // Consumers should be able to browse/poll every key before any input lands.
  EXPECT_EQ(node("mouse.x").value<int>(), 0);
  EXPECT_EQ(node("mouse.y").value<int>(), 0);
  EXPECT_FALSE(node("mouse.left").value<bool>());
  EXPECT_EQ(node("key.held_count").value<int>(), 0);
  EXPECT_EQ(node("modifiers").value<int>(), 0);
  EXPECT_EQ(node("event_seq").value<int>(), 0);
}

TEST_F(InputStateTest, MouseMovePublishesPositionDeltaAndButtons) {
  input->handleMouseMove(120, 45, 4, -3, Qt::LeftButton, Qt::ShiftModifier);

  EXPECT_EQ(node("mouse.x").value<int>(), 120);
  EXPECT_EQ(node("mouse.y").value<int>(), 45);
  EXPECT_EQ(node("mouse.dx").value<int>(), 4);
  EXPECT_EQ(node("mouse.dy").value<int>(), -3);
  EXPECT_TRUE(node("mouse.left").value<bool>());
  EXPECT_FALSE(node("mouse.right").value<bool>());
  EXPECT_TRUE(node("mouse.inside").value<bool>());
  EXPECT_TRUE(node("modifiers.shift").value<bool>());
  EXPECT_FALSE(node("modifiers.ctrl").value<bool>());
}

TEST_F(InputStateTest, ButtonFlagsTrackTheFullBitmask) {
  input->handleMousePress(Qt::RightButton, 1, 2, Qt::RightButton | Qt::MiddleButton, 0);

  EXPECT_TRUE(node("mouse.right").value<bool>());
  EXPECT_TRUE(node("mouse.middle").value<bool>());
  EXPECT_FALSE(node("mouse.left").value<bool>());

  input->handleMouseRelease(Qt::MiddleButton, 1, 2, Qt::RightButton, 0);
  EXPECT_TRUE(node("mouse.right").value<bool>());
  EXPECT_FALSE(node("mouse.middle").value<bool>());
}

TEST_F(InputStateTest, HeldKeysAccumulateAndClear) {
  input->handleKeyPress(Qt::Key_W, 0, "w");
  input->handleKeyPress(Qt::Key_A, 0, "a");

  EXPECT_EQ(node("key.held_count").value<int>(), 2);
  EXPECT_EQ(node("key.last_pressed").value<int>(), static_cast<int>(Qt::Key_A));
  EXPECT_EQ(node("key.last_text").value(), "a");

  input->handleKeyRelease(Qt::Key_W, 0);
  EXPECT_EQ(node("key.held_count").value<int>(), 1);
  EXPECT_EQ(node("key.last_released").value<int>(), static_cast<int>(Qt::Key_W));

  // Held list is a comma-separated set; only Key_A should remain.
  EXPECT_EQ(node("key.held").value(), std::to_string(static_cast<int>(Qt::Key_A)));
}

TEST_F(InputStateTest, ClearHeldDropsEverything) {
  // A key held while the window loses focus must not stay down forever.
  input->handleKeyPress(Qt::Key_Space, Qt::ControlModifier, " ");
  input->handleMousePress(Qt::LeftButton, 0, 0, Qt::LeftButton, Qt::ControlModifier);

  input->clearHeld();

  EXPECT_EQ(node("key.held_count").value<int>(), 0);
  EXPECT_EQ(node("key.held").value(), "");
  EXPECT_FALSE(node("mouse.left").value<bool>());
  EXPECT_EQ(node("modifiers").value<int>(), 0);
  EXPECT_FALSE(node("mouse.inside").value<bool>());
}

TEST_F(InputStateTest, WheelAccumulatesAcrossEvents) {
  input->handleWheel(0, 120, 10, 20, 0);
  input->handleWheel(0, 120, 10, 20, 0);
  input->handleWheel(0, -40, 10, 20, 0);

  EXPECT_EQ(node("wheel.dy").value<int>(), -40);
  EXPECT_EQ(node("wheel.accum_y").value<int>(), 200);
}

TEST_F(InputStateTest, EventSequenceAdvancesOnEveryEvent) {
  // Pollers use this to detect activity without diffing every field.
  const int start = node("event_seq").value<int>();
  input->handleMouseMove(1, 1, 1, 1, 0, 0);
  input->handleKeyPress(Qt::Key_X, 0, "x");
  input->handleWheel(0, 1, 0, 0, 0);

  EXPECT_EQ(node("event_seq").value<int>(), start + 3);
}

TEST_F(InputStateTest, StateIsReadableThroughTheSharedTreeRoot) {
  // Python reaches these via pycvc.state_set/-get on the same tree, so the
  // values must be visible from a plain rooted lookup, not just our handle.
  input->handleMouseMove(7, 9, 0, 0, 0, 0);

  auto &root = cvc::state::instance(ctx);
  EXPECT_EQ(root(prefix + ".mouse.x").value<int>(), 7);
  EXPECT_EQ(root(prefix + ".mouse.y").value<int>(), 9);
}
