// Live-scene bridge capstone: a script builds geometry with pycvc and adds it to
// the RUNNING scene via vrhost.scene() — and the SAME C++ SceneGraph the host
// owns gains the node. This proves the adopt-existing bridge:
//
//   C++ owns `app` + `scene` --> EmbeddedInterpreter injects BOTH as PyCapsules
//   ("cvc.app" + "cvc.scenegraph") --> vrhost.scene() feeds them to
//   pycvc_gl.scene_from_capsule, which builds a pycvc_gl.Scene that ADOPTS the
//   live SceneGraph (no new make_shared) --> s.add_geometry(...) mutates the SAME
//   graph, so the C++ side sees the node by name.
//
// If vrhost.scene() spun up a parallel scene (the pre-bridge behaviour), the C++
// scene->getGraphics(name) below would be null. Headless: builds the graph, never
// setRenderer()s a window — no display needed. Requires pycvc_gl (single mode).
//
// One interpreter for the whole suite (CPython won't re-Py_Initialize after
// finalize). Python home + module dir come from the CTest environment.

#include <cvc/core/app.h>
#include <cvc/gl/SceneGraph.h>
#include <gtest/gtest.h>
#include <volrover3/EmbeddedInterpreter.h>

#include <memory>
#include <string>

class VrHostSceneBridgeTest : public ::testing::Test {
protected:
  static std::shared_ptr<cvc::app> app;
  static std::shared_ptr<SceneGraph> scene;
  static std::unique_ptr<volrover3::EmbeddedInterpreter> interp;
  static void SetUpTestSuite() {
    app = std::make_shared<cvc::app>();
    // The live scene the host owns — same ctor MainWindow uses.
    scene = std::make_shared<SceneGraph>(*app, "volrover3");
    interp = std::make_unique<volrover3::EmbeddedInterpreter>(app, scene);
    ASSERT_TRUE(interp->ok());
  }
  static void TearDownTestSuite() {
    interp.reset();
    scene.reset();
    app.reset();
  }
};
std::shared_ptr<cvc::app> VrHostSceneBridgeTest::app;
std::shared_ptr<SceneGraph> VrHostSceneBridgeTest::scene;
std::unique_ptr<volrover3::EmbeddedInterpreter> VrHostSceneBridgeTest::interp;

// vrhost.scene() returns a pycvc_gl.Scene (the bridge is wired + importable).
TEST_F(VrHostSceneBridgeTest, SceneReturnsPycvcGlScene) {
  std::string out, err;
  ASSERT_TRUE(interp->run_string_capture("import vrhost; type(vrhost.scene()).__name__", out, err))
      << "vrhost.scene() raised (is pycvc_gl importable in single mode?)\n"
      << err;
  EXPECT_NE(out.find("Scene"), std::string::npos) << "out=[" << out << "] err=[" << err << "]";
}

// The capstone: geometry added via vrhost.scene() lands in the SAME live graph.
TEST_F(VrHostSceneBridgeTest, ScriptAddsToLiveScene) {
  ASSERT_FALSE(static_cast<bool>(scene->getGraphics("bridge_probe")))
      << "test precondition: node must not exist before the script runs";

  ASSERT_TRUE(interp->run_string("import vrhost, pycvc\n"
                                 "app = vrhost.app()\n"
                                 "g = pycvc.geometry(app)\n"
                                 "g.add_vertices([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0])\n"
                                 "g.add_triangles([0, 1, 2])\n"
                                 "s = vrhost.scene()\n"
                                 "s.add_geometry('bridge_probe', g)\n"
                                 "s.pump()\n"))
      << "the bridge script raised — see the traceback on stderr";

  // The node reached the host's OWN SceneGraph — vrhost.scene() adopted it, it did
  // not build a parallel scene.
  EXPECT_TRUE(static_cast<bool>(scene->getGraphics("bridge_probe")))
      << "add_geometry via vrhost.scene() did not reach the live SceneGraph "
         "(bridge adopted a separate graph?)";
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
