#include <QApplication>
#include <QCommandLineParser>
#include <QSurfaceFormat>
#include <QTimer>
#include <cvc/core/app.h>
#include <memory>
#include <volrover3/MainWindow.h>
#include <vtkAutoInit.h>

// VTK module initialization
VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);
VTK_MODULE_INIT(vtkRenderingFreeType);
VTK_MODULE_INIT(vtkRenderingVolumeOpenGL2);

int main(int argc, char *argv[]) {
  // Set up OpenGL format
  QSurfaceFormat format;
  format.setDepthBufferSize(24);
  format.setStencilBufferSize(8);
  format.setVersion(3, 3);
  format.setProfile(QSurfaceFormat::CoreProfile);
  QSurfaceFormat::setDefaultFormat(format);

  QApplication app(argc, argv);
  app.setApplicationName("VolRover3");
  app.setApplicationVersion("3.0.0");
  app.setOrganizationName("CVC");

  // Startup automation, so a demo can be launched non-interactively (and for
  // headless/CI capture): run a Python script in the embedded interpreter after
  // the window is up, optionally render the scene to a PNG, optionally quit.
  //   volrover3 --exec-script demo.py                 # run the demo, keep the window
  //   volrover3 --exec-script demo.py --screenshot out.png --exit-after  # headless
  QCommandLineParser parser;
  parser.setApplicationDescription("VolRover3 — volumetric viewer with an embedded Python interpreter");
  parser.addHelpOption();
  parser.addVersionOption();
  QCommandLineOption execOpt("exec-script",
                             "Run a Python script in the embedded interpreter at startup.", "file");
  QCommandLineOption shotOpt("screenshot",
                             "Render the scene to a PNG after the startup script runs.", "file");
  QCommandLineOption exitOpt("exit-after", "Quit once the startup script/screenshot is done.");
  parser.addOption(execOpt);
  parser.addOption(shotOpt);
  parser.addOption(exitOpt);
  parser.process(app);

  // Own the process-wide cvc::app explicitly (no singleton). This shared_ptr
  // is created before MainWindow and destroyed after it, so every reference
  // threaded down into MainWindow/AppState/widgets/SceneGraph stays valid for
  // the whole run.
  auto cvcApp = std::make_shared<cvc::app>();

  MainWindow mainWindow(cvcApp);
  mainWindow.show();

  if (parser.isSet(execOpt) || parser.isSet(shotOpt)) {
    const QString script = parser.value(execOpt);
    const QString shot = parser.value(shotOpt);
    const bool exitAfter = parser.isSet(exitOpt);
    // Let the window + the render widget's event timer come up first, then run
    // the script (which mutates the live scene), give the timer a couple of
    // ticks to pump + render the new nodes, then screenshot / quit.
    QTimer::singleShot(700, &mainWindow, [&mainWindow, script, shot, exitAfter]() {
      if (!script.isEmpty())
        mainWindow.execStartupScript(script);
      QTimer::singleShot(900, &mainWindow, [&mainWindow, shot, exitAfter]() {
        if (!shot.isEmpty())
          mainWindow.saveScreenshot(shot);
        if (exitAfter)
          QApplication::quit();
      });
    });
  }

  return app.exec();
}
