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
  QCommandLineOption jobOpt(
      "run-job",
      "Load a Python file as a JobScheduler job (defines step(dt)) — it appears in "
      "the Python Console Jobs tab and is ticked cooperatively, unlike --exec-script.",
      "file");
  parser.addOption(execOpt);
  parser.addOption(shotOpt);
  parser.addOption(exitOpt);
  parser.addOption(jobOpt);
  parser.process(app);

  // Own the process-wide cvc::app explicitly (no singleton). This shared_ptr
  // is created before MainWindow and destroyed after it, so every reference
  // threaded down into MainWindow/AppState/widgets/SceneGraph stays valid for
  // the whole run.
  auto cvcApp = std::make_shared<cvc::app>();

  MainWindow mainWindow(cvcApp);
  mainWindow.show();

  if (parser.isSet(execOpt) || parser.isSet(shotOpt) || parser.isSet(jobOpt)) {
    const QString script = parser.value(execOpt);
    const QString job = parser.value(jobOpt);
    const QString shot = parser.value(shotOpt);
    const bool exitAfter = parser.isSet(exitOpt);
    // Let the window + the render widget's event timer come up first, then run
    // the startup script and/or submit the job (which mutate the live scene),
    // give the scheduler several ticks to animate, then screenshot / quit.
    QTimer::singleShot(700, &mainWindow, [&mainWindow, script, job, shot, exitAfter]() {
      if (!script.isEmpty())
        mainWindow.execStartupScript(script);
      if (!job.isEmpty())
        mainWindow.runScriptAsJob(job, /*onWorker=*/false);
      // ~30 scheduler ticks (jobs run on the scheduler's QTimer) so the agent
      // visibly moves before the screenshot.
      QTimer::singleShot(3000, &mainWindow, [&mainWindow, shot, exitAfter]() {
        if (!shot.isEmpty())
          mainWindow.saveScreenshot(shot);
        if (exitAfter)
          QApplication::quit();
      });
    });
  }

  return app.exec();
}
