#include <QApplication>
#include <QSurfaceFormat>
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

  // Own the process-wide cvc::app explicitly (no singleton). This shared_ptr
  // is created before MainWindow and destroyed after it, so every reference
  // threaded down into MainWindow/AppState/widgets/SceneGraph stays valid for
  // the whole run.
  auto cvcApp = std::make_shared<cvc::app>();

  MainWindow mainWindow(cvcApp);
  mainWindow.show();

  return app.exec();
}
