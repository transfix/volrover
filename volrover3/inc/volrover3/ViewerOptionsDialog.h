/*
  Copyright 2026 CyberPC Angel, LLC

  This file is part of VolumeRover.

  VolumeRover is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  VolumeRover is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.
*/

#ifndef VIEWEROPTIONSDIALOG_H
#define VIEWEROPTIONSDIALOG_H

#include <QWidget>
#include <boost/signals2.hpp>
#include <memory>

class QCheckBox;
class QComboBox;
class QPushButton;
class VTKRenderWidget;
namespace cvc { namespace gl { class SceneGraph; } }
class AppState;

class ViewerOptionsDialog : public QWidget {
  Q_OBJECT

public:
  explicit ViewerOptionsDialog(AppState &appState, VTKRenderWidget *renderWidget,
                               std::shared_ptr<cvc::gl::SceneGraph> sceneGraph, QWidget *parent = nullptr);
  ~ViewerOptionsDialog() override;

protected:
  void showEvent(QShowEvent *event) override;
  void closeEvent(QCloseEvent *event) override;

private slots:
  void onShowFPSChanged(bool checked);
  void onGraphicsRootChanged(int index);
  void onCameraChanged(int index);
  void refreshGraphicsRoots();
  void refreshCameras();

private:
  void setupUI();
  void connectSignals();
  void loadFromState();

  AppState &m_appState;
  VTKRenderWidget *m_renderWidget;
  std::shared_ptr<cvc::gl::SceneGraph> m_sceneGraph;

  // Display options
  QCheckBox *m_showFPSCheckBox;

  // Scene selection
  QComboBox *m_graphicsRootComboBox;
  QPushButton *m_refreshRootsButton;

  // Camera selection
  QComboBox *m_cameraComboBox;
  QPushButton *m_refreshCamerasButton;

  std::vector<boost::signals2::connection> m_connections;
};

#endif // VIEWEROPTIONSDIALOG_H
