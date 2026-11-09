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

#ifndef ISOSURFACEDIALOG_H
#define ISOSURFACEDIALOG_H

#include <QDialog>
#include <boost/signals2.hpp>
#include <memory>
#include <string>
#include <vector>

class QComboBox;
class QDoubleSpinBox;
class QSpinBox;
class QCheckBox;
class QPushButton;
class QProgressBar;
class QLabel;
namespace cvc { namespace gl { class SceneGraph; } }
namespace cvc {
class app;
}

class IsosurfaceDialog : public QDialog {
  Q_OBJECT

public:
  explicit IsosurfaceDialog(cvc::app &app, std::shared_ptr<cvc::gl::SceneGraph> sceneGraph,
                            QWidget *parent = nullptr);
  ~IsosurfaceDialog() = default;

private slots:
  void onVolumeSelected(int index);
  void onComputeClicked();
  void onGraphicsChildrenChanged();

private:
  void setupUI();
  void connectSignals();
  void populateVolumeList();
  void updateProgress(int value);
  void onComputeFinished(bool success, const std::string &message);
  void setControlsEnabled(bool enabled);

  cvc::app &m_app;
  std::shared_ptr<cvc::gl::SceneGraph> m_sceneGraph;

  // UI elements
  QComboBox *m_volumeComboBox;
  QDoubleSpinBox *m_isovalueSpinBox;
  QComboBox *m_methodComboBox;
  QSpinBox *m_improveIterationsSpinBox;
  QComboBox *m_normalTypeComboBox;
  QPushButton *m_computeButton;
  QPushButton *m_cancelButton;
  QProgressBar *m_progressBar;
  QLabel *m_statusLabel;

  // Volume tracking
  std::vector<std::string> m_volumePaths; // Full state tree paths

  // Computation state
  bool m_computing;
  std::string m_activeThreadKey;

  // Signal connections
  boost::signals2::scoped_connection m_graphicsChangedConnection;
};

#endif // ISOSURFACEDIALOG_H
