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

#ifndef BOUNDINGBOXDIALOG_H
#define BOUNDINGBOXDIALOG_H

#include <QCheckBox>
#include <QComboBox>
#include <QDialog>
#include <QDoubleSpinBox>
#include <QLineEdit>
#include <QPushButton>
#include <QSpinBox>
#include <cvc/volume/bounding_box.h>
#include <memory>
#include <vector>

namespace cvc { namespace gl {
class GraphicsNode;
class SceneGraph;
} } // namespace cvc::gl

class BoundingBoxDialog : public QDialog {
  Q_OBJECT

public:
  explicit BoundingBoxDialog(std::shared_ptr<cvc::gl::SceneGraph> sceneGraph, QWidget *parent = nullptr);

private slots:
  void onGraphicsSelectionChanged(int index);
  void onResetToGraphics();
  void onBBoxVisibilityChanged(bool visible);
  void onBBoxColorChanged();
  void onApplyChanges();

private:
  void setupUI();
  void populateGraphicsComboBox();
  void loadGraphicsSettings();
  void updateColorButton();

  std::shared_ptr<cvc::gl::SceneGraph> m_sceneGraph;
  std::vector<std::shared_ptr<cvc::gl::GraphicsNode>> m_graphicsList;
  std::shared_ptr<cvc::gl::GraphicsNode> m_currentGraphics;

  QComboBox *m_graphicsComboBox;
  QLineEdit *m_minXEdit;
  QLineEdit *m_minYEdit;
  QLineEdit *m_minZEdit;
  QLineEdit *m_maxXEdit;
  QLineEdit *m_maxYEdit;
  QLineEdit *m_maxZEdit;

  QCheckBox *m_bboxVisibleCheckbox;
  QPushButton *m_bboxColorButton;
  double m_bboxColor[3];
};

#endif // BOUNDINGBOXDIALOG_H
