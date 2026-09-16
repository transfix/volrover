#ifndef GRAPHICSPARENTDIALOG_H
#define GRAPHICSPARENTDIALOG_H

#include <QDialog>
#include <memory>
#include <string>

class QComboBox;
class QPushButton;
namespace cvc { namespace gl {
class GraphicsNode;
class VolumeNode;
class SceneGraph;
} } // namespace cvc::gl

/**
 * @brief Dialog for selecting a parent graphics node for new geometry or volume
 *
 * Allows the user to select which graphics node should be the parent
 * for newly loaded geometry or volume files. Shows a hierarchical list of existing
 * graphics nodes and volume graphics nodes, with the root as the default option.
 */
class GraphicsParentDialog : public QDialog {
  Q_OBJECT

public:
  explicit GraphicsParentDialog(std::shared_ptr<cvc::gl::SceneGraph> sceneGraph, QWidget *parent = nullptr);
  ~GraphicsParentDialog() override;

  // Get the selected parent node name (empty string = root)
  std::string getSelectedParentName() const;

  // Get the selected parent node (nullptr = root)
  std::shared_ptr<cvc::gl::GraphicsNode> getSelectedParent() const;

  // Get the selected volume parent node (nullptr = root)
  std::shared_ptr<cvc::gl::VolumeNode> getSelectedVolumeParent() const;

private:
  void populateParentList();
  void addNodeToList(std::shared_ptr<cvc::gl::GraphicsNode> node, int depth = 0);
  void addVolumeNodeToList(std::shared_ptr<cvc::gl::GraphicsNode> node, int depth = 0);

  std::shared_ptr<cvc::gl::SceneGraph> m_sceneGraph;
  QComboBox *m_parentComboBox;
  QPushButton *m_okButton;
  QPushButton *m_cancelButton;
};

#endif // GRAPHICSPARENTDIALOG_H
