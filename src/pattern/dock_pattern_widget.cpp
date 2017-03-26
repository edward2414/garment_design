//
//  Filename         : dock_pattern_widget.cpp
//  Author           : Jamie Wither
//  Purpose          : The main widget of the dock for the Pattern module.
//  Date of creation : 13 Dec 2005
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include <qcheckbox.h>
#include <qspinbox.h>
#include <qlineedit.h>
#include <qlistview.h>
#include <qslider.h>
#include "floatspinbox.h"
#include "config_pattern.h"
#include "config.h"
#include "dock_window.h"
#include "dock_pattern_widget.h"

DockPatternWidget::DockPatternWidget(DockWindow* parent, const char* name, WFlags fl) :
  DockPatternWidgetBase(parent, name, fl)
{
  _dw = parent;
  
  _dw->config()->io()->setValue("display/pattern", true); // always start in pattern mode
  
  listView1->setColumnText( 0, "Atlas name" );
  listView1->addColumn( "Rows" );
  listView1->addColumn( "Cols" );
  listView1->clear();
  listView1->insertItem ( new QListViewItem( listView1, "Load pattern file" ) );
  
  updateConfig();
  show();
}

DockPatternWidget::~DockPatternWidget()
{

}

void DockPatternWidget::updateConfig()
{
  bool b1 = config::PATTERN_DISP_PATTERN;
  float twist_degree_val = 90.0;
  
  bool bmMesh = true;
  bool bcMesh = false;
  int atlas_cols = 2;
  int atlas_rows = 2;
  bool bSurface = true;
  bool twistEnabled = false;
  float ptf_val = 2.0/3.0;
  float hcompress_top_val = 1.0;
  float hcompress_bottom_val = 1.0;
  float twist_factor_val = 4.0;
  bool collision_enabled =false;
  bool display_model = true;
  int top_edges = 1;
  int bottom_edges = 1;
  int left_edges = 1;
  int right_edges = 1;
  bool e1,e2,e3,e4;
  e1 = e2 = e3 = e4 = false;
  int slide1 = 0;
  int slide2 = 0;
  
  
  if (_dw && _dw->config()) {
    _dw->config()->io()->getValue("display/pattern", b1);
    _dw->config()->io()->getValue("display/twist_degree_val", twist_degree_val);
    _dw->config()->io()->getValue("display/model_mesh", bmMesh);
    _dw->config()->io()->getValue("display/control_mesh", bcMesh);
    _dw->config()->io()->getValue("display/atlas_cols", atlas_cols);
    _dw->config()->io()->getValue("display/atlas_rows", atlas_rows);
    _dw->config()->io()->getValue("display/buckling_surface", bSurface);
    _dw->config()->io()->getValue("display/twist_enabled", twistEnabled);
    _dw->config()->io()->getValue("display/ptf_val", ptf_val);
    _dw->config()->io()->getValue("display/hcompress_top_val", hcompress_top_val);
    _dw->config()->io()->getValue("display/hcompress_bottom_val", hcompress_bottom_val);
    _dw->config()->io()->getValue("display/twist_factor_val", twist_factor_val);
    _dw->config()->io()->getValue("display/collision_detection_enabled", collision_enabled);
    _dw->config()->io()->getValue("display/model", display_model);
    _dw->config()->io()->getValue("display/top_edges", top_edges);
    _dw->config()->io()->getValue("display/bottom_edges", bottom_edges);
    _dw->config()->io()->getValue("display/left_edges", left_edges);
    _dw->config()->io()->getValue("display/right_edges", right_edges);
    _dw->config()->io()->getValue("display/effect1", e1);
    _dw->config()->io()->getValue("display/effect2", e2);
    _dw->config()->io()->getValue("display/effect3", e3);
    _dw->config()->io()->getValue("display/effect4", e4);
    _dw->config()->io()->getValue("display/slider1", slide1);
    _dw->config()->io()->getValue("display/slider2", slide2);
    
  }

  // Block signals in order not to trigger any callback.
  patternCheckBox->blockSignals(true);
  patternCheckBox->setChecked(b1);
  patternCheckBox->blockSignals(false);
  
  twistSpinBox->blockSignals(true);
  twistSpinBox->setValue(twist_degree_val);
  twistSpinBox->blockSignals(false);
  
  twistFactorFloatSpinBox->blockSignals(true);
  twistFactorFloatSpinBox->setValue(twist_factor_val);
  twistFactorFloatSpinBox->blockSignals(false);
  
  hCompressTopFloatSpinBox->blockSignals(true);
  hCompressTopFloatSpinBox->setValue(hcompress_top_val);
  hCompressTopFloatSpinBox->setLineStep(0.1);
  hCompressTopFloatSpinBox->blockSignals(false);
  
  hCompressBottomFloatSpinBox->blockSignals(true);
  hCompressBottomFloatSpinBox->setValue(hcompress_bottom_val);
  hCompressBottomFloatSpinBox->setLineStep(0.1);
  hCompressBottomFloatSpinBox->blockSignals(false);
  
  patchTangentFactorFloatSpinBox->blockSignals(true);
  patchTangentFactorFloatSpinBox->setValue(ptf_val);
  patchTangentFactorFloatSpinBox->setLineStep(0.05);
  patchTangentFactorFloatSpinBox->blockSignals(false);
  
  atlasColsSpinBox->blockSignals(true);
  atlasColsSpinBox->setValue(atlas_cols);
  atlasColsSpinBox->blockSignals(false);
  
  atlasRowsSpinBox->blockSignals(true);
  atlasRowsSpinBox->setValue(atlas_rows);
  atlasRowsSpinBox->blockSignals(false);
    
  leftEdgesSpinBox->blockSignals(true);
  leftEdgesSpinBox->setValue(left_edges);
  leftEdgesSpinBox->blockSignals(false);
    
  rightEdgesSpinBox->blockSignals(true);
  rightEdgesSpinBox->setValue(right_edges);
  rightEdgesSpinBox->blockSignals(false);
  
  topEdgesSpinBox->blockSignals(true);
  topEdgesSpinBox->setValue(top_edges);
  topEdgesSpinBox->blockSignals(false);
  
  bottomEdgesSpinBox->blockSignals(true);
  bottomEdgesSpinBox->setValue(bottom_edges);
  bottomEdgesSpinBox->blockSignals(false);
  
  modelMeshCheckBox->blockSignals(true);
  modelMeshCheckBox->setChecked(bmMesh);
  modelMeshCheckBox->blockSignals(false);
  
  controlMeshCheckBox->blockSignals(true);
  controlMeshCheckBox->setChecked(bcMesh);
  controlMeshCheckBox->blockSignals(false);
  
  
  bSurfaceCheckBox->blockSignals(true);
  bSurfaceCheckBox->setChecked(bSurface);
  bSurfaceCheckBox->blockSignals(false);
  
  twistEnabledCheckBox->blockSignals(true);
  twistEnabledCheckBox->setChecked(twistEnabled);
  twistEnabledCheckBox->blockSignals(false);
  
  collisionEnabledCheckBox->blockSignals(true);
  collisionEnabledCheckBox->setChecked(collision_enabled);
  collisionEnabledCheckBox->blockSignals(false);

  displayModelCheckBox->blockSignals(true);
  displayModelCheckBox->setChecked(display_model);
  displayModelCheckBox->blockSignals(false);
  
  effect1CheckBox->blockSignals(true);
  effect1CheckBox->setChecked(e1);
  effect1CheckBox->blockSignals(false);
  
  effect2CheckBox->blockSignals(true);
  effect2CheckBox->setChecked(e2);
  effect2CheckBox->blockSignals(false);
  
  effect3CheckBox->blockSignals(true);
  effect3CheckBox->setChecked(e3);
  effect3CheckBox->blockSignals(false);
  
  effect4CheckBox->blockSignals(true);
  effect4CheckBox->setChecked(e4);
  effect4CheckBox->blockSignals(false);
  
  slider1->blockSignals(true);
  slider1->setValue(slide1);
  slider1->blockSignals(false);
      
  slider2->blockSignals(true);
  slider2->setValue(slide2);
  slider2->blockSignals(false);
  
}

void DockPatternWidget::patternCheckBox_toggled(bool b)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/pattern", b);
  _dw->config()->io()->setValue("display/type_changed", true);
  _dw->configUpdated();
}

void DockPatternWidget::switchAtlas(QListViewItem* pItem)
{
  //(Pattern_dw->choose
  // FIXME tied in strings
  if (!_dw || !_dw->config())
    return;
  if(!pItem) return;
  _dw->config()->io()->setValue("display/atlas_name", pItem->text(0) );
  _dw->config()->io()->setValue("display/atlas_changed", true);
  _dw->configUpdated();
}

void DockPatternWidget::controlMeshCheckBox_toggled(bool b)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/control_mesh", b);
  _dw->config()->io()->setValue("display/redraw", true);
  _dw->configUpdated();
}

void DockPatternWidget::Redraw_clicked()
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/redraw", true);
  _dw->configUpdated();
}

void DockPatternWidget::twistSpinBox_valueChanged(float v)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/twist_degree_val", v);
}

void DockPatternWidget::twistFactorFloatSpinBox_valueChanged(float v)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/twist_factor_val", v);
}

void DockPatternWidget::hCompressTopFloatSpinBox_valueChanged(float v)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/hcompress_top_val", v);
}

void DockPatternWidget::hCompressBottomFloatSpinBox_valueChanged(float v)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/hcompress_bottom_val", v);
}

void DockPatternWidget::patchTangentFactorFloatSpinBox_valueChanged(float v)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/ptf_val", v);
}

void DockPatternWidget::modelMeshCheckBox_toggled(bool b)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/model_mesh", b);
  _dw->config()->io()->setValue("display/redraw", true);
  _dw->configUpdated();
}

void DockPatternWidget::atlasColsSpinBox_valueChanged(int val)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/atlas_cols", val);
  _dw->config()->io()->setValue("display/atlas_cols_changed", true);
  
  _dw->configUpdated();
}

void DockPatternWidget::atlasRowsSpinBox_valueChanged(int val)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/atlas_rows", val);
  _dw->config()->io()->setValue("display/atlas_rows_changed", true);
  
  _dw->configUpdated();
}

void DockPatternWidget::bSurfaceCheckBox_toggled(bool b)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/buckling_surface", b);
  _dw->config()->io()->setValue("display/redraw", true);
  _dw->configUpdated();
}

void DockPatternWidget::twistEnabledCheckBox_toggled(bool b)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/twist_enabled", b);
  _dw->config()->io()->setValue("display/redraw", true);
  _dw->configUpdated();
}

void DockPatternWidget::collisionEnabledCheckBox_toggled(bool b)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/collision_detection_enabled", b);
  _dw->config()->io()->setValue("display/redraw", true);
  _dw->configUpdated();
}

void DockPatternWidget::displayModelCheckBox_toggled(bool b)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/model", b);
  _dw->config()->io()->setValue("display/redraw", true);
  _dw->configUpdated();
}

void DockPatternWidget::axisAlignButton_clicked()
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/axis_align_clicked", true);
  _dw->configUpdated();
}

void DockPatternWidget::topEdgesSpinBox_valueChanged(int val)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/top_edges", val);
  _dw->config()->io()->setValue("display/edges_required_changed", true);
  
  _dw->configUpdated();
}

void DockPatternWidget::bottomEdgesSpinBox_valueChanged(int val)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/bottom_edges", val);
  _dw->config()->io()->setValue("display/edges_required_changed", true);
  
  _dw->configUpdated();
}

void DockPatternWidget::leftEdgesSpinBox_valueChanged(int val)
{  
  if (!_dw || !_dw->config())
  return;
  _dw->config()->io()->setValue("display/left_edges", val);
  _dw->config()->io()->setValue("display/edges_required_changed", true);
  
  _dw->configUpdated();
}

void DockPatternWidget::rightEdgesSpinBox_valueChanged(int val)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/right_edges", val);
  _dw->config()->io()->setValue("display/edges_required_changed", true);
  
  _dw->configUpdated();
}

void DockPatternWidget::effect1CheckBox_toggled(bool b)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/effect1", b);
  _dw->config()->io()->setValue("display/redraw", true);
  _dw->configUpdated();
}

void DockPatternWidget::effect2CheckBox_toggled(bool b)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/effect2", b);
  _dw->config()->io()->setValue("display/redraw", true);
  _dw->configUpdated();
}

void DockPatternWidget::effect3CheckBox_toggled(bool b)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/effect3", b);
  _dw->config()->io()->setValue("display/redraw", true);
  _dw->configUpdated();
}

void DockPatternWidget::effect4CheckBox_toggled(bool b)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/effect4", b);
  _dw->config()->io()->setValue("display/redraw", true);
  _dw->configUpdated();
}

void DockPatternWidget::slider1_valueChanged(int i)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/slider1", i);
  _dw->config()->io()->setValue("display/redraw", true);
  
  _dw->configUpdated();
}

void DockPatternWidget::slider2_valueChanged(int i)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/slider2", i);
  _dw->config()->io()->setValue("display/redraw", true);
  
  _dw->configUpdated();
}


























