//
//  Filename         : dock_pretreatment_widget.cpp
//  Author           : Emmanuel Turquin
//  Purpose          : The main widget of the dock for the Pretreatment module.
//  Date of creation : 04/18/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include <qradiobutton.h>
#include <qcheckbox.h>
#include <qspinbox.h>
#include <qslider.h>
#include <qlineedit.h>
#include "floatspinbox.h"
#include "config_pretreatment.h"
#include "config.h"
#include "dock_window.h"
#include "dock_pretreatment_widget.h"

DockPretreatmentWidget::DockPretreatmentWidget(DockWindow* parent, const char* name, WFlags fl) :
  DockPretreatmentWidgetBase(parent, name, fl)
{
  _dw = parent;

  textureSizeSpinBox->setMaxValue(config::PRETREATMENT_TEXTURE_MAX_SIZE);

  dfXSpinBox->setMaxValue(512);
  dfYSpinBox->setMaxValue(512);
  dfZSpinBox->setMaxValue(512);

  dfXMinSpinBox->setMaxValue(511);
  dfYMinSpinBox->setMaxValue(511);
  dfZMinSpinBox->setMaxValue(511);
  dfXMaxSpinBox->setMaxValue(511);
  dfYMaxSpinBox->setMaxValue(511);
  dfZMaxSpinBox->setMaxValue(511);
  
  dfSliceSpinBox->setMaxValue(99);

  updateConfig();
  show();
}

DockPretreatmentWidget::~DockPretreatmentWidget()
{
}

void DockPretreatmentWidget::updateConfig()
{
  bool b1 = config::PRETREATMENT_CAMERA_PERSPECTIVE;
  bool b2 = config::PRETREATMENT_DISP_BBOX;
  bool b3 = config::PRETREATMENT_DISP_MODEL;
  bool b4 = config::PRETREATMENT_DISP_OCTREE;
  bool b5 = config::PRETREATMENT_DISP_AXIS;
  bool b6 = config::PRETREATMENT_DISP_DF;
  int level = config::PRETREATMENT_OCTREE_LEVEL;
  float bbox_x = config::PRETREATMENT_BBOX_X;
  float bbox_y = config::PRETREATMENT_BBOX_Y;
  float bbox_z = config::PRETREATMENT_BBOX_Z;
  float bbox_x_step = config::PRETREATMENT_BBOX_X;
  float bbox_y_step = config::PRETREATMENT_BBOX_Y;
  float bbox_z_step = config::PRETREATMENT_BBOX_Z;
  unsigned texture_size = config::PRETREATMENT_TEXTURE_SIZE;
  unsigned df_x = config::PRETREATMENT_DF_X;
  unsigned df_y = config::PRETREATMENT_DF_Y;
  unsigned df_z = config::PRETREATMENT_DF_Z;
  unsigned df_x_min = 0;
  unsigned df_y_min = 0;
  unsigned df_z_min = 0;
  unsigned df_x_max = config::PRETREATMENT_DF_X - 1;
  unsigned df_y_max = config::PRETREATMENT_DF_Y - 1;
  unsigned df_z_max = config::PRETREATMENT_DF_Z - 1;
  if (_dw && _dw->config()) {
    _dw->config()->io()->getValue("camera/perspective", b1);
    _dw->config()->io()->getValue("display/bbox", b2);
    _dw->config()->io()->getValue("display/model", b3);
    _dw->config()->io()->getValue("display/octree", b4);
    _dw->config()->io()->getValue("display/bbox_axis", b5);
    _dw->config()->io()->getValue("display/df", b6);
    _dw->config()->io()->getValue("octree/level", level);
    _dw->config()->io()->getValue("bbox/size/x", bbox_x);
    _dw->config()->io()->getValue("bbox/size/y", bbox_y);
    _dw->config()->io()->getValue("bbox/size/z", bbox_z);
    _dw->config()->io()->getValue("bbox/size/x_step", bbox_x_step);
    _dw->config()->io()->getValue("bbox/size/y_step", bbox_y_step);
    _dw->config()->io()->getValue("bbox/size/z_step", bbox_z_step);
    _dw->config()->io()->getValue("texture/size", texture_size);
    _dw->config()->io()->getValue("df/size/x", df_x);
    _dw->config()->io()->getValue("df/size/y", df_y);
    _dw->config()->io()->getValue("df/size/z", df_z);
    _dw->config()->io()->getValue("display/dfsize/x_min", df_x_min);
    _dw->config()->io()->getValue("display/dfsize/x_max", df_x_max);
    _dw->config()->io()->getValue("display/dfsize/y_min", df_y_min);
    _dw->config()->io()->getValue("display/dfsize/y_max", df_y_max);
    _dw->config()->io()->getValue("display/dfsize/z_min", df_z_min);
    _dw->config()->io()->getValue("display/dfsize/z_max", df_z_max);
  }

 // Block signals in order not to trigger any callback.
  perspectiveButton->blockSignals(true);
  perspectiveButton->setChecked(b1);
  perspectiveButton->blockSignals(false);

  orthoButton->blockSignals(true);
  orthoButton->setChecked(!b1);
  orthoButton->blockSignals(false);

  bbCheckBox->blockSignals(true);
  bbCheckBox->setChecked(b2);
  bbCheckBox->blockSignals(false);

  modelCheckBox->blockSignals(true);
  modelCheckBox->setChecked(b3);
  modelCheckBox->blockSignals(false);

  octreeCheckBox->blockSignals(true);
  octreeCheckBox->setChecked(b4);
  octreeCheckBox->blockSignals(false);

  bbAxisCheckBox->blockSignals(true);
  bbAxisCheckBox->setChecked(b5);
  bbAxisCheckBox->blockSignals(false);

  dfCheckBox->blockSignals(true);
  dfCheckBox->setChecked(b6);
  dfCheckBox->blockSignals(false);

  octreeLevelSpinBox->blockSignals(true);
  octreeLevelSpinBox->setValue(level);
  octreeLevelSpinBox->blockSignals(false);

  bbXSpinBox->blockSignals(true);
  bbXSpinBox->setValue(bbox_x);
  bbXSpinBox->setLineStep(bbox_x_step);
  bbXSpinBox->blockSignals(false);

  bbYSpinBox->blockSignals(true);
  bbYSpinBox->setValue(bbox_y);
  bbYSpinBox->setLineStep(bbox_y_step);
  bbYSpinBox->blockSignals(false);

  bbZSpinBox->blockSignals(true);
  bbZSpinBox->setValue(bbox_z);
  bbZSpinBox->setLineStep(bbox_z_step);
  bbZSpinBox->blockSignals(false);

  textureSizeSpinBox->blockSignals(true);
  textureSizeSpinBox->setValue(texture_size);
  textureSizeSpinBox->blockSignals(false);

  dfXSpinBox->blockSignals(true);
  dfXSpinBox->setValue(df_x);
  dfXSpinBox->blockSignals(false);

  dfYSpinBox->blockSignals(true);
  dfYSpinBox->setValue(df_y);
  dfYSpinBox->blockSignals(false);

  dfZSpinBox->blockSignals(true);
  dfZSpinBox->setValue(df_z);
  dfZSpinBox->blockSignals(false);

  dfXMinSpinBox->blockSignals(true);
  dfXMinSpinBox->setMaxValue(df_x-1);
  dfXMinSpinBox->setValue(df_x_min);
  dfXMinSpinBox->blockSignals(false);

  dfXMaxSpinBox->blockSignals(true);
  dfXMaxSpinBox->setMaxValue(df_x-1);
  dfXMaxSpinBox->setValue(df_x_max);
  dfXMaxSpinBox->blockSignals(false);

  dfYMinSpinBox->blockSignals(true);
  dfYMinSpinBox->setMaxValue(df_y-1);
  dfYMinSpinBox->setValue(df_y_min);
  dfYMinSpinBox->blockSignals(false);

  dfYMaxSpinBox->blockSignals(true);
  dfYMaxSpinBox->setMaxValue(df_y-1);
  dfYMaxSpinBox->setValue(df_y_max);
  dfYMaxSpinBox->blockSignals(false);

  dfZMinSpinBox->blockSignals(true);
  dfZMinSpinBox->setMaxValue(df_z-1);
  dfZMinSpinBox->setValue(df_z_min);
  dfZMinSpinBox->blockSignals(false);

  dfZMaxSpinBox->blockSignals(true);
  dfZMaxSpinBox->setMaxValue(df_z-1);
  dfZMaxSpinBox->setValue(df_z_max);
  dfZMaxSpinBox->blockSignals(false);
}

void DockPretreatmentWidget::perspectiveButton_toggled(bool b)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("camera/perspective", b);
  _dw->configUpdated();
}

void DockPretreatmentWidget::orthoButton_toggled(bool b)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("camera/perspective", !b);
  _dw->configUpdated();
}

void DockPretreatmentWidget::bbCheckBox_toggled(bool b)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/bbox", b);
  _dw->configUpdated();
}

void DockPretreatmentWidget::bbAxisCheckBox_toggled(bool b)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/bbox_axis", b);
  _dw->configUpdated();
}

void DockPretreatmentWidget::modelCheckBox_toggled(bool b)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/model", b);
  _dw->configUpdated();
}

void DockPretreatmentWidget::octreeCheckBox_toggled(bool b)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/octree", b);
  _dw->configUpdated();
}

void DockPretreatmentWidget::dfCheckBox_toggled(bool b)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/df", b);
  _dw->configUpdated();
}

void DockPretreatmentWidget::computeButton_clicked()
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("compute/start", true);
  _dw->configUpdated();
}

void DockPretreatmentWidget::octreeLevelSpinBox_valueChanged(int val)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("octree/level", val);
  _dw->configUpdated();
}

void DockPretreatmentWidget::bbResetButton_clicked()
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("bbox/reset", true);
  _dw->config()->io()->setValue("bbox/modified", true);
  _dw->configUpdated();
}

void DockPretreatmentWidget::bbXSpinBox_valueChanged(float val)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("bbox/size/x", val);
  _dw->config()->io()->setValue("bbox/modified", true);
  _dw->configUpdated();
}

void DockPretreatmentWidget::bbYSpinBox_valueChanged(float val)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("bbox/size/y", val);
  _dw->config()->io()->setValue("bbox/modified", true);
  _dw->configUpdated();
}

void DockPretreatmentWidget::bbZSpinBox_valueChanged(float val)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("bbox/size/z", val);
  _dw->config()->io()->setValue("bbox/modified", true);
  _dw->configUpdated();
}

void DockPretreatmentWidget::textureSizeSpinBox_valueChanged(int val)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("texture/size", val);
  _dw->configUpdated();
}

void DockPretreatmentWidget::dfXSpinBox_valueChanged(int val)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("df/size/x", val);
  dfXMinSpinBox->blockSignals(true);
  if (dfXMinSpinBox->value() >= val)
    _dw->config()->io()->setValue("display/dfsize/x_min", val-1);
  dfXMinSpinBox->setMaxValue(val-1);
  dfXMinSpinBox->blockSignals(false);
  dfXMaxSpinBox->blockSignals(true);
  if (dfXMaxSpinBox->value() >= val)
    _dw->config()->io()->setValue("display/dfsize/x_max", val-1);
  dfXMaxSpinBox->setMaxValue(val-1);
  dfXMaxSpinBox->blockSignals(false);
  _dw->configUpdated();
}

void DockPretreatmentWidget::dfYSpinBox_valueChanged(int val)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("df/size/y", val);
  dfYMinSpinBox->blockSignals(true);
  if (dfYMinSpinBox->value() >= val)
    _dw->config()->io()->setValue("display/dfsize/y_min", val-1);
  dfYMinSpinBox->setMaxValue(val-1);
  dfYMinSpinBox->blockSignals(false);
  dfYMaxSpinBox->blockSignals(true);
  if (dfYMaxSpinBox->value() >= val)
    _dw->config()->io()->setValue("display/dfsize/y_max", val-1);
  dfYMaxSpinBox->setMaxValue(val-1);
  dfYMaxSpinBox->blockSignals(false);
  _dw->configUpdated();
}

void DockPretreatmentWidget::dfZSpinBox_valueChanged(int val)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("df/size/z", val);
  dfZMinSpinBox->blockSignals(true);
  if (dfZMinSpinBox->value() >= val)
    _dw->config()->io()->setValue("display/dfsize/z_min", val-1);
  dfZMinSpinBox->setMaxValue(val-1);
  dfZMinSpinBox->blockSignals(false);
  dfZMaxSpinBox->blockSignals(true);
  if (dfZMaxSpinBox->value() >= val)
    _dw->config()->io()->setValue("display/dfsize/z_max", val-1);
  dfZMaxSpinBox->setMaxValue(val-1);
  dfZMaxSpinBox->blockSignals(false);
  _dw->configUpdated();
}

void DockPretreatmentWidget::dfXMinSpinBox_valueChanged(int val)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/dfsize/x_min", val);
  dfXMaxSpinBox->blockSignals(true);
  if (dfXMaxSpinBox->value() < val) {
    dfXMaxSpinBox->setValue(val);
    _dw->config()->io()->setValue("display/dfsize/x_max", val);
  }
  dfXMaxSpinBox->blockSignals(false);
  _dw->config()->io()->setValue("display/dfsize/modified", true);
  _dw->configUpdated();
}

void DockPretreatmentWidget::dfXMaxSpinBox_valueChanged(int val)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/dfsize/x_max", val);
  dfXMinSpinBox->blockSignals(true);
  if (dfXMinSpinBox->value() > val) {
    dfXMinSpinBox->setValue(val);
    _dw->config()->io()->setValue("display/dfsize/x_min", val);
  }
  dfXMinSpinBox->blockSignals(false);
  _dw->config()->io()->setValue("display/dfsize/modified", true);
  _dw->configUpdated();
}

void DockPretreatmentWidget::dfYMinSpinBox_valueChanged(int val)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/dfsize/y_min", val);
  dfYMaxSpinBox->blockSignals(true);
  if (dfYMaxSpinBox->value() < val) {
    dfYMaxSpinBox->setValue(val);
    _dw->config()->io()->setValue("display/dfsize/y_max", val);
  }
  dfYMaxSpinBox->blockSignals(false);
  _dw->config()->io()->setValue("display/dfsize/modified", true);
  _dw->configUpdated();
}

void DockPretreatmentWidget::dfYMaxSpinBox_valueChanged(int val)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/dfsize/y_max", val);
  dfYMinSpinBox->blockSignals(true);
  if (dfYMinSpinBox->value() > val) {
    dfYMinSpinBox->setValue(val);
    _dw->config()->io()->setValue("display/dfsize/y_min", val);
  }
  dfYMinSpinBox->blockSignals(false);
  _dw->config()->io()->setValue("display/dfsize/modified", true);
  _dw->configUpdated();
}

void DockPretreatmentWidget::dfZMinSpinBox_valueChanged(int val)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/dfsize/z_min", val);
  dfZMaxSpinBox->blockSignals(true);
  if (dfZMaxSpinBox->value() < val) {
    dfZMaxSpinBox->setValue(val);
    _dw->config()->io()->setValue("display/dfsize/z_max", val);
  }
  dfZMaxSpinBox->blockSignals(false);
  _dw->config()->io()->setValue("display/dfsize/modified", true);
  _dw->configUpdated();
}

void DockPretreatmentWidget::dfZMaxSpinBox_valueChanged(int val)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/dfsize/z_max", val);
  dfZMinSpinBox->blockSignals(true);
  if (dfZMinSpinBox->value() > val) {
    dfZMinSpinBox->setValue(val);
    _dw->config()->io()->setValue("display/dfsize/z_min", val);
  }
  dfZMinSpinBox->blockSignals(false);
  _dw->config()->io()->setValue("display/dfsize/modified", true);
  _dw->configUpdated();
}

void DockPretreatmentWidget::dfSliceSlider_valueChanged(int val)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/dfslice/percent_z", val);
  
  dfSliceSlider->blockSignals(true);
  dfSliceSpinBox->blockSignals(true);
  dfSliceSpinBox->setValue(val);
  dfSliceSlider->blockSignals(false);
  dfSliceSpinBox->blockSignals(false);

  _dw->config()->io()->setValue("bbox/modified", true);

  _dw->configUpdated();
}

void DockPretreatmentWidget::dfSliceSpinBox_valueChanged(int val)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/dfslice/percent_z", val);
  
  dfSliceSlider->blockSignals(true);
  dfSliceSpinBox->blockSignals(true);
  dfSliceSlider->setValue(val);
  dfSliceSlider->blockSignals(false);
  dfSliceSpinBox->blockSignals(false);

  _dw->config()->io()->setValue("bbox/modified", true);

  _dw->configUpdated();
}


