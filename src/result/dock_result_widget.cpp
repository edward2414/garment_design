//
//  Filename         : dock_result_widget.cpp
//  Author           : Emmanuel Turquin
//  Purpose          : The main widget of the dock for the Pretreatment module.
//  Date of creation : 04/18/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include <qcheckbox.h>
#include <qcolordialog.h>
#include <qspinbox.h>
#include <qlineedit.h>
#include "floatspinbox.h"
#include "config_result.h"
#include "config.h"
#include "dock_window.h"
#include "dock_result_widget.h"

DockResultWidget::DockResultWidget(DockWindow* parent, const char* name, WFlags fl) :
  DockResultWidgetBase(parent, name, fl)
{
  _dw = parent;
  updateConfig();
  show();
}

DockResultWidget::~DockResultWidget()
{

}

void DockResultWidget::updateConfig()
{
  bool b1 = config::RESULT_DISP_MODEL;
  bool b2 = config::RESULT_DISP_SURFACE;
  bool b3 = config::RESULT_DISP_CONTOUR;
  bool b4 = config::RESULT_DISP_TYPE;
  bool b5 = false;
  if (_dw && _dw->config()) {
    _dw->config()->io()->getValue("display/model", b1);
    _dw->config()->io()->getValue("display/surface", b2);
    _dw->config()->io()->getValue("display/contour", b3);
    _dw->config()->io()->getValue("display/type", b4);
    _dw->config()->io()->getValue("display/surfaceRenderType", b5);
  }

  // Block signals in order not to trigger any callback.
  modelCheckBox->blockSignals(true);
  modelCheckBox->setChecked(b1);
  modelCheckBox->blockSignals(false);
  
  surfaceCheckBox->blockSignals(true);
  surfaceCheckBox->setChecked(b2);
  surfaceCheckBox->blockSignals(false);

  contourCheckBox->blockSignals(true);
  contourCheckBox->setChecked(b3);
  contourCheckBox->blockSignals(false);

  typeCheckBox->blockSignals(true);
  typeCheckBox->setChecked(b4);
  typeCheckBox->blockSignals(false);
  
  surfaceRenderTypeCheckBox->blockSignals(true);
  surfaceRenderTypeCheckBox->setChecked(b5);
  surfaceRenderTypeCheckBox->blockSignals(false);
}

void DockResultWidget::modelCheckBox_toggled(bool b)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/model", b);
  _dw->configUpdated();
}

void DockResultWidget::surfaceCheckBox_toggled(bool b)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/surface", b);
  _dw->configUpdated();
}

void DockResultWidget::surfaceRenderTypeCheckBox_toggled(bool b)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/surfaceRenderType", b);
  _dw->configUpdated();
}

void DockResultWidget::contourCheckBox_toggled(bool b)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/contour", b);
  _dw->configUpdated();
}

void DockResultWidget::typeCheckBox_toggled(bool b)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/type", b);
  _dw->config()->io()->setValue("display/type_changed", true);
  _dw->configUpdated();
}

void DockResultWidget::colorDialogButton_clicked()
{
  
  if (!_dw || !_dw->config())
    return;
  QColor col = QColorDialog::getColor();
  _dw->config()->io()->setValue("display/render/red", col.red());
  _dw->config()->io()->setValue("display/render/green", col.green());
  _dw->config()->io()->setValue("display/render/blue", col.blue());
  _dw->config()->io()->setValue("display/type_changed", true); // forces surfaces to be re-rendered
  
  _dw->configUpdated();
}

