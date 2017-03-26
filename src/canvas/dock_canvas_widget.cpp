//
//  Filename         : dock_canvas_widget.cpp
//  Author           : Emmanuel Turquin
//  Purpose          : The main widget of the dock for the Pretreatment module.
//  Date of creation : 04/18/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include <qcheckbox.h>
#include <qspinbox.h>
#include <qlineedit.h>
#include "floatspinbox.h"
#include "config_canvas.h"
#include "config.h"
#include "dock_window.h"
#include "dock_canvas_widget.h"

DockCanvasWidget::DockCanvasWidget(DockWindow* parent, const char* name, WFlags fl) :
  DockCanvasWidgetBase(parent, name, fl)
{
  _dw = parent;
  splittingSpinBox->setLineStep(0.1);
  splittingSpinBox->setRange(-1.0, 1.0);
  gmXSpinBox->setMaxValue(1024);
  gmYSpinBox->setMaxValue(1024);
  borderXSpinBox->setLineStep(0.1);
  borderYSpinBox->setLineStep(0.1);
  updateConfig();
  show();
}

DockCanvasWidget::~DockCanvasWidget()
{

}

void DockCanvasWidget::updateConfig()
{
  int sampling = config::CANVAS_SAMPLING_INT;
  int snapping = config::CANVAS_SNAPPING_INT;
  float threshold = config::CANVAS_SPLITTING_THRESHOLD;
  bool b1 = config::CANVAS_CHAIN_PTS;
  bool b2 = config::CANVAS_STROKE_PTS;
  bool b3 = config::CANVAS_TEXTURE;
  bool b4 = config::CANVAS_TYPE;
  bool b5 = config::CANVAS_ANTIALIASING;
  bool b6 = config::CANVAS_BBOX;
  int gmx = config::CANVAS_GM_X;
  int gmy = config::CANVAS_GM_Y;
  float borderx = config::CANVAS_BORDER_X;
  float bordery = config::CANVAS_BORDER_Y;
  if (_dw && _dw->config()) {
    _dw->config()->io()->getValue("intervals/sampling", sampling);
    _dw->config()->io()->getValue("intervals/snapping", snapping);
    _dw->config()->io()->getValue("splitting/threshold", threshold);
    _dw->config()->io()->getValue("display/chain_pts", b1);
    _dw->config()->io()->getValue("display/stroke_pts", b2);
    _dw->config()->io()->getValue("display/texture", b3);
    _dw->config()->io()->getValue("display/type", b4);
    _dw->config()->io()->getValue("display/antialiasing", b5);
    _dw->config()->io()->getValue("display/bbox", b6);
    _dw->config()->io()->getValue("gm/size/x", gmx);
    _dw->config()->io()->getValue("gm/size/y", gmy);
    _dw->config()->io()->getValue("border/factor/x", borderx);
    _dw->config()->io()->getValue("border/factor/y", bordery);
  }

  // Block signals in order not to trigger any callback.
  samplingSpinBox->blockSignals(true);
  samplingSpinBox->setValue(sampling);
  samplingSpinBox->blockSignals(false);

  snappingSpinBox->blockSignals(true);
  snappingSpinBox->setValue(snapping);
  snappingSpinBox->blockSignals(false);

  chainCheckBox->blockSignals(true);
  chainCheckBox->setChecked(b1);
  chainCheckBox->blockSignals(false);

  strokeCheckBox->blockSignals(true);
  strokeCheckBox->setChecked(b2);
  strokeCheckBox->blockSignals(false);

  textureCheckBox->blockSignals(true);
  textureCheckBox->setChecked(b3);
  textureCheckBox->blockSignals(false);

  typeCheckBox->blockSignals(true);
  typeCheckBox->setChecked(b4);
  typeCheckBox->blockSignals(false);

  antialiasingCheckBox->blockSignals(true);
  antialiasingCheckBox->setChecked(b5);
  antialiasingCheckBox->blockSignals(false);

  bbCheckBox->blockSignals(true);
  bbCheckBox->setChecked(b6);
  bbCheckBox->blockSignals(false);

  splittingSpinBox->blockSignals(true);
  splittingSpinBox->setValue(threshold);
  splittingSpinBox->blockSignals(false);

  gmXSpinBox->blockSignals(true);
  gmXSpinBox->setValue(gmx);
  gmXSpinBox->blockSignals(false);

  gmYSpinBox->blockSignals(true);
  gmYSpinBox->setValue(gmy);
  gmYSpinBox->blockSignals(false);

  borderXSpinBox->blockSignals(true);
  borderXSpinBox->setValue(borderx);
  borderXSpinBox->blockSignals(false);

  borderYSpinBox->blockSignals(true);
  borderYSpinBox->setValue(bordery);
  borderYSpinBox->blockSignals(false);
}

void DockCanvasWidget::samplingSpinBox_valueChanged(int val)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("intervals/sampling", val);
  _dw->configUpdated();
}

void DockCanvasWidget::snappingSpinBox_valueChanged(int val)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("intervals/snapping", val);
  _dw->configUpdated();
}

void DockCanvasWidget::splittingSpinBox_valueChanged(float val)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("splitting/threshold", val);
  _dw->configUpdated();
}

void DockCanvasWidget::chainCheckBox_toggled(bool b)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/chain_pts", b);
  _dw->configUpdated();
}

void DockCanvasWidget::strokeCheckBox_toggled(bool b)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/stroke_pts", b);
  _dw->configUpdated();
}

void DockCanvasWidget::textureCheckBox_toggled(bool b)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/texture", b);
  _dw->configUpdated();
}

void DockCanvasWidget::typeCheckBox_toggled(bool b)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/type", b);
  _dw->configUpdated();
}

void DockCanvasWidget::bbCheckBox_toggled(bool b)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/bbox", b);
  _dw->configUpdated();
}

void DockCanvasWidget::antialiasingCheckBox_toggled(bool b)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("display/antialiasing", b);
  _dw->configUpdated();
}

void DockCanvasWidget::gmXSpinBox_valueChanged(int val)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("gm/size/x", val);
  _dw->configUpdated();
}

void DockCanvasWidget::gmYSpinBox_valueChanged(int val)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("gm/size/y", val);
  _dw->configUpdated();
}


void DockCanvasWidget::borderXSpinBox_valueChanged(float val)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("border/factor/x", val);
  _dw->configUpdated();
}

void DockCanvasWidget::borderYSpinBox_valueChanged(float val)
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("border/factor/y", val);
  _dw->configUpdated();
}

void DockCanvasWidget::computeButton_clicked()
{
  if (!_dw || !_dw->config())
    return;
  _dw->config()->io()->setValue("compute/start", true);
  _dw->configUpdated();
}
