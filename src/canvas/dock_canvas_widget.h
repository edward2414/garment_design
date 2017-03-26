//
//  Filename         : dock_canvas_widget.h
//  Author           : Emmanuel Turquin
//  Purpose          : The main widget of the dock for the Pretreatment module.
//  Date of creation : 04/18/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  DOCK_CANVAS_WIDGET_H
# define DOCK_CANVAS_WIDGET_H

# include "dock_canvas_widget_base.h"

class DockWindow;

class DockCanvasWidget : public DockCanvasWidgetBase
{
  Q_OBJECT

public:

  DockCanvasWidget(DockWindow* parent = 0, const char* name = 0, WFlags fl = 0);
  ~DockCanvasWidget();

  void updateConfig();

public slots:

  virtual void samplingSpinBox_valueChanged(int val);
  virtual void snappingSpinBox_valueChanged(int val);
  virtual void splittingSpinBox_valueChanged(float val);
  virtual void chainCheckBox_toggled(bool);
  virtual void strokeCheckBox_toggled(bool);
  virtual void textureCheckBox_toggled(bool);
  virtual void typeCheckBox_toggled(bool);
  virtual void bbCheckBox_toggled(bool);
  virtual void antialiasingCheckBox_toggled(bool);
  virtual void gmXSpinBox_valueChanged(int val);
  virtual void gmYSpinBox_valueChanged(int val);
  virtual void borderXSpinBox_valueChanged(float val);
  virtual void borderYSpinBox_valueChanged(float val);
  virtual void computeButton_clicked();

private:

  DockWindow *_dw;
};

#endif // DOCK_CANVAS_WIDGET_H
