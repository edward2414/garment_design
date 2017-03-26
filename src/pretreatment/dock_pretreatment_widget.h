//
//  Filename         : dock_pretreatment_widget.h
//  Author           : Emmanuel Turquin
//  Purpose          : The main widget of the dock for the Pretreatment module.
//  Date of creation : 04/18/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  DOCK_PRETREATMENT_WIDGET_H
# define DOCK_PRETREATMENT_WIDGET_H

# include "dock_pretreatment_widget_base.h"

class DockWindow;

class DockPretreatmentWidget : public DockPretreatmentWidgetBase
{
  Q_OBJECT

public:

  DockPretreatmentWidget(DockWindow* parent = 0, const char* name = 0, WFlags fl = 0);
  ~DockPretreatmentWidget();

 void updateConfig();

public slots:

  virtual void perspectiveButton_toggled(bool);
  virtual void orthoButton_toggled(bool);
  virtual void bbCheckBox_toggled(bool);
  virtual void modelCheckBox_toggled(bool);
  virtual void bbZSpinBox_valueChanged(float);
  virtual void bbYSpinBox_valueChanged(float);
  virtual void bbXSpinBox_valueChanged(float);
  virtual void octreeCheckBox_toggled(bool);
  virtual void octreeLevelSpinBox_valueChanged(int);
  virtual void bbAxisCheckBox_toggled(bool);
  virtual void computeButton_clicked();
  virtual void bbResetButton_clicked();
  virtual void textureSizeSpinBox_valueChanged(int);
  virtual void dfXSpinBox_valueChanged(int);
  virtual void dfYSpinBox_valueChanged(int);
  virtual void dfZSpinBox_valueChanged(int);
  virtual void dfCheckBox_toggled(bool);
  virtual void dfXMinSpinBox_valueChanged(int);
  virtual void dfXMaxSpinBox_valueChanged(int);
  virtual void dfYMinSpinBox_valueChanged(int);
  virtual void dfYMaxSpinBox_valueChanged(int);
  virtual void dfZMinSpinBox_valueChanged(int);
  virtual void dfZMaxSpinBox_valueChanged(int);
    virtual void dfSliceSlider_valueChanged(int);
    virtual void dfSliceSpinBox_valueChanged(int);
    
private:

  DockWindow	*_dw;
};

#endif // DOCK_PRETREATMENT_WIDGET_H
