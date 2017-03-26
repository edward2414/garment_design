//
//  Filename         : dock_pattern_widget.h
//  Author           : Jamie
//  Purpose          : The main widget of the dock for the Pattern module.
//  Date of creation : 13 Dec 2005
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  DOCK_PATTERN_WIDGET_H
# define DOCK_PATTERN_WIDGET_H

# include "dock_pattern_widget_base.h"

class DockWindow;

class DockPatternWidget : public DockPatternWidgetBase
{
  Q_OBJECT

public:

  DockPatternWidget(DockWindow* parent = 0, const char* name = 0, WFlags fl = 0);
  ~DockPatternWidget();

  void updateConfig();


public slots:
    virtual void patternCheckBox_toggled(bool);
    virtual void switchAtlas(QListViewItem*);
    virtual void modelMeshCheckBox_toggled(bool);
    virtual void twistSpinBox_valueChanged(float);
    virtual void patchTangentFactorFloatSpinBox_valueChanged(float);
    virtual void twistFactorFloatSpinBox_valueChanged(float);
    virtual void hCompressTopFloatSpinBox_valueChanged(float);
    virtual void hCompressBottomFloatSpinBox_valueChanged(float);
    virtual void Redraw_clicked();
    virtual void controlMeshCheckBox_toggled(bool);
    virtual void atlasRowsSpinBox_valueChanged(int);
    virtual void atlasColsSpinBox_valueChanged(int);
    virtual void twistEnabledCheckBox_toggled(bool);
    virtual void bSurfaceCheckBox_toggled(bool);
    virtual void displayModelCheckBox_toggled(bool);
    virtual void collisionEnabledCheckBox_toggled(bool);
    
    virtual void axisAlignButton_clicked();
    virtual void rightEdgesSpinBox_valueChanged(int);
    virtual void leftEdgesSpinBox_valueChanged(int);
    virtual void bottomEdgesSpinBox_valueChanged(int);
    virtual void topEdgesSpinBox_valueChanged(int);
    
    virtual void slider2_valueChanged(int);
    virtual void slider1_valueChanged(int);
    virtual void effect4CheckBox_toggled(bool);
    virtual void effect3CheckBox_toggled(bool);
    virtual void effect2CheckBox_toggled(bool);
    virtual void effect1CheckBox_toggled(bool);

private:

  DockWindow *_dw;
};

#endif // DOCK_PATTERN_WIDGET_H
