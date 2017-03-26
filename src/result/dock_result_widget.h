//
//  Filename         : dock_result_widget.h
//  Author           : Emmanuel Turquin
//  Purpose          : The main widget of the dock for the Pretreatment module.
//  Date of creation : 04/18/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  DOCK_RESULT_WIDGET_H
# define DOCK_RESULT_WIDGET_H

# include "dock_result_widget_base.h"

class DockWindow;

class DockResultWidget : public DockResultWidgetBase
{
  Q_OBJECT

public:

  DockResultWidget(DockWindow* parent = 0, const char* name = 0, WFlags fl = 0);
  ~DockResultWidget();

  void updateConfig();


public slots:
    virtual void modelCheckBox_toggled(bool);
    virtual void surfaceCheckBox_toggled(bool);
    virtual void surfaceRenderTypeCheckBox_toggled(bool);
    virtual void contourCheckBox_toggled(bool);
    virtual void typeCheckBox_toggled(bool);
    virtual void colorDialogButton_clicked();

private:

  DockWindow *_dw;
};

#endif // DOCK_RESULT_WIDGET_H
