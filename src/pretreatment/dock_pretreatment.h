//
//  Filename         : dock_pretreatment.h
//  Author           : Emmanuel Turquin
//  Purpose          : The dock window for the Pretreatment module.
//  Date of creation : 04/18/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  DOCK_PRETREATMENT_H
# define DOCK_PRETREATMENT_H

# include "dock_window.h"

class DockPretreatmentWidget;

class DockPretreatment : public DockWindow
{
  Q_OBJECT

public:

  DockPretreatment(MainWindow* parent, const char* name = 0, WFlags f = 0, Config *conf = 0);
  virtual ~DockPretreatment();

  //  virtual void load();
  //  virtual void save();
  //  virtual void saveAs();
  //  void quit();

  virtual void updateConfig();

 private:

  DockPretreatmentWidget	*_widget;
};

#endif // DOCK_PRETREATMENT_H
