//
//  Filename         : dock_result.h
//  Author           : Emmanuel Turquin
//  Purpose          : The dock window for the Result module.
//  Date of creation : 04/18/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  DOCK_RESULT_H
# define DOCK_RESULT_H

# include "dock_window.h"

class DockResultWidget;

class DockResult : public DockWindow
{
  Q_OBJECT

public:

  DockResult(MainWindow* parent, const char* name = 0, WFlags f = 0, Config *conf = 0);
  virtual ~DockResult();

  //  virtual void load();
  //  virtual void save();
  //  virtual void saveAs();
  //  void quit();

  virtual void updateConfig();

 private:

  DockResultWidget	*_widget;
};

#endif // DOCK_RESULT_H
