//
//  Filename         : dock_pattern.h
//  Author           : Jamie Wither
//  Purpose          : The dock window for the Pattern module.
//  Date of creation : 13 Dec 2005
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  DOCK_PATTERN_H
# define DOCK_PATTERN_H

# include "dock_window.h"

class DockPatternWidget;

class DockPattern : public DockWindow
{
  Q_OBJECT

public:

  DockPattern(MainWindow* parent, const char* name = 0, WFlags f = 0, Config *conf = 0);
  virtual ~DockPattern();

  //  virtual void load();
  //  virtual void save();
  //  virtual void saveAs();
  //  void quit();

  virtual void updateConfig();

 private:

  DockPatternWidget	*_widget;
};

#endif // DOCK_PATTERN_H
