//
//  Filename         : dock_canvas.h
//  Author           : Emmanuel Turquin
//  Purpose          : The dock window for the Canvas module.
//  Date of creation : 04/18/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  DOCK_CANVAS_H
# define DOCK_CANVAS_H

# include "dock_window.h"

class DockCanvasWidget;

class DockCanvas : public DockWindow
{
  Q_OBJECT

public:

  DockCanvas(MainWindow* parent, const char* name = 0, WFlags f = 0, Config *conf = 0);
  virtual ~DockCanvas();

  //  virtual void load();
  //  virtual void save();
  //  virtual void saveAs();
  //  void quit();

  virtual void updateConfig();

 private:

  DockCanvasWidget	*_widget;
};

#endif // DOCK_CANVAS_H
