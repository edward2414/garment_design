//
//  Filename         : mdi_canvas.h
//  Author           : Emmanuel Turquin
//  Purpose          : An MDI window containing a 3DSViewer.
//  Date of creation : 04/20/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  MDI_CANVAS_H
# define MDI_CANVAS_H

# include "mdi_window.h"

class Canvas;

class MDICanvas : public MDIWindow
{
public:

  MDICanvas(QWidget* parent, const char* name = 0, int wflags = 0, Config *conf = 0);
  ~MDICanvas();

  void load();
  //  void save();
  //  void saveAs();

  void updateConfig();

private:

  Canvas *_canvas;
};

#endif // MDI_CANVAS_H
