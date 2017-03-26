//
//  Filename         : mdi_canvas.cpp
//  Author           : Emmanuel Turquin
//  Purpose          : An MDI window containing a 3DSViewer.
//  Date of creation : 04/20/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include <qpixmap.h>
#include "canvas.h"
#include "mdi_canvas.h"

#include "../data/pixmaps/canvas.xpm"

MDICanvas::MDICanvas(QWidget* parent, const char* name, int wflags, Config *conf)
  : MDIWindow(parent, name, wflags, conf)
{
  setIcon(QPixmap((const char**)canvas_xpm));
  _canvas = new Canvas(this, name);
  setFocusProxy(_canvas);
  setCentralWidget(_canvas);
}

MDICanvas::~MDICanvas()
{
  delete _canvas;
}

void MDICanvas::load()
{
  _canvas->loadFile();
}

void MDICanvas::updateConfig()
{
  _canvas->updateConfig();
  emit statusBarMessage(QString(name()) + " configuration updated.", 2000);
}
