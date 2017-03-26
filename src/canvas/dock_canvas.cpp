//
//  Filename         : dock_canvas.h
//  Author           : Emmanuel Turquin
//  Purpose          : The dock window for the Canvas module.
//  Date of creation : 04/18/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include "main_window.h"
#include "dock_canvas_widget.h"
#include "dock_canvas.h"

DockCanvas::DockCanvas(MainWindow* parent, const char* name, WFlags f, Config *conf) :
  DockWindow(parent, name, f, conf)
{
  _widget = new DockCanvasWidget(this);
  setWidget(_widget);

  if (parent) {
    parent->setDockEnabled(this, Qt::DockTop, false);
    parent->setDockEnabled(this, Qt::DockBottom, false);
    if (_place == UNDEFINED) {
      parent->moveDockWindow(this, Qt::DockRight);
      if (_conf)
	_conf->io()->setValue("dock_window/place", DOCK_RIGHT);
    }
  }
}

DockCanvas::~DockCanvas()
{
  delete _widget;
}

void DockCanvas::updateConfig()
{
  _widget->updateConfig();
}
