//
//  Filename         : dock_result.cpp
//  Author           : Emmanuel Turquin
//  Purpose          : The dock window for the Result module.
//  Date of creation : 04/18/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include "main_window.h"
#include "config.h"
#include "dock_result_widget.h"
#include "dock_result.h"

DockResult::DockResult(MainWindow* parent, const char* name, WFlags f, Config *conf) :
  DockWindow(parent, name, f, conf)
{
  _widget = new DockResultWidget(this);
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

DockResult::~DockResult()
{
  delete _widget;
}

void DockResult::updateConfig()
{
  _widget->updateConfig();
}
