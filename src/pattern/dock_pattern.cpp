//
//  Filename         : dock_pattern.cpp
//  Author           : Jamie Wither
//  Purpose          : The dock window for the Pattern module.
//  Date of creation : 13 Dec 2005
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include "main_window.h"
#include "config.h"
#include "dock_pattern_widget.h"
#include "dock_pattern.h"

DockPattern::DockPattern(MainWindow* parent, const char* name, WFlags f, Config *conf) :
  DockWindow(parent, name, f, conf)
{
  _widget = new DockPatternWidget(this);
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

DockPattern::~DockPattern()
{
  delete _widget;
}

void DockPattern::updateConfig()
{
  _widget->updateConfig();
}
