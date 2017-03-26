//
//  Filename         : dock_pretreatment.cpp
//  Author           : Emmanuel Turquin
//  Purpose          : The dock window for the Pretreatment module.
//  Date of creation : 04/18/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include "main_window.h"
#include "config.h"
#include "dock_pretreatment_widget.h"
#include "dock_pretreatment.h"

DockPretreatment::DockPretreatment(MainWindow* parent, const char* name, WFlags f, Config *conf) :
  DockWindow(parent, name, f, conf)
{
  _widget = new DockPretreatmentWidget(this);
  setWidget(_widget);

  if (parent) {
    parent->setDockEnabled(this, Qt::DockTop, false);
    parent->setDockEnabled(this, Qt::DockBottom, false);
    if (_place == UNDEFINED) {
      parent->moveDockWindow(this, Qt::DockLeft);
      if (_conf)
	_conf->io()->setValue("dock_window/place", DOCK_LEFT);
    }
  }
}

DockPretreatment::~DockPretreatment()
{
  delete _widget;
}

void DockPretreatment::updateConfig()
{
  _widget->updateConfig();
}
