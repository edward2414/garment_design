//
//  Filename         : config.cpp
//  Author(s)        : Emmanuel Turquin
//  Purpose          : Configuration management class.
//  Date of creation : 26/02/2003
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include "mdi_window.h"
#include "dock_window.h"
#include "config.h"

Config::Config(QString filename,
	       const QString& doc_type,
	       bool automatic,
	       MDIWindow *mdi,
	       DockWindow *dock) : _io(filename, doc_type, automatic) {
  _dock = dock;
  _mdi = mdi;
}

Config::~Config()
{

}

void Config::mdiUpdated() {
  if (_dock)
    _dock->updateConfig();
}

void Config::dockUpdated() {
  if (_mdi)
    _mdi->updateConfig();
}
