//
//  Filename         : config.h
//  Author(s)        : Emmanuel Turquin
//  Purpose          : Configuration management class.
//  Date of creation : 26/02/2003
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  CONFIG_H
# define CONFIG_H

# include <qstring.h>
# include "config_io.h"

class DockWindow;
class MDIWindow;

class Config
{
public:

  Config(QString filename = "",
	 const QString& doc_type = "config_file",
	 bool automatic = true,
	 MDIWindow *mdi = 0,
	 DockWindow *dock = 0);
  ~Config();

  void setMDIWindow(MDIWindow *mdi) {
    _mdi = mdi;
  }

  void setDockWindow(DockWindow *dock) {
    _dock = dock;
  }

  MDIWindow *mdiWindow() {
    return _mdi;
  }

  DockWindow *dockWindow() {
    return _dock;
  }

  void mdiUpdated();

  void dockUpdated();

  ConfigIO *io() {
    return &_io;
  }

private:

  ConfigIO	_io;
  DockWindow	*_dock;
  MDIWindow	*_mdi;
};

#endif // CONFIG_H
