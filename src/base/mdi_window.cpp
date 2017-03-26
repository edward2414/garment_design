//
//  Filename         : mdi_window.cpp
//  Author           : Emmanuel Turquin
//  Purpose          : The sub-window base class.
//  Date of creation : 04/18/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include <qfile.h>
#include <qfiledialog.h>
#include <qmessagebox.h>
#include <qmainwindow.h>
#include <qpixmap.h>
#include "config.h"
#include "config_mdi.h"
#include "config_app_info.h"
#include "mdi_window.h"

#include "../data/pixmaps/default.xpm"

MDIWindow::MDIWindow(QWidget* parent, const char* name, int wflags, Config *conf)
  : QMainWindow(parent, name, wflags)
{
  setCaption(name);
  setIcon(QPixmap((const char**)default_xpm));
  if (caption().isEmpty())
    setCaption(config::MDI_DEFAULT_CAPTION);
  setMinimumWidth(config::MDI_MIN_WIDTH);
  setMinimumHeight(config::MDI_MIN_HEIGHT);
  connect(this, SIGNAL(statusBarMessage(const QString&, int)),
	  topLevelWidget(), SLOT(statusBarMessage(const QString&, int)));
  connect(this, SIGNAL(message(const MsgHandler::MsgType&, const QString&, const QString&)),
	  topLevelWidget(), SIGNAL(message(const MsgHandler::MsgType&, const QString&, const QString&)));
  _conf = conf;

  // Load and initialize configuration.
  if (_conf) {
    _conf->setMDIWindow(this);

    int w = width();
    int h = height();
    _conf->io()->getValue("mdi_window/size/width", w);
    _conf->io()->getValue("mdi_window/size/height", h);

    int posx = x();
    int posy = y();
    _conf->io()->getValue("mdi_window/pos/x", posx);
    _conf->io()->getValue("mdi_window/pos/y", posy);

#ifdef WIN32
    resize(w, h);
    move(posx, posy);
#else
    setGeometry(posx, posy, w, h);
#endif // WIN32

    bool vis = true;
    _conf->io()->getValue("mdi_window/visible", vis);
    if (vis)
      show();
    else
      hide();

    bool max = false;
    _conf->io()->getValue("mdi_window/maximized", max);
    if (max)
      showMaximized();

    bool min = false;
    _conf->io()->getValue("mdi_window/minimized", min);
    if (min)
      showMinimized();
  }
}

MDIWindow::~MDIWindow()
{
}

void MDIWindow::closeEvent(QCloseEvent *)
{
  hide();
  emit statusBarMessage("Hide " + caption() + ".", 2000);
}

void MDIWindow::load()
{
  emit statusBarMessage("Warning: load() method not implemented.", 2000);
}

void MDIWindow::save()
{
  emit statusBarMessage("Warning: save() method not implemented.", 2000);
}

void MDIWindow::saveAs()
{
  emit statusBarMessage("Warning: saveAs() method not implemented.", 2000);
}

void MDIWindow::quit()
{
  if (_conf) {
    _conf->io()->setValue("mdi_window/maximized", isMaximized());
    _conf->io()->setValue("mdi_window/minimized", isMinimized());
    _conf->io()->setValue("mdi_window/visible", isVisible());
    if (isVisible())
      showNormal();
    _conf->io()->setValue("mdi_window/size/width", width());
    _conf->io()->setValue("mdi_window/size/height", height());
    _conf->io()->setValue("mdi_window/pos/x", ((QWidget*)parent())->x()); // Why parent() ???
    _conf->io()->setValue("mdi_window/pos/y", ((QWidget*)parent())->y()); // Why parent() ???
  }
}

void MDIWindow::updateConfig()
{
  emit statusBarMessage("Warning: updateConfig() method not implemented.", 2000);
}

void MDIWindow::configUpdated()
{
  if (_conf)
    _conf->mdiUpdated();
}

void MDIWindow::setConfig(Config *conf)
{
  _conf = conf;
}

Config *MDIWindow::config()
{
  return _conf;
}
