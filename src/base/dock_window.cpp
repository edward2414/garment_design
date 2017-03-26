//
//  Filename         : dock_window.cpp
//  Author           : Emmanuel Turquin
//  Purpose          : The dock-window base class.
//  Date of creation : 04/18/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include "config.h"
#include "config_dock.h"
#include "main_window.h"
#include "dock_window.h"

DockWindow::DockWindow(MainWindow* parent, const char* name, WFlags f, Config *conf)
  : QDockWindow((QWidget*)parent, name, f)
{
  setCaption(name);
  if (caption().isEmpty())
    setCaption(config::DOCK_DEFAULT_CAPTION);
  setResizeEnabled(true);
  setOpaqueMoving(true);
  setHorizontallyStretchable(true);
  setVerticallyStretchable(true);
  setCloseMode(QDockWindow::Always);
  setFixedExtentWidth(config::DOCK_MIN_WIDTH);
  setFixedExtentHeight(config::DOCK_MIN_HEIGHT);
  connect(this, SIGNAL(statusBarMessage(const QString&, int)),
	  topLevelWidget(), SLOT(statusBarMessage(const QString&, int)));
  connect(this, SIGNAL(message(const MsgHandler::MsgType&, const QString&, const QString&)),
	  topLevelWidget(), SIGNAL(message(const MsgHandler::MsgType&, const QString&, const QString&)));
  connect(this, SIGNAL(visibilityChanged(bool)),
	  this, SLOT(setVisibility(bool)));
  _place = UNDEFINED;
  _conf = conf;
  if (_conf) {
    connect(this, SIGNAL(placeChanged(QDockWindow::Place)),
	    this, SLOT(setPlace(QDockWindow::Place)));

    _conf->setDockWindow(this);

    _place = UNDEFINED;
    _conf->io()->getValue("dock_window/place", _place);
    switch (_place) {
    case DOCK_TOP:
      parent->moveDockWindow(this, Qt::DockTop);
      break;
    case DOCK_BOTTOM:
      parent->moveDockWindow(this, Qt::DockBottom);
      break;
    case DOCK_LEFT:
      parent->moveDockWindow(this, Qt::DockLeft);
      break;
    case DOCK_RIGHT:
      parent->moveDockWindow(this, Qt::DockRight);
      break;
    case UNDEFINED:
    default:
      ;
    }

    bool indock = true;
    _conf->io()->getValue("dock_window/indock", indock);
    if (!indock)
      undock();

    int w = width();
    int h = height();
    _conf->io()->getValue("dock_window/size/width", w);
    _conf->io()->getValue("dock_window/size/height", h);

    int posx = x();
    int posy = y();
    _conf->io()->getValue("dock_window/pos/x", posx);
    _conf->io()->getValue("dock_window/pos/y", posy);
 
#ifdef WIN32
    resize(w, h);
    move(posx, posy);
#else
    setGeometry(posx, posy, w, h);
#endif // WIN32

    bool vis = true;
    _conf->io()->getValue("dock_window/visible", vis);
    if (vis)
      show();
    else
      hide();
  }
}

DockWindow::~DockWindow()
{

}

void DockWindow::setVisibility(bool b)
{
  if (!b)
    emit statusBarMessage("Hide " + caption() + ".", 2000);
}

void DockWindow::load()
{
  emit statusBarMessage("Warning: load() method not implemented.", 2000);
}

void DockWindow::save()
{
  emit statusBarMessage("Warning: save() method not implemented.", 2000);
}

void DockWindow::saveAs()
{
  emit statusBarMessage("Warning: saveAs() method not implemented.", 2000);
}

void DockWindow::setPlace(QDockWindow::Place)
{
  if (_conf) {
    QDockArea *a = area();
    if (!a) {
      _conf->io()->setValue("dock_window/indock", false);
      return;
    }
    _conf->io()->setValue("dock_window/indock", true);
    MainWindow *mw = (MainWindow*)topLevelWidget();
    if (a == mw->topDock())
      _place = DOCK_TOP;
    else if (a == mw->bottomDock())
      _place = DOCK_BOTTOM;
    else if (a == mw->leftDock())
      _place = DOCK_LEFT;
    else if (a == mw->rightDock())
      _place = DOCK_RIGHT;
    else
      _place = UNDEFINED;
    _conf->io()->setValue("dock_window/place", _place);
  }
}

void DockWindow::quit()
{
  if (_conf) {
    _conf->io()->setValue("dock_window/size/width", width());
    _conf->io()->setValue("dock_window/size/height", height());
    _conf->io()->setValue("dock_window/pos/x", x());
    _conf->io()->setValue("dock_window/pos/y", y());
    _conf->io()->setValue("dock_window/visible", isVisible());
  }
}

void DockWindow::updateConfig()
{
  emit statusBarMessage("Warning: updateConfig() method not implemented.", 2000);
}

void DockWindow::configUpdated()
{
  if (_conf)
    _conf->dockUpdated();
}

void DockWindow::setConfig(Config *conf)
{
  _conf = conf;
}

Config *DockWindow::config()
{
  return _conf;
}
