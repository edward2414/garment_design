//
//  Filename         : dock_window.h
//  Author           : Emmanuel Turquin
//  Purpose          : The dock-window base class.
//  Date of creation : 04/18/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  DOCK_WINDOW_H
# define DOCK_WINDOW_H

# include <qdockwindow.h>
# include "msg_handler.h"

class MainWindow;
class Config;

class DockWindow : public QDockWindow
{
  Q_OBJECT

public:

  typedef enum {
    DOCK_TOP,
    DOCK_BOTTOM,
    DOCK_LEFT,
    DOCK_RIGHT,
    UNDEFINED
  } DockLocation;

  DockWindow(MainWindow* parent, const char* name = 0, WFlags f = 0, Config *conf = 0);
  virtual ~DockWindow();

  virtual void load();
  virtual void save();
  virtual void saveAs();
  void quit();

  virtual void updateConfig();
  void configUpdated();

  void setConfig(Config *conf);
  Config *config();

protected slots:

  void setVisibility(bool b);
  void setPlace(QDockWindow::Place p);

signals:

  void statusBarMessage(const QString&, int);
  void message(const MsgHandler::MsgType&, const QString&, const QString&);

protected:

  unsigned	_place;
  Config	*_conf;
};

#endif // DOCK_WINDOW_H
