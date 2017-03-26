//
//  Filename         : mdi_window.h
//  Author           : Emmanuel Turquin
//  Purpose          : The sub-window base class.
//  Date of creation : 04/18/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  MDI_WINDOW_H
# define MDI_WINDOW_H

# include <qmainwindow.h>
# include "msg_handler.h"

class Config;

class MDIWindow : public QMainWindow
{
  Q_OBJECT

public:

  MDIWindow(QWidget* parent, const char* name = 0, int wflags = 0, Config *conf = 0);
  virtual ~MDIWindow();

  virtual void load();
  virtual void save();
  virtual void saveAs();
  virtual void quit();

  virtual void updateConfig();
  void configUpdated();

  void setConfig(Config *conf);
  Config *config();

signals:

  void statusBarMessage(const QString&, int);
  void message(const MsgHandler::MsgType&, const QString&, const QString&);

protected:

  void closeEvent(QCloseEvent *);

  Config	*_conf;
};

#endif // MDI_WINDOW_H
