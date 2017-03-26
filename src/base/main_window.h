//
//  Filename         : main_window.h
//  Author           : Emmanuel Turquin
//  Purpose          : The application's main window, containing menubar,
//                     workspace, docks and msgbar.
//  Date of creation : 04/18/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  MAIN_WINDOW_H
# define MAIN_WINDOW_H

# include <qmainwindow.h>
# include "config_app_info.h"
# include "msg_handler.h"

class QToolBar;
class QPopupMenu;
class QWorkspace;
class ConfigIO;

class MainWindow: public QMainWindow
{
  Q_OBJECT

public:

  MainWindow(ConfigIO *config, const QString& caption = config::APP_NAME + " " + config::APP_VERSION);
  ~MainWindow();

  void initConfig();

  QWorkspace	*workspace() const;

protected:
  
  void closeEvent(QCloseEvent *);

public slots:

  void cascade();
  void tile();
  void tileHorizontal();
  void tileVertical();

  void statusBarMessage(const QString&, int);

private slots:
  
  void reset();
  void load();
  void save();
  void saveAs();
  void closeWindow();

  void contents();
  void about();

  void previousWindow();
  void nextWindow();

  void windowsMenuAboutToShow();
  void windowsMenuMDIActivated(int id);
  void windowsMenuDockActivated(int id);

signals:

  void message(const MsgHandler::MsgType&, const QString&, const QString&);

private:
  
  QWorkspace	*_ws;
  QPopupMenu	*_windows_menu;
  ConfigIO	*_config;
};

#endif // MAIN_WINDOW_H
