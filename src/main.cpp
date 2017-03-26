//
//  Filename         : main.cpp
//  Author           : Emmanuel Turquin
//  Purpose          : File containing the main function.
//  Date of creation : 04/18/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include <qapplication.h>
#include <qworkspace.h>
#include <qgl.h>
#include "config_os_dependant.h"
#include "config_app_info.h"
#include "config_pretreatment.h"
#include "config_canvas.h"
#include "config_pattern.h"
#include "config_result.h"
#include "config_msg.h"
#include "config_io.h"
#include "config.h"
#include "main_window.h"

#include "mdi_pretreatment.h"
#include "dock_pretreatment.h"
#include "mdi_canvas.h"
#include "mdi_pattern.h"
#include "dock_canvas.h"
#include "mdi_result.h"
#include "dock_result.h"
#include "dock_pattern.h"
#include "dock_msg_handler.h"

// Only required to test adding items to ListView in dock_pattern_widget
#include "dock_pattern_widget.h"
#include <qlistview.h>


//DockPattern *dock_pattern;
//DockPattern *dock_pattern2;
int other_model_dl;

int main(int argc, char **argv )
{
  
  other_model_dl=0;
  // Get the home directory from the HOME environment variable.
  QString home_dir(getenv("HOME"));
  if (home_dir.isEmpty()) {
    std::cerr << "Warning: You may want to set the $HOME environment variable to use "
	      << config::APP_NAME << "."
	      << std::endl
	      << "Otherwise, configuration files will be stored in the current directory."
	      << std::endl;
    home_dir = ".";
  }
  QString config_dir(home_dir + config::DIR_SEP + "." + config::APP_NAME + config::DIR_SEP);

   
  // Create the application.
  QApplication app(argc, argv);

  // Create the configuration managers.
  ConfigIO interface_config(config_dir + config::APP_INTERFACE_FILE,
			    "interface_configuration",
			    true);
  Config config_canvas(config_dir + config::CANVAS_CONFIG_FILE, "canvas_configuration");
  Config config_pattern(config_dir + config::PATTERN_CONFIG_FILE, "pattern_configuration");
  Config config_pattern2(config_dir + config::PATTERN_CONFIG_FILE2, "pattern_configuration2");
  Config config_pretreatment(config_dir + config::PRETREATMENT_CONFIG_FILE, "pretreatment_configuration");
  Config config_result(config_dir + config::RESULT_CONFIG_FILE, "result_configuration");
  Config config_msg(config_dir + config::MSG_CONFIG_FILE, "msg_handler_configuration");

  // Create the main window.
  MainWindow *mw = new MainWindow(&interface_config);
  app.setMainWidget(mw);

  // Add MDI windows.
  MDIWindow *mdi_pretreatment = new MDIPretreatment((QWidget*)mw->workspace(), "Pretreatment", 0, &config_pretreatment);
  MDIWindow *mdi_canvas = new MDICanvas((QWidget*)mw->workspace(), "Canvas", 0, &config_canvas);
  
  DockPattern* dockpattern = new DockPattern(mw, "Pattern Options", 0, &config_pattern);
  DockPattern* dockpattern2 = new DockPattern(mw, "Pattern Options2", 0, &config_pattern2);
  
  MDIWindow *mdi_pattern = new MDIPattern((QWidget*)mw->workspace(), "Pattern", 0, dynamic_cast<QGLWidget*>(mdi_pretreatment->centralWidget()),&config_pattern,dockpattern);
  
  MDIWindow *mdi_pattern2 = new MDIPattern((QWidget*)mw->workspace(), "Pattern2", 0, dynamic_cast<QGLWidget*>(mdi_pretreatment->centralWidget()),&config_pattern2,dockpattern2);
  
  
  MDIWindow *mdi_result = new MDIResult((QWidget*)mw->workspace(), "Result", 0, dynamic_cast<QGLWidget*>(mdi_pretreatment->centralWidget()), &config_result);

  // Add dock windows.
  DockPretreatment *dock_pretreatment = new DockPretreatment(mw, "Pretreatment Options", 0, &config_pretreatment);
  DockCanvas *dock_canvas = new DockCanvas(mw, "Canvas Options", 0, &config_canvas);
  DockResult *dock_result = new DockResult(mw, "Result Options", 0, &config_result);
  //DockPattern *dock_pattern = new DockPattern(mw, "Pattern Options", 0, &config_pattern);
  
  //dock_pattern2 = new DockPattern(mw, "Pattern Options", 0, &config_pattern);
  DockMsgHandler *dock_msg_handler = new DockMsgHandler(mw, "MsgHandler", "", &config_msg);

  // test adding an item

  // Show main window.
  mw->show();

  // Update tabs of the MsgHandler.
  dock_msg_handler->updateTabs();

  // If it is the first start, tile MDI windows.
  bool first_start = true;
  interface_config.getValue("app/first_start", first_start);
  if (first_start) {
    interface_config.setValue("app/first_start", false);
    mw->tile();
  }

  // Set the focus to the first visible MDI window.
  QWidgetList mdi_windows = mw->workspace()->windowList(QWorkspace::CreationOrder);
  for (int i = 0; i < int(mdi_windows.count()); ++i) {
    if (mdi_windows.at(i)->isVisible()) {
      mdi_windows.at(i)->setFocus();
      break;
    }
  }

  // A little welcome message. :)
  dock_msg_handler->message(MsgHandler::MSG_NORMAL, config::APP_NAME,
			    "Welcome to <b>" + config::APP_NAME + " " + config::APP_VERSION +
			    "</b>! If you need help, press <b>Ctrl+H</b> to read the manual.");

  // Start the application's main loop.
  int res = app.exec();
  return res;
}
