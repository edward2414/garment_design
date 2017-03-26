//
//  Filename         : main_window.cpp
//  Author           : Emmanuel Turquin
//  Purpose          : The application's main window, containing menubar,
//                     workspace, docks and msgbar.
//  Date of creation : 04/18/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include <qworkspace.h>
#include <qtoolbar.h>
#include <qpopupmenu.h>
#include <qmenubar.h>
#include <qstatusbar.h>
#include <qmessagebox.h>
#include <qapplication.h>
#include "mdi_window.h"
#include "dock_window.h"
#include "main_window.h"
#include "config_io.h"

MainWindow::MainWindow(ConfigIO *config, const QString& caption)
  : QMainWindow(0, "MainWindow", WDestructiveClose)
{
  setCaption(caption);

  QPopupMenu * file = new QPopupMenu(this);
  menuBar()->insertItem("&File", file);

  file->insertItem("&New", this, SLOT(reset()), CTRL+Key_N);
  file->insertItem("&Open...", this, SLOT(load()), CTRL+Key_O);
  file->insertItem("&Save", this, SLOT(save()), CTRL+Key_S);
  file->insertItem("Save &As...", this, SLOT(saveAs()));
  file->insertSeparator();
  file->insertItem("&Close", this, SLOT(closeWindow()), CTRL+Key_W);
  file->insertItem("&Quit", this, SLOT(close()), CTRL+Key_Q );

  _windows_menu = new QPopupMenu(this);
  _windows_menu->setCheckable(true);
  connect(_windows_menu, SIGNAL(aboutToShow()), this, SLOT(windowsMenuAboutToShow()));
  menuBar()->insertItem("&Windows", _windows_menu);

  // Initialize the windows menu for shortcuts (overwritten later).
  _windows_menu->insertItem("&Previous", this, SLOT(previousWindow()), Qt::Key_Prior);
  _windows_menu->insertItem("&Next", this, SLOT(nextWindow()), Qt::Key_Next);
  
  menuBar()->insertSeparator();

  QPopupMenu * help = new QPopupMenu(this);
  menuBar()->insertItem("&Help", help);
  
  help->insertItem("&Contents", this, SLOT(contents()), CTRL+Key_H);
  help->insertItem("&About", this, SLOT(about()));
  
  _ws = new QWorkspace(this);
  _ws->setScrollBarsEnabled(true);
  setCentralWidget(_ws);

  _config = config;
  if (_config)
    initConfig();
  
  statusBarMessage("Ready...", 2000);
}

void MainWindow::initConfig()
{
  int w = width();
  int h = height();
  _config->getValue("main_window/size/width", w);
  _config->getValue("main_window/size/height", h);
  
  int posx = x();
  int posy = y();
  _config->getValue("main_window/pos/x", posx);
  _config->getValue("main_window/pos/y", posy);
  
#ifdef WIN32
  resize(w, h);
  move(posx, posy);
#else
  setGeometry(posx, posy, w, h);
#endif // WIN32
}

MainWindow::~MainWindow()
{
  delete _ws;
  delete _windows_menu;
}


void MainWindow::reset()
{
  // FIXME
  emit statusBarMessage("Warning: reset() method not implemented.", 2000);
}

void MainWindow::load()
{
  MDIWindow* m = (MDIWindow*)_ws->activeWindow();
  if ( m )
    m->load();
}

void MainWindow::save()
{
  MDIWindow* m = (MDIWindow*)_ws->activeWindow();
  if ( m )
    m->save();
}


void MainWindow::saveAs()
{
  MDIWindow* m = (MDIWindow*)_ws->activeWindow();
  if ( m )
    m->saveAs();
}


void MainWindow::closeWindow()
{
  MDIWindow* m = (MDIWindow*)_ws->activeWindow();
  if ( m )
    m->close();
}


void MainWindow::contents()
{
  // FIXME:
  QMessageBox::warning(this, "Warning!!!",
		       "FIXME: no help yet");
}

void MainWindow::about()
{
  QMessageBox::about(this, "About " + config::APP_NAME + " " + config::APP_VERSION,
		     config::APP_ABOUT);
}

void MainWindow::previousWindow()
{
  _ws->activatePrevWindow();
  statusBarMessage("Focus to " + _ws->activeWindow()->caption() + ".", 2000);
}

void MainWindow::nextWindow()
{
  _ws->activateNextWindow();
  statusBarMessage("Focus to " + _ws->activeWindow()->caption() + ".", 2000);
}

void MainWindow::windowsMenuAboutToShow()
{
  _windows_menu->clear();
  int cascadeId = _windows_menu->insertItem("&Cascade", this, SLOT(cascade() ) );
  int tileId = _windows_menu->insertItem("&Tile", this, SLOT(tile() ) );
  int horTileId = _windows_menu->insertItem("Tile &Horizontally", this, SLOT(tileHorizontal() ) );
  int verTileId = _windows_menu->insertItem("Tile &Vertically", this, SLOT(tileVertical() ) );
  _windows_menu->insertSeparator();
  int prevId = _windows_menu->insertItem("&Previous", this, SLOT(previousWindow()), Qt::Key_Prior);
  int nextId = _windows_menu->insertItem("&Next", this, SLOT(nextWindow()), Qt::Key_Next);
  _windows_menu->insertSeparator();

  if (_ws->windowList(QWorkspace::CreationOrder).isEmpty()) {
    _windows_menu->setItemEnabled( cascadeId, FALSE );
    _windows_menu->setItemEnabled( tileId, FALSE );
    _windows_menu->setItemEnabled( horTileId, FALSE );
    _windows_menu->setItemEnabled( verTileId, FALSE );
    _windows_menu->setItemEnabled( prevId, FALSE );
    _windows_menu->setItemEnabled( nextId, FALSE );
  }

  QPopupMenu * windows_visible = new QPopupMenu(this);
  _windows_menu->insertItem("Visible", windows_visible);
  QPopupMenu * windows_hidden = new QPopupMenu(this);
  _windows_menu->insertItem("Hidden", windows_hidden);

  // MDIWindows
  QWidgetList mdi_windows = _ws->windowList(QWorkspace::CreationOrder);
  for ( int i = 0; i < int(mdi_windows.count()); ++i ) {
    if (mdi_windows.at(i)->isVisible()) {
      int id = windows_visible->insertItem(mdi_windows.at(i)->caption(),
					   this, SLOT( windowsMenuMDIActivated( int ) ) );
      windows_visible->setItemParameter( id, i );
      windows_visible->setItemChecked( id, _ws->activeWindow() == mdi_windows.at(i) );
    }
    else {
      int id = windows_hidden->insertItem(mdi_windows.at(i)->caption(),
					  this, SLOT( windowsMenuMDIActivated( int ) ) );
      windows_hidden->setItemParameter( id, i );
      windows_hidden->setItemChecked( id, _ws->activeWindow() == mdi_windows.at(i) );
    }
  }

  windows_visible->insertSeparator();
  windows_hidden->insertSeparator();

  // DockWindows
  QPtrList<QDockWindow> dock_windows = dockWindows();
  for (int i = 0; i < int(dock_windows.count()); ++i) {
    if (dock_windows.at(i)->isVisible()) {
      int id = windows_visible->insertItem(dock_windows.at(i)->caption(),
					   this, SLOT( windowsMenuDockActivated( int ) ) );
      windows_visible->setItemParameter( id, i );
    }
    else {
      int id = windows_hidden->insertItem(dock_windows.at(i)->caption(),
					  this, SLOT( windowsMenuDockActivated( int ) ) );
      windows_hidden->setItemParameter( id, i );
    }
  }
}

void MainWindow::windowsMenuMDIActivated( int id )
{
  QWidget* w = _ws->windowList(QWorkspace::CreationOrder).at( id );
  if (!w)
    return;
  w->showNormal();
  w->setFocus();
  statusBarMessage("Focus to " + w->caption() + ".", 2000);
}

void MainWindow::windowsMenuDockActivated( int id )
{
  QWidget* w = dockWindows().at( id );
  if (!w)
    return;
  w->showNormal();
  //  w->setFocus();
  statusBarMessage("Display " + w->caption() + ".", 2000);
}

void MainWindow::cascade()
{
  _ws->cascade();
}

void MainWindow::tile()
{
  _ws->tile();
}

void MainWindow::tileHorizontal()
{
  // primitive horizontal tiling
  QWidgetList windows = _ws->windowList(QWorkspace::CreationOrder);
  
  int visible_count = 0;
  int total_count = int(windows.count());
  for (int i = 0; i < total_count; ++i)
    if (windows.at(i)->isVisible())
      ++visible_count;
  
  int heightForEach = _ws->height() / visible_count;
  int y = 0;
  for (int i = 0; i < total_count; ++i) {
    QWidget *window = windows.at(i);
    if (window->isHidden())
      continue;
    if ( window->testWState( WState_Maximized ) ) {
      // prevent flicker
      window->hide();
      window->showNormal();
    }
    int preferredHeight = window->minimumHeight()+window->parentWidget()->baseSize().height();
    int actHeight = QMAX(heightForEach, preferredHeight);
    
    window->parentWidget()->setGeometry(0, y, _ws->width(), actHeight);
    y += actHeight;
  }
}

void MainWindow::tileVertical()
{
  // primitive vertical tiling
  QWidgetList windows = _ws->windowList(QWorkspace::CreationOrder);

  int visible_count = 0;
  int total_count = int(windows.count());
  for (int i = 0; i < total_count; ++i)
    if (windows.at(i)->isVisible())
      ++visible_count;
    
  int widthForEach = _ws->width() / visible_count;
  int x = 0;
  for ( int i = 0; i < total_count; ++i ) {
    QWidget *window = windows.at(i);
    if (window->isHidden())
      continue;
    if ( window->testWState( WState_Maximized ) ) {
      // prevent flicker
      window->hide();
      window->showNormal();
    }
    int preferredWidth = window->minimumWidth()+window->parentWidget()->baseSize().width();
    int actWidth = QMAX(widthForEach, preferredWidth);
	
    window->parentWidget()->setGeometry(x, 0, actWidth, _ws->height());
    x += actWidth;
  }
}

void MainWindow::statusBarMessage(const QString& msg, int t)
{
  statusBar()->message(msg, t);
}

QWorkspace *MainWindow::workspace() const
{
  return _ws;
}

void MainWindow::closeEvent( QCloseEvent *e )
{
  QWidgetList mdi_windows = _ws->windowList(QWorkspace::StackingOrder);
  for (int i = int(mdi_windows.count() - 1); i >= 0; --i) {
    MDIWindow* m = (MDIWindow*)mdi_windows.at(i);
    if (m)
      m->quit();
  }
  QPtrList<QDockWindow> dock_windows = dockWindows();
  for (int i = 0; i < int(dock_windows.count()); ++i) {
    DockWindow* d = (DockWindow*)dock_windows.at(i);
    if (d)
      d->quit();
  }

  if (_config) {
    _config->setValue("main_window/size/width", width());
    _config->setValue("main_window/size/height", height());
    _config->setValue("main_window/pos/x", x());
    _config->setValue("main_window/pos/y", y());
  }

  QMainWindow::closeEvent( e );
}
