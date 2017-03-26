//
//  Filename         : base_qglviewer.cpp
//  Author           : Emmanuel Turquin
//  Purpose          : An intermediate class between QGLViewer and
//                     the application's viewers.
//  Date of creation : 04/21/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include "base_qglviewer.h"

BaseQGLViewer::BaseQGLViewer(QWidget* parent, const char* name,
			     const QGLWidget* shareWidget, int wflags)
  : QGLViewer(parent, name, shareWidget, wflags)
{
  connect(this, SIGNAL(statusBarMessage(const QString&, int)),
	  topLevelWidget(), SLOT(statusBarMessage(const QString&, int)));
  connect(this, SIGNAL(message(const MsgHandler::MsgType&, const QString&, const QString&)),
	  topLevelWidget(), SIGNAL(message(const MsgHandler::MsgType&, const QString&, const QString&)));

  setShortcut(EXIT_VIEWER, 0);
  setShortcut(HELP, 0);
  //setShortcut(DISPLAY_Z_BUFFER, 0);
  setShortcut(DISPLAY_FPS, 0);
  setShortcut(DRAW_GRID, 0);
  setShortcut(DRAW_AXIS, 0);
  setShortcut(STEREO, 0);
  setShortcut(ANIMATION, 0);
  setShortcut(CAMERA_MODE, 0);
  setShortcut(EDIT_CAMERA, 0);
  setShortcut(SAVE_SCREENSHOT, Qt::ALT+Qt::Key_S);

  setStateFileName(QString::null);
};
