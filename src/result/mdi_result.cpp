//
//  Filename         : mdi_result.cpp
//  Author           : Emmanuel Turquin
//  Purpose          : An MDI window containing a Result widget.
//  Date of creation : 04/20/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include <qpixmap.h>
#include "result.h"
#include "mdi_result.h"

#include "../data/pixmaps/result.xpm"

MDIResult::MDIResult(QWidget* parent, const char* name, int wflags, QGLWidget *shared_qgl, Config *conf)
  : MDIWindow(parent, name, wflags, conf)
{
  setIcon(QPixmap((const char**)result_xpm));
  _result = new Result(this, name, shared_qgl);
  setFocusProxy(_result);
  setCentralWidget(_result);
}

MDIResult::~MDIResult()
{
  delete _result;
}

void MDIResult::save()
{
  _result->saveFile();
}

void MDIResult::saveAs()
{
  _result->saveFile();
}

void MDIResult::updateConfig()
{
  _result->updateConfig();
  emit statusBarMessage(QString(name()) + " configuration updated.", 2000);
}
