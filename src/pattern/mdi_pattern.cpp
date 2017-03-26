//
//  Filename         : mdi_pattern.cpp
//  Author           : Jamie Wither
//  Purpose          : An MDI window containing a Pattern widget.
//  Date of creation : 08 Dec 2005
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include <qpixmap.h>
#include "pattern.h"
#include "mdi_pattern.h"
#include "dock_pattern.h"

#include "../data/pixmaps/result.xpm"

MDIPattern::MDIPattern(QWidget* parent, const char* name, int wflags, QGLWidget *shared_qgl,  Config *conf, DockPattern* _dp)
  : MDIWindow(parent, name, wflags, conf)
{
  setIcon(QPixmap((const char**)result_xpm));
  _pattern = new Pattern(this, name, shared_qgl, 0, _dp);
  setFocusProxy(_pattern);
  setCentralWidget(_pattern);
}

MDIPattern::~MDIPattern()
{
  delete _pattern;
}


void MDIPattern::updateConfig()
{
  _pattern->updateConfig();
  emit statusBarMessage(QString(name()) + " configuration updated.", 2000);
}

void MDIPattern::load()
{
  
  _pattern->loadFile();
}
