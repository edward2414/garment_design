//
//  Filename         : mdi_pretreatment.cpp
//  Author           : Emmanuel Turquin
//  Purpose          : An MDI window containing a Pretreatment widget.
//  Date of creation : 04/20/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include <qpixmap.h>
#include "pretreatment.h"
#include "mdi_pretreatment.h"

#include "../data/pixmaps/pretreatment.xpm"

MDIPretreatment::MDIPretreatment(QWidget* parent, const char* name, int wflags, Config *conf)
  : MDIWindow(parent, name, wflags, conf)
{
  setIcon(QPixmap((const char**)pretreatment_xpm));
  _pretreatment = new Pretreatment(this, name);
  setFocusProxy(_pretreatment);
  setCentralWidget(_pretreatment);
}

MDIPretreatment::~MDIPretreatment()
{
  delete _pretreatment;
}

void MDIPretreatment::load()
{
  _pretreatment->loadFile();
}

void MDIPretreatment::updateConfig()
{
  _pretreatment->updateConfig();
  emit statusBarMessage(QString(name()) + " configuration updated.", 2000);
}
