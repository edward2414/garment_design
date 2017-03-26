//
//  Filename         : mdi_pretreatment.h
//  Author           : Emmanuel Turquin
//  Purpose          : An MDI window containing a Pretreatment widget.
//  Date of creation : 04/20/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  MDI_PRETREATMENT_H
# define MDI_PRETREATMENT_H

# include "mdi_window.h"

class Pretreatment;

class MDIPretreatment : public MDIWindow
{
public:

  MDIPretreatment(QWidget* parent, const char* name = 0, int wflags = 0, Config *conf = 0);
  ~MDIPretreatment();

  void load();
  //  void save();
  //  void saveAs();

  void updateConfig();

private:

  Pretreatment *_pretreatment;
};

#endif // MDI_PRETREATMENT_H
