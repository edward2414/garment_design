//
//  Filename         : mdi_pattern.h
//  Author           : Jamie Wither
//  Purpose          : An MDI window containing a Pattern widget.
//  Date of creation : 08 Dec 2005
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  MDI_PATTERN_H
# define MDI_PATTERN_H

# include "mdi_window.h"

class QGLWidget;
class Pattern;
class DockPattern;

class MDIPattern : public MDIWindow
{
public:

  MDIPattern(QWidget* parent, const char* name = 0, int wflags = 0, QGLWidget *shared_qgl = 0, Config *conf = 0, DockPattern* _dp=0);
  ~MDIPattern();

  void load();
  //void save();
  //void saveAs();

  void updateConfig();

private:

  Pattern *_pattern;
};

#endif // MDI_PATTERN_H
