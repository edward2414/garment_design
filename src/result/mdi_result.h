//
//  Filename         : mdi_result.h
//  Author           : Emmanuel Turquin
//  Purpose          : An MDI window containing a Result widget.
//  Date of creation : 04/20/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  MDI_RESULT_H
# define MDI_RESULT_H

# include "mdi_window.h"

class QGLWidget;
class Result;

class MDIResult : public MDIWindow
{
public:

  MDIResult(QWidget* parent, const char* name = 0, int wflags = 0, QGLWidget *shared_qgl = 0, Config *conf = 0);
  ~MDIResult();

  //  void load();
  void save();
  void saveAs();

  void updateConfig();

private:

  Result *_result;
};

#endif // MDI_RESULT_H
