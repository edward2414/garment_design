//
//  Filename         : result.h
//  Author           : Emmanuel Turquin
//  Purpose          : A viewer for the generated garments.
//  Date of creation : 04/20/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  RESULT_H
# define RESULT_H

# include "base_qglviewer.h"
#include <iostream>
#include <fstream>
#include <qcolor.h>

class MDIWindow;

class Result : public BaseQGLViewer
{
public:

  Result(MDIWindow* parent=NULL, const char* name=0,
	       const QGLWidget* shareWidget=0, int wflags=0);
  virtual ~Result();

  void updateConfig();

  void pretreatmentUpdated();
  void canvasUpdated();

  void saveFile();

private:

  virtual void draw();
  virtual void init();

  void renderModel();

  void renderSurface();
  void freeSurface();
  int prepareFrontSurface();
  int prepareBackSurface();

  void renderContour();
  void freeContour();
  int prepareFrontContour();
  int prepareFrontSeamContour();
  int prepareBackSeamContour();
  int prepareBackContour();

  int			_model_dl;
  int			_front_surface_dl;
  int			_back_surface_dl;
  int			_front_contour_dl;
  int         _front_seam_contour_dl;
  int         _back_seam_contour_dl;
  int			_back_contour_dl;

  bool			_disp_model;
  bool			_disp_surface;
  bool			_disp_contour;
  bool			_disp_type;
  bool          _disp_surface_wireframe;
  QColor        _renderColour;
  

  MDIWindow		*_mdi;
};

#endif // RESULT_H
