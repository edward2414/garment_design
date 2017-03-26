//
//  Filename         : pattern.h
//  Author           : Jamie Wither
//  Purpose          : A viewer for the supplied developable pattern
//  Date of creation : 08 Dec 2005
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  PATTERN_H
# define PATTERN_H

# include "base_qglviewer.h"
#include <iostream>
#include <fstream>
#include <qapplication.h>
# include "vectypes.h"
#include "garment_parameters.h"


#include <string>

class DockPattern;
class MDIWindow;

class LCModel;
class Garment;

class PushBackParams;


class Pattern : public BaseQGLViewer
{
public:

  Pattern(MDIWindow* parent=NULL, const char* name=0, const QGLWidget* shareWidget=0, int wflags=0, DockPattern* dp=0);
  virtual ~Pattern();

  void updateConfig();

  void loadFile();
  
  
  std::string currentAtlasName;
  std::string currentSectionName;
  
  static void setClothingModelDisplayList(int);
  static Vec3d movePointOutsideModel(const Vec3d& point, void *clientData);
  bool _garmentInitialised;

  protected:
    virtual void keyPressEvent(QKeyEvent *e);
    virtual void mousePressEvent(QMouseEvent *e);
    virtual void mouseReleaseEvent(QMouseEvent *e);
    virtual void mouseMoveEvent(QMouseEvent *e);
  

private:

  virtual void draw();
  virtual void init();
  int  prepareModel();
  void freeModel();
  void renderModel();
  void convertCoordinates(double x, double y, double& new_x, double& new_y, double& ratio);
  void renderPatternFaces();

  int			_model_dl;
  static int           _clothes_horse_dl;
  static bool _clothes_horse_available;
  LCModel       *_lcmodel;
  Vec3d       _bb_size;

  Garment* _garment;
  


  bool          _axis_select_mode;
  bool			_disp_pattern;
  bool          _disp_model;
  bool          _disp_control_mesh;
  bool          _disp_model_mesh;
  bool          _disp_buckling_surface;
  bool          _twist_enabled;
  bool          _collision_enabled;
  float        _twist_factor;
  float        _twist_degree;
  float        _ptf_val;
  float       _hcompress_bottom;
  float       _hcompress_top;
  

  MDIWindow		*_mdi;
  DockPattern* _dock_pattern;
  
};

#endif // PATTERN_H

