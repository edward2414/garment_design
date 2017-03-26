//
//  Filename         : canvas.h
//  Author           : Emmanuel Turquin
//  Purpose          : A 2D drawing canvas.
//  Date of creation : 04/23/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  _CANVAS_H
# define _CANVAS_H

# include <list>
# include "base_qglviewer.h"
#include <iostream>
#include "vectypes.h"


class MDIWindow;
class Layer;
class GarmentMap;
class Point;


class Canvas : public BaseQGLViewer
{
public:

  typedef std::list<GarmentMap*>	GarmentMapList;
  // JDW Added so I could generate the seam line from the Result Class
  typedef std::list < std::list<Point> > ListOfPointLists;
  typedef enum { OFF, READY, RIGHT, LEFT } MirrorType;
  
  
  ListOfPointLists listOfFrontSeamPointLists;
  ListOfPointLists listOfBackSeamPointLists;

  Canvas(MDIWindow* parent=NULL, const char* name=0,
	 const QGLWidget* shareWidget=0, int wflags=0);
  ~Canvas();

  void canvasMessage(const MsgHandler::MsgType&, const QString&);

  void updateConfig();
  
  void loadFile();

  void pretreatmentUpdated();
  void resultUpdated();
  ListOfPointLists* getFrontSeams() {
    return &listOfFrontSeamPointLists;
  }
  ListOfPointLists* getBackSeams() {
    return &listOfBackSeamPointLists;
  }
  
  
  void clearSeams();
  void saveLayers(std::ostream& out);
  void loadLayers(std::istream& in);
  void readToLayer(std::istream& in, Layer* layer);
    

protected:

  virtual void draw();
  virtual void init();
  virtual void keyPressEvent(QKeyEvent *e);
  virtual void mouseMoveEvent(QMouseEvent *e);
  virtual void mousePressEvent(QMouseEvent *e);
  virtual void mouseReleaseEvent(QMouseEvent *e);

private:

  void initLayers();
  void freeLayers();
  //  void copyLayer();

  void initGarmentMaps();
  void freeGarmentMaps();

  void preparePaperSheet();
  void freePaperSheet();

  void prepareTextures();
  void freeTextures();

  void prepareBBox();
  void freeBBox();

  void computeSkeleton();
  void computeSurface();
  void computeMesh();

  void convertCoordinates(double x, double y, double& new_x, double& new_y, double& ratio);
  
  // Draw at pos, dir is direction from the base to the top of the drawn gaussian profile
  void drawFoldEndParameters(const Vec3d& pos, const Vec3d& dir, const double& radius, const double& amplitude);
  

  inline double mirror(double x) {
    if (!_front)
      return -x;
    return x;
  }

  int			_paper_dl;
  int			_front_texture_dl;
  int			_back_texture_dl;
  int			_bb_dl;

  int			_gm_size_x;
  int			_gm_size_y;

  bool			_front;
  bool			_two_sides;
  bool			_left_button_pressed;

  Layer			*_front_layer;
  Layer			*_back_layer;
  Layer			*_front_seam_layer;
  Layer         *_back_seam_layer;
  Layer         *_front_gfold_layer; // gaussian fold technique layer
  Layer         *_back_gfold_layer; 
  

  GarmentMapList	_front_garment_maps;
  GarmentMapList	_back_garment_maps;

  bool			_chain_pts;
  bool			_stroke_pts;
  bool			_texture;
  bool			_type;
  bool			_bb;
  bool			_seam_mode;
  bool          _gfold_mode;
  MirrorType          _vertical_mirror;
  

  MDIWindow		*_mdi;
};

#endif // _CANVAS_H
