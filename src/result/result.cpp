//
//  Filename         : result.cpp
//  Author           : Emmanuel Turquin
//  Purpose          : A viewer for the generated garments.
//  Date of creation : 04/20/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include <vector>
#include <lib3ds/mesh.h>
#include <lib3ds/vector.h>
#include <lib3ds/quat.h>
#include <iostream>
#include <fstream>
#include <qfiledialog.h>
#include <qcursor.h>
#include "mdi_window.h"
#include "config.h"
#include "repository_pretreatment_result.h"
#include "repository_canvas_result.h"
#include "pretreatment.h"
#include "canvas.h"
#include "layer.h"
#include "garment_map.h"
#include "point_utils.h"
#include "utils.h"
#include "config_result.h"
#include "result.h"

using namespace qglviewer;

Result::Result(MDIWindow* parent, const char* name,
			   const QGLWidget* shareWidget, int wflags)
  : BaseQGLViewer(parent, name, shareWidget, wflags)
{
  _mdi = parent;
  _model_dl = 0;
  _front_surface_dl = 0;
  _back_surface_dl = 0;
  _front_contour_dl = 0;
  _front_seam_contour_dl = 0;
  _back_seam_contour_dl = 0;
  _back_contour_dl = 0;
  _disp_surface_wireframe = false;
  updateConfig();
  RepositoryPretreatmentResult::setResult(this);
  RepositoryCanvasResult::setResult(this);
  _renderColour = QColor( 0, 200, 0);
}

Result::~Result()
{
  freeSurface();
  freeContour();
  if (RepositoryPretreatmentResult::result() == this)
    RepositoryPretreatmentResult::setResult(0);
  if (RepositoryCanvasResult::result() == this)
    RepositoryCanvasResult::setResult(0);
}

void Result::updateConfig()
{
  makeCurrent();

  bool type_changed = false;
  _disp_model = config::RESULT_DISP_MODEL;
  _disp_surface = config::RESULT_DISP_SURFACE;
  _disp_surface_wireframe = false;
  _disp_contour = config::RESULT_DISP_CONTOUR;
  _disp_type = config::RESULT_DISP_TYPE;
  if (_mdi && _mdi->config()) {
    _mdi->config()->io()->getValue("display/model", _disp_model);
    _mdi->config()->io()->getValue("display/surface", _disp_surface);
    _mdi->config()->io()->getValue("display/contour", _disp_contour);
    _mdi->config()->io()->getValue("display/type", _disp_type);
    _mdi->config()->io()->getValue("display/type_changed", type_changed);
    _mdi->config()->io()->getValue("display/surfaceRenderType", _disp_surface_wireframe);
    int red = 200; int green = 0; int blue = 0;
    _mdi->config()->io()->getValue("display/render/red", red);
    _mdi->config()->io()->getValue("display/render/green", green);
    _mdi->config()->io()->getValue("display/render/blue", blue);
    _renderColour.setRgb( red, green, blue );
    
    if (type_changed) {
      _mdi->config()->io()->setValue("display/type_changed", false);
      freeSurface();
      freeContour();
      _front_surface_dl = prepareFrontSurface();
      _back_surface_dl = prepareBackSurface();
      _front_contour_dl = prepareFrontContour();
      _front_seam_contour_dl = prepareFrontSeamContour();
      _back_seam_contour_dl = prepareBackSeamContour();
      _back_contour_dl = prepareBackContour();
    }
    updateGL();
  }
}

void Result::init()
{
  glEnable(GL_BLEND);
#ifdef GL_FUNC_ADD    // Phil
  glBlendEquation(GL_FUNC_ADD);
#endif
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glShadeModel(GL_SMOOTH);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glDepthFunc(GL_LESS);
  glEnable(GL_DEPTH_TEST);
  glDisable(GL_COLOR_MATERIAL);
  glEnable(GL_CULL_FACE);
  glEnable(GL_NORMALIZE);
  glCullFace(GL_BACK);
}

// Internal use only.
inline void vertexCoord(const GarmentMap *gm, double dx, double dy, unsigned i, unsigned j, Vec3d& pt)
{
  if (gm->type(i, j) == GarmentMap::IN) {
    pt.x() = gm->bbMin()[0] + dx/2 + i*dx;
    pt.y() = gm->bbMin()[1] + dy/2 + j*dy;
  }
  else { // BORDER or BORDER_EXT
    pt.x() = gm->borderCoord(i, j)->x();
    pt.y() = gm->borderCoord(i, j)->y();
  }
  pt.z() = gm->z(i, j);
}

// Internal use only.
inline void vertexColor(double x, double y, double dist = 0)
{
  switch (point_utils::type(x, y, dist)) {
  case Point::IN:
    glColor3f(0.3, 0.9, 0.3);
    break;
  case Point::OUT:
    glColor3f(0.9, 0.3, 0.3);
    break;
  case Point::BORDER:
    glColor3f(0.9, 0.6, 0.3);
    break;
  case Point::UNKNOWN:
  default:
    glColor3f(0.3, 0.3, 0.3);
  }
}

// Internal use only.
inline void vertexNormal(const GarmentMap* gm, double dx, double dy, double x, double y, unsigned i, unsigned j, Vec3d& normal)
{
  double di = 0, dj = 0;
  if (gm->type(i-1, j) != GarmentMap::IN)
    di += x - gm->borderCoord(i-1, j)->x();
  else
    di += dx;
  if (gm->type(i+1, j) != GarmentMap::IN)
    di += gm->borderCoord(i+1, j)->x() - x;
  else
    di += dx;
  if (gm->type(i, j-1) != GarmentMap::IN)
    dj += y - gm->borderCoord(i, j-1)->y();
  else
    dj += dy;
  if (gm->type(i, j+1) != GarmentMap::IN)
    dj += gm->borderCoord(i, j+1)->y() - y;
  else
    dj += dy;

  double dzx = gm->z(i+1, j) - gm->z(i-1, j);
  double dzy = gm->z(i, j+1) - gm->z(i, j-1);
  normal.x() = dzy - dzx*dj;
  normal.y() = -dzy*di;
  normal.z() = di*dj;
  
  
}

// Internal use only.
inline void drawQuad(const GarmentMap* gm, unsigned i, unsigned j, bool disp_type, bool front)
{
  double dx, dy;
  dx = gm->bbSize()[0] / gm->sizeX();
  dy = gm->bbSize()[1] / gm->sizeY();

  Vec3d pt1, pt2, pt3, pt4;
  Vec3d normal;

  vertexCoord(gm, dx, dy, i, j, pt1);
  vertexCoord(gm, dx, dy, i+1, j, pt2);
  vertexCoord(gm, dx, dy, i, j+1, pt3);
  vertexCoord(gm, dx, dy, i+1, j+1, pt4);

  glBegin(GL_TRIANGLE_STRIP);

  if (disp_type)
    vertexColor(pt1.x(), pt1.y(), gm->dist(i, j));
  if (gm->type(i, j) == GarmentMap::IN)
    vertexNormal(gm, dx, dy, pt1.x(), pt1.y(), i, j, normal);
  else
    normal = (pt2 - pt1) ^ (pt3 - pt1);
  if (!front)
    normal = normal * -1;
  glNormal3dv(normal.address());
  glVertex3dv(pt1.address());

  if (disp_type)
    vertexColor(pt2.x(), pt2.y(), gm->dist(i+1, j));
  if (gm->type(i+1, j) == GarmentMap::IN)
    vertexNormal(gm, dx, dy, pt2.x(), pt2.y(), i+1, j, normal);
  else
    normal = (pt3 - pt2) ^ (pt1 - pt2);
  if (!front)
    normal = normal * -1;
  glNormal3dv(normal.address());
  glVertex3dv(pt2.address());

  if (disp_type)
    vertexColor(pt3.x(), pt3.y(), gm->dist(i, j+1));
  if (gm->type(i, j+1) == GarmentMap::IN)
    vertexNormal(gm, dx, dy, pt3.x(), pt3.y(), i, j+1, normal);
  else
    normal = (pt2 - pt3) ^ (pt4 - pt3);
  if (!front)
    normal = normal * -1;
  glNormal3dv(normal.address());
  glVertex3dv(pt3.address());

  if (disp_type)
    vertexColor(pt4.x(), pt4.y(), gm->dist(i+1, j+1));
  if (gm->type(i+1, j+1) == GarmentMap::IN)
    vertexNormal(gm, dx, dy, pt4.x(), pt4.y(), i+1, j+1, normal);
  else
    normal = (pt3 - pt4) ^ (pt2 - pt4);
  if (!front)
    normal = normal * -1;
  glNormal3dv(normal.address());
  
  glVertex3dv(pt4.address());

  glEnd();

}

// Internal use only.
inline void drawTriangleSW(const GarmentMap* gm, unsigned i, unsigned j, bool disp_type, bool front)
{
  double dx, dy;
  dx = gm->bbSize()[0] / gm->sizeX();
  dy = gm->bbSize()[1] / gm->sizeY();

  Vec3d pt1, pt2, pt3;
  Vec3d normal;

  vertexCoord(gm, dx, dy, i, j+1, pt1);
  vertexCoord(gm, dx, dy, i, j, pt2);
  vertexCoord(gm, dx, dy, i+1, j, pt3);

  glBegin(GL_TRIANGLES);

  if (disp_type)
    vertexColor(pt1.x(), pt1.y(), gm->dist(i, j+1));
  if (gm->type(i, j+1) == GarmentMap::IN)
    vertexNormal(gm, dx, dy, pt1.x(), pt1.y(), i, j+1, normal);
  else
    normal = (pt2 - pt1) ^ (pt3 - pt1);
  if (!front)
    normal = normal * -1;
  glNormal3dv(normal.address());
  glVertex3dv(pt1.address());

  if (disp_type)
    vertexColor(pt2.x(), pt2.y(), gm->dist(i, j));
  if (gm->type(i, j) == GarmentMap::IN)
    vertexNormal(gm, dx, dy, pt2.x(), pt2.y(), i, j, normal);
  else
    normal = (pt3 - pt2) ^ (pt1 - pt2);
  if (!front)
    normal = normal * -1;
  glNormal3dv(normal.address());
  glVertex3dv(pt2.address());

  if (disp_type)
    vertexColor(pt3.x(), pt3.y(), gm->dist(i+1, j));
  if (gm->type(i+1, j) == GarmentMap::IN)
    vertexNormal(gm, dx, dy, pt3.x(), pt3.y(), i+1, j, normal);
  else
    normal = (pt1 - pt3) ^ (pt2 - pt3);
  if (!front)
    normal = normal * -1;
  glNormal3dv(normal.address());
  glVertex3dv(pt3.address());

  glEnd();
}

// Internal use only.
inline void drawTriangleSE(const GarmentMap* gm, unsigned i, unsigned j, bool disp_type, bool front)
{
  double dx, dy;
  dx = gm->bbSize()[0] / gm->sizeX();
  dy = gm->bbSize()[1] / gm->sizeY();

  Vec3d pt1, pt2, pt3;
  Vec3d normal;

  vertexCoord(gm, dx, dy, i, j, pt1);
  vertexCoord(gm, dx, dy, i+1, j, pt2);
  vertexCoord(gm, dx, dy, i+1, j+1, pt3);

  glBegin(GL_TRIANGLES);

  if (disp_type)
    vertexColor(pt1.x(), pt1.y(), gm->dist(i, j));
  if (gm->type(i, j) == GarmentMap::IN)
    vertexNormal(gm, dx, dy, pt1.x(), pt1.y(), i, j, normal);
  else
    normal = (pt2 - pt1) ^ (pt3 - pt1);
  if (!front)
    normal = normal * -1;
  glNormal3dv(normal.address());
  glVertex3dv(pt1.address());

  if (disp_type)
    vertexColor(pt2.x(), pt2.y(), gm->dist(i+1, j));
  if (gm->type(i+1, j) == GarmentMap::IN)
    vertexNormal(gm, dx, dy, pt2.x(), pt2.y(), i+1, j, normal);
  else
    normal = (pt3 - pt2) ^ (pt1 - pt2);
  if (!front)
    normal = normal * -1;
  glNormal3dv(normal.address());
  glVertex3dv(pt2.address());

  if (disp_type)
    vertexColor(pt3.x(), pt3.y(), gm->dist(i+1, j+1));
  if (gm->type(i+1, j+1) == GarmentMap::IN)
    vertexNormal(gm, dx, dy, pt3.x(), pt3.y(), i+1, j+1, normal);
  else
    normal = (pt1 - pt3) ^ (pt2 - pt3);
  if (!front)
    normal = normal * -1;
  glNormal3dv(normal.address());
  glVertex3dv(pt3.address());

  glEnd();
}

// Internal use only.
inline void drawTriangleNE(const GarmentMap* gm, unsigned i, unsigned j, bool disp_type, bool front)
{
  double dx, dy;
  dx = gm->bbSize()[0] / gm->sizeX();
  dy = gm->bbSize()[1] / gm->sizeY();

  Vec3d pt1, pt2, pt3;
  Vec3d normal;

  vertexCoord(gm, dx, dy, i+1, j, pt1);
  vertexCoord(gm, dx, dy, i+1, j+1, pt2);
  vertexCoord(gm, dx, dy, i, j+1, pt3);

  glBegin(GL_TRIANGLES);

  if (disp_type)
    vertexColor(pt1.x(), pt1.y(), gm->dist(i+1, j));
  if (gm->type(i+1, j) == GarmentMap::IN)
    vertexNormal(gm, dx, dy, pt1.x(), pt1.y(), i+1, j, normal);
  else
    normal = (pt2 - pt1) ^ (pt3 - pt1);
  if (!front)
    normal = normal * -1;
  glNormal3dv(normal.address());
  glVertex3dv(pt1.address());

  if (disp_type)
    vertexColor(pt2.x(), pt2.y(), gm->dist(i+1, j+1));
  if (gm->type(i+1, j+1) == GarmentMap::IN)
    vertexNormal(gm, dx, dy, pt2.x(), pt2.y(), i+1, j+1, normal);
  else
    normal = (pt3 - pt2) ^ (pt1 - pt2);
  if (!front)
    normal = normal * -1;
  glNormal3dv(normal.address());
  glVertex3dv(pt2.address());

  if (disp_type)
    vertexColor(pt3.x(), pt3.y(), gm->dist(i, j+1));
  if (gm->type(i, j+1) == GarmentMap::IN)
    vertexNormal(gm, dx, dy, pt3.x(), pt3.y(), i, j+1, normal);
  else
    normal = (pt1 - pt3) ^ (pt2 - pt3);
  if (!front)
    normal = normal * -1;
  glNormal3dv(normal.address());
  glVertex3dv(pt3.address());

  glEnd();
}

// Internal use only.
inline void drawTriangleNW(const GarmentMap* gm, unsigned i, unsigned j, bool disp_type, bool front)
{
  double dx, dy;
  dx = gm->bbSize()[0] / gm->sizeX();
  dy = gm->bbSize()[1] / gm->sizeY();

  Vec3d pt1, pt2, pt3;
  Vec3d normal;

  vertexCoord(gm, dx, dy, i+1, j+1, pt1);
  vertexCoord(gm, dx, dy, i, j+1, pt2);
  vertexCoord(gm, dx, dy, i, j, pt3);

  glBegin(GL_TRIANGLES);

  if (disp_type)
    vertexColor(pt1.x(), pt1.y(), gm->dist(i+1, j+1));
  if (gm->type(i+1, j+1) == GarmentMap::IN)
    vertexNormal(gm, dx, dy, pt1.x(), pt1.y(), i+1, j+1, normal);
  else
    normal = (pt2 - pt1) ^ (pt3 - pt1);
  if (!front)
    normal = normal * -1;
  glNormal3dv(normal.address());
  glVertex3dv(pt1.address());

  if (disp_type)
    vertexColor(pt2.x(), pt2.y(), gm->dist(i, j+1));
  if (gm->type(i, j+1) == GarmentMap::IN)
    vertexNormal(gm, dx, dy, pt2.x(), pt2.y(), i, j+1, normal);
  else
    normal = (pt3 - pt2) ^ (pt1 - pt2);
  if (!front)
    normal = normal * -1;
  glNormal3dv(normal.address());
  glVertex3dv(pt2.address());

  if (disp_type)
    vertexColor(pt3.x(), pt3.y(), gm->dist(i, j));
  if (gm->type(i, j) == GarmentMap::IN)
    vertexNormal(gm, dx, dy, pt3.x(), pt3.y(), i, j, normal);
  else
    normal = (pt1 - pt3) ^ (pt2 - pt3);
  if (!front)
    normal = normal * -1;
  glNormal3dv(normal.address());
  glVertex3dv(pt3.address());

  glEnd();
}

void Result::draw()
{
  if (_disp_model)
    renderModel();


  
  if (_disp_surface)
  {
    glPushAttrib( GL_POLYGON_BIT );
    
    if(_disp_surface_wireframe) 
    {
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    }
    renderSurface();
    glPopAttrib();
  }
  
  glPopAttrib();
  
  if (_disp_contour)
    renderContour();
}

void Result::renderModel()
{
  if (_model_dl)
    glCallList(_model_dl);
}

void Result::renderContour()
{
  glPushMatrix();
  glMultMatrixd(RepositoryPretreatmentResult::matrix().address());
  if (_front_contour_dl)
    glCallList(_front_contour_dl);
  if (_front_seam_contour_dl)
    glCallList(_front_seam_contour_dl);
  if (_back_seam_contour_dl)
    glCallList(_back_seam_contour_dl);
  if (_back_contour_dl)
    glCallList(_back_contour_dl);
  glPopMatrix();
}

void Result::freeContour()
{
  if (_front_contour_dl) {
    glDeleteLists(_front_contour_dl, 1);
    _front_contour_dl = 0;
  }
  if (_front_seam_contour_dl) {
    glDeleteLists(_front_seam_contour_dl, 1);
    _front_seam_contour_dl = 0;
  }
  if (_back_seam_contour_dl) {
    glDeleteLists(_back_seam_contour_dl, 1);
    _back_seam_contour_dl = 0;
  }
  if (_back_contour_dl) {
    glDeleteLists(_back_contour_dl, 1);
    _back_contour_dl = 0;
  }
}

int Result::prepareFrontContour()
{
  if (!RepositoryCanvasResult::frontLayer()) // FIXME
  {
    std::cout << "Couldn't prepare front contour" << std::endl;
    return 0;
  }

  int contour_dl = glGenLists(1);
  if (!contour_dl)
    return 0;

  glNewList(contour_dl, GL_COMPILE);

  Layer::ChainList::const_iterator ch;
  Layer::ChainList::const_iterator ch_end;
  Chain::StrokeList::const_iterator st;
  Chain::StrokeList::const_iterator st_end;
  Stroke::SegmentList::const_iterator seg;
  Stroke::SegmentList::const_iterator seg_end;

  glDisable(GL_LIGHTING);

  glPointSize(5);
  
  ch = RepositoryCanvasResult::frontLayer()->chains()->begin(); // FIXME
  ch_end = RepositoryCanvasResult::frontLayer()->chains()->end(); // FIXME
  for (; ch != ch_end; ++ch) {
    st = (*ch)->strokes()->begin();
    st_end = (*ch)->strokes()->end();
    for (; st != st_end; ++st) {
      seg = (*st)->segments()->begin();
      seg_end = (*st)->segments()->end();
      glBegin(GL_LINES);
      for (; seg != seg_end; ++seg) {
	
	if (_disp_type) {
	  switch ((*seg)->type()) {
	  case Segment::IN:
	    glColor3f(0.3, 0.9, 0.3);
	    break;
	  case Segment::OUT:
	    glColor3f(0.9, 0.3, 0.3);
	    break;
	  case Segment::BORDER:
	    glColor3f(0.9, 0.6, 0.3);
	    break;
	  case Segment::MIXED:
	    glColor3f(0.9, 0.9, 0.3);
	    break;
	  case Segment::UNKNOWN:
	  default:
	    glColor3f(0.3, 0.3, 0.3);
	  }
	}
	
	glVertex3d((*seg)->pointA()->x(), (*seg)->pointA()->y(), (*seg)->pointA()->z());
	glVertex3d((*seg)->pointB()->x(), (*seg)->pointB()->y(), (*seg)->pointB()->z());
      }
      glEnd();
    }
  }
  
  counted_ptr<Point> p;
  
  ch = RepositoryCanvasResult::frontLayer()->chains()->begin(); // FIXME
  ch_end = RepositoryCanvasResult::frontLayer()->chains()->end(); // FIXME
  glBegin(GL_POINTS);
  for (; ch != ch_end; ++ch) {
    st = (*ch)->strokes()->begin();
    st_end = (*ch)->strokes()->end();
    for (; st != st_end; ++st) {
      vertexColor((*st)->pointA()->x(), (*st)->pointA()->y(), (*st)->pointA()->distance());
      glVertex3dv((*st)->pointA()->address());
      vertexColor((*st)->pointB()->x(), (*st)->pointB()->y(), (*st)->pointB()->distance());
      glVertex3dv((*st)->pointB()->address());
    }
  }
  glEnd();
    
  glEnable(GL_LIGHTING);

  glEndList();
  return contour_dl;
}

int Result::prepareFrontSeamContour()
{
  if (!RepositoryCanvasResult::frontSeamLayer()) {
    std::cout << "Couldn't prepare front seam contour" << std::endl;
    return 0;
  }
  
  int contour_dl = glGenLists(1);
  if (!contour_dl)
    return 0;

  glNewList(contour_dl, GL_COMPILE);

  Layer::ChainList::const_iterator ch;
  Layer::ChainList::const_iterator ch_end;
  Chain::StrokeList::const_iterator st;
  Chain::StrokeList::const_iterator st_end;
  Stroke::SegmentList::const_iterator seg;
  Stroke::SegmentList::const_iterator seg_end;

  glDisable(GL_LIGHTING);

  glPointSize(5);
  
  ch = RepositoryCanvasResult::frontSeamLayer()->chains()->begin(); 
  ch_end = RepositoryCanvasResult::frontSeamLayer()->chains()->end(); 
  for (; ch != ch_end; ++ch) {
    if((*ch)->cycle()) // shouldn't be any cyclic chains
      continue;
    
    st = (*ch)->strokes()->begin();
    st_end = (*ch)->strokes()->end();
    for (; st != st_end; ++st) {
      seg = (*st)->segments()->begin();
      seg_end = (*st)->segments()->end();
      glBegin(GL_LINES);
      for (; seg != seg_end; ++seg) {
        glColor3f(0.8, 0.0, 0.0);
    
        glVertex3d((*seg)->pointA()->x(), (*seg)->pointA()->y(), (*seg)->pointA()->z());
        glVertex3d((*seg)->pointB()->x(), (*seg)->pointB()->y(), (*seg)->pointB()->z());
      }
      glEnd();
    }
  }
  
  glEnd();
    
  glEnable(GL_LIGHTING);

  glEndList();
  return contour_dl;
}


int Result::prepareBackSeamContour()
{
  if (!RepositoryCanvasResult::backSeamLayer()) {
    std::cout << "Couldn't prepare back seam contour" << std::endl;
    return 0;
  }
  
  int contour_dl = glGenLists(1);
  if (!contour_dl)
    return 0;

  glNewList(contour_dl, GL_COMPILE);

  Layer::ChainList::const_iterator ch;
  Layer::ChainList::const_iterator ch_end;
  Chain::StrokeList::const_iterator st;
  Chain::StrokeList::const_iterator st_end;
  Stroke::SegmentList::const_iterator seg;
  Stroke::SegmentList::const_iterator seg_end;

  glDisable(GL_LIGHTING);

  glPointSize(5);
  
  ch = RepositoryCanvasResult::backSeamLayer()->chains()->begin(); 
  ch_end = RepositoryCanvasResult::backSeamLayer()->chains()->end(); 
  for (; ch != ch_end; ++ch) {
    if((*ch)->cycle()) // shouldn't be any cyclic chains
      continue;
    
    st = (*ch)->strokes()->begin();
    st_end = (*ch)->strokes()->end();
    for (; st != st_end; ++st) {
      seg = (*st)->segments()->begin();
      seg_end = (*st)->segments()->end();
      glBegin(GL_LINES);
      for (; seg != seg_end; ++seg) {
        glColor3f(0.8, 0.0, 0.0);
    
        glVertex3d((*seg)->pointA()->x(), (*seg)->pointA()->y(), (*seg)->pointA()->z());
        glVertex3d((*seg)->pointB()->x(), (*seg)->pointB()->y(), (*seg)->pointB()->z());
      }
      glEnd();
    }
  }
  
  glEnd();
    
  glEnable(GL_LIGHTING);

  glEndList();
  return contour_dl;
}

int Result::prepareBackContour()
{
  if (!RepositoryCanvasResult::backLayer()) {// FIXME
    std::cout << "Couldn't prepare back contour" << std::endl;
    return 0;
  }
  int contour_dl = glGenLists(1);
  if (!contour_dl)
    return 0;

  glNewList(contour_dl, GL_COMPILE);

  Layer::ChainList::const_iterator ch;
  Layer::ChainList::const_iterator ch_end;
  Chain::StrokeList::const_iterator st;
  Chain::StrokeList::const_iterator st_end;
  Stroke::SegmentList::const_iterator seg;
  Stroke::SegmentList::const_iterator seg_end;

  glDisable(GL_LIGHTING);

  glPointSize(5);
  
  ch = RepositoryCanvasResult::backLayer()->chains()->begin(); // FIXME
  ch_end = RepositoryCanvasResult::backLayer()->chains()->end(); // FIXME
  for (; ch != ch_end; ++ch) {
    st = (*ch)->strokes()->begin();
    st_end = (*ch)->strokes()->end();
    for (; st != st_end; ++st) {
      seg = (*st)->segments()->begin();
      seg_end = (*st)->segments()->end();
      glBegin(GL_LINES);
      for (; seg != seg_end; ++seg) {
	
	if (_disp_type) {
	  switch ((*seg)->type()) {
	  case Segment::IN:
	    glColor3f(0.3, 0.9, 0.3);
	    break;
	  case Segment::OUT:
	    glColor3f(0.9, 0.3, 0.3);
	    break;
	  case Segment::BORDER:
	    glColor3f(0.9, 0.6, 0.3);
	    break;
	  case Segment::MIXED:
	    glColor3f(0.9, 0.9, 0.3);
	    break;
	  case Segment::UNKNOWN:
	  default:
	    glColor3f(0.3, 0.3, 0.3);
	  }
	}
	
	glVertex3d((*seg)->pointA()->x(), (*seg)->pointA()->y(), (*seg)->pointA()->z());
	glVertex3d((*seg)->pointB()->x(), (*seg)->pointB()->y(), (*seg)->pointB()->z());
      }
      glEnd();
    }
  }
  
  counted_ptr<Point> p;
  
  ch = RepositoryCanvasResult::backLayer()->chains()->begin(); // FIXME
  ch_end = RepositoryCanvasResult::backLayer()->chains()->end(); // FIXME
  glBegin(GL_POINTS);
  for (; ch != ch_end; ++ch) {
    st = (*ch)->strokes()->begin();
    st_end = (*ch)->strokes()->end();
    for (; st != st_end; ++st) {
      vertexColor((*st)->pointA()->x(), (*st)->pointA()->y(), (*st)->pointA()->distance());
      glVertex3dv((*st)->pointA()->address());
      vertexColor((*st)->pointB()->x(), (*st)->pointB()->y(), (*st)->pointB()->distance());
      glVertex3dv((*st)->pointB()->address());
    }
  }
  glEnd();
    
  glEnable(GL_LIGHTING);

  glEndList();
  return contour_dl;
}

void Result::renderSurface()
{
  glPushMatrix();
  glMultMatrixd(RepositoryPretreatmentResult::matrix().address());
  if (_front_surface_dl)
    glCallList(_front_surface_dl);
  if (_back_surface_dl)
    glCallList(_back_surface_dl);
  glPopMatrix();
}

void Result::freeSurface()
{
  if (_front_surface_dl) {
    glDeleteLists(_front_surface_dl, 1);
    _front_surface_dl = 0;
  }
  if (_back_surface_dl) {
    glDeleteLists(_back_surface_dl, 1);
    _back_surface_dl = 0;
  }
}

int Result::prepareFrontSurface()
{
  if (!RepositoryCanvasResult::frontLayer()) // FIXME
    return 0;

  int surface_dl = glGenLists(1);
  if (!surface_dl)
    return 0;


  
  glNewList(surface_dl, GL_COMPILE);
  glDisable(GL_CULL_FACE);
  //float no_shininess[] = {0.0, 0.0, 0.0, 0.0}; // FIXME
  //glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, no_shininess);
  
  //GLfloat light_ambient[] = { 0.5, 0.5, 0.5, 1.0 };
  //GLfloat light_diffuse[] = { 0.5, 0.5, 0.5, 1.0 };
  //GLfloat no_mat[] = { 0.0, 0.0, 0.0, 1.0 };
  //GLfloat light_position[] = { 1.0, 1.0, 1.0, 0.0 };

  //glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
  //glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
  //glLightfv(GL_LIGHT0, GL_SPECULAR, no_mat);
  //glLightfv(GL_LIGHT0, GL_POSITION, light_position);
  
  
  GLfloat mat_specular[] = { 0.0, 0.0, 0.0,  1.0 };
  //GLfloat mat_diffuse[] = { 0.0, 0.0, 0.0, 1.0 };
  //GLfloat mat_shininess[] = { 0.0 };
  
  //mat_diffuse[0] = _renderColour.red() / 255.;
  //mat_diffuse[1] = _renderColour.green() / 255.;
  //mat_diffuse[2] = _renderColour.blue() / 255.;

  glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
  //glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
  //glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mat_diffuse);
  glEnable(GL_COLOR_MATERIAL);
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  
  //glColor3f(0.9, 0.3, 0.3);
  
  glEnable(GL_LIGHT0);
  drawLight(GL_LIGHT0);
  glDisable(GL_LIGHT1);
  glDisable(GL_LIGHT2);
  
  glColor3f( _renderColour.red() / 255., _renderColour.green() / 255., _renderColour.blue() / 255.);
  
  if (RepositoryCanvasResult::frontGarmentMaps()) {
    Canvas::GarmentMapList::const_iterator gm = RepositoryCanvasResult::frontGarmentMaps()->begin();
    Canvas::GarmentMapList::const_iterator gm_end = RepositoryCanvasResult::frontGarmentMaps()->end();
    bool sw, se, nw, ne;
    unsigned count, count_in;

    
    for (; gm != gm_end; ++gm) {
      for (unsigned i = 0; i < (*gm)->sizeX() - 1; ++i) {
	for (unsigned j = 0; j < (*gm)->sizeY() - 1; ++j) {
	  count = 0;
	  count_in = 0;
	  if ((*gm)->type(i, j) != GarmentMap::OUT && (*gm)->type(i, j) != GarmentMap::UNKNOWN) {
	    sw = true;
	    ++count;
	    if ((*gm)->type(i, j) == GarmentMap::IN)
	      ++count_in;
	  }
	  else {
	    sw = false;
	  }
	  if((*gm)->type(i+1, j) != GarmentMap::OUT && (*gm)->type(i, j) != GarmentMap::UNKNOWN) {
	    se = true;
	    ++count;
	    if ((*gm)->type(i+1, j) == GarmentMap::IN)
	      ++count_in;
	  }
	  else {
	    se = false;
	  }
	  if ((*gm)->type(i, j+1) != GarmentMap::OUT && (*gm)->type(i, j) != GarmentMap::UNKNOWN) {
	    nw = true;
	    ++count;
	    if ((*gm)->type(i, j+1) == GarmentMap::IN)
	      ++count_in;
	  }
	  else {
	    nw = false;
	  }
	  if ((*gm)->type(i+1, j+1) != GarmentMap::OUT && (*gm)->type(i, j) != GarmentMap::UNKNOWN) {
	    ne = true;
	    ++count;
	    if ((*gm)->type(i+1, j+1) == GarmentMap::IN)
	      ++count_in;
	  }
	  else {
	    ne = false;
	  }

	  switch (count) {
	  case 4:
	    if (count_in == 4)
	      drawQuad(*gm, i, j, _disp_type, true);
	    else {
	      drawTriangleSW(*gm, i, j, _disp_type, true);
	      drawTriangleNE(*gm, i, j, _disp_type, true);
	    }
	    break;
	  case 3:
	    if (!count_in) {
	      if ((*gm)->type(i, j) == GarmentMap::OUT && (*gm)->type(i+1, j+1) == GarmentMap::BORDER_EXT ||
		  (*gm)->type(i+1, j+1) == GarmentMap::OUT && (*gm)->type(i, j) == GarmentMap::BORDER_EXT ||
		  (*gm)->type(i, j+1) == GarmentMap::OUT && (*gm)->type(i+1, j) == GarmentMap::BORDER_EXT ||
		  (*gm)->type(i+1, j) == GarmentMap::OUT && (*gm)->type(i, j+1) == GarmentMap::BORDER_EXT)
		continue;
	    }
	    if (nw && sw && se) // SW
	      drawTriangleSW(*gm, i, j, _disp_type, true);
	    else if (sw && se && ne) // SE
	      drawTriangleSE(*gm, i, j, _disp_type, true);
	    else if (se && ne && nw) // NE
	      drawTriangleNE(*gm, i, j, _disp_type, true);
	    else // NW
	      drawTriangleNW(*gm, i, j, _disp_type, true);
	    break;
	  default:
	    ;
	  }
	}
      }
    }
  }
  glDisable(GL_COLOR_MATERIAL);
  glEnable(GL_CULL_FACE);
  


  glEndList();
  return surface_dl;
}

int Result::prepareBackSurface()
{
  if (!RepositoryCanvasResult::backLayer()) // FIXME
    return 0;

  int surface_dl = glGenLists(1);
  if (!surface_dl)
    return 0;
  


  glNewList(surface_dl, GL_COMPILE);
  glDisable(GL_CULL_FACE);
  float no_shininess[] = {0.0, 0.0, 0.0, 1.0}; // FIXME
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, no_shininess);
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  //glColor3f(0.9, 0.3, 0.3);
  glColor3f( _renderColour.red() / 255., _renderColour.green() / 255., _renderColour.blue() / 255.);
  glEnable(GL_COLOR_MATERIAL);
  if (RepositoryCanvasResult::backGarmentMaps()) {
    Canvas::GarmentMapList::const_iterator gm = RepositoryCanvasResult::backGarmentMaps()->begin();
    Canvas::GarmentMapList::const_iterator gm_end = RepositoryCanvasResult::backGarmentMaps()->end();
    bool sw, se, nw, ne;
    unsigned count, count_in;
    for (; gm != gm_end; ++gm) {
      for (unsigned i = 0; i < (*gm)->sizeX() - 1; ++i) {
	for (unsigned j = 0; j < (*gm)->sizeY() - 1; ++j) {
	  count = 0;
	  count_in = 0;
	  if ((*gm)->type(i, j) != GarmentMap::OUT && (*gm)->type(i, j) != GarmentMap::UNKNOWN) {
	    sw = true;
	    ++count;
	    if ((*gm)->type(i, j) == GarmentMap::IN)
	      ++count_in;
	  }
	  else {
	    sw = false;
	  }
	  if((*gm)->type(i+1, j) != GarmentMap::OUT && (*gm)->type(i, j) != GarmentMap::UNKNOWN) {
	    se = true;
	    ++count;
	    if ((*gm)->type(i+1, j) == GarmentMap::IN)
	      ++count_in;
	  }
	  else {
	    se = false;
	  }
	  if ((*gm)->type(i, j+1) != GarmentMap::OUT && (*gm)->type(i, j) != GarmentMap::UNKNOWN) {
	    nw = true;
	    ++count;
	    if ((*gm)->type(i, j+1) == GarmentMap::IN)
	      ++count_in;
	  }
	  else {
	    nw = false;
	  }
	  if ((*gm)->type(i+1, j+1) != GarmentMap::OUT && (*gm)->type(i, j) != GarmentMap::UNKNOWN) {
	    ne = true;
	    ++count;
	    if ((*gm)->type(i+1, j+1) == GarmentMap::IN)
	      ++count_in;
	  }
	  else {
	    ne = false;
	  }

	  switch (count) {
	  case 4:
	    if (count_in == 4)
	      drawQuad(*gm, i, j, _disp_type, false);
	    else {
	      drawTriangleSW(*gm, i, j, _disp_type, false);
	      drawTriangleNE(*gm, i, j, _disp_type, false);
	    }
	    break;
	  case 3:
	    if (!count_in) {
	      if ((*gm)->type(i, j) == GarmentMap::OUT && (*gm)->type(i+1, j+1) == GarmentMap::BORDER_EXT ||
		  (*gm)->type(i+1, j+1) == GarmentMap::OUT && (*gm)->type(i, j) == GarmentMap::BORDER_EXT ||
		  (*gm)->type(i, j+1) == GarmentMap::OUT && (*gm)->type(i+1, j) == GarmentMap::BORDER_EXT ||
		  (*gm)->type(i+1, j) == GarmentMap::OUT && (*gm)->type(i, j+1) == GarmentMap::BORDER_EXT)
		continue;
	    }
	    if (nw && sw && se) // SW
	      drawTriangleSW(*gm, i, j, _disp_type, false);
	    else if (sw && se && ne) // SE
	      drawTriangleSE(*gm, i, j, _disp_type, false);
	    else if (se && ne && nw) // NE
	      drawTriangleNE(*gm, i, j, _disp_type, false);
	    else // NW
	      drawTriangleNW(*gm, i, j, _disp_type, false);
	    break;
	  default:
	    ;
	  }
	}
      }
    }
  }
  glDisable(GL_COLOR_MATERIAL);
  glEnable(GL_CULL_FACE);
 

  glEndList();
  return surface_dl;
}

void Result::pretreatmentUpdated()
{
  unsigned modified = RepositoryPretreatmentResult::modified();
  
  if (modified & RepositoryPretreatmentResult::MODEL_DL) {
    freeSurface();
    freeContour();
    _model_dl = RepositoryPretreatmentResult::modelDisplayList();
    if (_model_dl) {
      setSceneBoundingBox(RepositoryPretreatmentResult::modelBBMin().address(),
			  RepositoryPretreatmentResult::modelBBMax().address());
      camera()->showEntireScene();
      
      message(MsgHandler::MSG_NORMAL, name(), QString("Model received from ") +
	      RepositoryPretreatmentResult::pretreatment()->name() + ".");
    }
  }
  updateGL();
}

void Result::canvasUpdated()
{
  freeSurface();
  freeContour();
  _front_surface_dl = prepareFrontSurface();
  _back_surface_dl = prepareBackSurface();
  _front_contour_dl = prepareFrontContour();
  _front_seam_contour_dl = prepareFrontSeamContour();
  _back_seam_contour_dl = prepareBackSeamContour();
  _back_contour_dl = prepareBackContour();
  updateGL();
}

void Result::saveFile() // FIXME
{
  
  if (
      (!RepositoryCanvasResult::frontGarmentMaps() ||
       RepositoryCanvasResult::frontGarmentMaps()->empty())
      &&
      (!RepositoryCanvasResult::backGarmentMaps() ||
       RepositoryCanvasResult::backGarmentMaps()->empty())
      ) {
    emit statusBarMessage("No mesh to save.", 2000);
    emit message(MsgHandler::MSG_ERROR, name(), "No mesh to save.");
    return;
  }

  QString path;
  
  if (_mdi && _mdi->config())
    _mdi->config()->io()->getValue("paths/file_3ds", path);
  QString filename = QFileDialog::getSaveFileName(path, "3DS files (*.3ds *.3DS);;All files (*)", this);
  
  // In case of Cancel.
  if (filename.isEmpty())
    return;

  QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
  
  // Build vertices and triangles lists.
  double dx, dy;
  unsigned *indices;
  unsigned current_index;
  Vec3d pt;
  std::vector<Vec3d> vertices;
  std::vector<Vec3ui> triangles;
  bool sw, se, nw, ne;
  unsigned count, count_in;

  // Front part.
  Canvas::GarmentMapList::const_iterator gm = RepositoryCanvasResult::frontGarmentMaps()->begin();
  Canvas::GarmentMapList::const_iterator gm_end = RepositoryCanvasResult::frontGarmentMaps()->end();

  current_index = 0;
  for (; gm != gm_end; ++gm) {
    dx = (*gm)->bbSize()[0] / (*gm)->sizeX();
    dy = (*gm)->bbSize()[1] / (*gm)->sizeY();
    indices = new unsigned[(*gm)->sizeX() * (*gm)->sizeY()];
    for (unsigned i = 0; i < (*gm)->sizeX(); ++i)
    {
      for (unsigned j = 0; j < (*gm)->sizeY(); ++j)
      {
        if ((*gm)->type(i, j) != GarmentMap::OUT && (*gm)->type(i, j) != GarmentMap::UNKNOWN)
        {
          indices[i * (*gm)->sizeY() + j] = current_index++;
          vertexCoord(*gm, dx, dy, i, j, pt);
          vertices.push_back(pt);
        }
      }
    }

    for (unsigned i = 0; i < (*gm)->sizeX() - 1; ++i) {
      for (unsigned j = 0; j < (*gm)->sizeY() - 1; ++j) {
	count = 0;
	count_in = 0;
	if ((*gm)->type(i, j) != GarmentMap::OUT && (*gm)->type(i, j) != GarmentMap::UNKNOWN) {
	  sw = true;
	  ++count;
	  if ((*gm)->type(i, j) == GarmentMap::IN)
	    ++count_in;
	}
	else {
	  sw = false;
	}
	if((*gm)->type(i+1, j) != GarmentMap::OUT && (*gm)->type(i, j) != GarmentMap::UNKNOWN) {
	  se = true;
	  ++count;
	  if ((*gm)->type(i+1, j) == GarmentMap::IN)
	    ++count_in;
	}
	else {
	  se = false;
	}
	if ((*gm)->type(i, j+1) != GarmentMap::OUT && (*gm)->type(i, j) != GarmentMap::UNKNOWN) {
	  nw = true;
	  ++count;
	  if ((*gm)->type(i, j+1) == GarmentMap::IN)
	    ++count_in;
	}
	else {
	  nw = false;
	}
	if ((*gm)->type(i+1, j+1) != GarmentMap::OUT && (*gm)->type(i, j) != GarmentMap::UNKNOWN) {
	  ne = true;
	  ++count;
    if ((*gm)->type(i+1, j+1) == GarmentMap::IN)
	    ++count_in;
	}
	else {
	  ne = false;
	}

	switch (count) {
	case 4:
	  triangles.push_back(Vec3ui(indices[i * (*gm)->sizeY() + j],
				     indices[(i+1) * (*gm)->sizeY() + j],
				     indices[i * (*gm)->sizeY() + (j+1)]));
	  triangles.push_back(Vec3ui(indices[(i+1) * (*gm)->sizeY() + j],
				     indices[(i+1) * (*gm)->sizeY() + (j+1)],
				     indices[i * (*gm)->sizeY() + (j+1)]));
	  break;
	case 3:
	  if (!count_in) {
	    if ((*gm)->type(i, j) == GarmentMap::OUT && (*gm)->type(i+1, j+1) == GarmentMap::BORDER_EXT ||
		(*gm)->type(i+1, j+1) == GarmentMap::OUT && (*gm)->type(i, j) == GarmentMap::BORDER_EXT ||
		(*gm)->type(i, j+1) == GarmentMap::OUT && (*gm)->type(i+1, j) == GarmentMap::BORDER_EXT ||
		(*gm)->type(i+1, j) == GarmentMap::OUT && (*gm)->type(i, j+1) == GarmentMap::BORDER_EXT)
	      continue;
	  }
	  if (nw && sw && se) // SW
	    triangles.push_back(Vec3ui(indices[i * (*gm)->sizeY() + j],
				       indices[(i+1) * (*gm)->sizeY() + j],
				       indices[i * (*gm)->sizeY() + (j+1)]));
	  else if (sw && se && ne) // SE
	    triangles.push_back(Vec3ui(indices[i * (*gm)->sizeY() + j],
				       indices[(i+1) * (*gm)->sizeY() + j],
				       indices[(i+1) * (*gm)->sizeY() + (j+1)]));
	  else if (se && ne && nw) // NE
	    triangles.push_back(Vec3ui(indices[(i+1) * (*gm)->sizeY() + j],
				       indices[(i+1) * (*gm)->sizeY() + (j+1)],
				       indices[i * (*gm)->sizeY() + (j+1)]));
	  else // NW
	    triangles.push_back(Vec3ui(indices[i * (*gm)->sizeY() + j],
				       indices[(i+1) * (*gm)->sizeY() + (j+1)],
				       indices[i * (*gm)->sizeY() + (j+1)]));
	  break;
	default:
	  ;
	}
      }
    }

    delete indices;
  }
  
  // JDW - write an obj file for the front mesh
  unsigned previousVertexCount = vertices.size();
  unsigned previousTriangleCount = triangles.size();
  QString obj_file_name = filename;
  obj_file_name.replace(".3ds",".obj");
  if(obj_file_name == filename) { // there wasn't a 3ds extension
    obj_file_name += QString(".obj");
  }
  
  
  
  std::ofstream _obj_file;
  _obj_file.open(obj_file_name.latin1());
  _obj_file << "g front" << std::endl;
  for(unsigned i=0;i<vertices.size();i++) {
    _obj_file  << "v " 
        << vertices.at(i).x() << " "
        << vertices.at(i).y() << " "
        << vertices.at(i).z() << std::endl;
  }
  for(unsigned i=0;i<triangles.size();i++) {
    _obj_file  << "f " 
        << triangles.at(i).x()+1 << " " // plus one because indices in OBJ files are 1 based, not 0 based
        << triangles.at(i).y()+1 << " "
        << triangles.at(i).z()+1 << std::endl;
  }
  

  //

  // Back part.
  gm = RepositoryCanvasResult::backGarmentMaps()->begin();
  gm_end = RepositoryCanvasResult::backGarmentMaps()->end();

  for (; gm != gm_end; ++gm) {
    dx = (*gm)->bbSize()[0] / (*gm)->sizeX();
    dy = (*gm)->bbSize()[1] / (*gm)->sizeY();
    indices = new unsigned[(*gm)->sizeX() * (*gm)->sizeY()];
    for (unsigned i = 0; i < (*gm)->sizeX(); ++i) {
      for (unsigned j = 0; j < (*gm)->sizeY(); ++j) {
	if ((*gm)->type(i, j) != GarmentMap::OUT && (*gm)->type(i, j) != GarmentMap::UNKNOWN) {
	  indices[i * (*gm)->sizeY() + j] = current_index++;
	  vertexCoord(*gm, dx, dy, i, j, pt);
	  vertices.push_back(pt);
	}
      }
    }

    for (unsigned i = 0; i < (*gm)->sizeX() - 1; ++i) {
      for (unsigned j = 0; j < (*gm)->sizeY() - 1; ++j) {
	count = 0;
	count_in = 0;
	if ((*gm)->type(i, j) != GarmentMap::OUT && (*gm)->type(i, j) != GarmentMap::UNKNOWN) {
	  sw = true;
	  ++count;
	  if ((*gm)->type(i, j) == GarmentMap::IN)
	    ++count_in;
	}
	else {
	  sw = false;
	}
	if((*gm)->type(i+1, j) != GarmentMap::OUT && (*gm)->type(i, j) != GarmentMap::UNKNOWN) {
	  se = true;
	  ++count;
	  if ((*gm)->type(i+1, j) == GarmentMap::IN)
	    ++count_in;
	}
	else {
	  se = false;
	}
	if ((*gm)->type(i, j+1) != GarmentMap::OUT && (*gm)->type(i, j) != GarmentMap::UNKNOWN) {
	  nw = true;
	  ++count;
	  if ((*gm)->type(i, j+1) == GarmentMap::IN)
	    ++count_in;
	}
	else {
	  nw = false;
	}
	if ((*gm)->type(i+1, j+1) != GarmentMap::OUT && (*gm)->type(i, j) != GarmentMap::UNKNOWN) {
	  ne = true;
	  ++count;
	  if ((*gm)->type(i+1, j+1) == GarmentMap::IN)
	    ++count_in;
	}
	else {
	  ne = false;
	}

	switch (count) {
	case 4:
	  triangles.push_back(Vec3ui(indices[i * (*gm)->sizeY() + (j+1)],
				     indices[(i+1) * (*gm)->sizeY() + j],
				     indices[i * (*gm)->sizeY() + j]));
	  triangles.push_back(Vec3ui(indices[i * (*gm)->sizeY() + (j+1)],
				     indices[(i+1) * (*gm)->sizeY() + (j+1)],
				     indices[(i+1) * (*gm)->sizeY() + j]));
	  break;
	case 3:
	  if (!count_in) {
	    if ((*gm)->type(i, j) == GarmentMap::OUT && (*gm)->type(i+1, j+1) == GarmentMap::BORDER_EXT ||
		(*gm)->type(i+1, j+1) == GarmentMap::OUT && (*gm)->type(i, j) == GarmentMap::BORDER_EXT ||
		(*gm)->type(i, j+1) == GarmentMap::OUT && (*gm)->type(i+1, j) == GarmentMap::BORDER_EXT ||
		(*gm)->type(i+1, j) == GarmentMap::OUT && (*gm)->type(i, j+1) == GarmentMap::BORDER_EXT)
	      continue;
	  }
	  if (nw && sw && se) // SW
	    triangles.push_back(Vec3ui(indices[i * (*gm)->sizeY() + (j+1)],
				       indices[(i+1) * (*gm)->sizeY() + j],
				       indices[i * (*gm)->sizeY() + j]));
	  else if (sw && se && ne) // SE
	    triangles.push_back(Vec3ui(indices[(i+1) * (*gm)->sizeY() + (j+1)],
				       indices[(i+1) * (*gm)->sizeY() + j],
				       indices[i * (*gm)->sizeY() + j]));
	  else if (se && ne && nw) // NE
	    triangles.push_back(Vec3ui(indices[i * (*gm)->sizeY() + (j+1)],
				       indices[(i+1) * (*gm)->sizeY() + (j+1)],
				       indices[(i+1) * (*gm)->sizeY() + j]));
	  else // NW
	    triangles.push_back(Vec3ui(indices[i * (*gm)->sizeY() + (j+1)],
				       indices[(i+1) * (*gm)->sizeY() + (j+1)],
				       indices[i * (*gm)->sizeY() + j]));
	  break;
	default:
	  ;
	}
      }
    }

    delete indices;
  }
  
  _obj_file << "g back" << std::endl;
  for(unsigned i=previousVertexCount;i<vertices.size();i++) {
    _obj_file  << "v " 
        << vertices.at(i).x() << " "
        << vertices.at(i).y() << " "
        << vertices.at(i).z() << std::endl;
  }
  for(unsigned i=previousTriangleCount;i<triangles.size();i++) {
    _obj_file  << "f " 
        << triangles.at(i).x()+1 << " "
        << triangles.at(i).y()+1 << " "
        << triangles.at(i).z()+1 << std::endl;
  }
  previousVertexCount = vertices.size();
  

  std::cout << "Testing seam layer" << std::endl;
  
  Canvas::ListOfPointLists::iterator lol_it = RepositoryCanvasResult::canvas()->getFrontSeams()->begin();
  Canvas::ListOfPointLists::iterator lol_it_end = RepositoryCanvasResult::canvas()->getFrontSeams()->end();

  int seamNumber=0;
  for(; lol_it != lol_it_end; ++lol_it) {
    std::list<Point>::iterator spl = lol_it->begin();
    std::list<Point>::iterator spl_end = lol_it->end();
    
    _obj_file<<"g fseam"<<++seamNumber << std::endl;
    for(; spl != spl_end; ++spl) {
      _obj_file << "v " // OBJ file format for vertices
          << spl->x() << " "
          << spl->y() << " "
          << spl->z() << std::endl;
    }
  }
  
  lol_it = RepositoryCanvasResult::canvas()->getBackSeams()->begin();
  lol_it_end = RepositoryCanvasResult::canvas()->getBackSeams()->end();

  seamNumber=0;
  for(; lol_it != lol_it_end; ++lol_it) {
    std::list<Point>::iterator spl = lol_it->begin();
    std::list<Point>::iterator spl_end = lol_it->end();
    
    _obj_file<<"g bseam"<<++seamNumber << std::endl;
    for(; spl != spl_end; ++spl) {
      _obj_file << "v " // OBJ file format for vertices
          << spl->x() << " "
          << spl->y() << " "
          << spl->z() << std::endl;
    }
  }
  
  
  _obj_file.close();
  
  // Save strokes to file
  QString stroke_file_name = filename;
  stroke_file_name.replace(".3ds", ".strokes");
  if(stroke_file_name == filename) { // there was no .3ds extension
    stroke_file_name += QString(".strokes");
  }
  
  std::ofstream stroke_file;
  stroke_file.open(stroke_file_name.latin1());
  RepositoryCanvasResult::canvas()->saveLayers(stroke_file);
  stroke_file.close();

  // Save triangles and vertices to disk.
  char str[] = "garment";
  Lib3dsFile *file = lib3ds_file_new();
  Lib3dsMesh *mesh = lib3ds_mesh_new(str);
  Lib3dsNode *node = lib3ds_node_new_object();

  strcpy(node->name, str);
  node->parent_id = LIB3DS_NO_PARENT;

  Lib3dsLin3Key *key_pos = new Lib3dsLin3Key();
  lib3ds_vector_zero(key_pos->value);
  lib3ds_lin3_track_insert(&node->data.object.pos_track, key_pos);

  Lib3dsQuatKey *key_rot = new Lib3dsQuatKey();
  lib3ds_quat_identity(key_rot->q);
  lib3ds_quat_track_insert(&node->data.object.rot_track, key_rot);

  Lib3dsLin3Key *key_scl = new Lib3dsLin3Key();
  key_scl->value[0] = 1.0;
  key_scl->value[1] = 1.0;
  key_scl->value[2] = 1.0;
  lib3ds_lin3_track_insert(&node->data.object.scl_track, key_scl);

  lib3ds_file_insert_node(file, node);
  lib3ds_file_insert_mesh(file, mesh);
  lib3ds_mesh_new_point_list(mesh, vertices.size());
  lib3ds_mesh_new_face_list(mesh, triangles.size());

  for (unsigned i = 0; i < mesh->points; ++i) {
    mesh->pointL[i].pos[0] = vertices[i].x();
    mesh->pointL[i].pos[1] = vertices[i].y();
    mesh->pointL[i].pos[2] = vertices[i].z();
  }
  
  for (unsigned i = 0; i < mesh->faces; ++i) {
    mesh->faceL[i].points[0] = triangles[i].x();
    mesh->faceL[i].points[1] = triangles[i].y();
    mesh->faceL[i].points[2] = triangles[i].z();
  }

  triangles.clear();
  vertices.clear();

  if (!lib3ds_file_save(file, filename)) {
    emit statusBarMessage("Unable to save file \"" + filename + "\".", 2000);
    emit message(MsgHandler::MSG_ERROR, name(), "Unable to save file \"" + filename + "\".");
  } else {
    emit statusBarMessage("File \"" + filename + "\" successfully saved.", 2000);
    emit message(MsgHandler::MSG_NORMAL, name(), "File \"" +
		 filename + "\" successfully saved.");
    if (_mdi && _mdi->config())
      _mdi->config()->io()->setValue("paths/file_3ds", filename);
  }
  lib3ds_file_free(file);
  QApplication::restoreOverrideCursor();
}
