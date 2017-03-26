//
//  Filename         : canvas.cpp
//  Author           : Emmanuel Turquin
//  Purpose          : A 2D drawing canvas.
//  Date of creation : 04/23/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include <qimage.h>
#include <qbitmap.h>
#include <qcursor.h>
#include <qfiledialog.h>
#include "config.h"
#include "config_app_info.h"
#include "config_canvas.h"
#include "config_pretreatment.h"
#include "utils.h"
#include "glutils.h"
#include "mdi_window.h"
#include "layer.h"
#include "repository_pretreatment_canvas.h"
#include "repository_canvas_result.h"
#include "pretreatment.h"
#include "point_utils.h"
#include "garment_map.h"
#include "canvas.h"
#include "point.h"
#include <string>
#include <istream>
#include <fstream>

#include "../data/pixmaps/paintbrush_small.xbm"
#include "../data/pixmaps/paintbrush_small_mask.xbm"
#include "../data/textures/paper/whitepaper.xpm"

using namespace qglviewer;


class CanvasConstraint : public WorldConstraint
{
public:

  CanvasConstraint(double size = 1.0) : WorldConstraint() {
    _size = size;
  }
  ~CanvasConstraint() {}

  void setSize(double size) {
    _size = size;
  }

  virtual void constrainTranslation(Vec& t, Frame * const fr)
  {
    if (fr->position().z + t.z < 0.1 * _size)
      t.z = 0.1 * _size - fr->position().z;
    else if (fr->position().z + t.z > 10.0 * _size)
      t.z = 10 * _size - fr->position().z;
    if (fr->position().x + t.x < -1.0 * _size)
      t.x = -1.0 * _size - fr->position().x;
    else if (fr->position().x + t.x > 1.0 * _size)
      t.x = 1.0 * _size - fr->position().x;
    if (fr->position().y + t.y < -1.0 * _size)
      t.y = -1.0 * _size - fr->position().y;
    else if (fr->position().y + t.y > 1.0 * _size)
      t.y = 1.0 * _size - fr->position().y;
  }

private:

  double _size;
};

Canvas::Canvas(MDIWindow* parent, const char* name,
	       const QGLWidget* shareWidget, int wflags)
  : BaseQGLViewer(parent, name, shareWidget, wflags)
{
  setCursor(QCursor(QBitmap(paintbrush_small_width, paintbrush_small_height, paintbrush_small_bits, true),
		    QBitmap(paintbrush_small_mask_width, paintbrush_small_mask_height, paintbrush_small_mask_bits, true),
		    4, 3));
  _mdi = parent;
  _front_texture_dl = 0;
  _back_texture_dl = 0;
  _bb_dl = 0;
  _front = true;
  _two_sides = true;
  _seam_mode = false;
  _gfold_mode = false;
  _vertical_mirror = Canvas::READY;
  initLayers();
  initGarmentMaps();
  updateConfig();
  RepositoryPretreatmentCanvas::setCanvas(this);
  RepositoryCanvasResult::setCanvas(this);
}

Canvas::~Canvas()
{
  freeLayers();
  freeGarmentMaps();

  if (RepositoryPretreatmentCanvas::canvas() == this)
    RepositoryPretreatmentCanvas::setCanvas(0);
  if (RepositoryCanvasResult::canvas() == this)
    RepositoryCanvasResult::setCanvas(0);
}

void Canvas::canvasMessage(const MsgHandler::MsgType& type, const QString& msg) {
  emit message(type, name(), msg);
}

void Canvas::updateConfig()
{
  makeCurrent();

  int sampling = config::CANVAS_SAMPLING_INT;
  int snapping = config::CANVAS_SNAPPING_INT;
  float threshold = config::CANVAS_SPLITTING_THRESHOLD;
  _chain_pts = config::CANVAS_CHAIN_PTS;
  _stroke_pts = config::CANVAS_STROKE_PTS;
  _texture = config::CANVAS_TEXTURE;
  _type = config::CANVAS_TYPE;
  _bb = config::CANVAS_BBOX;
  _gm_size_x = config::CANVAS_GM_X;
  _gm_size_y =  config::CANVAS_GM_Y;
  double border_factor_x = config::CANVAS_BORDER_X;
  double border_factor_y = config::CANVAS_BORDER_Y;
  bool antialiasing = config::CANVAS_ANTIALIASING;
  bool compute = false;
  if (_mdi && _mdi->config()) {
    _mdi->config()->io()->getValue("intervals/sampling", sampling);
    _mdi->config()->io()->getValue("intervals/snapping", snapping);
    _mdi->config()->io()->getValue("gm/size/x", _gm_size_x);
    _mdi->config()->io()->getValue("gm/size/y", _gm_size_y);
    _mdi->config()->io()->getValue("border/factor/x", border_factor_x);
    _mdi->config()->io()->getValue("border/factor/y", border_factor_y);
    _mdi->config()->io()->getValue("splitting/threshold", threshold);
    _mdi->config()->io()->getValue("display/chain_pts", _chain_pts);
    _mdi->config()->io()->getValue("display/stroke_pts", _stroke_pts);
    _mdi->config()->io()->getValue("display/texture", _texture);
    _mdi->config()->io()->getValue("display/type", _type);
    _mdi->config()->io()->getValue("display/bbox", _bb);
    _mdi->config()->io()->getValue("display/antialiasing", antialiasing);
    _mdi->config()->io()->getValue("compute/start", compute);
  }
  RepositoryCanvasResult::setBorderFactor(Vec2d(border_factor_x, border_factor_x));
  _front_layer->setSampling(sampling);
  _front_layer->setSnapping(snapping);
  _front_layer->setSplittingThreshold(threshold);
  _back_layer->setSampling(sampling);
  _back_layer->setSnapping(snapping);
  _back_layer->setSplittingThreshold(threshold);
  if (antialiasing) {
    glEnable(GL_POINT_SMOOTH); // FIXME: buggy
    glEnable(GL_LINE_SMOOTH); // FIXME: buggy
  }
  else {
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_LINE_SMOOTH);
  }
  if (compute) { // Start computation.
    QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
    _mdi->config()->io()->setValue("compute/start", false);
    computeMesh();
    RepositoryCanvasResult::canvasUpdated();
    QApplication::restoreOverrideCursor();
  }
  updateGL();
}

void Canvas::init()
{
  // Snippet of code copied from updateConfig(), needed to initialize GL antialiasing.
  bool antialiasing = config::CANVAS_ANTIALIASING;
  if (_mdi && _mdi->config())
    _mdi->config()->io()->getValue("display/antialiasing", antialiasing);
  if (antialiasing) {
    glEnable(GL_POINT_SMOOTH); // FIXME: buggy
    glEnable(GL_LINE_SMOOTH); // FIXME: buggy
  }
  else {
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_LINE_SMOOTH);
  }

  // GL parameters.
  glEnable(GL_BLEND);
#ifdef GL_FUNC_ADD    // Phil
  glBlendEquation(GL_FUNC_ADD);
#endif
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);
  glEnable(GL_MAP1_VERTEX_3);

  // Camera initialization.
  camera()->setType(Camera::ORTHOGRAPHIC);

  // Forbid rotation.
  CanvasConstraint* constraint = new CanvasConstraint(utils::max(RepositoryPretreatmentCanvas::bbMax()[0], RepositoryPretreatmentCanvas::bbMax()[1]));
  constraint->setRotationConstraintType(AxisPlaneConstraint::FORBIDDEN);
  camera()->frame()->setConstraint(constraint);

  // Enable mouse tracking.
  setMouseTracking(true);
  _left_button_pressed = false;
  _seam_mode = false;
  _gfold_mode = false;

  // Initialize the display list for the sheet of paper.
  preparePaperSheet();
}

void Canvas::initLayers()
{
  _front_layer = new Layer(this, RepositoryPretreatmentCanvas::bbMin(), RepositoryPretreatmentCanvas::bbMax());
  _back_layer = new Layer(this, RepositoryPretreatmentCanvas::bbMin(), RepositoryPretreatmentCanvas::bbMax());
  _front_seam_layer = new Layer(this, RepositoryPretreatmentCanvas::bbMin(), RepositoryPretreatmentCanvas::bbMax());
  _back_seam_layer = new Layer(this, RepositoryPretreatmentCanvas::bbMin(), RepositoryPretreatmentCanvas::bbMax());
  _front_gfold_layer= new Layer(this, RepositoryPretreatmentCanvas::bbMin(), RepositoryPretreatmentCanvas::bbMax());
  _front_gfold_layer->_ltype = Layer::GFOLD;
  _back_gfold_layer= new Layer(this, RepositoryPretreatmentCanvas::bbMin(), RepositoryPretreatmentCanvas::bbMax());
  _back_gfold_layer->_ltype = Layer::GFOLD;
  
  RepositoryCanvasResult::setFrontLayer(_front_layer);
  RepositoryCanvasResult::setBackLayer(_back_layer);
  RepositoryCanvasResult::setFrontSeamLayer(_front_seam_layer);
  RepositoryCanvasResult::setBackSeamLayer(_back_seam_layer);
}

void Canvas::freeLayers()
{
  delete _front_layer;
  delete _back_layer;
  delete _front_seam_layer;
  delete _back_seam_layer;
  delete _front_gfold_layer;
  delete _back_gfold_layer;
  RepositoryCanvasResult::setFrontLayer(0);
  RepositoryCanvasResult::setFrontSeamLayer(0);
  RepositoryCanvasResult::setBackSeamLayer(0);
  RepositoryCanvasResult::setBackLayer(0);
}

void Canvas::initGarmentMaps()
{
  RepositoryCanvasResult::setFrontGarmentMaps(&_front_garment_maps);
  RepositoryCanvasResult::setBackGarmentMaps(&_back_garment_maps);
}

void Canvas::freeGarmentMaps()
{
  // Front garment maps.
  GarmentMapList::iterator it = _front_garment_maps.begin();
  GarmentMapList::iterator it_end = _front_garment_maps.end();
  for ( ; it != it_end; ++it)
    delete *it;
  _front_garment_maps.clear();
  RepositoryCanvasResult::setFrontGarmentMaps(0);
  // Back garment maps.
  it = _back_garment_maps.begin();
  it_end = _back_garment_maps.end();
  for ( ; it != it_end; ++it)
    delete *it;
  _back_garment_maps.clear();
  RepositoryCanvasResult::setBackGarmentMaps(0);
}

void Canvas::preparePaperSheet()
{
  _paper_dl = glGenLists(1);
  if (!_paper_dl)
    return;

  QImage paper = convertToGLFormat(QImage(whitepaper_xpm));

  unsigned tex_id;
  glGenTextures(1, &tex_id);
  glBindTexture(GL_TEXTURE_2D, tex_id);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, paper.width(),
	       paper.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE,
	       paper.bits());

  glNewList(_paper_dl, GL_COMPILE);

  double size = 100.0;
  double max = utils::max(RepositoryPretreatmentCanvas::bbMax()[0], RepositoryPretreatmentCanvas::bbMax()[1]);
  if (max)
    size *= max;

  glEnable(GL_TEXTURE_2D);
  glColor3f(0.95, 0.95, 0.95);
  glBindTexture(GL_TEXTURE_2D, tex_id);
  glBegin(GL_TRIANGLE_STRIP);
  glTexCoord2f(-100.0, -100.0); // 1
  glVertex2f(-size, size);
  glTexCoord2f(100.0, -100.0); // 2
  glVertex2f(-size, -size);
  glTexCoord2f(-100.0, 100.0); // 3
  glVertex2f(size, size);
  glTexCoord2f(100.0, 100.0); // 4
  glVertex2f(size, -size);
  glEnd();
  glDisable(GL_TEXTURE_2D);

  glEndList();

  dynamic_cast<CanvasConstraint*>(camera()->frame()->constraint())->setSize(size / 100.0);
  //camera()->fitSphere(0.0, 0.0, 0.0, size / 100.0);
  // JDW method interface has changed. Needs a Vec and a radius
  camera()->fitSphere(qglviewer::Vec(0.0, 0.0, 0.0), size / 100.0);
}

void Canvas::freePaperSheet()
{
  if (_paper_dl)
    glDeleteLists(_paper_dl, 1);
}

void Canvas::keyPressEvent(QKeyEvent* e)
{
  switch (e->key())
    {
    case Qt::Key_M :
      _two_sides = !_two_sides;
      if (_two_sides)
	message(MsgHandler::MSG_NORMAL, name(), "Drawing mode set to \"two sides\".");
      else
	message(MsgHandler::MSG_NORMAL, name(), "Drawing mode set to \"one side\".");
      updateGL();
      break;
    case Qt::Key_I :
      
      if (_vertical_mirror == Canvas::OFF) {
        _vertical_mirror = Canvas::READY;
        message(MsgHandler::MSG_NORMAL, name(), "Vertical mirror set to \"On\".");
      }
      else
      {
        _vertical_mirror = Canvas::OFF;
        message(MsgHandler::MSG_NORMAL, name(), "Vertical mirror set to \"Off\".");
      }
      updateGL();
      break;
    case Qt::Key_V :
      _front = !_front;
      if (_front)
	message(MsgHandler::MSG_NORMAL, name(), "View set to \"front\".");
      else
	message(MsgHandler::MSG_NORMAL, name(), "View set to \"back\".");
      updateGL();
      break;
    case Qt::Key_S :
	 _seam_mode = !_seam_mode;
      _gfold_mode = false; // mutally exclusive
	if(_seam_mode)
		message(MsgHandler::MSG_NORMAL, name(), "Seam mode set to \"on\".");
	else
		message(MsgHandler::MSG_NORMAL, name(), "Seam mode set to \"off\".");
	break;
      case Qt::Key_F :
        _gfold_mode = !_gfold_mode;
        _seam_mode = false;
        if(_gfold_mode)
          message(MsgHandler::MSG_NORMAL, name(), "Gauss fold mode set to \"on\".");
        else
          message(MsgHandler::MSG_NORMAL, name(), "Gauss fold mode set to \"off\".");
        break;
    default:
      QGLViewer::keyPressEvent(e);
    }
}

void Canvas::convertCoordinates(double x, double y, double& new_x, double& new_y, double& ratio)
{
  float camera_pos[3];
  double w, h;
  //JDW - method depreciated in later versions of QGLViewer. FIXME needs more tidying up. Should use Vec everywhere
  //camera()->getPosition(camera_pos[0], camera_pos[1], camera_pos[2]);
  qglviewer::Vec pos=qglviewer::Vec();
  pos = camera()->position();
  camera_pos[0] = pos[0];
  camera_pos[1] = pos[1];
  camera_pos[2] = pos[2];


  camera()->getOrthoWidthHeight(w, h);
  ratio = 2.0 * w / width();
  new_x = camera_pos[0] + w * (2.0 * x / width() - 1);
  new_y = camera_pos[1] - h * (2.0 * y / height() - 1);

  if (_left_button_pressed) {
    if (new_x > RepositoryPretreatmentCanvas::bbMax()[0])
      new_x = RepositoryPretreatmentCanvas::bbMax()[0];
    else if (new_x < -RepositoryPretreatmentCanvas::bbMax()[0])
      new_x = -RepositoryPretreatmentCanvas::bbMax()[0];
    if (new_y > RepositoryPretreatmentCanvas::bbMax()[1])
      new_y = RepositoryPretreatmentCanvas::bbMax()[1];
    else if (new_y < -RepositoryPretreatmentCanvas::bbMax()[1])
      new_y = -RepositoryPretreatmentCanvas::bbMax()[1];
  }
}

void Canvas::mouseMoveEvent(QMouseEvent *e)
{
  double coord[2];
  double ratio;
  convertCoordinates(e->x(), e->y(), coord[0], coord[1], ratio);

  emit statusBarMessage("Mouse coordinates: [" +
			QString::number(coord[0]) + "; " +
			QString::number(coord[1]) + "]", 1000);

  if (_left_button_pressed) {
    if(_seam_mode) {
      // create seams
      if (_front)
      {
        _front_seam_layer->intermediatePoint(coord[0], coord[1], ratio);
        if (_two_sides)
          _back_seam_layer->intermediatePoint(coord[0], coord[1], ratio);
      }
      else
      {
        _back_seam_layer->intermediatePoint(-coord[0], coord[1], ratio);
        if (_two_sides)
          _front_seam_layer->intermediatePoint(-coord[0], coord[1], ratio);
      }
    } 
    else if (_gfold_mode) {
      if(_front) {
        _front_gfold_layer->intermediatePoint(coord[0], coord[1], ratio);
      }
      else {
        _back_gfold_layer->intermediatePoint(-coord[0], coord[1], ratio);
      }
    }
    else {
      // normal sketching
      if (_front)
      {
        // disallow mirror strokes on other side of centre line
        if(_vertical_mirror != OFF) {
          // terminate stroke when vertical line is crossed
          if(_vertical_mirror == RIGHT && coord[0] < 0.0) {
            Canvas::mouseReleaseEvent(e);
            return;
          }
          if(_vertical_mirror == LEFT && coord[0] >= 0.0) {
            Canvas::mouseReleaseEvent(e);
            return;
          }
        }
        
        _front_layer->intermediatePoint(coord[0], coord[1], ratio);
        if (_two_sides)
          _back_layer->intermediatePoint(coord[0], coord[1], ratio);
      }
      else
      {
                // disallow mirror strokes on other side of centre line
        if(_vertical_mirror != OFF) {
          // terminate stroke when vertical line is crossed
          if(_vertical_mirror == RIGHT && coord[0] < 0.0) {
            Canvas::mouseReleaseEvent(e);
            return;
          }
          if(_vertical_mirror == LEFT && coord[0] >= 0.0) {
            Canvas::mouseReleaseEvent(e);
            return;
          }
        }
        
        _back_layer->intermediatePoint(-coord[0], coord[1], ratio);
        if (_two_sides)
          _front_layer->intermediatePoint(-coord[0], coord[1], ratio);
      }
    }
    updateGL();
    return;
  }

  BaseQGLViewer::mouseMoveEvent(e);
}


void Canvas::saveLayers(std::ostream& out) {
  out << *_front_layer << std::endl;
  out << *_back_layer << std::endl;
  out << *_front_seam_layer << std::endl;
  out << *_back_seam_layer << std::endl;
  out << *_front_gfold_layer << std::endl;
  out << *_back_gfold_layer << std::endl;
}

void Canvas::loadLayers(std::istream& in) {
  _front_layer->setSnapping(1);
  _front_layer->setSampling(1);
  readToLayer(in, _front_layer);
  _front_layer->setSnapping(10);
  _front_layer->setSampling(5);
  
  _back_layer->setSnapping(1);
  _back_layer->setSampling(1);
  readToLayer(in, _back_layer);
  _back_layer->setSnapping(10);
  _back_layer->setSampling(5);
  
  _front_seam_layer->setSnapping(1);
  _front_seam_layer->setSampling(1);
  readToLayer(in, _front_seam_layer);
  _front_seam_layer->setSnapping(10);
  _front_seam_layer->setSampling(5);
  
  
  
  _back_seam_layer->setSnapping(1);
  _back_seam_layer->setSampling(1);
  readToLayer(in, _back_seam_layer);
  _back_seam_layer->setSnapping(10);
  _back_seam_layer->setSampling(5);
  
  
  _front_gfold_layer->setSnapping(1);
  _front_gfold_layer->setSampling(1);
  readToLayer(in, _front_gfold_layer);
  _front_gfold_layer->setSnapping(10);
  _front_gfold_layer->setSampling(5);
  
  
  _back_gfold_layer->setSnapping(1);
  _back_gfold_layer->setSampling(1);
  readToLayer(in, _back_gfold_layer);
  _back_gfold_layer->setSnapping(10);
  _back_gfold_layer->setSampling(5);

}


void Canvas::loadFile()
{
  if (!RepositoryPretreatmentCanvas::frontTexture()) {
    message(MsgHandler::MSG_ERROR, name(), "Cannot load strokes until model has been loaded into pretreatment");
    return;
  }
  
  QString path("/tmp");
  if (_mdi && _mdi->config())
    _mdi->config()->io()->getValue("paths/file_strokes", path);
  QString filename = QFileDialog::getOpenFileName(path, "stroke files (*.strokes);;All files (*)", this);
  
  // In case of Cancel.
  if (filename.isEmpty())
    return;

  QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));

  

  // Load a strokes definition file.
  std::ifstream _file(filename.latin1());
  if (!_file)
  {
    emit statusBarMessage("Unable to open file \"" + filename + "\".", 2000);
    emit message(MsgHandler::MSG_ERROR, name(), "Unable to open file \"" + filename + "\".");
    return;
  }

  if (_mdi && _mdi->config())
    _mdi->config()->io()->setValue("paths/file_strokes", filename);
    
  loadLayers( _file );
  emit statusBarMessage("Loaded stroke file file \"" + filename + "\".", 2000);
  _file.close();
  
  QApplication::restoreOverrideCursor();
}


void Canvas::readToLayer(std::istream& in, Layer* layer) {
  layer->clear();
  
  double w, h;
  
  std::string s;
  Point p;
  
  camera()->getOrthoWidthHeight(w, h);
  double ratio = 2.0 * w / width();
  
  in >> s;
  if(s == "Layer{")
  {
    int numOfChainsLeft;
    in >> numOfChainsLeft;
    
    while(numOfChainsLeft>0) 
    {
      in >> s;
      if(s == "Chain{")
      {
        int numStrokesLeft;
        in >> numStrokesLeft;
        
        while(numStrokesLeft>0 && !in.eof()) 
        { // loop over strokes
                
          in >> s; // Assume this is a "Stroke{" token
          int numPointsLeft;
          in >> numPointsLeft;
          
          // nasty hack to support gaussian folds. We need more classes when we have the time
          Vec2d paramsA, paramsB;
          
          in >> paramsA[0];
          in >> paramsA[1];
          in >> paramsB[0];
          in >> paramsB[1];
          
          // assume at least 2 points (one Segment)
          in >> p; numPointsLeft--;
          
          layer->firstPoint(p.x(),p.y(),ratio);
          
          layer->currentStroke()->setAParams( paramsA );
          layer->currentStroke()->setBParams( paramsB );
          
          while(numPointsLeft > 1 && !in.eof())
          {
            in >> p; numPointsLeft--;
            layer->intermediatePoint(p.x(),p.y(),ratio);
          }
          in >> p; numPointsLeft--;
          layer->lastPoint(p.x(), p.y(), ratio);
          
          in >> s; // Assume closing stroke token "}"
          numStrokesLeft--;
        }
        in >> s; // Assume closing chain token "}"
        numOfChainsLeft--;
      }
    }
  }
  in >> s; // Assume closing layer token "}"
}


void Canvas::mousePressEvent(QMouseEvent *e)
{
  
  if ((e->button() == Qt::LeftButton) && (e->state() == Qt::NoButton)) {
    if (!RepositoryPretreatmentCanvas::frontTexture()) {
      message(MsgHandler::MSG_ERROR, name(), "A texture is required to draw onto.");
      return;
    }
    _left_button_pressed = true;

    double coord[2];
    double ratio;
    convertCoordinates(e->x(), e->y(), coord[0], coord[1], ratio);
    
    if(_vertical_mirror != Canvas::OFF) {
    if(coord[0] >0.0)
      _vertical_mirror = Canvas::RIGHT;
    else
      _vertical_mirror = Canvas::LEFT;
    }
    
	// sketch seams
    if(_seam_mode) {
      if (_front) {
        _front_seam_layer->firstPoint(coord[0], coord[1], ratio);
        if (_two_sides)
          _back_seam_layer->firstPoint(coord[0], coord[1], ratio);
      }
      else {
        _back_seam_layer->firstPoint(-coord[0], coord[1], ratio);
        if (_two_sides)
          _front_seam_layer->firstPoint(-coord[0], coord[1], ratio);
      }
		
    }
     // gaussian folds
    else if (_gfold_mode){
      if(_front) {
        _front_gfold_layer->firstPoint(coord[0], coord[1], ratio);
      }
      else {
        _back_gfold_layer->firstPoint(-coord[0], coord[1], ratio);
      }
    }
    // normal sketching
    else {
		if (_front) {
			_front_layer->firstPoint(coord[0], coord[1], ratio);
			if (_two_sides)
				_back_layer->firstPoint(coord[0], coord[1], ratio);
		}
		else {
			_back_layer->firstPoint(-coord[0], coord[1], ratio);
			if (_two_sides)
				_front_layer->firstPoint(-coord[0], coord[1], ratio);
		}
    }

    updateGL();
  }
  else
    BaseQGLViewer::mousePressEvent(e);
}

void Canvas::mouseReleaseEvent(QMouseEvent *e)
{
  if (!RepositoryPretreatmentCanvas::frontTexture()) {
    message(MsgHandler::MSG_ERROR, name(), "A texture is required to draw onto.");
    return;
  }
  
  if ((e->button() == Qt::LeftButton) && _left_button_pressed) {
    double coord[2];
    double ratio;
    convertCoordinates(e->x(), e->y(), coord[0], coord[1], ratio);
    
    // Save points to create mirrored strokes if required
    std::list< counted_ptr<Point> > pList;

    
    if(_seam_mode) {
      if (_front) {
        _front_seam_layer->lastPoint(coord[0], coord[1], ratio);
        if (_two_sides)
          _back_seam_layer->lastPoint(coord[0], coord[1], ratio);
      }
      else {
        _back_seam_layer->lastPoint(-coord[0], coord[1], ratio);
        if (_two_sides)
          _front_seam_layer->lastPoint(-coord[0], coord[1], ratio);
      }
    }
    else if (_gfold_mode) {
      if(_front) {
        _front_gfold_layer->lastPoint(coord[0], coord[1], ratio);
      }
      else {
        _back_gfold_layer->lastPoint(-coord[0], coord[1], ratio);
      }
    }
    // normal sketching
    else {
      
      if(_vertical_mirror != Canvas::OFF) {
        
        std::list< Segment* >::const_iterator it, it_end;
      
        if(_front) {
          it = _front_layer->currentStroke()->segments()->begin();
          it_end = _front_layer->currentStroke()->segments()->end();
        } else {
          it = _back_layer->currentStroke()->segments()->begin();
          it_end = _back_layer->currentStroke()->segments()->end();
        }
        if(it!=it_end) {
          pList.push_back( (*it)->pointA() );
        }
        for(;it!=it_end;it++) {
          pList.push_back( (*it)->pointB() );
        }
        
      }
      
      if (_front) {
        
        // disallow mirror strokes on other side of centre line
        if(_vertical_mirror != OFF) {
          if(_vertical_mirror == RIGHT && coord[0] < 0.0) {
            coord[0] = 0.0;
          }
          if(_vertical_mirror == LEFT && coord[0] >= 0.0) {
            coord[0] = 0.0;
          }
        }
        
        _front_layer->lastPoint(coord[0], coord[1], ratio);
        
        // Create the mirror stroke also
        if(_vertical_mirror != Canvas::OFF) {
          std::list< counted_ptr< Point > >::const_iterator it, it_end;
          
          if(pList.size() >= 2) 
          {
            it = pList.begin(); it_end=pList.end();
            
            _front_layer->firstPoint( -(*it)->x(), (*it)->y(), ratio );
            
            it++;
            for(;it!=it_end;it++) {
              _front_layer->intermediatePoint( -(*it)->x(), (*it)->y(), ratio );
            }
            _front_layer->lastPoint( -coord[0], coord[1], ratio );
          }
        }
        
        if (_two_sides)
        {
          _back_layer->lastPoint(coord[0], coord[1], ratio);
          
                  
        // Create the mirror stroke also
          if(_vertical_mirror != Canvas::OFF) {
            if(pList.size() >= 2) 
            {
              std::list< counted_ptr< Point > >::const_iterator it, it_end;
              it = pList.begin(); it_end=pList.end();
              _back_layer->firstPoint( -(*it)->x(), (*it)->y(), ratio );
              it++;
              for(;it!=it_end;it++) {
                _back_layer->intermediatePoint( -(*it)->x(), (*it)->y(), ratio );
              }
              _back_layer->lastPoint( -coord[0], coord[1], ratio );
            }
          }
          
        }
      }
      else { // viewing from back
        // disallow mirror strokes on other side of centre line
        if(_vertical_mirror != OFF) {
          if(_vertical_mirror == RIGHT && coord[0] < 0.0) {
            coord[0] = 0.0;
          }
          if(_vertical_mirror == LEFT && coord[0] >= 0.0) {
            coord[0] = 0.0;
          }
        }
        
        _back_layer->lastPoint(-coord[0], coord[1], ratio);
        
        // Create the mirror stroke also
        if(_vertical_mirror != Canvas::OFF) {
          if(pList.size() >= 2) 
          {
            std::list< counted_ptr< Point > >::const_iterator it, it_end;
            it = pList.begin(); it_end=pList.end();
            _back_layer->firstPoint( -(*it)->x(), (*it)->y(), ratio );
            it++;
            for(;it!=it_end;it++) {
              _back_layer->intermediatePoint( -(*it)->x(), (*it)->y(), ratio );
            }
            _back_layer->lastPoint( coord[0], coord[1], ratio );
          }
        }
        
        if (_two_sides) {
            _front_layer->lastPoint(-coord[0], coord[1], ratio);
            
                    // Create the mirror stroke also
            if(_vertical_mirror != Canvas::OFF) 
            {
              if(pList.size() >= 2) 
              {
                std::list< counted_ptr< Point > >::const_iterator it, it_end;
                it = pList.begin(); it_end=pList.end();
                _front_layer->firstPoint( -(*it)->x(), (*it)->y(), ratio );
                it++;
                for(;it!=it_end;it++) {
                  _front_layer->intermediatePoint( -(*it)->x(), (*it)->y(), ratio );
                }
              }
              _front_layer->lastPoint( coord[0], coord[1], ratio );
            }
            
        }
      }
    }
    updateGL();
    _left_button_pressed = false;
  }
  else
    BaseQGLViewer::mouseReleaseEvent(e);
}

// void Canvas::copyLayer()
// {
//   Layer *current_front_layer;
//   Layer *current_back_layer;
//   Layer::ChainList::const_iterator ch;
//   Layer::ChainList::const_iterator ch_end;
//   Chain::StrokeList::const_iterator st;
//   Chain::StrokeList::const_iterator st_end;
//   Stroke::SegmentList::const_iterator seg;
//   Stroke::SegmentList::const_iterator seg_end;

//   if (_front) {
//     current_back_layer = _back_layer;
//     current_front_layer = _front_layer;
//   }
//   else {
//     current_back_layer = _front_layer;
//     current_front_layer = _back_layer;
//   }

//   double ratio, tmp;
//   camera()->getOrthoWidthHeight(ratio, tmp);
//   ratio = 2.0 * ratio / width();

//   current_front_layer->clear();

//   ch = current_back_layer->chains()->begin();
//   ch_end = current_back_layer->chains()->end();
//   for (; ch != ch_end; ++ch) {
//     st = (*ch)->strokes()->begin();
//     st_end = (*ch)->strokes()->end();
//     for (; st != st_end; ++st) {
//       current_front_layer->firstPoint((*st)->pointA()->x(), (*st)->pointA()->y(), ratio);
//       seg = (*st)->segments()->begin();
//       seg_end = (*st)->segments()->end();
//       for (++seg; seg != seg_end; ++seg)
// 	current_front_layer->intermediatePoint((*seg)->pointA()->x(), (*seg)->pointA()->y(), ratio);
//       current_front_layer->lastPoint((*st)->pointB()->x(), (*st)->pointB()->y(), ratio);
//     }
//   }
// }

void Canvas::draw()
{
  if (_paper_dl)
    glCallList(_paper_dl);

  if (_bb && _bb_dl)
    glCallList(_bb_dl);

  if (_texture) {
    if (_front && _front_texture_dl)
      glCallList(_front_texture_dl);
    else if (!_front && _back_texture_dl)
      glCallList(_back_texture_dl);
  }

  float camera_pos[3];
  //JDW depreciated
  //camera()->getPosition(camera_pos[0], camera_pos[1], camera_pos[2]);
  qglviewer::Vec pos = qglviewer::Vec();
  camera()->position();
  camera_pos[0] = pos[0];
  camera_pos[1] = pos[1];
  camera_pos[2] = pos[2];

  double w, h;
  //JDW not sure why it was here twice
  //camera()->getPosition(camera_pos[0], camera_pos[1], camera_pos[2]);
  camera()->getOrthoWidthHeight(w, h);
//  std::cout << "Ortho w: " << w << std::endl;
  double line_width = utils::min(0.5 * utils::min(RepositoryPretreatmentCanvas::bbMax()[0], RepositoryPretreatmentCanvas::bbMax() [1]) / w, 5.0); 
//  std::cout << "Line width: " << line_width << std::endl;

  Layer *current_front_layer;
  Layer *current_back_layer;
  
  Layer *current_front_seam_layer;
  Layer *current_back_seam_layer;
  
  Layer *current_front_gfold_layer;
  Layer *current_back_gfold_layer;
  
  Layer::ChainList::const_iterator ch;
  Layer::ChainList::const_iterator ch_end;
  Chain::StrokeList::const_iterator st;
  Chain::StrokeList::const_iterator st_end;
  Stroke::SegmentList::const_iterator seg;
  Stroke::SegmentList::const_iterator seg_end;

  float *ctrl_pts;
  int pts_nb;
  int i, j;

  if (_front) {
    current_back_layer = _back_layer;
    current_front_layer = _front_layer;
    
    current_front_seam_layer = _front_seam_layer;
    current_back_seam_layer = _back_seam_layer;
    
    current_front_gfold_layer = _front_gfold_layer;
    current_back_gfold_layer = _back_gfold_layer;
  }
  else {
    current_back_layer = _front_layer;
    current_front_layer = _back_layer;
    
    current_front_seam_layer = _back_seam_layer;
    current_back_seam_layer = _front_seam_layer;
    current_front_gfold_layer = _back_gfold_layer;
    current_back_gfold_layer = _front_gfold_layer;
  }

  glEnable(GL_DEPTH_TEST);

  // Current front layer.

  glPolygonOffset(-1.0, -1.0);
  glEnable(GL_POLYGON_OFFSET_LINE);

  glColor3f(0.3, 0.3, 0.3);
  glLineWidth(line_width);

  /////////////// CUBIC SPLINES ///////////////

  ch = current_front_layer->chains()->begin();
  ch_end = current_front_layer->chains()->end();

  for (; ch != ch_end; ++ch) {
    st = (*ch)->strokes()->begin();
    st_end = (*ch)->strokes()->end();
    for (; st != st_end; ++st) {

      if (_type) {
	switch ((*st)->type()) {
	case Stroke::IN:
	  glColor3f(0.3, 0.9, 0.3);
	  break;
	case Stroke::OUT:
	  glColor3f(0.9, 0.3, 0.3);
	  break;
	case Stroke::MIXED:
	  glColor3f(0.9, 0.9, 0.3);
	  break;
	case Stroke::UNKNOWN:
	default:
	  glColor3f(0.3, 0.3, 0.3);
	}
      }

      pts_nb = (*st)->segments()->size() + 1;
      ctrl_pts = new float[3 * pts_nb];
      seg = (*st)->segments()->begin();
      seg_end = (*st)->segments()->end();
      if (seg != seg_end) {
	ctrl_pts[0] = mirror((*seg)->pointA()->x());
	ctrl_pts[1] = (*seg)->pointA()->y();
	ctrl_pts[2] = 0;
      }
      for (i = 3; seg != seg_end; ++seg, i += 3) {
	ctrl_pts[i] = mirror((*seg)->pointB()->x());
	ctrl_pts[i + 1] = (*seg)->pointB()->y();
	ctrl_pts[i + 2] = 0;
      }
      for (i = 0; i < pts_nb - 3; i += 3) {
	glMap1f(GL_MAP1_VERTEX_3, 0, 1, 3, 4, &ctrl_pts[3*i]);
	glBegin(GL_LINE_STRIP);
	for (j = 0; j <= 8; ++j)
	  glEvalCoord1f(float(j) / 8);
	glEnd();
      }
      glMap1f(GL_MAP1_VERTEX_3, 0, 1, 3, pts_nb - i, &ctrl_pts[3*i]);
      glBegin(GL_LINE_STRIP);
      for (j = 0; j <= 8; ++j)
	glEvalCoord1f(float(j) / 8);
      glEnd();
      delete[] ctrl_pts;
    }
  }

  /////////////// POLYLINES ///////////////

//   ch = current_front_layer->chains()->begin();
//   ch_end = current_front_layer->chains()->end();
//   for (; ch != ch_end; ++ch) {
//     st = (*ch)->strokes()->begin();
//     st_end = (*ch)->strokes()->end();
//     for (; st != st_end; ++st) {
//       seg = (*st)->segments()->begin();
//       seg_end = (*st)->segments()->end();
//       glBegin(GL_LINE_STRIP);
//       if (seg != seg_end)
// 	glVertex2d(mirror((*seg)->pointA()->x()), (*seg)->pointA()->y());
//       for (; seg != seg_end; ++seg)
// 	glVertex2d(mirror((*seg)->pointB()->x()), (*seg)->pointB()->y());
//       glEnd();
//     }
//   }

  glDisable(GL_POLYGON_OFFSET_LINE);

  // Current back layer.

  glColor3f(0.3, 0.3, 0.3);
  glLineWidth(line_width / 2);

  glEnable(GL_LINE_STIPPLE);
  glLineStipple (1, 0xF0F0); // FIXME

  ch = current_back_layer->chains()->begin();
  ch_end = current_back_layer->chains()->end();

  for (; ch != ch_end; ++ch) {
    st = (*ch)->strokes()->begin();
    st_end = (*ch)->strokes()->end();
    for (; st != st_end; ++st) {

      if (_type) {
	switch ((*st)->type()) {
	case Stroke::IN:
	  glColor3f(0.3, 0.9, 0.3);
	  break;
	case Stroke::OUT:
	  glColor3f(0.9, 0.3, 0.3);
	  break;
	case Stroke::MIXED:
	  glColor3f(0.9, 0.9, 0.3);
	  break;
	case Stroke::UNKNOWN:
	default:
	  glColor3f(0.3, 0.3, 0.3);
	}
      }

      pts_nb = (*st)->segments()->size() + 1;
      ctrl_pts = new float[3 * pts_nb];
      seg = (*st)->segments()->begin();
      seg_end = (*st)->segments()->end();
      if (seg != seg_end) {
	ctrl_pts[0] = mirror((*seg)->pointA()->x());
	ctrl_pts[1] = (*seg)->pointA()->y();
	ctrl_pts[2] = 0;
      }
      for (i = 3; seg != seg_end; ++seg, i += 3) {
	ctrl_pts[i] = mirror((*seg)->pointB()->x());
	ctrl_pts[i + 1] = (*seg)->pointB()->y();
	ctrl_pts[i + 2] = 0;
      }
      for (i = 0; i < pts_nb - 3; i += 3) {
	glMap1f(GL_MAP1_VERTEX_3, 0, 1, 3, 4, &ctrl_pts[3*i]);
	glBegin(GL_LINE_STRIP);
	for (j = 0; j <= 8; ++j)
	  glEvalCoord1f(float(j) / 8);
	glEnd();
      }
      glMap1f(GL_MAP1_VERTEX_3, 0, 1, 3, pts_nb - i, &ctrl_pts[3*i]);
      glBegin(GL_LINE_STRIP);
      for (j = 0; j <= 8; ++j)
	glEvalCoord1f(float(j) / 8);
      glEnd();
      delete[] ctrl_pts;
    }
  }

  glDisable(GL_LINE_STIPPLE);
  
  // JDW Start drawing the seam layer 
  
  /////////////// CURRENT FRONT SEAM LAYER ///////////////
  
  //std::cout << "CurrFrontSeamLayer: " << current_front_seam_layer << std::endl;
  //std::cout << "CurrBackSeamLayer: " << current_back_seam_layer << std::endl;
  
  

  ch = current_front_seam_layer->chains()->begin();
  ch_end = current_front_seam_layer->chains()->end();

  for (; ch != ch_end; ++ch) {

    
    st = (*ch)->strokes()->begin();
    st_end = (*ch)->strokes()->end();
    for (; st != st_end; ++st) {

      // Make obvious red for development
      glColor3f(0.5, 0.0, 0.0);

      pts_nb = (*st)->segments()->size() + 1;
      ctrl_pts = new float[3 * pts_nb];
      seg = (*st)->segments()->begin();
      seg_end = (*st)->segments()->end();
      if (seg != seg_end) {
        ctrl_pts[0] = mirror((*seg)->pointA()->x());
        ctrl_pts[1] = (*seg)->pointA()->y();
        ctrl_pts[2] = 0;
      }
      for (i = 3; seg != seg_end; ++seg, i += 3) {
        ctrl_pts[i] = mirror((*seg)->pointB()->x());
        ctrl_pts[i + 1] = (*seg)->pointB()->y();
        ctrl_pts[i + 2] = 0;
      }
      for (i = 0; i < pts_nb - 3; i += 3) {
        glMap1f(GL_MAP1_VERTEX_3, 0, 1, 3, 4, &ctrl_pts[3*i]);
        glBegin(GL_LINE_STRIP);
        for (j = 0; j <= 8; ++j)
          glEvalCoord1f(float(j) / 8);
        glEnd();
      }
      glMap1f(GL_MAP1_VERTEX_3, 0, 1, 3, pts_nb - i, &ctrl_pts[3*i]);
      glBegin(GL_LINE_STRIP);
      for (j = 0; j <= 8; ++j)
        glEvalCoord1f(float(j) / 8);
      glEnd();
      delete[] ctrl_pts;
    }
  }
  
  // JDW Show any seam currently being drawn
  glColor3f(0.8, 0.0, 0.0);
  if (current_front_seam_layer->currentStroke()) {
    seg = current_front_seam_layer->currentStroke()->segments()->begin();
    seg_end = current_front_seam_layer->currentStroke()->segments()->end();
    glBegin(GL_LINE_STRIP);
    if (seg != seg_end)
      glVertex2d(mirror((*seg)->pointA()->x()), (*seg)->pointA()->y());
    for (; seg != seg_end; ++seg)
      glVertex2d(mirror((*seg)->pointB()->x()), (*seg)->pointB()->y());
    glEnd();
  }
  // JDW Finish drawing the seam layer

  // GAUSSIAN FOLD DRAWING
  
  /////////////// CURRENT FRONT FOLD LAYER ///////////////
  
  ch = current_front_gfold_layer->chains()->begin();
  ch_end = current_front_gfold_layer->chains()->end();

  for (; ch != ch_end; ++ch) {

    
    st = (*ch)->strokes()->begin();
    st_end = (*ch)->strokes()->end();
    for (; st != st_end; ++st) {

      // Make obvious green for development
      glColor3f(0.0, 0.7, 0.0);

      pts_nb = (*st)->segments()->size() + 1;
      ctrl_pts = new float[3 * pts_nb];
      seg = (*st)->segments()->begin();
      seg_end = (*st)->segments()->end();
      if (seg != seg_end) {
        ctrl_pts[0] = mirror((*seg)->pointA()->x());
        ctrl_pts[1] = (*seg)->pointA()->y();
        ctrl_pts[2] = 0;
      }
      for (i = 3; seg != seg_end; ++seg, i += 3) {
        ctrl_pts[i] = mirror((*seg)->pointB()->x());
        ctrl_pts[i + 1] = (*seg)->pointB()->y();
        ctrl_pts[i + 2] = 0;
      }
      for (i = 0; i < pts_nb - 3; i += 3) {
        glMap1f(GL_MAP1_VERTEX_3, 0, 1, 3, 4, &ctrl_pts[3*i]);
        glBegin(GL_LINE_STRIP);
        for (j = 0; j <= 8; ++j)
          glEvalCoord1f(float(j) / 8);
        glEnd();
      }
      glMap1f(GL_MAP1_VERTEX_3, 0, 1, 3, pts_nb - i, &ctrl_pts[3*i]);
      glBegin(GL_LINE_STRIP);
      for (j = 0; j <= 8; ++j)
        glEvalCoord1f(float(j) / 8);
      glEnd();
      delete[] ctrl_pts;
      
      // Draw the parameters controlling the gaussian
      Stroke *s = (*st);
      double heightA = s->getAParams()[0];
      double heightB = s->getBParams()[0];
      
      double radiusA = s->getAParams()[1];
      double radiusB = s->getBParams()[1];
      
      /*
      glBegin(GL_LINES);
      glVertex3d( s->pointA()->x() - radiusA, s->pointA()->y(), 0.0);
      glVertex3d( s->pointA()->x() + radiusA, s->pointA()->y(), 0.0);
      glEnd();
      glBegin(GL_LINES);
      glVertex3d( s->pointB()->x() - radiusB, s->pointB()->y(), 0.0);
      glVertex3d( s->pointB()->x() + radiusB, s->pointB()->y(), 0.0);
      glEnd();
      */
      
      
      
      Vec3d posA( mirror(s->pointA()->x()), s->pointA()->y(), 0.0);
      Vec3d posB( mirror(s->pointB()->x()), s->pointB()->y(), 0.0);
      
      /*
      glColor3f(1.0,0.0,0.0);
      glBegin(GL_LINES);
      glVertex3d( mirror(s->pointA()->x()), s->pointA()->y()+heightA, 0.0);
      glVertex3d( mirror(s->pointA()->x()), s->pointA()->y(), 0.0);
      glEnd();
      glBegin(GL_LINES);
      glVertex3d( mirror(s->pointB()->x()), s->pointB()->y()+heightB, 0.0);
      glVertex3d( mirror(s->pointB()->x()), s->pointB()->y(), 0.0);
      glEnd();
      */
      
      Vec3d Apa(
          ( *( (*st)->segments()->begin() ) )->pointA()->x(),
          (*((*st)->segments()->begin()))->pointA()->y(),
      0.0
               );
      
      Vec3d Apb(
          (*((*st)->segments()->begin()))->pointB()->x(),
      (*((*st)->segments()->begin()))->pointB()->y(),
      0.0
               );
      
      Vec3d dirA( Apb - Apa );
      
      Vec3d Bpa(
          ( *( (*st)->segments()->rbegin() ) )->pointA()->x(),
      (*((*st)->segments()->rbegin()))->pointA()->y(),
      0.0
               );
      
      Vec3d Bpb(
          (*((*st)->segments()->rbegin()))->pointB()->x(),
      (*((*st)->segments()->rbegin()))->pointB()->y(),
      0.0
               );
      
      Vec3d dirB( Bpa - Bpb );
      
      
      drawFoldEndParameters(posA, dirA,  radiusA,  heightA);
      drawFoldEndParameters(posB, dirB,  radiusB,  heightB);
      
    }
    
    
  }
  
  /////////////// CURRENT BACK FOLD LAYER ///////////////
  
  ch = current_back_gfold_layer->chains()->begin();
  ch_end = current_back_gfold_layer->chains()->end();

  for (; ch != ch_end; ++ch) {

    
    st = (*ch)->strokes()->begin();
    st_end = (*ch)->strokes()->end();
    for (; st != st_end; ++st) {

      // Make obvious green for development
      glColor3f(0.0, 0.7, 0.0);
      glLineWidth(line_width / 2);

      glEnable(GL_LINE_STIPPLE);
      glLineStipple (1, 0xF0F0); // FIXME

      pts_nb = (*st)->segments()->size() + 1;
      ctrl_pts = new float[3 * pts_nb];
      seg = (*st)->segments()->begin();
      seg_end = (*st)->segments()->end();
      if (seg != seg_end) {
        ctrl_pts[0] = mirror((*seg)->pointA()->x());
        ctrl_pts[1] = (*seg)->pointA()->y();
        ctrl_pts[2] = 0;
      }
      for (i = 3; seg != seg_end; ++seg, i += 3) {
        ctrl_pts[i] = mirror((*seg)->pointB()->x());
        ctrl_pts[i + 1] = (*seg)->pointB()->y();
        ctrl_pts[i + 2] = 0;
      }
      for (i = 0; i < pts_nb - 3; i += 3) {
        glMap1f(GL_MAP1_VERTEX_3, 0, 1, 3, 4, &ctrl_pts[3*i]);
        glBegin(GL_LINE_STRIP);
        for (j = 0; j <= 8; ++j)
          glEvalCoord1f(float(j) / 8);
        glEnd();
      }
      glMap1f(GL_MAP1_VERTEX_3, 0, 1, 3, pts_nb - i, &ctrl_pts[3*i]);
      glBegin(GL_LINE_STRIP);
      for (j = 0; j <= 8; ++j)
        glEvalCoord1f(float(j) / 8);
      glEnd();
      delete[] ctrl_pts;
      
      // Draw the parameters controlling the gaussian
      Stroke *s = (*st);
      double heightA = s->getAParams()[0];
      double heightB = s->getBParams()[0];
      
      double radiusA = s->getAParams()[1];
      double radiusB = s->getBParams()[1];
      
      /*
      glBegin(GL_LINES);
      glVertex3d( s->pointA()->x() - radiusA, s->pointA()->y(), 0.0);
      glVertex3d( s->pointA()->x() + radiusA, s->pointA()->y(), 0.0);
      glEnd();
      glBegin(GL_LINES);
      glVertex3d( s->pointB()->x() - radiusB, s->pointB()->y(), 0.0);
      glVertex3d( s->pointB()->x() + radiusB, s->pointB()->y(), 0.0);
      glEnd();
      */
      
      glutils::circle(mirror(s->pointA()->x()), s->pointA()->y(), radiusA);
      glutils::circle(mirror(s->pointB()->x()), s->pointB()->y(), radiusB);
      
      glColor3f(1.0,0.0,0.0);
      glBegin(GL_LINES);
      glVertex3d( mirror(s->pointA()->x()), s->pointA()->y()+heightA, 0.0);
      glVertex3d( mirror(s->pointA()->x()), s->pointA()->y(), 0.0);
      glEnd();
      glBegin(GL_LINES);
      glVertex3d( mirror(s->pointB()->x()), s->pointB()->y()+heightB, 0.0);
      glVertex3d( mirror(s->pointB()->x()), s->pointB()->y(), 0.0);
      glEnd();
      
      glDisable(GL_LINE_STIPPLE);
      
    }
    
    
  }
  
  // JDW Show any fold currently being drawn
  glColor3f(0.0, 1.0, 0.0);
  if (current_front_gfold_layer->currentStroke()) {
    seg = current_front_gfold_layer->currentStroke()->segments()->begin();
    seg_end = current_front_gfold_layer->currentStroke()->segments()->end();
    glBegin(GL_LINE_STRIP);
    if (seg != seg_end)
      glVertex2d(mirror((*seg)->pointA()->x()), (*seg)->pointA()->y());
    for (; seg != seg_end; ++seg)
      glVertex2d(mirror((*seg)->pointB()->x()), (*seg)->pointB()->y());
    glEnd();
  }  
  // JDW Finish drawing the seam layer
  
  // END GAUSSIAN FOLD DRAWING
  
  glColor3f(0.3, 0.3, 0.3);

  if (current_front_layer->currentStroke()) {
    seg = current_front_layer->currentStroke()->segments()->begin();
    seg_end = current_front_layer->currentStroke()->segments()->end();
    glBegin(GL_LINE_STRIP);
    if (seg != seg_end) {
      glVertex2d(mirror((*seg)->pointA()->x()), (*seg)->pointA()->y());
      
    }
    for (; seg != seg_end; ++seg) {
      glVertex2d(mirror((*seg)->pointB()->x()), (*seg)->pointB()->y());
    }
    glEnd();
  }
  // Drawing any mirroring 
  if (_vertical_mirror!= Canvas::OFF && current_front_layer->currentStroke()) {
    seg = current_front_layer->currentStroke()->segments()->begin();
    seg_end = current_front_layer->currentStroke()->segments()->end();
    glBegin(GL_LINE_STRIP);
    if (seg != seg_end) {
      
      glVertex2d(-mirror((*seg)->pointA()->x()), (*seg)->pointA()->y()); // Y-Axis mirroring JDW
    }
    for (; seg != seg_end; ++seg) {
      
      glVertex2d(-mirror((*seg)->pointB()->x()), (*seg)->pointB()->y()); // Y-Axis mirroring JDW
    }
    glEnd();
  }
  
  

  counted_ptr<Point> p;

  if (_chain_pts) {
    ch = current_front_layer->chains()->begin();
    ch_end = current_front_layer->chains()->end();
    glColor4f(0.1, 0.1, 0.1, 0.05);
    glPointSize(2 * current_front_layer->snapping());
    glBegin(GL_POINTS);
    for (; ch != ch_end; ++ch) {
      if (!(*ch)->cycle()) {
	p = (*ch)->pointA();
	if (p.get())
	  glVertex2d(mirror(p->x()), p->y());
	p = (*ch)->pointB();
	if (p.get())
	  glVertex2d(mirror(p->x()), p->y());
      }
    }
    glEnd();
  }

  if (_stroke_pts) {
    ch = current_front_layer->chains()->begin();
    ch_end = current_front_layer->chains()->end();
    glColor4f(0.1, 0.1, 0.1, 0.2);
    glPointSize(utils::min(utils::max(6.0, line_width * 6), 2.0 * current_front_layer->snapping()));
    glBegin(GL_POINTS);
    for (; ch != ch_end; ++ch) {
      st = (*ch)->strokes()->begin();
      st_end = (*ch)->strokes()->end();
      for (; st != st_end; ++st) {
	p = (*st)->pointA();
	if (p.get())
	  glVertex2d(mirror(p->x()), p->y());
	p = (*st)->pointB();
	if (p.get())
	  glVertex2d(mirror(p->x()), p->y());
      }
    }
    glEnd();
  }

  glDisable(GL_DEPTH_TEST);
  
  
  // JDW TEST
  // DRAW GarmentMap Grid
  /*
  GarmentMapList::iterator gmit = _front_garment_maps.begin();
  GarmentMapList::iterator gmit_end = _front_garment_maps.end();
  for ( ; gmit != gmit_end; ++gmit)
  {
    (*gmit)->renderGL(); // testing only
  }
  */

  // Draw textual information.
  glColor3f(0.3, 0.3, 0.3);
  unsigned voffset = 2*((QApplication::font().pixelSize()>0)?QApplication::font().pixelSize():QApplication::font().pointSize());
  if (_front)
    drawText(10, voffset, "<V>iew: front");
  else
    drawText(10, voffset, "<V>iew: back");
  if (_two_sides)
    drawText(10, 2*voffset, "Drawing <M>ode: two sides");
  else
    drawText(10, 2*voffset, "Drawing <M>ode: one side");
  if (_vertical_mirror != Canvas::OFF)
    drawText(10, 3*voffset, "Left/Right M<i>rror: On");
  else
    drawText(10, 3*voffset, "Left/Right M<i>rror: Off");
}

void Canvas::prepareTextures()
{
  // Front texture.
  unsigned char *texture = RepositoryPretreatmentCanvas::frontTexture();
  if (!texture)
    return;

  unsigned texture_size = RepositoryPretreatmentCanvas::frontTextureSize();

  unsigned tex_id;
  glGenTextures(1, &tex_id);
  glBindTexture(GL_TEXTURE_2D, tex_id);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texture_size,
	       texture_size, 0, GL_RGBA, GL_UNSIGNED_BYTE,
	       texture);

  _front_texture_dl = glGenLists(1);
  glNewList(_front_texture_dl, GL_COMPILE);
    
  glEnable(GL_TEXTURE_2D);
  glColor4f(0.95, 0.95, 0.95, 0.1);
  glBindTexture(GL_TEXTURE_2D, tex_id);
  glBegin(GL_TRIANGLE_STRIP);
  glTexCoord2f(0.0, 1.0); // 1
  glVertex2d(-RepositoryPretreatmentCanvas::bbMax()[0], RepositoryPretreatmentCanvas::bbMax()[1]);
  glTexCoord2f(0.0, 0.0); // 2
  glVertex2d(-RepositoryPretreatmentCanvas::bbMax()[0], -RepositoryPretreatmentCanvas::bbMax()[1]);
  glTexCoord2f(1.0, 1.0); // 3
  glVertex2d(RepositoryPretreatmentCanvas::bbMax()[0], RepositoryPretreatmentCanvas::bbMax()[1]);
  glTexCoord2f(1.0, 0.0); // 4
  glVertex2d(RepositoryPretreatmentCanvas::bbMax()[0], -RepositoryPretreatmentCanvas::bbMax()[1]);
  glEnd();
  glDisable(GL_TEXTURE_2D);
    
  glEndList();

  // Back texture.
  texture = RepositoryPretreatmentCanvas::backTexture();
  if (!texture)
    return;

  texture_size = RepositoryPretreatmentCanvas::backTextureSize();

  glGenTextures(1, &tex_id);
  glBindTexture(GL_TEXTURE_2D, tex_id);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texture_size,
	       texture_size, 0, GL_RGBA, GL_UNSIGNED_BYTE,
	       texture);

  _back_texture_dl = glGenLists(1);
  glNewList(_back_texture_dl, GL_COMPILE);
    
  glEnable(GL_TEXTURE_2D);
  glColor4f(0.95, 0.95, 0.95, 0.1);
  glBindTexture(GL_TEXTURE_2D, tex_id);
  glBegin(GL_TRIANGLE_STRIP);
  glTexCoord2f(0.0, 1.0); // 1
  glVertex2d(-RepositoryPretreatmentCanvas::bbMax()[0], RepositoryPretreatmentCanvas::bbMax()[1]);
  glTexCoord2f(0.0, 0.0); // 2
  glVertex2d(-RepositoryPretreatmentCanvas::bbMax()[0], -RepositoryPretreatmentCanvas::bbMax()[1]);
  glTexCoord2f(1.0, 1.0); // 3
  glVertex2d(RepositoryPretreatmentCanvas::bbMax()[0], RepositoryPretreatmentCanvas::bbMax()[1]);
  glTexCoord2f(1.0, 0.0); // 4
  glVertex2d(RepositoryPretreatmentCanvas::bbMax()[0], -RepositoryPretreatmentCanvas::bbMax()[1]);
  glEnd();
  glDisable(GL_TEXTURE_2D);
    
  glEndList();
}

void Canvas::freeTextures()
{
  if (_front_texture_dl)
    glDeleteLists(_front_texture_dl, 1);
  _front_texture_dl = 0;
  if (_back_texture_dl)
    glDeleteLists(_back_texture_dl, 1);
  _back_texture_dl = 0;
}

void Canvas::prepareBBox()
{
  if (!RepositoryPretreatmentCanvas::frontTexture())
    return;

  _bb_dl = glGenLists(1);
  glNewList(_bb_dl, GL_COMPILE);

  glLineWidth(1);
  glutils::fadedBox2D(-RepositoryPretreatmentCanvas::bbMax()[0], -RepositoryPretreatmentCanvas::bbMax()[1], RepositoryPretreatmentCanvas::bbMax()[0], RepositoryPretreatmentCanvas::bbMax()[1]);

  glEndList();
}

void Canvas::freeBBox()
{
  if (_bb_dl)
    glDeleteLists(_bb_dl, 1);
  _bb_dl = 0;
}


void Canvas::pretreatmentUpdated()
{
  unsigned modified = RepositoryPretreatmentCanvas::modified();

  if (modified & RepositoryPretreatmentCanvas::TEXTURES) {
    freeTextures();
    prepareTextures();

    freeBBox();
    prepareBBox();

    freePaperSheet();
    preparePaperSheet();

    freeLayers();
    initLayers();

    freeGarmentMaps();
    initGarmentMaps();

    if (RepositoryPretreatmentCanvas::frontTexture())
      message(MsgHandler::MSG_NORMAL, name(), QString("Textures received from ") +
	      RepositoryPretreatmentCanvas::pretreatment()->name() + ".");
  }
  updateGL();
}

void Canvas::resultUpdated()
{

}

// Internal use only.
inline double slowFastSlow(double t) {
  return sin(M_PI*(t-0.5))*0.5+0.5;
}

// Internal use only.
void strokeAdjust(Stroke *st)
{
  if (!st)
    return;
  DistanceField *df = RepositoryPretreatmentCanvas::distanceField();
  if (!df)
    return;

  unsigned count, i;
  double z1, z2, z3, z4, tmp, alpha;
  Stroke::SegmentList::iterator seg;
  Stroke::SegmentList::iterator seg_begin;
  Stroke::SegmentList::iterator seg_end;
  Stroke::SegmentList::iterator seg1;
  Stroke::SegmentList::iterator seg2;
  seg_begin = st->segments()->begin();
  seg_end = st->segments()->end();
  seg = seg_begin;
  while (seg != seg_end) {
    while (seg != seg_end && (*seg)->type() == Segment::IN)
      ++seg;
    if (seg == seg_end)
      break;
    seg1 = seg;
    count = 0;
    do {
      ++count;
      ++seg;
    } while (seg != seg_end && (*seg)->type() != Segment::IN);
    seg2 = seg;
    --seg2;
    if (seg1 == seg2)
      continue;
    if (seg1-- == seg_begin || ++seg2 == seg_end)
      continue;
    z1 = (*seg1)->pointB()->z();
    z2 = (*seg1)->pointB()->z() + ((*seg1)->pointB()->z() - (*seg1)->pointA()->z()) * count / 3.0;
    z3 = (*seg2)->pointA()->z() + ((*seg2)->pointA()->z() - (*seg2)->pointB()->z()) * count / 3.0;
    z4 = (*seg2)->pointA()->z();
    ++seg1;
    --seg2;
    for (i = 0; seg1 != seg2; ++seg1, ++i) {
      alpha = static_cast<double>(i+1)/(count+1);
      utils::bezier3(tmp, z1, z2, z3, z4, alpha); // FIXME
      (*seg1)->pointB()->setZ(tmp);

      point_utils::distFromZ((*seg1)->pointB()->x(), (*seg1)->pointB()->y(), (*seg1)->pointB()->z(), tmp);
      
      (*seg1)->pointB()->setDistance(tmp);
    }
  }
}

//     z1 = (*seg1)->pointB()->z();
//     z2 = (*seg2)->pointA()->z();
//     dz1 = z1 - (*seg1)->pointA()->z();
//     dz2 = (*seg2)->pointB()->z() - z2;
//     dz1_tmp = 2 * (z2 - z1) / count - dz2;
//     dz2_tmp = 2 * (z2 - z1) / count - dz1;
//     dz1 = (dz1 + dz1_tmp) / 2;
//     dz2 = (dz2 + dz2_tmp) / 2;
//     ++seg1;
//     --seg2;
//     dz = new double[count-1];
//     for (i = 0; i < count-1; ++i) {
//       alpha = static_cast<double>(i+1)/(count+1);
//       dz[i] = (1 - alpha) * dz1 + alpha * dz2;
//     }
//     for (i = 0; seg1 != seg2; ++seg1, ++i) {
//       (*seg1)->pointB()->setZ((*seg1)->pointA()->z() + dz[i]);
//       point_utils::distFromZ((*seg1)->pointB()->x(), (*seg1)->pointB()->y(), (*seg1)->pointB()->z(), dz1_tmp);
//       (*seg1)->pointB()->setDistance(dz1_tmp);
//     }
//     delete[] dz;
//   }
// }

// Internal use only.
void strokeSmooth(Stroke *st)
{
  if (!st)
    return;
  DistanceField *df = RepositoryPretreatmentCanvas::distanceField();
  if (!df)
    return;

  Vec3d mask(1.0, 1.0, 1.0);
  unsigned i, j;
  unsigned size = st->segments()->size();
  double *z = new double[size];
  double *new_z = new double[size];
  double *dist = new double[size];
  double alpha, tmp;
  Stroke::SegmentList::iterator seg;
  Stroke::SegmentList::iterator seg_end;

  seg = st->segments()->begin();
  seg_end = st->segments()->end();
  for (i = 0; seg != seg_end; ++seg, ++i) {
    z[i] = (*seg)->pointA()->z();
    dist[i] = (*seg)->pointA()->distance();
  }

  for (i = 0; i < 100; ++i) { // FIXME: nb of iterations
    for (j = 1; j < size-1; ++j) {
      alpha = dist[j] / df->maxDist();
      tmp = mask[0] * z[j-1] +
	mask[1] * z[j] +
	mask[2] * z[j+1];
      tmp /= mask[0] + mask[1] + mask[2];
      new_z[j] = alpha * tmp + (1 - alpha) * z[j];
    }
    for (j = 1; j < size-1; ++j) {
      z[j] = new_z[j];
    }
  }
  seg = st->segments()->begin();
  seg_end = st->segments()->end();
  for (++seg, i = 1; seg != seg_end; ++seg, ++i) {
    (*seg)->pointA()->setZ(z[i]);

    point_utils::distFromZ((*seg)->pointA()->x(), (*seg)->pointA()->y(), (*seg)->pointA()->z(), tmp);
    // FIXME:    
//     if ((*seg)->pointA()->type() != Point::OUT && tmp <= point_utils::zEpsilon()) {
//       (*seg)->pointA()->setDistance(point_utils::zEpsilon());
//       point_utils::zFromDist((*seg)->pointA()->x(), (*seg)->pointA()->y(), point_utils::zEpsilon(), tmp);
//       (*seg)->pointA()->setZ(tmp);
//     }
//     else
    (*seg)->pointA()->setDistance(tmp);
  }
  delete[] z;
  delete[] new_z;
  delete[] dist;
}

void Canvas::computeSkeleton()
{
  Layer::ChainList::iterator ch;
  Layer::ChainList::iterator ch_end;
  Chain::StrokeList::iterator st;
  Chain::StrokeList::iterator st_end;
  Stroke::SegmentList::iterator seg;
  Stroke::SegmentList::iterator seg_end;
  double dist, z, partial_length, alpha;

  // Front layer.

  ch = _front_layer->chains()->begin();
  ch_end = _front_layer->chains()->end();
  for (; ch != ch_end; ++ch) {

    // Discard non-cyclic chains.
    if (!(*ch)->cycle())
      continue;

    // First pass on the strokes, to set up length and type, and extremities.
    st = (*ch)->strokes()->begin();
    st_end = (*ch)->strokes()->end();
    for (; st != st_end; ++st) {
      seg = (*st)->segments()->begin();
      seg_end = (*st)->segments()->end();
      for (; seg != seg_end; ++seg) {
        // Set the point type from the canvas texture (IN, OUT, BORDER, UNKNOWN)
	(*seg)->pointA()->setType(point_utils::type((*seg)->pointA()->x(), (*seg)->pointA()->y()));
	(*seg)->pointB()->setType(point_utils::type((*seg)->pointB()->x(), (*seg)->pointB()->y()));
	(*seg)->updateLength();
	(*seg)->updateType();
      }
      // find minimum distance to surface along z direction for this x and y
      point_utils::minDist((*st)->pointA()->x(), (*st)->pointA()->y(), dist, z, true);
      (*st)->pointA()->setDistance(dist);
      //      point_utils::zFromDist((*st)->pointA()->x(), (*st)->pointA()->y(), dist, z, true);
      (*st)->pointA()->setZ(z);
      (*st)->updateLength();
      (*st)->updateType();
//      std::cout << "Set Distance for stroke point: " << dist << std::endl;
//      std::cout << "Set Z for stroke point: " << z << std::endl;
    }

    // Second pass on the strokes.
    //  Interpolate z dist along strokes - use to set distance
    // 
    st = (*ch)->strokes()->begin();
    st_end = (*ch)->strokes()->end();
    for (; st != st_end; ++st) {
      if ((*st)->type() == Stroke::OUT) {
	partial_length = 0;
	seg = (*st)->segments()->begin();
	seg_end = (*st)->segments()->end();
	for (; seg != seg_end; ++seg) {
	  partial_length += (*seg)->length();
	  alpha = partial_length / (*st)->length();
	  (*seg)->pointB()->setZ((1 - alpha) * (*st)->pointA()->z() +
				 alpha * (*st)->pointB()->z());
   
   point_utils::distFromZ((*seg)->pointB()->x(), (*seg)->pointB()->y(), (*seg)->pointB()->z(), dist);
	  (*seg)->pointB()->setDistance(dist);
	}
      }
      else {
	partial_length = 0;
	seg = (*st)->segments()->begin();
	seg_end = (*st)->segments()->end();
	for (; (*seg)->pointB().get() != (*st)->pointB().get(); ++seg) {
	  partial_length += (*seg)->length();
	  alpha = partial_length / (*st)->length();
	  (*seg)->pointB()->setDistance((1 - alpha) * (*st)->pointA()->distance() +
					alpha * (*st)->pointB()->distance());
	  point_utils::zFromDist((*seg)->pointB()->x(), (*seg)->pointB()->y(), (*seg)->pointB()->distance(), z, true);
	  (*seg)->pointB()->setZ(z);
   
	}
	if ((*st)->type() == Stroke::MIXED) {
	  seg = (*st)->segments()->begin();
	  seg_end = (*st)->segments()->end();
	  for (; seg != seg_end; ++seg) {
	    (*seg)->pointA()->setType(point_utils::type((*seg)->pointA()->x(), (*seg)->pointA()->y(), (*seg)->pointA()->distance()));
	    (*seg)->pointB()->setType(point_utils::type((*seg)->pointB()->x(), (*seg)->pointB()->y(), (*seg)->pointB()->distance()));
	    (*seg)->updateType();
	  }
	  strokeAdjust(*st);
	  strokeSmooth(*st);
	}
      }
    }
  }

  // Back layer.

  ch = _back_layer->chains()->begin();
  ch_end = _back_layer->chains()->end();
  for (; ch != ch_end; ++ch) {

    // Discard non-cyclic chains.
    if (!(*ch)->cycle())
      continue;

    // First pass on the strokes, to set up length and type, and extremities.
    st = (*ch)->strokes()->begin();
    st_end = (*ch)->strokes()->end();
    for (; st != st_end; ++st) {
      seg = (*st)->segments()->begin();
      seg_end = (*st)->segments()->end();
      for (; seg != seg_end; ++seg) {
	(*seg)->pointA()->setType(point_utils::type((*seg)->pointA()->x(), (*seg)->pointA()->y()));
	(*seg)->pointB()->setType(point_utils::type((*seg)->pointB()->x(), (*seg)->pointB()->y()));
	(*seg)->updateLength();
	(*seg)->updateType();
      }
      point_utils::minDist((*st)->pointA()->x(), (*st)->pointA()->y(), dist, z, false);
      (*st)->pointA()->setDistance(dist);
      //      point_utils::zFromDist((*st)->pointA()->x(), (*st)->pointA()->y(), dist, z, false);
      (*st)->pointA()->setZ(z);
      (*st)->updateLength();
      (*st)->updateType();
    }

    // Second pass on the strokes.
    st = (*ch)->strokes()->begin();
    st_end = (*ch)->strokes()->end();
    for (; st != st_end; ++st) {
      if ((*st)->type() == Stroke::OUT) {
	partial_length = 0;
	seg = (*st)->segments()->begin();
	seg_end = (*st)->segments()->end();
	for (; seg != seg_end; ++seg) {
	  partial_length += (*seg)->length();
	  alpha = partial_length / (*st)->length();
	  (*seg)->pointB()->setZ((1 - alpha) * (*st)->pointA()->z() +
				 alpha * (*st)->pointB()->z());
   
	  point_utils::distFromZ((*seg)->pointB()->x(), (*seg)->pointB()->y(), (*seg)->pointB()->z(), dist);
	  (*seg)->pointB()->setDistance(dist);
	}
      }
      else {
	partial_length = 0;
	seg = (*st)->segments()->begin();
	seg_end = (*st)->segments()->end();
	for (; (*seg)->pointB().get() != (*st)->pointB().get(); ++seg) {
	  partial_length += (*seg)->length();
	  alpha = partial_length / (*st)->length();
	  (*seg)->pointB()->setDistance((1 - alpha) * (*st)->pointA()->distance() +
					alpha * (*st)->pointB()->distance());
	  point_utils::zFromDist((*seg)->pointB()->x(), (*seg)->pointB()->y(), (*seg)->pointB()->distance(), z, false);
	  (*seg)->pointB()->setZ(z);
	}
	if ((*st)->type() == Stroke::MIXED) {
	  seg = (*st)->segments()->begin();
	  seg_end = (*st)->segments()->end();
	  for (; seg != seg_end; ++seg) {
	    (*seg)->pointA()->setType(point_utils::type((*seg)->pointA()->x(), (*seg)->pointA()->y(), (*seg)->pointA()->distance()));
	    (*seg)->pointB()->setType(point_utils::type((*seg)->pointB()->x(), (*seg)->pointB()->y(), (*seg)->pointB()->distance()));
	    (*seg)->updateType();
	  }
	  strokeAdjust(*st);
	  strokeSmooth(*st);
	}
      }
    }
  }
}

void Canvas::clearSeams() {
  listOfFrontSeamPointLists.clear();
  listOfBackSeamPointLists.clear();
}

void Canvas::computeSurface()
{
  Layer::ChainList::iterator ch;
  Layer::ChainList::iterator ch_end;
  Chain::StrokeList::iterator st;
  Chain::StrokeList::iterator st_end;
  Stroke::SegmentList::iterator seg;
  Stroke::SegmentList::iterator seg_end;

  Mat33d mask;
  mask(0, 0) = 1;
  mask(0, 1) = 1;
  mask(0, 2) = 1;
  mask(1, 0) = 1;
  mask(1, 1) = 1;
  mask(1, 2) = 1;
  mask(2, 0) = 1;
  mask(2, 1) = 1;
  mask(2, 2) = 1;

  GarmentMap *db;

  freeGarmentMaps();
  initGarmentMaps();
  clearSeams();

  // Front garment maps.

  ch = _front_layer->chains()->begin();
  ch_end = _front_layer->chains()->end();
  for (; ch != ch_end; ++ch) {
    // Discard non-cyclic chains.
    if (!(*ch)->cycle())
      continue;
    db = new GarmentMap(_gm_size_x, _gm_size_y);
    db->compute(true, *ch, mask, mask, 0, 0.01);
    _front_garment_maps.push_back(db);
  }

  // Back garment maps.

  ch = _back_layer->chains()->begin();
  ch_end = _back_layer->chains()->end();
  for (; ch != ch_end; ++ch) {
    // Discard non-cyclic chains.
    if (!(*ch)->cycle())
      continue;
    db = new GarmentMap(_gm_size_x, _gm_size_y);
    db->compute(false, *ch, mask, mask, 0, 0.01);
    _back_garment_maps.push_back(db);
  }
  
  // JDW Compute seam lines from any chains in each seam layer
  // using the z buffer from the appropriate garment map.
  // Simple implementation determines which garment map the seam belongs to
  // using a boundary check.
  
  // First determine extents of any seams in front layer, and then work out
  // which is the first (for this simple implementation) garment map which it 
  // could belong to
  
    // Front garment maps.

  ch = _front_seam_layer->chains()->begin();
  ch_end = _front_seam_layer->chains()->end();
  
  
  // For each chain in the front seam layer
  for (; ch != ch_end; ++ch) {


    // FIXME Each chain is a seam. 
    Vec3d seamMinBB(999999,999999,0), seamMaxBB(-9999999,-999999,0); // store maximum extents of the seam
    seamMinBB = (*ch)->bbMin();
    seamMaxBB = (*ch)->bbMax();
   
    // Compare seams to available GarmentMaps. Pick first overlapping set
    if (RepositoryCanvasResult::frontGarmentMaps()) {
      Canvas::GarmentMapList::const_iterator gm = RepositoryCanvasResult::frontGarmentMaps()->begin();
      Canvas::GarmentMapList::const_iterator gm_end = RepositoryCanvasResult::frontGarmentMaps()->end();
      int gmCount = 0;
      bool foundOverlappingSeam = false;
      for (; gm != gm_end; ++gm) {
        std::cout << "GM: " << gmCount++ << std::endl;
        std::cout << "GMBBMin: " << (*gm)->bbMin() << std::endl;
        std::cout << "GMBBMax: " << (*gm)->bbMax() << std::endl;
        if( (   seamMaxBB[0]  > (*gm)->bbMin()[0]
            && seamMinBB[0]  < (*gm)->bbMax()[0]) // overlap in X plane
                &&                                 // AND
            (   seamMaxBB[1]  > (*gm)->bbMin()[1]
            && seamMinBB[1]  < (*gm)->bbMax()[1]) // overlap in Y plane
        ) 
        {
          // We have a seam for this garment map. End here
            foundOverlappingSeam = true;
            break;
        }
      } // end loop over GarmentMaps
      
      if(foundOverlappingSeam) {
        // Act on this Seam and GarmentMap
               
          st     = (*ch)->strokes()->begin();
          st_end = (*ch)->strokes()->end();
          
          std::list<Point> pointList; // to be added to the list of lists for later file storage
          
          for (; st != st_end; ++st) {
            seg     = (*st)->segments()->begin();
            seg_end = (*st)->segments()->end();
            
            for (; seg != seg_end; ++seg) {
              //std::cout << "Pt[x]: " << (*seg)->pointA()->x() << std::endl;
              //std::cout << "Pt[y]: " << (*seg)->pointA()->y() << std::endl;
              
              // Use integer lookup to find Z from depth buffer
              // FIXME add range checks to ensure seam outside of garment doesn't index
              // outside zbuffer
              unsigned i = static_cast<unsigned>((((*seg)->pointA()->x() - (*gm)->bbMin()[0]) / ((*gm)->bbMax()[0] - (*gm)->bbMin()[0])) * (*gm)->sizeX());
              unsigned j = static_cast<unsigned>((((*seg)->pointA()->y() - (*gm)->bbMin()[1]) / ((*gm)->bbMax()[1] - (*gm)->bbMin()[1])) * (*gm)->sizeY());
              //std::cout << "Pt[z]: " << (*gm)->z(i, j) << std::endl;
              
              (*seg)->pointA()->setZ( (*gm)->z(i, j) );
              
              
              pointList.push_back( (*(*seg)->pointA()) );
              
        
              
              //  do pointB only for last segment in chain...
              
              if((*seg) == (*st)->segments()->back())
              {
                unsigned i = static_cast<unsigned>((((*seg)->pointB()->x() - (*gm)->bbMin()[0]) / ((*gm)->bbMax()[0] - (*gm)->bbMin()[0])) * (*gm)->sizeX());
                unsigned j = static_cast<unsigned>((((*seg)->pointB()->y() - (*gm)->bbMin()[1]) / ((*gm)->bbMax()[1] - (*gm)->bbMin()[1])) * (*gm)->sizeY());
                (*seg)->pointB()->setZ( (*gm)->z(i,j) );
                pointList.push_back( (*(*seg)->pointB()) );
              }

            } // end seg loop

          } // end stroke loop
          listOfFrontSeamPointLists.push_back( pointList );
        } // end chain loop
        

        
    } // end if GarmentMap exists block
  } // end loop over chains in seam layer

  
  // Determine extents of any seams in back layer, and then work out
  // which is the first (for this simple implementation) garment map which it 
  // could belong to
  
  // Back garment maps.

  ch = _back_seam_layer->chains()->begin();
  ch_end = _back_seam_layer->chains()->end();
  
  
  // For each chain in the back seam layer
  for (; ch != ch_end; ++ch) {


    // FIXME Each chain is a seam. 
    Vec3d seamMinBB(999999,999999), seamMaxBB(-9999999,-999999); // store maximum extents of the seam
    seamMinBB = (*ch)->bbMin();
    seamMaxBB = (*ch)->bbMax();
 
 
    // Compare seams to available GarmentMaps. Pick first overlapping set
    if (RepositoryCanvasResult::backGarmentMaps()) {
      Canvas::GarmentMapList::const_iterator gm = RepositoryCanvasResult::backGarmentMaps()->begin();
      Canvas::GarmentMapList::const_iterator gm_end = RepositoryCanvasResult::backGarmentMaps()->end();
      int gmCount = 0;
      bool foundOverlappingSeam = false;
      for (; gm != gm_end; ++gm) {
        std::cout << "GM: " << gmCount++ << std::endl;
        std::cout << "GMBBMin: " << (*gm)->bbMin() << std::endl;
        std::cout << "GMBBMax: " << (*gm)->bbMax() << std::endl;
        if( (   seamMaxBB[0]  > (*gm)->bbMin()[0]
                && seamMinBB[0]  < (*gm)->bbMax()[0]) // overlap in X plane
                &&                                 // AND
                (   seamMaxBB[1]  > (*gm)->bbMin()[1]
                && seamMinBB[1]  < (*gm)->bbMax()[1]) // overlap in Y plane
          ) 
        {
          // We have a seam for this garment map. End here
          foundOverlappingSeam = true;
          break;
        }
      } // end loop over GarmentMaps
      
      if(foundOverlappingSeam) {
        // Act on this Seam and GarmentMap
               
        st     = (*ch)->strokes()->begin();
        st_end = (*ch)->strokes()->end();
          
        std::list<Point> pointList; // to be added to the list of lists for later file storage
          
        for (; st != st_end; ++st) {
          seg     = (*st)->segments()->begin();
          seg_end = (*st)->segments()->end();
            
          for (; seg != seg_end; ++seg) {
              //std::cout << "Pt[x]: " << (*seg)->pointA()->x() << std::endl;
              //std::cout << "Pt[y]: " << (*seg)->pointA()->y() << std::endl;
              
              // Use integer lookup to find Z from depth buffer
              // FIXME add range checks to ensure seam outside of garment doesn't index
              // outside zbuffer
            unsigned i = static_cast<unsigned>((((*seg)->pointA()->x() - (*gm)->bbMin()[0]) / ((*gm)->bbMax()[0] - (*gm)->bbMin()[0])) * (*gm)->sizeX());
            unsigned j = static_cast<unsigned>((((*seg)->pointA()->y() - (*gm)->bbMin()[1]) / ((*gm)->bbMax()[1] - (*gm)->bbMin()[1])) * (*gm)->sizeY());
              //std::cout << "Pt[z]: " << (*gm)->z(i, j) << std::endl;
              
            (*seg)->pointA()->setZ( (*gm)->z(i, j) );
              
              
            pointList.push_back( (*(*seg)->pointA()) );
              
        
              
              //  do pointB only for last segment in chain...
              
            if((*seg) == (*st)->segments()->back())
            {
              unsigned i = static_cast<unsigned>((((*seg)->pointB()->x() - (*gm)->bbMin()[0]) / ((*gm)->bbMax()[0] - (*gm)->bbMin()[0])) * (*gm)->sizeX());
              unsigned j = static_cast<unsigned>((((*seg)->pointB()->y() - (*gm)->bbMin()[1]) / ((*gm)->bbMax()[1] - (*gm)->bbMin()[1])) * (*gm)->sizeY());
              (*seg)->pointB()->setZ( (*gm)->z(i,j) );
              pointList.push_back( (*(*seg)->pointB()) );
            }

          } // end seg loop

        } // end stroke loop
        listOfBackSeamPointLists.push_back( pointList );
      } // end chain loop
        

        
    } // end if GarmentMap exists block
  }
  
  //////// Modify front GarmentMaps according to fold information
  
  ch = _front_gfold_layer->chains()->begin();
  ch_end = _front_gfold_layer->chains()->end();
  
  
  // For each chain in the front gfold layer
  for (; ch != ch_end; ++ch) {

    Vec3d chMinBB, chMaxBB; 
    (*ch)->recalculateBB();
    chMinBB = (*ch)->bbMin();
    chMaxBB = (*ch)->bbMax();
   
    // Compare folds to available GarmentMaps. Pick first overlapping set
    if (RepositoryCanvasResult::frontGarmentMaps()) {
      Canvas::GarmentMapList::const_iterator gm  =
          RepositoryCanvasResult::frontGarmentMaps()->begin();
      Canvas::GarmentMapList::const_iterator gm_end =
          RepositoryCanvasResult::frontGarmentMaps()->end();
      int gmCount = 0;
      bool foundOverlapping = false;
      for (; gm != gm_end; ++gm) {
        //std::cout << "GM: " << gmCount++ << std::endl;
        //std::cout << "GMBBMin: " << (*gm)->bbMin() << std::endl;
        //std::cout << "GMBBMax: " << (*gm)->bbMax() << std::endl;
        if( (   chMaxBB[0]  > (*gm)->bbMin()[0]
                && chMinBB[0]  < (*gm)->bbMax()[0]) // overlap in X plane
                &&                                 // AND
                (   chMaxBB[1]  > (*gm)->bbMin()[1]
                && chMinBB[1]  < (*gm)->bbMax()[1]) // overlap in Y plane
          ) 
        {
          // We have a seam for this garment map. End here
          foundOverlapping = true;
          break;
        }
      } // end loop over GarmentMaps
      
      if(foundOverlapping) {
        // Act on this fold and GarmentMap
               
        st     = (*ch)->strokes()->begin();
        st_end = (*ch)->strokes()->end();
          
        
        for (; st != st_end; ++st) {

            
          (*gm)->gfold( (*st), true ); // concave fold on the front


        } // end stroke loop
        
      } // end chain loop
        

        
    } // end if GarmentMap exists block
  } // end loop over chains in front gfold layer

  //////// Modify back GarmentMaps according to fold information
  
  ch = _back_gfold_layer->chains()->begin();
  ch_end = _back_gfold_layer->chains()->end();
  
  
  // For each chain in the back gfold layer
  for (; ch != ch_end; ++ch) {

    Vec3d chMinBB, chMaxBB; 
    (*ch)->recalculateBB();
    chMinBB = (*ch)->bbMin();
    chMaxBB = (*ch)->bbMax();
   
    // Compare folds to available GarmentMaps. Pick first overlapping set
    if (RepositoryCanvasResult::backGarmentMaps()) {
      Canvas::GarmentMapList::const_iterator gm  =
          RepositoryCanvasResult::backGarmentMaps()->begin();
      Canvas::GarmentMapList::const_iterator gm_end =
          RepositoryCanvasResult::backGarmentMaps()->end();
      int gmCount = 0;
      bool foundOverlapping = false;
      for (; gm != gm_end; ++gm) {
        //std::cout << "GM: " << gmCount++ << std::endl;
        //std::cout << "GMBBMin: " << (*gm)->bbMin() << std::endl;
        //std::cout << "GMBBMax: " << (*gm)->bbMax() << std::endl;
        if( (   chMaxBB[0]  > (*gm)->bbMin()[0]
                && chMinBB[0]  < (*gm)->bbMax()[0]) // overlap in X plane
                &&                                 // AND
                (   chMaxBB[1]  > (*gm)->bbMin()[1]
                && chMinBB[1]  < (*gm)->bbMax()[1]) // overlap in Y plane
          ) 
        {
          // We have a seam for this garment map. End here
          foundOverlapping = true;
          break;
        }
      } // end loop over GarmentMaps
      
      if(foundOverlapping) {
        // Act on this fold and GarmentMap
               
        st     = (*ch)->strokes()->begin();
        st_end = (*ch)->strokes()->end();
          
        
        for (; st != st_end; ++st) {

            
          (*gm)->gfold( (*st), false ); // fold on the back


        } // end stroke loop
        //listOfFrontSeamPointLists.push_back( pointList );
      } // end chain loop
        

        
    } // end if GarmentMap exists block
  } // end loop over chains in back gfold layer
  
} //end computeSurface()

void Canvas::computeMesh()
{
  computeSkeleton();
  computeSurface();
}

void Canvas::drawFoldEndParameters(const Vec3d& pos, const Vec3d& dir, const double& radius, const double& height)
{
  //glMatrixMode(GL_MODELVIEW);
  Vec3d direction = dir;
  direction.normalize();
  
  glPushMatrix();
  glutils::circle(pos.x(), pos.y(), radius);
  Vec3d yAxis(0.0,1.0,0.0);
  double ctheta = (yAxis * direction ) / (direction.norm() * yAxis.norm());
  //std::cout<< "Theta: " << ctheta << std::endl;
  
  double rotDir = 1.0;
  if(dir.x() < 0.0) 
    rotDir = -1.0;
    
  
  glTranslatef( pos.x(), pos.y(), 0.0 );
  glRotatef( rotDir * -acos(ctheta) / (M_PI * 2.0) * 360., 0.0,0.0,1.0); 
  
  // indicate height
  /*
  glColor3f(1.0,0.0,1.0);
  glBegin(GL_LINES);
  glVertex3f( 0.0,0.0,0.0 );
  glVertex3f( 0.0,-height,0.0 );
  glEnd();
  */
  glColor3f(1.0,0.0,1.0);
  // sample gaussian along x direction
  // 3*sigma is equal to radius
  int nbSamples = 10;
  double sigma = radius / 3.;
  double dx = radius / static_cast<double>(nbSamples-1);
  
  double x,y;
  glBegin(GL_LINE_STRIP);
  for(int i=-(nbSamples-1);i<nbSamples;i++) 
  {
    x = i * dx;
    y = utils::gauss2d( -height, sigma, x, 0.0);
    glVertex3f(x,y,0.0); 
  }
  glEnd();
  
  glPopMatrix();
  

}