//
//  Filename         : pretreatment.cpp
//  Author           : Emmanuel Turquin
//  Purpose          : A generator of df and textures needed by the algorithm.
//  Date of creation : 04/20/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include <cmath>
#include <lib3ds/mesh.h>
#include <lib3ds/material.h>
#include <lib3ds/matrix.h>
#include <lib3ds/vector.h>
#include <qfiledialog.h>
#include <qcursor.h>
#include "chrono.h"
#include "config.h"
#include "config_app_info.h"
#include "config_pretreatment.h"
#include "utils.h"
#include "glutils.h"
#include "mdi_window.h"
#include "lcmodel.h"
#include "octree.h"
#include "distance_field.h"
#include "repository_pretreatment_canvas.h"
#include "repository_pretreatment_result.h"
#include "canvas.h"
#include "result.h"
#include "offscreen.h"
#include "pretreatment.h"
#include "pattern.h"
#include "vecfield.h"

#include <qimage.h> // DEBUG

static offscreen::OffScreenArea	*_offscreen;

using namespace qglviewer;

Pretreatment::Pretreatment(MDIWindow* parent, const char* name,
			   const QGLWidget* shareWidget, int wflags)
  : BaseQGLViewer(parent, name, shareWidget, wflags),
    _file(NULL)
{
  _mdi = parent;
  _bb_dl = 0;
  _model_dl = 0;
  _df_dl = 0;
  _octree = 0;
  _octree_dl = 0;
  _lcmodel = 0;
  _offscreen = 0;
  updateConfig();
  RepositoryPretreatmentCanvas::setPretreatment(this);
  RepositoryPretreatmentResult::setPretreatment(this);
}

Pretreatment::~Pretreatment()
{
  delete _offscreen;
  freeModel();
  freeBBox();
  freeOctree();
  freeDistanceField();
  freeTextures();
  if (_mdi && _mdi->config()) {
    _mdi->config()->io()->setValue("bbox/size/x", config::PRETREATMENT_BBOX_X);
    _mdi->config()->io()->setValue("bbox/size/y", config::PRETREATMENT_BBOX_Y);
    _mdi->config()->io()->setValue("bbox/size/z", config::PRETREATMENT_BBOX_Z);
  }
  if (RepositoryPretreatmentCanvas::pretreatment() == this)
    RepositoryPretreatmentCanvas::setPretreatment(0);
  if (RepositoryPretreatmentResult::pretreatment() == this)
    RepositoryPretreatmentResult::setPretreatment(0);
}

void Pretreatment::updateConfig()
{
  makeCurrent();

  bool persp = config::PRETREATMENT_CAMERA_PERSPECTIVE;
  bool reset = false;
  bool bb_modified = false;
  bool df_modified = false;
  bool compute = false;
  int octree_level = config::PRETREATMENT_OCTREE_LEVEL;
  unsigned df_x = config::PRETREATMENT_DF_X;
  unsigned df_y = config::PRETREATMENT_DF_Y;
  unsigned df_z = config::PRETREATMENT_DF_Z;
  unsigned df_x_min = 0;
  unsigned df_y_min = 0;
  unsigned df_z_min = 0;
  unsigned df_x_max = config::PRETREATMENT_DF_X - 1;
  unsigned df_y_max = config::PRETREATMENT_DF_Y - 1;
  unsigned df_z_max = config::PRETREATMENT_DF_Z - 1;
  unsigned texture_size = config::PRETREATMENT_TEXTURE_SIZE;
  _disp_df = config::PRETREATMENT_DISP_DF;
  _disp_bb = config::PRETREATMENT_DISP_BBOX;
  _disp_axis = config::PRETREATMENT_DISP_AXIS;
  _disp_model = config::PRETREATMENT_DISP_MODEL;
  _disp_octree = config::PRETREATMENT_DISP_OCTREE;
  _bb_size = RepositoryPretreatmentResult::modelBBMax() - RepositoryPretreatmentResult::modelBBMin();
  int df_slice_z=33;
  if (_mdi && _mdi->config()) {
    _mdi->config()->io()->getValue("camera/perspective", persp);
    _mdi->config()->io()->getValue("display/bbox", _disp_bb);
    _mdi->config()->io()->getValue("display/bbox_axis", _disp_axis);
    _mdi->config()->io()->getValue("display/model", _disp_model);
    _mdi->config()->io()->getValue("display/octree", _disp_octree);
    _mdi->config()->io()->getValue("display/df", _disp_df);
    _mdi->config()->io()->getValue("octree/level", octree_level);
    _mdi->config()->io()->getValue("compute/start", compute);
    _mdi->config()->io()->getValue("bbox/modified", bb_modified);
    _mdi->config()->io()->getValue("bbox/reset", reset);
    _mdi->config()->io()->getValue("df/size/x", df_x);
    _mdi->config()->io()->getValue("df/size/y", df_y);
    _mdi->config()->io()->getValue("df/size/z", df_z);
    _mdi->config()->io()->getValue("display/dfsize/x_min", df_x_min);
    _mdi->config()->io()->getValue("display/dfsize/x_max", df_x_max);
    _mdi->config()->io()->getValue("display/dfsize/y_min", df_y_min);
    _mdi->config()->io()->getValue("display/dfsize/y_max", df_y_max);
    _mdi->config()->io()->getValue("display/dfsize/z_min", df_z_min);
    _mdi->config()->io()->getValue("display/dfsize/z_max", df_z_max);
    _mdi->config()->io()->getValue("display/dfsize/modified", df_modified);
    _mdi->config()->io()->getValue("texture/size", texture_size);
    _mdi->config()->io()->getValue("display/dfslice/percent_z", df_slice_z);
    _df_slice_z = static_cast<float>((99-df_slice_z) /99.0); // slider range is 0-99, invert value to match orientation of model

    if (reset) {
      _mdi->config()->io()->setValue("bbox/reset", false);
      _mdi->config()->io()->setValue("bbox/size/x", _bb_size[0]);
      _mdi->config()->io()->setValue("bbox/size/y", _bb_size[1]);
      _mdi->config()->io()->setValue("bbox/size/z", _bb_size[2]);
      _mdi->configUpdated();
      manipulatedFrame()->setTranslationAndRotation(sceneCenter(), Quaternion());
    }
    else {
      _mdi->config()->io()->getValue("bbox/size/x", _bb_size[0]);
      _mdi->config()->io()->getValue("bbox/size/y", _bb_size[1]);
      _mdi->config()->io()->getValue("bbox/size/z", _bb_size[2]);
    }

    if (bb_modified) {
      _mdi->config()->io()->setValue("bbox/modified", false);
      freeBBox();
      _bb_dl = prepareBBox();
      /* JDW to display slice of distance field
      if(RepositoryPretreatmentCanvas::distanceField()) 
        DistanceField::initSlice(RepositoryPretreatmentCanvas::distanceField(),_df_slice_z);
      */
    }

    if (df_modified) {
      _mdi->config()->io()->setValue("display/dfsize/modified", false);
      if (_df_dl)
	glDeleteLists(_df_dl, 1);
      _df_dl = prepareDistanceField(df_x_min, df_x_max, df_y_min, df_y_max, df_z_min, df_z_max);
    }

    if (persp)
	camera()->setType(Camera::PERSPECTIVE);
    else
	camera()->setType(Camera::ORTHOGRAPHIC);

    if (compute) { // Start computation.
      QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
      _mdi->config()->io()->setValue("compute/start", false);
      freeOctree();
      freeDistanceField();
      freeTextures();
      computeTextures(texture_size); // Has to be computed before the distance field. 
      computeOctree(octree_level);
      computeDistanceField(df_x, df_y, df_z);
      _octree_dl = prepareOctree();
      _df_dl = prepareDistanceField(df_x_min, df_x_max, df_y_min, df_y_max, df_z_min, df_z_max);
      QApplication::restoreOverrideCursor();
      /* JDW to display slice of distance field
      if(RepositoryPretreatmentCanvas::distanceField()) 
        DistanceField::initSlice(RepositoryPretreatmentCanvas::distanceField(),_df_slice_z);
      */
    }
    updateGL();
  }
}

void Pretreatment::loadFile()
{
  QString path;
  if (_mdi && _mdi->config())
    _mdi->config()->io()->getValue("paths/file_3ds", path);
  QString filename = QFileDialog::getOpenFileName(path, "3DS files (*.3ds *.3DS);;All files (*)", this);
  
  // In case of Cancel.
  if (filename.isEmpty())
    return;

  QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));

  // If a scene is already loaded, it has to be freed.
  if (_lcmodel) {
    freeModel();
    freeBBox();
    freeOctree();
    freeDistanceField();
    freeTextures();
  }

  // Load a 3DS model.
  _file = lib3ds_file_load(filename.latin1());
  if (!_file)
    {
      emit statusBarMessage("Unable to open file \"" + filename + "\".", 2000);
      emit message(MsgHandler::MSG_ERROR, name(), "Unable to open file \"" + filename + "\".");
      return;
    }

  if (_mdi && _mdi->config())
    _mdi->config()->io()->setValue("paths/file_3ds", filename);
    
  lib3ds_file_eval(_file, 0);
  
  // Compute a world coordinates representation of the model.
  unsigned tri_count;
  computeModel(&tri_count);

  // Set the bounding boxes.
  _bb_size = _lcmodel->bbMax() - _lcmodel->bbMin();

  // note often we need lots more depth (z) for garment to occupy
  _bb_size[2] *= 3.;

  if (_mdi && _mdi->config()) {
    _mdi->config()->io()->setValue("bbox/size/x", _bb_size[0]);
    _mdi->config()->io()->setValue("bbox/size/y", _bb_size[1]);
    _mdi->config()->io()->setValue("bbox/size/z", _bb_size[2]); 
    _mdi->config()->io()->setValue("bbox/size/x_step", _bb_size[0] / 100);
    _mdi->config()->io()->setValue("bbox/size/y_step", _bb_size[1] / 100);
    _mdi->config()->io()->setValue("bbox/size/z_step", _bb_size[2] / 100);
    _mdi->configUpdated();
  }

  setSceneBoundingBox(_lcmodel->bbMin().address(), _lcmodel->bbMax().address());
  camera()->showEntireScene();
  setSceneRadius(_bb_size.norm()); // Voluntarily too big (factor 2).
  
  std::cout << "Doll BBMin: " << _lcmodel->bbMin() << std::endl;
  std::cout << "Doll BBMax: " << _lcmodel->bbMax() << std::endl;

  manipulatedFrame()->setTranslationAndRotation(sceneCenter(), Quaternion());

  // Prepare GL display lists for both the model and the BBox.
  _model_dl = prepareModel();
  _bb_dl = prepareBBox();

  // 3DS data isn't needed anymore, delete it.
  lib3ds_file_free(_file);

  RepositoryPretreatmentResult::setModelDisplayList(_model_dl);
  Pattern::setClothingModelDisplayList(_model_dl);
  RepositoryPretreatmentResult::setModelBBMin(_lcmodel->bbMin());
  RepositoryPretreatmentResult::setModelBBMax(_lcmodel->bbMax());

  emit statusBarMessage("File \"" + filename + "\" successfully loaded.", 2000);
  emit message(MsgHandler::MSG_NORMAL, name(), "File \"" +
	       filename + "\" successfully loaded.");
  emit message(MsgHandler::MSG_NORMAL, name(), "The model comprises " +
	       QString::number(tri_count) + " triangles.");

  if (!RepositoryPretreatmentResult::result())
    message(MsgHandler::MSG_ERROR, name(), "Could not insert model into the workflow.");
  else {
    message(MsgHandler::MSG_NORMAL, name(), QString("Model sent to ") +
	    RepositoryPretreatmentResult::result()->name() + ".");
    RepositoryPretreatmentResult::pretreatmentUpdated();
  }

  QApplication::restoreOverrideCursor();
}


void Pretreatment::init()
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
  glCullFace(GL_BACK);

  _offscreen = new offscreen::OffScreenArea(offscreen::OffScreenArea::PBUFFER_OFFSCREEN_TYPE,
					    glXGetCurrentContext());
  if (!_offscreen->AllocateOffScreenArea(config::PRETREATMENT_TEXTURE_MAX_SIZE, config::PRETREATMENT_TEXTURE_MAX_SIZE)) {
    message(MsgHandler::MSG_ERROR, name(), "Could not allocate a " +
	    QString::number(config::PRETREATMENT_TEXTURE_MAX_SIZE) + "x" +
	    QString::number(config::PRETREATMENT_TEXTURE_MAX_SIZE) + " offscreen area.");
    delete _offscreen;
    _offscreen = 0;
  }

  setManipulatedFrame(new ManipulatedFrame());
}

void Pretreatment::computeModel(unsigned *tri_count)
{
  if (!_file)
    return;

  _lcmodel = new LCModel();

  if (tri_count)
    *tri_count = 0;
  for (Lib3dsNode* p=_file->nodes; p!=0; p=p->next)
    computeModelRec(p, tri_count);
}

void Pretreatment::computeModelRec(Lib3dsNode *node, unsigned *tri_count)
{
  for (Lib3dsNode* p=node->childs; p!=0; p=p->next)
    computeModelRec(p, tri_count);

  if (node->type == LIB3DS_OBJECT_NODE) {
    if (strcmp(node->name,"$$$DUMMY")==0)
      return;
    
    Lib3dsMesh *mesh=lib3ds_file_mesh_by_name(_file, node->name);
    if (!mesh)
      return;

    Lib3dsObjectData* d = &node->data.object;

    Lib3dsMatrix m;
    lib3ds_matrix_copy(m, mesh->matrix);
    lib3ds_matrix_inv(m);
    
    Lib3dsVector v, vtmp;

    unsigned pts_nb = _lcmodel->points().size();

    Vec3d *pt;
    unsigned p = 0;
    for (; p < mesh->points; ++p) {
      lib3ds_vector_transform(vtmp, m, mesh->pointL[p].pos);
      lib3ds_vector_sub(vtmp, vtmp, d->pivot);
      lib3ds_vector_transform(v, node->matrix, vtmp);

      pt = new Vec3d(v[0], v[1], v[2]);
      _lcmodel->addPoint(pt);
    }

    Triangle *tr;
    p = 0;
    for (; p<mesh->faces; ++p) {
      Lib3dsFace *f=&mesh->faceL[p];
      tr = new Triangle(_lcmodel->points()[pts_nb + f->points[0]],
			_lcmodel->points()[pts_nb + f->points[1]],
			_lcmodel->points()[pts_nb + f->points[2]]);
      _lcmodel->addTriangle(tr);
    }
    if (tri_count)
      *tri_count += p;
  }
}

void Pretreatment::freeModel()
{
  delete _lcmodel;
  _lcmodel = 0;
  if (_model_dl) {
    RepositoryPretreatmentResult::setModelDisplayList(0);
    RepositoryPretreatmentResult::pretreatmentUpdated();
    glDeleteLists(_model_dl, 1);
    _model_dl = 0;
  }
}

int Pretreatment::prepareModel()
{
  int model_dl = glGenLists(1);
  if (!model_dl)
    return 0;
  glNewList(model_dl, GL_COMPILE);
  for (Lib3dsNode* p=_file->nodes; p!=0; p=p->next)
    prepareModelRec(p);
  glEndList();
  
  return model_dl;
}

void Pretreatment::prepareModelRec(Lib3dsNode *node)
{
  for (Lib3dsNode* p=node->childs; p!=0; p=p->next)
    prepareModelRec(p);
    
  if (node->type == LIB3DS_OBJECT_NODE)
    {
      if (strcmp(node->name,"$$$DUMMY")==0)
	return;

      Lib3dsMesh *mesh=lib3ds_file_mesh_by_name(_file, node->name);
      if (!mesh)
	return;

      glPushMatrix();
      Lib3dsObjectData* d = &node->data.object;
      glMultMatrixf(&node->matrix[0][0]);
      glTranslatef(-d->pivot[0], -d->pivot[1], -d->pivot[2]);

      Lib3dsVector *normalL = new Lib3dsVector[3*mesh->faces];

      Lib3dsMatrix M;
      lib3ds_matrix_copy(M, mesh->matrix);
      lib3ds_matrix_inv(M);
      glMultMatrixf(&M[0][0]);

      lib3ds_mesh_calculate_normals(mesh, normalL);

      for (unsigned int p=0; p<mesh->faces; ++p)
	{
	  Lib3dsFace *f=&mesh->faceL[p];
	  Lib3dsMaterial *mat=0;
	  if (f->material[0])
	    mat=lib3ds_file_material_by_name(_file, f->material);

	  if (mat)
	    {
	      static GLfloat a[4]={0,0,0,1};
	      float s;
	      glMaterialfv(GL_FRONT, GL_AMBIENT, a);
	      glMaterialfv(GL_FRONT, GL_DIFFUSE, mat->diffuse);
	      glMaterialfv(GL_FRONT, GL_SPECULAR, mat->specular);
	      s = pow(2, 10.0*mat->shininess);
	      if (s>128.0)
		s=128.0;
	      glMaterialf(GL_FRONT, GL_SHININESS, s);
	    }
	  else
	    {
	      Lib3dsRgba a={0.2, 0.2, 0.2, 1.0};
	      Lib3dsRgba d={0.8, 0.8, 0.8, 1.0};
	      Lib3dsRgba s={0.0, 0.0, 0.0, 1.0};
	      glMaterialfv(GL_FRONT, GL_AMBIENT, a);
	      glMaterialfv(GL_FRONT, GL_DIFFUSE, d);
	      glMaterialfv(GL_FRONT, GL_SPECULAR, s);
	    }

	  glBegin(GL_TRIANGLES);
	  glNormal3fv(f->normal);
	  for (int i=0; i<3; ++i)
	    {
	      glNormal3fv(normalL[3*p+i]);
	      glVertex3fv(mesh->pointL[f->points[i]].pos);
	    }
	  glEnd();
	}

      delete[] normalL;

      glPopMatrix();
    }
}

void Pretreatment::renderModel()
{
  if (!_model_dl)
    return;

  Vec3d size(_bb_size / 2);

  double eq0[4] = {-1, 0, 0, 0};
  eq0[3] = size[0];
  double eq1[4] = {1, 0, 0, 0};
  eq1[3] = size[0];
  double eq2[4] = {0, -1, 0, 0};
  eq2[3] = size[1];
  double eq3[4] = {0, 1, 0, 0};
  eq3[3] = size[1];
  double eq4[4] = {0, 0, -1, 0};
  eq4[3] = size[2];
  double eq5[4] = {0, 0, 1, 0};
  eq5[3] = size[2];

  glPushMatrix();
  glMultMatrixd(manipulatedFrame()->matrix());
  glClipPlane(GL_CLIP_PLANE0, eq0);
  glClipPlane(GL_CLIP_PLANE1, eq1);
  glClipPlane(GL_CLIP_PLANE2, eq2);
  glClipPlane(GL_CLIP_PLANE3, eq3);
  glClipPlane(GL_CLIP_PLANE4, eq4);
  glClipPlane(GL_CLIP_PLANE5, eq5);
  glPopMatrix();

  glEnable(GL_CLIP_PLANE0);
  glEnable(GL_CLIP_PLANE1);
  glEnable(GL_CLIP_PLANE2);
  glEnable(GL_CLIP_PLANE3);
  glEnable(GL_CLIP_PLANE4);
  glEnable(GL_CLIP_PLANE5);

  glCallList(_model_dl);

  glDisable(GL_CLIP_PLANE0);
  glDisable(GL_CLIP_PLANE1);
  glDisable(GL_CLIP_PLANE2);
  glDisable(GL_CLIP_PLANE3);
  glDisable(GL_CLIP_PLANE4);
  glDisable(GL_CLIP_PLANE5);

  glEnable(GL_COLOR_MATERIAL);
  glColor4f(0.4, 0.4, 0.4, 1); // FIXME
  glPushMatrix();
  const float *tr = (camera()->viewDirection() * (sceneRadius() * 0.001)).address();
  glTranslatef(tr[0], tr[1], tr[2]);
  glCallList(_model_dl);

  glPopMatrix();
  glDisable(GL_COLOR_MATERIAL);
}

int Pretreatment::prepareBBox()
{
  int bb_dl = glGenLists(1);
  if (!bb_dl)
    return 0;
  Vec3d bb_min(_bb_size/(-2));
  Vec3d bb_max(_bb_size/2);
  glNewList(bb_dl, GL_COMPILE);
  glDisable(GL_LIGHTING);
  glDepthMask(false);
  glColor4f(0.6, 0.6, 0.6, 0.1); // FIXME
  glutils::cube(bb_min[0], bb_min[1], bb_min[2],
		bb_max[0], bb_max[1], bb_max[2]);
  
    // JDW draw df slice location
  /*
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glColor4f(0.7, 0.3, 0.3, 0.6); 
  glDisable(GL_CULL_FACE);
  glutils::squareY(bb_min[0], bb_min[2],
                   bb_max[0], bb_max[2],
                    bb_min[1] + _df_slice_z*_bb_size[1] );
  glEnable(GL_CULL_FACE);
  */
  
  glDepthMask(true);
  
  glColor4f(0.6, 0.6, 0.6, 0.1); // FIXME
  glutils::fadedBox3D(bb_min[0], bb_min[1], bb_min[2],
		      bb_max[0], bb_max[1], bb_max[2]);
  
  
  
  
  glEnable(GL_LIGHTING);
  glEndList();
  return bb_dl;
}

void Pretreatment::renderBBox()
{
  glPushMatrix();
  glMultMatrixd(manipulatedFrame()->matrix());
  if (_bb_dl)
    glCallList(_bb_dl);
  glPopMatrix();
}

void Pretreatment::freeBBox()
{
  if (_bb_dl) {
    glDeleteLists(_bb_dl, 1);
    _bb_dl = 0;
  }
}

void Pretreatment::renderAxis()
{
  glPushMatrix();
  glMultMatrixd(manipulatedFrame()->matrix());
  drawAxis(sceneRadius()/2);
  glPopMatrix();
}

void Pretreatment::draw()
{
  if (!_lcmodel)
    return;

  if (_disp_model)
    renderModel();

  if (_disp_axis)
    renderAxis();

  if (_disp_df)
    renderDistanceField();

  if (_disp_octree)
    renderOctree();
  
  if (_disp_bb)
    renderBBox();
  
  // jdw to display a slice of the distance field
  /*
  
  if(RepositoryPretreatmentCanvas::distanceField()!=NULL)
  {
  	//DistanceField::initSlice(RepositoryPretreatmentCanvas::distanceField());
  	DistanceField::drawSlice();
  }
  */
  
}

void Pretreatment::computeOctree(unsigned level)
{
  Chronometer chrono;
  chrono.start();
  if (!_lcmodel) {
    message(MsgHandler::MSG_ERROR, name(), "Cannot compute octree: no 3D model loaded yet.");
    return;
  }
  Mat44d mat(manipulatedFrame()->worldMatrix());
  _lcmodel->transform(mat);
  RepositoryPretreatmentResult::setMatrix(mat);
  //  RepositoryPretreatmentResult::pretreatmentUpdated();
  _octree = new Octree(_lcmodel, _bb_size / (-2), _bb_size / 2, level);
  message(MsgHandler::MSG_NORMAL, name(), "Octree computed in " +
	  QString::number(chrono.stop()) + " seconds.");
  emit message(MsgHandler::MSG_NORMAL, name(), "The octree comprises " +
	       QString::number(_octree->trianglesNumber()) + " triangle pointers.");
}

void Pretreatment::freeOctree()
{
  if (!_octree)
    return;
  delete _octree;
  _octree = 0;
  if (_octree_dl) {
    glDeleteLists(_octree_dl, 1);
    _octree_dl = 0;
  }
  message(MsgHandler::MSG_NORMAL, name(), "Octree deleted.");
}

void Pretreatment::prepareOctreeRec(OctreeNode *node, const Vec3d& min, const Vec3d& max)
{
  if (!node || node->level() > _octree->level())
    return;

  glColor4f(0.6, 0.6, 0.6, 0.6); // FIXME
  glutils::box3D(min[0], min[1], min[2],
		 max[0], max[1], max[2]);
  if (node->level() == _octree->level() && !node->triangles().empty()) {
    glColor4f(0.6, 0.6, 0.6, 0.2); // FIXME
    glutils::cube(min[0], min[1], min[2],
		  max[0], max[1], max[2]);
  }

  double size[3] = { utils::max(0.0000001, (max[0] - min[0]) / 2),
		     utils::max(0.0000001, (max[1] - min[1]) / 2),
		     utils::max(0.0000001, (max[2] - min[2]) / 2) };

  Vec3d new_min, new_max;
  for (unsigned location = 0; location < 8; ++location) {
    new_min[0] = !(location & 1) ? min[0] : min[0] + size[0];
    new_min[1] = !(location & 2) ? min[1] : min[1] + size[1];
    new_min[2] = !(location & 4) ? min[2] : min[2] + size[2];
    new_max[0] = !(location & 1) ? max[0] - size[0] : max[0];
    new_max[1] = !(location & 2) ? max[1] - size[1] : max[1];
    new_max[2] = !(location & 4) ? max[2] - size[2] : max[2];
    prepareOctreeRec(node->child(location), new_min, new_max);
  }
}

int Pretreatment::prepareOctree()
{
  if (!_octree)
    return 0;
  int octree_dl = glGenLists(1);
  if (!octree_dl)
    return 0;
  glNewList(octree_dl, GL_COMPILE);
  glPushMatrix();
  glMultMatrixd(manipulatedFrame()->matrix());
  glDisable(GL_LIGHTING);
  glDepthMask(false);
  prepareOctreeRec(_octree->root(), _bb_size / (-2), _bb_size / 2);
  glDepthMask(true);
  glEnable(GL_LIGHTING);
  glPopMatrix();
  glEndList();
  return octree_dl;
}

void Pretreatment::renderOctree()
{
  if (_octree_dl)
    glCallList(_octree_dl);
}

void Pretreatment::computeDistanceField(unsigned x, unsigned y, unsigned z)
{
  Chronometer chrono;
  chrono.start();
  if (!_octree) {
    message(MsgHandler::MSG_ERROR, name(), "Cannot compute distance field: no octree computed yet.");
    return;
  }
  RepositoryPretreatmentCanvas::setDistanceField(new DistanceField(x, y, z));
  
  RepositoryPretreatmentCanvas::distanceField()->compute(_octree);
  
  using std::cout; using std::endl;
  cout << "Euclidean: " << endl;
  
  /*
  std::ofstream df_file("dist_field.vol");
  //df_file.open();
  RepositoryPretreatmentCanvas::distanceField()->write(df_file);
  df_file.close();
  */
  
  /*
  std::ifstream df_file("dist_field.vol");
  //df_file.open();
  RepositoryPretreatmentCanvas::distanceField()->read(df_file,32,32,32);
  df_file.close();
  */
  RepositoryPretreatmentCanvas::distanceField()->dumpStats();
  
  
  message(MsgHandler::MSG_NORMAL, name(), "Distance field computed in " +
      QString::number(chrono.stop()) + " seconds.");
  
  // JDW Testing normalised distance field
  //VecData<Vec3d> *gradField = DistanceField::giveGradientField(RepositoryPretreatmentCanvas::distanceField());
  
  //DistanceField *newdf = giveNormalisedField(RepositoryPretreatmentCanvas::distanceField(),gradField);
  
  //DistanceField *newdf = DistanceField::giveXNormalField(RepositoryPretreatmentCanvas::distanceField(),gradField,0);
  
  //RepositoryPretreatmentCanvas::setDistanceField(newdf);
  
  //cout << "2nd degree normal: " << endl;
  
//  newdf->dumpStats();
  
  
  //  RepositoryPretreatmentCanvas::pretreatmentUpdated();

}

int Pretreatment::prepareDistanceField(unsigned x_min, unsigned x_max,
				       unsigned y_min, unsigned y_max,
				       unsigned z_min, unsigned z_max)
{
  if (!RepositoryPretreatmentCanvas::distanceField())
    return 0;
  int df_dl = glGenLists(1);
  if (!df_dl)
    return 0;
  glNewList(df_dl, GL_COMPILE);
  glPushMatrix();
  glMultMatrixd(manipulatedFrame()->matrix());
  glDisable(GL_LIGHTING);
  glDepthMask(false);
  glPointSize(5);
  glBegin(GL_POINTS);
  double dx = _bb_size[0] / RepositoryPretreatmentCanvas::distanceField()->sizeX();
  double dy = _bb_size[1] / RepositoryPretreatmentCanvas::distanceField()->sizeY();
  double dz = _bb_size[2] / RepositoryPretreatmentCanvas::distanceField()->sizeZ();
  double val;
  Vec3d point;
  point[0] = (-_bb_size[0] + dx) / 2 + x_min * dx;
  for (unsigned i = x_min; i <= x_max; ++i) {
    point[1] = (-_bb_size[1] + dy) / 2 + y_min * dy;
    for (unsigned j = y_min; j <= y_max; ++j) {
      point[2] = (-_bb_size[2] + dz) / 2 + z_min * dz;
      for (unsigned k = z_min; k <= z_max; ++k) {
	val = (RepositoryPretreatmentCanvas::distanceField()->maxDist() - utils::abs(RepositoryPretreatmentCanvas::distanceField()->dist(i, j, k))) / RepositoryPretreatmentCanvas::distanceField()->maxDist();
	val *= val;
	val *= val;
	val *= val;
	val *= val;
	if (RepositoryPretreatmentCanvas::distanceField()->dist(i, j, k) >= 0)
	  glColor4d(0.8, 0.8, 0.8, val);
	else
	  glColor4d(0.8, 0.0, 0.0, val);
	glVertex3d(point[0], point[1], point[2]);
	point[2] += dz;
      }
      point[1] += dy;
    }
    point[0] += dx;
  }
  glEnd();
  
  
  glDepthMask(true);
  glEnable(GL_LIGHTING);
  glPopMatrix();
  glEndList();
  return df_dl;
}

void Pretreatment::renderDistanceField()
{
  if (_df_dl)
    glCallList(_df_dl);
}

void Pretreatment::freeDistanceField()
{
  if (!RepositoryPretreatmentCanvas::distanceField())
    return;
  RepositoryPretreatmentCanvas::freeDistanceField();
  //    RepositoryPretreatmentCanvas::pretreatmentUpdated();
  if (_df_dl) {
    glDeleteLists(_df_dl, 1);
    _df_dl = 0;
  }
  message(MsgHandler::MSG_NORMAL, name(), "Distance field deleted.");
}

void Pretreatment::computeTextures(unsigned texture_size)
{
  if (!_offscreen) {
    message(MsgHandler::MSG_ERROR, name(), "Cannot compute textures: offscreen area could not be allocated.");
    return;
  }

  if (!_lcmodel) {
    message(MsgHandler::MSG_ERROR, name(), "Cannot compute textures: no 3D model loaded yet.");
    return;
  }

  unsigned n = static_cast<unsigned>(log(texture_size) / log(2));
  texture_size = 1 << n;

  if (texture_size > config::PRETREATMENT_TEXTURE_MAX_SIZE) {
    message(MsgHandler::MSG_ERROR, name(), "Cannot compute textures: size too big.");
    return;
  }

  Chronometer chrono;
  chrono.start();

  unsigned char *front_texture = new unsigned char[texture_size * texture_size * 4];
  unsigned char *back_texture = new unsigned char[texture_size * texture_size * 4];

  // Offscreen area grabs the context.
  _offscreen->MakeCurrent();
  glPixelStorei(GL_PACK_ALIGNMENT, 1);

  GLint d;
  glGetIntegerv(GL_DRAW_BUFFER, &d);
  glReadBuffer(d);

  Vec3d size((_bb_size / 2).address());

  // Projection.
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-size[0], size[0], -size[1], size[1], -size[2], size[2]);

  // ModelView.
  glMatrixMode(GL_MODELVIEW);
  //JDW depreciated in v2.0.0 of QGLViewer use inverse().matrix() instead
  //glLoadMatrixd(manipulatedFrame()->worldInverseMatrix());
  glLoadMatrixd(manipulatedFrame()->inverse().matrix());

  // Viewport.
  glViewport(0, 0, texture_size, texture_size);

  // Various GL parameters.
  glShadeModel(GL_SMOOTH);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glDepthFunc(GL_LEQUAL);
  glEnable(GL_DEPTH_TEST);
  glDisable(GL_COLOR_MATERIAL);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);

  // Draw the front texture image.
  glClearColor(0, 0, 0, 0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  if (_model_dl)
    glCallList(_model_dl);

  // Read the front texture image from GL draw buffer.
  glReadPixels(0, 0, texture_size, texture_size, GL_RGBA, GL_UNSIGNED_BYTE, front_texture);
  RepositoryPretreatmentCanvas::setFrontTexture(front_texture, texture_size);

  //     QImage texture(texture_size, texture_size, 32);
  //     texture.setAlphaBuffer(true);
  //     for (int j = 0; j < texture_size; ++j)
  //       for (int i = 0; i < texture_size; ++i) {
  // 	texture.setPixel(i, texture_size - 1 - j,
  // 			 qRgba(front_texture[4*(i+texture_size*j)],
  // 			       front_texture[4*(i+texture_size*j)+1],
  // 			       front_texture[4*(i+texture_size*j)+2],
  // 			       front_texture[4*(i+texture_size*j)+3]));
  //       }
  //     texture.save("texture.png", "PNG");

  // Draw the back texture image.
  double roty[] = {
    -1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, -1, 0,
    0, 0, 0, 1
  };
  glLoadMatrixd(roty); // FIXME
  // JDW Depreciated
  //glMultMatrixd(manipulatedFrame()->worldInverseMatrix());
  glMultMatrixd(manipulatedFrame()->inverse().matrix());

  glClearColor(0, 0, 0, 0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  if (_model_dl)
    glCallList(_model_dl);

  // Read the back texture image from GL draw buffer.
  glReadPixels(0, 0, texture_size, texture_size, GL_RGBA, GL_UNSIGNED_BYTE, back_texture);
  RepositoryPretreatmentCanvas::setBackTexture(back_texture, texture_size);

  Vec2d size2d(size.address());
  RepositoryPretreatmentCanvas::setBBMin(size2d * -1);
  RepositoryPretreatmentCanvas::setBBMax(size2d);
  message(MsgHandler::MSG_NORMAL, name(), "Textures (" +
	  QString::number(texture_size) + "x" + QString::number(texture_size) +
	  ") computed in " + QString::number(chrono.stop()) + " seconds.");

  if (!RepositoryPretreatmentCanvas::canvas())
    message(MsgHandler::MSG_ERROR, name(), "Could not insert textures into the workflow.");
  else {
    message(MsgHandler::MSG_NORMAL, name(), QString("Textures sent to ") +
	    RepositoryPretreatmentCanvas::canvas()->name() + ".");
    RepositoryPretreatmentCanvas::pretreatmentUpdated();
  }

  // Restore the context
  makeCurrent();
}

void Pretreatment::freeTextures()
{
  if (!RepositoryPretreatmentCanvas::frontTexture() && !RepositoryPretreatmentCanvas::backTexture())
    return;
  RepositoryPretreatmentCanvas::freeFrontTexture();
  RepositoryPretreatmentCanvas::freeBackTexture();
  Vec2d null;
  RepositoryPretreatmentCanvas::setBBMin(null);
  RepositoryPretreatmentCanvas::setBBMax(null);
  RepositoryPretreatmentCanvas::pretreatmentUpdated();
  message(MsgHandler::MSG_NORMAL, name(), "Textures deleted.");
}

void Pretreatment::canvasUpdated()
{

}

void Pretreatment::resultUpdated()
{

}
