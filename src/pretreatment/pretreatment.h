//
//  Filename         : pretreatment.h
//  Author           : Emmanuel Turquin
//  Purpose          : A generator of df and textures needed by the algorithm.
//  Date of creation : 04/20/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  PRETREATMENT_H
# define PRETREATMENT_H

# include <lib3ds/file.h>
# include <lib3ds/node.h>
# include "vectypes.h"
# include "base_qglviewer.h"
# include <qapplication.h>

class MDIWindow;
class LCModel;
class OctreeNode;
class Octree;
class DistanceField;

class Pretreatment : public BaseQGLViewer
{
public:

  Pretreatment(MDIWindow* parent=NULL, const char* name=0,
	       const QGLWidget* shareWidget=0, int wflags=0);
  virtual ~Pretreatment();

  void loadFile();
  
  void updateConfig();

  void canvasUpdated();
  void resultUpdated();

private:

  virtual void draw();
  virtual void init();

  int  prepareBBox();
  void renderBBox();
  void freeBBox();

  void renderAxis();

  void computeModel(unsigned *tri_count = 0);
  int  prepareModel();
  void renderModel();
  void freeModel();

  void computeOctree(unsigned level);
  int  prepareOctree();
  void renderOctree();
  void freeOctree();

  void computeDistanceField(unsigned x, unsigned y, unsigned z);
  int  prepareDistanceField(unsigned x_min, unsigned x_max,
			    unsigned y_min, unsigned y_max,
			    unsigned z_min, unsigned z_max);
  void renderDistanceField();
  void freeDistanceField();

  void computeTextures(unsigned texture_size);
  void freeTextures();

  void computeModelRec(Lib3dsNode *node, unsigned *tri_count = 0);
  void prepareModelRec(Lib3dsNode *node);
  void prepareOctreeRec(OctreeNode *node, const Vec3d& min, const Vec3d& max);

  int			_model_dl;
  int			_df_dl;
  int			_octree_dl;
  int			_bb_dl;

  Lib3dsFile		*_file;

  LCModel		*_lcmodel;

  Octree		*_octree;

  Vec3d			_bb_size;

  bool			_disp_model;
  bool			_disp_octree;
  bool			_disp_bb;
  bool			_disp_axis;
  bool			_disp_df;
  float           _df_slice_z; // percentage of z extent to use when drawing df slice

  MDIWindow		*_mdi;
};

#endif // PRETREATMENT_H
