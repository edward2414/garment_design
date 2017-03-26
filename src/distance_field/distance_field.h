//
//  Filename         : distance_field.h
//  Author           : Emmanuel Turquin
//  Purpose          : Distance Field.
//  Date of creation : 04/20/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  DISTANCE_FIELD_H
# define DISTANCE_FIELD_H

# include "point.h"
# include "voldata.h"
# include "vectypes.h"
# include "vecfield.h"

class Octree;
class OctreeNode;
class LCModel;
class Triangle;



class DistanceField : public VolumeData<short>
{
public:

  DistanceField(unsigned x = 0, unsigned y = 0, unsigned z = 0);
  DistanceField(DistanceField* df); // creates a copy
  ~DistanceField();

  double actualMaxDist(); // JDW required in garment map smoothing. Not sure why maxDist is deliberately factor 2 too small...
  
  void compute(const Octree *octree);

  const LCModel *model() const {
    return _model;
  }

  const Vec3d& min() const {
    return _min;
  }

  const Vec3d& max() const {
    return _max;
  }

  const Vec3d size() const {
    return _max - _min;
  }

  // i,j,k are integer indexes into distance field structure. 
  double dist(unsigned i, unsigned j, unsigned k) const {
    return _max_dist * val(i, j, k) / std::numeric_limits<valuetype>::max();
  }

  // i,j,k are double indexes into distance field structure. 
  // use point_utils::distFromZ if you wish to specify point in world co-ordinates
  double dist(double i, double j, double k) const {
    return _max_dist * valTrilinear(i, j, k) / std::numeric_limits<valuetype>::max();
  }

  Point::Type type(unsigned i, unsigned j) const {
    return _type[i * _y + j];
  }

  double maxDist() const {
    return _max_dist;
  }
  
  // Return the distance gradient vector for a voxel in the distance field
  Vec3d gradCentralDifference(unsigned i, unsigned j, unsigned k) const;
  
  // Return the interpolated distance gradient vector at a point in space
  Vec3d gradTrilinear(double x, double y, double z) const;
  Vec3d bilinearVectorInterpolation(const Vec3d &v1, const Vec3d &v2, const double &alpha) const;
  
  static void drawSlice();
  static void initSlice(DistanceField *df,float percentage);
  
  // return a new field. The gradient field of the supplied field.
  friend DistanceField *giveNormalisedField(DistanceField *df, VecData<Vec3d> *gf);
  static VecData<Vec3d> *giveGradientField(DistanceField *df);
  // numerically differentiate the field wrt the vector field supplied
  static DistanceField *giveDirectionallyDifferentiatedField(DistanceField *df, VecData<Vec3d> *gf);
  static DistanceField *give2NormalField(DistanceField *df, VecData<Vec3d> *gf);
  static DistanceField *giveXNormalField(DistanceField *df, VecData<Vec3d> *gf,  int order);
  // print out stats
  void dumpStats();

private:

  void computeRec(double& dist, const Triangle*& tr, const Vec3d& point, const OctreeNode *octree, const Vec3d& min, const Vec3d& max, unsigned level);

  const LCModel			*_model;
  Vec3d				_min;
  Vec3d				_max;
  double			_max_dist;
  std::vector<Point::Type>	_type;
};

#endif // DISTANCE_FIELD_H
