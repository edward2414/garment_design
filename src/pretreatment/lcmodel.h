//
//  Filename         : lcmodel.h
//  Author           : Emmanuel Turquin
//  Purpose          : 3D model in world coordinates.
//  Date of creation : 05/24/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  LCMODEL_H
# define LCMODEL_H

# include <vector>
# include <limits>
#include <math.h>
# include "utils.h"
# include "vectypes.h"
# include "mattypes.h"
# include "triangle.h"
# include "indexed_triangle.h"


#include "GarmentAtlas.h"

class LCModel
{
public:

  typedef std::vector<Vec3d*>		PointList;
  typedef std::vector<const Triangle*>	TriangleList;
  typedef std::vector<IndexedTriangle*>    IndexedTriangleList;
  typedef std::vector<Vec3d*>  TextureCoordList; // only x and y valid. z is zero in all cases
  typedef std::vector<GarmentAtlas*> AtlasList; // stores the various atlases in the model - will be changed for a full hierarchy later
  
  LCModel(const Mat44d& m = Mat44d::identity()) {
    _mat = m;
    _bb_min[0] = std::numeric_limits<Vec3d::value_type>::max();
    _bb_min[1] = std::numeric_limits<Vec3d::value_type>::max();
    _bb_min[2] = std::numeric_limits<Vec3d::value_type>::max();
    _bb_max[0] = -std::numeric_limits<Vec3d::value_type>::max();
    _bb_max[1] = -std::numeric_limits<Vec3d::value_type>::max();
    _bb_max[2] = -std::numeric_limits<Vec3d::value_type>::max();
    _tc_bb_min[0] = std::numeric_limits<Vec3d::value_type>::max();
    _tc_bb_min[1] = std::numeric_limits<Vec3d::value_type>::max();
    _tc_bb_min[2] = std::numeric_limits<Vec3d::value_type>::max();
    _tc_bb_max[0] = -std::numeric_limits<Vec3d::value_type>::max();
    _tc_bb_max[1] = -std::numeric_limits<Vec3d::value_type>::max();
    _tc_bb_max[2] = -std::numeric_limits<Vec3d::value_type>::max();
  }

  ~LCModel() {
    TriangleList::iterator t = _tl.begin();
    TriangleList::iterator t_end = _tl.end();
    for (; t != t_end; ++t)
      delete *t;
    
    PointList::iterator p = _pl.begin();
    PointList::iterator p_end = _pl.end();
    for (; p != p_end; ++p)
      delete *p;
    
    TextureCoordList::iterator tc = _tcl.begin();
    TextureCoordList::iterator tc_end = _tcl.end();
    for (; tc != tc_end; ++tc)
      delete *tc;
    
    IndexedTriangleList::iterator itl = _itl.begin();
    IndexedTriangleList::iterator itl_end = _itl.end();
    for (; itl != itl_end; ++itl)
      delete *itl;
    
  }

  void translateAtlas(Vec3d v)
  {
    std::vector<Vec3d*>::iterator vi, vi_end;
    vi = _tcl.begin();
    vi_end = _tcl.end();
    for(;vi!=vi_end;vi++)
    {
      **vi += v;
    }
  }
  
  void correctAtlasRotation() {
    // Chose a 2D (X,Y) vector from the 3D mesh (mesh has the 'correct' vertical orientation)
    // Then find the corresponding vector in the 2D pattern and determine the 
    // rotation required to align the pattern to the mesh and apply it.
    
    // Determine whether this is a front or back panel by using face normals. This will
    // determine which way to rotate the pattern
    
    // Find face normal for mesh
    Vec3d meshFaceNormal = Vec3d();
    { // block to avoid naming clashes
      Vec3d* p1 = points()[ indexedTriangles()[0]->indexA() ];
      Vec3d* p2 = points()[ indexedTriangles()[0]->indexB() ];
      Vec3d* p3 = points()[ indexedTriangles()[0]->indexC() ];
      
      // Two vectors in the plane of the face
      Vec3d v1( p2->x() - p1->x(),
                p2->y() - p1->y(),
                p2->z() - p1->z());
      Vec3d v2( p3->x() - p2->x(),
                p3->y() - p2->y(),
                p3->z() - p2->z());
      
      // cross product to give faceNormal
      meshFaceNormal = v1^v2;
      std::cout << "mNormal: " << meshFaceNormal << std::endl;
      
    }
    
        // Find face normal for pattern
    Vec3d patternFaceNormal = Vec3d();
    { // block to avoid naming clashes
      Vec3d* p1 = tex_coords()[ indexedTriangles()[0]->indexA() ];
      Vec3d* p2 = tex_coords()[ indexedTriangles()[0]->indexB() ];
      Vec3d* p3 = tex_coords()[ indexedTriangles()[0]->indexC() ];
      
      // Two vectors in the plane of the face
      Vec3d v1( p2->x() - p1->x(),
                p2->y() - p1->y(),
                p2->z() - p1->z());
      Vec3d v2( p3->x() - p2->x(),
                p3->y() - p2->y(),
                p3->z() - p2->z());
      
      // cross product to give normal
      patternFaceNormal = v1^v2;
    }
    
    // If Z sign is different then flip the X comp of vector for rotation correction.
    // FIXME Determining rotation correction this way could be too simple if assumptions change
    double flipX = 1.0;
    { // avoid name clash

      if( (patternFaceNormal.z() * meshFaceNormal.z()) < 0.0 ) { // angle greater than 90 deg, so facing away
        flipX = -1.0;
      }
      
    }
    
    // Points from 3d mesh
    // FIXME - need a robust method of picking points to form a long vertical vector close
    // to the centre along the Y axis direction. This will do for now as it's not a core problem
    unsigned ia = (points().size() - 1);
    unsigned ib = ia - (points().size() / 3);
    Vec3d* pa1 = points()[ indexedTriangles()[ia]->indexA() ];
    Vec3d* pa2 = points()[ indexedTriangles()[ib]->indexA() ];
    
    // Points from 2d mesh
    Vec3d* pb1 = tex_coords()[ indexedTriangles()[ia]->indexA() ];
    Vec3d* pb2 = tex_coords()[ indexedTriangles()[ib]->indexA() ];
    
    // Create vectors in the XY plane from these
    Vec3d vMesh(flipX * ( pa2->x() - pa1->x() ), // mirror the x component if required as this is the back of the clothing.
                pa2->y() - pa1->y(),
                0);
    std::cout <<"Flip is: " << flipX << std::endl;
    Vec3d vPattern(pb2->x() - pb1->x(),
                   pb2->y() - pb1->y(),
                   0);

    double dp = utils::dotProduct(vMesh.x(),vMesh.y(),vPattern.x(),vPattern.y());
    double productOfNorms = vMesh.norm() * vPattern.norm();
    double theta = acos(dp / productOfNorms); // FIXME check for division by zero - could happen
    
    
    
    // mat[M][N] M - row, N = col
    double mat[4][4]; // rotation matrix in homogeneous co-ords. - use to initialise Mat44d
    for(int i = 0; i < 4; i++)
      for(int j = 0; j < 4; j++)
        mat[i][j] = 0.0;
    
    mat[2][2] = 1.0;
    mat[3][3] = 1.0;
    mat[0][0] = cos(theta);
    mat[0][1] = sin(theta);
    mat[1][0] = -1.0 * mat[0][1];
    mat[1][1] = mat[0][0];
    
    Mat44d rotationMatrixAroundZ = Mat44d(mat);
    transformAtlas(rotationMatrixAroundZ);
  }
  
  
  void correctAtlasRotationUser(Vec3d pointA, Vec3d pointB) {
    // Users has chosen the central axis of the pattern.
    // Determine the rotation (about Z) required to align the pattern to Y-axis and apply it.
    
    // Determine whether this is a front or back panel by using face normals. This will
    // determine which way to rotate the pattern
    
    // Find face normal for mesh
    /*
    Vec3d meshFaceNormal = Vec3d();
    { // block to avoid naming clashes
      Vec3d* p1 = points()[ indexedTriangles()[0]->indexA() ];
      Vec3d* p2 = points()[ indexedTriangles()[0]->indexB() ];
      Vec3d* p3 = points()[ indexedTriangles()[0]->indexC() ];
      
      // Two vectors in the plane of the face
      Vec3d v1( p2->x() - p1->x(),
                p2->y() - p1->y(),
                p2->z() - p1->z());
      Vec3d v2( p3->x() - p2->x(),
                p3->y() - p2->y(),
                p3->z() - p2->z());
      
      // cross product to give faceNormal
      meshFaceNormal = v1^v2;
      std::cout << "mNormal: " << meshFaceNormal << std::endl;
      
    }
    */
    
        // Find face normal for pattern
    /*
    Vec3d patternFaceNormal = Vec3d();
    { // block to avoid naming clashes
      Vec3d* p1 = tex_coords()[ indexedTriangles()[0]->indexA() ];
      Vec3d* p2 = tex_coords()[ indexedTriangles()[0]->indexB() ];
      Vec3d* p3 = tex_coords()[ indexedTriangles()[0]->indexC() ];
      
      // Two vectors in the plane of the face
      Vec3d v1( p2->x() - p1->x(),
                p2->y() - p1->y(),
                p2->z() - p1->z());
      Vec3d v2( p3->x() - p2->x(),
                p3->y() - p2->y(),
                p3->z() - p2->z());
      
      // cross product to give normal
      patternFaceNormal = v1^v2;
    }
    */
    /*
    // If Z sign is different then flip the X comp of vector for rotation correction.
    // FIXME Determining rotation correction this way could be too simple if assumptions change
    double flipX = 1.0;
    { // avoid name clash

      if( (patternFaceNormal.z() * meshFaceNormal.z()) < 0.0 ) { // angle greater than 90 deg, so facing away
        flipX = -1.0;
      }
      
    }
    */
    /*
    // Points from 3d mesh
    // FIXME - need a robust method of picking points to form a long vertical vector close
    // to the centre along the Y axis direction. This will do for now as it's not a core problem
    unsigned ia = (points().size() - 1);
    unsigned ib = ia - (points().size() / 3);
    Vec3d* pa1 = points()[ indexedTriangles()[ia]->indexA() ];
    Vec3d* pa2 = points()[ indexedTriangles()[ib]->indexA() ];
    
    // Points from 2d mesh
    Vec3d* pb1 = tex_coords()[ indexedTriangles()[ia]->indexA() ];
    Vec3d* pb2 = tex_coords()[ indexedTriangles()[ib]->indexA() ];
    */
    // Create vectors in the XY plane from these
    /*
    Vec3d vMesh(flipX * ( pa2->x() - pa1->x() ), // mirror the x component if required as this is the back of the clothing.
                pa2->y() - pa1->y(),
                0);
    std::cout <<"Flip is: " << flipX << std::endl;
    Vec3d vPattern(pb2->x() - pb1->x(),
                   pb2->y() - pb1->y(),
                   0);
    */
    Vec3d verticalAxis(0.0,1.0,0.0);
    Vec3d centralAxis = pointA - pointB;
    
    std::cout << "Central axis is: " << centralAxis << std::endl;

    double dp = utils::dotProduct(verticalAxis.x(),verticalAxis.y(),centralAxis.x(),centralAxis.y());
    double productOfNorms = verticalAxis.norm() * centralAxis.norm();
    double theta = acos(dp / productOfNorms); // FIXME check for division by zero - could happen
    
    std::cout << "Theta is: " << theta << std::endl;
    // Determine cross product for sign of theta
    Vec3d cross = centralAxis^verticalAxis;
    std::cout << "Cross is: " << cross << std::endl;
    if(cross.z() > 0) theta *= -1.0;
    
    // We need to rotate about the centre of the pattern. So find transform which does this
    //tcbbMax() - tcbbMin()
    
    // mat[M][N] M - row, N = col
    double mat[4][4]; // rotation matrix in homogeneous co-ords. - use to initialise Mat44d
    for(int i = 0; i < 4; i++)
      for(int j = 0; j < 4; j++)
        mat[i][j] = 0.0;
    
    mat[2][2] = 1.0;
    mat[3][3] = 1.0;
    mat[0][0] = cos(theta);
    mat[0][1] = sin(theta);
    mat[1][0] = -1.0 * mat[0][1];
    mat[1][1] = mat[0][0];
    
    Mat44d rotationMatrixAroundZ = Mat44d(mat);
    transformAtlas(rotationMatrixAroundZ);
    
    recalculateTCBB();
    
  }
  
  void correctAtlasRotationUser(Vec3d axis) {

    Vec3d verticalAxis(0.0,1.0,0.0);
    Vec3d centralAxis = axis;
    
    std::cout << "Central axis supplied is: " << centralAxis << std::endl;

    double dp = utils::dotProduct(verticalAxis.x(),verticalAxis.y(),centralAxis.x(),centralAxis.y());
    double productOfNorms = verticalAxis.norm() * centralAxis.norm();
    double theta = acos(dp / productOfNorms); // FIXME check for division by zero - could happen
    
    std::cout << "Theta is: " << theta << std::endl;
    // Determine cross product for sign of theta
    Vec3d cross = centralAxis^verticalAxis;
    std::cout << "Cross is: " << cross << std::endl;
    if(cross.z() > 0) theta *= -1.0;
    
    // We need to rotate about the centre of the pattern. So find transform which does this
    //tcbbMax() - tcbbMin()
    
    // mat[M][N] M - row, N = col
    double mat[4][4]; // rotation matrix in homogeneous co-ords. - use to initialise Mat44d
    for(int i = 0; i < 4; i++)
      for(int j = 0; j < 4; j++)
        mat[i][j] = 0.0;
    
    mat[2][2] = 1.0;
    mat[3][3] = 1.0;
    mat[0][0] = cos(theta);
    mat[0][1] = sin(theta);
    mat[1][0] = -1.0 * mat[0][1];
    mat[1][1] = mat[0][0];
    
    Mat44d rotationMatrixAroundZ = Mat44d(mat);
    transformAtlas(rotationMatrixAroundZ);
    
    recalculateTCBB();
  }
  
  const PointList& points() const {
    return _pl;
  }
  
  const TextureCoordList& tex_coords() const {
    return _tcl;
  }

  void addPoint(Vec3d *p) {
    if (!p)
      return;
    _pl.push_back(p);
    _bb_min.min(*p);
    _bb_max.max(*p);
  }
  
  void addTextureCoord(Vec3d *p) {
    if (!p)
      return;
    _tcl.push_back(p);
    _tc_bb_min.min(*p);
    _tc_bb_max.max(*p);
  }
  
  void recalculateTCBB() {
    
    _tc_bb_min[0] = std::numeric_limits<Vec3d::value_type>::max();
    _tc_bb_min[1] = std::numeric_limits<Vec3d::value_type>::max();
    _tc_bb_min[2] = std::numeric_limits<Vec3d::value_type>::max();
    _tc_bb_max[0] = -std::numeric_limits<Vec3d::value_type>::max();
    _tc_bb_max[1] = -std::numeric_limits<Vec3d::value_type>::max();
    _tc_bb_max[2] = -std::numeric_limits<Vec3d::value_type>::max();
    
    std::vector<Vec3d*>::const_iterator vi = _tcl.begin();
    std::vector<Vec3d*>::const_iterator vi_end = _tcl.end();
    for(;vi != vi_end;vi++) 
    {
      _tc_bb_min.min(**vi);
      _tc_bb_max.max(**vi);
    }
  }

  const TriangleList& triangles() const {
    return _tl;
  }
  
  IndexedTriangleList& indexedTriangles()  {
    return _itl;
  }

  void addTriangle(const Triangle *t) {
    if (!t)
      return;
    _tl.push_back(t);
  }
  
  void addIndexedTriangle( IndexedTriangle *it) {
    if (!it)
      return;
    _itl.push_back(it);
  }

  const Vec3d& bbMin() const {
    return _bb_min;
  }

  Vec3d& bbMin() {
    return _bb_min;
  }
  
  Vec3d& tcbbMin() {
    return _tc_bb_min;
  }

  void setBBMin(const Vec3d& v) {
    _bb_min = v;
  }

  const Vec3d& bbMax() const {
    return _bb_max;
  }

  Vec3d& bbMax() {
    return _bb_max;
  }
  
  Vec3d& tcbbMax() {
    return _tc_bb_max;
  }

  void setBBMax(const Vec3d& v) {
    _bb_max = v;
  }

  const Mat44d& matrix() const {
    return _mat;
  }

  Mat44d& matrix() {
    return _mat;
  }
  


  void transform(const Mat44d& m = Mat44d::identity()) {
    double m_inv[4][4];
    m.getArray(m_inv);
    utils::invertMatrix(m_inv);
    Mat44d trans = (_mat * Mat44d(m_inv)).transpose();
    PointList::iterator p = _pl.begin();
    PointList::iterator p_end = _pl.end();
    HVec3d ptmp;
    _bb_min[0] = std::numeric_limits<Vec3d::value_type>::max();
    _bb_min[1] = std::numeric_limits<Vec3d::value_type>::max();
    _bb_min[2] = std::numeric_limits<Vec3d::value_type>::max();
    _bb_max[0] = -std::numeric_limits<Vec3d::value_type>::max();
    _bb_max[1] = -std::numeric_limits<Vec3d::value_type>::max();
    _bb_max[2] = -std::numeric_limits<Vec3d::value_type>::max();
    for (; p != p_end; ++p) {
      ptmp = **p;
      **p = HVec3d(trans * ptmp);
      _bb_min.min(**p);
      _bb_max.max(**p);
    }
    _mat = m;
  }

  void transformAtlas(const Mat44d& m = Mat44d::identity()) {
    double m_inv[4][4];
    m.getArray(m_inv);
    utils::invertMatrix(m_inv);
    Mat44d trans = (_mat * Mat44d(m_inv)).transpose();
    TextureCoordList::iterator p = _tcl.begin();
    TextureCoordList::iterator p_end = _tcl.end();
    HVec3d ptmp;
    _tc_bb_min[0] = std::numeric_limits<Vec3d::value_type>::max();
    _tc_bb_min[1] = std::numeric_limits<Vec3d::value_type>::max();
    _tc_bb_min[2] = std::numeric_limits<Vec3d::value_type>::max();
    _tc_bb_max[0] = -std::numeric_limits<Vec3d::value_type>::max();
    _tc_bb_max[1] = -std::numeric_limits<Vec3d::value_type>::max();
    _tc_bb_max[2] = -std::numeric_limits<Vec3d::value_type>::max();
    for (; p != p_end; ++p) {
      ptmp = **p;
      **p = HVec3d(trans * ptmp);
      _tc_bb_min.min(**p);
      _tc_bb_max.max(**p);
    }
    _mat = m;
  }
private:

  PointList	_pl;
  TriangleList	_tl;
  IndexedTriangleList  _itl;
  TextureCoordList  _tcl;

  Vec3d		_bb_min;
  Vec3d		_bb_max;
  Vec3d     _tc_bb_min;
  Vec3d     _tc_bb_max;
  Mat44d	_mat;
};

#endif // LCMODEL_H


