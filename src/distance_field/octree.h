//
//  Filename         : octree.h
//  Author           : Emmanuel Turquin
//  Purpose          : Octree containing triangles (use lib3ds).
//  Date of creation : 05/18/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  OCTREE_H
# define OCTREE_H

# include <vector>
# include "vectypes.h"

class Triangle;
class LCModel;

//////////////////// OctreeNode ////////////////////

class OctreeNode
{
public:

  typedef std::vector<const Triangle*>	TriangleList;

  OctreeNode(OctreeNode *parent, unsigned location);
  ~OctreeNode();

  OctreeNode *child(unsigned location) {
    return _children[location];
  }

  const OctreeNode *child(unsigned location) const {
    return _children[location];
  }

  OctreeNode *createChild(unsigned location);

  unsigned level() const {
    return _level;
  }

  unsigned location() const {
    return _location;
  }

  void addTriangle(const Triangle *t) {
    if (!t)
      return;
    _triangles.push_back(t);
  }

  const TriangleList& triangles() const {
    return _triangles;
  }

private:

  OctreeNode		*_parent;
  OctreeNode		*_children[8];
  unsigned		_level;
  unsigned		_location;
  TriangleList		_triangles;
};

//////////////////// Octree ////////////////////

class Octree
{
public:

  Octree(const LCModel *model, const Vec3d& min, const Vec3d& max, unsigned subdiv_level);
  ~Octree();

  
  OctreeNode *root() {
    return _root;
  }

  const OctreeNode *root() const {
    return _root;
  }

  unsigned level() const {
    return _level;
  }

  const LCModel *model() const {
    return _model;
  }

  const Vec3d& min() const {
    return _min;
  }

  const Vec3d& max() const {
    return _max;
  }

  unsigned trianglesNumber() const {
    return _triangles_nb;
  }

private:

  void addTriangle(const Triangle* tr);

  void addTriangleRec(const Triangle *tr, OctreeNode *node, const Vec3d& min, const Vec3d& max);

  OctreeNode		*_root;
  unsigned		_level;
  unsigned		_triangles_nb;
  const LCModel		*_model;
  Vec3d			_min;
  Vec3d			_max;
};

#endif // OCTREE_H
