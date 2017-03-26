//
//  Filename         : octree.cpp
//  Author           : Emmanuel Turquin
//  Purpose          : Octree containing triangles (use lib3ds).
//  Date of creation : 05/18/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include <cmath>
#include "utils.h"
#include "lcmodel.h"
#include "octree.h"

//////////////////// OctreeNode ////////////////////

OctreeNode::OctreeNode(OctreeNode *parent, unsigned location)
{
  if (!parent)
    _level = 0;
  else
    _level = parent->_level + 1;
  for (unsigned i = 0; i < 8; ++i)
    _children[i] = 0;
  _location = location;
}

OctreeNode::~OctreeNode()
{
  for (unsigned i = 0; i < 8; ++i)
    delete _children[i];
}

OctreeNode *OctreeNode::createChild(unsigned location)
{
  if (!_children[location])
    _children[location] = new OctreeNode(this, location);
  return _children[location];
}

//////////////////// Octree ////////////////////

Octree::Octree(const LCModel *model, const Vec3d& min, const Vec3d& max, unsigned subdiv_level)
{
  _model = model;
  _min = min;
  _max = max;
  _level = subdiv_level;
  _triangles_nb = 0;
  _root = new OctreeNode(0, 0);
  LCModel::TriangleList::const_iterator it = _model->triangles().begin();
  LCModel::TriangleList::const_iterator it_end = _model->triangles().end();
  for (; it != it_end; ++it)
    addTriangle(*it);
}

Octree::~Octree()
{
  delete _root;
}

void Octree::addTriangle(const Triangle* tr)
{
  if (!_model || !tr)
    return;

  addTriangleRec(tr, _root, _min, _max);
}

// We presume that (tr != 0) and (node != 0), since they have been tested before.
void Octree::addTriangleRec(const Triangle *tr, OctreeNode *node, const Vec3d& min, const Vec3d& max)
{
  if (node->level() == _level) {
    node->addTriangle(tr);
    ++_triangles_nb;
    return;
  }

  Vec3d center((min + max) / 2);
  Vec3d size((max - center) / 2);
  Vec3d new_min, new_max;
  for (unsigned location = 0; location < 8; ++location) {
    new_min[0] = location & 1 ? center[0] : min[0];
    new_min[1] = location & 2 ? center[1] : min[1];
    new_min[2] = location & 4 ? center[2] : min[2];
    new_max[0] = location & 1 ? max[0] : center[0];
    new_max[1] = location & 2 ? max[1] : center[1];
    new_max[2] = location & 4 ? max[2] : center[2];
    if (utils::overlapTriangleBox((new_min + new_max) / 2, size, *tr)) {
      addTriangleRec(tr, node->createChild(location), new_min, new_max);
    }
  }
}
