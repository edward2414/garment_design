#include "Edge.h"
#include "Face.h"
#include <iostream>
#include "lcmodel.h"
#include <cassert>


namespace polymesh {



Edge::Edge(long i1, long i2) : isDartEdge(false), isDartBoundary(false)
{
  // Assume indexes are >= 0
  // Assume pFace is not null
  // Assume i1 != i2
  // maintain the ordering of indexes as it gives us the ability to walk around the border of a mesh easily
  
  
    _vertexIndex[0] = i1;
    _vertexIndex[1] = i2;
    
    // this is messy but too many parts of the code needed to make calculations directly on the Edge.
    // Note this is a pointer to the edge in the 2d atlas (pattern / texcoords) not the 3d mesh
    _vertexPtr[0] = 0;
    _vertexPtr[1] = 0;

  _type = UNKNOWN;
  _borderEdgeType = MIDDLE;  // FIXME default, should be in a subclass really
  _next = 0;
  _previous=0;
  _lcmodel=0;
}

double Edge::length()
{
  Vec3d v1( 
      _lcmodel->tex_coords()[ _vertexIndex[1] ]->x() 
      - _lcmodel->tex_coords()[ _vertexIndex[0] ]->x(),
  _lcmodel->tex_coords()[ _vertexIndex[1] ]->y() 
      - _lcmodel->tex_coords()[ _vertexIndex[0] ]->y(),
  0 );
  return v1.norm();
}

Vec3d Edge::normal()
{
  // calculate the normal using the face to supply the other vector...
  assert(_type == polymesh::Edge::BOUNDARY); // FIXME messy
  
  //  calculate the edge normal 
      
  //    using the ordered vertices to form vector v1
  Vec3d v1( 
      _lcmodel->tex_coords()[ _vertexIndex[1] ]->x() 
      - _lcmodel->tex_coords()[ _vertexIndex[0] ]->x(),
  _lcmodel->tex_coords()[ _vertexIndex[1] ]->y() 
      - _lcmodel->tex_coords()[ _vertexIndex[0] ]->y(),
  0 );
      
      //std::cout << "V1 created" << std::endl;
   //   and taking the faceNormal as v2
      // faceNormal is crossProduct of v1 and next vector in the face
  polymesh::Face* pFace = getFace( 0 ); // Only one face belonging to BOUNDARY edges
      
      
      // get a vector from the face which isn't this Edge
      
      
      // FIXME We need the edge vector which follows on... 
  Vec3d vOther;
      
  for(int a=0;a<3;a++) {
        
    if(pFace->_vertex[ a ] != _vertexIndex[1]) { // doesn't follow on
      continue;
    }
    int b=(a+1)%3;
    vOther= Vec3d( 
        _lcmodel->tex_coords()[ pFace->_vertex[b] ]->x() 
        - _lcmodel->tex_coords()[ pFace->_vertex[a] ]->x(),
    _lcmodel->tex_coords()[ pFace->_vertex[b] ]->y() 
        - _lcmodel->tex_coords()[ pFace->_vertex[a] ]->y(),
    0 );
  }

      
      //std::cout << "vOther created" << std::endl;
      
  Vec3d faceNormal = v1^vOther;
  Vec3d edgeNormal = v1^faceNormal;
  edgeNormal.normalize();
      
  // 
  //normal( edgeNormal );
  return edgeNormal;
}



void Edge::addFace(polymesh::Face* pFace) {
  if(!pFace)
    return;
  if(_faceList.size() >= 2) {
    std::cerr << "Warning: Invalid edge. Edge shared by more than 2 faces!" << std::endl;
  }
  _faceList.push_back(pFace);
}

int Edge::faceCount() {
  return _faceList.size();
}

void Edge::flip() 
{
  long _tmp = _vertexIndex[0];
  _vertexIndex[0] = _vertexIndex[1];
  _vertexIndex[1] = _tmp;
}

} // end of namespace polymesh
