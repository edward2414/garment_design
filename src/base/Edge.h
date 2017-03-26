#ifndef  POLYMESH_EDGE_H
#define POLYMESH_EDGE_H

#include <vector>
#include "vectypes.h"

class LCModel;

namespace polymesh {

class Face;

class Edge
{
public:
  typedef enum { BOUNDARY, INTERIOR, UNKNOWN } EdgeType; 
  typedef enum { BEGIN, MIDDLE, END } BorderEdgeType; 
  
  Edge(long i1, long i2);
  ~Edge(void) {};
	
 EdgeType _type;
 
 // FIXME Would be neater in a subclass as only used for boundary edges
 polymesh::Edge* _next;
 polymesh::Edge* _previous;
 BorderEdgeType _borderEdgeType;
 bool isDartEdge;
 bool isDartBoundary;
 
 
 
 // The order here is maintained from the vertices in the Face. However only 
 // one Edge is created per pair of vertices. So this direction is only 'correct' for boundary
 // edges which only belong to one Face. Which is fine because we want to traverse the boundary 
 // in one direction later.
 long _vertexIndex[2]; // //FIXME perhaps should be a friend of GarmentAtlas 
 // pointers into the 2d pattern
 Vec3d* _vertexPtr[2];
 

 
 
 void addFace(polymesh::Face*);
 double length();
 int faceCount();
 
 // Normal only assigned for BOUNDARY edges
 Vec3d normal();
 //void normal( Vec3d n ) { _normal = n;};
 polymesh::Face* getFace( int faceIndex ) { return _faceList[ faceIndex ]; };
 //Vec3d* vector(); // FIXME Edge doesn't know where to get it's values from so don't implement yet
 //void recalculateVector(); // For use after resizing pattern
 
 void flip();
 std::vector<polymesh::Face*> faces() { return _faceList; };
 
 
 LCModel *_lcmodel;  // required to recalculate normals as this is an indexed object
 
  private:
   //Vec3d _normal;	
   std::vector<polymesh::Face*> _faceList;
   
};

}
#endif // POLYMESH_EDGE_H
