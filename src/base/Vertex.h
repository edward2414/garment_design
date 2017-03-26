#ifndef  POLYMESH_VERTEX_H
#define POLYMESH_VERTEX_H

#include <vector>
#include "Edge.h"
#include "Face.h"
#include "vectypes.h"


namespace polymesh {
	class Vertex
	{
	public:
   std::vector<polymesh::Edge*> _edgeList;
   std::vector<polymesh::Vertex*> _coincidentPointList;
   
   Vertex(long index) : isDartPoint(false), _index(index) { _point=0;};
   ~Vertex(void) {};
   void addFace(polymesh::Face*);
   void addEdge(polymesh::Edge*);
   void point(Vec3d* p) {_point = p;}; // set the point information
   Vec3d* point() { return _point; };
   void addCoincidentPoint(polymesh::Vertex*);
   int numberOfCoincidentPoints();
   
   bool isDartPoint;
   //std::vector<polymesh::Face*> faces() { return _faceList; };   
   std::vector<polymesh::Face*> _faceList;
	private:
		long _index;
		
		
        Vec3d* _point; // pointer to co-ordinate information
	};
}
#endif //POLYMESH_VERTEX_H
