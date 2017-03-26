#ifndef  POLYMESH_FACE_H
#define POLYMESH_FACE_H

#include <vector>
#include <stack>
class GarmentAtlas;

namespace polymesh {

class Edge;

class Face
{
public:
  Face(long a, long b, long c);
  ~Face(void) {};
   
  bool visited; // flag used when determining connectivity
  
  void addEdge(Edge*);
  std::vector<polymesh::Edge*> edges() { return _edgeList;};
  long _vertex[3];
  
  // reverse the vertex ordering so as to flip the normal, will also flip associated edges
  void flip();
  void clearEdgeList() { _edgeList.clear(); };
  
  // Passed a reference Face make all connected Faces use the same ordering
  static void makeMeshConsistentWith( polymesh::Face* , GarmentAtlas*  );
  static void makeAdjacentFacesConsistent( std::stack<polymesh::Face*>*,                                                                             GarmentAtlas* );

  private:
	
    std::vector<polymesh::Edge*> _edgeList;
 
 
};
}

#endif // POLYMESH_FACE_H

