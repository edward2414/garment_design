#include "Vertex.h"
#include <stdlib.h>
#include <vector>
#include "assert.h"

using namespace polymesh;


void Vertex::addEdge(Edge* pEdge) {
  
  if(!pEdge)
    return;
  
  // Due to edge construction procedure there will be many attempts to add a duplicate
  // edge pointer to this vertex. So check by lookup
  
  // FIXME inefficient but doesn't matter as GarmentAtlas construction is once only upon loading
  for(int i = 0; i < _edgeList.size(); i++) {
    if(_edgeList[i] == pEdge) { // compare addresses, i.e. both point to same Edge
      return;
    }
  }
  _edgeList.push_back(pEdge);
}

void Vertex::addFace(polymesh::Face* pFace) {
  if(!pFace)
    return;
  _faceList.push_back(pFace);
}

void Vertex::addCoincidentPoint(polymesh::Vertex* pV)
{
  _coincidentPointList.push_back(pV);
  assert(!(_coincidentPointList.size() >3)); // shouldn't happen when only four atlases can join at one point
}

int Vertex::numberOfCoincidentPoints()
{
  return _coincidentPointList.size();
}

