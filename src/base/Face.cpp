#include "Face.h"
#include "Edge.h"
#include <iostream>
#include <stack>
#include <vector>
#include "GarmentAtlas.h"



namespace polymesh {



Face::Face(long a, long b, long c)
{
  _vertex[0] = a; 
  _vertex[1] = b; 
  _vertex[2] = c; 
  visited = false;
};

void Face::addEdge(Edge* pEdge) {
  if(!pEdge) {
    std::cerr << "attempt to add null edge pointer to face" << std::endl;
    return;
  }
  if(_edgeList.size() >= 3) {
    std::cerr << "attempt to add more than 3 edges to a face" << std::endl;
    return;
  }
  _edgeList.push_back(pEdge);
  
}

void Face::makeMeshConsistentWith( polymesh::Face* pF, GarmentAtlas* pGarmentAtlas )
{
  std::stack<polymesh::Face*> _consistentButUnvisitedFaceStack;
  _consistentButUnvisitedFaceStack.push(pF);
  
  //std::cout << "Made stack and pushed face" << std::endl;
  
  while(_consistentButUnvisitedFaceStack.size() > 0) {
    makeAdjacentFacesConsistent( &_consistentButUnvisitedFaceStack, pGarmentAtlas );
  }
  

}

void Face::makeAdjacentFacesConsistent( std::stack<polymesh::Face*> *pStack, 
                                               GarmentAtlas *pGarmentAtlas )
{
  polymesh::Face* pF = pStack->top();
  pStack->pop();
  pF->visited = true;
  
  //std::cout << "Pushed top OK" << std::endl;
  
  // FIXME Messy method of ensuring consistency
  // We can't trust the internal ordering of Edges so we need to lookup edges using 
  // lookup on the GarmentAtlas. This Edge then gives us the attached Face. We also have the correct
  // vertex ordering to make the comparison.
  
  // loop over the three vertex pairs
  for( int a = 0; a < 3; a++) 
  {
     int b = (a+1) %3;
     long v1 = pF->_vertex[a];
     long v2 = pF->_vertex[b];
     
     
     polymesh::Edge* pE = pGarmentAtlas->lookup(v1,v2); // Assume always exists
     //std::cout << "E1: " << pE->_vertexIndex[0] << std::endl;
     //std::cout << "E2: " << pE->_vertexIndex[1] << std::endl;
     //std::cout << "SIZE: " << pE->faces().size() <<std::endl;
     //std::cout << "Past lookup" << std::endl;
     // get the adjacent faces, if any
     //std::vector<polymesh::Face*>::const_iterator fi = pE->faces().begin();
     //std::vector<polymesh::Face*>::const_iterator fi_end = pE->faces().end();
     //std::cout << "Edge pointer was OK" << std::endl;
     //for(;fi!= fi_end;fi++)
     for(int ii=0; ii < pE->faceCount(); ii++)
     {
       polymesh::Face* pAdjacentFace = pE->getFace( ii );
       
       //std::cout << "Face verts:" << pF->_vertex[0] << std::endl;
       //std::cout << pF->_vertex[1] << std::endl;
       //std::cout << pF->_vertex[2] << std::endl;
       
       //std::cout << "AFace verts:" << pAdjacentFace->_vertex[0] << std::endl;
       //std::cout << pAdjacentFace->_vertex[1] << std::endl;
       //std::cout << pAdjacentFace->_vertex[2] << std::endl;
       
       if(pAdjacentFace == pF) continue; // Don't process the current face
       if(pAdjacentFace->visited == true) continue;// Don't reprocess visited faces
       //std::cout << "Processing adj face" << std::endl;
       bool flipFace = true; // default is to flip if we don't find a consistent ordering
       // compare ordering in this adjacent face, vertices must be in reverse order to be consistent
       for( int aa = 0; aa < 3; aa++)
       {
         long av1 = pAdjacentFace->_vertex[aa];
         //std::cout << "got first vertex: " << av1 << std::endl;
         if(av1 == v2) {
          int bb = (aa+1) % 3;
          long av2 = pAdjacentFace->_vertex[bb];
          if(av2 == v1) {
            // face direction within this adjacent Face is consistent so no need to flip
            flipFace = false;
            
          }
         }
       }
       // if we didn't find the vertex ordering we expected then flip the face
       if(flipFace)  { 
         pAdjacentFace->flip(); 
         std::cout << "Flipped" << std::endl;
       }
       // now it's consistent we can put it on the stack to be visited
       pStack->push( pAdjacentFace );
       //std::cout << "Pushed adjacent " << std::endl;
     }
     
  } // end of vertex pair iteration
  //std::cout << "Leaving makeConsi..." << std::endl;
  

}

void Face::flip() 
{
  long _tmp = _vertex[0];
  _vertex[0] = _vertex[2];
  _vertex[2] = _tmp;
  
  //std::vector<polymesh::Edge*>::const_iterator ei = _edgeList.begin();
  //std::vector<polymesh::Edge*>::const_iterator ei_end = _edgeList.end();
  //for(; ei != ei_end; ei++) 
  //{
  //  (*ei)->flip();
  //}
}

} // end of namespace polymesh
