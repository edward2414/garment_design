//
// C++ Interface: BoundaryEdge
//
// Description: 
//
//
// Author: Jamie Wither <wither@stalactite>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef  POLYMESH_BOUNDARYEDGE_H
#define POLYMESH_BOUNDARYEDGE_H

// Encapsulates an Edge and adds pointer information to track an ordered list of edges.

namespace polymesh {
  class Edge;

  class BoundaryEdge {
    public:
      BoundaryEdge() {};
      ~BoundaryEdge() {};
      polymesh::Edge* _edge;
      polymesh::Edge* _next;
      polymesh::Edge* _previous;

  };

};

#endif //POLYMESH_BOUNDARYEDGE_H