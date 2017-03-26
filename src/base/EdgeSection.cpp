//
// C++ Implementation: EdgeSection
//
// Description: 
//
//
// Author: Jamie Wither <wither@stalactite>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "EdgeSection.h"
#include "vectypes.h"
#include "Edge.h"
#include "lcmodel.h"

using namespace std;
using namespace polymesh;

long EdgeSection::closestTCIndexTo(double x, double y) const
{
  
  double minDist = 9999999.9;
  vector<Edge*>::const_iterator ei;
  long closestIndex;
  for(ei = _edges.begin();
      ei != _edges.end();
      ei++)
  {
    Edge* pEdge = *ei;
    
    Vec3d* p = pEdge->_vertexPtr[0];
    double distSq = (p->x() - x) * (p->x() - x) + (p->y() - y) * (p->y() - y);
    if(distSq < minDist) 
    {
      closestIndex = pEdge->_vertexIndex[0];
      minDist = distSq; 
    }
    
    p = pEdge->_vertexPtr[1];
    distSq = (p->x() - x) * (p->x() - x) + (p->y() - y) * (p->y() - y);
    if(distSq < minDist)
    {
      closestIndex = pEdge->_vertexIndex[1];
      minDist = distSq; 
    }
    
  }
  return closestIndex;

}

EdgeSection::EdgeSection() : _orientation(UNKNOWN), _classification(UNCLASSIFIED) {
  next = 0;
  previous = 0;
}

void EdgeSection::merge(EdgeSection* pEStoMerge) 
{
  std::vector<polymesh::Edge*>::iterator ei;
  for(ei = pEStoMerge->_edges.begin();
      ei!=pEStoMerge->_edges.end();
      ei++)
  {
    _edges.push_back( *ei );
  }
}

void EdgeSection::recalculateAttributes(LCModel *_lcmodel) 
{
  
    EdgeSection* pEdgeSection = this;
    
     
    Vec3d accumulatePosition(0.0,0.0,0.0);
    Vec3d accumulateNormal(0.0,0.0,0.0);
    double accumulateLength=0.0;
    vector<Edge*>::const_iterator ei, ei_end;
    ei= _edges.begin();
    ei_end = _edges.end();
    for(;ei!=ei_end;ei++) // each edge in this section
    {
      polymesh::Edge* pEdge = *ei;
      Vec3d centrePoint = *(_lcmodel->tex_coords()[ pEdge->_vertexIndex[0] ]);
      centrePoint +=      *(_lcmodel->tex_coords()[ pEdge->_vertexIndex[1] ]);
      centrePoint /= 2.0;
      
      Vec3d edgeVec = *(_lcmodel->tex_coords()[ pEdge->_vertexIndex[0] ]) - 
          *(_lcmodel->tex_coords()[ pEdge->_vertexIndex[1] ]);
      accumulateLength += edgeVec.norm();
      
      accumulatePosition += centrePoint;
      accumulateNormal += pEdge->normal();
    }
    
    
    Vec3d aN = accumulateNormal.normalize();
    pEdgeSection->averageNormal = aN;
    pEdgeSection->averagePos = accumulatePosition / (double)_edges.size();
    pEdgeSection->length = accumulateLength;
    
    // Determine orientation of edge
    if(aN.y() > 0) 
    {
      
      if( aN.y() > fabs(aN.x() ) )
        pEdgeSection->_orientation = EdgeSection::UP;
      else if( aN.x() < 0 )
        pEdgeSection->_orientation = EdgeSection::LEFT;
      else 
        pEdgeSection->_orientation = EdgeSection::RIGHT;
    }
    else
    {
      if( fabs(aN.y()) > fabs(aN.x() ) )
        pEdgeSection->_orientation = EdgeSection::DOWN;
      else if( aN.x() < 0 )
        pEdgeSection->_orientation = EdgeSection::LEFT;
      else 
        pEdgeSection->_orientation = EdgeSection::RIGHT;
    }
    
    
  
}