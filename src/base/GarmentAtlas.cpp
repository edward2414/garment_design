//
// C++ Implementation: GarmentAtlas
//
// Description: 
//
//
// Author: Jamie Wither <wither@stalactite>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "GarmentAtlas.h"
#include "GarmentSection.h"
#include "Garment.h"
#include <iostream>
#include "lcmodel.h"
#include "indexed_triangle.h"
#include "Edge.h"
#include "Face.h"
#include "Vertex.h"
#include <stdlib.h>
#include <qstring.h>
#include <utility>
#include "vectypes.h"
#include "vecmat.h"
#include <math.h>
#include "EdgeSection.h"
#include <set>
#include <algorithm>
#include "Parabola.h"
#include <cassert>
#include <utility>
#include "garment_parameters.h"
#include <string>
#include "PushBackParams.h"

Vec3d topTorsoLeft;
Vec3d middleTorsoLeft;
Vec3d topTorsoRight;
Vec3d middleTorsoRight;
Vec3d rightCuffTop;
Vec3d rightCuffBottom;
Vec3d leftCuffTop;
Vec3d leftCuffBottom;

  struct pointDistance { 
    double dist;
    long index;
    Vec3d* pPoint;
  };
  
    
  struct closestPointFunctor
  {
    bool operator()(const pointDistance pd1, const pointDistance pd2)
    {
      return( pd1.dist < pd2.dist );
    }
  };

struct esPtrYLess
{
  bool operator()(const EdgeSection* pEs1, const EdgeSection* pEs2)
  {
    return( pEs1->averagePos.y() < pEs2->averagePos.y() );
  }
};

struct esPtrXLess
{
  bool operator()(const EdgeSection* pEs1, const EdgeSection* pEs2)
  {
    return( pEs1->averagePos.x() < pEs2->averagePos.x() );
  }
};

// Allows ordering of Edges in a container. Assumes all Edges have same orientation.
struct edgePtrYLess
{
  bool operator()(const polymesh::Edge* pE1, const polymesh::Edge* pE2)
  {
    return( pE1->_vertexPtr[0]->y() < pE2->_vertexPtr[0]->y() );
  }
};

 GarmentAtlas::~GarmentAtlas() {
      // release Vertex objects
  std::vector<polymesh::Vertex*>::const_iterator iv = _vertexList.begin();
  std::vector<polymesh::Vertex*>::const_iterator iv_end = _vertexList.end();
  for(; iv != iv_end; iv++) 
  {
    delete (*iv);
  }
      
      // release Edges
  std::vector<polymesh::Edge*>::const_iterator eiv = _edgeList.begin();
  std::vector<polymesh::Edge*>::const_iterator eiv_end = _edgeList.end();
  for(; eiv != eiv_end; eiv++) 
  {
    delete (*eiv);
  }
      
      // release Faces
  std::vector<polymesh::Face*>::const_iterator fiv = _faceList.begin();
  std::vector<polymesh::Face*>::const_iterator fiv_end = _faceList.end();
  for(; fiv != fiv_end; fiv++) 
  {
    delete (*fiv);
  }
      
      // release model
  if(_lcmodel)
    delete (_lcmodel);
            
  if(_listOfSeamEdgeSections.size())
  {
    std::vector<EdgeSection*>::iterator esi = _listOfSeamEdgeSections.begin();
    for(;esi!=_listOfSeamEdgeSections.end();esi++)
    {
      delete (*esi);
    }
  }
      
          
      
}

void GarmentAtlas::linkBoundaryEdges()
{
  assert(_unsortedBoundaryEdges.size() > 0);
  //for each BOUNDARY edge
  //  determine which edges are associated with begin and end vertex (near and far end of edge)
  //    and are also BOUNDARY
  //    and are not this edge
  // point to as the 'next' or 'previous' edge

  std::vector<polymesh::Edge*>::const_iterator ei,ei_end;
  ei = _unsortedBoundaryEdges.begin();
  ei_end = _unsortedBoundaryEdges.end();
  for(;ei!=ei_end;ei++)
  {
    polymesh::Edge* pEdge = *ei;
    long endIndex = pEdge->_vertexIndex[1]; // End point of edge;
    polymesh::Vertex* pVertex = _vertexList[endIndex];
    std::vector<polymesh::Edge*>::const_iterator vei = pVertex->_edgeList.begin();
    std::vector<polymesh::Edge*>::const_iterator vei_end = pVertex->_edgeList.end();
    for(;vei!=vei_end;vei++)  // iterate over edges associated with this vertex
    {
      polymesh::Edge* pAssocEdge = *vei;
      if(pAssocEdge == pEdge) continue; // Ignore this edge
      if(pAssocEdge->_type == polymesh::Edge::BOUNDARY)
      { // we only want the next BOUNDARY Edge

        pEdge->_next = pAssocEdge;

      }

    }
    long beginIndex = pEdge->_vertexIndex[0]; // Start point of edge;
    pVertex = _vertexList[beginIndex];
    vei = pVertex->_edgeList.begin();
    vei_end = pVertex->_edgeList.end();
    for(;vei!=vei_end;vei++)  // iterate over edges associated with this vertex
    {
      polymesh::Edge* pAssocEdge = *vei;
      if(pAssocEdge == pEdge) continue; // Ignore this edge
      if(pAssocEdge->_type == polymesh::Edge::BOUNDARY)
      { // we only want the next BOUNDARY Edge

        pEdge->_previous = pAssocEdge;

      }
    }
  } // End loop over unsorted edges
  
}
  
void GarmentAtlas::fitParabola()
{

  assert(_edgeSectionsCategorised == true);
  assert( _topEdgeSections.size() > 0);
  assert( _bottomEdgeSections.size() > 0);
  
  // Run through all candidate Edges and select 3 positions leftmost, centre and rightmost
  // These points used to fit a parabola to the top edge
  
  // Vec3d compares X co-ord first 
  // can't use a set here as we want to access the middle of the range
  // so sort a vector after populating it.
  
  std::vector<EdgeSection*>::const_iterator esi, esi_end;
  std::vector<polymesh::Edge*>::const_iterator ei, ei_end;
  std::vector<Vec3d> sortedPoints;
  
  
  esi = _bottomEdgeSections.begin();
  esi_end = _bottomEdgeSections.end();
  
  for(;esi!=esi_end;esi++)
  {
    EdgeSection* pEdgeSection = *esi;
    ei = pEdgeSection->_edges.begin();
    ei_end = pEdgeSection->_edges.end();
    for(;ei!=ei_end;ei++)
    {
      polymesh::Edge* pEdge = *ei;
      sortedPoints.push_back( *(_lcmodel->tex_coords()[ pEdge->_vertexIndex[0] ]) );
      sortedPoints.push_back( *(_lcmodel->tex_coords()[ pEdge->_vertexIndex[1] ]) );
    }
  }
  
  // sort on X
  sort(sortedPoints.begin(),sortedPoints.end());
  Parabola* pBottomParabola = new Parabola;
  parabolas.push_back(pBottomParabola);
  
  pBottomParabola->fit_points[0] = sortedPoints.front();
  pBottomParabola->fit_points[1] = sortedPoints.at( (sortedPoints.size() - 1) / 2 );
  pBottomParabola->fit_points[2] = sortedPoints.back();
  
  pBottomParabola->makeFit();
  
  
  esi = _topEdgeSections.begin();
  esi_end = _topEdgeSections.end();
  sortedPoints.clear();
  
  for(;esi!=esi_end;esi++)
  {
    EdgeSection* pEdgeSection = *esi;
    ei = pEdgeSection->_edges.begin();
    ei_end = pEdgeSection->_edges.end();
    for(;ei!=ei_end;ei++)
    {
      polymesh::Edge* pEdge = *ei;
      sortedPoints.push_back( *(_lcmodel->tex_coords()[ pEdge->_vertexIndex[0] ]) );
      sortedPoints.push_back( *(_lcmodel->tex_coords()[ pEdge->_vertexIndex[1] ]) );
    }
  }
  
  // sort on X
  sort(sortedPoints.begin(),sortedPoints.end());

  Parabola* pTopParabola = new Parabola;
  parabolas.push_back(pTopParabola);
  
  pTopParabola->fit_points[0] = sortedPoints.front();
  pTopParabola->fit_points[1] = sortedPoints.at( (sortedPoints.size() - 1) / 2 );
  pTopParabola->fit_points[2] = sortedPoints.back();
  
  pTopParabola->makeFit();
  
    // If we have more than one edge section, assume sleeves and calculate an extra parabola below the topmost
  // edge on either side.
  // FIXME too many hard coded assumptions about parabola numbers, make more general
  if ( _leftEdgeSections.size() >1 && _rightEdgeSections.size() > 1) // FIXME enable once we have an example
  {
    std::cout << "We have a TShirt" << std::endl; // FIXME remove debug printing
    // we need to calculate an extra parabola.
    sort( _leftEdgeSections.begin(), _leftEdgeSections.end(), esPtrYLess() );
    // now sort the edges in top EdgeSection, from least to greatest Y
    EdgeSection *pEdgeSection = _leftEdgeSections[ _leftEdgeSections.size() -1 ];
    sort( pEdgeSection->_edges.begin(), pEdgeSection->_edges.end(), edgePtrYLess() );
    // So we should have the bottom edge in the top EdgeSection
    

    // So take point with least Y from leastY edge
    polymesh::Edge* pEdge = pEdgeSection->_edges[0];
    Vec3d* pLeftIntercept = pEdge->_vertexPtr[0];
    if(pEdge->_vertexPtr[1]->y() < pLeftIntercept->y())
      pLeftIntercept = pEdge->_vertexPtr[1];
    
    // Determine the right hand intercept
    
    sort( _rightEdgeSections.begin(), _rightEdgeSections.end(), esPtrYLess() );
    // now sort the edges in top EdgeSection, from least to greatest Y
    pEdgeSection = _rightEdgeSections[ _rightEdgeSections.size() -1 ];
    sort( pEdgeSection->_edges.begin(), pEdgeSection->_edges.end(), edgePtrYLess() );
    // So we should have the bottom edge in the top EdgeSection
    

    
    // So take point with least Y from leastY edge
    pEdge = pEdgeSection->_edges[0];
    Vec3d* pRightIntercept = pEdge->_vertexPtr[0];
    if(pEdge->_vertexPtr[1]->y() < pRightIntercept->y())
      pRightIntercept = pEdge->_vertexPtr[1];
    
    // now fit a parabola to these two points plus a centre point formed from interpolating the existing 
    // two parabola
    
    // Find the percentage height of the average Y of these intercept points. 
    double aveY = (pLeftIntercept->y() + pRightIntercept->y())/2.0;
    double aveX = (pLeftIntercept->x() + pRightIntercept->x())/2.0;
    // Convert to alpha in the range between top and bottom parabola
    double yRange = pTopParabola->fit_points[0].y() - pBottomParabola->fit_points[0].y();
    double alpha = (aveY - pBottomParabola->fit_points[0].y()) / yRange;
    std::cout << " alpha: " << alpha << std::endl;
    
    double centreY = Parabola::yInterpolate( *pBottomParabola, *pTopParabola, alpha, 0.5 );
    long index = findClosestTCTo( aveX, centreY );
    
    Parabola* pMiddleParabola = new Parabola;
    pMiddleParabola->fit_points[0] = *pLeftIntercept;
    pMiddleParabola->fit_points[1] = *(_lcmodel->tex_coords()[ index ]);
    pMiddleParabola->fit_points[2] = *pRightIntercept;
    pMiddleParabola->makeFit();
    std::vector<Parabola*>::iterator p = parabolas.begin();
    // insert after bottom parabola (and before top).
    parabolas.insert( ++p, pMiddleParabola );

    
  }
  
  // END PASTE
}

  
  
void GarmentAtlas::determineEdgesUsingSeams() 
{

  if(!_boundaryEdgesLinked) linkBoundaryEdges();
  _listOfSeamEdgeSections.clear(); // FIXME memory leak
  assert(_unsortedBoundaryEdges.size() >0);
  
      
    // Darts also need splitting into their own edge sections.
    // mark edges which are darts
  determineDartPoints();
  
  std::cout<< "Determine Edges for: " << name << " using co-incident point information" << std::endl;
  
  std::vector<polymesh::Edge*>::const_iterator ei = _unsortedBoundaryEdges.begin();
  std::vector<polymesh::Edge*>::const_iterator ei_end = _unsortedBoundaryEdges.end();
  for(;ei!=ei_end;ei++)
  {
    (*ei)->_borderEdgeType = polymesh::Edge::MIDDLE; // clear any previous attempts
  }
  
  // use the boundary edge connectivity information to determine edge sections.
  // Consider the coincident vertex count as the boundary is followed. The gradient
  // of this determines where to split edges into separate sections. The max possible
  // given no more than three atlases can join at one point is 2 other coincident points.
  // begin an EdgeSection where gradient is negative, unless first vertex co-incidence count is 2. (already accounted for by next case)
  // begin an EdgeSection at the next edge where gradient is positive
  
  int countOfEdgesLeftToProcess = _unsortedBoundaryEdges.size();
  polymesh::Edge* pStartEdge = _unsortedBoundaryEdges.at(0);
  polymesh::Edge* pCurrentEdge = pStartEdge;
  polymesh::Edge* pNextEdge =0;
  std::vector<polymesh::Edge*> edgeSectionBeginnings;
  int edgeCount=0;
    do {
      pNextEdge = pCurrentEdge->_next;
      
      countOfEdgesLeftToProcess--;
      //std::cout << "[" << countOfEdgesLeftToProcess << "]CurrentEdge: " << pCurrentEdge << std::endl;
      
      polymesh::Vertex* pV0 = _vertexList.at( pCurrentEdge->_vertexIndex[0] );
      polymesh::Vertex* pV1 = _vertexList.at( pCurrentEdge->_vertexIndex[1] );
      
      //std::cout << "v0 count: " << pV0->_coincidentPointList.size() << std::endl;
      
      if(pV0->_coincidentPointList.size() != pV1->_coincidentPointList.size())
      {
        std::cout << "[" << edgeCount << "] v0 coincd point count: " << pV0->_coincidentPointList.size()
            << " v1 coincd point count: " << pV1->_coincidentPointList.size() << std::endl;
        if(pV0->_coincidentPointList.size() < 2) // 2 shouldn't occur in ideal case of one point for 4 atlas join (should be 3 other co-incident maximum), but our data is not ideal
        { // 'falling' edges already accounted for
          if(pV0->_coincidentPointList.size() < pV1->_coincidentPointList.size()) // +ve gradient
          {
            pNextEdge->_borderEdgeType = polymesh::Edge::BEGIN;
            //edgeSectionBeginnings.push_back(pNextEdge);
            std::cout << "Begin (+ve grad) index is: " << pNextEdge->_vertexIndex[0] << std::endl;
          }
          else // -ve gradient
          {
            pCurrentEdge->_borderEdgeType = polymesh::Edge::BEGIN;
            //edgeSectionBeginnings.push_back(pCurrentEdge);
            std::cout << "Begin (-ve grad) index is: " << pCurrentEdge->_vertexIndex[0] << std::endl;
          }
        }
      }
      
      // Darts have separate begins placed inline before making edge sections
      // This is to avoid normal edge sections being lumped with dart edge sections
      if(pNextEdge->isDartEdge && pCurrentEdge->isDartBoundary)
      {
        pNextEdge->_borderEdgeType = polymesh::Edge::BEGIN;
        //edgeSectionBeginnings.push_back(pNextEdge);
        std::cout << "Begin dart (+ve grad) index is: " << pNextEdge->_vertexIndex[0] << std::endl;
      }
      else if(pNextEdge->isDartBoundary && pCurrentEdge->isDartEdge)
      {
        pNextEdge->_borderEdgeType = polymesh::Edge::BEGIN;
        //edgeSectionBeginnings.push_back(pNextEdge);
        std::cout << "Begin dart (+ve grad) index is: " << pNextEdge->_vertexIndex[0] << std::endl;
      }
      
      pCurrentEdge = pNextEdge;
      edgeCount++;
      
    } while ((pCurrentEdge != pStartEdge) && (countOfEdgesLeftToProcess>-5));
    if(countOfEdgesLeftToProcess<0) 
    {
      std::cout << "Error with edge linking, didn't arrive back at start!" << std::endl;
      std::cout << "Count left was: " << countOfEdgesLeftToProcess << std::endl;
    }

    // find the (unique) beginning sections
    ei = _unsortedBoundaryEdges.begin();
    ei_end = _unsortedBoundaryEdges.end();
    
    pStartEdge = _unsortedBoundaryEdges.at(0);
    pCurrentEdge = pStartEdge;
    do
    {
      if(pCurrentEdge->_borderEdgeType == polymesh::Edge::BEGIN)
      {
        edgeSectionBeginnings.push_back(pCurrentEdge);
      }
      pCurrentEdge = pCurrentEdge->_next;
    }
    while(pCurrentEdge != pStartEdge);
    
    
  
  // For each edgeSection BEGIN take the next edge and add follow chain adding edges to an edgeSection list
  // until the next end edge is encountered.
    ei = edgeSectionBeginnings.begin();
    ei_end = edgeSectionBeginnings.end();
    std::cout << "Begin seam edge sections count for:" << name << "is: " << edgeSectionBeginnings.size() << std::endl;
  
  for(;ei!=ei_end;ei++)
  {
    pCurrentEdge = *ei;
    pNextEdge = pCurrentEdge->_next;
    EdgeSection* edgeSection = new EdgeSection();
    _listOfSeamEdgeSections.push_back(edgeSection);
    edgeSection->_edges.push_back( pCurrentEdge );
    while(pNextEdge->_borderEdgeType != polymesh::Edge::BEGIN)
    {
      pCurrentEdge = pNextEdge;
      pNextEdge = pCurrentEdge->_next;
      edgeSection->_edges.push_back(pCurrentEdge);
    }
  }
  
  // Add connectivity information to the EdgeSections. The list is in order so trivial.
  edgeCount=0;
  for(int i = 0;i<_listOfSeamEdgeSections.size();i++) 
  {
    int indexPrevious = (i-1);
    if(indexPrevious<0) indexPrevious = _listOfSeamEdgeSections.size() - 1;
    int indexNext = (i+1)%_listOfSeamEdgeSections.size();
    EdgeSection* pES = _listOfSeamEdgeSections.at(i);
    pES->next = _listOfSeamEdgeSections.at(indexNext);
    pES->previous = _listOfSeamEdgeSections.at(indexPrevious);
    
    // also calculate the average attributes for each edge section
    pES->recalculateAttributes(_lcmodel);
    std::cout << "Edge[" << edgeCount++ << "] length: " << pES->length <<std::endl;
  }
  
  _seamEdgesDetermined = true;
  

  
  std::vector<EdgeSection*> removeSections;
  edgeCount=0;
  for(int i = 0;i<_listOfSeamEdgeSections.size();i++) 
  {
    EdgeSection* pES = _listOfSeamEdgeSections.at(i);
    std::vector<polymesh::Edge*>::const_iterator ei = pES->_edges.begin();
    int dartEdgeCount=0;
    for(;ei!=pES->_edges.end();ei++)
    {
      polymesh::Edge* pEdge = *ei;
      
      if(pEdge->isDartEdge)
      {
        dartEdgeCount++;
      }
    }
    if(dartEdgeCount>0) 
    {
      pES->isDartSeam = true;
      std::cout << name << ". Edge[" << i << "] length: " << pES->length << " IS A SEAM" << std::endl;
      removeSections.push_back(pES);
    }
  }
  
  // remove dart sections from consideration
  
  // remove the apexes of dart sections from consideration
  
  std::vector<EdgeSection*>::iterator esi = _listOfSeamEdgeSections.begin();
  for(;esi != _listOfSeamEdgeSections.end();
       esi++)
  {
    EdgeSection* pES = *esi;
    if(pES->_edges.size() == 2)
      // Check both ends are dart boundaries
      if((*(pES->_edges.begin()))->isDartBoundary)
        if((*(pES->_edges.end()-1))->isDartBoundary)
          removeSections.push_back(pES);
  }
  
  esi = removeSections.begin();
  for(;esi != removeSections.end();
       esi++)
  {
    EdgeSection* pES = *esi;
    pES->previous->next = pES->next;
    pES->next->previous = pES->previous;

    
    std::vector<EdgeSection*>::iterator es_it = std::find(_listOfSeamEdgeSections.begin(),_listOfSeamEdgeSections.end(),*esi);
    _listOfSeamEdgeSections.erase( es_it );


  }
  

  

}


bool GarmentAtlas::addVertex(polymesh::Vertex* pVertex) {
  if(!pVertex) {
    std::cerr << "addVertex: passed null vertex pointer" << std::endl;
    return false;
  }
  _vertexList.push_back(pVertex);
  
  return true;
}

bool GarmentAtlas::addEdge(polymesh::Edge* pEdge) {
  if(!pEdge) {
    std::cerr << "addEdge: passed null Edge pointer" << std::endl;
    return false;
  }
  _edgeList.push_back(pEdge);
  
  long vi1 = pEdge->_vertexIndex[0];
  long vi2 = pEdge->_vertexIndex[1];
  // Order indexes
  if(vi1 > vi2) { long _tmp = vi1; vi1 = vi2; vi2 = _tmp; }
  
  
  std::pair<long, long> key(vi1, vi2);
  //std::cout << "Adding edge with key: " << key << std::endl;
  
  _edgeLookupMap[ key ] = pEdge; // Creates a simple lookup for existing edges
  return true;
}

void GarmentAtlas::determineCoincidentPoints(Garment* pGarment)
{
  // For every 3d point, in every Boundary Edge, in every GarmentAtlas, in every GarmentSection in the garment
  //  compare to every other point in every other boundary edge in every other GarmentAtlas and Section
  //  if distance between them is below small threshold. Then they are co-incident. Record the number 
  // of other points each point is co-incident with. This information is useful in determining edges and seams.
  
  assert(_unsortedBoundaryEdges.size() > 0);
  
  double distThresholdSq = 0.10 * (_averageMeshVectorLength * _averageMeshVectorLength); // 10 percent
  
  // loop over all boundary edges in this atlas
  std::vector<polymesh::Edge*>::const_iterator ei;
  for(ei = _unsortedBoundaryEdges.begin();
      ei != _unsortedBoundaryEdges.end();
      ei++)
  {
    polymesh::Edge* pCurrentEdge = *ei;
    // loop over all other garment sections
    std::vector<GarmentSection*>::const_iterator gsi = pGarment->sectionList.begin();
    for(;gsi!=pGarment->sectionList.end();gsi++)
    {
      
      // loop over all other garment atlases
      std::vector<GarmentAtlas*>::const_iterator gai = (*gsi)->atlasList.begin();
      for(;gai!=(*gsi)->atlasList.end();gai++)
      {
        if(*gai == this) continue; // these would be dart points, determined separately
        if((*gai)->_coincidentPointSearchVisited == true) continue;
      
        int compCount = 0;
        Vec3d* pCurrentVec3d = _lcmodel->points()[ pCurrentEdge->_vertexIndex[0] ];
    
        std::vector<polymesh::Edge*>::const_iterator oei;
        for(oei = (*gai)->_unsortedBoundaryEdges.begin();
            oei != (*gai)->_unsortedBoundaryEdges.end();
            oei++)
        {
          polymesh::Edge* pOtherEdge = *oei; 
          
          if(pOtherEdge == pCurrentEdge) continue; // don't compare the same point
          
          compCount++;
          
          Vec3d* pOtherVec3d = (*gai)->_lcmodel->points()[ pOtherEdge->_vertexIndex[0] ];
          double xDistSq = (pCurrentVec3d->x() - pOtherVec3d->x())*(pCurrentVec3d->x() - pOtherVec3d->x());
          double yDistSq = (pCurrentVec3d->y() - pOtherVec3d->y())*(pCurrentVec3d->y() - pOtherVec3d->y());
          double zDistSq = (pCurrentVec3d->z() - pOtherVec3d->z())*(pCurrentVec3d->z() - pOtherVec3d->z());
          double distSq = xDistSq + yDistSq + zDistSq;
          if(distSq < distThresholdSq)
          {
            polymesh::Vertex* pThisPoint = this->_vertexList[ pCurrentEdge->_vertexIndex[0] ];
            polymesh::Vertex* pOtherPoint = (*gai)->_vertexList[ pOtherEdge->_vertexIndex[0] ];
            pThisPoint->addCoincidentPoint( pOtherPoint );
            pOtherPoint->addCoincidentPoint( pThisPoint );
            
          }
          
        } // end loop over other edges

        
      } // end loop over atlases
    } // end of loop over garment sections
  }
  this->_coincidentPointSearchVisited = true; // no need to check this atlas again.
  std::cout << "Coincident points determined for atlas: " << name << std::endl;
}

void GarmentAtlas::determineDartPoints()
{
  // For every 3d point, in every Edge, in every EdgeSection
  //  compare to every other point in everyother edge in every other EdgeSection
  //  if distance between them is below small threshold. Then they are co-incident and therefore darts.
  double distThresholdSq = 0.05 * (_averageMeshVectorLength * _averageMeshVectorLength); // 5 percent
  
    std::vector<polymesh::Edge*>::const_iterator ei;
    for(ei = _unsortedBoundaryEdges.begin();
        ei != _unsortedBoundaryEdges.end();
        ei++)
    {
      polymesh::Edge* pCurrentEdge = *ei;
      
      Vec3d* pCurrentVec3d = _lcmodel->points()[ pCurrentEdge->_vertexIndex[0] ];

    
      
        std::vector<polymesh::Edge*>::const_iterator oei;
        for(oei = _unsortedBoundaryEdges.begin();
            oei != _unsortedBoundaryEdges.end();
            oei++)
        {
          polymesh::Edge* pOtherEdge = *oei; 
          
          if(pOtherEdge == pCurrentEdge) continue; // don't compare the same point
          
          Vec3d* pOtherVec3d = _lcmodel->points()[ pOtherEdge->_vertexIndex[0] ];

          
          double xDistSq = (pCurrentVec3d->x() - pOtherVec3d->x())*(pCurrentVec3d->x() - pOtherVec3d->x());
          double yDistSq = (pCurrentVec3d->y() - pOtherVec3d->y())*(pCurrentVec3d->y() - pOtherVec3d->y());
          double zDistSq = (pCurrentVec3d->z() - pOtherVec3d->z())*(pCurrentVec3d->z() - pOtherVec3d->z());
          double distSq = xDistSq + yDistSq + zDistSq;
          if(distSq < distThresholdSq)
          {
            _dartPoints.push_back(pCurrentEdge->_vertexIndex[0]); // FIXME Optimize: We have actually identified two candidates here so can eliminate half the work with a more elegant method
            _dartPointsSet.insert( pCurrentEdge->_vertexIndex[0] );
            _vertexList.at( pCurrentEdge->_vertexIndex[0] )->isDartPoint = true;
          }
          
        } // end 2nd loop over boundary edges
      
      
    } // end 1st loop over boundary edges

    // now classify edges as dart or dart boundary for later segmentation
    for(ei = _unsortedBoundaryEdges.begin();
        ei != _unsortedBoundaryEdges.end();
        ei++)
    {
      polymesh::Edge* pEdge = *ei;
      polymesh::Vertex* vertA = _vertexList.at( pEdge->_vertexIndex[0] );
      polymesh::Vertex* vertB = _vertexList.at( pEdge->_vertexIndex[1] );
      if(vertA->isDartPoint && vertB->isDartPoint)
        pEdge->isDartEdge = true;
      else if(vertA->isDartPoint || vertB->isDartPoint)
        pEdge->isDartBoundary = true;
    }
  
  std::cout << "Dart Point Count: " << _dartPoints.size() << std::endl;
  _dartPointsDetermined = true;
}


double GarmentAtlas::lengthOfInterpolatedParabola(Parabola *ppa, Parabola *ppb, double alpha, int sections) const
{
  // take sections number of lengths of a parabola interpolated at alpha between pa and pb. Where pa is below (on Y axis) pb. Returns total length
  
  double yLeft = Parabola::yInterpolate( *ppa, *ppb, alpha, 0.0 );
  double yRight = Parabola::yInterpolate( *ppa, *ppb, alpha, 1.0 );
      
  Vec3d leftInterceptPoint, rightInterceptPoint;
      
  double minDiff = 999999.9;
      
      // Find left intercept
  std::vector<EdgeSection*>::const_iterator esi, esi_end;
  esi = _leftEdgeSections.begin();
  esi_end = _leftEdgeSections.end();
  for(;esi!=esi_end;esi++)
  {
    std::vector<polymesh::Edge*>::const_iterator ei, ei_end;
        // iterate through each Edge in the left EdgeSection
    ei = (*esi)->_edges.begin();
    ei_end = (*esi)->_edges.end();
    for(;ei!=ei_end;ei++)
    {
          // Only examining one of the points in each edge to avoid duplicating effort
          // Will miss the termination point at one or other end of edgeSection but this
          // is the same as the left or right most point of either top or bottom curve
          // so shouldn't matter
      polymesh::Edge* pEdge = *ei;
      double edgeY = pEdge->_vertexPtr[0]->y();
      double yDiff = fabs(edgeY - yLeft);
      if(yDiff < minDiff) {
        leftInterceptPoint = *(pEdge->_vertexPtr[0]);
        minDiff = yDiff;
      }
          
    }
  }
      
      // find right intercept
  minDiff = 999999.9;
  esi = _rightEdgeSections.begin();
  esi_end = _rightEdgeSections.end();
  for(;esi!=esi_end;esi++)
  {
    std::vector<polymesh::Edge*>::const_iterator ei, ei_end;
    ei = (*esi)->_edges.begin();
    ei_end = (*esi)->_edges.end();
        // iterate through each Edge in the right EdgeSection
    for(;ei!=ei_end;ei++)
    {
          // Only examining one of the points in each edge to avoid duplicating effort
          // Will miss the termination point at one or other end of edgeSection but this
          // is the same as the left or right most point of either top or bottom curve
          // so shouldn't matter
      polymesh::Edge* pEdge = *ei;
      double edgeY = pEdge->_vertexPtr[0]->y();
      double yDiff = fabs(edgeY - yRight);
      if(yDiff < minDiff) {
        rightInterceptPoint = *(pEdge->_vertexPtr[0]);
        minDiff = yDiff;
      }
    }
  }
      
//      std::cout << "LeftIPoint: " << leftInterceptPoint <<std::endl;
  //    std::cout << "RightIPoint: " << rightInterceptPoint << std::endl <<std::endl;
      
      
  double xRange = rightInterceptPoint.x() - leftInterceptPoint.x();
  
  Vec3d lastPoint = leftInterceptPoint;
  double totalLength = 0.0;
  // sample along the range
  for(int j=1; j<sections+1;j++) 
  {
    double beta = j/sections;
    double y = Parabola::yInterpolate( *ppa, *ppb, alpha, beta );
    double x = leftInterceptPoint.x() + beta*xRange;
    // FIXME cut and paste wrong
    Vec3d currentPoint = Vec3d( x, y, 0.0 );
    Vec3d vectorBetween = currentPoint - lastPoint;
    totalLength += vectorBetween.norm();
    lastPoint = currentPoint;
  }
  return totalLength;
}

void GarmentAtlas::centreModelOnOrigin() {
  // Place centre of BB on origin
  std::cout << "Old BBMin: " << _lcmodel->tcbbMin() << std::endl;
  std::cout << "Old BBMax: " << _lcmodel->tcbbMax() << std::endl;
  
  Vec3d size = _lcmodel->tcbbMax() - _lcmodel->tcbbMin();
  Vec3d centre = _lcmodel->tcbbMin() + 0.5 * size;
  
  //std::cout << "Centre: " << centre << std::endl;
  // FIXME Homogeneous transform doesn't work as expected. Do manually
  _lcmodel->translateAtlas(-1.0*centre);
    
  /*
  // mat[M][N] M - row, N = col
  double mat[4][4]; // matrix in homogeneous co-ords. - use to initialise Mat44d
  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++)
      mat[i][j] = 0.0;
    
  mat[0][3] =  -1.0*centre.x();  // x transform
  mat[1][3] =  -1.0*centre.y();  // y transform
  mat[2][3] =  -1.0*centre.z();  // z transform
  
  mat[0][0] = 1.0;
  mat[1][1] = 1.0;
  mat[2][2] = 1.0;
  mat[3][3] = 1.0;
  
    
  Mat44d translation = Mat44d(mat);
  std::cout << "Trans...";
  _lcmodel->transformAtlas(translation);
  std::cout << "form" << std::endl;
  */
  _lcmodel->recalculateTCBB();
  std::cout << "New BBMin: " << _lcmodel->tcbbMin() << std::endl;
  std::cout << "New BBMax: " << _lcmodel->tcbbMax() << std::endl;
}



void GarmentAtlas::scaleAtlasFromMesh() {
  
  // For each Edge
  double ratio = 0;
  assert((_edgeList.size() > 0) && (_edgeList.size() < 20000)); // sanity check
  
  std::vector<polymesh::Edge*>::const_iterator ei = _edgeList.begin();
  std::vector<polymesh::Edge*>::const_iterator ei_end = _edgeList.end();
  
  for(;ei != ei_end;ei++) {
    
    Vec3d edgeVector( 
        _lcmodel->tex_coords()[ (*ei)->_vertexIndex[1] ]->x() 
        - _lcmodel->tex_coords()[ (*ei)->_vertexIndex[0] ]->x(),
    _lcmodel->tex_coords()[ (*ei)->_vertexIndex[1] ]->y() 
        - _lcmodel->tex_coords()[ (*ei)->_vertexIndex[0] ]->y(),
    0 );
    
    Vec3d meshVector( 
        _lcmodel->points()[ (*ei)->_vertexIndex[1] ]->x() 
        - _lcmodel->points()[ (*ei)->_vertexIndex[0] ]->x(),
    _lcmodel->points()[ (*ei)->_vertexIndex[1] ]->y() 
        - _lcmodel->points()[ (*ei)->_vertexIndex[0] ]->y(),
    _lcmodel->points()[ (*ei)->_vertexIndex[1] ]->z() 
        - _lcmodel->points()[ (*ei)->_vertexIndex[0] ]->z() );
    
    // FIXME ignore rare tiny values for stability (indicates some upstream problem)
    bool smallEdgeVectorNorm = false;
    if(edgeVector.norm() < 0.001f)
    {
      std::cout << "small pattern norm: " << edgeVector.norm() << std::endl;
      smallEdgeVectorNorm = true;
      continue;
    }

    if(meshVector.norm() > 100.0f)
    {
      std::cout << "ignoring large mesh norm: " << meshVector.norm() << std::endl;
      continue;
    }
    
    ratio += (meshVector.norm() / edgeVector.norm());
    _averageMeshVectorLength +=meshVector.norm();
  }
  
  
  
  _averageMeshVectorLength /= _edgeList.size();
  
  assert(!(_averageMeshVectorLength > 100.0)); // sanity check
  
  double scaleFactor = ratio / _edgeList.size();
  std::vector<Vec3d*>::const_iterator vi = _lcmodel->tex_coords().begin();
  std::vector<Vec3d*>::const_iterator vi_end = _lcmodel->tex_coords().end();
  
  
  // scale up
  for(;vi != vi_end; vi++) {
    **vi *= scaleFactor;
  }
  
  
  _lcmodel->recalculateTCBB();
}


// Populate the Edge, Vertex, Face information and determine information about edges
void GarmentAtlas::populateAtlasFromModel() 
{
  // Do this first to maintain same Vertex indexing as LCModel
  // For each vertex in the LCModel
  //  add new Vertex to GarmentAtlas (same index)
  
  _vertexList.clear();
  _faceList.clear();
  _edgeList.clear();
  _edgeLookupMap.clear();
  _unsortedBoundaryEdges.clear();
  _edgeSectionsCategorised = false;
  _controlMeshAssigned = false;
  _dartPointsDetermined =false;
  _atlas_populated = false;
  _topEdgeSections.clear();
  _bottomEdgeSections.clear();
  _leftEdgeSections.clear();
  _rightEdgeSections.clear();
  parabolas.clear();
  _dartPoints.clear();
  
  std::vector<Vec3d*>::const_iterator pl = _lcmodel->points().begin();
  std::vector<Vec3d*>::const_iterator pl_end = _lcmodel->points().end();
  long vCount=0;
  for(;pl!=pl_end;pl++) {
    polymesh::Vertex* pV = new polymesh::Vertex(vCount++);
    pV->point( *pl );
    if(!addVertex( pV )) {
      std::cerr << "ERROR: Could not add vertex" << std::endl;
    }
  }
  
  // For each face in LCModel
  std::vector<IndexedTriangle*>::const_iterator itl = _lcmodel->indexedTriangles().begin();
  std::vector<IndexedTriangle*>::const_iterator itl_end = _lcmodel->indexedTriangles().end();
  
  for(;itl!=itl_end;itl++) {
//    std::cout << "FCount: " << faceCount++ << std::endl;
    polymesh::Face* pF = new polymesh::Face( (*itl)->indexA(),
                                             (*itl)->indexB(),
                                             (*itl)->indexC()
                                           );
    
    if(!addFace( pF )) {
      std::cerr << "ERROR: Could not add Face" << std::endl;
    }
    
    // assign this face to the relevant Vertices
    vertices()[ pF->_vertex[0] ]->addFace( pF );
    vertices()[ pF->_vertex[1] ]->addFace( pF );
    vertices()[ pF->_vertex[2] ]->addFace( pF );
    
    // First pass of Edge creation. Some edges will be inconsistent but we want Face connectivity
    // So build edges. But later we will recreate edges once Face ordering is consistent.
    // for each pair of vertices in order, check if edge exists, if not create one
    polymesh::Edge* pEdge = 0;
    for(int vi = 0; vi < 3; vi++) {
      int vi_next = (vi + 1) % 3; // cycle 0-1, 1-2, 2-0
      pEdge = lookup( pF->_vertex[ vi ], pF->_vertex[ vi_next ] );
      if(!pEdge) {
        pEdge = new polymesh::Edge(pF->_vertex[ vi ], pF->_vertex[ vi_next ] );
        addEdge( pEdge );
      }
      // Add reference to this face to found or created edge
      pEdge->addFace( pF );
      pF->addEdge( pEdge );
      
     
    }
  } // End of iterator over faces
  std::cout << "Initial Vertex/Edge/Face construction complete" << std::endl;
 
  
  // Make faces consistent with initial Face
  polymesh::Face::makeMeshConsistentWith( _faceList[ 0 ], this );
  std::cout << "Faces are now consistent" << std::endl;
  
  // Second pass of Edge creation now Faces are consistent
  // first delete old Edges
  _edgeLookupMap.clear();
  std::vector<polymesh::Edge*>::const_iterator ei = _edgeList.begin();
  std::vector<polymesh::Edge*>::const_iterator ei_end = _edgeList.end();
  for(;ei != ei_end;ei++) {
    if(*ei) {
      delete (*ei);
    }
  }
  _edgeList.clear();
  
  // now create new edges and assign to faces and vertices
  std::vector<polymesh::Face*>::const_iterator fi = _faceList.begin();
  std::vector<polymesh::Face*>::const_iterator fi_end = _faceList.end();
  for(;fi!=fi_end;fi++) {
    polymesh::Face* pF = *fi;
    pF->clearEdgeList();
    // Lookup or create an Edge for each vertex pair inorder
    polymesh::Edge* pEdge = 0;
    for(int vi = 0; vi < 3; vi++) {
      int vi_next = (vi + 1) % 3; // cycle 0-1, 1-2, 2-0
      pEdge = lookup( pF->_vertex[ vi ], pF->_vertex[ vi_next ] );
      if(!pEdge) {
        pEdge = new polymesh::Edge(pF->_vertex[ vi ], pF->_vertex[ vi_next ] );
        pEdge->_vertexPtr[0] = _lcmodel->tex_coords()[ pF->_vertex[ vi ] ];
        pEdge->_vertexPtr[1] = _lcmodel->tex_coords()[ pF->_vertex[ vi_next ] ];
        pEdge->_lcmodel = _lcmodel;
        addEdge( pEdge );
      }
      // Add reference to this face to found or created edge
      pEdge->addFace( pF );
      // Add edge to face
      pF->addEdge( pEdge );
      // add edge to vertices
            // Add reference to this Edge to both vertices. Vertices ignore duplicate Edges.
      vertices()[ pF->_vertex[ vi ] ]->addEdge( pEdge );
      vertices()[ pF->_vertex[ vi_next ] ]->addEdge( pEdge );
    }
    
  }
  
  scaleAtlasFromMesh();
  
  // Mark exterior edges
   ei = _edgeList.begin();
   ei_end = _edgeList.end();
   for(;ei != ei_end;ei++) { // Begin iteration over edges for purposes of marking
    if( (*ei)->faceCount() == 1) {
      
      (*ei)->_type = polymesh::Edge::BOUNDARY;
      _unsortedBoundaryEdges.push_back( *ei );

    }
    else {
      (*ei)->_type = polymesh::Edge::INTERIOR;
    }
  } // End of iteration over edges for purposes of marking
  
      
} // end populateAtlasFromModel

/* OBSOLETE - now use seam edge sections. but keep in case required again
// Populate _left, _right, _top, _bottom edge section lists for use by parabola and left/right intercept point determination
// Use constraints for required number of edges supplied by user
void GarmentAtlas::categoriseEdgeSections()
{
  std::cout << "Categorising edges for: " << name << std::endl;
    // Categorise the edges. We want to determine connectivity and then split into edge sections

  linkBoundaryEdges();
  
  // loop again, this time compare this->edgeNormal with next->edgeNormal using dotProduct
  //  if angle between is larger than threshold than mark this edge as end of edge section 
  //  and add to a list of edgeSection ends.
  
  
  polymesh::Edge* pStartEdge = _unsortedBoundaryEdges.at(0);
  polymesh::Edge* pCurrentEdge = pStartEdge;
  polymesh::Edge* pNextEdge;
  std::vector<polymesh::Edge*> edgeSectionBeginnings;
  double angle_threshold = M_PI / 3.0;
  do {
    std::cout << "Finding edges for: " << name << ", with threshold at: " << angle_threshold << std::endl;
    edgeSectionBeginnings.clear();
    do {
      pNextEdge = pCurrentEdge->_next;
      Vec3d n1=pCurrentEdge->normal();
      Vec3d n2=pNextEdge->normal();
      double dp = utils::dotProduct( n1.x(), n1.y(), n2.x(), n2.y() );
      double productOfNorms = n1.norm() * n2.norm();
      double theta = acos(dp / productOfNorms);
      //std::cout << "Theta: " << theta << std::endl;
      if(theta > angle_threshold)
      {
        //std::cout << "EDGE" <<  std::endl;
        //pCurrentEdge->_borderEdgeType = polymesh::Edge::END;
        pNextEdge->_borderEdgeType = polymesh::Edge::BEGIN;
        edgeSectionBeginnings.push_back(pNextEdge);
      }
      pCurrentEdge = pNextEdge;
    } while (pCurrentEdge != pStartEdge);
    angle_threshold /= 2.0; // reduce threshold if not enough edges detected
    pCurrentEdge = pStartEdge;
  } while (edgeSectionBeginnings.size() < 4); // We need at least 4 edges for classification
  
  
  
  
  // For each edgeSection BEGIN take the next edge and add follow chain adding edges to an edgeSection list
  // until the next end edge is encountered.
  std::vector<polymesh::Edge*>::const_iterator ei = edgeSectionBeginnings.begin();
  std::vector<polymesh::Edge*>::const_iterator ei_end = edgeSectionBeginnings.end();
  std::cout << "Begin sections count: " << edgeSectionBeginnings.size() << std::endl;
  
  for(;ei!=ei_end;ei++) 
  {
    pCurrentEdge = *ei;
    pNextEdge = pCurrentEdge->_next;
    EdgeSection* edgeSection = new EdgeSection();
    _listOfBoundaryEdgeSections.push_back(edgeSection);
    edgeSection->_edges.push_back( pCurrentEdge );
    while(pNextEdge->_borderEdgeType != polymesh::Edge::BEGIN)
    {
      pCurrentEdge = pNextEdge;
      pNextEdge = pCurrentEdge->_next;
      edgeSection->_edges.push_back(pCurrentEdge);
    }
  }
  
  // Add connectivity information to the EdgeSections. The list is in order so trivial.
  for(int i = 0;i<_listOfBoundaryEdgeSections.size();i++) 
  {
    int indexPrevious = (i-1);
    if(indexPrevious<0) indexPrevious = _listOfBoundaryEdgeSections.size() - 1;
    int indexNext = (i+1)%_listOfBoundaryEdgeSections.size();
    EdgeSection* pES = _listOfBoundaryEdgeSections.at(i);
    pES->next = _listOfBoundaryEdgeSections.at(indexNext);
    pES->previous = _listOfBoundaryEdgeSections.at(indexPrevious);
  }
  
  
  
  std::cout << "Separate edge count: " << _listOfBoundaryEdgeSections.size() << std::endl;
  
  // Now classify EdgeSections
  std::vector<EdgeSection*>::const_iterator esi = _listOfBoundaryEdgeSections.begin();
  std::vector<EdgeSection*>::const_iterator esi_end = _listOfBoundaryEdgeSections.end();
  // Look at edges in each section. Calculate average position and normal

  for(;esi!=esi_end;esi++)  
  {
    EdgeSection* pEdgeSection = *esi;
    
     
    Vec3d accumulatePosition;
    Vec3d accumulateNormal;
    double accumulateLength;
    ei = pEdgeSection->_edges.begin();
    ei_end = pEdgeSection->_edges.end();
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
    pEdgeSection->averagePos = accumulatePosition / (double)pEdgeSection->_edges.size();
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
    
    
  } // end iterator over edge sections
  std::cout << "Edges assigned orientations" << std::endl;
  
  
  
  // Find top and bottom edges
  //FIXME Use constraints set by user.

  
  std::set<EdgeSection*,esPtrYLess> setESTBptr;
  esi = _listOfBoundaryEdgeSections.begin();
  esi_end = _listOfBoundaryEdgeSections.end();
  for(;esi!=esi_end;esi++)
  {
    setESTBptr.insert( *esi );
  }
  
  std::cout << "Inserted EdgeSections" << std::endl;
  
  // FIXME use a vector for post-sorting instead of a set
  std::vector<EdgeSection*> sortedESTB;
  for(std::set<EdgeSection*,esPtrYLess>::const_iterator p = setESTBptr.begin();
      p != setESTBptr.end();
      p++
     )
  {
    //std::cout << "Ordered: " << (*p)->averagePos << std::endl;
    sortedESTB.push_back( *p );
  }
  
  std::cout << "Sorted EdgeSections" << std::endl;
  
  // Setup range
  EdgeSection* pLeastY = sortedESTB.at(0);
  
  EdgeSection* pGreatestY = sortedESTB.at( sortedESTB.size()-1 );
  double yRange = pGreatestY->averagePos.y() - pLeastY->averagePos.y();
  double yWindow = yRange / 10.0; // Find edge sections with nearest 10% of y().

  
  esi = _listOfBoundaryEdgeSections.begin();
  esi_end = _listOfBoundaryEdgeSections.end();
  for(;esi!=esi_end;esi++)
  {
    EdgeSection* pEdgeSection=*esi;
    if(pEdgeSection->_orientation == EdgeSection::UP) {
      if(pEdgeSection->averagePos.y() <= (pGreatestY->averagePos.y() + yWindow/2.0)
         &&
         pEdgeSection->averagePos.y() >= (pGreatestY->averagePos.y() - yWindow/2.0))
      { // close to the edge at the top, so include it
        
        _topEdgeSections.push_back(pEdgeSection);
        pEdgeSection->_classification = EdgeSection::TOP;
      }
    }
    else if(pEdgeSection->_orientation == EdgeSection::DOWN) {
      if(pEdgeSection->averagePos.y() <= (pLeastY->averagePos.y() + yWindow/2.0)
         &&
         pEdgeSection->averagePos.y() >= (pLeastY->averagePos.y() - yWindow/2.0))
      { // close to the edge at the top, so include it
        
        _bottomEdgeSections.push_back(pEdgeSection);
        pEdgeSection->_classification = EdgeSection::BOTTOM;
      }
    }
  }
  
  std::cout << "Initial Top EdgeSection count: " <<_topEdgeSections.size() << std::endl;
  std::cout << "Initial Bottom EdgeSection count: " <<_bottomEdgeSections.size() << std::endl;
  
  // Now determine left and right edge sections
  
  std::set<EdgeSection*,esPtrXLess> setESLRptr;
  esi = _listOfBoundaryEdgeSections.begin();
  esi_end = _listOfBoundaryEdgeSections.end();
  for(;esi!=esi_end;esi++)
  {
    setESLRptr.insert( *esi );
  }
  
  // Sort edge sections on X
  std::vector<EdgeSection*> sortedESLR;
  for(std::set<EdgeSection*,esPtrXLess>::const_iterator p = setESLRptr.begin();
      p != setESLRptr.end();
      p++
     )
  {
    //std::cout << "Ordered: " << (*p)->averagePos << std::endl;
    sortedESLR.push_back( *p );
  }
  
  //std::cout << "Sorted EdgeSections" << std::endl;
  
  // Setup range
  EdgeSection* pLeastX = sortedESLR.at(0);
  
  EdgeSection* pGreatestX = sortedESLR.at( sortedESLR.size()-1 );
  double xRange = pGreatestX->averagePos.x() - pLeastX->averagePos.x();
  double xWindow = xRange / 10.0; // Find edge sections with nearest 10% of x().
  
  
  
  esi = _listOfBoundaryEdgeSections.begin();
  esi_end = _listOfBoundaryEdgeSections.end();
  for(;esi!=esi_end;esi++)
  {
    EdgeSection* pEdgeSection=*esi;
    if(pEdgeSection->_orientation == EdgeSection::RIGHT) {
      if(pEdgeSection->averagePos.x() <= (pGreatestX->averagePos.x() + xWindow/2.0)
         &&
         pEdgeSection->averagePos.x() >= (pGreatestX->averagePos.x() - xWindow/2.0))
      { // close to the edge at the top, so include it
        
        _rightEdgeSections.push_back(pEdgeSection);
        pEdgeSection->_classification = EdgeSection::RIGHT_SIDE;
      }
    }
    else if(pEdgeSection->_orientation == EdgeSection::LEFT) {
      if(pEdgeSection->averagePos.x() <= (pLeastX->averagePos.x() + xWindow/2.0)
         &&
         pEdgeSection->averagePos.x() >= (pLeastX->averagePos.x() - xWindow/2.0))
      { // close to the edge at the top, so include it
        
        _leftEdgeSections.push_back(pEdgeSection);
        pEdgeSection->_classification = EdgeSection::LEFT_SIDE;
      }
    }
  }
  
  std::cout << "Number of left edge sections for finding intercepts: " << _leftEdgeSections.size() << std::endl;
  std::cout << "Number of right edge sections for finding intercepts: " << _rightEdgeSections.size() << std::endl;
  
  // Constrain
  if(_leftEdgeSections.size() > leftEdgesReq)
  {
    std::cout << "Merging left edgesections..." << std::endl;    
    do {
      std::cout << "Iteration begins." << std::endl;
      // find neighbours to the edge sections.
      for(int i=0;i<_leftEdgeSections.size();i++)
      {
        EdgeSection* pES = _leftEdgeSections.at(i);
        EdgeSection* pNeighbouringES = 0; // will point to ES we want to merge
        
        // Look for shortest neighbouring unclassified edge section
        
        if(pES->next->_classification == EdgeSection::UNCLASSIFIED
          && pES->previous->_classification == EdgeSection::UNCLASSIFIED)
        {
          pNeighbouringES = pES->next->length < pES->previous->length ? pES->next : pES->previous;
        }
        else if(pES->next->_classification == EdgeSection::UNCLASSIFIED)
        {
          pNeighbouringES = pES->next;
        }
        else if(pES->previous->_classification == EdgeSection::UNCLASSIFIED)
        {
          pNeighbouringES = pES->previous;
        }
        else
        {
          std::cout << "No unclassified ES neighbours, continuing" << std::endl;
          continue; // No unclassified edges, move on
        }
        
        // Examine the normals of the neighbours of this neighbour, and then merge this neighbour with
        // the one with closest normal.
        double anglePrevious;
        double angleNext;
        EdgeSection *pESMergeCandidate = 0;
        {
          double dp = utils::dotProduct(pNeighbouringES->previous->averageNormal.x(),pNeighbouringES->previous->averageNormal.y()
                ,pNeighbouringES->averageNormal.x(), pNeighbouringES->averageNormal.y());
          double productOfNorms = pNeighbouringES->previous->averageNormal.norm() *  
                    pNeighbouringES->averageNormal.norm();
          anglePrevious = acos(dp / productOfNorms); // FIXME check for division by zero - could happen
        }
        {
          double dp = utils::dotProduct(
                pNeighbouringES->averageNormal.x(), pNeighbouringES->averageNormal.y()
                ,pNeighbouringES->next->averageNormal.x(),pNeighbouringES->next->averageNormal.y());
          double productOfNorms = pNeighbouringES->next->averageNormal.norm() *  
                pNeighbouringES->averageNormal.norm();
          angleNext = acos(dp / productOfNorms); // FIXME check for division by zero - could happen
        }
        pESMergeCandidate = anglePrevious < angleNext ? pNeighbouringES->previous : pNeighbouringES->next;
        
        pESMergeCandidate->merge( pNeighbouringES );
        pESMergeCandidate->recalculateAttributes( _lcmodel );
        std::vector<EdgeSection*>::iterator esi = _listOfBoundaryEdgeSections.begin();
        
        // Delete merged edge section from vector
        for(;esi!= _listOfBoundaryEdgeSections.end();
            esi++)
        {
          if(*esi == pNeighbouringES)
          {
            _listOfBoundaryEdgeSections.erase ( esi);
            delete pNeighbouringES;
            break;
          }
          
        } // end of search to delete merged edgesection
      } // end of loop over left edgesections
    } while (_leftEdgeSections.size() > leftEdgesReq);
    
  } // end of correcting left edge sections
  
  // repeat for RHS
  
  if(_rightEdgeSections.size() > rightEdgesReq)
  {
    std::cout << "Merging right edgesections..." << std::endl;    
    do {
      std::cout << "Iteration begins." << std::endl;
      // find neighbours to the edge sections.
      for(int i=0;i<_rightEdgeSections.size();i++)
      {
        EdgeSection* pES = _rightEdgeSections.at(i);
        EdgeSection* pNeighbouringES = 0; // will point to ES we want to merge
        
        // Look for shortest neighbouring unclassified edge section
        
        if(pES->next->_classification == EdgeSection::UNCLASSIFIED
           && pES->previous->_classification == EdgeSection::UNCLASSIFIED)
        {
          pNeighbouringES = pES->next->length < pES->previous->length ? pES->next : pES->previous;
        }
        else if(pES->next->_classification == EdgeSection::UNCLASSIFIED)
        {
          pNeighbouringES = pES->next;
        }
        else if(pES->previous->_classification == EdgeSection::UNCLASSIFIED)
        {
          pNeighbouringES = pES->previous;
        }
        else
        {
          std::cout << "No unclassified ES neighbours, continuing" << std::endl;
          continue; // No unclassified edges, move on
        }
        
        // Examine the normals of the neighbours of this neighbour, and then merge this neighbour with
        // the one with closest normal.
        double anglePrevious;
        double angleNext;
        EdgeSection *pESMergeCandidate = 0;
        {
          double dp = utils::dotProduct(pNeighbouringES->previous->averageNormal.x(),pNeighbouringES->previous->averageNormal.y()
                ,pNeighbouringES->averageNormal.x(), pNeighbouringES->averageNormal.y());
          double productOfNorms = pNeighbouringES->previous->averageNormal.norm() *  
                pNeighbouringES->averageNormal.norm();
          anglePrevious = acos(dp / productOfNorms); // FIXME check for division by zero - could happen
        }
        {
          double dp = utils::dotProduct(
                pNeighbouringES->averageNormal.x(), pNeighbouringES->averageNormal.y()
                ,pNeighbouringES->next->averageNormal.x(),pNeighbouringES->next->averageNormal.y());
          double productOfNorms = pNeighbouringES->next->averageNormal.norm() *  
                pNeighbouringES->averageNormal.norm();
          angleNext = acos(dp / productOfNorms); // FIXME check for division by zero - could happen
        }
        pESMergeCandidate = anglePrevious < angleNext ? pNeighbouringES->previous : pNeighbouringES->next;
        
        pESMergeCandidate->merge( pNeighbouringES );
        pESMergeCandidate->recalculateAttributes( _lcmodel );
        std::vector<EdgeSection*>::iterator esi = _listOfBoundaryEdgeSections.begin();
        
        // Delete merged edge section from vector
        for(;esi!= _listOfBoundaryEdgeSections.end();
             esi++)
        {
          if(*esi == pNeighbouringES)
          {
            _listOfBoundaryEdgeSections.erase ( esi);
            delete pNeighbouringES;
            break;
          }
          
        } // end of search to delete merged edgesection
      } // end of loop over right edgesections
    } while (_rightEdgeSections.size() > rightEdgesReq);
    
  } // end of correcting right edge sections
  
  _edgeSectionsCategorised = true;
      
}
// END of obsolete categoriseedgesections*/

void GarmentAtlas::updateEdgeSectionInformation()
{
  assert(_listOfSeamEdgeSections.size()>0);
  std::vector<EdgeSection*>::iterator esi = _listOfSeamEdgeSections.begin();
  unsigned edgeCount=0;
  for(;esi!=_listOfSeamEdgeSections.end();esi++)
  {
    EdgeSection* pES = *esi;
    std::cout <<"Edge["<<edgeCount<<"] had length: " <<pES->length;
    pES->recalculateAttributes(_lcmodel);
    std::cout <<" now has length: "<<pES->length << std::endl;
  }
}



void GarmentAtlas::reduceEdgeSectionsToRequiredNumber(int numReq)
{
  
  std::cout << "reduceEdgeSectionsToRequiredNumber called with: " << numReq << std::endl;
  std::cout << "Size is: " << _listOfSeamEdgeSections.size() <<std::endl;
  assert(_listOfSeamEdgeSections.size() > 0);
  // iterate through the edges merging the shortest edges until the required number is acheived.
  while(_listOfSeamEdgeSections.size() > numReq)
  {
    std::vector<EdgeSection*>::iterator esi = _listOfSeamEdgeSections.begin();
    
    EdgeSection* pShortestEdgeSection=0;
    double minLength = 9999999.9;

    int count =0;
    for(;esi!= _listOfSeamEdgeSections.end();
         esi++)
    {
      EdgeSection* pES = (*esi);
      std::cout << "ES[" << count << "] has length: " << pES->length << std::endl;
      if( (pES->length < minLength) )
      {
        pShortestEdgeSection = pES;
        minLength = pES->length;
        std::cout << "ES[" << count << "] is shortest" << std::endl;
      }
      count++;
    }
    
    std::cout << "Found shortest edge section" << std::endl;
    
            // Examine the normals of the neighbours of this edgesection, and then merge this neighbour with
        // the one with closest normal.
    double anglePrevious;
    double angleNext;
    EdgeSection *pEStoMergeWith = 0;
    {
      double dp = utils::dotProduct(pShortestEdgeSection->previous->averageNormal.x(),pShortestEdgeSection->previous->averageNormal.y()
            ,pShortestEdgeSection->averageNormal.x(), pShortestEdgeSection->averageNormal.y());
      double productOfNorms = pShortestEdgeSection->previous->averageNormal.norm() *  
            pShortestEdgeSection->averageNormal.norm();
      anglePrevious = acos(dp / productOfNorms); // FIXME check for division by zero - could happen
    }
    {
      double dp = utils::dotProduct(
            pShortestEdgeSection->averageNormal.x(), pShortestEdgeSection->averageNormal.y()
            ,pShortestEdgeSection->next->averageNormal.x(),pShortestEdgeSection->next->averageNormal.y());
      double productOfNorms = pShortestEdgeSection->next->averageNormal.norm() *  
            pShortestEdgeSection->averageNormal.norm();
      angleNext = acos(dp / productOfNorms); // FIXME check for division by zero - could happen
    }
    pEStoMergeWith = anglePrevious < angleNext ? pShortestEdgeSection->previous : pShortestEdgeSection->next;
    

    
    pEStoMergeWith->merge( pShortestEdgeSection );
        
    
  
    pEStoMergeWith->recalculateAttributes( _lcmodel );
    
    
        
    std::cout << "Deleting merged edge section" << std::endl;
    

    // relink
    if(pEStoMergeWith == pShortestEdgeSection->previous) // we merged with previous edge
    {
      pEStoMergeWith->next = pShortestEdgeSection->next;
    }
    else if(pEStoMergeWith == pShortestEdgeSection->next)// we merged with next edge
    {
      pEStoMergeWith->previous = pShortestEdgeSection->previous;
    }
    else
    {
      std::cout << "Relink error" << std::endl;
      assert(1==2);
    }
    
    // delete
    std::vector<EdgeSection*>::iterator p = std::find(_listOfSeamEdgeSections.begin(),_listOfSeamEdgeSections.end(),pShortestEdgeSection);
    _listOfSeamEdgeSections.erase ( p );
    delete pShortestEdgeSection;

  } // end while too many edge sections
    
}

long GarmentAtlas::findClosestTCTo(double x, double y)
{
  double minDist = 99999999.9;
  long closestIndex=0;
  long indexCount=0;
  // FIXME Uh oh. Brute force first. Optimize later
//  std::cout << "Searching for a texture co-ordinate..." << std::endl;
  std::vector<Vec3d*>::const_iterator p;
  for(p = _lcmodel->tex_coords().begin();
      p!= _lcmodel->tex_coords().end();
      p++)
  {
    double dist = ((*p)->x() - x) * ((*p)->x() - x) + ((*p)->y() - y) * ((*p)->y() - y);
    if(dist < minDist) 
    {
      minDist = dist;
      closestIndex = indexCount;
    }
    indexCount++;
  }
  return closestIndex;
}


// If closest point is a dart, keep looking...
long GarmentAtlas::findClosestTCToIgnoringDarts(double x, double y)
{
  double minDist = 99999999.9;
  long closestIndex=0;
  long indexCount=0;
  // FIXME Uh oh. Brute force first. Optimize later
//  std::cout << "Searching for a texture co-ordinate..." << std::endl;
  std::vector<Vec3d*>::const_iterator p;
  for(p = _lcmodel->tex_coords().begin();
      p!= _lcmodel->tex_coords().end();
      p++)
  {
    double dist = ((*p)->x() - x) * ((*p)->x() - x) + ((*p)->y() - y) * ((*p)->y() - y);
    if(dist < minDist) 
    {
      if( _vertexList.at( indexCount )->isDartPoint)
      {
        std::cout << "Ignoring dart point" << std::endl;
      }
      else
      {
        minDist = dist;
        closestIndex = indexCount;
      }
    }
    indexCount++;
  }
  return closestIndex;
}

long GarmentAtlas::findClosest3TCTo(double x, double y, pointDistance* npArray)
{
  double minDist = 99999999.9;
  long indexCount=0;
  long closestIndex=0;
  // FIXME Uh oh. Brute force first. Optimize later
//  std::cout << "Searching for a texture co-ordinate..." << std::endl;
  
  assert(_lcmodel->tex_coords().size() > 0);

  std::vector<pointDistance> nearestPoints;
  
  std::vector<Vec3d*>::const_iterator p;
  for(p = _lcmodel->tex_coords().begin();
      p!= _lcmodel->tex_coords().end();
      p++)
  {
    
    double dist = ((*p)->x() - x) * ((*p)->x() - x) + ((*p)->y() - y) * ((*p)->y() - y);

    pointDistance pd;
    
    pd.index = indexCount; 
    pd.dist = dist;
    pd.pPoint = _lcmodel->tex_coords()[ pd.index ];
    
    nearestPoints.push_back(pd);
    
    indexCount++;
  }
  
  // now sort closest first

  sort(nearestPoints.begin(),nearestPoints.end(),closestPointFunctor());
  assert(nearestPoints.size() >= 3);
  
  npArray[0] = nearestPoints.at(0);
  npArray[1] = nearestPoints.at(1);
  npArray[2] = nearestPoints.at(2);
}

void GarmentAtlas::categoriseSeamEdgeSections()
{
  
  
  // Find top and bottom edges
  //FIXME Use constraints set by user.

  std::vector<EdgeSection*>::const_iterator esi, esi_end;
  
  std::set<EdgeSection*,esPtrYLess> setESTBptr;
  esi = _listOfSeamEdgeSections.begin();
  esi_end = _listOfSeamEdgeSections.end();
  for(;esi!=esi_end;esi++)
  {
    setESTBptr.insert( *esi );
  }
  
  std::cout << "Inserted EdgeSections" << std::endl;
  
  // FIXME use a vector for post-sorting instead of a set
  std::vector<EdgeSection*> sortedESTB;
  for(std::set<EdgeSection*,esPtrYLess>::const_iterator p = setESTBptr.begin();
      p != setESTBptr.end();
      p++
     )
  {
    //std::cout << "Ordered: " << (*p)->averagePos << std::endl;
    sortedESTB.push_back( *p );
  }
  
  std::cout << "Sorted EdgeSections" << std::endl;
  
  // Setup range
  EdgeSection* pLeastY = sortedESTB.at(0);
  
  EdgeSection* pGreatestY = sortedESTB.at( sortedESTB.size()-1 );
  double yRange = pGreatestY->averagePos.y() - pLeastY->averagePos.y();
  double yWindow = yRange / 10.0; // Find edge sections with nearest 10% of y().

  
  esi = _listOfSeamEdgeSections.begin();
  esi_end = _listOfSeamEdgeSections.end();
  for(;esi!=esi_end;esi++)
  {
    EdgeSection* pEdgeSection=*esi;
    if(pEdgeSection->_orientation == EdgeSection::UP) {
      if(pEdgeSection->averagePos.y() <= (pGreatestY->averagePos.y() + yWindow/2.0)
         &&
         pEdgeSection->averagePos.y() >= (pGreatestY->averagePos.y() - yWindow/2.0))
      { // close to the edge at the top, so include it
        
        _topEdgeSections.push_back(pEdgeSection);
        pEdgeSection->_classification = EdgeSection::TOP;
      }
    }
    else if(pEdgeSection->_orientation == EdgeSection::DOWN) {
      if(pEdgeSection->averagePos.y() <= (pLeastY->averagePos.y() + yWindow/2.0)
         &&
         pEdgeSection->averagePos.y() >= (pLeastY->averagePos.y() - yWindow/2.0))
      { // close to the edge at the top, so include it
        
        _bottomEdgeSections.push_back(pEdgeSection);
        pEdgeSection->_classification = EdgeSection::BOTTOM;
      }
    }
  }
  
  std::cout << "Initial Top EdgeSection count: " <<_topEdgeSections.size() << std::endl;
  std::cout << "Initial Bottom EdgeSection count: " <<_bottomEdgeSections.size() << std::endl;
  
  // Now determine left and right edge sections
  
  std::set<EdgeSection*,esPtrXLess> setESLRptr;
  esi = _listOfSeamEdgeSections.begin();
  esi_end = _listOfSeamEdgeSections.end();
  for(;esi!=esi_end;esi++)
  {
    setESLRptr.insert( *esi );
  }
  
  // Sort edge sections on X
  std::vector<EdgeSection*> sortedESLR;
  for(std::set<EdgeSection*,esPtrXLess>::const_iterator p = setESLRptr.begin();
      p != setESLRptr.end();
      p++
     )
  {
    //std::cout << "Ordered: " << (*p)->averagePos << std::endl;
    sortedESLR.push_back( *p );
  }
  
  //std::cout << "Sorted EdgeSections" << std::endl;
  
  // Setup range
  EdgeSection* pLeastX = sortedESLR.at(0);
  
  EdgeSection* pGreatestX = sortedESLR.at( sortedESLR.size()-1 );
  double xRange = pGreatestX->averagePos.x() - pLeastX->averagePos.x();
  double xWindow = xRange / 2.0; // Find edge sections with nearest 50% of x(). // FIXME made larger (obsolete? with new edge categorisation using seams
  
  
  
  esi = _listOfSeamEdgeSections.begin();
  esi_end = _listOfSeamEdgeSections.end();
  for(;esi!=esi_end;esi++)
  {
    EdgeSection* pEdgeSection=*esi;
    if(pEdgeSection->_orientation == EdgeSection::RIGHT) {
      if(pEdgeSection->averagePos.x() <= (pGreatestX->averagePos.x() + xWindow/2.0)
         &&
         pEdgeSection->averagePos.x() >= (pGreatestX->averagePos.x() - xWindow/2.0))
      { // close to the edge at the top, so include it
        
        _rightEdgeSections.push_back(pEdgeSection);
        pEdgeSection->_classification = EdgeSection::RIGHT_SIDE;
      }
    }
    else if(pEdgeSection->_orientation == EdgeSection::LEFT) {
      if(pEdgeSection->averagePos.x() <= (pLeastX->averagePos.x() + xWindow/2.0)
         &&
         pEdgeSection->averagePos.x() >= (pLeastX->averagePos.x() - xWindow/2.0))
      { // close to the edge at the top, so include it
        
        _leftEdgeSections.push_back(pEdgeSection);
        pEdgeSection->_classification = EdgeSection::LEFT_SIDE;
      }
    }
  }
  
  std::cout << "Number of left edge sections for finding intercepts: " << _leftEdgeSections.size() << std::endl;
  std::cout << "Number of right edge sections for finding intercepts: " << _rightEdgeSections.size() << std::endl;
  
  
  _edgeSectionsCategorised = true;
}


void GarmentAtlas::findBarycentricCoefficients(const Vec3d P, const Vec3d A, const Vec3d B, const Vec3d C, double *b, double *c)
{
  
    double top = (C.y() - A.y())*(P.x() - A.x()) - ((C.x()-A.x())*(P.y()-A.y()));
    double bottom = ((B.x()-A.x())*(C.y()-A.y())) - ((C.x()-A.x())*(B.y()-A.y()));
    if(fabs(bottom) < 0.00001)
    { // unstable - default to A as it is the nearest point.
      *b = 0.0; *c = 0.0; 
      return;
    }
    *b = top / bottom;
    
    top = (-1.0)*(B.y() - A.y())*(P.x() - A.x()) + ((B.x()-A.x())*(P.y()-A.y()));
    // bottom is the same
    
    *c = top / bottom;

    //std::cout << "Xp - Xa: " << P.x() - A.x() << " == " << *b*(B.x() - A.x()) + *c*(C.x() - A.x()) << std::endl;
    //std::cout << "Yp - Ya: " << P.y() - A.y() << " == " << *b*(B.y() - A.y()) + *c*(C.y() - A.y()) << std::endl;
    //std::cout << "b: " << *b << ", c: " << *c << std::endl;
    //assert((*b+*c)<=1.0);

}



void GarmentAtlas::populateControlMesh(int rows, int cols)
{
  if(!_edgeSectionsCategorised) 
  {
    std::cout << "Warning: can't populate control mesh for: " << name << " until edges categorised" << std::endl; 
    return;
  }
std::cout<<"inside POP ATLAS"<<std::endl;  
  
  int i, j;
  Vec3d cp;
  
  _controlMesh.rows = rows * (parabolas.size()-1); // not correct for multiple parabola tshirt torso, it is set again later
  _controlMesh.cols = cols;
  
  // We need to find (rows) number of interpolating parabolas between bottom and top parabola.
  // To do this we need the left and right intersection points for each interpolated parabola.
  // Once we have these ranges we need to calculate y at each beta along curve. Then calculate
  // Actual x from beta and range. This gives us a x,y co-ord which we can use to find the nearest 
  // texture coordinate.
  
  // We get problems with welding if we use b-values and the vertical height calculation
  // when spacing interpolated parabolas in cases where the edge is a strong curve (such
  // as arm holes of t-shirt. So for now distribute points regularly along boundary of
  // top edge when we have 3 parabolas
  
  std::vector<Vec3d> controlPoints; // Phil
  std::vector<bool> controlPointIsSeam; // JDW
  
  std::vector<Vec3d> patternPoints; // JDW
  std::vector<int> patternIndices; // JDW
  for(int lowerParabolaIndex=0;lowerParabolaIndex<parabolas.size()-1;lowerParabolaIndex++)
  {
    int upperParabolaIndex=lowerParabolaIndex+1;
    
    bool topSectionOfTshirt = false; // FIXME hack for Tshirt
    
    if((name.find("torso_") != std::string::npos))
    {
      if((upperParabolaIndex == parabolas.size()-1) && (upperParabolaIndex > 1))
      {
        topSectionOfTshirt=true;
        rows = TSHIRT_ARM_COLS;
      }
      else
      {
        // assume only two sections in Tshirt length
        rows = TSHIRT_TORSO_BOTTOM_ROWS;
      }
      _controlMesh.rows = TSHIRT_ARM_COLS + TSHIRT_TORSO_BOTTOM_ROWS;
      
    }
    
    std::vector<EdgeSection*>::const_iterator esi, esi_end;
    Parabola *pBottomParabola, *pTopParabola;
    pBottomParabola = parabolas.at(lowerParabolaIndex);
    pTopParabola = parabolas.at(upperParabolaIndex); 
    
    

    
    const int NUM_RADII = 100;
    double radii[NUM_RADII];
    for( int rIdx=0; rIdx<NUM_RADII; ++rIdx )
      radii[rIdx] = lengthOfInterpolatedParabola(pBottomParabola,pTopParabola, (double)rIdx/(double)(NUM_RADII-1), 50);
    double *b = new double[2*rows];
    _controlMesh.computeBValues(b, 2*rows, radii, NUM_RADII);
  
    
    // Special case for first and last parabola. Use edge sections to supply the co-ordinates to
    // ensure control mesh runs right to the edges of the garment
    
    // First calculate control points for the bottom row (easiest as we have parabola) (alpha is 0 here)
    if(lowerParabolaIndex==0)
    {
      for(int i=0; i<(2*cols)+1;i++) // so for 2 columns (2 patchs across) would want 5 points
      {
        double beta = (double) i / (double) (2*cols);
        double y = pBottomParabola->yBeta( beta );
        double x = pBottomParabola->xFromBeta( beta );
        long closestIndex=0;
        double minDist = 999999.9;
        // FIXME this isn't working yet...
        // Select closest from edge points // FIXME Add check for e.g. neck section
        std::vector<EdgeSection*>::const_iterator esi = _bottomEdgeSections.begin();
        std::vector<EdgeSection*>::const_iterator esi_end = _bottomEdgeSections.end();
        for(;esi!=esi_end;esi++) 
        {
          long index = (*esi)->closestTCIndexTo( x, y );
          
          Vec3d* pV = _lcmodel->tex_coords()[ index ];
          double distSq = (pV->x() - x)*(pV->x() - x) + (pV->y() - y) * (pV->y() - y);
          if(distSq < minDist)
          {
            minDist = distSq;
            closestIndex = index;
          }
        }
        
        //long index = findClosestTCTo( x, y );
        
        cp = *(_lcmodel->points()[closestIndex]);
        controlPoints.push_back( cp );
        patternPoints.push_back( *(_lcmodel->tex_coords()[closestIndex]) );
        interpolatedPatternPoints.push_back( *(_lcmodel->tex_coords()[closestIndex]) ); // not interpolated for edge
        patternIndices.push_back( closestIndex );
        if(_vertexList.at(closestIndex)->numberOfCoincidentPoints()>0)
        { // This point is a seam and should be welded in the control mesh
          controlPointIsSeam.push_back(true);
        }
        else
        {
          controlPointIsSeam.push_back(false);
        }
      }
    }
    else
    {/*intentionally empty*/} // This lower bound is somewhere in middle of garment 
      // However this row is generated as the top row of the previous iteration.
      // so do nothing
    
    
    // Now interpolate between bottom and top parabolas 
    
    for(int i = 1; i<(rows*2); i++)
    {
      
      
        double percentageOfHeight = 0.0;
        if(i > 0)
        {
          for(int z=i-1; z>=0;z--)
          {
            percentageOfHeight += b[z];
          }
        }
        
        if(topSectionOfTshirt)
        { // Force an even distribution in Y to avoid the effects of highly curved left and right arm seams
          percentageOfHeight = (double)i/(double)(rows*2);
        }
        
        
        
        double yLeft = Parabola::yInterpolate( *pBottomParabola, *pTopParabola, percentageOfHeight, 0.0 );
        double yRight = Parabola::yInterpolate( *pBottomParabola, *pTopParabola, percentageOfHeight, 1.0 );
        //double yRight = Parabola::yInterpolate( *pBottomParabola, *pTopParabola, (double) i / (rows*2.0), 1.0 );
        
        Vec3d leftInterceptPoint, rightInterceptPoint;
        
        double minDiff = 999999.9;
        
        if(topSectionOfTshirt)
        {
          std::cout << "Finding intersection points along top left of TShirt" << std::endl;
          // find left intercept
          // take the top edge section from the left
          sort( _leftEdgeSections.begin(), _leftEdgeSections.end(), esPtrYLess() );
          // now sort the edges in top EdgeSection, from least to greatest Y
          EdgeSection *pEdgeSection = _leftEdgeSections[ _leftEdgeSections.size() -1 ];
          sort( pEdgeSection->_edges.begin(), pEdgeSection->_edges.end(), edgePtrYLess() );
          
          // walk along the length of the edge section and pick points which are
          // evenly distributed
          int totalSegments = rows*2;
          double segmentLength = pEdgeSection->length / (double) totalSegments;
          double travelledLength = segmentLength * i;
          
          std::vector<polymesh::Edge*>::const_iterator ei;
          ei = pEdgeSection->_edges.begin();
          double lengthAccumulator = 0.0;
          for(;ei!= pEdgeSection->_edges.end();ei++)
          {
            polymesh::Edge* pEdge = *ei;
            lengthAccumulator += pEdge->length();
            if(lengthAccumulator >= travelledLength)
            {
                // we are travelling up the y direction
                // find the nearest edge point to the travelledLength
                
                // take lengthAccumulator, subtract this edge length, this gives length to first point in the edge
                // take the difference between this and the travelledLength. This is the distance to walk along this edge
                // Find the edge direction vector (lower point to higher point), find the point this distance along this vector
                // from the lower point. 
                
                // divide by length of edge. If < 0.5 then choose lower ( in Y ) point, otherwise higher
                
              double len_a = lengthAccumulator - pEdge->length();
              double diff = travelledLength - len_a;
              double ratio = diff / pEdge->length();
              Vec3d edgePointA = *(_lcmodel->tex_coords()[ pEdge->_vertexIndex[0] ]);
              Vec3d edgePointB = *(_lcmodel->tex_coords()[ pEdge->_vertexIndex[1] ]);
              if(edgePointA.y() < edgePointB.y())
              {
                if(ratio<=0.5) {
                  leftInterceptPoint = edgePointA;
                }
                else
                  leftInterceptPoint = edgePointB;
              }
              else
              {
                if(ratio<=0.5) {
                  leftInterceptPoint = edgePointB;
                }
                else
                  leftInterceptPoint = edgePointA;
              }
              break;
            }
          }

        } // end if topsection of tshirt
        else // bottom section
        {
          // Find left intercept
          esi = _leftEdgeSections.begin();
          esi_end = _leftEdgeSections.end();
          for(;esi!=esi_end;esi++)
          {
            std::vector<polymesh::Edge*>::const_iterator ei, ei_end;
            // iterate through each Edge in the left EdgeSection
            ei = (*esi)->_edges.begin();
            ei_end = (*esi)->_edges.end();
            for(;ei!=ei_end;ei++)
            {
              // Only examining one of the points in each edge to avoid duplicating effort
              // Will miss the termination point at one or other end of edgeSection but this
              // is the same as the left or right most point of either top or bottom curve
              // so shouldn't matter
              polymesh::Edge* pEdge = *ei;
              double edgeY = _lcmodel->tex_coords()[ pEdge->_vertexIndex[0] ]->y();
              double yDiff = fabs(edgeY - yLeft);
              if(yDiff < minDiff) {
                leftInterceptPoint = *(_lcmodel->tex_coords()[ pEdge->_vertexIndex[0] ]);
                minDiff = yDiff;
              }
            }
          } // end loop over all edge sections to find left intercept for non-tshirt top section
        } // end test for top tshirt section
        
        if(topSectionOfTshirt)
        {
          
          // find left intercept
          // take the top edge section from the left
            sort( _rightEdgeSections.begin(), _rightEdgeSections.end(), esPtrYLess() );
          // now sort the edges in top EdgeSection, from least to greatest Y
            EdgeSection *pEdgeSection = _rightEdgeSections[ _rightEdgeSections.size() -1 ];
            sort( pEdgeSection->_edges.begin(), pEdgeSection->_edges.end(), edgePtrYLess() );
          
          // walk along the length of the edge section and pick points which are
          // evenly distributed
            int totalSegments = rows*2;
            double segmentLength = pEdgeSection->length / (double) totalSegments;
            double travelledLength = segmentLength * i;
          
            std::vector<polymesh::Edge*>::const_iterator ei;
            ei = pEdgeSection->_edges.begin();
            double lengthAccumulator = 0.0;
            for(;ei!= pEdgeSection->_edges.end();ei++)
            {
              polymesh::Edge* pEdge = *ei;
              lengthAccumulator += pEdge->length();
              if(lengthAccumulator >= travelledLength)
              {
                // we are travelling up the y direction
                // find the nearest edge point to the travelledLength
                
                // take lengthAccumulator, subtract this edge length, this gives length to first point in the edge
                // take the difference between this and the travelledLength. This is the distance to walk along this edge
                // Find the edge direction vector (lower point to higher point), find the point this distance along this vector
                // from the lower point. 
                
                // divide by length of edge. If < 0.5 then choose lower ( in Y ) point, otherwise higher
                
                double len_a = lengthAccumulator - pEdge->length();
                double diff = travelledLength - len_a;
                double ratio = diff / pEdge->length();
                Vec3d edgePointA = *(_lcmodel->tex_coords()[ pEdge->_vertexIndex[0] ]);
                Vec3d edgePointB = *(_lcmodel->tex_coords()[ pEdge->_vertexIndex[1] ]);
                if(edgePointA.y() < edgePointB.y())
                {
                  if(ratio<=0.5) {
                    rightInterceptPoint = edgePointA;
                  }
                  else
                    rightInterceptPoint = edgePointB;
                }
                else
                {
                  if(ratio<=0.5) {
                    rightInterceptPoint = edgePointB;
                  }
                  else
                    rightInterceptPoint = edgePointA;
                }
                
                break;
              }
            }
        }
        else {
          // find right intercept
          minDiff = 999999.9;
          esi = _rightEdgeSections.begin();
          esi_end = _rightEdgeSections.end();
          for(;esi!=esi_end;esi++)
          {
            std::vector<polymesh::Edge*>::const_iterator ei, ei_end;
            ei = (*esi)->_edges.begin();
            ei_end = (*esi)->_edges.end();
            // iterate through each Edge in the right EdgeSection
            for(;ei!=ei_end;ei++)
            {
              // Only examining one of the points in each edge to avoid duplicating effort
              // Will miss the termination point at one or other end of edgeSection but this
              // is the same as the left or right most point of either top or bottom curve
              // so shouldn't matter
              polymesh::Edge* pEdge = *ei;
              double edgeY = _lcmodel->tex_coords()[ pEdge->_vertexIndex[0] ]->y();
              double yDiff = fabs(edgeY - yRight);
              if(yDiff < minDiff) {
                rightInterceptPoint = *(_lcmodel->tex_coords()[ pEdge->_vertexIndex[0] ]);
                minDiff = yDiff;
              }
              
            }
          }
          
    //      std::cout << "LeftIPoint: " << leftInterceptPoint <<std::endl;
      //    std::cout << "RightIPoint: " << rightInterceptPoint << std::endl <<std::endl;
        } // end finding intercept points
        
        double xRange = rightInterceptPoint.x() - leftInterceptPoint.x();
        for(int j=0; j<(2*cols)+1;j++) // so for 2 columns (2 patchs across) would want 5 points
        {
          double beta = (double) j / (double) (2*cols);
          //double y = Parabola::yInterpolate( *pBottomParabola, *pTopParabola, (double) i / (rows*2.0), beta );
          double y = Parabola::yInterpolate( *pBottomParabola, *pTopParabola, percentageOfHeight, beta );
          double x = leftInterceptPoint.x() + beta*xRange;
          
          
          
          //long index = findClosestTCTo( x, y );
          long index = findClosestTCToIgnoringDarts( x, y );
          // Find the triangles this vertex belongs to and calculate barycentric co-ordinates
          // b and c and then validate against a+b+c==1 (i.e. a=1-b-c).
          // If Each of a, b and c are within the range 0<=x<=1 then
          // assume point is inside triangle.
          
          Vec3d pA2d;
          Vec3d pB2d;
          Vec3d pC2d;
          double b,c;
          Vec3d P(x,y,0);
          Vec3d pA3d; 
          Vec3d pB3d; 
          Vec3d pC3d; 
          Vec3d interpolated3dPoint;
          
          polymesh::Vertex* pV = _vertexList.at( index );
          std::vector<polymesh::Face*>::const_iterator fi = pV->_faceList.begin();
          bool found = false;
          
          for(;fi!=pV->_faceList.end();fi++)
          {
            pA2d = *_lcmodel->tex_coords()[ (*fi)->_vertex[0] ]; 
            pB2d = *_lcmodel->tex_coords()[ (*fi)->_vertex[1] ];
            pC2d = *_lcmodel->tex_coords()[ (*fi)->_vertex[2] ];
            
            findBarycentricCoefficients(P, pA2d, pB2d, pC2d, &b, &c);
            if(b < 0.0 || b > 1.0) continue;
            if(c < 0.0 || c > 1.0) continue;
            
            double a = 1.0 - (b+c);
            
            if(a < 0.0 || a > 1.0) continue; // Some room for rounding
            found=true;
            
            
            pA3d = *_lcmodel->points()[ (*fi)->_vertex[0] ]; 
            pB3d = *_lcmodel->points()[ (*fi)->_vertex[1] ];
            pC3d = *_lcmodel->points()[ (*fi)->_vertex[2] ];
            
             // These three points are used to provide an interpolated position in 3d.
             interpolated3dPoint = pA3d + b*(pB3d - pA3d) + c*(pC3d - pA3d);
          
             cp = interpolated3dPoint;
             
            break;
          }
          
          if(!found) 
          {
            // it's outside the mesh. Just pick the nearest point
            cp = *(_lcmodel->points()[ index ]);
          }

          Vec3d interpPatternPoint = Vec3d(x,y,0); 
          // FIXME bit messy, keep index from first time around in tshirt case
          if(topSectionOfTshirt && j==0)
          {
            index = findClosestTCTo( leftInterceptPoint.x(), leftInterceptPoint.y() );
            cp = *(_lcmodel->points()[ index ]);
            interpPatternPoint = *(_lcmodel->tex_coords()[index]);
          }
          if(topSectionOfTshirt && j==(2*cols))
          {
            index = findClosestTCTo( rightInterceptPoint.x(), rightInterceptPoint.y() );
            cp = *(_lcmodel->points()[ index ]);
            interpPatternPoint = *(_lcmodel->tex_coords()[index]);
          }
          
          
          controlPoints.push_back( cp );
          patternPoints.push_back( *(_lcmodel->tex_coords()[index]) );
          interpolatedPatternPoints.push_back( interpPatternPoint );
          patternIndices.push_back( index );
          if(_vertexList.at(index)->numberOfCoincidentPoints()>0)
          { // This point is a seam and should be welded in the control mesh
            controlPointIsSeam.push_back(true);
          }
          else
          {
            controlPointIsSeam.push_back(false);
          }
        }
    } // end interpolation between bottom and top parabolas
    
      // Finally calculate control points for the top row (we have parabola) (alpha is 1 here)
    if(upperParabolaIndex==(parabolas.size()-1)) // We have EdgeSections to choose points from
    {
      for(int i=0; i<(2*cols)+1;i++) // so for 2 columns (2 patchs across) would want 5 points
      {
        double beta = (double) i / (double) (2*cols);
        double y = pTopParabola->yBeta( beta );
        double x = pTopParabola->xFromBeta( beta );
        // FIXME Don't snap to edges if in a neck hole for example. Compare with dartpoints.
        double minDist =9999999.9;
        long closestIndex=0;
        
        std::vector<EdgeSection*>::const_iterator esi = _topEdgeSections.begin();
        std::vector<EdgeSection*>::const_iterator esi_end = _topEdgeSections.end();
        for(;esi!=esi_end;esi++) 
        {
          long index = (*esi)->closestTCIndexTo( x, y );
          Vec3d* pV = _lcmodel->tex_coords()[ index ];
          double distSq = (pV->x() - x)*(pV->x() - x) + (pV->y() - y) * (pV->y() - y);
          if(distSq < minDist)
          {
            minDist = distSq;
            closestIndex = index;
          }
        }
        
        //long index = findClosestTCTo( x, y );
        cp = *(_lcmodel->points()[closestIndex]);
        controlPoints.push_back( cp );
        patternPoints.push_back( *(_lcmodel->tex_coords()[closestIndex]) );
        interpolatedPatternPoints.push_back( *(_lcmodel->tex_coords()[closestIndex]) );
        patternIndices.push_back( closestIndex );
        if(_vertexList.at(closestIndex)->numberOfCoincidentPoints()>0)
        { // This point is a seam and should be welded in the control mesh
          controlPointIsSeam.push_back(true);
        }
        else
        {
          controlPointIsSeam.push_back(false);
        }
      }
    }
    else // this top row is in the middle of the garment
      // we need to choose points from the whole garment and we have the left-right intercept values from the parabola
    {
      for(int j=0; j<(2*cols)+1;j++) // so for 2 columns (2 patchs across) would want 5 points
      {
        double beta = (double) j / (double) (2*cols);
        double y = pTopParabola->yBeta( beta );
        double x = pTopParabola->xFromBeta( beta );
          
        //long index = findClosestTCTo( x, y );
        
        // start paste
        
        long index = findClosestTCToIgnoringDarts( x, y );
          // Find the triangles this vertex belongs to and calculate barycentric co-ordinates
          // b and c and then validate against a+b+c==1 (i.e. a=1-b-c).
          // If Each of a, b and c are within the range 0<=x<=1 then
          // assume point is inside triangle.
          
        Vec3d pA2d;
        Vec3d pB2d;
        Vec3d pC2d;
        double b,c;
        Vec3d P(x,y,0);
        Vec3d pA3d; 
        Vec3d pB3d; 
        Vec3d pC3d; 
        Vec3d interpolated3dPoint;
        Vec3d interpolated2dPoint = Vec3d(x,y,0);
          
        polymesh::Vertex* pV = _vertexList.at( index );
        std::vector<polymesh::Face*>::const_iterator fi = pV->_faceList.begin();
        bool found = false;
          
        for(;fi!=pV->_faceList.end();fi++)
        {
          pA2d = *_lcmodel->tex_coords()[ (*fi)->_vertex[0] ]; 
          pB2d = *_lcmodel->tex_coords()[ (*fi)->_vertex[1] ];
          pC2d = *_lcmodel->tex_coords()[ (*fi)->_vertex[2] ];
            
          findBarycentricCoefficients(P, pA2d, pB2d, pC2d, &b, &c);
          if(b < 0.0 || b > 1.0) continue;
          if(c < 0.0 || c > 1.0) continue;
            
          double a = 1.0 - (b+c);
            
          if(a < 0.0 || a > 1.0) continue; // Some room for rounding
          found=true;
            
            
          pA3d = *_lcmodel->points()[ (*fi)->_vertex[0] ]; 
          pB3d = *_lcmodel->points()[ (*fi)->_vertex[1] ];
          pC3d = *_lcmodel->points()[ (*fi)->_vertex[2] ];
            
             // These three points are used to provide an interpolated position in 3d.
          interpolated3dPoint = pA3d + b*(pB3d - pA3d) + c*(pC3d - pA3d);
          
          cp = interpolated3dPoint;
             
          break;
        }
          
        if(!found) 
        {
            // it's outside the mesh. Just pick the nearest point
          cp = *(_lcmodel->points()[ index ]);
          interpolated2dPoint = *(_lcmodel->tex_coords()[ index ]);
        }
        
        
        
        
        // end paste

        controlPoints.push_back( cp );
        patternPoints.push_back( *(_lcmodel->tex_coords()[index]) );
        interpolatedPatternPoints.push_back( interpolated2dPoint );
        patternIndices.push_back( index);
        if(_vertexList.at(index)->numberOfCoincidentPoints()>0)
        { // This point is a seam and should be welded in the control mesh
          controlPointIsSeam.push_back(true);
        }
        else
        {
          controlPointIsSeam.push_back(false);
        }
      }
    }
    
    
    
  
  } // end loop over parabolas
  
  assert( (int)controlPoints.size()== ((cols*2+1) * (_controlMesh.rows*2+1)) );
    
    /* Phil
    _controlMesh.patchs.clear();
    _controlMesh.patchs.resize(0);
    ControlMesh::Patch p;
    */
    _controlMesh.reset();	// Phil
    Vec3d cpPos[9];		// Phil
    double cpEps = 0.00001 * (_lcmodel->bbMax()-_lcmodel->bbMin()).norm();	// Phil
    
    int initialWeldIndex=2;
    int currentWeldIndex[2];
    currentWeldIndex[0]=initialWeldIndex; // front left edge
    currentWeldIndex[1]=initialWeldIndex + (_controlMesh.rows*2+1); // front right edge
    
    // for tshirt the seams bifurcate as you approach the armpits to accommodate the sleeves 
    // so we need another 4 sets of unique indices
    
    int initialArmPitsIndices[4];
    
    // Note, Tshirt body and upper chest can have different row counts. Correct total is in _controlMesh.rows
    int torsoEdgePointCount;
    int backSleeveAdditionalPointCount;
    int frontLeftMiddleIndex;
    int frontRightMiddleIndex;
    if((name.find("torso_") != std::string::npos) )
    {
      torsoEdgePointCount = (_controlMesh.rows *2) +1; 
      backSleeveAdditionalPointCount = TSHIRT_ARM_COLS*2 -1; // arm cols is top section rows
      frontLeftMiddleIndex = TSHIRT_TORSO_BOTTOM_ROWS*2;
      frontRightMiddleIndex= frontLeftMiddleIndex + torsoEdgePointCount;
    }
    if(name.find("_arm_") != std::string::npos)
    {
      torsoEdgePointCount = ((TSHIRT_ARM_COLS+TSHIRT_TORSO_BOTTOM_ROWS) *2) +1; // needed in the arm atlases (remember arm cols == torso rows)
      backSleeveAdditionalPointCount = TSHIRT_ARM_COLS*2 -1; // arm cols is top section rows
      frontLeftMiddleIndex = TSHIRT_TORSO_BOTTOM_ROWS*2;
      frontRightMiddleIndex = frontLeftMiddleIndex + torsoEdgePointCount;
    }

    
    initialArmPitsIndices[0] = initialWeldIndex + frontLeftMiddleIndex; // left front
    initialArmPitsIndices[1] = initialWeldIndex + frontRightMiddleIndex; // right front
    initialArmPitsIndices[2] = initialWeldIndex + 2*torsoEdgePointCount + backSleeveAdditionalPointCount; // left torso back (same as LAB)
    initialArmPitsIndices[3] = initialWeldIndex + 2*torsoEdgePointCount ; // right back

    
    for( j=0; j<_controlMesh.rows; ++j )
    {
      for( i=0; i<cols; ++i )
      {
        double cpLength[6];
        int index[9];
        int atlasIndices[9];
        Vec3d interpPatternPoints[9];
        index[0] = 2*i             + (j*2)*(cols*2+1);
        index[1] = 2*i+1           + (j*2)*(cols*2+1);
        index[2] = 2*i+2          + (j*2)*(cols*2+1);
        index[3] = 2*i             + (j*2+1)*(cols*2+1);
        index[4] = 2*i+1           + (j*2+1)*(cols*2+1);
        index[5] = 2*i+2           + (j*2+1)*(cols*2+1);
        index[6] = 2*i             + (j*2+2)*(cols*2+1);
        index[7] = 2*i+1           + (j*2+2)*(cols*2+1);
        index[8] = 2*i+2           + (j*2+2)*(cols*2+1);
        
        for(int ii=0;ii<9;ii++)
          atlasIndices[ii] = patternIndices[ index[ii] ];
        
        for(int ii=0;ii<9;ii++)
          interpPatternPoints[ii] = interpolatedPatternPoints[ index[ii] ];
        
  
        // NOTE this are ignored later as the controlmesh recalculates them
        /*
        cpLength[0] =  (patternPoints[index[1]]-patternPoints[index[0]]).norm()
            + (patternPoints[index[2]]-patternPoints[index[1]]).norm();
        cpLength[1] =  (patternPoints[index[6]]-patternPoints[index[7]]).norm()
            + (patternPoints[index[8]]-patternPoints[index[7]]).norm();
        cpLength[2] =  (patternPoints[index[6]]-patternPoints[index[3]]).norm()
            + (patternPoints[index[3]]-patternPoints[index[0]]).norm();
        cpLength[3] =  (patternPoints[index[8]]-patternPoints[index[5]]).norm()
            + (patternPoints[index[5]]-patternPoints[index[2]]).norm();
        cpLength[4] =  (patternPoints[index[6]]-patternPoints[index[4]]).norm()
            + (patternPoints[index[4]]-patternPoints[index[2]]).norm();
        cpLength[5] =  (patternPoints[index[8]]-patternPoints[index[4]]).norm()
            + (patternPoints[index[4]]-patternPoints[index[0]]).norm();
        */

        
        cpPos[0] = controlPoints[ 2*i             + (j*2)*(cols*2+1) ];
        cpPos[1] = controlPoints[ 2*i+1           + (j*2)*(cols*2+1) ];
        cpPos[2] = controlPoints[ 2*i+2           + (j*2)*(cols*2+1) ];
              
        cpPos[3] = controlPoints[ 2*i             + (j*2+1)*(cols*2+1) ];
        cpPos[4] = controlPoints[ 2*i+1           + (j*2+1)*(cols*2+1) ];
        cpPos[5] = controlPoints[ 2*i+2           + (j*2+1)*(cols*2+1) ];
              
        cpPos[6] = controlPoints[ 2*i             + (j*2+2)*(cols*2+1) ];
        cpPos[7] = controlPoints[ 2*i+1           + (j*2+2)*(cols*2+1) ];
        cpPos[8] = controlPoints[ 2*i+2           + (j*2+2)*(cols*2+1) ];
        
        int cpIsSeam[9] = { 0,0,0, 0,0,0, 0,0,0 };
        
        cpIsSeam[0] = controlPointIsSeam[ 2*i             + (j*2)*(cols*2+1) ] ? 1:0;
        cpIsSeam[1] = controlPointIsSeam[ 2*i+1           + (j*2)*(cols*2+1) ]? 1:0;
        cpIsSeam[2] = controlPointIsSeam[ 2*i+2           + (j*2)*(cols*2+1) ]? 1:0;
              
        cpIsSeam[3] = controlPointIsSeam[ 2*i             + (j*2+1)*(cols*2+1) ]? 1:0;
        cpIsSeam[4] = controlPointIsSeam[ 2*i+1           + (j*2+1)*(cols*2+1) ]? 1:0;
        cpIsSeam[5] = controlPointIsSeam[ 2*i+2           + (j*2+1)*(cols*2+1) ]? 1:0;
              
        cpIsSeam[6] = controlPointIsSeam[ 2*i             + (j*2+2)*(cols*2+1) ]? 1:0;
        cpIsSeam[7] = controlPointIsSeam[ 2*i+1           + (j*2+2)*(cols*2+1) ]? 1:0;
        cpIsSeam[8] = controlPointIsSeam[ 2*i+2           + (j*2+2)*(cols*2+1) ]? 1:0;
        
        // FIXME Now force seam welding for particular design.
        if(name.find("skirt_front") != std::string::npos)
        { // SKIRT FRONT. 
          // Front right edge to back left edge
          
          if( 0 == i ) // Left edge
          {
            cpIsSeam[0] = currentWeldIndex[0]++;
            cpIsSeam[3] = currentWeldIndex[0]++;
            cpIsSeam[6] = currentWeldIndex[0];// current patch top left  index is next bottom right index
            
          }
          
          if( (cols-1) == i ) // Right edge
          {
            cpIsSeam[2] = currentWeldIndex[1]++;
            cpIsSeam[5] = currentWeldIndex[1]++;
            cpIsSeam[8] = currentWeldIndex[1]; // current patch top right index is next bottom right index
            
          }
        }
        
        if(name.find("skirt_back") != std::string::npos)
        { // SKIRT FRONT. 
         
          
          if( 0 == i ) // Back left edge corresponds to front right edge
          {
            cpIsSeam[0] = currentWeldIndex[1]++;
            cpIsSeam[3] = currentWeldIndex[1]++;
            cpIsSeam[6] = currentWeldIndex[1];// current patch top left index is next bottom right index
            
          }
          
          if( (cols-1) == i ) // Back right edge corresponds to front left edge
          {
            cpIsSeam[2] = currentWeldIndex[0]++;
            cpIsSeam[5] = currentWeldIndex[0]++;
            cpIsSeam[8] = currentWeldIndex[0]; // current patch top right index is next bottom right index
            
          }
        }
        

        if(name.find("torso_front") != std::string::npos)
        { 
          //if(j < (parabolas.size()-1)*rows/2) // only the bottom half
          {
            // Front right edge to back left edge
            if( 0 == i ) // Left edge
            {
              //std::cout << "tf l weld index: " << currentWeldIndex[0] << std::endl;
              cpIsSeam[0] = currentWeldIndex[0]++;
              //std::cout << "tf l weld index: " << currentWeldIndex[0] << std::endl;
              cpIsSeam[3] = currentWeldIndex[0]++;
              //std::cout << "tf l weld index: " << currentWeldIndex[0] << std::endl;
              cpIsSeam[6] = currentWeldIndex[0];// current patch top left  index is next bottom right index
              
              // Hack set params for push back
              // This is the top most row
              if(j == _controlMesh.rows-1){
                topTorsoLeft = cpPos[6];
              }
            }
            
            if( (cols-1) == i ) // Right edge
            {
              //std::cout << "tf r weld index: " << currentWeldIndex[1] << std::endl;
              cpIsSeam[2] = currentWeldIndex[1]++;
              //std::cout << "tf r weld index: " << currentWeldIndex[1] << std::endl;
              cpIsSeam[5] = currentWeldIndex[1]++;
              //std::cout << "tf r weld index: " << currentWeldIndex[1] << std::endl;
              cpIsSeam[8] = currentWeldIndex[1]; // current patch top right index is next bottom right index
            }
          }
          
        }
        
        if(name.find("torso_back") != std::string::npos)
        {  
          
          if(j < TSHIRT_TORSO_BOTTOM_ROWS) // only the bottom half
          {
            if( 0 == i ) // Back left edge corresponds to front right edge
            {
              //std::cout << "tb l weld index: " << currentWeldIndex[1] << std::endl;
              cpIsSeam[0] = currentWeldIndex[1]++;
              //std::cout << "tb l weld index: " << currentWeldIndex[1] << std::endl;
              cpIsSeam[3] = currentWeldIndex[1]++;
              //std::cout << "tb l weld index: " << currentWeldIndex[1] << std::endl;
              cpIsSeam[6] = currentWeldIndex[1];// current patch top left index is next bottom right index
              
            }
            
            if( (cols-1) == i ) // Back right edge corresponds to front left edge
            {
              //std::cout << "tb r weld index: " << currentWeldIndex[0] << std::endl;
              cpIsSeam[2] = currentWeldIndex[0]++;
              //std::cout << "tb r weld index: " << currentWeldIndex[0] << std::endl;
              cpIsSeam[5] = currentWeldIndex[0]++;
              //std::cout << "tb r weld index: " << currentWeldIndex[0] << std::endl;
              cpIsSeam[8] = currentWeldIndex[0]; // current patch top right index is next bottom right index
            }
            
          }
          else // Top half of TShirt
          {
            if( 0 == i ) // Left edge
            {
              if(j == TSHIRT_TORSO_BOTTOM_ROWS) // the critical middle row
              {
                cpIsSeam[0] = currentWeldIndex[1]++;
              }
              else
              {
                //std::cout << "tb l topsect weld index: " << initialArmPitsIndices[2] << std::endl;
                cpIsSeam[0] = initialArmPitsIndices[2]++;
              }
              
              //std::cout << "tb l topsect weld index: " << initialArmPitsIndices[2] << std::endl;
              cpIsSeam[3] = initialArmPitsIndices[2]++;
              
              
              if(j==_controlMesh.rows-1)
              { // top most point must weld front atlas again
                cpIsSeam[6] = initialArmPitsIndices[2] + torsoEdgePointCount-1;
              }
              else
              {
                //std::cout << "tb l topsect weld index: " << initialArmPitsIndices[2] << std::endl;
                cpIsSeam[6] = initialArmPitsIndices[2]; // don't inc
              }
              
            }
            
            if( (cols-1) == i ) // Right edge
            {
              if(j == TSHIRT_TORSO_BOTTOM_ROWS) // the critical middle row
              {
                cpIsSeam[2] = currentWeldIndex[0]++;
              }
              else
                cpIsSeam[2] = initialArmPitsIndices[3]++;
              
              
              //std::cout << "tb r topsect weld index: " << initialArmPitsIndices[3] << std::endl;
              cpIsSeam[5] = initialArmPitsIndices[3]++;
              //std::cout << "tb r topsect weld index: " << initialArmPitsIndices[3] << std::endl;
              
              if(j==_controlMesh.rows-1)
              {// top most point must weld front atlas again
                initialArmPitsIndices[0] + torsoEdgePointCount-1;
              }
              else
              {
                cpIsSeam[8] =initialArmPitsIndices[3];// don't inc
              }
              
            }
          }
        }
        
        if(name.find("left_arm_front") != std::string::npos)
        {
          int raOffset = 1000; // FIXME Define this according to the indices used up already
          if(i==0 && j==0)
          {
            // FIXME messy, only offset once
            currentWeldIndex[0] += raOffset;
            currentWeldIndex[1] += raOffset;
          }
          // top row is critical as it welds to front torso atlas as well as back right arm atlas in corners
          if(j==_controlMesh.rows-1) // TOP ROW
          {
            //we need an array corresponding to the indices used for the top right section of the torso front
            int *indices = new int[ cols*2+1 ];
            
            //int torsoVerticalEdgeCount = ((pointsAlongVerticalEdge-1)*2) + 1; // double because for tshirt two sections in torso
            
            //int halfwayIndex = ((torsoVerticalEdgeCount-1)/2) + 1;
            int ii = initialWeldIndex + frontRightMiddleIndex;
            for(int aa=0;aa<cols*2+1;aa++)
            {
              indices[aa] = ii++;
              //std::cout << "LAF indices[" << aa << "]: " << indices[aa] << std::endl;
            }
            
            
            if( 0 == i ) // Left edge
            {
              cpIsSeam[0] = currentWeldIndex[0]++;
              cpIsSeam[3] = currentWeldIndex[0]++;
              
            }
            
            cpIsSeam[6] = indices[ (i*2) ];
            cpIsSeam[7] = indices[ (i*2+1) ];
            cpIsSeam[8] = indices[ (i*2+2) ];
            
            
            if( (cols-1) == i ) // Right edge
            {
              cpIsSeam[2] = currentWeldIndex[1]++;
              cpIsSeam[5] = currentWeldIndex[1]++;
            }
          }
          else
          {
            
            if( 0 == i ) // Left edge
            {
              cpIsSeam[0] = currentWeldIndex[0]++;
              cpIsSeam[3] = currentWeldIndex[0]++;
              cpIsSeam[6] = currentWeldIndex[0];// current patch top left  index is next bottom right index
            }
            
            if( (cols-1) == i ) // Right edge
            {
              cpIsSeam[2] = currentWeldIndex[1]++;
              cpIsSeam[5] = currentWeldIndex[1]++;
              cpIsSeam[8] = currentWeldIndex[1]; // current patch top right index is next bottom right index
            }
          }
        }
        
        if(name.find("left_arm_back") != std::string::npos)
        { 
          int raOffset = 1000; // FIXME Define this according to the indices used up already
          if(i==0 && j==0)
          {
            currentWeldIndex[0] += raOffset;
            currentWeldIndex[1] += raOffset;
          }
          
          
          int *indices = new int[ cols*2+1 ];
            
          //int torsoVerticalEdgeCount = ((pointsAlongVerticalEdge-1)*2) + 1; // double because for tshirt two sections in torso
          //int halfwayIndex = ((torsoVerticalEdgeCount-1)/2) + 1;
          indices[cols*2] = initialWeldIndex + frontRightMiddleIndex;
          indices[0] = initialWeldIndex + 2*torsoEdgePointCount -1;
          
          
             
          //std::cout << "LAB indices[cols*2]: " << indices[cols*2] << std::endl;
          int ii = initialWeldIndex + 2*torsoEdgePointCount + backSleeveAdditionalPointCount; 
          // just the middle values
          for(int aa=cols*2-1;aa>0;aa--)
          {
            indices[aa] = ii++;
            //std::cout << "LAB indices[" << aa << "]: " << indices[aa] << std::endl;
          }
          //std::cout << "LAB indices[0]: " << indices[0] << std::endl;
          
          // top row is critical as it welds to front torso atlas as well as back right arm atlas in corners
          if(j==_controlMesh.rows-1) // TOP ROW
          {
            if( 0 == i ) // Left edge
            {
              cpIsSeam[0] = currentWeldIndex[1]++;
              cpIsSeam[3] = currentWeldIndex[1]++;
              
              topTorsoRight = cpPos[6];
            }
            
            cpIsSeam[6] = indices[ (i*2) ];
            cpIsSeam[7] = indices[ (i*2+1) ];
            cpIsSeam[8] = indices[ (i*2+2) ];
            
            
            if( (cols-1) == i ) // Right edge
            {
              // link to left edge of arm front
              cpIsSeam[2] = currentWeldIndex[0]++;
              cpIsSeam[5] = currentWeldIndex[0]++;
              
              
              middleTorsoRight = cpPos[8];
              
              
            }
          }
          else
          {
            if( 0 == i ) // Left edge
            {
              cpIsSeam[0] = currentWeldIndex[1]++;
              cpIsSeam[3] = currentWeldIndex[1]++;
              cpIsSeam[6] = currentWeldIndex[1];// current patch top left  index is next bottom right index
              
              if(j==0) // bottom row
                leftCuffTop = cpPos[0];
            }
            
            if( (cols-1) == i ) // Right edge
            {
              // link to left edge of arm front
              cpIsSeam[2] = currentWeldIndex[0]++;
              cpIsSeam[5] = currentWeldIndex[0]++;
              cpIsSeam[8] = currentWeldIndex[0]; // current patch top right index is next bottom right index
              
              if(j==0) // bottom row
                leftCuffBottom = cpPos[2];
            }
          }
        }
        
        if(name.find("right_arm_front") != std::string::npos)
        {
          int raOffset = 2000; // FIXME Define this according to the indices used up already
          if(i==0 && j==0)
          {
            currentWeldIndex[0] += raOffset;
            currentWeldIndex[1] += raOffset;
            
            // hack for pushback
              rightCuffTop = cpPos[0];
            
            
          }
          
          if( ((cols-1) == i) && (j==0))
            // bottom of front right cuff
          {
            rightCuffBottom = cpPos[2];
          }
          
          int *indices = new int[ cols*2+1 ];
            
          //int torsoVerticalEdgeCount = ((pointsAlongVerticalEdge-1)*2) + 1; // double because for tshirt two sections in torso
          
          //indices[0] = initialWeldIndex + torsoVerticalEdgeCount -1; // TL of front torso
          //indices[cols*2] = initialWeldIndex + halfwayIndex; // ML of front torso
          
          
          int ii = initialWeldIndex + frontLeftMiddleIndex; // torso back right side
          for(int aa=cols*2;aa>=0;aa--)
          {
            indices[aa] = ii++;
            //std::cout << "RAF indices[" << aa << "]: " << indices[aa] << std::endl;
          }
          
          // top row is critical as it welds to front torso atlas as well as back right arm atlas in corners
          if(j==_controlMesh.rows-1) // TOP ROW
          {
            if( 0 == i ) // Left edge
            {
              cpIsSeam[0] = currentWeldIndex[0]++;
              cpIsSeam[3] = currentWeldIndex[0]++;

            }
            
            
            
            cpIsSeam[6] = indices[ (i*2) ];
            cpIsSeam[7] = indices[ (i*2+1) ];
            cpIsSeam[8] = indices[ (i*2+2) ];
                
            if( (cols-1) == i ) // Right edge
            {
              cpIsSeam[2] = currentWeldIndex[1]++;
              cpIsSeam[5] = currentWeldIndex[1]++;
              
              // Hack for push back
              middleTorsoLeft = cpPos[8];

            }
          }
          else
          {
            // Front right edge to back left edge
            
            if( 0 == i ) // Left edge
            {
              cpIsSeam[0] = currentWeldIndex[0]++;
              cpIsSeam[3] = currentWeldIndex[0]++;
              cpIsSeam[6] = currentWeldIndex[0];// current patch top left  index is next bottom right index
            }
            
            if( (cols-1) == i ) // Right edge
            {
              cpIsSeam[2] = currentWeldIndex[1]++;
              cpIsSeam[5] = currentWeldIndex[1]++;
              cpIsSeam[8] = currentWeldIndex[1]; // current patch top right index is next bottom right index
            }
          }
        }
        
        if(name.find("right_arm_back") != std::string::npos)
        {
          int raOffset = 2000; // FIXME Define this according to the indices used up already
          if(i==0 && j==0)
          {
            currentWeldIndex[0] += raOffset;
            currentWeldIndex[1] += raOffset;
          }
          
          int *indices = new int[ cols*2+1 ];
            
          //int torsoVerticalEdgeCount = ((pointsAlongVerticalEdge-1)*2) + 1; // double because for tshirt two sections in torso
          //int halfwayIndex = ((torsoVerticalEdgeCount-1)/2) + 1;
        
          //indices[0] = initialWeldIndex + halfwayIndex - 1;
          //indices[cols*2] = initialWeldIndex + torsoVerticalEdgeCount -1;
          indices[0] = initialWeldIndex + frontLeftMiddleIndex;
          indices[cols*2] = initialWeldIndex + torsoEdgePointCount -1;
          
          int ii = initialWeldIndex + 2*torsoEdgePointCount;
          
          //std::cout << "RAB indices[0]: " << indices[0] << std::endl;
          for(int aa=1;aa<cols*2;aa++)
          {
            indices[aa] = ii++;
            //std::cout << "RAB indices[" << aa << "]: " << indices[aa] << std::endl;
          }
          //std::cout << "RAB indices[cols*2]: " << indices[cols*2] << std::endl;
          
          // top row is critical as it welds to front torso atlas as well as back right arm atlas in corners
          if(j==_controlMesh.rows-1) // TOP ROW
          {
            if( 0 == i ) // Left edge
            {
              cpIsSeam[0] = currentWeldIndex[1]++;
              cpIsSeam[3] = currentWeldIndex[1]++;
            }
            
            cpIsSeam[6] = indices[ (i*2) ];
            cpIsSeam[7] = indices[ (i*2+1) ];
            cpIsSeam[8] = indices[ (i*2+2) ];
            
            if( (cols-1) == i ) // Right edge
            {
              // link to left edge of arm front
              cpIsSeam[2] = currentWeldIndex[0]++;
              cpIsSeam[5] = currentWeldIndex[0]++;
            }
          }
          else
          {
            if( 0 == i ) // Left edge
            {
              cpIsSeam[0] = currentWeldIndex[1]++;
              cpIsSeam[3] = currentWeldIndex[1]++;
              cpIsSeam[6] = currentWeldIndex[1];// current patch top left  index is next bottom right index
            }
            
            if( (cols-1) == i ) // Right edge
            {
              // link to left edge of arm front
              cpIsSeam[2] = currentWeldIndex[0]++;
              cpIsSeam[5] = currentWeldIndex[0]++;
              cpIsSeam[8] = currentWeldIndex[0]; // current patch top right index is next bottom right index
            }
          }
        }
        
        // 2D length calculation WITHOUT accounting for dart space 
        //_controlMesh.addPatch(cpPos, cpEps, cpIsSeam, 0, cpLength, NULL, NULL, atlasIndices, name);
        // let ControlMesh calculate the lengths from 3d
        double weldRadius=(cpPos[0]-cpPos[8]).norm() *4;
              //_averageMeshVectorLength*4; // set large for Tshirt (default) as it has unique weld indices
        if((name.find("skirt_front") != std::string::npos) || (name.find("skirt_back") != std::string::npos))
        {
          std::cout << "Using skirt settings for atlas WELD RADIUS" << std::endl;
          weldRadius = (cpPos[0]-cpPos[8]).norm();
              //_averageMeshVectorLength*2; // lower for the non-unique skirt seams
        }
        
        _controlMesh.addPatch(cpPos, cpEps, cpIsSeam, weldRadius, NULL, NULL, -1, atlasIndices, name, interpPatternPoints);
      }
    }
   
  assert( _controlMesh.getNumPatchs()==_controlMesh.rows*cols );
  
  PushBackParams *p = &(this->parent->parent->_push_back_params_right);
  PushBackParams *pleft = &(this->parent->parent->_push_back_params_left);
  
  p->startPoint = middleTorsoLeft + (topTorsoLeft - middleTorsoLeft) /2;
  p->endPoint = rightCuffBottom + (rightCuffTop- rightCuffBottom) /2;
  
  p->endPoint += 0.2 * (p->endPoint - p->startPoint);  
  double len = (p->startPoint - p->endPoint).norm();
  
  p->shrunkLength = 0.9 * len;
  p->a = 0.3 * len;
  p->b = 0.4*len;
  p->influenceRadius = (topTorsoLeft - middleTorsoLeft).norm()*2/3; // test value
  p->valid = true;
  
   
  // LEFT ARM PUSH BACK PARAMS
  
  pleft->startPoint = middleTorsoRight + (topTorsoRight - middleTorsoRight) /2;
  pleft->endPoint = leftCuffBottom + (leftCuffTop- leftCuffBottom) /2;
  
   pleft->endPoint += 0.2 * (pleft->endPoint - pleft->startPoint);  
   len = (pleft->startPoint - pleft->endPoint).norm();
  
   pleft->shrunkLength = 0.8 * len;
   pleft->a = 0.05 * len;
   pleft->b = 0.25 * len;
   pleft->influenceRadius = (topTorsoRight - middleTorsoRight).norm()*2/3; // test value
   pleft->valid = true;
   
   
  _controlMesh.update();
  _controlMeshAssigned = true;
  std::cout << "Control mesh updated for " << name << std::endl;
}

bool GarmentAtlas::addFace(polymesh::Face* pFace) {
  if(!pFace) {
    std::cerr << "addFace: passed null Face pointer" << std::endl;
    return false;
  }
  _faceList.push_back(pFace);
  return true;
}

// Lookup an Edge by it's two endpoints. Order doesn't matter.
polymesh::Edge* GarmentAtlas::lookup(long vi1,long vi2) {
  // Assert not the same or negative
  if(vi1 < 0 || vi2 < 0) {
    std::cerr << "Passed negative index" << std::endl;
    return (polymesh::Edge*)0;
  }
  if(vi1 == vi2) {
    std::cerr << "Cannot lookup edge with same start and end vertex" << std::endl;
    return (polymesh::Edge*)0;
  }
  // Order indexes
  if(vi1 > vi2) { long _tmp = vi1; vi1 = vi2; vi2 = _tmp; }
  
   
  //std::cout <<"Lookup using: " << key << std::endl;
  std::pair<long, long> key(vi1, vi2);
  
  std::map<std::pair<long, long>, polymesh::Edge*>::iterator p = _edgeLookupMap.find( key );
  if(p == _edgeLookupMap.end()) {
    //std::cout << "Could not find edge" << std::endl;
    return (polymesh::Edge*) 0;
  }
  //std::cout << "Found edge: " << (*p).first << std::endl;
  return (*p).second;
}

GarmentAtlas::VertexPList GarmentAtlas::vertices() {
  return _vertexList;
}

GarmentAtlas::FacePList GarmentAtlas::faces() {
  return _faceList;
}
