//
// C++ Interface: EdgeSection
//
// Description: 
//
//
// Author: Jamie Wither <wither@stalactite>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef  EDGE_SECTION_H
#define EDGE_SECTION_H

#include <vector>
#include "vectypes.h"
#include "Edge.h"

class LCModel;

class EdgeSection {
  
  public:
    
    typedef enum { UP, DOWN, LEFT, RIGHT, UNKNOWN } Orientation;
    typedef enum { TOP, BOTTOM, LEFT_SIDE, RIGHT_SIDE, UNCLASSIFIED } Classification;
    EdgeSection() ;
    ~EdgeSection() {};
    Orientation _orientation;
    Classification _classification;
    std::vector<polymesh::Edge*> _edges;
    Vec3d averagePos;
    Vec3d averageNormal;
    double length;
    bool isDartSeam;
    
    EdgeSection* next;
    EdgeSection* previous;
    
    long closestTCIndexTo( double x, double y ) const;
    
    void merge(EdgeSection*); // Merge supplied ES with this edgesection
    void recalculateAttributes(LCModel*);

};

#endif // EDGE_SECTION_H
