//
//  Filename         : segment.cpp
//  Author           : Emmanuel Turquin
//  Purpose          : Class representing a 2d segment for the Canvas.
//  Date of creation : 04/26/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include "grid.h"
#include "segment.h"

Segment::~Segment() {
  CellList::iterator it = _cells.begin();
  CellList::iterator it_end = _cells.end();
  for (; it != it_end; ++it) {
    (*it)->removeSegment(this);
  }
}

Segment Segment::clippedToRect(const Vec2d tl, const Vec2d br)
{
  // use type to indictate whether they intersect at all
  Segment result;
  result.setType( OUT );
  counted_ptr<Point> pa(new Point(this->pointA()->x(), this->pointA()->y()));
  counted_ptr<Point> pb(new Point(this->pointB()->x(), this->pointB()->y()));
  result.setPointA(pa);
  result.setPointB(pb);
  
  //cout << result.pointA()->x() << endl;
  
  Vec2d dir( this->pointB()->x() - this->pointA()->x(),
             this->pointB()->y() - this->pointA()->y() );
  dir.normalizeSafe();
  
  bool reclip = true;
  int clipCounter = 12; // avoid problems when we don't get convergence //TODO better clip algorithm here
  
  while(reclip && clipCounter) {
    reclip = false;
    clipCounter--;
  // clip PointA to tl
  if( result.pointA()->x() < tl[0]) 
  {
    if( result.pointB()->x() < tl[0])
    {
      return result;
    }
    else {
      // clip point A.x
      double dx = tl[0] - result.pointA()->x();
      double ratio = dx / dir.x();
      double dy = ratio * dir.y();
      result.pointA()->setX( tl[0] );
      result.pointA()->setY( result.pointA()->y() + dy );
      reclip = true;
    }
  }
  
  if( result.pointA()->y() > tl[1]) 
  {
    if( result.pointB()->y() > tl[1])
    {
      return result;
    }
    else {
      // clip point A.y
      
      double dy = tl[1] - result.pointA()->y(); 
      double ratio = dy / dir.y(); 
      double dx = ratio * dir.x();
      
      result.pointA()->setX( result.pointA()->x() + dx );
      result.pointA()->setY( tl[1] );
      reclip = true;
      
      
    }
  }
  
  // clip PointA to br
  if( result.pointA()->x() > br[0]) 
  {
    if( result.pointB()->x() > br[0])
    {
      return result;
    }
    else {
      // clip point A.x
      double dx = br[0] - result.pointA()->x();
      double ratio = dx / dir.x();
      double dy = ratio * dir.y();
      result.pointA()->setX( br[0] );
      result.pointA()->setY( result.pointA()->y() + dy );
      reclip = true;
     
    }
  }
  
  if(result.pointA()->y() < br[1]) 
  {
    if( result.pointB()->y() < br[1])
    {
      return result;
    }
    else {
      // clip point A.y
      result.pointA()->setY( br[1] );
      
      double dy = br[1] - result.pointA()->y(); 
      double ratio = dy / dir.y(); 
      double dx = ratio * dir.x();
      
      result.pointA()->setX( result.pointA()->x() + dx );
      result.pointA()->setY( br[1] );
      reclip = true;
    }
  }
  
  result.setType( IN );
   
  /// Clip for pointB

  // clip PointB to tl
  
  if( result.pointB()->x() < tl[0])
    {
      
      
      double dx = tl[0] - result.pointB()->x();
      double ratio = dx / dir.x();
      double dy = ratio * dir.y();
      result.pointB()->setX( tl[0] );
      result.pointB()->setY( result.pointB()->y() + dy );
      reclip = true;
    }
  
    if( result.pointB()->y() > tl[1])
    {
      
      double dy = tl[1] - result.pointB()->y(); 
      double ratio = dy / dir.y(); 
      double dx = ratio * dir.x();
      
      result.pointB()->setX( result.pointB()->x() + dx );
      result.pointB()->setY( tl[1] );
      reclip = true;
    }
  
  
  // clip PointB to br
  
    if( result.pointB()->x() > br[0])
    {
      
      
      double dx = br[0] - result.pointB()->x();
      double ratio = dx / dir.x();
      double dy = ratio * dir.y();
      result.pointB()->setX( br[0] );
      result.pointB()->setY( result.pointB()->y() + dy );
      reclip = true;
    }
  
    if( result.pointB()->y() < br[1])
    {
      
      double dy = br[1] - result.pointB()->y(); 
      double ratio = dy / dir.y(); 
      double dx = ratio * dir.x();
      
      result.pointB()->setX( result.pointB()->x() + dx );
      result.pointB()->setY( br[1] );
      reclip = true;
    }
  
    //cout << "PA x:" << result.pointA()->x() << ", y:" << result.pointA()->y() << endl;
    //cout << "PB x:" << result.pointB()->x() << ", y:" << result.pointB()->y() << endl;
    //cout <<  endl;
    
  } // end while reclip
  
  return result;
}

void Segment::addCell(Cell *cell) {
  if (!cell)
    return;
  _cells.push_back(cell);
}

void Segment::removeCell(Cell *cell) {
  if (!cell)
    return;
  _cells.remove(cell);
}

/*
int main(int argc, char** argv)
{
  Segment test;
  
  test.setType( Segment::OUT );
  counted_ptr<Point> pa(new Point(1.0,-3.0));
  counted_ptr<Point> pb(new Point(0.0,-1.0));
  test.setPointA(pa);
  test.setPointB(pb);
  
  using namespace std;
  
  Segment clipped = test.clippedToRect(Vec2d(0.0,0.0),Vec2d(1.0,1.0));
  
  cout << "Type: " << clipped.type() << endl;
  
  cout << clipped.pointA()->x() << endl;
  cout << clipped.pointA()->y() << endl;
  
  cout << clipped.pointB()->x() << endl;
  cout << clipped.pointB()->y() << endl;
}
*/
