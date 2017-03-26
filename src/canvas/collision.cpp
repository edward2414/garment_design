//
//  Filename         : collision.cpp
//  Author           : Emmanuel Turquin
//  Purpose          : Functions to detect collisions between segments.
//  Date of creation : 05/04/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include "utils.h"
#include "stroke.h"
#include "collision.h"

bool Collision::collisionTest(Segment *seg1, Segment *seg2, List& l)
{
  if (!seg1 || !seg2)
    return false;

  if (seg1->pointA().get() == seg2->pointB().get() ||
      seg1->pointB().get() == seg2->pointB().get() ||
      seg1->pointA().get() == seg2->pointA().get() ||
      seg1->pointB().get() == seg2->pointA().get())
    return false;
    
  if (utils::collisionSegmentSegment(*seg1, *seg2)) {
    List::iterator it = l.begin();
    List::iterator it_end = l.end();
    for (; it != it_end; ++it)
      // Collision already registered.
      if (seg1 == it->segment1() && seg2 == it->segment2())
	return false;
    // New collision detected!
    l.push_back(Collision(seg1, seg2));
    return true;
  }

  return false;
}

bool Collision::self() const
{
  return (_seg1 && _seg2) && (_seg1->stroke() == _seg2->stroke());
}
