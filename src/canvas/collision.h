//
//  Filename         : collision.h
//  Author           : Emmanuel Turquin
//  Purpose          : Functions to detect collisions between segments.
//  Date of creation : 05/04/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  _COLLISION_H
# define _COLLISION_H

# include <list>

class Segment;

class Collision
{
public:

  typedef std::list<Collision>	List;

  Collision(Segment *seg1 = 0, Segment *seg2 = 0)
  {
    _seg1 = seg1;
    _seg2 = seg2;
  }

  Collision(const Collision& col)
  {
    _seg1 = col._seg1;
    _seg2 = col._seg2;
  }

  ~Collision()
  {

  }

  Segment *segment1() {
    return _seg1;
  }

  Segment *segment2() {
    return _seg2;
  }

  bool self() const;

  static bool collisionTest(Segment *seg1, Segment *seg2, List& l);

private:

  Segment *_seg1;
  Segment *_seg2;
};



#endif // _COLLISION_H
