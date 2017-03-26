//
//  Filename         : stroke.h
//  Author           : Emmanuel Turquin
//  Purpose          : Class representing a 2d stroke for the Canvas.
//  Date of creation : 04/26/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  _STROKE_H
# define _STROKE_H

# include <list>
# include "segment.h"
#include "vectypes.h"

class Chain;

class Stroke
{
public:

  typedef std::list<Segment *> SegmentList;

  typedef enum {
    UNKNOWN,
    IN,
    OUT,
    BORDER,
    MIXED
  } Type;

  Stroke(Segment *seg = 0) {
    _chain = 0;
    _pt_a_params = Vec2d(-0.3,1.0);
    _pt_b_params = Vec2d(-0.3,1.0);
    if (!seg) {
      _length = 0;
      _type = UNKNOWN;
      return;
    }
    seg->setStroke(this);
    _pt_a = seg->pointA();
    _pt_b = seg->pointB();
    _segments.push_back(seg);

  }

  Stroke(const Stroke& stroke) : _segments(stroke._segments) {
    _length = stroke._length;
    _pt_a = stroke._pt_a;
    _pt_b = stroke._pt_b;
    _pt_a_params = stroke._pt_a_params;
    _pt_b_params = stroke._pt_b_params;
    _chain = stroke._chain;
    _type = stroke._type;
  }

  ~Stroke() {
    SegmentList::const_iterator it = _segments.begin();
    SegmentList::const_iterator it_end = _segments.end();
    for (; it != it_end; ++it)
      delete *it;
  }

  // Points.
  
  counted_ptr<Point> pointA() const {
    return _pt_a;
  }

  counted_ptr<Point> pointB() const {
    return _pt_b;
  }

  void setPointA(const counted_ptr<Point>& pt_a) {
    _pt_a = pt_a;
  }

  void setPointB(const counted_ptr<Point>& pt_b) {
    _pt_b = pt_b;
  }

  void reverse() {
    SegmentList::iterator it = _segments.begin();
    SegmentList::iterator it_end = _segments.end();
    for (; it != it_end; ++it)
      (*it)->reverse();
    _segments.reverse();
    counted_ptr<Point>	pt = _pt_a;
    _pt_a = _pt_b;
    _pt_b = pt;
  }

  void updateLength() {
    SegmentList::const_iterator it = _segments.begin();
    SegmentList::const_iterator it_end = _segments.end();
    _length = 0;
    for ( ; it != it_end; ++it)
    {
      (*it)->updateLength();
      _length += (*it)->length();
    }
  }

  double length() const {
    return _length;
  }

  void updateType() {
    unsigned in_count = 0;
    unsigned out_count = 0;
    unsigned mixed_count = 0;
    SegmentList::const_iterator it = _segments.begin();
    SegmentList::const_iterator it_end = _segments.end();
    for ( ; it != it_end; ++it) {
      switch ((*it)->type()) {
      case UNKNOWN:
	_type = UNKNOWN;
	return;
      case IN:
      case BORDER:
	++in_count;
	break;
      case OUT:
	++out_count;
	break;
      case MIXED:
      default:
	++mixed_count;
      }
    }
    if (in_count && !out_count && !mixed_count)
      _type = IN;
    else if (out_count && !in_count && !mixed_count)
      _type = OUT;
    else
      _type = MIXED;
  }

  Type type() const {
    return _type;
  }

  void addSegment(Segment *seg);

  bool join(Stroke *stroke);

  Stroke *split(SegmentList::iterator seg_it);
  Stroke *split(Segment *seg);

  Stroke *cleverSplit(SegmentList::iterator seg_it);
  Stroke *cleverSplit(Segment *seg);

  SegmentList *segments() {
    return &_segments;
  }

  const SegmentList *segments() const {
    return &_segments;
  }

  void setChain(Chain *chain) {
    _chain = chain;
  }

  Chain *chain() {
    return _chain;
  }
  // Used for gaussian fold parameters
  void setAParams(Vec2d params);
  void setBParams(Vec2d params);
  Vec2d getAParams();
  Vec2d getBParams();
  
  void recalculateBB(); 
  const Vec3d bbMin() const { return _bbMin; } // warning - these are not auto recalculated
  const Vec3d bbMax() const { return _bbMax; } // warning - these are not auto recalculated

private:

  double		_length;
  counted_ptr<Point>	_pt_a;
  counted_ptr<Point>	_pt_b;
  SegmentList		_segments;
  Chain			*_chain;
  Type			_type;
  Vec3d _bbMin;
  Vec3d _bbMax;
  
  // JDW store additional parameters with the end points, use for gaussian folds radius and amplitude
  Vec2d _pt_a_params;
  Vec2d _pt_b_params;
};

#endif // _STROKE_H
