//
//  Filename         : chain.h
//  Author           : Emmanuel Turquin
//  Purpose          : Class representing a chain of strokes for the Canvas.
//  Date of creation : 04/26/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  _CHAIN_H
# define _CHAIN_H

# include <list>
# include "stroke.h"
#include "vectypes.h"

class Chain
{
public:

  typedef std::list<Stroke *> StrokeList;

  Chain(Stroke *stroke = 0) {
    if (!stroke)
      return;

    stroke->setChain(this);
    _pt_a = stroke->pointA();
    _pt_b = stroke->pointB();
    _strokes.push_back(stroke);
  }

  Chain(const Chain& chain) : _strokes(chain._strokes) {
    _pt_a = chain._pt_a;
    _pt_b = chain._pt_b;
  }

  ~Chain() {
    StrokeList::const_iterator it = _strokes.begin();
    StrokeList::const_iterator it_end = _strokes.end();
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

  bool cycle() const {
    return _pt_a.get() == _pt_b.get();
  }

  void reverse() {
    StrokeList::iterator it = _strokes.begin();
    StrokeList::iterator it_end = _strokes.end();
    for ( ; it != it_end; ++it)
      (*it)->reverse();
    _strokes.reverse();
    counted_ptr<Point>	pt = _pt_a;
    _pt_a = _pt_b;
    _pt_b = pt;
  }

  void addStroke(Stroke *stroke);

  bool join(Chain *chain);

  Chain *split(StrokeList::iterator str_it);
  Chain *split(Stroke *stroke);

  void removeStroke(StrokeList::iterator str_it);
  void removeStroke(Stroke *stroke);

  StrokeList *strokes() {
    return &_strokes;
  }

  const StrokeList *strokes() const {
    return &_strokes;
  }
  
  const Vec3d bbMin() const {
    return _bbMin; // warning - these are not auto recalculated
  }
  
  const Vec3d bbMax() const {
    return _bbMax; // warning - these are not auto recalculated
  }
  
  void recalculateBB();
  

private:

  counted_ptr<Point>	_pt_a;
  counted_ptr<Point>	_pt_b;
  StrokeList	_strokes;
  Vec3d _bbMin;
  Vec3d _bbMax;
};

#endif // _CHAIN_H
