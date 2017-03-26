//
//  Filename         : segment.h
//  Author           : Emmanuel Turquin
//  Purpose          : Class representing a 2d segment for the Canvas.
//  Date of creation : 04/26/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  _SEGMENT_H
# define _SEGMENT_H

# include <list>
# include <cmath>
# include "counted_ptr.h"
# include "point.h"
#include "vectypes.h"

class Stroke;
class Cell;


class Segment
{
public:

  typedef std::list<Cell*> CellList;

  typedef enum {
    UNKNOWN,
    IN,
    OUT,
    BORDER,
    MIXED
  } Type;

  Segment()
  {
    _stroke = 0;
    _length = 0;
    _type = UNKNOWN;
  }

  Segment(const counted_ptr<Point>& pt_a, const counted_ptr<Point>& pt_b) : _pt_a(pt_a), _pt_b(pt_b)
  {
    _stroke = 0;
    _length = 0;
    _type = UNKNOWN;
  }

  Segment(const Segment& seg) : _cells(seg._cells) {
    _length = seg._length;
    _pt_a = seg._pt_a;
    _pt_b = seg._pt_b;
    _stroke = seg._stroke;
    _type = seg._type;
  }

  ~Segment();

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
    counted_ptr<Point>	pt = _pt_a;
    _pt_a = _pt_b;
    _pt_b = pt;
  }

  void updateLength() {
    if (!_pt_a.get() || !_pt_b.get())
      return;
    double dx = _pt_a->x() - _pt_b->x();
    double dy = _pt_a->y() - _pt_b->y();
    _length = sqrt(dx * dx + dy * dy);
  }

  double length() const {
    return _length;
  }

  void updateType() {
    if (!_pt_a.get() || !_pt_b.get())
      return;
    Point::Type type_a = _pt_a->type();
    Point::Type type_b = _pt_b->type();
    if (type_a == Point::UNKNOWN || type_b == Point::UNKNOWN)
      _type = UNKNOWN;
    else if (type_a == Point::IN && type_b != Point::OUT ||
	     type_a != Point::OUT && type_b == Point::IN)
      _type = IN;
    else if (type_a == Point::OUT && type_b == Point::OUT)
      _type = OUT;
    else if (type_a == Point::BORDER && type_b == Point::BORDER)
      _type = BORDER;
    else
      _type = MIXED;
  }

  void setType(Type type) {
    _type = type;
  }

  Type type() const {
    return _type;
  }

  void addCell(Cell *cell);

  void removeCell(Cell *cell);

  void setStroke(Stroke *stroke) {
    _stroke = stroke;
  }

  Stroke *stroke() {
    return _stroke;
  }
  
  Segment clippedToRect(const Vec2d topLeft, const Vec2d bottomRight);

private:

  double		_length;
  counted_ptr<Point>	_pt_a;
  counted_ptr<Point>	_pt_b;
  Stroke		*_stroke;
  CellList		_cells;
  Type			_type;
};

#endif // _SEGMENT_H
