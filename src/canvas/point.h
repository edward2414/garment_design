//
//  Filename         : point.h
//  Author           : Emmanuel Turquin
//  Purpose          : Class representing a 2d point for the Canvas.
//  Date of creation : 04/26/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  _POINT_H
# define _POINT_H

#include <iostream>
using namespace std;

class Point
{
  
  friend ostream& operator<<
      (ostream& out, const Point& point);
  friend istream& operator>>
      (istream& in, Point& point);
public:
  
  typedef enum {
    UNKNOWN,
    IN,
    OUT,
    BORDER
  } Type;

  Point(double x = 0, double y = 0, const Type& type = UNKNOWN) {
    _coord[0] = x;
    _coord[1] = y;
    _coord[2] = 0;
    _dist = 0;
    _type = type;
  }

  Point(const Point& point) {
    _coord[0] = point._coord[0];
    _coord[1] = point._coord[1];
    _coord[2] = point._coord[2];
    _dist = point._dist;
    _type = point._type;
  }
  
  ~Point() {}

  // Coordinates.
  
  double x() const {
    return _coord[0];
  }

  double y() const {
    return _coord[1];
  }

  double z() const {
    return _coord[2];
  }

  void setX(double x) {
    _coord[0] = x;
  }

  void setY(double y) {
    _coord[1] = y;
  }

  void setZ(double z) {
    _coord[2] = z;
  }

  const double *address() const {
    return _coord;
  }

  // Return (and set) the type of the point.

  Type type() const {
    return _type;
  }

  void setType(Type type) {
    _type = type;
  }

  // Return (and set) the distance from the surface of the point.

  double distance() const {
    return _dist;
  }

  void setDistance(double dist) {
    _dist = dist;
  }

private:

  double	_coord[3];
  double	_dist;
  Type		_type;
};

#endif // _POINT_H
