//
//  Filename         : triangle.h
//  Author(s)        : Emmanuel Turquin
//  Purpose          : Triangle (in 3D).
//  Date of creation : 05/2'/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  TRIANGLE_H
# define TRIANGLE_H

# include "vectypes.h"

class Triangle
{
public:

  Triangle(const Vec3d *a = 0, const Vec3d *b = 0, const Vec3d *c = 0) {
    _a = a;
    _b = b;
    _c = c;
  }

  Triangle(const Triangle& t) {
    _a = t._a;
    _b = t._b;
    _c = t._c;
  }

  ~Triangle() {}

  const Vec3d *pointA() const {
    return _a;
  }

  const Vec3d *pointB() const {
    return _b;
  }

  const Vec3d *pointC() const {
    return _c;
  }

  void setPointA(const Vec3d *a) {
     _a = a;
  }

  void setPointB(const Vec3d *b) {
    _b = b;
  }

  void setPointC(const Vec3d *c) {
    _c = c;
  }

  Vec3d normal() const {
    return (*_b - *_a) ^ (*_c - *_a);
  }

private:

  const Vec3d *_a;
  const Vec3d *_b;
  const Vec3d *_c;
};

#endif // TRIANGLE_H
