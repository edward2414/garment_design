//
//  Filename         : utils.h
//  Author           : Emmanuel Turquin
//  Purpose          : Misc. 2D utility functions.
//  Date of creation : 04/23/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  UTILS_H
# define UTILS_H

# include "vectypes.h"
#include <string>

class Segment;
class Triangle;
class LCModel;
class Garment;
class DockPattern;

namespace utils {

  // General.
  
  double gauss2d(double amp,double sigma,double x,double y);
  double distPointLine(const Vec3d p1, const Vec3d p2, const Vec3d r);
  Vec3d vecFromPointToLine(const Vec3d p1, const Vec3d p2, const Vec3d r);

  double normalizedRand();

  template <class T>
    T max(const T& a, const T& b)
    {
      return a >= b ? a : b;
    }

  template <class T>
    T min(const T& a, const T& b)
    {
      return a <= b ? a : b;
    }

  template <class T>
    T abs(const T& val)
    {
      return val < 0 ? -val : val;
    }

  template <class T>
    void swap(T& a, T& b)
    {
      T tmp(a);
      a = b;
      b = tmp;
    }

  // 2D.

  unsigned segmentsInSphere(double a, double b, double x, double y, double max_length);

  bool isInSphere(double a, double b, double x, double y, double max_length);

  bool isInWindow(double a, double b, double x, double y, double max_dx, double max_dy);

  double dotProduct(double x1, double y1, double x2, double y2);

  double normalizedDotProduct(double x1, double y1, double x2, double y2);

  double direction(double x1, double y1, double x2, double y2, double x3, double y3);

  bool collisionSegmentSegment(const Segment& seg1, const Segment& seg2);

  bool overlapSegmentBox(const Vec2d& boxcenter,
			 const Vec2d& boxhalfsize,
			 const Segment& seg);

  bool collisionSegmentX(const Segment& seg, double x, double& y);
  bool collisionSegmentY(const Segment& seg, double y, double& x);

  // 3D.

  // Adapted from Tomas Akenine-M�ler code.
  bool overlapTriangleBox(const Vec3d& boxcenter,
			  const Vec3d& boxhalfsize,
			  const Triangle& tr);
  bool invertMatrix(double m[4][4]);

  double squareDistPointBox(const Vec3d& point, const Vec3d& min, const Vec3d& max);
  double distPointBox(const Vec3d& point, const Vec3d& min, const Vec3d& max);

  int factorial(int n);
  // Adapted from David Eberly code.
  double squareDistPointTriangle(const Vec3d& point, const Triangle& tr);
  double distPointTriangle(const Vec3d& point, const Triangle& tr);

  void bezier3(double& z,
	       const double z1,
	       const double z2,
	       const double z3,
	       const double z4,
	       double t);
  bool importOBJ(LCModel* model, const char* filename);
  bool importOBJAtlas(LCModel* model, const char* filename);
  bool importOBJGarment(Garment* pGarment, const char* filename, DockPattern* dockpattern);
} // end of namespace utils

#endif // UTILS_H
