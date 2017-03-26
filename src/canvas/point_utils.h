//
//  Filename         : point_utils.h
//  Author           : Emmanuel Turquin
//  Purpose          : Various utility functions for points.
//  Date of creation : 04/23/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef POINT_UTILS_H
# define POINT_UTILS_H

# include "repository_pretreatment_canvas.h"
# include "point.h"

namespace point_utils {

  double zEpsilon();
  // Look at the body texture map to determine whether point is IN, OUT or BORDER
  Point::Type type(double x, double y, double dist = -1);
  void minDist(double x, double y, double& dist, double& z, bool front);
  void zFromDist(double x, double y, double dist, double& z, bool front);
  void distFromZ(double x, double y, double z, double& dist);
  void distNormalisedTo1stOrder(double x, double y, double z, double& dist);

}; // end of namespace point_utils

#endif
