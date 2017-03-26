//
//  Filename         : point_utils.cpp
//  Author           : Emmanuel Turquin
//  Purpose          : Various utility functions for points.
//  Date of creation : 04/23/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include "config_pretreatment.h"
#include "repository_pretreatment_canvas.h"
#include "repository_pretreatment_result.h"
#include "repository_canvas_result.h"
#include "utils.h"
#include "point_utils.h"

#include "math.h"

namespace point_utils {

//   Point::Type type(double x, double y)
//   {
//     if (!RepositoryPretreatmentCanvas::texture())
//       return Point::UNKNOWN;
//     unsigned half_size = RepositoryPretreatmentCanvas::textureSize() / 2;
//     unsigned u = static_cast<unsigned>(half_size + (x * half_size) / RepositoryPretreatmentCanvas::bbMax()[0]);
//     unsigned v = static_cast<unsigned>(half_size + (y * half_size) / RepositoryPretreatmentCanvas::bbMax()[1]);
//     if (u && u >= RepositoryPretreatmentCanvas::textureSize())
//       u = RepositoryPretreatmentCanvas::textureSize() - 1;
//     if (v && v >= RepositoryPretreatmentCanvas::textureSize())
//       v = RepositoryPretreatmentCanvas::textureSize() - 1;
//     if (RepositoryPretreatmentCanvas::texture()[(v * RepositoryPretreatmentCanvas::textureSize() + u) * 4 + 3])
//       return Point::IN;
//     return Point::OUT;
//   }

  Point::Type type(double x, double y, double dist)
  {
    unsigned char *current_texture = 0;
    unsigned current_texture_size = 0;

    current_texture = RepositoryPretreatmentCanvas::frontTexture();
    current_texture_size = RepositoryPretreatmentCanvas::frontTextureSize();

    if (!current_texture)
      return Point::UNKNOWN;

    unsigned du, dv;
    if (dist < 0) {
      du = current_texture_size / 100;
      dv = du;
    }
    else {
      if (!RepositoryPretreatmentCanvas::distanceField())
	return Point::UNKNOWN;
      double factor = (dist * current_texture_size) / (RepositoryPretreatmentCanvas::distanceField()->maxDist() * 10);
      du = static_cast<unsigned>(factor * RepositoryCanvasResult::borderFactor()[0] * RepositoryPretreatmentResult::modelBBSize()[0] / RepositoryPretreatmentCanvas::distanceField()->size()[0]);
      dv = static_cast<unsigned>(factor * RepositoryCanvasResult::borderFactor()[1] * RepositoryPretreatmentResult::modelBBSize()[1] / RepositoryPretreatmentCanvas::distanceField()->size()[1]);
    }
    unsigned half_size = current_texture_size / 2;
    unsigned u = static_cast<unsigned>(half_size + (x * half_size) / RepositoryPretreatmentCanvas::bbMax()[0]);
    unsigned v = static_cast<unsigned>(half_size + (y * half_size) / RepositoryPretreatmentCanvas::bbMax()[1]);
    if (u && u >= current_texture_size)
      u = current_texture_size - 1;
    if (v && v >= current_texture_size)
      v = current_texture_size - 1;
    unsigned min_u = u >= du ? u - du : u;
    unsigned max_u = u + du < current_texture_size ? u + du : u;
    unsigned min_v = v >= dv ? v - dv : v;
    unsigned max_v = v + dv < current_texture_size ? v + dv : v;
    unsigned in = 0;
    unsigned total = 0;
    for (unsigned i = min_u; i <= max_u; ++i)
      for (unsigned j = min_v; j <= max_v; ++j) {
	if (current_texture[(j * current_texture_size + i) * 4 + 3])
	  ++in;
	++total;
      }
    if (in == total)
      return Point::IN;
    if (!current_texture[(v * current_texture_size + u) * 4 + 3])
      return Point::OUT;
    return Point::BORDER;
  }

  double zEpsilon() {
    DistanceField *df = RepositoryPretreatmentCanvas::distanceField();
    if (!df)
      return 0;
    return df->size()[2] / (2 * df->sizeZ());
  }

  // Internal use only.
  inline void minDistOut(unsigned i, unsigned j, double& dist, double& z)
  {
    // Note: df != 0 already checked.
    DistanceField *df = RepositoryPretreatmentCanvas::distanceField();
    dist = df->maxDist();
    double dz = df->size()[2] / df->sizeZ();
    double dist_tmp;
    for (int k = df->sizeZ() - 1; k >= 0; --k) {
      dist_tmp = df->dist(i, j, k);
      if (dist_tmp <= zEpsilon()) {
	dist = zEpsilon();
	z = k + (dist - dist_tmp) / (df->dist(i, j, k+1) - dist_tmp);
	for (unsigned k = 0; k < df->sizeZ(); ++k) {
	  dist_tmp = df->dist(i, j, k);
	  if (dist_tmp <= zEpsilon()) {
	    z += k - (dist - dist_tmp) / (df->dist(i, j, k-1) - dist_tmp);
	    break;
	  }
	}
	z /= 2;
	break;
      }
      if (dist_tmp < dist) {
	dist = dist_tmp;
	z = k;
      }
    }
    z = df->min()[2] + dz / 2 + z * dz;
  }

  // Internal use only.
  inline void minDistFront(unsigned i, unsigned j, double& dist, double& z)
  {
    // Note: df != 0 already checked.
    DistanceField *df = RepositoryPretreatmentCanvas::distanceField();
    dist = df->maxDist();
    double dz = df->size()[2] / df->sizeZ();
    double dist_tmp;
    for (int k = df->sizeZ() - 1; k >= 0; --k) {
      dist_tmp = df->dist(i, j, k);
      if (dist_tmp <= zEpsilon()) {
	dist = zEpsilon();
	z = k + (dist - dist_tmp) / (df->dist(i, j, k+1) - dist_tmp);
	break;
      }
      if (dist_tmp < dist) {
	dist = dist_tmp;
	z = k;
      }
    }
    z = df->min()[2] + dz / 2 + z * dz;
  }

  // Internal use only.
  inline void minDistBack(unsigned i, unsigned j, double& dist, double& z)
  {
    // Note: df != 0 already checked.
    DistanceField *df = RepositoryPretreatmentCanvas::distanceField();
    dist = df->maxDist();
    double dz = df->size()[2] / df->sizeZ();
    double dist_tmp;
    for (unsigned k = 0; k < df->sizeZ(); ++k) {
      dist_tmp = df->dist(i, j, k);
      if (dist_tmp <= zEpsilon()) {
	dist = zEpsilon();
	z = k - (dist - dist_tmp) / (df->dist(i, j, k-1) - dist_tmp);
	break;
      }
      if (dist_tmp < dist) {
	dist = dist_tmp;
	z = k;
      }
    }
    z = df->min()[2] + dz / 2 + z * dz;
  }

  void minDist(double x, double y, double& dist, double& z, bool front)
  {
    DistanceField *df = RepositoryPretreatmentCanvas::distanceField();
    if (!df)
      return;

    void (*minDist)(unsigned, unsigned, double&, double&);
    if (type(x, y) == Point::OUT)
      minDist = minDistOut;
    else if (front)
      minDist = minDistFront;
    else
      minDist = minDistBack;

    x = df->sizeX() * (x - df->min()[0]) / (df->size()[0]) - 0.5;
    y = df->sizeY() * (y - df->min()[1]) / (df->size()[1]) - 0.5;

    unsigned char region = 0;
    if (x <= 0 || x >= df->sizeX() - 1)
      region += 1;
    if (y <= 0 || y >= df->sizeY() - 1)
      region += 2;

    unsigned i, j;
    i = static_cast<unsigned>(utils::max(0.0, x));
    j = static_cast<unsigned>(utils::max(0.0, y));
    if (i && i >= df->sizeX())
      i = df->sizeX() - 1;
    if (j && j >= df->sizeY())
      j = df->sizeY() - 1;

    if (region == 3) {
      minDist(i, j, dist, z);
    }
    else if (region == 2) {
      double z1, z2;
      double dist1, dist2;
      double alpha = x - i;
      minDist(i, j, dist1, z1);
      minDist(i+1, j, dist2, z2);
      dist = (1 - alpha) * dist1 + alpha * dist2;
      z = (1 - alpha) * z1 + alpha * z2;
    }
    else if (region == 1) {
      double z1, z2;
      double dist1, dist2;
      double beta = y - j;
      minDist(i, j, dist1, z1);
      minDist(i, j+1, dist2, z2);
      dist = (1 - beta) * dist1 + beta * dist2;
      z = (1 - beta) * z1 + beta * z2;
    }
    else { // region == 0
      double z1, z2, z3, z4;;
      double dist1, dist2, dist3, dist4;
      double alpha = x - i;
      double beta = y - j;
      minDist(i, j, dist1, z1);
      minDist(i+1, j, dist2, z2);
      minDist(i, j+1, dist3, z3);
      minDist(i+1, j+1, dist4, z4);
      dist1 = (1 - alpha) * dist1 + alpha * dist2;
      dist2 = (1 - alpha) * dist3 + alpha * dist4;
      dist = (1 - beta) * dist1 + beta * dist2;
      z1 = (1 - alpha) * z1 + alpha * z2;
      z2 = (1 - alpha) * z3 + alpha * z4;
      z = (1 - beta) * z1 + beta * z2;
    }
  }

  // Internal use only.
  inline void zFromDistFront(unsigned i, unsigned j, double dist, double& z)
  {
    // Note: df != 0 already checked.
    DistanceField *df = RepositoryPretreatmentCanvas::distanceField();
    double dz = df->size()[2] / df->sizeZ();
    double smaller_diff = df->maxDist();
    if (dist < zEpsilon())
      dist = zEpsilon();
    double dist_tmp;
    for (int k = df->sizeZ() - 1; k >= 0; --k) {
      dist_tmp = df->dist(i, j, k);
      if (dist_tmp <= dist) {
        z = k + (dist - dist_tmp) / (df->dist(i, j, k+1) - dist_tmp);
        break;
      }
      if (dist_tmp - dist < smaller_diff) {
        smaller_diff = dist_tmp - dist;
        z = k;
      }
    }
    z = df->min()[2] + dz / 2 + z * dz;
  }

  // Internal use only.
  inline void zFromDistBack(unsigned i, unsigned j, double dist, double& z)
  {
    // Note: df != 0 already checked.
    DistanceField *df = RepositoryPretreatmentCanvas::distanceField();
    double dz = df->size()[2] / df->sizeZ();
    double smaller_diff = df->maxDist();
    if (dist < zEpsilon())
      dist = zEpsilon();
    double dist_tmp;
    for (unsigned k = 0; k < df->sizeZ(); ++k) {
      dist_tmp = df->dist(i, j, k);
      if (dist_tmp <= dist) {
	z = k - (dist - dist_tmp) / (df->dist(i, j, k-1) - dist_tmp);
	break;
      }
      if (dist_tmp - dist < smaller_diff) {
	smaller_diff = dist_tmp - dist;
	z = k;
      }
    }
    z = df->min()[2] + dz / 2 + z * dz;
  }

  void zFromDist(double x, double y, double dist, double& z, bool front)
  {
    DistanceField *df = RepositoryPretreatmentCanvas::distanceField();
    if (!df)
      return;

    void (*zFromDist)(unsigned, unsigned, double, double&);
    if (front)
      zFromDist = zFromDistFront;
    else
      zFromDist = zFromDistBack;

    x = df->sizeX() * (x - df->min()[0]) / (df->size()[0]) - 0.5;
    y = df->sizeY() * (y - df->min()[1]) / (df->size()[1]) - 0.5;

    unsigned char region = 0;
    if (x <= 0 || x >= df->sizeX() - 1)
      region += 1;
    if (y <= 0 || y >= df->sizeY() - 1)
      region += 2;

    unsigned i, j;
    i = static_cast<unsigned>(utils::max(0.0, x));
    j = static_cast<unsigned>(utils::max(0.0, y));
    if (i && i >= df->sizeX())
      i = df->sizeX() - 1;
    if (j && j >= df->sizeY())
      j = df->sizeY() - 1;

    if (region == 3) {
      zFromDist(i, j, dist, z);
    }
    else if (region == 2) {
      double z1, z2;
      double alpha = x - i;
      zFromDist(i, j, dist, z1);
      zFromDist(i+1, j, dist, z2);
      z = (1 - alpha) * z1 + alpha * z2;
    }
    else if (region == 1) {
      double z1, z2;
      double beta = y - j;
      zFromDist(i, j, dist, z1);
      zFromDist(i, j+1, dist, z2);
      z = (1 - beta) * z1 + beta * z2;
    }
    else { // region == 0
      double z1, z2, z3, z4;
      double alpha = x - i;
      double beta = y - j;
      zFromDist(i, j, dist, z1);
      zFromDist(i+1, j, dist, z2);
      zFromDist(i, j+1, dist, z3);
      zFromDist(i+1, j+1, dist, z4);
      z1 = (1 - alpha) * z1 + alpha * z2;
      z2 = (1 - alpha) * z3 + alpha * z4;
      z = (1 - beta) * z1 + beta * z2;
    }
     // TODO Remove this
    if(-26.0 > z || z > 27.0)
      std::cout << "zFromDist[" << x << ", " << y << ", " << z << "]: " << dist << std::endl;
    
  }

  void distFromZ(double x, double y, double z, double& dist)
  {
    DistanceField *df = RepositoryPretreatmentCanvas::distanceField();
    if (!df)
      return;

    x = df->sizeX() * (x - df->min()[0]) / (df->size()[0]) - 0.5;
    y = df->sizeY() * (y - df->min()[1]) / (df->size()[1]) - 0.5;
    z = df->sizeZ() * (z - df->min()[2]) / (df->size()[2]) - 0.5;

    dist = df->dist(x, y, z);
    /* //TODO remove this
    if(-26.0 > dist || dist > 27.0)
      std::cout << "Dist[" << x << ", " << y << ", " << z << "]: " << dist << std::endl;
    */
  }
  
  void distNormalisedTo1stOrder(double x, double y, double z, double& dist)
  {
    DistanceField *df = RepositoryPretreatmentCanvas::distanceField();
    if (!df)
      return;

    //x = df->sizeX() * (x - df->min()[0]) / (df->size()[0]) - 0.5;
    //y = df->sizeY() * (y - df->min()[1]) / (df->size()[1]) - 0.5;
    //z = df->sizeZ() * (z - df->min()[2]) / (df->size()[2]) - 0.5;
    
    double d = 0.0;
    Vec3d grad_d;
    distFromZ(x, y, z, d);
    grad_d = df->gradTrilinear(x,y,z);
    grad_d.normalizeSafe();
    double grad_d_sq = grad_d[0]*grad_d[0]+
                       grad_d[1]*grad_d[1]+
                       grad_d[2]*grad_d[2];
    double d_sq = d * d;
          

    dist =  d / sqrt(d_sq + grad_d_sq);
    // scale back into DF range
    //dist *= df->maxDist();
  }

}; // end of namespace point_utils
