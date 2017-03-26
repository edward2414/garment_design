//
//  Filename         : utils.cpp
//  Author           : Emmanuel Turquin
//  Purpose          : Misc. 2D utility functions.
//  Date of creation : 04/23/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include <stdlib.h>
#include <math.h>
#include <string>
#include <map>
#include "GarmentSection.h"
#include "GarmentAtlas.h"
#include "Garment.h"
#include "segment.h"
#include "triangle.h"
#include "indexed_triangle.h"
#include "utils.h"
#include "lcmodel.h"
#include <cassert>

// Only required to test adding items to ListView in dock_pattern_widget
#include "dock_pattern.h"
#include "dock_pattern_widget.h"
#include <qlistview.h>

#include<iostream>
#include<boost/tokenizer.hpp>

using namespace std;

// defined in main.cpp in global namespace
//extern DockPattern *dock_pattern;
//extern DockPattern *dock_pattern2;

namespace utils {

  // General.
  
  double gauss2d(double amp, double sigma, double x, double y)
  {
    double exponent = -1.0f*((x * x) + (y * y))/(2.0*sigma*sigma) ;
    return (amp * pow(M_E, exponent));
  }

  double normalizedRand()
  {
    return (double)rand() / (double)(RAND_MAX);
  }

  // 2D.

  unsigned segmentsInSphere(double a, double b, double x, double y, double length)
  {
    if (!length)
      return 0;
    double dx = a - x;
    double dy = b - y;
    double radius = sqrt(dx * dx + dy * dy);
    return static_cast<unsigned>(radius / length);
  }

  bool isInSphere(double a, double b, double x, double y, double length)
  {
    double dx = a - x;
    double dy = b - y;
    double radius = sqrt(dx * dx + dy * dy);
    return radius <= length;
  }

  bool isInWindow(double a, double b, double x, double y, double max_dx, double max_dy)
  {
    double dx = a - x;
    double dy = b - y;
    if (dx < 0)
      dx = -dx;
    if (dy < 0)
      dy = -dy;
    return (dx <= max_dx && dy <= max_dy);
  }

  double dotProduct(double x1, double y1, double x2, double y2)
  {
    return x1 * x2 + y1 * y2;
  }

  double normalizedDotProduct(double x1, double y1, double x2, double y2)
  {
    double len1 = sqrt(x1 * x1 + y1 * y1);
    x1 /= len1;
    y1 /= len1;
    double len2 = sqrt(x2 * x2 + y2 * y2);
    x2 /= len2;
    y2 /= len2;
    return dotProduct(x1, y1, x2, y2);
  }

  double direction(double x1, double y1, double x2, double y2, double x3, double y3)
  {
    double u_x = x2 - x1;
    double u_y = y2 - y1;
    double v_x = x3 - x1;
    double v_y = y3 - y1;
    return (u_x * v_y - u_y * v_x);
  }

  //Internal use only.
  inline bool collisionSegmentSegment(double x1a, double y1a, double x1b, double y1b,
				      double x2a, double y2a, double x2b, double y2b)
  {
    // A simple test to discard most of the segments.
    if (utils::max(x1a, x1b) < utils::min(x2a, x2b) ||
	utils::max(x2a, x2b) < utils::min(x1a, x1b) ||
	utils::max(y1a, y1b) < utils::min(y2a, y2b) ||
	utils::max(y2a, y2b) < utils::min(y1a, y1b))
      return false;

    double d1 = utils::direction(x1a, y1a, x1b, y1b, x2a, y2a);
    double d2 = utils::direction(x1a, y1a, x1b, y1b, x2b, y2b);
    double d3 = utils::direction(x2a, y2a, x2b, y2b, x1a, y1a);
    double d4 = utils::direction(x2a, y2a, x2b, y2b, x1b, y1b);
    if (d1 != 0 && d2 != 0 && d1 * d2 > 0)
      return false;
    if (d3 != 0 && d4 != 0 && d3 * d4 > 0)
      return false;

    return true;
  }

  bool collisionSegmentSegment(const Segment& seg1, const Segment& seg2)
  {
    return collisionSegmentSegment(seg1.pointA()->x(),
				   seg1.pointA()->y(),
				   seg1.pointB()->x(),
				   seg1.pointB()->y(),
				   seg2.pointA()->x(),
				   seg2.pointA()->y(),
				   seg2.pointB()->x(),
				   seg2.pointB()->y());
  }

  bool overlapSegmentBox(const Vec2d& boxcenter,
			 const Vec2d& boxhalfsize,
			 const Segment& seg)
  {
    if (seg.pointA()->x() >= boxcenter[0] - boxhalfsize[0] &&
	seg.pointA()->x() <= boxcenter[0] + boxhalfsize[0] &&
	seg.pointA()->y() >= boxcenter[1] - boxhalfsize[1] &&
	seg.pointA()->y() <= boxcenter[1] + boxhalfsize[1])
      return true;
    if (seg.pointB()->x() >= boxcenter[0] - boxhalfsize[0] &&
	seg.pointB()->x() <= boxcenter[0] + boxhalfsize[0] &&
	seg.pointB()->y() >= boxcenter[1] - boxhalfsize[1] &&
	seg.pointB()->y() <= boxcenter[1] + boxhalfsize[1])
      return true;
    if (collisionSegmentSegment(seg.pointA()->x(),
				seg.pointA()->y(),
				seg.pointB()->x(),
				seg.pointB()->y(),
				boxcenter[0] - boxhalfsize[0],
				boxcenter[1] - boxhalfsize[1],
				boxcenter[0] - boxhalfsize[0],
				boxcenter[1] + boxhalfsize[1]))
      return true;
    if (collisionSegmentSegment(seg.pointA()->x(),
				seg.pointA()->y(),
				seg.pointB()->x(),
				seg.pointB()->y(),
				boxcenter[0] - boxhalfsize[0],
				boxcenter[1] + boxhalfsize[1],
				boxcenter[0] + boxhalfsize[0],
				boxcenter[1] + boxhalfsize[1]))
      return true;
    if (collisionSegmentSegment(seg.pointA()->x(),
				seg.pointA()->y(),
				seg.pointB()->x(),
				seg.pointB()->y(),
				boxcenter[0] + boxhalfsize[0],
				boxcenter[1] + boxhalfsize[1],
				boxcenter[0] + boxhalfsize[0],
				boxcenter[1] - boxhalfsize[1]))
      return true;
    if (collisionSegmentSegment(seg.pointA()->x(),
				seg.pointA()->y(),
				seg.pointB()->x(),
				seg.pointB()->y(),
				boxcenter[0] + boxhalfsize[0],
				boxcenter[1] - boxhalfsize[1],
				boxcenter[0] - boxhalfsize[0],
				boxcenter[1] - boxhalfsize[1]))
      return true;
    return false;
  }

  bool collisionSegmentX(const Segment& seg, double x, double& y)
  {
    double x1 = seg.pointA()->x();
    double y1 = seg.pointA()->y();
    double x2 = seg.pointB()->x();
    double y2 = seg.pointB()->y();
    if (x2 < x1) {
      utils::swap(x1, x2);
      utils::swap(y1, y2);
    }
    if (x < x1 || x > x2)
      return false;
    y = y1 + (x - x1) * (y2 - y1) / (x2 - x1);
    return true;
  }

  bool collisionSegmentY(const Segment& seg, double y, double& x)
  {
    double x1 = seg.pointA()->x();
    double y1 = seg.pointA()->y();
    double x2 = seg.pointB()->x();
    double y2 = seg.pointB()->y();
    if (y2 < y1) {
      utils::swap(x1, x2);
      utils::swap(y1, y2);
    }
    if (y < y1 || y > y2)
      return false;
    x = x1 + (y - y1) * (x2 - x1) / (y2 - y1);
    return true;
  }

  // 3D.

  // AABB-triangle overlap test code
  // by Tomas Akenine-M�ler
  // Function: int triBoxOverlap(real boxcenter[3],
  //          real boxhalfsize[3],real triverts[3][3]);
  // History:
  //   2001-03-05: released the code in its first version
  //   2001-06-18: changed the order of the tests, faster
  //
  // Acknowledgement: Many thanks to Pierre Terdiman for
  // suggestions and discussions on how to optimize code.
  // Thanks to David Hunt for finding a ">="-bug!

#define X 0
#define Y 1
#define Z 2

#define FINDMINMAX(x0, x1, x2, min, max) \
  min = max = x0;    \
  if(x1<min) min=x1; \
  if(x1>max) max=x1; \
  if(x2<min) min=x2; \
  if(x2>max) max=x2;

  //======================== X-tests ========================//
#define AXISTEST_X01(a, b, fa, fb)                         \
        p0 = a*v0[Y] - b*v0[Z];                            \
        p2 = a*v2[Y] - b*v2[Z];                            \
        if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
        rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];   \
        if(min>rad || max<-rad) return 0;

#define AXISTEST_X2(a, b, fa, fb)                          \
        p0 = a*v0[Y] - b*v0[Z];                            \
        p1 = a*v1[Y] - b*v1[Z];                            \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
        rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];   \
        if(min>rad || max<-rad) return 0;

  //======================== Y-tests ========================//
#define AXISTEST_Y02(a, b, fa, fb)                         \
        p0 = -a*v0[X] + b*v0[Z];                           \
        p2 = -a*v2[X] + b*v2[Z];                           \
        if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
        rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];   \
        if(min>rad || max<-rad) return 0;

#define AXISTEST_Y1(a, b, fa, fb)                          \
        p0 = -a*v0[X] + b*v0[Z];                           \
        p1 = -a*v1[X] + b*v1[Z];                           \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
        rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];   \
        if(min>rad || max<-rad) return 0;

  //======================== Z-tests ========================//
#define AXISTEST_Z12(a, b, fa, fb)                         \
        p1 = a*v1[X] - b*v1[Y];                            \
        p2 = a*v2[X] - b*v2[Y];                            \
        if(p2<p1) {min=p2; max=p1;} else {min=p1; max=p2;} \
        rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];   \
        if(min>rad || max<-rad) return 0;

#define AXISTEST_Z0(a, b, fa, fb)                          \
        p0 = a*v0[X] - b*v0[Y];                            \
        p1 = a*v1[X] - b*v1[Y];                            \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
        rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];   \
        if(min>rad || max<-rad) return 0;
 
  // Internal procedure.
  inline bool overlapPlaneBox(const Vec3d& normal, double d, const Vec3d& maxbox) {
    Vec3d vmin, vmax;

    for(unsigned q = X; q <= Z; q++) {
      if(normal[q] > 0.0f) {
	vmin[q] = -maxbox[q];
	vmax[q] = maxbox[q];
      }
      else {
	vmin[q] = maxbox[q];
	vmax[q] = -maxbox[q];
      }
    }
    if((normal * vmin) + d > 0.0f)
      return false;
    if((normal * vmax) + d >= 0.0f)
      return true;
    return false;
  }

  // Use separating axis theorem to test overlap between triangle and box
  // need to test for overlap in these directions:
  // 1) the {x,y,z}-directions (actually, since we use the AABB of the triangle
  //    we do not even need to test these)
  // 2) normal of the triangle
  // 3) crossproduct(edge from tri, {x,y,z}-directin)
  //    this gives 3x3=9 more tests
  bool overlapTriangleBox(const Vec3d& boxcenter,
			  const Vec3d& boxhalfsize,
			  const Triangle& tr)
  {
    double min, max, d, p0, p1, p2, rad, fex, fey, fez;  

    // This is the fastest branch on Sun
    // move everything so that the boxcenter is in (0, 0, 0)
    Vec3d v0(*(tr.pointA()) - boxcenter);
    Vec3d v1(*(tr.pointB()) - boxcenter);
    Vec3d v2(*(tr.pointC()) - boxcenter);

    // compute triangle edges
    Vec3d e0(v1 - v0);
    Vec3d e1(v2 - v1);
    Vec3d e2(v0 - v2);

    // Bullet 3:
    // Do the 9 tests first (this was faster)
    fex = fabs(e0[X]);
    fey = fabs(e0[Y]);
    fez = fabs(e0[Z]);
    AXISTEST_X01(e0[Z], e0[Y], fez, fey);
    AXISTEST_Y02(e0[Z], e0[X], fez, fex);
    AXISTEST_Z12(e0[Y], e0[X], fey, fex);

    fex = fabs(e1[X]);
    fey = fabs(e1[Y]);
    fez = fabs(e1[Z]);
    AXISTEST_X01(e1[Z], e1[Y], fez, fey);
    AXISTEST_Y02(e1[Z], e1[X], fez, fex);
    AXISTEST_Z0(e1[Y], e1[X], fey, fex);
  
    fex = fabs(e2[X]);
    fey = fabs(e2[Y]);
    fez = fabs(e2[Z]);
    AXISTEST_X2(e2[Z], e2[Y], fez, fey);
    AXISTEST_Y1(e2[Z], e2[X], fez, fex);
    AXISTEST_Z12(e2[Y], e2[X], fey, fex);

    // Bullet 1:
    // first test overlap in the {x,y,z}-directions
    // find min, max of the triangle each direction, and test for overlap in
    // that direction -- this is equivalent to testing a minimal AABB around
    // the triangle against the AABB

    // test in X-direction
    FINDMINMAX(v0[X], v1[X], v2[X], min, max);
    if (min > boxhalfsize[X] || max < -boxhalfsize[X])
      return false;

    // test in Y-direction
    FINDMINMAX(v0[Y], v1[Y], v2[Y], min, max);
    if (min > boxhalfsize[Y] || max < -boxhalfsize[Y])
      return false;

    // test in Z-direction
    FINDMINMAX(v0[Z], v1[Z], v2[Z], min, max);
    if(min > boxhalfsize[Z] || max < -boxhalfsize[Z])
      return false;

    // Bullet 2:
    // test if the box intersects the plane of the triangle
    // compute plane equation of triangle: normal * x + d = 0
    Vec3d normal(e0 ^ e1);
    d = -(normal * v0); // plane eq: normal.x + d = 0
    if (!overlapPlaneBox(normal, d, boxhalfsize))
      return false;
  
    return true; // box and triangle overlaps
  }

#define EPSILON 0.00000001

  // Matrix inversion (in place).
  bool invertMatrix(double m[4][4])
  {                          
    int i,j,k;               
    int pvt_i[4], pvt_j[4];       /* Locations of pivot elements */
    double pvt_val;               /* Value of current pivot element */
    double hold;                  /* Temporary storage */
    double determinat;            

    determinat = 1.0f;
    for (k=0; k<4; k++)  {
      /* Locate k'th pivot element */
      pvt_val=m[k][k];            /* Initialize for search */
      pvt_i[k]=k;
      pvt_j[k]=k;
      for (i=k; i<4; i++) {
	for (j=k; j<4; j++) {
	  if (fabs(m[i][j]) > fabs(pvt_val)) {
	    pvt_i[k]=i;
	    pvt_j[k]=j;
	    pvt_val=m[i][j];
	  }
	}
      }

      /* Product of pivots, gives determinant when finished */
      determinat*=pvt_val;
      if (fabs(determinat) < EPSILON) {    
	return false;           /* Matrix is singular (zero determinant) */
      }

      /* "Interchange" rows (with sign change stuff) */
      i=pvt_i[k];
      if (i!=k) {               /* If rows are different */
	for (j=0; j<4; j++) {
	  hold=-m[k][j];
	  m[k][j]=m[i][j];
	  m[i][j]=hold;
	}
      }

      /* "Interchange" columns */
      j=pvt_j[k];
      if (j!=k) {              /* If columns are different */
	for (i=0; i<4; i++) {
	  hold=-m[i][k];
	  m[i][k]=m[i][j];
	  m[i][j]=hold;
	}
      }
    
      /* Divide column by minus pivot value */
      for (i=0; i<4; i++) {
	if (i!=k) m[i][k]/=( -pvt_val) ; 
      }

      /* Reduce the matrix */
      for (i=0; i<4; i++) {
	hold = m[i][k];
	for (j=0; j<4; j++) {
	  if (i!=k && j!=k) m[i][j]+=hold*m[k][j];
	}
      }

      /* Divide row by pivot */
      for (j=0; j<4; j++) {
	if (j!=k) m[k][j]/=pvt_val;
      }

      /* Replace pivot by reciprocal (at last we can touch it). */
      m[k][k] = 1.0f/pvt_val;
    }

    /* That was most of the work, one final pass of row/column interchange */
    /* to finish */
    for (k=4-2; k>=0; k--) { /* Don't need to work with 1 by 1 corner*/
      i=pvt_j[k];            /* Rows to swap correspond to pivot COLUMN */
      if (i!=k) {            /* If rows are different */
	for(j=0; j<4; j++) {
	  hold = m[k][j];
	  m[k][j]=-m[i][j];
	  m[i][j]=hold;
	}
      }

      j=pvt_i[k];           /* Columns to swap correspond to pivot ROW */
      if (j!=k)             /* If columns are different */
	for (i=0; i<4; i++) {
	  hold=m[i][k];
	  m[i][k]=-m[i][j];
	  m[i][j]=hold;
	}
    }
    return true;                          
  }

  double distPointLine(const Vec3d p1, const Vec3d p2, const Vec3d r)
  {
    // line from p1 to p2, find distance of r from line.
    Vec3d w = r - p1; // vector from first point in line to required point r
    Vec3d v = p2 -p1; // the line

    return (v^w).norm() / v.norm();
  }
  
  Vec3d vecFromPointToLine(const Vec3d p1, const Vec3d p2, const Vec3d r)
  {
    // line from p1 to p2, find Vector from r to the line.
    Vec3d w = r - p1; // vector from first point in line to required point r
    Vec3d v = p2 -p1; // the line
    
    // want point q, which is p + tv, find t

    double t = (v * w) / v.squareNorm();
    Vec3d q = p1 + t * v; // point q on the original line;
    return (q-r);
  }
  
  double squareDistPointBox(const Vec3d& point, const Vec3d& min, const Vec3d& max)
  {
    double dist = 0;
    double delta;
    if (point[0] < min[0]) {
      delta = point[0] - min[0];
      dist += delta * delta;
    }
    else if (point[0] > max[0]) {
      delta = point[0] - max[0];
      dist += delta * delta;
    }
    if (point[1] < min[1]) {
      delta = point[1] - min[1];
      dist += delta * delta;
    }
    else if (point[1] > max[1]) {
      delta = point[1] - max[1];
      dist += delta * delta;
    }
    if (point[2] < min[2]) {
      delta = point[2] - min[2];
      dist += delta * delta;
    }
    else if (point[2] > max[2]) {
      delta = point[2] - max[2];
      dist += delta * delta;
    }
    return dist;
  }

  double distPointBox(const Vec3d& point, const Vec3d& min, const Vec3d& max) {
    return sqrt(squareDistPointBox(point, min, max));
  }

  double squareDistPointTriangle(const Vec3d& point, const Triangle& tr)
  {
    Vec3d kDiff = *tr.pointA() - point;
    Vec3d edge0 = *tr.pointB() - *tr.pointA();
    Vec3d edge1 = *tr.pointC() - *tr.pointA();
    double fA00 = edge0.squareNorm();
    double fA01 = edge0 * edge1;
    double fA11 = edge1.squareNorm();
    double fB0 = kDiff * edge0;
    double fB1 = kDiff * edge1;
    double fC = kDiff.squareNorm();
    double fDet = utils::abs(fA00*fA11-fA01*fA01);
    double fS = fA01*fB1-fA11*fB0;
    double fT = fA01*fB0-fA00*fB1;
    double fSqrDist;

    if ( fS + fT <= fDet )
      {
        if ( fS < 0.0 )
	  {
            if ( fT < 0.0 )  // region 4
	      {
                if ( fB0 < 0.0 )
		  {
                    fT = 0.0;
                    if ( -fB0 >= fA00 )
		      {
                        fS = 1.0;
                        fSqrDist = fA00+(2.0)*fB0+fC;
		      }
                    else
		      {
                        fS = -fB0/fA00;
                        fSqrDist = fB0*fS+fC;
		      }
		  }
                else
		  {
                    fS = 0.0;
                    if ( fB1 >= 0.0 )
		      {
                        fT = 0.0;
                        fSqrDist = fC;
		      }
                    else if ( -fB1 >= fA11 )
		      {
                        fT = 1.0;
                        fSqrDist = fA11+(2.0)*fB1+fC;
		      }
                    else
		      {
                        fT = -fB1/fA11;
                        fSqrDist = fB1*fT+fC;
		      }
		  }
	      }
            else  // region 3
	      {
                fS = 0.0;
                if ( fB1 >= 0.0 )
		  {
                    fT = 0.0;
                    fSqrDist = fC;
		  }
                else if ( -fB1 >= fA11 )
		  {
                    fT = 1.0;
                    fSqrDist = fA11+(2.0)*fB1+fC;
		  }
                else
		  {
                    fT = -fB1/fA11;
                    fSqrDist = fB1*fT+fC;
		  }
	      }
	  }
        else if ( fT < 0.0 )  // region 5
	  {
            fT = 0.0;
            if ( fB0 >= 0.0 )
	      {
                fS = 0.0;
                fSqrDist = fC;
	      }
            else if ( -fB0 >= fA00 )
	      {
                fS = 1.0;
                fSqrDist = fA00+(2.0)*fB0+fC;
	      }
            else
	      {
                fS = -fB0/fA00;
                fSqrDist = fB0*fS+fC;
	      }
	  }
        else  // region 0
	  {
            // minimum at interior point
            double fInvDet = (1.0)/fDet;
            fS *= fInvDet;
            fT *= fInvDet;
            fSqrDist = fS*(fA00*fS+fA01*fT+(2.0)*fB0) +
	      fT*(fA01*fS+fA11*fT+(2.0)*fB1)+fC;
	  }
      }
    else
      {
        double fTmp0, fTmp1, fNumer, fDenom;

        if ( fS < 0.0 )  // region 2
	  {
            fTmp0 = fA01 + fB0;
            fTmp1 = fA11 + fB1;
            if ( fTmp1 > fTmp0 )
	      {
                fNumer = fTmp1 - fTmp0;
                fDenom = fA00-2.0f*fA01+fA11;
                if ( fNumer >= fDenom )
		  {
                    fS = 1.0;
                    fT = 0.0;
                    fSqrDist = fA00+(2.0)*fB0+fC;
		  }
                else
		  {
                    fS = fNumer/fDenom;
                    fT = 1.0 - fS;
                    fSqrDist = fS*(fA00*fS+fA01*fT+2.0f*fB0) +
		      fT*(fA01*fS+fA11*fT+(2.0)*fB1)+fC;
		  }
	      }
            else
	      {
                fS = 0.0;
                if ( fTmp1 <= 0.0 )
		  {
                    fT = 1.0;
                    fSqrDist = fA11+(2.0)*fB1+fC;
		  }
                else if ( fB1 >= 0.0 )
		  {
                    fT = 0.0;
                    fSqrDist = fC;
		  }
                else
		  {
                    fT = -fB1/fA11;
                    fSqrDist = fB1*fT+fC;
		  }
	      }
	  }
        else if ( fT < 0.0 )  // region 6
	  {
            fTmp0 = fA01 + fB1;
            fTmp1 = fA00 + fB0;
            if ( fTmp1 > fTmp0 )
	      {
                fNumer = fTmp1 - fTmp0;
                fDenom = fA00-(2.0)*fA01+fA11;
                if ( fNumer >= fDenom )
		  {
                    fT = 1.0;
                    fS = 0.0;
                    fSqrDist = fA11+(2.0)*fB1+fC;
		  }
                else
		  {
                    fT = fNumer/fDenom;
                    fS = 1.0 - fT;
                    fSqrDist = fS*(fA00*fS+fA01*fT+(2.0)*fB0) +
		      fT*(fA01*fS+fA11*fT+(2.0)*fB1)+fC;
		  }
	      }
            else
	      {
                fT = 0.0;
                if ( fTmp1 <= 0.0 )
		  {
                    fS = 1.0;
                    fSqrDist = fA00+(2.0)*fB0+fC;
		  }
                else if ( fB0 >= 0.0 )
		  {
                    fS = 0.0;
                    fSqrDist = fC;
		  }
                else
		  {
                    fS = -fB0/fA00;
                    fSqrDist = fB0*fS+fC;
		  }
	      }
	  }
        else  // region 1
	  {
            fNumer = fA11 + fB1 - fA01 - fB0;
            if ( fNumer <= 0.0 )
	      {
                fS = 0.0;
                fT = 1.0;
                fSqrDist = fA11+(2.0)*fB1+fC;
	      }
            else
	      {
                fDenom = fA00-2.0f*fA01+fA11;
                if ( fNumer >= fDenom )
		  {
                    fS = 1.0;
                    fT = 0.0;
                    fSqrDist = fA00+(2.0)*fB0+fC;
		  }
                else
		  {
                    fS = fNumer/fDenom;
                    fT = 1.0 - fS;
                    fSqrDist = fS*(fA00*fS+fA01*fT+(2.0)*fB0) +
		      fT*(fA01*fS+fA11*fT+(2.0)*fB1)+fC;
		  }
	      }
	  }
      }

    return fSqrDist; // FIXME: abs() or not ?
  }

  double distPointTriangle(const Vec3d& point,  const Triangle& tr) {
    return sqrt(squareDistPointTriangle(point, tr));
  }

  int factorial(int n)
  {
    if(n==1) return 1;
    if(n==2) return 2;
    return n*factorial(n-1);
  }

  
  void bezier3(double& z,
	       const double z1,
	       const double z2,
	       const double z3,
	       const double z4,
	       double t)
  {
    double tm1 = 1 - t;
    double tm13 = tm1 * tm1 * tm1;
    double t3 = t * t * t;
    
    z = tm13*z1 + 3*t*tm1*tm1*z2+ 3*t*t*tm1*z3 + t3*z4;
  }

  // load an obj file into an LCModel class. Supports texture coordinates
  bool importOBJ(LCModel* pModel, const char* pFilename) {
    
    
    FILE *pFile = fopen(pFilename,"rt");
    if(pFile == NULL) {
      return false;
    }
    fseek(pFile,0,SEEK_SET);
    char pLine[512];
    Vec3d* pt;
    
    
    while(fgets(pLine,512,pFile))
    {
      if(pLine[0] == 'v' && pLine[1] == ' ')
      {
        float x,y,z;
        if(sscanf(pLine,"v %f %f %f",&x,&y,&z) == 3)
        {
          pt = new Vec3d(x,y,z);
          pModel->addPoint(pt);
        }
      }
      else if(pLine[0] == 'v' && pLine[1] == 't') {
        // load texture co-ordinates, assume mapping and ordering is 1:1 with vertices
        float x,y,z;
        z = 0.0f;
        if(sscanf(pLine,"vt %f %f",&x,&y) == 2)
        {
          pt = new Vec3d(1.0 - x,y,z); // flip x because TC seem to be mirrored with respect to 3d garment
          pModel->addTextureCoord(pt);
        }
        
      }
    } // end while
    
    fclose(pFile);
    pFile = fopen(pFilename,"rt");
    
    Triangle *tr;
    IndexedTriangle *itr;
    fseek(pFile,0,SEEK_SET);
    int triCount=0;
    int fLineCount=0;
    
    while(fgets(pLine,512,pFile))
    {
      if(pLine[0] == 'f')
      {
        fLineCount++;
        int v1,v2,v3;
        int vt1,vt2,vt3;
        // The case for no texture co-ordinates
        if(sscanf(pLine,"f %d %d %d",&v1,&v2,&v3) == 3)
        {
          // the indexes in the file begin at one, not zero, hence the adjustment below
          tr = new Triangle(pModel->points()[v1 - 1],
                            pModel->points()[v2 - 1],
                            pModel->points()[v3 - 1]
                          );
            pModel->addTriangle(tr);
            
            // used indexed triangles to track the mapping between vertices and tc
            itr = new IndexedTriangle(v1 - 1,
                              v2 - 1,
                              v3 - 1
                             );
            pModel->addIndexedTriangle(itr);
            triCount++;
        }
        else if(sscanf(pLine,"f %d/%d %d/%d %d/%d",&v1,&vt1,&v2,&vt2,&v3,&vt3) == 6)
        {
          // the indexes in the file begin at one, not zero, hence the adjustment below
          
          tr = new Triangle(pModel->points()[v1 - 1],
                            pModel->points()[v2 - 1],
                            pModel->points()[v3 - 1]
                          );
          pModel->addTriangle(tr);
          
          // used indexed triangles to track the mapping between vertices and tc
          // the same index can be used on the vertex and texcord arrays
          itr = new IndexedTriangle(v1 - 1,
                                    v2 - 1,
                                    v3 - 1
                                   );
          pModel->addIndexedTriangle(itr);
          triCount++;

          
        }
      }
    } // end while
    
    std::cout << "F Line count: " << fLineCount << std::endl;
    std::cout << "Triangle count: " << triCount << std::endl;
    fclose(pFile);
    return true;
  }
  
   
  
      // load an obj file into an LCModel class supporting atlases. 
  // We are reading a file produced by Dan Julius.
  // We have to handle: Arbitrary rotation of patterns
  //                    TC indexing different to V indexing, but there is a 1:1 mapping to be exploited for lookups
  //                    The TCs are scaled between 0 and 1, we want a 'real' sized pattern so it needs to be scaled.
  // FIXME need to separate out the components of the garment
  bool importOBJGarment(Garment* pGarment, const char* pFilename, DockPattern* dockpattern)
  {
    
    
    FILE *pFile = fopen(pFilename,"rt");
    if(pFile == NULL) {
      return false;
    }
    fseek(pFile,0,SEEK_SET);
    char pLine[512];
    Vec3d* pt;
    GarmentAtlas* currentAtlas;
    std::map<std::string, int> groupCounter; // store the number of vertices in each group
    std::vector<std::string> groups; // store the group names in the order they are found
    std::string currentGroupName; // Assume there are always group names in file
    std::vector<Vec3d*> tmpV; // hold temporarily so we can split into groups
    std::vector<Vec3d*> tmpTC; // temporarily store the texture co-ordinates in original ordering - so we can finally add them to the model in the correct ordering.
    std::map<unsigned, unsigned> tcLookupFromVI; // Use vertex index as key and tc index as value.
    
    // Read file in two passes. First pass: Determine Garment structure, group names and temporarily store vertices and texture co-ords
    // in given order. Use a groupCounter to determine which vertex index ranges belong to which groups.
    
    while(fgets(pLine,512,pFile))
    {
      // Files should start with hierarchy information encoded in the comments
      // Create a Garment using these. First hierarchy line encountered defines the root of the hierarchy.
      // Every polygon in a garment control mesh should be connected.
      if(pLine[0] == '#' && pLine[1] == 'h' && pLine[2] == ' ')
      {
          std::string sname(pLine);
          sname.erase(0, 3); // remove first 3 characters
          
          
          
          typedef boost::tokenizer<boost::char_separator<char> > 
              tokenizer;
          
          boost::char_separator<char> sep(" \n");

          tokenizer tok(sname, sep);
          
          
          
          
          std::cout << "TOKENS" << std::endl;
          tokenizer::iterator beg=tok.begin();

          
          if(*beg == "garment") {
            ++beg; // advance to garment name
            pGarment->name = *beg;
            continue; // no more info required from this line
          }
          else if(*beg == "garment_section")
          {
            ++beg; // advance to section name
            GarmentSection* pGarmentSection = new GarmentSection();
            pGarmentSection->name = *beg;
            pGarment->addSection( pGarmentSection->name, pGarmentSection );
            std::cout << "Added GarmentSection: " << pGarmentSection->name << std::endl;
            
            for(tokenizer::iterator p=++beg; p!=tok.end();++p) // rest of tokens are atlas names, create them
            {
              GarmentAtlas *pAtlas = new GarmentAtlas();
              LCModel *pModel = new LCModel();
              pAtlas->_lcmodel = pModel;
              pAtlas->name = *p;
              pGarmentSection->addAtlas( pAtlas->name, pAtlas );
              std::cout << "Added GarmentAtlas: " << pAtlas->name << std::endl;
            }
          }
          
      } // end of processing hierarchy information line
      
      // Look to see which group we are in. 
      if(pLine[0] == 'g' && pLine[1] == ' ')
      {
        char name[30];
        if(sscanf(pLine,"g %s",name) == 1) {
          std::string sname(name);
          std::cout << "Reading section: " << sname << std::endl;
          groups.push_back(sname);
          currentGroupName = sname;
          groupCounter[sname]=0; // initialise the counter
        }
      }
      
      
      
      if(pLine[0] == 'v' && pLine[1] == ' ')
      {
        float x,y,z;
        if(sscanf(pLine,"v %f %f %f",&x,&y,&z) == 3)
        {
          pt = new Vec3d(x,y,z);
          tmpV.push_back(pt);
          groupCounter[currentGroupName] += 1;
        }
      }
      else if(pLine[0] == 'v' && pLine[1] == 't') {
        // load texture co-ordinates, assume mapping and ordering is 1:1 with vertices
        float x,y,z;
        z = 0.0f;
        if(sscanf(pLine,"vt %f %f",&x,&y) == 2)
        {
          pt = new Vec3d(1.0 - x,y,z); // flip x because TC seem to be mirrored with respect to 3d garment
          tmpTC.push_back(pt); // store temporarily in original order
        }
        
      }
    } // end while
    
    fclose(pFile);
    
    
    // Second pass through file. Add face information to Atlases (which correspond to groups in the file). 
    // Also build lookup table for finding correct texture co-ordinates for each vertex.
    
    
    
    std::cout << "Whole file vertex count: " << tmpV.size() << std::endl;
    std::cout << "Whole file texture coord count: " << tmpTC.size() << std::endl;
    
    pFile = fopen(pFilename,"rt");
    
    Triangle *tr;
    IndexedTriangle *itr;
    fseek(pFile,0,SEEK_SET);
    int triCount=0;
    int fLineCount=0;
    int currentGroupNumber = -1; // Will be incremented to 0 on first line of file
    int offset = 0;
    
    while(fgets(pLine,512,pFile))
    {
      int v1, v2, v3;
      int vt1, vt2, vt3;
      
      
      // Identify groups again, determine offsets into vertex indices per group
      if(pLine[0] == 'g' && pLine[1] == ' ')
      {
        char name[30];
        if(sscanf(pLine,"g %s",name) == 1) 
        {
          std::string sname(name);
          std::cout << "Reading group: " << sname << std::endl;
          
          currentGroupName = sname;
          currentGroupNumber++;
          // Determine correct section and atlas to populate
          
          
          
          GarmentAtlas *pAtlas = pGarment->getAtlas ( sname );
          assert(pAtlas!=0); // should be found unless hierarchy information was corrupt - Do group names in file correspond to hierarchy labels supplied at the beginning of the file you loaded?
          
          currentAtlas = pAtlas;
          
          if(currentGroupNumber > 0)
          {
          // offset by the number of vertices in the previous groups to get correct index
            offset += groupCounter[groups[ currentGroupNumber - 1] ];
          }
          
        }
      }
      
      
      if(sscanf(pLine,"f %d/%d %d/%d %d/%d",&v1,&vt1,&v2,&vt2,&v3,&vt3) == 6)
      {
        fLineCount++;

        v1 -= 1; // correct for starting index
        v2 -= 1;
        v3 -= 1;
                  
        vt1 -= 1; // correct for starting index
        vt2 -= 1;
        vt3 -= 1;

          
          // The range of indexes in tmpV cover all the groups in the file. But the indexes we
          // add to the atlas will only span the number of vertices in the group we require
          // hence shifting the index ranges below
          
        tr = new Triangle(tmpV[v1 ] ,
                          tmpV[v2 ] ,
                          tmpV[v3 ] 
                         );
        
        currentAtlas->_lcmodel->addTriangle(tr);

          // used indexed triangles to track the mapping between vertices and tc
          // the same index can be used on the vertex and texcord arrays
          
          // The range of indexes read from file cover all the groups in the file. But the indexes we
          // add to the atlas will only span the number of vertices in the group we require
          // hence shifting the index ranges below

        
        itr = new IndexedTriangle(v1- offset,
                                  v2- offset,
                                  v3- offset
                                 );
          
          
          
        currentAtlas->_lcmodel->addIndexedTriangle(itr);
        //pModel->addIndexedTriangle(itr);
         
          
          // Use vt1, vt2, vt3 to build a map which can be used to populate the tex-coords
          // in vertex order
          
        tcLookupFromVI[ v1 ] = vt1;
        tcLookupFromVI[ v2 ] = vt2;
        tcLookupFromVI[ v3 ] = vt3;
        triCount++;
          //std::cout << "Written lookup" << std::endl;
          
      }
      
    } // end while
    fclose(pFile);
    std::cout <<"Closed file" << std::endl;
    
    

    unsigned indexOffset = 0;
    // Populate the models with the Vs and TCs in the same order as the Vs. 
    for(int gi = 0; gi < groups.size(); gi++) 
    {
      
      
      unsigned upperBound = 0;
      if(gi>0) {
        indexOffset += groupCounter[ groups[ gi - 1] ];
      }
      
      upperBound = groupCounter[ groups[gi] ] + indexOffset;
    
      for(unsigned i = indexOffset; i < upperBound; i++) 
      {
        pGarment->getAtlas( groups[gi] )->_lcmodel->addPoint( tmpV.at(i)  );
        unsigned tcIndex = tcLookupFromVI[ i ];
        pGarment->getAtlas( groups[gi] )->_lcmodel->addTextureCoord( tmpTC.at( tcIndex ) );
      }
    
    }
    
    pGarment->printHierarchy();
    
    for(int gi = 0; gi < groups.size(); gi++) {
      // Now done via user interaction in Pattern
      //pGarment->getAtlas( groups[gi] )->_lcmodel->correctAtlasRotation();
        // Scale so edge Normals are correct when calculated next
      //std::cout << "Min: " << pGarment->getAtlas( groups[gi] )->_lcmodel->tcbbMin() << std::endl;
      //pGarment->getAtlas( groups[gi] )->scaleAtlasFromMesh();
  // So rotations are performed about the centre of the pattern
      //std::cout << "Min: " << pGarment->getAtlas( groups[gi] )->_lcmodel->tcbbMin() << std::endl;
      pGarment->getAtlas( groups[gi] )->centreModelOnOrigin();
      //pGarment->getAtlas( groups[gi] )->populateAtlasFromModel(); // Can't call this until user has corrected rotations
    }
    
    
    
    std::cout << "F Line count: " << fLineCount << std::endl;
    std::cout << "Triangle count: " << triCount << std::endl;
    
    QListView* qList = ((DockPatternWidget*) dockpattern->widget())->listView1;
    qList->clear();
    
    
    for(int i=0; i<groups.size();i++)
    {
      std::cout << "Group: " << groups[i] <<
          ", VertexCount: " << groupCounter[ groups[i] ] << std::endl;
      
      
      new QListViewItem( qList, groups[i],
                         QString::number(pGarment->getAtlas(groups[i])->_controlMesh.rows),
                         QString::number(pGarment->getAtlas(groups[i])->_controlMesh.cols)
                       );
    }
    
    std::vector<GarmentSection*>::iterator gsi, gsi_end;
    gsi= pGarment->sectionList.begin();
    gsi_end = pGarment->sectionList.end();
    for(;gsi!=gsi_end;gsi++) 
    {
      //(*gsi)->populateControlMesh(); Can't call until rotation corrected
      
    }
    
    for(int gi = 0; gi < groups.size(); gi++) {
      pGarment->getAtlas( groups[gi] )->populateAtlasFromModel();
    }
    pGarment->determineCoincidentPoints(); // If possible

    return true;
  } // end of importOBJGarment
  

  
} // end of namespace utils
