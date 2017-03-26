//
// C++ Interface: Parabola
//
// Description: 
//
//
// Author: Jamie Wither <wither@stalactite>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef  PARABOLA_H
#define PARABOLA_H

#include "vectypes.h"

class Parabola {

      
  public:
    Parabola() {};
    ~Parabola() {};
    // The left most, middle and right most points of a horizontal edgesection to use when fitting a parabola
    Vec3d fit_points[3];
    // parabola coefficients in reverse order a+bx+cx^2
    double coefficients[3]; // a, b, c
    
    void makeFit();
    
    // calculate y(x)
    double y(double x);
    // calculate y where beta is a proportion of the defined range. beta in [0-1]
    double yBeta(double beta);
    // given a proportion of the range return the corresponding x value;
    double xFromBeta(double beta);
    // given a proportion value for y (alpha between 0 and 1) and x (beta [0-1]) produce the point from the parabola which interpolates the two
    static double yInterpolate( Parabola y0, Parabola y1, double alpha, double beta );
};

#endif // PARABOLA_H
