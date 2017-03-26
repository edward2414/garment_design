//
// C++ Implementation: Parabola
//
// Description: 
//
//
// Author: Jamie Wither <wither@stalactite>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "Parabola.h"
#include "Matrix.h"


extern void pollsf( Matrix x, Matrix y, int M, Matrix& a_fit);

double Parabola::yInterpolate( Parabola y0, Parabola y1, double alpha, double beta )
{
  return (1-alpha)*y0.yBeta(beta) + alpha*y1.yBeta(beta);
}

double Parabola::xFromBeta( double beta ) 
{
  double xMin = fit_points[0].x();
  double xRange = fit_points[2].x() - xMin;
  return (xMin + beta*xRange);
}

double Parabola::yBeta( double beta ) 
{
  // calculate y within the the defined range of the parabola
  double xMin = fit_points[0].x();
  double xRange = fit_points[2].x() - xMin;
  return y( xFromBeta( beta) );
}

// calculate y given x
double Parabola::y(double x) 
{
  return ( (x*x)*coefficients[2] + x*coefficients[1] + coefficients[0] );
}

// Fit a parabola to the fit_points
void Parabola::makeFit() {

  Matrix xVal=Matrix(3,1);
  Matrix yVal=Matrix(3,1);
  for(int i=0;i<3;i++) {
    xVal(i+1,1) = fit_points[i].x();
    yVal(i+1,1) = fit_points[i].y();
  }

  // result in here
  Matrix fit(3);
  
  pollsf(xVal, yVal, 3, fit);
  for(int i=0;i<3;i++) {
    coefficients[i] = fit(i+1);
  }
}