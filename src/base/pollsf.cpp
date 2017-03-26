#include "NumMeth.h"

// JDW Cut down version of numerical methods pollsf

double inv(Matrix a, Matrix& aInv);

void pollsf( Matrix x, Matrix y, int M, Matrix& a_fit) {
// Function to fit a polynomial to data
// Inputs 
//   x       Independent variable
//   y       Dependent variable

//   M       Number of parameters used to fit data
// Outputs
//   a_fit   Fit parameters; a(1) is intercept, a(2) is slope

  //* Form the vector b and design matrix A
  int i, j, k, N = x.nRow();
  Matrix b(N), A(N,M);
  for( i=1; i<=N; i++ ) {
   b(i) = y(i);
   for( j=1; j<=M; j++ )
     A(i,j) = pow(x(i),(double)(j-1));  
  }

  //* Compute the correlation matrix C 
  Matrix C(M,M), Cinv(M,M);
  for( i=1; i<=M; i++ ) {   // (C inverse) = (A transpose) * A
    for( j=1; j<=M; j++ ) {   
      Cinv(i,j) = 0.0;
      for( k=1; k<=N; k++ )
        Cinv(i,j) += A(k,i)*A(k,j);
    }
  }
  inv( Cinv, C );  // C = ( (C inverse) inverse)
  
  //* Compute the least squares polynomial coefficients a_fit
  for( k=1; k<=M; k++ ) {
    a_fit(k) = 0.0;
    for( j=1; j<=M; j++ )
     for( i=1; i<=N; i++ )
       a_fit(k) += C(k,j) * A(i,j) * b(i);
  }

}
