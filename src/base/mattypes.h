//
//  Filename         : mattypes.h
//  Author(s)        : Emmanuel Turquin
//  Purpose          : Matrix instantiated types for easier use.
//  Date of creation : 12/06/2003
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  MATTYPES_H
# define MATTYPES_H

# include <vecmat.h>

// 2D.

typedef vecmat::SquareMatrix<char, 3>		Mat33c;
typedef vecmat::SquareMatrix<unsigned char, 3>	Mat33uc;
typedef vecmat::SquareMatrix<short, 3>		Mat33s;
typedef vecmat::SquareMatrix<unsigned short, 3>	Mat33us;
typedef vecmat::SquareMatrix<int, 3>		Mat33i;
typedef vecmat::SquareMatrix<unsigned int, 3>	Mat33ui;
typedef vecmat::SquareMatrix<float, 3>		Mat33f;
typedef vecmat::SquareMatrix<double, 3>		Mat33d;

// 3D.

typedef vecmat::SquareMatrix<char, 4>		Mat44c;
typedef vecmat::SquareMatrix<unsigned char, 4>	Mat44uc;
typedef vecmat::SquareMatrix<short, 4>		Mat44s;
typedef vecmat::SquareMatrix<unsigned short, 4>	Mat44us;
typedef vecmat::SquareMatrix<int, 4>		Mat44i;
typedef vecmat::SquareMatrix<unsigned int, 4>	Mat44ui;
typedef vecmat::SquareMatrix<float, 4>		Mat44f;
typedef vecmat::SquareMatrix<double, 4>		Mat44d;

#endif // MATTYPES_H
