//
//  Filename         : indexed_triangle.h
//  Author(s)        : Jamie Wither
//  Purpose          : Triangle which can be used to lookup vertices or texcoords via indexing (in 3D).
//                   : indexes into LCModel point/tc lists
//  Date of creation : 13 Dec 2005
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  INDEXED_TRIANGLE_H
# define INDEXED_TRIANGLE_H


class IndexedTriangle
{
public:

  IndexedTriangle( int a,  int b,  int c) {
    _a = a;
    _b = b;
    _c = c;
  }

  ~IndexedTriangle() {}

   int indexA()  {
    return _a;
  }
    
   int indexB()  {
    return _b;
  }
  
   int indexC()  {
    return _c;
  }

  void setIndexA( int a) {
     _a = a;
  }

  void setIndexB( int b) {
    _b = b;
  }

  void setIndexC( int c) {
    _c = c;
  }

private:

   int _a;
   int _b;
   int _c;
};

#endif // INDEXED_TRIANGLE_H
