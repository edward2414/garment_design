//
//  Filename         : vecmat.h
//  Author(s)        : Emmanuel Turquin
//  Purpose          : Vectors and Matrices definition and manipulation.
//  Date of creation : 12/06/2003
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  VECMAT_H
# define VECMAT_H

# include <cmath>
# include <vector>
# include <iostream>

namespace vecmat {

  namespace internal {

    template <bool B>
    struct is_false {};

    template <>
    struct is_false<false> {
      static inline void ensure() {}
    };

  } // end of namespace internal

  //
  //  Vector class
  //    - T: value type
  //    - N: dimension
  //
  /////////////////////////////////////////////////////////////////////////////

  template <class T, unsigned N>
  class Vector
  {
  public:

    typedef T value_type;
    
    // constructors

    inline Vector() {
      for (unsigned i = 0; i < N; i++)
	_coord[i] = 0;
    }

    ~Vector() {
      internal::is_false<(N == 0)>::ensure();
    }

    template <class U>
    explicit inline Vector(const U tab[N]) {
      for (unsigned i = 0; i < N; i++)
	_coord[i] = (T)tab[i];
    }

    template <class U>
    explicit inline Vector(const std::vector<U>& tab) {
      for (unsigned i = 0; i < N; i++)
	_coord[i] = (T)tab[i];
    }

    template <class U>
    explicit inline Vector(const Vector<U, N>& v) {
      for (unsigned i = 0; i < N; i++)
        _coord[i] = (T)v[i];
    }

    // accessors

    inline value_type  operator[](const unsigned i) const {
      return _coord[i];
    }

    inline value_type& operator[](const unsigned i) {
      return _coord[i];
    }

    static inline unsigned dim() {
      return N;
    }
    
    // various useful methods

    inline value_type norm() const {
      return (T)sqrt(squareNorm());
    }

    inline value_type squareNorm() const {
      return (*this) * (*this);
    }

    inline Vector<T, N>& normalize() {
      value_type n = norm();
      for (unsigned i = 0; i < N; i++)
	_coord[i] /= n;
      return *this;
    }

    inline Vector<T, N>& normalizeSafe() {
      value_type n = norm();
      if (n)
	for (unsigned i=0; i < N; i++)
	  _coord[i] /= n;
      return *this;
    }

    inline Vector<T, N>& min(const Vector<T, N>& v) {
      for (unsigned i=0; i < N; i++)
	if (_coord[i]  > v._coord[i])
	  _coord[i] = v._coord[i];
      return *this;
    }

    inline Vector<T, N>& max(const Vector<T, N>& v) {
      for (unsigned i=0; i < N; i++)
	if (_coord[i]  < v._coord[i])
	  _coord[i] = v._coord[i];
      return *this;
    }

    inline const value_type* address() const {
      return _coord;
    }

    // classical operators

    template <class U>
    inline Vector<T, N>& operator=(const Vector<U, N>& v) {
      if (this != &v)
	for (unsigned i = 0; i < N; i++)
	  _coord[i] = (T)v[i];
      return *this;
    }

    template <class U>
    inline Vector<T, N>& operator+=(const Vector<U, N>& v) {
      for (unsigned i = 0 ; i < N; i++)
	_coord[i] += (T)v[i];
      return *this;
    }

    template <class U>
    inline Vector<T, N>& operator-=(const Vector<U, N>& v) {
      for (unsigned i = 0 ; i < N; i++)
	_coord[i] -= (T)v[i];
      return *this;
    }

    template <class U>
    inline Vector<T, N>& operator*=(const U r) {
      for (unsigned i = 0 ; i < N; i++)
	_coord[i] *= r;
      return *this;
    }

    template <class U>
    inline Vector<T, N>& operator/=(const U r) {
      if (r)
	for (unsigned i = 0 ; i < N; i++)
	  _coord[i] /= r;
      return *this;
    }


    inline bool operator==(const Vector<T, N>& v) const {
      for(unsigned i = 0; i < N; i++)
	if (_coord[i] != v[i])
	  return false;
      return true;
    }

    inline bool operator!=(const Vector<T, N>& v) const {
      for(unsigned i = 0; i < N; i++)
	if (_coord[i] != v[i])
	  return true;
      return false;
    }

    inline bool operator<(const Vector<T, N>& v) const {
      for (unsigned i = 0; i<N; i++) {
	if (_coord[i] < v[i])
	  return true;
	if (_coord[i] > v[i])
	  return false;
	if (_coord[i] == v[i])
	  continue;
      }
      return false;  
    }

    inline bool operator>(const Vector<T, N>& v) const {
      for (unsigned i=0; i<N; i++) {
	if(_coord[i] > v[i])
	  return true;
	if(_coord[i] < v[i])
	  return false;
	if(_coord[i] == v[i])
	  continue;
      }
      return false;  
    }

  protected:

    value_type _coord[N];
    enum {
      _dim = N,
    };
  };


  //
  //  Vec2 class (2D Vector)
  //    - T: value type
  //
  /////////////////////////////////////////////////////////////////////////////

  template <class T>
  class Vec2 : public Vector<T, 2>
  {
  public:

    typedef typename Vector<T, 2>::value_type value_type;

    inline Vec2() : Vector<T, 2>() {}

    template <class U>
    explicit inline Vec2(const U tab[2]) : Vector<T, 2>(tab) {}

    template <class U>
    explicit inline Vec2(const std::vector<U>& tab) : Vector<T, 2>(tab) {}

    template <class U>
    inline Vec2(const Vector<U, 2>& v) : Vector<T, 2>(v) {}

    inline Vec2(const value_type x,
		const value_type y = 0) : Vector<T, 2>() {
      this->_coord[0] = (T)x;
      this->_coord[1] = (T)y;
    }
    
    inline value_type x() const {
      return this->_coord[0];
    }

    inline value_type& x() {
      return this->_coord[0];
    }

    inline value_type y() const {
      return this->_coord[1];
    }

    inline value_type& y() {
      return this->_coord[1];
    }
  };


  //
  //  HVec3 class (3D Vector in homogeneous coordinates)
  //    - T: value type
  //
  /////////////////////////////////////////////////////////////////////////////

  template <class T>
  class HVec3 : public Vector<T, 4>
  {
  public:

    typedef typename Vector<T, 4>::value_type value_type;

    inline HVec3() : Vector<T, 4>() {}

    template <class U>
    explicit inline HVec3(const U tab[4]) : Vector<T, 4>(tab) {}

    template <class U>
    explicit inline HVec3(const std::vector<U>& tab) : Vector<T, 4>(tab) {}

    template<class U>
    inline HVec3(const Vector<U, 4>& v) : Vector<T, 4>(v) {}
  
    inline HVec3(const value_type sx,
		 const value_type sy = 0,
		 const value_type sz = 0,
		 const value_type s = 1) {
      this->_coord[0] = sx;
      this->_coord[1] = sy;
      this->_coord[2] = sz;
      this->_coord[3] = s;
    }
    
    template <class U>
    inline HVec3(const Vector<U, 3>& sv) {
      this->_coord[0] = (T)sv[0];
      this->_coord[1] = (T)sv[1];
      this->_coord[2] = (T)sv[2];
      this->_coord[3] = (T)1;
    }
    

    template <class U>
    inline HVec3(const Vector<U, 3>& sv,
		 const U) {
      this->_coord[0] = (T)sv[0];
      this->_coord[1] = (T)sv[1];
      this->_coord[2] = (T)sv[2];
      this->_coord[3] = (T)s;
    }
    
    inline value_type sx() const {
      return this->_coord[0];
    }

    inline value_type& sx() {
      return this->_coord[0];
    }

    inline value_type sy() const {
      return this->_coord[1];
    }

    inline value_type& sy() {
      return this->_coord[1];
    }

    inline value_type sz() const {
      return this->_coord[2];
    }

    inline value_type& sz() {
      return this->_coord[2];
    }

    inline value_type s() const {
      return this->_coord[3];
    }

    inline value_type& s() {
      return this->_coord[3];
    }

    // Acces to non-homogeneous coordinates in 3D

    inline value_type x() const {
      return this->_coord[0] / this->_coord[3];
    }

    inline value_type y() const {
      return this->_coord[1] / this->_coord[3];
    }

    inline value_type z() const {
      return this->_coord[2] / this->_coord[3];
    }
  };


  //
  //  Vec3 class (3D Vector)
  //    - T: value type
  //
  /////////////////////////////////////////////////////////////////////////////

  template <class T>
  class Vec3 : public Vector<T, 3>
  {
  public:
    
    typedef typename Vector<T, 3>::value_type value_type;

    inline Vec3() : Vector<T, 3>() {}

    template <class U>
    explicit inline Vec3(const U tab[3]) : Vector<T, 3>(tab) {}

    template <class U>
    explicit inline Vec3(const std::vector<U>& tab) : Vector<T, 3>(tab) {}

    template<class U>
    inline Vec3(const Vector<U, 3>& v) : Vector<T, 3>(v) {}
  
    template<class U>
      inline Vec3(const HVec3<U>& v) {
      this->_coord[0] = (T)v.x();
      this->_coord[1] = (T)v.y();
      this->_coord[2] = (T)v.z();
    }

    inline Vec3(const value_type x,
		const value_type y = 0,
		const value_type z = 0) : Vector<T, 3>() {
      this->_coord[0] = x;
      this->_coord[1] = y;
      this->_coord[2] = z;
    }
    
    inline value_type x() const {
      return this->_coord[0];
    }

    inline value_type& x() {
      return this->_coord[0];
    }

    inline value_type y() const {
      return this->_coord[1];
    }

    inline value_type& y() {
      return this->_coord[1];
    }

    inline value_type z() const {
      return this->_coord[2];
    }

    inline value_type& z() {
      return this->_coord[2];
    }
  };


  //
  //  Matrix class
  //    - T: value type
  //    - M: rows
  //    - N: cols
  //
  /////////////////////////////////////////////////////////////////////////////
  
  // Dirty, but icc under Windows needs this
# define _SIZE (M * N)

  template <class T, unsigned M, unsigned N>
  class Matrix
  {
  public:

    typedef T value_type;
    
    inline Matrix() {
      for (unsigned i = 0; i < _SIZE; i++)
	this->_coord[i] = 0;
    }

    ~Matrix() {
      internal::is_false<(M == 0)>::ensure();
      internal::is_false<(N == 0)>::ensure();
    }

    template <class U>
    explicit inline Matrix(const U tab[M][N]) {
      for (unsigned i = 0; i < M; i++)
	for (unsigned j = 0; j < N; j++)
	  this->_coord[i * N + j] = tab[i][j];
    }

   template <class U>
    explicit inline Matrix(const U tab[_SIZE]) {
      for (unsigned i = 0; i < _SIZE; i++)
	this->_coord[i] = tab[i];
    }

    template <class U>
    explicit inline Matrix(const std::vector<U>& tab) {
      for (unsigned i = 0; i < _SIZE; i++)
	this->_coord[i] = tab[i];
    }

    template <class U>
    inline Matrix(const Matrix<U, M, N>& m) {
      for (unsigned i = 0; i < M; i++)
	for (unsigned j = 0; j < N; j++)
	  this->_coord[i * N + j] = (T)m(i, j);
    }

    inline value_type operator()(const unsigned i, const unsigned j) const {
      return this->_coord[i * N + j];
    }

    inline value_type& operator()(const unsigned i, const unsigned j) {
      return this->_coord[i * N + j];
    }

    static inline unsigned rows() {
      return M;
    }

    static inline unsigned cols() {
      return N;
    }

    inline Matrix<T, M, N> transpose() const {
      Matrix<T, N, M> res;
      for (unsigned i = 0; i < M; i++)
	for (unsigned j = 0; j < N; j++)
	  res(j,i) = this->_coord[i * N + j];
      return res;
    }

    inline void getArray(value_type res[M][N]) const {
      for (unsigned i = 0; i < M; i++)
	for (unsigned j = 0; j < N; j++)
	  res[i][j] = this->_coord[i * N + j];
    }

    inline void getArray(value_type res[_SIZE]) const {
      for (unsigned i = 0; i < _SIZE; i++)
	res[i] = this->_coord[i];
    }

    inline const value_type* address() const {
      return this->_coord;
    }

    template <class U>
    inline Matrix<T, M, N>& operator=(const Matrix<U, M, N>& m) {
      if (this != &m)
	for (unsigned i = 0; i < M; i++)
	  for (unsigned j = 0; j < N; j++)
	    this->_coord[i * N + j] = (T)m(i, j);
      return *this;
    }

    template <class U>
    inline Matrix<T, M, N>& operator+=(const Matrix<U, M, N>& m) {
      for (unsigned i = 0; i < M; i++)
	for (unsigned j = 0; j < N; j++)
	  this->_coord[i * N + j] += (T)m(i, j);
      return *this;
    }

    template <class U>
    inline Matrix<T, M, N>& operator-=(const Matrix<U, M, N>& m) {
      for (unsigned i = 0; i < M; i++)
	for (unsigned j = 0; j < N; j++)
	  this->_coord[i * N + j] -= (T)m(i, j);
      return *this;
    }

    template <class U>
    inline Matrix<T, M, N>& operator*=(const U lambda) {
      for (unsigned i = 0; i < M; i++)
	for (unsigned j = 0; j < N; j++)
	  this->_coord[i * N + j] *= lambda;
      return *this;
    }

    template <class U>
    inline Matrix<T, M, N>& operator/=(const U lambda) {
      if (lambda)
	for (unsigned i = 0; i < M; i++)
	  for (unsigned j = 0; j < N; j++)
	    this->_coord[i * N + j] /= lambda;
      return *this;
    }

  protected:

    value_type _coord[_SIZE];
  };


  //
  //  SquareMatrix class
  //    - T: value type
  //    - N: rows & cols
  //
  /////////////////////////////////////////////////////////////////////////////

  // Dirty, but icc under Windows needs this
# define __SIZE (N * N)

  template <class T, unsigned N>
  class SquareMatrix : public Matrix<T, N, N>
  {
  public:

    typedef T value_type;
  
    inline SquareMatrix() : Matrix<T, N, N>() {}

    template <class U>
    explicit inline SquareMatrix(const U tab[__SIZE]) : Matrix<T, N, N>(tab) {}

    template <class U>
    explicit inline SquareMatrix(const std::vector<U>& tab) : Matrix<T, N, N>(tab) {}

    template <class U>
    inline SquareMatrix(const Matrix<U, N, N>& m) : Matrix<T, N, N>(m) {}

    static inline SquareMatrix<T, N> identity() {
      SquareMatrix<T, N> res;
      for (unsigned i = 0; i < N; i++)
	res(i, i) = 1;
      return res;
    }
  };


  //
  // Vector external functions
  //
  /////////////////////////////////////////////////////////////////////////////

  template <class T, unsigned N>
  inline Vector<T, N> operator+(const Vector<T, N>& v1,
				const Vector<T, N>& v2) {
    Vector<T, N> res(v1);
    res += v2;
    return res;
  }

  template <class T, unsigned N>
  inline Vector<T, N> operator-(const Vector<T, N>& v1,
				const Vector<T, N>& v2) {
    Vector<T, N> res(v1);
    res -= v2;
    return res;
  }
  template <class T, unsigned N>
  inline Vector<T, N> operator*(const Vector<T, N>& v,
				const typename Vector<T, N>::value_type r) {
    Vector<T, N> res(v);
    res *= r;
    return res;
  }

  template <class T, unsigned N>
  inline Vector<T, N> operator*(const typename Vector<T, N>::value_type r,
				const Vector<T, N>& v) {
    Vector<T, N> res(v);
    res *= r;
    return res;
  }

  template <class T, unsigned N>
  inline Vector<T, N> operator/(const Vector<T, N>& v,
				const typename Vector<T, N>::value_type r) {
    Vector<T, N> res(v);
    if (r)
      res /= r;
    return res;
  }

  // dot product
  template <class T, unsigned N>
  inline typename Vector<T, N>::value_type operator*(const Vector<T, N>& v1,
						     const Vector<T, N>& v2) {
    typename Vector<T, N>::value_type sum = 0;
    for (unsigned i = 0; i < N; i++)
      sum += v1[i] * v2[i];
    return sum;
  }

  // cross product for 3D Vectors
  template <typename T>
  inline Vec3<T> operator^(const Vector<T, 3>& v1,
			   const Vector<T, 3>& v2) {
    Vec3<T> res(v1[1] * v2[2] - v1[2] * v2[1],
		v1[2] * v2[0] - v1[0] * v2[2],
		v1[0] * v2[1] - v1[1] * v2[0]);
    return res;
  }

  // stream operator
  template <class T, unsigned N>
  inline std::ostream& operator<<(std::ostream& s,
				  const Vector<T, N>& v) {
    unsigned i;
    s << "[";
    for (i = 0; i < N - 1; i++)
      s << v[i] << ", ";
    s << v[i] << "]";
    return s;
  }


  //
  // Matrix external functions
  //
  /////////////////////////////////////////////////////////////////////////////

  template <class T, unsigned M, unsigned N>
  inline Matrix<T, M, N>
  operator+(const Matrix<T, M, N>& m1,
	    const Matrix<T, M, N>& m2) {
    Matrix<T, M, N> res(m1);
    res += m2;
    return res;
  }

  template <class T, unsigned M, unsigned N>
  inline Matrix<T, M, N>
  operator-(const Matrix<T, M, N>& m1,
	    const Matrix<T, M, N>& m2) {
    Matrix<T, M, N> res(m1);
    res -= m2;
    return res;
  }

  template <class T, unsigned M, unsigned N>
  inline Matrix<T, M, N>
  operator*(const Matrix<T, M, N>& m1,
	    const typename Matrix<T, M, N>::value_type lambda) {
    Matrix<T, M, N> res(m1);
    res *= lambda;
    return res;
  }

  template <class T, unsigned M, unsigned N>
  inline Matrix<T, M, N>
  operator*(const typename Matrix<T, M, N>::value_type lambda,
	    const Matrix<T, M, N>& m1) {
    Matrix<T, M, N> res(m1);
    res *= lambda;
    return res;
  }

  template <class T, unsigned M, unsigned N>
  inline Matrix<T, M, N>
  operator/(const Matrix<T, M, N>& m1,
	    const typename Matrix<T, M, N>::value_type lambda) {
    Matrix<T, M, N> res(m1);
    res /= lambda;
    return res;
  }

  template <class T, unsigned M, unsigned N, unsigned P>
  inline Matrix<T, M, P>
  operator*(const Matrix<T, M, N>& m1,
	    const Matrix<T, N, P>& m2) {
    unsigned i, j, k;
    Matrix<T, M, P> res;
    typename  Matrix<T, N, P>::value_type scale;
  
    for (j = 0; j < P; j++) {
      for (k = 0; k < N; k++) {
	scale = m2(k, j);
	for (i = 0; i < N; i++)
	  res(i, j) += m1(i, k) * scale;
      }
    }
    return res;
  }

  template <class T, unsigned M, unsigned N>
  inline Vector<T, M>
  operator*(const Matrix<T, M, N>& m,
	    const Vector<T, N>& v) {
  
    Vector<T, M> res;
    typename Matrix<T, M, N>::value_type scale;
  
    for (unsigned j = 0; j < M; j++) {
      scale = v[j];
      for (unsigned i = 0; i < N; i++)
	res[i] += m(i, j) * scale;
    }
    return res;
  }

  // stream operator
  template <class T, unsigned M, unsigned N>
  inline std::ostream& operator<<(std::ostream& s,
				  const Matrix<T, M, N>& m) {
    unsigned i, j;
    for (i = 0; i < M; i++) {
      s << "[";
      for (j = 0; j < N - 1; j++)
	s << m(i, j) << ", ";
      s << m(i, j) << "]" << std::endl;
    }
    return s;
  }

} // end of namespace vecmat

#endif // VECMAT_H
