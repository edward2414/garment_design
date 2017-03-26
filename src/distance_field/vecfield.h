//
//  Filename         : vec_field.h
//  Author           : Jamie Wither
//  Purpose          : Store a vector field along the lines of VolData
//  Date of creation : 18/4/2006
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  VECDATA_H
# define VECDATA_H

# include <fstream>
# include <vector>
# include <stdexcept>


// =========================== VolumeData class ===============================

template <typename T>
class VecData
{
public:
  
  typedef T valuetype;

  VecData(unsigned x = 0, unsigned y = 0, unsigned z = 0)
  {
    _x = x;
    _y = y;
    _z = z;
    if (_x && _y && _z)
    {
      _val.resize(_x * _y * _z);
    }
  }

  inline unsigned sizeX() const
  {
    return(_x);
  }

  inline unsigned sizeY() const
  {
    return(_y);
  }

  inline unsigned sizeZ() const
  {
    return(_z);
  }

  inline unsigned size() const
  {
    return(_x * _y * _z);
  }

  const T getVal(unsigned x, unsigned y, unsigned z) const
  {
	  return val(x,y,z);
  }

  const T getVal(double x, double y, double z) const
  {
	  return valTrilinear(x,y,z);
  }

  void setVal(unsigned x, unsigned y, unsigned z, T value) {
	  val(x,y,z) =  value;
  }

protected:

  T valNearest(double x, double y, double z) const
    {
      return val(round(x), round(y), round(z));
    }

  T valTrilinear(double x, double y, double z) const
    {
      int i = static_cast<int>(x), j = static_cast<int>(y), k = static_cast<int>(z);

      if (i < 0)
	i = 0;
      else if (i >= static_cast<int>(_x) - 1)
	i = _x - 2;
      if (j < 0)
	j = 0;
      else if (j >= static_cast<int>(_y) - 1)
	j = _y - 2;
      if (k < 0)
	k = 0;
      else if (k >= static_cast<int>(_z) - 1)
	k = _z - 2;
      
      T a11, a12, a21, a22, b11, b12, b21, b22;
      T a, a1, a2, b, b1, b2, c;

      a11 = val(i, j, k);
      a12 = val(i + 1, j, k);
      a21 = val(i, j + 1, k);
      a22 = val(i + 1, j + 1, k);
      b11 = val(i, j, k + 1);
      b12 = val(i + 1, j, k + 1);
      b21 = val(i, j + 1, k + 1);
      b22 = val(i + 1, j + 1, k + 1);

      a1 = bilinear(a11, a12, x - i);
      a2 = bilinear(a21, a22, x - i);
      a = bilinear(a1, a2, y - j);
      b1 = bilinear(b11, b12, x - i);
      b2 = bilinear(b21, b22, x - i);
      b = bilinear(b1, b2, y - j);
      c = bilinear(a, b, z - k);

      return c;
    }

  inline int round(double f) const
  {
    return ((int)(f + 0.5));
  }

  inline T bilinear(T a, T b, double alpha) const
  {
    return a + alpha * (b - a);
  }

  T& val(int x, int y, int z)
  {
    if (_val.empty() ||
	x < 0 || y < 0 || z < 0 || x >= (int)_x || y >= (int)_y || z >= (int)_z)
      throw(std::domain_error("This volume data does not contain such value."));
    return _val[z * _y * _x + x * _y + y];
  }

  T val(int x, int y, int z) const
  {
    if (_val.empty())
      throw(std::domain_error("This volume data does not contain such value."));
    if (x < 0 || y < 0 || z < 0 || x >= (int)_x || y >= (int)_y || z >= (int)_z)
      return 0;
    return _val[z * _y * _x + x * _y + y];
  }

  unsigned _x;
  unsigned _y;
  unsigned _z;
  std::vector<T> _val;
};

#endif // VECDATA_H
