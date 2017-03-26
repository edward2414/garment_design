//
//  Filename         : voldata.h
//  Author           : Emmanuel Turquin
//  Purpose          : Volume data class.
//  Date of creation : 04/20/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  VOLDATA_H
# define VOLDATA_H

# include <fstream>
# include <vector>
# include <stdexcept>

// =========================== VolumeData class ===============================

template <typename T>
class VolumeData
{
public:
  
  typedef T valuetype;

  VolumeData(unsigned x = 0, unsigned y = 0, unsigned z = 0)
  {
    _x = x;
    _y = y;
    _z = z;
    _max_val = 0;
    if (_x && _y && _z)
      _val.resize(_x * _y * _z);
  }

  VolumeData(std::istream& is, unsigned x, unsigned y, unsigned z)
  {
    _max_val = 0;
    read(is, x, y, z);
  }

  template <class T1>
  VolumeData(std::istream& is, unsigned x, unsigned y, unsigned z)
  {
    _max_val = 0;
    read<T1>(is, x, y, z);
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

  inline T max(bool compute = false)
  {
    if (compute)
      computeMax();
    return(_max_val);
  }

  void computeMax()
  {
    _max_val = 0;
    for (typename std::vector<T>::const_iterator it = _val.begin(); it != _val.end(); ++it)
      if (*it > _max_val)
	_max_val = *it;
  }

  void classify(T i1, T i2, T j1, T j2)
  {
    if (i1 > i2)
      swap (i1, i2);
    for (typename std::vector<T>::iterator it = _val.begin(); it != _val.end(); ++it)
      if (*it)
	*it = (T)((double)(j2 - j1) / (i2 - i1) * *it) + j1;
    computeMax();
  }

  inline unsigned size() const
  {
    return(_x * _y * _z);
  }

  T& val(unsigned int i)
  {
    if (_val.empty() || i >= size())
      throw(std::domain_error("This volume data does not contain such value."));
    return _val[i];
  }

  T val(unsigned int i) const
  {
    if (_val.empty())
      throw(std::domain_error("This volume data does not contain such value."));
    if (i >= size())
      return 0;
    return _val[i];
  }

  T valNearest(double x, double y, double z) const
    {
      return val(round(x), round(y), round(z));
    }

  double valTrilinear(double x, double y, double z) const
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
      double a, a1, a2, b, b1, b2, c;

      a11 = val(i, j, k);
      a12 = val(i + 1, j, k);
      a21 = val(i, j + 1, k);
      a22 = val(i + 1, j + 1, k);
      b11 = val(i, j, k + 1);
      b12 = val(i + 1, j, k + 1);
      b21 = val(i, j + 1, k + 1);
      b22 = val(i + 1, j + 1, k + 1);
	
      if (!a11 && !a12 && !a21 && !a22 && !b11 && !b12 && !b21 && !b22)
	return 0;
      
      a1 = bilinear(a11, a12, x - i);
      a2 = bilinear(a21, a22, x - i);
      a = bilinear(a1, a2, y - j);
      b1 = bilinear(b11, b12, x - i);
      b2 = bilinear(b21, b22, x - i);
      b = bilinear(b1, b2, y - j);
      c = bilinear(a, b, z - k);

      return c;
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

  template<class T1>
  void read(std::istream& is, unsigned x, unsigned y, unsigned z)
  {
    _x = x;
    _y = y;
    _z = z;

    // Read data
    unsigned size = x * y * z;
    
    T1* val = new T1[size];

    is.read((char*)val, size * sizeof(T1));
    if (is.eof())
      throw(std::domain_error("Wrong data size (volume actually smaller)."));
    is.get();
    if (!is.eof())
      throw(std::domain_error("Wrong data size (volume actually bigger)."));
    _val.resize(size);

    
    for (unsigned i = 0; i < size; ++i)
      {
	if (val[i] > _max_val)
	  _max_val = val[i];
	_val[i] = val[i];
      }
    delete val;
  }

  void read(std::istream& is, unsigned x, unsigned y, unsigned z)
  {
    read<T>(is, x, y, z);
  }

  void write(std::ostream& os) const
  {
    unsigned size = _x * _y * _z;

    for (unsigned i = 0; i < size; ++i)
      os << _val[i] << " ";
  }

protected:
  
  inline int round(double f) const
  {
    return ((int)(f + 0.5));
  }

  inline double bilinear(double a, double b, double alpha) const
  {
    return a + alpha * (b - a);
  }

  unsigned _x;
  unsigned _y;
  unsigned _z;
  std::vector<T> _val;
  T _max_val;
};

#endif // VOLDATA_H
