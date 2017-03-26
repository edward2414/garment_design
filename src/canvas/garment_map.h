//
//  Filename         : garment_map.h
//  Author           : Emmanuel Turquin
//  Purpose          : Store various maps to be used for the cloth surface.
//  Date of creation : 06/01/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  DIST_BUFFER_H
# define DIST_BUFFER_H

# include "vectypes.h"
# include "mattypes.h"
# include "config_canvas.h"

class Chain;

class GarmentMap
{
public:

  GarmentMap(unsigned x = config::CANVAS_GM_X,
	     unsigned y = config::CANVAS_GM_Y);

  ~GarmentMap();

  typedef enum {
    UNKNOWN,
    IN,
    IN_OUT,
    OUT,
    BORDER,
    BORDER_EXT
  } Type;

  // Test altering the garment using a gaussian along a stroke
  void gfold(Stroke *stroke, bool front);
  
  
  inline unsigned sizeX() const
  {
    return(_size_x);
  }

  inline unsigned sizeY() const
  {
    return(_size_y);
  }

  const Vec2d& bbMin() const {
    return _bb_min;
  }

  const Vec2d& bbMax() const {
    return _bb_max;
  }

  const Vec2d bbSize() const {
    return _bb_max - _bb_min;
  }

  double minDist() const {
    return _min_dist;
  }

  double maxDist() const {
    return _max_dist;
  }

  double meanDist() const {
    return _mean_dist;
  }

  double minZ() const {
    return _min_z;
  }

  double maxZ() const {
    return _max_z;
  }

  double z(unsigned i, unsigned j) const {
    if (!_z)
      return 0;
    return _z[i * _size_y + j];
  }

  double dist(unsigned i, unsigned j) const {
    if (!_dist)
      return 0;
    return _dist[i * _size_y + j];
  }

  Type type(unsigned i, unsigned j) const {
    if (!_type)
      return UNKNOWN;
    return _type[i * _size_y + j];
  }

  const Vec2d* borderCoord(unsigned i, unsigned j) const {
    if (!_border_coord)
      return 0;
    return _border_coord[i * _size_y + j];
  }

  void compute(bool front, const Chain* chain, const Mat33d& diffusion_mask, const Mat33d& smoothing_mask,
	       unsigned iter_nb = 100, double stop_threshold = 0);
  
  void renderGL();

private:

  double applyMask(unsigned i, unsigned j, const double *array, const Mat33d& mask,
		   Type in = IN, Type out = OUT) const;

  void propagateType(unsigned i, unsigned j);

  unsigned	_size_x;
  unsigned	_size_y;
  Vec2d		_bb_min;
  Vec2d		_bb_max;
  double	_min_dist;
  double	_max_dist;
  double	_mean_dist;
  double	_min_z;
  double	_max_z;
  double	*_dist;
  double	*_z;
  Type		*_type;
  bool      *_clamp;
  Vec2d		**_border_coord;
  vecmat::Vector<double, 5> **pointGrid; // held here so we can render the result
  int pointGridSize;
};

#endif // DIST_BUFFER_H
