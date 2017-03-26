//
//  Filename         : layer.h
//  Author           : Emmanuel Turquin
//  Purpose          : Class representing a layer of chains for the Canvas.
//  Date of creation : 04/26/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  _LAYER_H
# define _LAYER_H

# include "vectypes.h"
# include "chain.h"
# include "collision.h"
# include "config_canvas.h"
#include <iostream>

class Grid;
class Canvas;

class Layer
{
  friend std::ostream& operator<<(std::ostream& out, const Layer& layer);
  
public:

  typedef std::list<Chain *>	ChainList;

  Layer(Canvas *parent,
	const Vec2d& min,
	const Vec2d& max,
	int sampling = config::CANVAS_SAMPLING_INT,
	int snapping = config::CANVAS_SNAPPING_INT,
	double splitting_threshold = config::CANVAS_SPLITTING_THRESHOLD);
  Layer(const Layer& layer);
  ~Layer();

  void clear();

  void firstPoint(double x, double y, double ratio);

  void intermediatePoint(double x, double y, double ratio);

  void lastPoint(double x, double y, double ratio);
  
  typedef enum { GARMENT, GFOLD, SEAM } LayerType;
  LayerType _ltype; // This warrents another class but time is of the essence 

  Stroke *currentStroke() {
    return _current_stroke;
  }

  ChainList *chains() {
    return &_chains;
  }

  const ChainList *chains() const {
    return &_chains;
  }

  void setSampling(int sampling) {
    _sampling = sampling;
  }

  int sampling() const {
    return _sampling;
  }

  void setSnapping(int snapping) {
    _snapping = snapping;
  }

  int snapping() const {
    return _snapping;
  }

  void setSplittingThreshold(double splitting_threshold) {
    _splitting_threshold = splitting_threshold;
  }

  double splittingThreshold() const {
    return _splitting_threshold;
  }

private:

  void curvatureTreatment(Stroke *stroke);

  bool collisionTreatment();

  ChainList		_chains;

  Chain			*_current_chain_begin;
  Chain			*_current_chain_end;
  Stroke		*_current_stroke;
  counted_ptr<Point>	_current_point;

  Grid			*_grid;
  Collision::List	_collisions;

  int			_sampling;
  int			_snapping;
  double		_splitting_threshold;

  Canvas		*_parent;
};

#endif // _LAYER_H
