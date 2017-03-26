//
//  Filename         : garment_map.cpp
//  Author           : Emmanuel Turquin
//  Purpose          : Store various maps to be used for the cloth surface.
//  Date of creation : 06/01/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include <limits>
#include <stack>
#include "utils.h"
#include "chain.h"
#include "point_utils.h"
#include "repository_pretreatment_canvas.h"
#include "repository_pretreatment_result.h"
#include "garment_map.h"

#include <GL/gl.h>
#include <GL/glu.h>
#include <math.h>
#include <algorithm>

GarmentMap::GarmentMap(unsigned x, unsigned y)
{
  _size_x = x;
  _size_y = y;
  _min_dist = 0;
  _max_dist = 0;
  _mean_dist = 0;
  _min_z = 0;
  _max_z = 0;
  _dist = 0;
  _z = 0;
  _type = 0;
  _border_coord = 0;
  _clamp = 0;
  pointGrid = 0;
  pointGridSize =0;
}

GarmentMap::~GarmentMap()
{
  if (_border_coord) {
    for (unsigned i = 0; i < _size_x; ++i)
      for (unsigned j = 0; j < _size_y; ++j)
        if (_border_coord[i * _size_y + j])
          delete _border_coord[i * _size_y + j];
    delete[] _border_coord;
  }
  delete[] _dist;
  delete[] _z;
  delete[] _type;
  delete[] _clamp;
  
  if (pointGrid) {
    for (unsigned i = 0; i < pointGridSize; ++i)
        if (pointGrid[i])
          delete pointGrid[i];
    delete[] pointGrid;
  }
}

double GarmentMap::applyMask(unsigned i, unsigned j, const double *array, const Mat33d& mask, Type in, Type out) const
{
  if (_type[i * _size_y + j] != in)
    return array[i * _size_y + j];

  unsigned min_i = i >= 1 ? i-1 : i;
  unsigned max_i = i < _size_x-1 ? i+1 : i;
  unsigned min_j = j >= 1 ? j-1 : j;
  unsigned max_j = j < _size_y-1 ? j+1 : j;

  double sum = 0;
  double denom = 0;
  for (i = min_i; i <= max_i; ++i)
    for (j = min_j; j <= max_j; ++j)
      if (_type[i * _size_y + j] != out) {
        sum += mask(2 - j + min_j, i - min_i) * array[i * _size_y + j];
        denom += mask(2 - j + min_j, i - min_i);
      }

  if (!denom) {
    return 0;
  }

  return sum / denom;
}

void GarmentMap::propagateType(unsigned i, unsigned j)
{
  std::stack<Vec2ui> coord_stack;
  coord_stack.push(Vec2ui(i, j));
  while (!coord_stack.empty()) {
    i = coord_stack.top()[0];
    j = coord_stack.top()[1];
    coord_stack.pop();    
    if (i > _size_x - 1 || j > _size_y - 1 || _type[i * _size_y + j] != UNKNOWN)
      continue;
    if (i == 0 || _type[(i-1) * _size_y + j] == OUT)
      _type[i * _size_y + j] = OUT;
    else if (j == 0 || _type[i * _size_y + (j-1)] == OUT)
      _type[i * _size_y + j] = OUT;
    else if (i == _size_x - 1 || _type[(i+1) * _size_y + j] == OUT)
      _type[i * _size_y + j] = OUT;
    else if (j == _size_y - 1 || _type[i * _size_y + (j+1)] == OUT)
      _type[i * _size_y + j] = OUT;
    else
      _type[i * _size_y + j] = IN;
    coord_stack.push(Vec2ui(i-1, j));
    coord_stack.push(Vec2ui(i, j-1));
    coord_stack.push(Vec2ui(i+1, j));
    coord_stack.push(Vec2ui(i, j+1));
  }
}

// Internal use only.
// inline double slowFastSlow(double t) {
//   return sin(M_PI*(t-0.5))*0.5+0.5;
// }

void GarmentMap::compute(bool front, const Chain* chain, const Mat33d& diffusion_mask, const Mat33d& smoothing_mask,
        		 unsigned iter_nb, double stop_threshold)
{
  if (!chain || !chain->cycle())
    return;

  Chain::StrokeList::const_iterator st;
  Chain::StrokeList::const_iterator st_end;
  Stroke::SegmentList::const_iterator seg;
  Stroke::SegmentList::const_iterator seg_end;

  // First iteration on points.
  // - compute the bbox of the chain
  // - find _min_dist and _max_dist
  // - compute _mean_dist
  
  _bb_min[0] = std::numeric_limits<Vec2d::value_type>::max();
  _bb_min[1] = std::numeric_limits<Vec2d::value_type>::max();
  _bb_max[0] = -std::numeric_limits<Vec2d::value_type>::max();
  _bb_max[1] = -std::numeric_limits<Vec2d::value_type>::max();
  _mean_dist = 0;
  _min_dist = std::numeric_limits<double>::max();
  _max_dist = 0;
  unsigned pts_nb = 0;

  st = chain->strokes()->begin();
  st_end = chain->strokes()->end();
  for (; st != st_end; ++st) {
    seg = (*st)->segments()->begin();
    seg_end = (*st)->segments()->end();
    for (; seg != seg_end; ++seg) {
      if (_bb_min[0] > (*seg)->pointA()->x())
        _bb_min[0] = (*seg)->pointA()->x();
      if (_bb_max[0] < (*seg)->pointA()->x())
        _bb_max[0] = (*seg)->pointA()->x();
      if (_bb_min[1] > (*seg)->pointA()->y())
        _bb_min[1] = (*seg)->pointA()->y();
      if (_bb_max[1] < (*seg)->pointA()->y())
        _bb_max[1] = (*seg)->pointA()->y();
      if (_min_dist > (*seg)->pointA()->distance())
        _min_dist = (*seg)->pointA()->distance();
      if (_max_dist < (*seg)->pointA()->distance())
        _max_dist = (*seg)->pointA()->distance();
      _mean_dist += (*seg)->pointA()->distance();
      ++pts_nb;
    }
  }
  


  // Ensure that all the segments are contained in the box.

  double dx = (_bb_max[0] - _bb_min[0]) / 100;
  double dy = (_bb_max[1] - _bb_min[1]) / 100;
  _bb_min[0] -= dx;
  _bb_min[1] -= dy;
  _bb_max[0] += dx;
  _bb_max[1] += dy;

  _mean_dist /= pts_nb;
  
  //std::cout << "GM: Min dist: " << _min_dist << std::endl;
  //std::cout << "GM: Max dist: " << _max_dist << std::endl;
  //std::cout << "GM: Mean dist: " << _mean_dist << std::endl;
  
  // Compute the size of the maps and allocate them.

  _size_x = static_cast<unsigned>(((_bb_max[0] - _bb_min[0]) / RepositoryPretreatmentResult::modelBBSize()[0]) * _size_x);
  _size_y = static_cast<unsigned>(((_bb_max[1] - _bb_min[1]) / RepositoryPretreatmentResult::modelBBSize()[1]) * _size_y);
  _z = new double[_size_x * _size_y];
  _dist = new double[_size_x * _size_y];
  _type   = new Type[_size_x * _size_y];
  
  // JDW when altering using folding tools which points need to be clamped? (Seams)
  _clamp = new bool[_size_x * _size_y];
  
  // JDW Used to store the shifted co-ordinates for the border pixels. Corresponds to Sect 5.3 in the paper
  _border_coord = new Vec2d*[_size_x * _size_y];

  // First iteration on the arrays:
  // - initialize the distance to _mean_dist
  // - initialize the type to UNKNOWN
  
  unsigned i, j;
  
  for (i = 0; i < _size_x; ++i)
    for (j = 0; j < _size_y; ++j) {
      _dist[i * _size_y + j] = _mean_dist;
      _type[i * _size_y + j] = UNKNOWN;
      _border_coord[i * _size_y + j] = 0;
      _clamp[i * _size_y + j] = false;
    }

  // Second iteration on points:
  // - add all the points of the chain in the array, and set their type to BORDER or BORDER_EXT
    // - shift the co-ordinates of border points and record them in an array

  Type border_type;
  st = chain->strokes()->begin();
  st_end = chain->strokes()->end();
  for (; st != st_end; ++st)
  {
    seg = (*st)->segments()->begin();
    seg_end = (*st)->segments()->end();
    for (; seg != seg_end; ++seg)
    {
      i = static_cast<unsigned>((((*seg)->pointA()->x() - _bb_min[0]) / (_bb_max[0] - _bb_min[0])) * _size_x);
      j = static_cast<unsigned>((((*seg)->pointA()->y() - _bb_min[1]) / (_bb_max[1] - _bb_min[1])) * _size_y);
      if (i && i >= _size_x)
        i = _size_x - 1;
      if (j && j >= _size_y)
        j = _size_y - 1;
      if ((*seg)->pointA().get() == (*seg)->stroke()->pointA().get())
        border_type = BORDER_EXT;
      else
        border_type = BORDER;
      
      // JDW Set cells to be clamped (first pass)
      if( (*st)->type() == Stroke::OUT)
        _clamp[i * _size_y + j] = true;
      
      // If more than one point it assigned to a cell, average the various contributions
      
      if (_type[i * _size_y + j] == BORDER_EXT)
      {
        if (border_type == BORDER)
          continue;
        _border_coord[i * _size_y + j]->x() = (_border_coord[i * _size_y + j]->x() + (*seg)->pointA()->x()) / 2;
        _border_coord[i * _size_y + j]->y() = (_border_coord[i * _size_y + j]->y() + (*seg)->pointA()->y()) / 2;
        _z[i * _size_y + j] = (_z[i * _size_y + j] + (*seg)->pointA()->z()) / 2;
        _dist[i * _size_y + j] = (_dist[i * _size_y + j] + (*seg)->pointA()->distance()) / 2;
      }
      else if (_type[i * _size_y + j] == BORDER)
      {
        if (border_type == BORDER_EXT)
        {
          _border_coord[i * _size_y + j]->x() = (*seg)->pointA()->x();
          _border_coord[i * _size_y + j]->y() = (*seg)->pointA()->y();
          _z[i * _size_y + j] = (*seg)->pointA()->z();
          _dist[i * _size_y + j] = (*seg)->pointA()->distance();
          _type[i * _size_y + j] = border_type;
        }
        else
        {
          _border_coord[i * _size_y + j]->x() = (_border_coord[i * _size_y + j]->x() + (*seg)->pointA()->x()) / 2;
          _border_coord[i * _size_y + j]->y() = (_border_coord[i * _size_y + j]->y() + (*seg)->pointA()->y()) / 2;
          _z[i * _size_y + j] = (_z[i * _size_y + j] + (*seg)->pointA()->z()) / 2;
          _dist[i * _size_y + j] = (_dist[i * _size_y + j] + (*seg)->pointA()->distance()) / 2;
        }
      }
      else
        // this is the first assignment for this cell
      {
        _border_coord[i * _size_y + j] = new Vec2d((*seg)->pointA()->x(), (*seg)->pointA()->y());
        _z[i * _size_y + j] = (*seg)->pointA()->z();
        _dist[i * _size_y + j] = (*seg)->pointA()->distance();
        _type[i * _size_y + j] = border_type;
      }
    }
  }

  // Third iteration on points:
  // - detect other border points with segment box overlapping tests

  unsigned k, l, x, y;
  double val1, val2, val3, val4, alpha;
  dx = (_bb_max[0] - _bb_min[0]) / _size_x;
  dy = (_bb_max[1] - _bb_min[1]) / _size_y;
  Vec2d halfsize(dx/2, dy/2);
  Vec2d center;
  Vec2d coord;

  st = chain->strokes()->begin();
  st_end = chain->strokes()->end();
  for (; st != st_end; ++st) {
    seg = (*st)->segments()->begin();
    seg_end = (*st)->segments()->end();
    for (; seg != seg_end; ++seg) {
      i = static_cast<unsigned>((((*seg)->pointA()->x() - _bb_min[0]) / (_bb_max[0] - _bb_min[0])) * _size_x);
      j = static_cast<unsigned>((((*seg)->pointA()->y() - _bb_min[1]) / (_bb_max[1] - _bb_min[1])) * _size_y);
      k = static_cast<unsigned>((((*seg)->pointB()->x() - _bb_min[0]) / (_bb_max[0] - _bb_min[0])) * _size_x);
      l = static_cast<unsigned>((((*seg)->pointB()->y() - _bb_min[1]) / (_bb_max[1] - _bb_min[1])) * _size_y);
      if (i && i >= _size_x)
        i = _size_x - 1;
      if (j && j >= _size_y)
        j = _size_y - 1;
      if (k && k >= _size_x)
        k = _size_x - 1;
      if (l && l >= _size_y)
        l = _size_y - 1;
      if (i > k)
        utils::swap(i, k);
      if (j > l)
        utils::swap(j, l);
      for (x = i; x <= k; ++x) {
        center[0] = _bb_min[0] + halfsize[0] + x * dx;
        for (y = j; y <= l; ++y) {
   
         // JDW Set cells to be clamped (second pass)
      if( (*st)->type() == Stroke::OUT)
        _clamp[x * _size_y + y] = true;
   
          if (_type[x * _size_y + y] != UNKNOWN)
            continue;
          center[1] = _bb_min[1] + halfsize[1] + y * dy;
          pts_nb = 0;
          if (utils::collisionSegmentX(**seg, center[0], coord[1]) &&
              coord[1] >= center[1] - halfsize[1] && coord[1] <= center[1] + halfsize[1])
            pts_nb += 1;
          if (utils::collisionSegmentY(**seg, center[1], coord[0]) &&
              coord[0] >= center[0] - halfsize[0] && coord[0] <= center[0] + halfsize[0])
            pts_nb += 2;
          if (!pts_nb)
            continue;
          switch (pts_nb) {
          case 3:
            _border_coord[x * _size_y + y] = new Vec2d((center[0] + coord[0]) / 2,
        					       (center[1] + coord[1]) / 2);
            break;
          case 2:
            _border_coord[x * _size_y + y] = new Vec2d(coord[0], center[1]);
            break;
          case 1:
            _border_coord[x * _size_y + y] = new Vec2d(center[0], coord[1]);
            break;
          default:
            ;
          }
          val1 = _border_coord[x * _size_y + y]->x() - (*seg)->pointA()->x();
          val2 = _border_coord[x * _size_y + y]->y() - (*seg)->pointA()->y();
          alpha = sqrt(val1*val1 + val2*val2) / (*seg)->length();
          _z[x * _size_y + y] = (1 - alpha) * (*seg)->pointA()->z() + alpha * (*seg)->pointB()->z();
          _dist[x * _size_y + y] = (1 - alpha) * (*seg)->pointA()->distance() + alpha * (*seg)->pointB()->distance();
          _type[x * _size_y + y] = BORDER;
        }
      }
    }
  }

  // Second iteration (recursive) on the arrays:
  // - determine IN and OUT
  // JDW Propagates type OUT from the 4 outer corners of the garment map
  // by examining the immediate 4-neighbours. If none of these are in the
  // corner or already labelled OUT and current type is not UNKNOWN then
  // it must be IN.

  k = _size_x - 1;
  l = _size_y - 1;
  for (j = 0; j <= l/2; ++j)
    for (i = 0; i <= k/2; ++i) {
      for (x = i; x <= k-i; ++x) {
        propagateType(x, j);
        propagateType(k-x, l-j);
      }
      for (y = j+1; y <= l-j-1; ++y) {
        propagateType(i, y);
        propagateType(k-i, l-y);
      }
    }

  // Third iteration on the arrays:
  // - diffuse the distance values

  double *array_tmp = new double[_size_x * _size_y];
  double dmax;
  double dmax_tmp;
  if (!iter_nb)
    iter_nb = 10000;
  for (; iter_nb > 0; --iter_nb) {
    dmax = 0;
    for (i = 0; i < _size_x; ++i)
      for (j = 0; j < _size_y; ++j) {
        array_tmp[i * _size_y + j] = applyMask(i, j, _dist, diffusion_mask);
        dmax_tmp = utils::abs(array_tmp[i * _size_y + j] - _dist[i * _size_y + j]);
        if (dmax_tmp > dmax)
          dmax = dmax_tmp;
      }
    for (i = 0; i < _size_x; ++i)
      for (j = 0; j < _size_y; ++j)
        _dist[i * _size_y + j] = array_tmp[i * _size_y + j];
    if (stop_threshold && dmax <= stop_threshold)
      break;
  }

  // Fourth iteration to compute the z buffer.
  // And marks points inside the garment but outside the body (e.g. between legs)
  double zMin = numeric_limits<double>::max();
  double zMax = -numeric_limits<double>::max();
  for (i = 0; i < _size_x; ++i) {
    center[0] = _bb_min[0] + halfsize[0] + i * dx;
    for (j = 0; j < _size_y; ++j) {
      if  (_type[i * _size_y + j] == IN) {
           center[1] = _bb_min[1] + halfsize[1] + j * dy;
           point_utils::zFromDist(center[0], center[1], _dist[i * _size_y + j], _z[i * _size_y + j], front);
 	   if (point_utils::type(center[0], center[1], _dist[i * _size_y + j]) != Point::IN)
       _type[i * _size_y + j] = IN_OUT; // JDW Inside garment but outside body
     
     // JDW stats for debug
     zMin = utils::min(zMin, _z[i * _size_y + j]);
     if(zMax < _z[i * _size_y + j]) zMax = _z[i * _size_y + j];
     if(_z[i * _size_y + j] > 30) {
       std::cout << "Problem z: " << _z[i * _size_y + j] << ", dist: " << _dist[i * _size_y + j] << std::endl;
     }
      }
    }
  }
  
  std::cout << "Garment zMin: " << zMin << std::endl;
  std::cout << "Garment zMax: " << zMax << std::endl;
  
  // Fifth iteration to fill regions that are between 2 IN segments horizontally.
 
  for (j = 0; j < _size_y; ++j) {
    center[1] = _bb_min[1] + halfsize[1] + j * dy;
    i = 0;
    while (i < _size_x) {
      // Find next IN_OUT cell in this row and store location in k
      while (_type[i * _size_y + j] != IN_OUT) {
        if (i >= _size_x)
          break;
        ++i;
      }
      k = i;
      // Find end of this IN_OUT section and store location in l
      while (_type[i * _size_y + j] == IN_OUT) {
        if (i >= _size_x)
          break;
        ++i;
      }
      l = i;
      // Check that this stretch of IN_OUT is fully inside garment
      if (_type[(k-1) * _size_y + j] != IN || _type[l * _size_y + j] != IN)
        continue;
      
      // use a bezier curve to interpolate the depth values horizontally
      y = l - k;
      val1 = _z[(k-1) * _size_y + j];
      val2 = _z[(k-1) * _size_y + j] + (_z[(k-1) * _size_y + j] - _z[(k-2) * _size_y + j]) * (y+1) / 3.0;
      val3 = _z[l * _size_y + j] + (_z[l * _size_y + j] - _z[(l+1) * _size_y + j]) * (y+1) / 3.0;
      val4 = _z[l * _size_y + j];
      
      // fill in the values
      for (x = 0; x < y; ++x) {
        alpha = static_cast<double>(x+1) / (y+1);
        utils::bezier3(_z[(k+x) * _size_y + j], val1, val2, val3, val4, alpha); // FIXME
        center[0] = _bb_min[0] + halfsize[0] + (k+x) * dx;
        point_utils::distFromZ(center[0], center[1], _z[(k+x) * _size_y + j], _dist[(k+x) * _size_y + j]);

 
        _type[(k+x) * _size_y + j] = IN;
      }
    }
  }
  
  

  // Sixth iteration on the arrays:
  // - diffuse the distance values in the IN_OUT regions
  // TODO JDW why?
  
   // JDW OFF

  if (!iter_nb)
    iter_nb = 10000;
  for (; iter_nb > 0; --iter_nb) {
    dmax = 0;
    for (i = 0; i < _size_x; ++i)
      for (j = 0; j < _size_y; ++j) {
        array_tmp[i * _size_y + j] = applyMask(i, j, _dist, diffusion_mask, IN_OUT);
        dmax_tmp = utils::abs(array_tmp[i * _size_y + j] - _dist[i * _size_y + j]);
        if (dmax_tmp > dmax)
          dmax = dmax_tmp;
      }
    for (i = 0; i < _size_x; ++i)
      for (j = 0; j < _size_y; ++j)
        _dist[i * _size_y + j] = array_tmp[i * _size_y + j];
    if (stop_threshold && dmax <= stop_threshold)
      break;
  }
  
  
  // JDW END OFF
  // Seventh iteration to compute the z buffer in the IN_OUT regions.

  // JDW Debug
  zMax = -numeric_limits<double>::max();
  zMin = numeric_limits<double>::max();
  
  for (i = 0; i < _size_x; ++i) {
    center[0] = _bb_min[0] + halfsize[0] + i * dx;
    for (j = 0; j < _size_y; ++j) {
      center[1] = _bb_min[1] + halfsize[1] + j * dy;
      if (_type[i * _size_y + j] == IN_OUT) {
        point_utils::zFromDist(center[0], center[1], _dist[i * _size_y + j], _z[i * _size_y + j], front);
        _type[i * _size_y + j] = IN;
 // JDW
 zMin = utils::min(_z[i * _size_y + j],zMin);
 zMax = utils::max(_z[i * _size_y + j],zMax);
      }
    }
  }
  
  std::cout << "7th It Garment zMin: " << zMin << std::endl;
  std::cout << "7th It Garment zMax: " << zMax << std::endl;
  

  
  
  // Last iterations to smooth the result and set _z_min and _z_max.
  
   // TODO JDW turned off to test problem with zdepth in normalised fields
  
  // smooth the zdepth surface, but interpolate between the original and smoothed using
  // an alpha factor which increases as the distance from the surface increases.
  // closer to the surface use the actual value more
  // further from the surface use the smoothed value more.
  
  iter_nb = 100; // FIXME
  for (; iter_nb > 0; --iter_nb) {
    for (i = 0; i < _size_x; ++i) {
      for (j = 0; j < _size_y; ++j) {
            alpha = _dist[i * _size_y + j] / RepositoryPretreatmentCanvas::distanceField()->actualMaxDist();
       // alpha = 0.95; // JDW set constant for now while testing new fields.
             array_tmp[i * _size_y + j] = alpha * applyMask(i, j, _z, smoothing_mask) + (1 - alpha) * _z[i * _size_y + j];
      } //std::cout << "Alpha " << alpha << std::endl;
    }
    for (i = 0; i < _size_x; ++i)
      for (j = 0; j < _size_y; ++j)
        _z[i * _size_y + j] = array_tmp[i * _size_y + j];
  }
  
  
  delete[] array_tmp;
  
  

  _min_z = std::numeric_limits<double>::max();
  _max_z = -std::numeric_limits<double>::max();
  for (i = 0; i < _size_x; ++i) {
    for (j = 0; j < _size_y; ++j) {
      if (_type[i * _size_y + j] != OUT) {
        if (_z[i * _size_y + j] < _min_z)
          _min_z = _z[i * _size_y + j];
        if (_z[i * _size_y + j] > _max_z)
          _max_z = _z[i * _size_y + j];
      }
    }
  }
  std::cout << "FinalMaxZ: " << _max_z << std::endl;
  std::cout << "FinalMinZ: " << _min_z << std::endl;

}

// Just renders the grid
void GarmentMap::renderGL()
{
  float dx = bbSize()[0] / sizeX();
  float dy = bbSize()[1] / sizeY();
  
  glPushAttrib(GL_POLYGON_BIT | GL_LIGHTING_BIT);
  
  glColor3f(1.0,1.0,1.0);
  glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
  
  
  for(int j=0;j<sizeY();j++)
  {
    glBegin(GL_QUAD_STRIP);
    for(int i=0;i<sizeX()+1;i++)
    {
      
        glVertex3d( _bb_min[0] + (i * dx), _bb_min[1] + (j * dy) + dy, 0);
        glVertex3d( _bb_min[0] + (i * dx), _bb_min[1] + (j * dy) , 0);

    }
    glEnd();
  } // end loop over increasing y
  
  glPopAttrib();
  
  if (pointGrid) {
    glBegin(GL_POINTS);
    for (unsigned i = 0; i < pointGridSize; ++i)
    {
      if(pointGrid[i]) 
      {
        vecmat::Vector<double,5>* p = pointGrid[i];
        glVertex3d((*p)[0], (*p)[1], 0.0);
      }
    }
    glEnd();
  }
  
}

// Alter the garment z-depth using a 2d gaussian along a stroke
void GarmentMap::gfold(Stroke *stroke, bool front)
{
  // sense means whether the fold is concave or convex. 
  // This needs to depends on the gestures supplied but currently
  // both front and back garment folds are generated concave.
  std::cout << "Applying fold alterations" << std::endl;
  stroke->recalculateBB();
  const double gmCellWidth = this->bbSize()[0] / this->sizeX();
  const double gmCellHeight = this->bbSize()[1] / this->sizeY();
  
  /*
  std::cout << "GarmentMap Details:" << std::endl;
  std::cout << "BBMin: " << bbMin() << std::endl;
  std::cout << "BBMax: " << bbMax() << std::endl;
  std::cout << "CellWidth: " << gmCellWidth << std::endl;
  std::cout << "CellHeight: " << gmCellHeight << std::endl;
  std::cout << "SizeX: " << sizeX() << std::endl;
  std::cout << "SizeY: " << sizeX() << std::endl;
  */
  

  
  // The fold doesn't exist outside the garment - so clip it
  Vec2d clippedStrokeBBMin = Vec2d( stroke->bbMin()[0], stroke->bbMin()[1] );
  clippedStrokeBBMin.max( this->bbMin() );
  Vec2d clippedStrokeBBMax = Vec2d( stroke->bbMax()[0], stroke->bbMax()[1] );
  clippedStrokeBBMax.min( this->bbMax() );
  
  // Define a new grid to use in combining the gaussian mask contributions to the offset
  double xOffset = abs(this->bbMin()[0] - clippedStrokeBBMin[0]);
  double yOffset = abs(this->bbMin()[1] - clippedStrokeBBMin[1]);
  
  
  std::cout << "XOff: " << xOffset << std::endl;
  std::cout << "YOff: " << yOffset << std::endl;
  
  Vec2d strokeBBSize( clippedStrokeBBMax[0] - clippedStrokeBBMin[0], 
                      clippedStrokeBBMax[1] - clippedStrokeBBMin[1] );
  int strokeBBNbCellsWide = static_cast<int>( floor(strokeBBSize[0] / gmCellWidth) + 1 ); // add one so we have integer width definitely wider than the original BB
  int strokeBBNbCellsHigh = static_cast<int>( floor(strokeBBSize[1] / gmCellHeight) + 1 ); 
  
  std::cout << "Initial tmpGrid width: " << strokeBBNbCellsWide << std::endl;
  
  // align grid to GarmentMap (move borders to nearest cell borders)
  double xRem = fmod(xOffset, gmCellWidth);
  double yRem = fmod(yOffset, gmCellHeight);
  
  std::cout << "XREM: " << xRem << std::endl;
  std::cout << "YREM: " << yRem << std::endl;
  
  // Now add another cell to right hand and upper boundaries to cover any shortfall due to alignment.
  strokeBBNbCellsWide+=1;
  strokeBBNbCellsHigh+=1;
  
  
  Vec2d tmpGridBBMin = Vec2d( clippedStrokeBBMin[0] - xRem, clippedStrokeBBMin[1] - yRem );
  std::cout << "Initial BBMin: " << tmpGridBBMin << std::endl;
  
  
  double startGaussianHeight = stroke->getAParams()[0];
  double endGaussianHeight   = stroke->getBParams()[0];
  
  
  double startGaussianRadius  = stroke->getAParams()[1];
  double endGaussianRadius    = stroke->getBParams()[1];
  
  double maxRadius = utils::max(startGaussianRadius,endGaussianRadius);
  
  // expand grid to allow for edge effects of mask

  
  int expansionX = static_cast<int>((maxRadius) / gmCellWidth) +1; // make larger to account for rounding down
  //if(!(expansionX%2)) expansionX += 1; // make odd
  int expansionY = static_cast<int>((maxRadius) / gmCellHeight) +1;
  //if(!(expansionY%2)) expansionY += 1; // make odd
  
  tmpGridBBMin[0] -= (expansionX * gmCellWidth);
  tmpGridBBMin[1] -= (expansionY * gmCellHeight);
  
  strokeBBNbCellsWide += (expansionX  * 2);
  strokeBBNbCellsHigh += (expansionY * 2);
  
  
  
  std::cout << "BBMin after edge effect: " << tmpGridBBMin << std::endl;
  std::cout << "New tmpGrid width: " << strokeBBNbCellsWide << std::endl;
  
  // TODO boundary checks against garment map - don't want to index outside of it
  
  Vec2d tmpGridBBMax = Vec2d( tmpGridBBMin[0] + (strokeBBNbCellsWide*gmCellWidth) 
                             ,tmpGridBBMin[1] + (strokeBBNbCellsHigh*gmCellHeight) );
  
  Vec2d tmpGridBBSize = tmpGridBBMax - tmpGridBBMin;
  
  int gridOffsetX = abs(bbMin()[0] - tmpGridBBMin[0]) / gmCellWidth;
  int gridOffsetY = abs(bbMin()[1] - tmpGridBBMin[1]) / gmCellHeight;
  if(tmpGridBBMin[0] < bbMin()[0]) gridOffsetX *= -1;
  if(tmpGridBBMin[1] < bbMin()[1]) gridOffsetY *= -1;
  
  std::cout << "tmpGrid Details:" << std::endl;
  std::cout << "BBMin: " << tmpGridBBMin << std::endl;
  std::cout << "BBMax: " << tmpGridBBMax << std::endl;
  std::cout << "CellWidth: " << gmCellWidth << std::endl;
  std::cout << "CellHeight: " << gmCellHeight << std::endl;
  std::cout << "SizeX: " << strokeBBNbCellsWide << std::endl;
  std::cout << "SizeY: " << strokeBBNbCellsHigh << std::endl;
  
  std::cout << "gridOffsetX: " << gridOffsetX << std::endl;
  std::cout << "gridOffsetY: " << gridOffsetY << std::endl;
  
  
  double *tmpGrid = new double[ strokeBBNbCellsWide * strokeBBNbCellsHigh ];
  // initialise grid to 0.0;
  for(int i =0; i<strokeBBNbCellsWide;i++)
    for(int j =0; j<strokeBBNbCellsHigh;j++)
      tmpGrid[ i*strokeBBNbCellsHigh + j ] = 0.0;
  
  // Sample gaussian correctly. 
  //  For each map cell overlapped by fold strokes:
  //   Determine position at which gaussian centre should be placed (average when >1 seg in cell)
  //   Sample gaussian at centre of this and all surrounding cells within a significant radius
 
  
  if (pointGrid) {
    for (unsigned i = 0; i < pointGridSize; ++i)
        if (pointGrid[i])
          delete pointGrid[i];
    delete[] pointGrid;
  }
  
  pointGrid = new vecmat::Vector<double, 5>* [ strokeBBNbCellsWide * strokeBBNbCellsHigh ];
  pointGridSize = strokeBBNbCellsWide * strokeBBNbCellsHigh;
  // initialise grid to null pointers;
  for(int i =0; i<strokeBBNbCellsWide;i++)
    for(int j =0; j<strokeBBNbCellsHigh;j++)
      pointGrid[ i*strokeBBNbCellsHigh + j ] = static_cast<vecmat::Vector<double,5>*>(0);
  
  
  std::list<Segment *>::const_iterator sgit, sgit_end;
  sgit = stroke->segments()->begin();
  sgit_end = stroke->segments()->end();
  int segTotal = stroke->segments()->size();
  //std::cout << "Nb of Segments: " << segTotal << std::endl;
  int scount=1;
  
  stroke->updateLength();
  //std::cout << "Stroke length: " << stroke->length() << std::endl;
  double currentLength = 0.0;
  
  for(;sgit!=sgit_end;sgit++)
  { 
    //std:cout << "Processing seg: " << scount++ << "/" << segTotal << std::endl;
    Segment *seg = (*sgit);
    seg->updateLength();
//    std::cout << " Seg length: " << seg->length() << std::endl;
    currentLength += seg->length();
    double percentageLength = currentLength / stroke->length();

    // ignore points outside of clipped boundary
    /*
    if(seg->pointA()->x() < clippedStrokeBBMin[0] 
       || seg->pointA()->x() > clippedStrokeBBMax[0] 
       || seg->pointA()->y() < clippedStrokeBBMin[1] 
       || seg->pointA()->y() > clippedStrokeBBMax[1] 
      ) continue;
    */
    Segment clippedSeg = seg->clippedToRect( Vec2d(clippedStrokeBBMin[0], clippedStrokeBBMax[1]),
                                             Vec2d(clippedStrokeBBMax[0], clippedStrokeBBMin[1])
                                           );
    if(clippedSeg.type() == Segment::OUT) continue;
    //std::cout << "Fold stroke overlaps garment boundary" << std::endl;
    
    
    double percentageOfXDist = ((seg->pointA()->x() - tmpGridBBMin[0]) / tmpGridBBSize[0]);
    double percentageOfYDist = ((seg->pointA()->y() - tmpGridBBMin[1]) / tmpGridBBSize[1]);
    //std::cout << "PX: " << percentageOfXDist << std::endl;
    //std::cout << "PY: " << percentageOfYDist << std::endl;
    int AXindex = static_cast<int> (percentageOfXDist * (strokeBBNbCellsWide));
    if(AXindex >= strokeBBNbCellsWide) AXindex = strokeBBNbCellsWide -1; // rare case, shouldn't exceed index range 
    int AYindex = static_cast<int> (percentageOfYDist * (strokeBBNbCellsHigh));
    if(AYindex >= strokeBBNbCellsHigh) AYindex = strokeBBNbCellsHigh -1;// rare case 
    
    percentageOfXDist = ((seg->pointB()->x() - tmpGridBBMin[0]) / tmpGridBBSize[0]);
    percentageOfYDist = ((seg->pointB()->y() - tmpGridBBMin[1]) / tmpGridBBSize[1]);
    
    int BXindex = static_cast<int> (percentageOfXDist * (strokeBBNbCellsWide));
    if(BXindex >= strokeBBNbCellsWide) BXindex = strokeBBNbCellsWide -1; // rare case, shouldn't exceed index range 
    int BYindex = static_cast<int> (percentageOfYDist * (strokeBBNbCellsHigh));
    if(BYindex >= strokeBBNbCellsHigh) BYindex = strokeBBNbCellsHigh -1;// rare case 
    
    // We have range of cells potentially covered by this straight line segment
    // For each cell, clip and determine centre point - then store. If a measurement
    // already exists, then average together.
    
    if(AXindex > BXindex) std::swap( AXindex, BXindex );
    if(AYindex > BYindex) std::swap( AYindex, BYindex );
    
    for(int i=AXindex;i<=BXindex;i++)
    {
      for(int j=AYindex;j<=BYindex;j++)
      {
        //std::cout << "Range loop [" << i << "/" << BXindex  
        //         << ", " << j << "/" << BYindex << "]" << std::endl;
        Vec2d tl, br;
        tl[0] = tmpGridBBMin[0] + i * gmCellWidth;
        tl[1] = tmpGridBBMin[1] + j * gmCellHeight + gmCellHeight;
        
        br[0] = tl[0] + gmCellWidth;
        br[1] = tl[1] - gmCellHeight;
        
        Segment cSeg = clippedSeg.clippedToRect(tl,br);
        if(cSeg.type() == Segment::OUT) {
          //std::cout << "segment doesn't overlap fold stroke" << std::endl;
          continue; // no intersection
        }
        //std::cout << "Overlap found" << std::endl;
        
        int indexPointGrid1D = i*strokeBBNbCellsHigh + j;
        
        vecmat::Vector<double,5>* averageOfThisSeg = new vecmat::Vector<double,5>;
        (*averageOfThisSeg)[0] = (cSeg.pointA()->x() + cSeg.pointB()->x()) / 2.;
        (*averageOfThisSeg)[1] = (cSeg.pointA()->y() + cSeg.pointB()->y()) / 2.;
        (*averageOfThisSeg)[2] = 1.0; // count of points in this cell
        //(*averageOfThisSeg)[3] = 2.0; // radius of gaussian
        
        
        double radius = startGaussianRadius * (1.0 - percentageLength) + 
                       endGaussianRadius * percentageLength; 
    //    std::cout << "Percentlength: " << percentageLength << std::endl;
        (*averageOfThisSeg)[3] = radius;
        
        double gaussheight = startGaussianHeight * (1.0 - percentageLength) +
              endGaussianHeight * percentageLength;
        
        (*averageOfThisSeg)[4] = gaussheight; // height of gaussian

        
        vecmat::Vector<double,5> *p = pointGrid[indexPointGrid1D];
        if(p)
        {
          // combine with existing point, add count of points to z co-ord for averaging
          (*p)[0] = ( (*averageOfThisSeg)[0] + (*p)[0] ); // x
          (*p)[1] = ( (*averageOfThisSeg)[1] + (*p)[1] ); // y
          (*p)[2] += 1.0; // pb points contributing to this average
          (*p)[3] = ( (*averageOfThisSeg)[3] + (*p)[3] ); // radius of gaussian at this point
          (*p)[4] = ( (*averageOfThisSeg)[4] + (*p)[4] ); // height of gaussian at this point
  //        std::cout << "Made average point" << std::endl;
          delete averageOfThisSeg;
        }
        else
        {
          pointGrid[indexPointGrid1D] = averageOfThisSeg;
//          std::cout << "Made new point" << std::endl;
        }
         
      }
    } // end loop over affected area of pointGrid[]
  } // end loop over all segments
 // std::cout << "End seg loop" << std::endl;
  
  // We have now got a grid of sample positions with which to sample the gaussian.
  // loop over the grid and each time a sample point exists determine the sample value for 
  // this cell and cells within an affected radius.
     
  double maxInGaussGrid = 0.0;
        
  for(int i=0;i<strokeBBNbCellsWide;i++)
    for(int j=0;j<strokeBBNbCellsHigh;j++)
  {
    vecmat::Vector<double,5> *p; // check if we sample this cell
    if((p = pointGrid[i * strokeBBNbCellsHigh + j]))
    {
      // now average points sampled for multiple cell samples.
      if((*p)[2] > 1.0) {
        (*p)[0] /= (*p)[2];
        (*p)[1] /= (*p)[2];
        
        (*p)[3] /= (*p)[2];
        (*p)[4] /= (*p)[2];
        (*p)[2] = 1.0; // reset
      }
      // loop over affected cells
      
      // We want to sample ~99 of the area under the gaussian. 
      // and ensure gaussian falls to near zero at edges of sample area.
      // So we want the radius to be 3 sigma.
      double radius = (*p)[3];
      double dHalfWidth = radius;
      int iHalfWidth = static_cast<int>(dHalfWidth / gmCellWidth)+1;
      int iHalfHeight = static_cast<int>(dHalfWidth / gmCellHeight)+1;
      
      
      
      //for(int ii=-hwidth;ii<=hwidth;ii++) 
      for(int ii=-iHalfWidth;ii<=iHalfWidth;ii++) 
      {
        if( 0 > (i+ii) || (i+ii) > strokeBBNbCellsWide-1) continue;
        //for(int jj=-hheight;jj<=hheight;jj++) 
        for(int jj=-iHalfHeight;jj<=iHalfHeight;jj++) 
        {
          if( 0 > (j+jj) || (j+jj) > strokeBBNbCellsHigh-1) continue;
          Vec2d cellCentre;
          cellCentre.x() = tmpGridBBMin[0] + (i+ii+0.5) * gmCellWidth;
          cellCentre.y() = tmpGridBBMin[1] + (j+jj+0.5) * gmCellHeight;
          
          double ddx = cellCentre.x() - (*p)[0];
          double ddy = cellCentre.y() - (*p)[1];
          int indextmpGrid1D = (ii+i)*strokeBBNbCellsHigh + (jj+j);
          
          // modify ddx and ddy so we sample the gaussian over the required range.
          // roughly one cell width = 2.0 in gaussian space
          
          //ddx; // /= gmCellWidth;
          //ddx *= dx;
          
          //ddy; // /= gmCellHeight;
          //ddy *= dy;
          
          //std::cout << "At [" << ii+i << ", " << jj+j << "]" << std::endl;
          //std::cout << "DX: " << ddx << ", DY: " << ddy << std::endl;
          //std::cout << "Gauss: " << utils::gauss2d(1.0,2.0,ddx,ddy) << std::endl << std::endl;
          double gaussianAmplitude = (*p)[4];
          double guassianSigma = (*p)[3] / 3.0; // 99.7% of the power within 3sigma
          
          //tmpGrid[ indextmpGrid1D  ] = 
              //max( utils::gauss2d( gaussianAmplitude,guassianSigma,ddx,ddy) , tmpGrid[ indextmpGrid1D ]);
          
          tmpGrid[ indextmpGrid1D  ] += 
              utils::gauss2d( gaussianAmplitude,guassianSigma,ddx,ddy);
          
          maxInGaussGrid = utils::max(abs(tmpGrid[ indextmpGrid1D  ]), abs(maxInGaussGrid));
          
          //std::cout << "Radius: " << (*p)[3] << std::endl;
        }
      }
    }
  }

  
  // Copy back to z buffer
  double zRange = _max_z - _min_z;
  const double percentOffset = 0.05; // offset some percentage of maximum z range
  
  for(int i =0; i<strokeBBNbCellsWide;i++)
  {
    for(int j =0; j<strokeBBNbCellsHigh;j++)
    {
      if(   i+gridOffsetX < 0 || i+gridOffsetX >= _size_x
         || j+gridOffsetY < 0 || j+gridOffsetY >= _size_y
        ) continue;
      
      int zindex = (i+gridOffsetX) * _size_y + j + gridOffsetY;
      if(zindex < 0 || zindex >= (_size_y*_size_x) )
      {
        // The fold doesn't exist outside the garment so don't try to write there
        std::cerr << "Warning: indexing outside _z on copy back fold offset: " << zindex << std::endl;
        continue;
      }
      
      /*
      if(_type[zindex] == BORDER || _type[zindex] == BORDER_EXT) {
        std::out << "Clamping borders during smoothing pass" << std::endl;
      } 
      */
      
      if(_clamp[zindex]) // don't alter - this is a seam that should co-incide with the back garment
        continue;
      
      // Scale back into 0-1 range;
      tmpGrid[ i * strokeBBNbCellsHigh + j ] /= maxInGaussGrid;
      // scale back to user specified range 
      tmpGrid[ i * strokeBBNbCellsHigh + j ] *= utils::max( abs( startGaussianHeight ), abs(endGaussianHeight) );
      
      if(front) {
        
        _z[zindex] += 
            tmpGrid[ i * strokeBBNbCellsHigh + j]; 
      } 
      else {
          _z[zindex] -= 
              tmpGrid[ i * strokeBBNbCellsHigh + j]; 
      } 

          //(percentOffset * zRange) * tmpGrid[ i * strokeBBNbCellsHigh + j]; // offset some percentage of maximum z range
      
    }
  }
  
  // smooth altered region 
  // reuse tmpGrid for result of smoothing
  
  Mat33d smask;
  smask(0, 0) = 1;
  smask(0, 1) = 1;
  smask(0, 2) = 1;
  smask(1, 0) = 1;
  smask(1, 1) = 1;
  smask(1, 2) = 1;
  smask(2, 0) = 1;
  smask(2, 1) = 1;
  smask(2, 2) = 1;
  
  for(int i =0; i<strokeBBNbCellsWide;i++)
  {
    for(int j =0; j<strokeBBNbCellsHigh;j++)
    {
      double sum = 0.0;
      double denom = 0.0;
      
      for(int mj = 0; mj < 3; mj++)
      {
        for(int mi = 0; mi < 3; mi++)
        {
          
          if(
             i+gridOffsetX+(mi-1) < 0 || i+gridOffsetX+(mi-1) >= _size_x
             || j + (mj-1) + gridOffsetY < 0 || j + (mj-1) + gridOffsetY >= _size_y
            ) continue;
          
          int zindex = (i+gridOffsetX+(mi-1)) * _size_y + j + (mj-1) + gridOffsetY;
          if(zindex < 0 || zindex >= (_size_y*_size_x) )
          {
            // fold doesn't exist outside of garment so don't write there
            std::cerr << "Warning, Smoothing Error: indexing outside _z: " << zindex << std::endl;
            continue;
          }
          sum += smask(mi, mj) * _z[zindex];
          denom += smask(mi, mj);
        }
      }
      if(denom != 0)
        tmpGrid[ i * strokeBBNbCellsHigh + j] = sum / denom;
      else 
        tmpGrid[ i * strokeBBNbCellsHigh + j] = _z[ (i+gridOffsetX) * _size_y + j +  gridOffsetY ];
      
    }
  }

  // copy back smoothed result
  for(int i =0; i<strokeBBNbCellsWide;i++)
  {
    for(int j =0; j<strokeBBNbCellsHigh;j++)
    {
      if(   i+gridOffsetX < 0 || i+gridOffsetX >= _size_x
            || j+gridOffsetY < 0 || j+gridOffsetY >= _size_y
        ) continue;
      
      int zindex = (i+gridOffsetX) * _size_y + j + gridOffsetY;
      if(zindex < 0 || zindex >= (_size_y*_size_x) )
      {
        // don't write outside the garment
        std::cerr << "Warning: indexing outside _z on z buffer write back, ignoring: " << zindex << std::endl;
        continue;
      }
      _z[zindex] = tmpGrid[ i * strokeBBNbCellsHigh + j];
    }
  }
  
  
  delete[] tmpGrid;
  /*
  if (pointGrid) {
    for (unsigned i = 0; i < strokeBBNbCellsWide; ++i)
      for (unsigned j = 0; j < strokeBBNbCellsHigh; ++j)
        if (pointGrid[i * strokeBBNbCellsHigh + j])
          delete pointGrid[i * strokeBBNbCellsHigh + j];
    delete[] pointGrid;
  }
  */
  
}


