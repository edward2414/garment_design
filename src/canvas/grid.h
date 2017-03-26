//
//  Filename         : grid.h
//  Author           : Emmanuel Turquin
//  Purpose          : Class representing a grid of cells containing segments.
//  Date of creation : 05/04/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  _GRID_H
# define _GRID_H

# include <list>

class Segment;

class Cell
{
public:

  typedef std::list<Segment *>	SegmentList;

  Cell() {}
  ~Cell() {}

  void addSegment(Segment *seg);
  void removeSegment(Segment *seg);

  SegmentList* getSegments();

  void clear();

private:

  SegmentList	_segments;
};

class Grid
{
public:

  Grid(unsigned x_cells = 20, unsigned y_cells = 20,
       double x_min = -1, double y_min = -1,
       double x_max = 1, double y_max = 1);

  ~Grid();

  void addSegment(Segment *seg);
  void removeSegment(Segment *seg);

  Cell::SegmentList getSegments(double x, double y);
  Cell::SegmentList getSegments(Segment *seg);

  void clear();

private:

  void convertCoordinates(double x, double y, unsigned& i, unsigned& j)
  {
    i = static_cast<unsigned>(((x - _x_min) * (_x_cells - 1)) / (_x_max - _x_min));
    j = static_cast<unsigned>(((y - _y_min) * (_y_cells - 1)) / (_y_max - _y_min));
  }

  unsigned _x_cells;
  unsigned _y_cells;
  double _x_min;
  double _x_max;
  double _y_min;
  double _y_max;
  Cell*	_cells;
};

#endif // _GRID_H
