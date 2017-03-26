//
//  Filename         : grid.cpp
//  Author           : Emmanuel Turquin
//  Purpose          : Class representing a grid of cells containing segments.
//  Date of creation : 05/04/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include "utils.h"
#include "segment.h"
#include "grid.h"

void Cell::addSegment(Segment *seg)
{
  if (!seg)
    return;
  _segments.push_back(seg);
}

void Cell::removeSegment(Segment *seg)
{
  if (!seg)
    return;
  _segments.remove(seg);
}

void Cell::clear()
{
  _segments.clear();
}

Cell::SegmentList* Cell::getSegments()
{
  return &_segments;
}

Grid::Grid(unsigned x_cells, unsigned y_cells,
	   double x_min, double y_min,
	   double x_max, double y_max)
{
  _x_cells = x_cells;
  _y_cells = y_cells;
  _x_min = x_min;
  _y_min = y_min;
  _x_max = x_max;
  _y_max = y_max;
  _cells = new Cell[x_cells * y_cells];
}

Grid::~Grid()
{
  delete[] _cells;
}

void Grid::clear()
{
  for (unsigned i = 0; i < _x_cells * _y_cells; ++i)
    _cells[i].clear();
}

void Grid::addSegment(Segment *seg)
{
  unsigned ai, aj, bi, bj;
  convertCoordinates(seg->pointA()->x(), seg->pointA()->y(), ai, aj);
  convertCoordinates(seg->pointB()->x(), seg->pointB()->y(), bi, bj);
  if (ai > bi)
    utils::swap(ai, bi);
  if (aj > bj)
    utils::swap(aj, bj);
  for (unsigned i = ai; i <= bi; ++i)
    for (unsigned j = aj; j <= bj; ++j) {
      _cells[i * _y_cells + j].addSegment(seg);
      seg->addCell(&_cells[i * _y_cells + j]);
    }
}

void Grid::removeSegment(Segment *seg)
{
  unsigned ai, aj, bi, bj;
  convertCoordinates(seg->pointA()->x(), seg->pointA()->y(), ai, aj);
  convertCoordinates(seg->pointB()->x(), seg->pointB()->y(), bi, bj);
  if (ai > bi)
    utils::swap(ai, bi);
  if (aj > bj)
    utils::swap(aj, bj);
  for (unsigned i = ai; i <= bi; ++i)
    for (unsigned j = aj; j <= bj; ++j) {
      _cells[i * _y_cells + j].removeSegment(seg);
      seg->removeCell(&_cells[i * _y_cells + j]);
    }
}

Cell::SegmentList Grid::getSegments(double x, double y)
{
  unsigned i, j;
  convertCoordinates(x, y, i, j);
  return *_cells[i * _y_cells + j].getSegments();
}

Cell::SegmentList Grid::getSegments(Segment *seg)
{
  Cell::SegmentList l;
  unsigned ai, aj, bi, bj;
  convertCoordinates(seg->pointA()->x(), seg->pointA()->y(), ai, aj);
  convertCoordinates(seg->pointB()->x(), seg->pointB()->y(), bi, bj);
  if (ai > bi)
    utils::swap(ai, bi);
  if (aj > bj)
    utils::swap(aj, bj);
  for (unsigned i = ai; i <= bi; ++i)
    for (unsigned j = aj; j <= bj; ++j)
      l.insert(l.end(), _cells[i * _y_cells + j].getSegments()->begin(), _cells[i * _y_cells + j].getSegments()->end());
  return l;
}
