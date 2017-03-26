//
//  Filename         : stroke.cpp
//  Author           : Emmanuel Turquin
//  Purpose          : Class representing a 2d stroke for the Canvas.
//  Date of creation : 04/26/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include <algorithm>
#include "utils.h"
#include "chain.h"
#include "stroke.h"
#include "vectypes.h"


void Stroke::setAParams(Vec2d params)
{
  _pt_a_params = params;
}
void Stroke::setBParams(Vec2d params)
{
  _pt_b_params = params;
}
Vec2d Stroke::getAParams()
{
  return _pt_a_params;
}
Vec2d Stroke::getBParams()
{
  return _pt_b_params;
}

void Stroke::addSegment(Segment *seg) {
  if (!seg)
    return;

  seg->setStroke(this);

  if (!_pt_a.get())
    _pt_a = seg->pointA();
  _pt_b = seg->pointB();
  _segments.push_back(seg);
}

bool Stroke::join(Stroke *stroke)
{
  if (!stroke || stroke == this)
    return false;

  if (!_pt_a.get() || !_pt_b.get()) {
    _pt_a = stroke->_pt_a;
    _pt_b = stroke->_pt_b;
    _segments = stroke->_segments;
  }
  else if (_pt_a.get() == stroke->_pt_a.get()) {
    stroke->reverse();
    _pt_a = stroke->_pt_a;
    _segments.splice(_segments.begin(), stroke->_segments);
  }
  else if (_pt_a.get() == stroke->_pt_b.get()) {
    _pt_a = stroke->_pt_a;
    _segments.splice(_segments.begin(), stroke->_segments);
  }
  else if (_pt_b.get() == stroke->_pt_a.get()) {
    _pt_b = stroke->_pt_b;
    _segments.splice(_segments.end(), stroke->_segments);
  }
  else if (_pt_b.get() == stroke->_pt_b.get()) {
    stroke->reverse();
    _pt_b = stroke->_pt_b;
    _segments.splice(_segments.end(), stroke->_segments);
  }
  else // These strokes shouldn't be joined.
    return false;

  SegmentList::iterator it = _segments.begin();
  SegmentList::iterator it_end = _segments.end();
  for (; it != it_end; ++it)
    (*it)->setStroke(this);

  chain()->strokes()->remove(stroke);
  delete(stroke);

  if (chain()->cycle()) {
    chain()->setPointA(chain()->strokes()->front()->pointA());
    chain()->setPointB(chain()->strokes()->front()->pointA());
  }    

  return true;
}

Stroke *Stroke::split(SegmentList::iterator seg_it)
{
  if (seg_it == _segments.begin())
    return 0;

  Stroke *new_stroke = new Stroke();
  new_stroke->setPointA((*seg_it)->pointA());
  new_stroke->setPointB(_pt_b);
  new_stroke->setChain(_chain);
  new_stroke->_segments.splice(new_stroke->_segments.end(), _segments, seg_it, _segments.end());
  SegmentList::iterator it = new_stroke->_segments.begin();
  SegmentList::iterator it_end = new_stroke->_segments.end();
  for (; it != it_end; ++it)
    (*it)->setStroke(new_stroke);

  Chain::StrokeList::iterator cur = find(_chain->strokes()->begin(), _chain->strokes()->end(), this);
  _chain->strokes()->insert(++cur, new_stroke);

  _pt_b = (*seg_it)->pointA();

  return new_stroke;
}

Stroke *Stroke::split(Segment *seg) {
  if (!seg)
    return 0;

  return split(std::find(_segments.begin(), _segments.end(), seg));
}

Stroke *Stroke::cleverSplit(SegmentList::iterator seg_it)
{
  if (seg_it == _segments.begin())
    return 0;

  SegmentList::iterator seg_it_next(seg_it);
  ++seg_it_next;

  if (seg_it_next == _segments.end())
    return split(seg_it);

  SegmentList::iterator seg_it_prev(seg_it);
  --seg_it_prev;

  double dotp1 = utils::normalizedDotProduct((*seg_it_prev)->pointA()->x() - (*seg_it_prev)->pointB()->x(),
					     (*seg_it_prev)->pointA()->y() - (*seg_it_prev)->pointB()->y(),
					     (*seg_it)->pointA()->x() - (*seg_it)->pointB()->x(),
					     (*seg_it)->pointA()->y() - (*seg_it)->pointB()->y());
  double dotp2 = utils::normalizedDotProduct((*seg_it)->pointA()->x() - (*seg_it)->pointB()->x(),
					     (*seg_it)->pointA()->y() - (*seg_it)->pointB()->y(),
					     (*seg_it_next)->pointA()->x() - (*seg_it_next)->pointB()->x(),
					     (*seg_it_next)->pointA()->y() - (*seg_it_next)->pointB()->y());
  if (dotp2 < dotp1)
    return split(seg_it_next);
  return split(seg_it);
}

Stroke *Stroke::cleverSplit(Segment *seg) {
  if (!seg)
    return 0;

  return cleverSplit(std::find(_segments.begin(), _segments.end(), seg));
}

void Stroke::recalculateBB() 
{
  Vec3d minBB(999999,999999,0), maxBB(-9999999,-999999,0); // store maximum extents
    

  Stroke::SegmentList::iterator seg;
  Stroke::SegmentList::iterator seg_end;
  
  seg = this->segments()->begin();
  seg_end = this->segments()->end();
    
  for (; seg != seg_end; ++seg) {
      // Update X extents if required
    if( (*seg)->pointA()->x() < minBB.x()) {
      minBB[0] = (*seg)->pointA()->x() ;
    }
    if( (*seg)->pointB()->x() < minBB.x()) {
      minBB[0] = (*seg)->pointB()->x();
    }
    if( (*seg)->pointA()->x() > maxBB.x()) {
      maxBB[0] = (*seg)->pointA()->x();
    }
    if( (*seg)->pointB()->x() > maxBB.x()) {
      maxBB[0] = (*seg)->pointB()->x() ;
    }
      
      // Update Y extents if required
    if( (*seg)->pointA()->y() < minBB.y()) {
      minBB[1] = (*seg)->pointA()->y() ;
    }
    if( (*seg)->pointB()->y() < minBB.y()) {
      minBB[1] = (*seg)->pointB()->y() ;
    }
    if( (*seg)->pointA()->y() > maxBB.y()) {
      maxBB[1] = (*seg)->pointA()->y();
    }
    if( (*seg)->pointB()->y() > maxBB.y()) {
      maxBB[1] = (*seg)->pointB()->y();
    }
  }

  _bbMin = minBB;
  _bbMax = maxBB;
}
