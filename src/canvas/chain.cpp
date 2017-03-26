//
//  Filename         : chain.cpp
//  Author           : Emmanuel Turquin
//  Purpose          : Class representing a chain of strokes for the Canvas.
//  Date of creation : 04/26/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include <algorithm>
#include "chain.h"
#include "vectypes.h"

void Chain::addStroke(Stroke *stroke) {
  if (!stroke)
    return;

  stroke->setChain(this);

  if (!_pt_a.get() || !_pt_b.get()) {
    _pt_a = stroke->pointA();
    _pt_b = stroke->pointB();
    _strokes.push_back(stroke);
  }
  else if (_pt_a.get() == stroke->pointA().get()) {
    stroke->reverse();
    _pt_a = stroke->pointA();
    _strokes.push_front(stroke);
  }
  else if (_pt_a.get() == stroke->pointB().get()) {
    _pt_a = stroke->pointA();
    _strokes.push_front(stroke);
  }

  else if (_pt_b.get() == stroke->pointA().get()) {
    _pt_b = stroke->pointB();
    _strokes.push_back(stroke);
  }
  else { // _pt_b.get() == stroke->pointB().get())
    stroke->reverse();
    _pt_b = stroke->pointB();
    _strokes.push_back(stroke);
  }
}

bool Chain::join(Chain *chain)
{
  if (!chain || chain == this)
    return false;

  if (!_pt_a.get() || !_pt_b.get()) {
    _pt_a = chain->_pt_a;
    _pt_b = chain->_pt_b;
    _strokes = chain->_strokes;
  }
  else if (_pt_a.get() == chain->_pt_a.get()) {
    chain->reverse();
    _pt_a = chain->_pt_a;
    _strokes.splice(_strokes.begin(), chain->_strokes);
  }
  else if (_pt_a.get() == chain->_pt_b.get()) {
    _pt_a = chain->_pt_a;
    _strokes.splice(_strokes.begin(), chain->_strokes);
  }
  else if (_pt_b.get() == chain->_pt_a.get()) {
    _pt_b = chain->_pt_b;
    _strokes.splice(_strokes.end(), chain->_strokes);
  }
  else if (_pt_b.get() == chain->_pt_b.get()) {
    chain->reverse();
    _pt_b = chain->_pt_b;
    _strokes.splice(_strokes.end(), chain->_strokes);
  }
  else // These chains shouldn't be joined.
    return false;

  StrokeList::iterator it = _strokes.begin();
  StrokeList::iterator it_end = _strokes.end();
  for (; it != it_end; ++it)
    (*it)->setChain(this);

  delete(chain);

  return true;
}

Chain *Chain::split(StrokeList::iterator str_it)
{
  if (str_it == _strokes.begin())
    return 0;

  if (cycle()) {
    _pt_a = (*str_it)->pointA();
    _pt_b = (*str_it)->pointA();
    _strokes.splice(_strokes.end(), _strokes, _strokes.begin(), str_it);
    return 0;
  }
  
  StrokeList::iterator str_it_prev = str_it;
  --str_it_prev;
    
  Chain *new_chain = new Chain();
  new_chain->setPointA((*str_it)->pointA());
  new_chain->setPointB(pointB());
  setPointB((*str_it)->pointA());
  new_chain->_strokes.splice(new_chain->_strokes.end(), _strokes, str_it, _strokes.end());
  StrokeList::iterator it = new_chain->_strokes.begin();
  StrokeList::iterator it_end = new_chain->_strokes.end();
  for (; it != it_end; ++it)
    (*it)->setChain(new_chain);

  return new_chain;
}

Chain *Chain::split(Stroke *stroke)
{
  if (!stroke)
    return 0;

  return split(std::find(_strokes.begin(), _strokes.end(), stroke));
}

void Chain::removeStroke(StrokeList::iterator str_it)
{
  if (_strokes.empty())
    return;

  if (_strokes.size() == 1) {
    _pt_a = counted_ptr<Point>(0);
    _pt_b = counted_ptr<Point>(0);
  }
  else {
    if (_pt_a.get() == (*str_it)->pointA().get())
      _pt_a = (*str_it)->pointB();
    else if (_pt_a.get() == (*str_it)->pointB().get())
      _pt_a = (*str_it)->pointA();
    else if (_pt_b.get() == (*str_it)->pointA().get())
      _pt_b = (*str_it)->pointB();
    else if (_pt_b.get() == (*str_it)->pointB().get())
      _pt_b = (*str_it)->pointA();
    else // Cannot delete this stroke, which is not an extremity of the chain.
      return;
  }
  
  delete *str_it;
  _strokes.erase(str_it);
}

void Chain::removeStroke(Stroke *stroke)
{
  if (!stroke)
    return;

  return removeStroke(std::find(_strokes.begin(), _strokes.end(), stroke));
}

void Chain::recalculateBB() 
{
  Vec3d minBB(999999,999999,0), maxBB(-9999999,-999999,0); // store maximum extents of the seam
    
  Chain::StrokeList::iterator st;
  Chain::StrokeList::iterator st_end;
  
  Stroke::SegmentList::iterator seg;
  Stroke::SegmentList::iterator seg_end;
  
  st = this->strokes()->begin();
  st_end = this->strokes()->end();
    
    // calculate the extents of the seam
  for (; st != st_end; ++st) {
    
      (*st)->recalculateBB();
        // Update X extents if required
      if( (*st)->bbMin()[0] < minBB.x()) {
        minBB[0] = (*st)->bbMin()[0] ;
      }
      if( (*st)->bbMin()[1] < minBB.y()) {
        minBB[1] = (*st)->bbMin()[1];
      }
      if( (*st)->bbMax()[0] > maxBB.x()) {
        maxBB[0] = (*st)->bbMax()[0];
      }
      if( (*st)->bbMax()[1] > maxBB.y()) {
        maxBB[1] = (*st)->bbMax()[1] ;
      }
    
  }
  _bbMin = minBB;
  _bbMax = maxBB;
}
