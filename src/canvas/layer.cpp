//
//  Filename         : layer.cpp
//  Author           : Emmanuel Turquin
//  Purpose          : Class representing a layer of chains for the Canvas.
//  Date of creation : 04/26/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include <set>
#include "utils.h"
#include "grid.h"
#include "canvas.h"
#include "layer.h"
#include "point.h"


#define MSG(type, msg) if (_parent) _parent->canvasMessage(type, msg);

ostream& operator<<(std::ostream& out, const Layer& layer) 
{
  // print out the information required to reconstruct the strokes in the layer
  std::list<Chain *>::const_iterator it = layer.chains()->begin();
  std::list<Chain *>::const_iterator it_end = layer.chains()->end();
  
  out << "Layer{" << std::endl;

  out << layer._chains.size() << std::endl; // number of chains in this layer
  
  for(;it != it_end;it++)
  {
    Chain* ptrChain = *it;
    std::list<Stroke *>::const_iterator sit = ptrChain->strokes()->begin();
    std::list<Stroke *>::const_iterator sit_end = ptrChain->strokes()->end();
    
    out << "Chain{" << std::endl;
    out << ptrChain->strokes()->size() << std::endl; 
    
    for(;sit != sit_end;sit++)
    {
      Stroke* ptrStroke = *sit;
      std::list<Segment *>::const_iterator git = ptrStroke->segments()->begin();
      std::list<Segment *>::const_iterator git_end = ptrStroke->segments()->end();
      
      out << "Stroke{" << std::endl;
      out << ptrStroke->segments()->size()+1 << std::endl; // plus one because we want number of points
      out << ptrStroke->getAParams()[0] << std::endl;
      out << ptrStroke->getAParams()[1] << std::endl;
      out << ptrStroke->getBParams()[0] << std::endl;
      out << ptrStroke->getBParams()[1] << std::endl;
      
      for(;git!=git_end;git++)
      {
        Segment* ptrSegment = *git;
        if(git == ptrStroke->segments()->begin())
        {
          out << *(ptrSegment->pointA()) << std::endl; // only print point A for the beginning of the stroke
        }
        out << *(ptrSegment->pointB()) << std::endl;
      }
      
      out << "}" << std::endl; // end Stroke
    }
    
    out << "}" << std::endl; // end Chain
    
  }
  out << "}" << std::endl; // end Layer
  return(out);
}

Layer::Layer(Canvas *parent, const Vec2d& min, const Vec2d& max, int sampling, int snapping, double splitting_threshold) {
  _parent = parent;
  _current_chain_begin = 0;
  _current_chain_end = 0;
  _current_stroke = 0;
  _sampling = sampling;
  _snapping = snapping;
  _splitting_threshold = splitting_threshold;
  _grid = new Grid(config::GRID_WIDTH, config::GRID_HEIGHT, min[0], min[1], max[0], max[1]);
  _ltype = GARMENT;
}

Layer::Layer(const Layer& layer) :
  _chains(layer._chains),
  _grid(layer._grid),
  _collisions(layer._collisions)
{
  _parent = layer._parent;
  _current_chain_begin = layer._current_chain_begin;
  _current_chain_end = layer._current_chain_end;
  _current_stroke = layer._current_stroke;
  _current_point = layer._current_point;
  _sampling = layer._sampling;
  _snapping = layer._snapping;
  _splitting_threshold = layer._splitting_threshold;
}

Layer::~Layer() {
  clear();
  delete _grid;
}

void Layer::clear() {
  // Delete chains.
  ChainList::iterator it1 = _chains.begin();
  ChainList::iterator it_end1 = _chains.end();
  for ( ; it1 != it_end1; ++it1)
    delete *it1;
  _chains.clear();

  // Delete current stroke.
  delete _current_stroke;

  // Clears the grid and collisions list.
  _grid->clear();
  _collisions.clear();

  // reinitialize all pointers.
  _current_chain_begin = 0;
  _current_chain_end = 0;
  _current_stroke = 0;
  _current_point = counted_ptr<Point>(0);
}

void Layer::firstPoint(double x, double y, double ratio) {
  _current_chain_begin = 0;
  ChainList::iterator it = _chains.begin();
  ChainList::iterator it_end = _chains.end();
  counted_ptr<Point> p;

  // JDW Cycle through existing (non-cyclic) chains in this layer. If this point is close to an existing
  // start or end point then snap to the existing point
  
  for ( ; it != it_end; ++it) {
    if ((*it)->cycle())
      continue;
    if ((*it)->pointA().get() && utils::isInSphere(x, y, (*it)->pointA()->x(), (*it)->pointA()->y(), ratio * _snapping)) {
      p = (*it)->pointA();
      _current_chain_begin = *it;
      break;
    }
    if ((*it)->pointB().get() && utils::isInSphere(x, y, (*it)->pointB()->x(), (*it)->pointB()->y(), ratio * _snapping)) {
      p = (*it)->pointB();
      _current_chain_begin = *it;
      break;
    }
  }

  if (!p.get())
    p = counted_ptr<Point>(new Point(x, y));
  _current_point = p;

  Stroke *s = new Stroke();
  _current_stroke = s;
}

void Layer::intermediatePoint(double x, double y, double ratio) {
  if (!_current_stroke)
    return;

  // JDW Don't add a new point unless it's sufficiently far from the last one. 
  unsigned seg_nb = utils::segmentsInSphere(x, y, _current_point->x(), _current_point->y(), ratio * _sampling);

  if (!seg_nb)
    return;
  
  double dx = (x - _current_point->x()) / seg_nb;
  double dy = (y - _current_point->y()) / seg_nb;

  for (; seg_nb > 0; --seg_nb) {
    counted_ptr<Point> p(new Point(_current_point->x() + dx, _current_point->y() + dy));
    Segment *s = new Segment(_current_point, p);

    _current_stroke->addSegment(s);

    // Collision test.
    Cell::SegmentList l = _grid->getSegments(s);
    Cell::SegmentList::const_iterator it = l.begin();
    Cell::SegmentList::const_iterator it_end = l.end();
    for (; it != it_end; ++it)
      Collision::collisionTest(s, (*it), _collisions);

    _grid->addSegment(s);
    _current_point = p;
  }
}

void Layer::lastPoint(double x, double y, double ratio) {
  if (!_current_stroke)
    return;

  if (_current_stroke->segments()->empty()) {
    delete _current_stroke;
    _current_stroke = 0;
    return;
  }

  _current_chain_end = 0;
  ChainList::iterator it = _chains.begin();
  ChainList::iterator it_end = _chains.end();
  counted_ptr<Point> p;
  if (!_current_chain_begin &&
      utils::isInSphere(x, y, _current_stroke->pointA()->x(), _current_stroke->pointA()->y(), ratio * _snapping))
    p = _current_stroke->pointA();
  else {
    for ( ; it != it_end; ++it) {
      if ((*it)->cycle())
	continue;
      if (utils::isInSphere(x, y, (*it)->pointA()->x(), (*it)->pointA()->y(), ratio * _snapping)) {
	p = (*it)->pointA();
	_current_chain_end = *it;
	break;
      }
      if (utils::isInSphere(x, y, (*it)->pointB()->x(), (*it)->pointB()->y(), ratio * _snapping)) {
	p = (*it)->pointB();
	_current_chain_end = *it;
	break;
      }
    }
  }

  // Create last segment if needed
  if (p.get() || !utils::isInSphere(x, y, _current_point->x(), _current_point->y(), ratio * _sampling)) {
    if (!p.get())
      p = counted_ptr<Point>(new Point(x, y));
    Segment *s = new Segment(_current_point, p);
    _current_point = p;
    _current_stroke->addSegment(s);

    // Collision test.
    Cell::SegmentList l = _grid->getSegments(s);
    Cell::SegmentList::const_iterator it = l.begin();
    Cell::SegmentList::const_iterator it_end = l.end();
    for (; it != it_end; ++it)
      Collision::collisionTest(s, (*it), _collisions);

    _grid->addSegment(s);
  }

  if (collisionTreatment())
    return;

  MSG(MsgHandler::MSG_NORMAL, "New stroke created.");

  _current_stroke->updateType();

  // No chain present at either end of stroke, create a new one.
  if (!_current_chain_begin && !_current_chain_end) {
    Chain *c = new Chain(_current_stroke);
    _chains.push_back(c);
    MSG(MsgHandler::MSG_NORMAL, "New chain created.");
  }

  // A chain already exists at the beginning.
  else if (!_current_chain_end) {
    _current_chain_begin->addStroke(_current_stroke);
    MSG(MsgHandler::MSG_NORMAL, "Stroke added to an existing chain.");
  }

  // A chain already exists at the end.
  else if (!_current_chain_begin) {
    _current_chain_end->addStroke(_current_stroke);
    MSG(MsgHandler::MSG_NORMAL, "Stroke added to an existing chain.");
  }

  // A chain is closed (a cycle is created).
  else if (_current_chain_begin == _current_chain_end) {
    _current_chain_begin->addStroke(_current_stroke);
    MSG(MsgHandler::MSG_NORMAL, "Chain closed (cycle present).");
  }

  // Two chains are joined by the stroke, merge them.
  else {
    _current_chain_begin->addStroke(_current_stroke);
    _current_chain_begin->join(_current_chain_end);
    _chains.remove(_current_chain_end);
    MSG(MsgHandler::MSG_NORMAL, "Chain joined with another one.");
  }

  curvatureTreatment(_current_stroke);
  _current_stroke = 0;
}

void Layer::curvatureTreatment(Stroke *stroke)
{
  if (!stroke)
    return;

  Stroke::SegmentList::iterator it = stroke->segments()->begin();
  Stroke::SegmentList::iterator it_end = stroke->segments()->end();
  Segment *seg = *it;
  double dotp;
  for (++it; it != it_end; ++it) {
    dotp = utils::normalizedDotProduct(seg->pointA()->x() - seg->pointB()->x(),
				       seg->pointA()->y() - seg->pointB()->y(),
				       (*it)->pointA()->x() - (*it)->pointB()->x(),
				       (*it)->pointA()->y() - (*it)->pointB()->y());
    if (dotp < _splitting_threshold) {
      // In this situation, split will always work, no need to test for a null result.
      Stroke *new_stroke = stroke->split(it);
      MSG(MsgHandler::MSG_NORMAL, "Stroke automatically splitted.");
      curvatureTreatment(new_stroke);
      return;
    }
    seg = *it;
  }
}

bool Layer::collisionTreatment()
{
  if (_collisions.empty())
    return false;

  Collision::List::iterator it = _collisions.begin();
  Collision::List::iterator it_end = _collisions.end();
  Collision::List self_collisions, ext_collisions;
  for (; it != it_end; ++it)
    if (it->self())
      self_collisions.push_back(*it);
    else
      ext_collisions.push_back(*it);
  unsigned self_count = self_collisions.size();
  unsigned ext_count = ext_collisions.size();
  unsigned total_count = self_count + ext_count;
  Segment *seg;
  Chain *c1, *c2;

  // Rough gestures (self and external collisions counts serve as borns)
  if (total_count >= 40 && ext_count) {
    MSG(MsgHandler::MSG_NORMAL, "Current layer cleared.");
    clear();
  }
  else if (total_count >= 5 && ext_count) {
    std::set<Stroke*> strokes;
    std::set<Chain*> chains;
    it = ext_collisions.begin();
    it_end = ext_collisions.end();
    for (; it != it_end; ++it) {
      strokes.insert(it->segment2()->stroke());
      chains.insert(it->segment2()->stroke()->chain());
    }
    if (chains.size() > 1) {
      std::set<Chain*>::iterator ch = chains.begin();
      std::set<Chain*>::iterator ch_end = chains.end();
      for (; ch != ch_end; ++ch) {
	_chains.remove(*ch);
	delete *ch;
      }
      MSG(MsgHandler::MSG_NORMAL, QString::number(chains.size()) + " chains deleted.");
    }
    else {
      MSG(MsgHandler::MSG_NORMAL, QString::number(strokes.size()) + " stroke(s) deleted.");
      std::set<Stroke*>::iterator st = strokes.begin();
      std::set<Stroke*>::iterator st_end = strokes.end();
      for (; st != st_end; ++st) {
	c1 = (*st)->chain();
	c2 = c1->split(*st);
	if (c2) {
	  _chains.push_back(c2);
	  MSG(MsgHandler::MSG_NORMAL, "Chain splitted.");
	}
	else
	  c2 = c1;
	c2->removeStroke(*st);
	if (c2->strokes()->empty()) {
	  _chains.remove(c2);
	  delete c2;
	  MSG(MsgHandler::MSG_NORMAL, "Chain deleted.");
	}
      }
    }
  }

  // Exact count (self and external collisions)
  else {
    switch (self_count) {
    case 0: // self_count = 0
      switch (ext_count) {
        case 1: // ext_count = 1
          if(_ltype != Layer::GFOLD) 
          { // JDW usual case 
            it = ext_collisions.begin();
            seg =  it->segment2();
            seg->stroke()->cleverSplit(seg);
            MSG(MsgHandler::MSG_NORMAL, "Stroke splitted.");
            break;
          }
          else
          { // gfold case
            it = ext_collisions.begin();
            seg =  it->segment2();
            seg->stroke()->updateLength();
            // Iterate over segments in intersected stroke to determine how far along length of stroke we are
            
            Stroke::SegmentList::iterator sit = seg->stroke()->segments()->begin();
            Stroke::SegmentList::iterator sit_end = seg->stroke()->segments()->end();
            double currentLength = 0.0;
            for(;sit!= sit_end;sit++) 
            {
              Segment *s = *sit;
              s->updateLength();
              currentLength += s->length();
              if (s == seg) {
                break;
              }
            }
            double percentageLength = currentLength / seg->stroke()->length();
            double dx = currentStroke()->pointA()->x() - currentStroke()->pointB()->x();
            double dy = currentStroke()->pointA()->y() - currentStroke()->pointB()->y();
            double dist = sqrt ( dx*dx + dy*dy ) /2.0; // use as radius of gaussian
            
            
            
            // find gaussian height parameter
            double maxHeight = -10.0;
            Vec3d pa = Vec3d(currentStroke()->pointA()->x(), currentStroke()->pointA()->y(), 0.0);
            Vec3d pb = Vec3d(currentStroke()->pointB()->x(), currentStroke()->pointB()->y(), 0.0);
            sit = currentStroke()->segments()->begin();
            sit_end = currentStroke()->segments()->end();
            Vec3d maxPoint;
            if(sit != sit_end) {
              double d = utils::distPointLine( pa, pb, 
                                               Vec3d((*sit)->pointA()->x(), (*sit)->pointA()->y(), 0.0));
              if(d >= maxHeight) {
                maxHeight = d;
                maxPoint.x() = (*sit)->pointA()->x();
                maxPoint.y() = (*sit)->pointA()->y();
              }
            }
            
            for(;sit != sit_end;sit++)
            {
              double d = utils::distPointLine( pa, pb, 
                                               Vec3d((*sit)->pointB()->x(), (*sit)->pointB()->y(), 0.0));
              if(d >= maxHeight) {
                maxHeight = d;
                maxPoint.x() = (*sit)->pointB()->x();
                maxPoint.y() = (*sit)->pointB()->y();
              }
            }
            
            Vec3d gaussDir = utils::vecFromPointToLine( pa, pb, maxPoint );
            
            Vec3d a( seg->pointA()->x(), seg->pointA()->y(), 0.0);
            Vec3d b( seg->pointB()->x(), seg->pointB()->y(), 0.0);

            if(percentageLength <= 0.5) {
              
              // We are closer to the beginning of the fold stroke
              Vec3d segDir = b - a; 
              double dp = gaussDir * segDir;
              double sense = 1.0;
              if(dp >= 0.0) // same direction, ridge fold
                sense = 1.0;
              else // different direction, valley fold (default)
                sense = -1.0;

              seg->stroke()->setAParams( Vec2d(maxHeight*sense,dist) );
              MSG(MsgHandler::MSG_NORMAL, "Radius at A  set to: " + QString::number(dist) + ", Height at A set to: " + QString::number(maxHeight));
            }
            else {
              
              // We are closer to the end of the fold stroke
              Vec3d segDir = a - b; 
              double dp = gaussDir * segDir;
              double sense = 1.0;
              if(dp >= 0.0) // same direction, ridge fold
                sense = 1.0;
              else // different direction, valley fold (previous default)
                sense = -1.0;
            
              seg->stroke()->setBParams( Vec2d(maxHeight*sense,dist) );
              MSG(MsgHandler::MSG_NORMAL, "Radius at B  set to: " + QString::number(dist) + ", Height at B set to: " + QString::number(maxHeight));
            }
            
            break;
          }
        default:
            MSG(MsgHandler::MSG_WARNING, "No action specified for a " + QString::number(ext_count) +
                " external and " + QString::number(self_count) + " self collision(s) gesture.");
      }
      break;
    case 1: // self_count = 1
      switch (ext_count) {
      case 2: // ext_count = 2
	if (ext_collisions.front().segment2()->stroke()->join(ext_collisions.back().segment2()->stroke()))
	  MSG(MsgHandler::MSG_NORMAL, "Strokes joined.");
	break;
      default:
	MSG(MsgHandler::MSG_WARNING, "No action specified for a " + QString::number(ext_count) +
	    " external and " + QString::number(self_count) + " self collision(s) gesture.");
      }
      break;
      break;
    default:
      MSG(MsgHandler::MSG_WARNING, "No action specified for a " + QString::number(ext_count) +
	  " external and " + QString::number(self_count) + " self collision(s) gesture.");
    }
  }

  delete _current_stroke;
  _current_stroke = 0;
  _collisions.clear();
  return true;
}
