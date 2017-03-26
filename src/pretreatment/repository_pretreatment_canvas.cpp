//
//  Filename         : repository_pretreatment_canvas.cpp
//  Author(s)        : Emmanuel Turquin
//  Purpose          : Common repositoty for Pretreatment and Canvas modules.
//  Date of creation : 05/22/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include "canvas.h"
#include "pretreatment.h"
#include "repository_pretreatment_canvas.h"

Pretreatment	*RepositoryPretreatmentCanvas::_pre = 0;
Canvas		*RepositoryPretreatmentCanvas::_can = 0;
DistanceField	*RepositoryPretreatmentCanvas::_df = 0;
unsigned char	*RepositoryPretreatmentCanvas::_front_texture = 0;
unsigned	RepositoryPretreatmentCanvas::_front_texture_size = 0;
unsigned char	*RepositoryPretreatmentCanvas::_back_texture = 0;
unsigned	RepositoryPretreatmentCanvas::_back_texture_size = 0;
Vec2d		RepositoryPretreatmentCanvas::_bb_min;
Vec2d		RepositoryPretreatmentCanvas::_bb_max;
unsigned char	RepositoryPretreatmentCanvas::_modified = 0;

void RepositoryPretreatmentCanvas::pretreatmentUpdated()
{
  if (_can) {
    _can->makeCurrent();
    _can->pretreatmentUpdated();
    _modified = 0;
    if (_pre)
      _pre->makeCurrent();
  }
}

void RepositoryPretreatmentCanvas::canvasUpdated()
{
  if (_pre) {
    _pre->makeCurrent();
    _pre->canvasUpdated();
    _modified = 0;
    if (_can)
      _can->makeCurrent();
  }
}
