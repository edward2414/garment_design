//
//  Filename         : repository_canvas_result.cpp
//  Author(s)        : Emmanuel Turquin
//  Purpose          : Common repositoty for Canvas and Result modules.
//  Date of creation : 05/22/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include "result.h"
#include "canvas.h"
#include "config_canvas.h"
#include "repository_canvas_result.h"

Canvas		*RepositoryCanvasResult::_can = 0;
Result		*RepositoryCanvasResult::_res = 0;
unsigned char	RepositoryCanvasResult::_modified = 0;
Vec2d		RepositoryCanvasResult::_border_factor(config::CANVAS_BORDER_X, config::CANVAS_BORDER_Y);
const Layer		*RepositoryCanvasResult::_front_layer = 0;
const Layer     *RepositoryCanvasResult::_back_layer = 0;
const Layer     *RepositoryCanvasResult::_front_seam_layer = 0;
const Layer     *RepositoryCanvasResult::_back_seam_layer = 0;

const RepositoryCanvasResult::GarmentMapList *RepositoryCanvasResult::_front_garment_maps = 0;
const RepositoryCanvasResult::GarmentMapList *RepositoryCanvasResult::_back_garment_maps = 0;

void RepositoryCanvasResult::canvasUpdated()
{
  if (_res) {
    _res->makeCurrent();
    _res->canvasUpdated();
    _modified = 0;
    if (_can)
      _can->makeCurrent();
  }
}

void RepositoryCanvasResult::resultUpdated()
{
  if (_can) {
    _can->makeCurrent();
    _can->resultUpdated();
    _modified = 0;
    if (_res)
      _res->makeCurrent();
  }
}
