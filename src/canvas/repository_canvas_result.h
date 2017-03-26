//
//  Filename         : repository_canvas_result.h
//  Author(s)        : Emmanuel Turquin
//  Purpose          : Common repositoty for Canvas and Result modules.
//  Date of creation : 05/22/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  REPOSITORY_CANVAS_RESULT_H
# define REPOSITORY_CANVAS_RESULT_H

# include <list> // FIXME DEBUG
# include "vectypes.h"

class Canvas;
class Result;
class GarmentMap;
class Layer;

class RepositoryCanvasResult
{
public:

  typedef std::list<GarmentMap*>	GarmentMapList;

  static const unsigned char BORDER_FACTOR =	1 << 0;
  static const unsigned char LAYERS =		1 << 1;
  static const unsigned char GARMENT_MAPS =	1 << 2;

  static inline void setCanvas(Canvas *can) {
    _can = can;
  }

  static inline void setResult(Result *res) {
    _res = res;
  }

  static inline Canvas *canvas() {
    return _can;
  }

  static inline Result *result() {
    return _res;
  }

  static inline void setBorderFactor(const Vec2d& v) {
    _border_factor = v;
   _modified |= BORDER_FACTOR;
  }

  static inline const Vec2d& borderFactor() {
    return _border_factor;
  }

  static inline void setFrontLayer(const Layer* layer) {
    _front_layer = layer;
    _modified |= LAYERS;
  }

  static inline void setFrontSeamLayer(const Layer* layer) {
    _front_seam_layer = layer;
    _modified |= LAYERS;
  }
  
  static inline void setBackSeamLayer(const Layer* layer) {
    _back_seam_layer = layer;
    _modified |= LAYERS;
  }
  
  static inline const Layer *frontLayer() {
    return _front_layer;
  }

  static inline const Layer *frontSeamLayer() {
    return _front_seam_layer;
  }

  static inline const Layer *backSeamLayer() {
    return _back_seam_layer;
  }
  
  static inline void setBackLayer(const Layer* layer) {
    _back_layer = layer;
    _modified |= LAYERS;
  }

  static inline const Layer *backLayer() {
    return _back_layer;
  }

  static inline void setFrontGarmentMaps(const GarmentMapList *gm) {
    _front_garment_maps = gm;
    _modified |= GARMENT_MAPS;
  }

  static inline const GarmentMapList *frontGarmentMaps() {
    return _front_garment_maps;
  }

  static inline void setBackGarmentMaps(const GarmentMapList *gm) {
    _back_garment_maps = gm;
    _modified |= GARMENT_MAPS;
  }

  static inline const GarmentMapList *backGarmentMaps() {
    return _back_garment_maps;
  }

  static void canvasUpdated();
  static void resultUpdated();

  static inline unsigned char modified() {
    return _modified;
  }

private:

  static Canvas		*_can;
  static Result		*_res;
  static unsigned char	_modified;
  static Vec2d		_border_factor;
  static const GarmentMapList	*_front_garment_maps;
  static const GarmentMapList	*_back_garment_maps;
  static const Layer		*_front_layer;
  static const Layer        *_back_layer;
  static const Layer      *_front_seam_layer;
  static const Layer      *_back_seam_layer;
};

#endif // REPOSITORY_CANVAS_RESULT_H
