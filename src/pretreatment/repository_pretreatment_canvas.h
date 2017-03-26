//
//  Filename         : repository_pretreatment_canvas.h
//  Author(s)        : Emmanuel Turquin
//  Purpose          : Common repositoty for Pretreatment and Canvas modules.
//  Date of creation : 05/22/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  REPOSITORY_PRETREATMENT_CANVAS_H
# define REPOSITORY_PRETREATMENT_CANVAS_H

# include "vectypes.h"
# include "distance_field.h"

class Pretreatment;
class Canvas;

class RepositoryPretreatmentCanvas
{
public:

  static const unsigned char TEXTURES =		1 << 0;
  static const unsigned char BOUNDING_BOX =	1 << 1;
  static const unsigned char DISTANCE_FIELD =	1 << 2;

  static inline void setPretreatment(Pretreatment *pre) {
    _pre = pre;
  }

  static inline void setCanvas(Canvas *can) {
    _can = can;
  }

  static inline Pretreatment *pretreatment() {
    return _pre;
  }

  static inline Canvas *canvas() {
    return _can;
  }

  static void pretreatmentUpdated();
  static void canvasUpdated();

  static inline void setBBMin(const Vec2d& bb_min) {
    _bb_min = bb_min;
    _modified |= BOUNDING_BOX;
  }

  static inline const Vec2d& bbMin() {
    return _bb_min;
  }

  static inline void setBBMax(const Vec2d& bb_max) {
    _bb_max = bb_max;
    _modified |= BOUNDING_BOX;
  }

  static inline const Vec2d& bbMax() {
    return _bb_max;
  }

  static inline void setFrontTexture(unsigned char *texture, unsigned texture_size) {
    _front_texture = texture;
    _front_texture_size = texture_size;
    _modified |= TEXTURES;
  }

  static inline unsigned frontTextureSize() {
    return _front_texture_size;
  }

  static inline void freeFrontTexture() {
    delete[] _front_texture;
    _front_texture = 0;
    _front_texture_size = 0;
    _modified |= TEXTURES;
  }

  static inline unsigned char* frontTexture() {
    return _front_texture;
  }

  static inline void setBackTexture(unsigned char *texture, unsigned texture_size) {
    _back_texture = texture;
    _back_texture_size = texture_size;
    _modified |= TEXTURES;
  }

  static inline unsigned backTextureSize() {
    return _back_texture_size;
  }

  static inline void freeBackTexture() {
    delete[] _back_texture;
    _back_texture = 0;
    _back_texture_size = 0;
    _modified |= TEXTURES;
  }

  static inline unsigned char* backTexture() {
    return _back_texture;
  }

  static inline void setDistanceField(DistanceField *df) {
    _df = df;
    _modified |= DISTANCE_FIELD;
  }

  static inline DistanceField *distanceField() {
    return _df;
  }

  static inline void freeDistanceField() {
    delete _df;
    _df = 0;
    _modified |= DISTANCE_FIELD;
  }

  static inline unsigned char modified() {
    return _modified;
  }

private:

  static Pretreatment	*_pre;
  static Canvas		*_can;
  static Vec2d		_bb_min;
  static Vec2d		_bb_max;
  static unsigned char	*_front_texture;
  static unsigned	_front_texture_size;
  static unsigned char	*_back_texture;
  static unsigned	_back_texture_size;
  static DistanceField	*_df;
  static unsigned char	_modified;
};

#endif // REPOSITORY_PRETREATMENT_CANVAS_H
