//
//  Filename         : repository_pretreatment_result.h
//  Author(s)        : Emmanuel Turquin
//  Purpose          : Common repositoty for Pretreatment and Result modules.
//  Date of creation : 05/22/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#ifndef  REPOSITORY_PRETREATMENT_RESULT_H
# define REPOSITORY_PRETREATMENT_RESULT_H

# include "vectypes.h"
# include "mattypes.h"

class Pretreatment;
class Result;

class RepositoryPretreatmentResult
{
public:

  static const unsigned char MODEL_BB =	1 << 0;
  static const unsigned char MODEL_DL =	1 << 1;
  static const unsigned char MATRIX =	1 << 2;

  static inline void setPretreatment(Pretreatment *pre) {
    _pre = pre;
  }

  static inline void setResult(Result *res) {
    _res = res;
  }

  static inline Pretreatment *pretreatment() {
    return _pre;
  }

  static inline Result *result() {
    return _res;
  }

  static void pretreatmentUpdated();
  static void resultUpdated();

  static inline void setModelBBMin(const Vec3d& min) {
    _model_bb_min = min;
   _modified |= MODEL_BB;
  }

  static inline const Vec3d& modelBBMin() {
    return _model_bb_min;
  }

  static inline void setModelBBMax(const Vec3d& max) {
    _model_bb_max = max;
    _modified |= MODEL_BB;
  }

  static inline const Vec3d& modelBBMax() {
    return _model_bb_max;
  }

  static inline Vec3d modelBBSize() {
    return _model_bb_max - _model_bb_min;
  }

  static inline void setMatrix(const Mat44d& matrix) {
    _matrix = matrix;
    _modified |= MATRIX;
  }

  static inline const Mat44d& matrix() {
    return _matrix;
  }

  static inline void setModelDisplayList(int dl) {
    _model_dl = dl;
   _modified |= MODEL_DL;
  }

  static inline int modelDisplayList() {
    return _model_dl;
  }

  static inline unsigned char modified() {
    return _modified;
  }

private:

  static Pretreatment	*_pre;
  static Result		*_res;
  static Vec3d		_model_bb_min;
  static Vec3d		_model_bb_max;
  static Mat44d		_matrix;
  static int		_model_dl;
  static unsigned char	_modified;
};

#endif // REPOSITORY_PRETREATMENT_RESULT_H
