//
//  Filename         : repository_pretreatment_result.cpp
//  Author(s)        : Emmanuel Turquin
//  Purpose          : Common repositoty for Pretreatment and Result modules.
//  Date of creation : 05/22/2004
//
///////////////////////////////////////////////////////////////////////////////

// LICENSE

#include "result.h"
#include "pretreatment.h"
#include "repository_pretreatment_result.h"

Pretreatment	*RepositoryPretreatmentResult::_pre = 0;
Result		*RepositoryPretreatmentResult::_res = 0;
Vec3d		RepositoryPretreatmentResult::_model_bb_min;
Vec3d		RepositoryPretreatmentResult::_model_bb_max;
Mat44d		RepositoryPretreatmentResult::_matrix = Mat44d::identity();
int		RepositoryPretreatmentResult::_model_dl = 0;
unsigned char	RepositoryPretreatmentResult::_modified = 0;

void RepositoryPretreatmentResult::pretreatmentUpdated()
{
  if (_res) {
    _res->makeCurrent();
    _res->pretreatmentUpdated();
    _modified = 0;
    if (_pre)
      _pre->makeCurrent();
  }
}

void RepositoryPretreatmentResult::resultUpdated()
{
  if (_pre) {
    _pre->makeCurrent();
    _pre->resultUpdated();
    _modified = 0;
    if (_res)
      _res->makeCurrent();
  }
}
