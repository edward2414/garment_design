#ifndef  PUSH_BACK_PARAMS_H
# define PUSH_BACK_PARAMS_H
#include "vectypes.h"

class PushBackParams
{
  public:
  Vec3d startPoint;
  Vec3d endPoint;
  double shrunkLength;
  double a;
  double b;
  double influenceRadius;
  bool valid;
};

#endif // PUSH_BACK_PARAMS_H

