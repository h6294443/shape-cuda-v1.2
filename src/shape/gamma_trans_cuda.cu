/***************************************************************************
                                                              gamma_trans.c

Applies a gamma transformation to a data point.

Modified 2004 Feb 13 by CM:
    Removed obsolete "sdev" argument and the commented-out code
    which used to make use of that argument
***************************************************************************/
extern "C" {
#include "head.h"
}

__device__ int dev_gamma_trans(float *datum, double gamma)
{
  if ((*datum) <= 0.0)
    return 0;
  (*datum) = pow( (double)(*datum), 1/gamma);
  return 1;
}
__device__ int dev_gamma_trans_f(float *datum, float gamma)
{
  if ((*datum) <= 0.0)
    return 0;
  (*datum) = __powf( (*datum), 1/gamma);
  return 1;
}
