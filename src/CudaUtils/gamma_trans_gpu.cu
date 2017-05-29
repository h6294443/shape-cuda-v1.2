/***************************************************************************
                                                              gamma_trans.c

Applies a gamma transformation to a data point.

Modified 2004 Feb 13 by CM:
    Removed obsolete "sdev" argument and the commented-out code
    which used to make use of that argument
***************************************************************************/
extern "C" {
#include "../shape/head.h"
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
__global__ void cf_gamma_trans_streams_krnl(struct par_t *dpar, struct dat_t *ddat,
		int s, int f, int nThreads, unsigned char type, int flt) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	if (offset < nThreads) {
		/*  Carry out a gamma transformation on the fit image if requested  */
		if (dpar->dd_gamma != 1.0) {
			switch (type) {
			case DELAY:
				if (flt)
					dev_gamma_trans_f(&ddat->set[s].desc.deldop.frame[f].fit_s[offset],
							__double2float_rn(dpar->dd_gamma));
				else
					dev_gamma_trans(&ddat->set[s].desc.deldop.frame[f].fit_s[offset],
						dpar->dd_gamma);
				break;
			case DOPPLER:
				//cf_dop_frame->fit[offset] = fit[offset];
				break;
			}
		}
	}
}
