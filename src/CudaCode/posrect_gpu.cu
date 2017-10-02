
extern "C" {
#include "../shape/head.h"
#include <limits.h>
}

__device__ static float atomicMinf(float* address, float val) {
	int* address_as_i = (int*) address;
	int old = *address_as_i, assumed;
	do {
		assumed = old;
		old = ::atomicCAS(address_as_i, assumed,
				__float_as_int(::fminf(val, __int_as_float(assumed))));
	} while (assumed != old);
	return __int_as_float(old);
}
__device__ static float atomicMin64(double* address, double val)
{
	unsigned long long* address_as_i = (unsigned long long*) address;
	unsigned long long old = *address_as_i, assumed;
	do {
		assumed = old;
		old = ::atomicCAS(address_as_i, assumed,
				__double_as_longlong(::fminf(val, __longlong_as_double(assumed))));
	} while (assumed != old);
	return __longlong_as_double(old);
}
__device__ static float atomicMaxf(float* address, float val) {
	int* address_as_i = (int*) address;
	int old = *address_as_i, assumed;
	do {
		assumed = old;
		old = ::atomicCAS(address_as_i, assumed,
				__float_as_int(::fmaxf(val, __int_as_float(assumed))));
	} while (assumed != old);
	return __int_as_float(old);
}
__device__ static float atomicMax64(double* address, double val)
{
	unsigned long long* address_as_i = (unsigned long long*) address;
	unsigned long long old = *address_as_i, assumed;
	do {
		assumed = old;
		old = ::atomicCAS(address_as_i, assumed,
				__double_as_longlong(::fmaxf(val, __longlong_as_double(assumed))));
	} while (assumed != old);
	return __longlong_as_double(old);
}
__device__ void dev_POSrect_gpu32(
		struct pos_t **pos,
		int src,
		float imin_dbl,
		float imax_dbl,
		float jmin_dbl,
		float jmax_dbl,
		float4 *ijminmax_overall,
		int frm)	{
	int n, imin, imax, jmin, jmax;
	n = pos[frm]->n;

	/* Update the POS region that contains the target without
	 * regard to whether or not it extends beyond the POS frame */
	atomicMinf(&ijminmax_overall[frm].w, imin_dbl);
	atomicMaxf(&ijminmax_overall[frm].x, imax_dbl);
	atomicMinf(&ijminmax_overall[frm].y, jmin_dbl);
	atomicMaxf(&ijminmax_overall[frm].z, jmax_dbl);

	/*  Update the subset of the POS frame that contains the target  */
	imin = (imin_dbl < INT_MIN) ? INT_MIN : (int) imin_dbl;
	imax = (imax_dbl > INT_MAX) ? INT_MAX : (int) imax_dbl;
	jmin = (jmin_dbl < INT_MIN) ? INT_MIN : (int) jmin_dbl;
	jmax = (jmax_dbl > INT_MAX) ? INT_MAX : (int) jmax_dbl;

	/* Make sure it's smaller than n */
	imin = MAX(imin,-n);
	imax = MIN(imax, n);
	jmin = MAX(jmin,-n);
	jmax = MIN(jmax, n);

	if (src) {
		atomicMin(&pos[frm]->xlim2[0], imin);
		atomicMax(&pos[frm]->xlim2[1], imax);
		atomicMin(&pos[frm]->ylim2[0], jmin);
		atomicMax(&pos[frm]->ylim2[1], jmax);
	} else {
		atomicMin(&pos[frm]->xlim[0], imin);
		atomicMax(&pos[frm]->xlim[1], imax);
		atomicMin(&pos[frm]->ylim[0], jmin);
		atomicMax(&pos[frm]->ylim[1], jmax);
	}
}
__device__ void dev_POSrect_gpu64(
		struct pos_t **pos,
		int src,
		double imin_dbl,
		double imax_dbl,
		double jmin_dbl,
		double jmax_dbl,
		double4 *ijminmax_overall,
		int frm)	{
	int n, imin, imax, jmin, jmax;
	n = pos[frm]->n;

	/* Update the POS region that contains the target without
	 * regard to whether or not it extends beyond the POS frame */
	atomicMin64(&ijminmax_overall[frm].w, imin_dbl);
	atomicMax64(&ijminmax_overall[frm].x, imax_dbl);
	atomicMin64(&ijminmax_overall[frm].y, jmin_dbl);
	atomicMax64(&ijminmax_overall[frm].z, jmax_dbl);

	/*  Update the subset of the POS frame that contains the target  */
	imin = (imin_dbl < INT_MIN) ? INT_MIN : (int) imin_dbl;
	imax = (imax_dbl > INT_MAX) ? INT_MAX : (int) imax_dbl;
	jmin = (jmin_dbl < INT_MIN) ? INT_MIN : (int) jmin_dbl;
	jmax = (jmax_dbl > INT_MAX) ? INT_MAX : (int) jmax_dbl;

	/* Make sure it's smaller than n */
	imin = MAX(imin,-n);
	imax = MIN(imax, n);
	jmin = MAX(jmin,-n);
	jmax = MIN(jmax, n);

	if (src) {
		atomicMin(&pos[frm]->xlim2[0], imin);
		atomicMax(&pos[frm]->xlim2[1], imax);
		atomicMin(&pos[frm]->ylim2[0], jmin);
		atomicMax(&pos[frm]->ylim2[1], jmax);
	} else {
		atomicMin(&pos[frm]->xlim[0], imin);
		atomicMax(&pos[frm]->xlim[1], imax);
		atomicMin(&pos[frm]->ylim[0], jmin);
		atomicMax(&pos[frm]->ylim[1], jmax);
	}
}
