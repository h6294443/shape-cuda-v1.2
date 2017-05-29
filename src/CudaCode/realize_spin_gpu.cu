/*****************************************************************************************
                                                                           realize_spin.c

Takes the initial spin state described in mod, computes the spin state at the epoch of
each data frame, and produces the various coordinate transformation matrices needed in
dat.  Also computes the total apparent spin vector at the epoch of each data frame.

Modified 2015 June 3 by CM:
    Implement smearing for the "fit" and "write" actions

Modified 2011 August 15 by CM:
    Determine which spin impulses must be applied to each frame or
        lightcurve point
    Pass the "int_abstol" parameter to the inteuler routine

Modified 2006 June 18 by CM:
    Eliminate range datasets

Modified 2005 January 20 by CM:
    For POS and range datasets, save the intrisic spin vector and total
        (intrinsic plus orbital) spin vector

Modified 2004 March 22 by CM:
    For lightcurve points, save the intrisic spin vector and total
        (intrinsic plus orbital) spin vector

Modified 2004 Feb 5 by CM:
    Implement "=" state for angle and spin offsets by creating
    routines realize_angleoff and realize_omegaoff

Modified 2003 May 4 by CM:
    Apply angle offsets to Doppler datasets, not just delay-Doppler
 *****************************************************************************************/
extern "C" {
#include "../shape/head.h"
}

__global__ void add_offsets_to_euler_krnl(struct mod_t *dmod,
		struct dat_t *ddat, double3 *angle_omega_save, int s)
{
	/* Single-threaded kernel */
	/*	angle_omega_save[0].x,y,z = original anglesave[3]
	 * 	angle_omega_save[1].x,y,z = original omegasave[3]
	 * 		 */

	if (threadIdx.x == 0) {

		angle_omega_save[0].x = dmod->spin.angle[0].val;
		angle_omega_save[0].y = dmod->spin.angle[1].val;
		angle_omega_save[0].z = dmod->spin.angle[2].val;
		angle_omega_save[1].x = dmod->spin.omega[0].val;
		angle_omega_save[1].y = dmod->spin.omega[1].val;
		angle_omega_save[1].z = dmod->spin.omega[2].val;
//		Original code:
//		anglesave[j] = dmod->spin.angle[j].val;
//		omegasave[j] = dmod->spin.omega[j].val;
		for (int j=0; j<=2; j++)
			dmod->spin.angle[j].val += ddat->set[s].angleoff[j].val;
	}
}
__global__ void add_offsets_to_euler_f_krnl(struct mod_t *dmod,
		struct dat_t *ddat, float3 *angle_omega_save, int s)
{
	/* Single-threaded kernel */
	/*	angle_omega_save[0].x,y,z = original anglesave[3]
	 * 	angle_omega_save[1].x,y,z = original omegasave[3]
	 * 		 */

	if (threadIdx.x == 0) {

		angle_omega_save[0].x = __double2float_rn(dmod->spin.angle[0].val);
		angle_omega_save[0].y = __double2float_rn(dmod->spin.angle[1].val);
		angle_omega_save[0].z = __double2float_rn(dmod->spin.angle[2].val);
		angle_omega_save[1].x = __double2float_rn(dmod->spin.omega[0].val);
		angle_omega_save[1].y = __double2float_rn(dmod->spin.omega[1].val);
		angle_omega_save[1].z = __double2float_rn(dmod->spin.omega[2].val);
//		Original code:
//		anglesave[j] = dmod->spin.angle[j].val;
//		omegasave[j] = dmod->spin.omega[j].val;
		for (int j=0; j<=2; j++)
			dmod->spin.angle[j].val += ddat->set[s].angleoff[j].val;
	}
}
__global__ void realize_spin_dop_krnl(struct mod_t *dmod, struct dat_t *ddat,
		struct par_t *dpar, int nviews, int s, int f)
{
	/* nview-threaded kernel */
	int k = blockIdx.x * blockDim.x + threadIdx.x;
	int j;
	if (k < nviews) {

		dev_realize_impulse(dmod->spin,
				ddat->set[s].desc.doppler.frame[f].view[k].t,
				ddat->set[s].desc.doppler.frame[f].t_integrate,
				ddat->set[s].desc.doppler.frame[f].impulse,
				&ddat->set[s].desc.doppler.frame[f].n_integrate,s,f,k);

		dev_inteuler(dmod->spin,
				ddat->set[s].desc.doppler.frame[f].t_integrate,
				ddat->set[s].desc.doppler.frame[f].impulse,
				ddat->set[s].desc.doppler.frame[f].n_integrate,
				ddat->set[s].desc.doppler.frame[f].view[k].intspin,
				ddat->set[s].desc.doppler.frame[f].view[k].ae,
				dmod->spin.pa, dpar->int_method, dpar->int_abstol);

		for (j=0; j<=2; j++)
			ddat->set[s].desc.doppler.frame[f].view[k].intspin[j] += ddat->set[s].omegaoff[j].val;

		dev_cotrans2(ddat->set[s].desc.doppler.frame[f].view[k].intspin,
				ddat->set[s].desc.doppler.frame[f].view[k].ae,
				ddat->set[s].desc.doppler.frame[f].view[k].intspin, -1);

		for (j=0; j<=2; j++)
			ddat->set[s].desc.doppler.frame[f].view[k].spin[j] = ddat->set[s].desc.doppler.frame[f].view[k].orbspin[j] +
			ddat->set[s].desc.doppler.frame[f].view[k].intspin[j];
	}
}
__global__ void realize_spin_dop_f_krnl(struct mod_t *dmod, struct dat_t *ddat,
		struct par_t *dpar, int nframes, int s, int k)
{
	/* nview-threaded kernel */
	int f = blockIdx.x * blockDim.x + threadIdx.x;
	int j;
	if (f < nframes) {

		dev_realize_impulse(dmod->spin,
				ddat->set[s].desc.doppler.frame[f].view[k].t,
				ddat->set[s].desc.doppler.frame[f].t_integrate,
				ddat->set[s].desc.doppler.frame[f].impulse,
				&ddat->set[s].desc.doppler.frame[f].n_integrate,s,f,k);

		dev_inteuler(dmod->spin,
				ddat->set[s].desc.doppler.frame[f].t_integrate,
				ddat->set[s].desc.doppler.frame[f].impulse,
				ddat->set[s].desc.doppler.frame[f].n_integrate,
				ddat->set[s].desc.doppler.frame[f].view[k].intspin,
				ddat->set[s].desc.doppler.frame[f].view[k].ae,
				dmod->spin.pa, dpar->int_method, dpar->int_abstol);

		for (j=0; j<=2; j++)
			ddat->set[s].desc.doppler.frame[f].view[k].intspin[j] += ddat->set[s].omegaoff[j].val;

		dev_cotrans2(ddat->set[s].desc.doppler.frame[f].view[k].intspin,
				ddat->set[s].desc.doppler.frame[f].view[k].ae,
				ddat->set[s].desc.doppler.frame[f].view[k].intspin, -1);

		for (j=0; j<=2; j++)
			ddat->set[s].desc.doppler.frame[f].view[k].spin[j] = ddat->set[s].desc.doppler.frame[f].view[k].orbspin[j] +
			ddat->set[s].desc.doppler.frame[f].view[k].intspin[j];
	}
}
__global__ void realize_spin_deldop_krnl(struct mod_t *dmod, struct dat_t *ddat,
		struct par_t *dpar, int nviews, int s, int f)
{
	/* nview-threaded kernel */
	int k = blockIdx.x * blockDim.x + threadIdx.x;
	int j;
	if (k < nviews) {
		dev_realize_impulse(dmod->spin,
				ddat->set[s].desc.deldop.frame[f].view[k].t,
				ddat->set[s].desc.deldop.frame[f].t_integrate,
				ddat->set[s].desc.deldop.frame[f].impulse,
				&ddat->set[s].desc.deldop.frame[f].n_integrate,
				s, f, k);

		dev_inteuler(dmod->spin,
				ddat->set[s].desc.deldop.frame[f].t_integrate,
				ddat->set[s].desc.deldop.frame[f].impulse,
				ddat->set[s].desc.deldop.frame[f].n_integrate,
				ddat->set[s].desc.deldop.frame[f].view[k].intspin,
				ddat->set[s].desc.deldop.frame[f].view[k].ae,
				dmod->spin.pa, dpar->int_method, dpar->int_abstol);

		for (j=0; j<=2; j++)
			ddat->set[s].desc.deldop.frame[f].view[k].intspin[j] += ddat->set[s].omegaoff[j].val;

		dev_cotrans2(ddat->set[s].desc.deldop.frame[f].view[k].intspin,
				ddat->set[s].desc.deldop.frame[f].view[k].ae,
				ddat->set[s].desc.deldop.frame[f].view[k].intspin, -1);

		for (j=0; j<=2; j++)
			ddat->set[s].desc.deldop.frame[f].view[k].spin[j] = ddat->set[s].desc.deldop.frame[f].view[k].orbspin[j] +
			ddat->set[s].desc.deldop.frame[f].view[k].intspin[j];
	}
}
__global__ void realize_spin_deldop_f_krnl(struct mod_t *dmod, struct dat_t *ddat,
		struct par_t *dpar, int nframes, int s, int k)
{
	/* nview-threaded kernel */
	int f = blockIdx.x * blockDim.x + threadIdx.x;
	int j;
	if (f < nframes) {
		dev_realize_impulse(dmod->spin,
				ddat->set[s].desc.deldop.frame[f].view[k].t,
				ddat->set[s].desc.deldop.frame[f].t_integrate,
				ddat->set[s].desc.deldop.frame[f].impulse,
				&ddat->set[s].desc.deldop.frame[f].n_integrate,
				s, f, k);

		dev_inteuler(dmod->spin,
				ddat->set[s].desc.deldop.frame[f].t_integrate,
				ddat->set[s].desc.deldop.frame[f].impulse,
				ddat->set[s].desc.deldop.frame[f].n_integrate,
				ddat->set[s].desc.deldop.frame[f].view[k].intspin,
				ddat->set[s].desc.deldop.frame[f].view[k].ae,
				dmod->spin.pa, dpar->int_method, dpar->int_abstol);

		for (j=0; j<=2; j++)
			ddat->set[s].desc.deldop.frame[f].view[k].intspin[j] += ddat->set[s].omegaoff[j].val;

		dev_cotrans2(ddat->set[s].desc.deldop.frame[f].view[k].intspin,
				ddat->set[s].desc.deldop.frame[f].view[k].ae,
				ddat->set[s].desc.deldop.frame[f].view[k].intspin, -1);

		for (j=0; j<=2; j++)
			ddat->set[s].desc.deldop.frame[f].view[k].spin[j] = ddat->set[s].desc.deldop.frame[f].view[k].orbspin[j] +
			ddat->set[s].desc.deldop.frame[f].view[k].intspin[j];
	}
}
__global__ void realize_spin_poset_krnl(struct mod_t *dmod, struct dat_t *ddat,
		struct par_t *dpar, int nviews, int s, int f)
{
	/* nview-threaded kernel */
	int k = blockIdx.x * blockDim.x + threadIdx.x;
	int j;
	if (k < nviews)
	{
		dev_realize_impulse(dmod->spin,
				ddat->set[s].desc.poset.frame[f].view[k].t,
				ddat->set[s].desc.poset.frame[f].t_integrate,
				ddat->set[s].desc.poset.frame[f].impulse,
				&ddat->set[s].desc.poset.frame[f].n_integrate,s,f,k);

		dev_inteuler(dmod->spin,
				ddat->set[s].desc.poset.frame[f].t_integrate,
				ddat->set[s].desc.poset.frame[f].impulse,
				ddat->set[s].desc.poset.frame[f].n_integrate,
				ddat->set[s].desc.poset.frame[f].view[k].intspin,
				ddat->set[s].desc.poset.frame[f].view[k].ae,
				dmod->spin.pa, dpar->int_method, dpar->int_abstol);

		for (j=0; j<=2; j++)
			ddat->set[s].desc.poset.frame[f].view[k].intspin[j] += ddat->set[s].omegaoff[j].val;

		dev_cotrans2(ddat->set[s].desc.poset.frame[f].view[k].intspin,
				ddat->set[s].desc.poset.frame[f].view[k].ae,
				ddat->set[s].desc.poset.frame[f].view[k].intspin, -1);

		for (j=0; j<=2; j++)
			ddat->set[s].desc.poset.frame[f].view[k].spin[j] = ddat->set[s].desc.poset.frame[f].view[k].orbspin[j] +
			ddat->set[s].desc.poset.frame[f].view[k].intspin[j];
	}
}
__global__ void realize_spin_poset_f_krnl(struct mod_t *dmod, struct dat_t *ddat,
		struct par_t *dpar, int nframes, int s, int k)
{
	/* nview-threaded kernel */
	int f = blockIdx.x * blockDim.x + threadIdx.x;
	int j;
	if (f <= nframes)
	{
		dev_realize_impulse(dmod->spin,
				ddat->set[s].desc.poset.frame[f].view[k].t,
				ddat->set[s].desc.poset.frame[f].t_integrate,
				ddat->set[s].desc.poset.frame[f].impulse,
				&ddat->set[s].desc.poset.frame[f].n_integrate,s,f,k);

		dev_inteuler(dmod->spin,
				ddat->set[s].desc.poset.frame[f].t_integrate,
				ddat->set[s].desc.poset.frame[f].impulse,
				ddat->set[s].desc.poset.frame[f].n_integrate,
				ddat->set[s].desc.poset.frame[f].view[k].intspin,
				ddat->set[s].desc.poset.frame[f].view[k].ae,
				dmod->spin.pa, dpar->int_method, dpar->int_abstol);

		for (j=0; j<=2; j++)
			ddat->set[s].desc.poset.frame[f].view[k].intspin[j] += ddat->set[s].omegaoff[j].val;

		dev_cotrans2(ddat->set[s].desc.poset.frame[f].view[k].intspin,
				ddat->set[s].desc.poset.frame[f].view[k].ae,
				ddat->set[s].desc.poset.frame[f].view[k].intspin, -1);

		for (j=0; j<=2; j++)
			ddat->set[s].desc.poset.frame[f].view[k].spin[j] = ddat->set[s].desc.poset.frame[f].view[k].orbspin[j] +
			ddat->set[s].desc.poset.frame[f].view[k].intspin[j];
	}
}
__global__ void realize_spin_lghtcrv_krnl(struct mod_t *dmod, struct dat_t *ddat,
		struct par_t *dpar, int s, int i)
{
	/* Single-threaded kernel */
	int j;
	i = i+1;	/* This fixes the offset for rend[i] */
	if (threadIdx.x == 0)
	{
		dev_realize_impulse(dmod->spin,
				ddat->set[s].desc.lghtcrv.x[i],
				ddat->set[s].desc.lghtcrv.rend[i].t_integrate,
				ddat->set[s].desc.lghtcrv.rend[i].impulse,
				&ddat->set[s].desc.lghtcrv.rend[i].n_integrate,
				s,i,0);	// s = s,  f = i, k = 0

		dev_inteuler(dmod->spin,
				ddat->set[s].desc.lghtcrv.rend[i].t_integrate,
				ddat->set[s].desc.lghtcrv.rend[i].impulse,
				ddat->set[s].desc.lghtcrv.rend[i].n_integrate,
				ddat->set[s].desc.lghtcrv.rend[i].intspin,
				ddat->set[s].desc.lghtcrv.rend[i].ae,
				dmod->spin.pa, dpar->int_method, dpar->int_abstol);

		for (j=0; j<=2; j++)
			ddat->set[s].desc.lghtcrv.rend[i].intspin[j] += ddat->set[s].omegaoff[j].val;

		dev_cotrans2(ddat->set[s].desc.lghtcrv.rend[i].intspin,
				ddat->set[s].desc.lghtcrv.rend[i].ae,
				ddat->set[s].desc.lghtcrv.rend[i].intspin, -1);

		for (j=0; j<=2; j++)
			ddat->set[s].desc.lghtcrv.rend[i].spin[j] = ddat->set[s].desc.lghtcrv.rend[i].orbspin[j] +
			ddat->set[s].desc.lghtcrv.rend[i].intspin[j];
	}
}
__global__ void realize_spin_lghtcrv_f_krnl(struct mod_t *dmod, struct dat_t *ddat,
		struct par_t *dpar, int s, int size)
{
	/* ncalc-threaded kernel */
	int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
	int j;
	if (i <= size)
	{
		dev_realize_impulse(dmod->spin,
				ddat->set[s].desc.lghtcrv.x[i],
				ddat->set[s].desc.lghtcrv.rend[i].t_integrate,
				ddat->set[s].desc.lghtcrv.rend[i].impulse,
				&ddat->set[s].desc.lghtcrv.rend[i].n_integrate,
				s,i,0);	// s = s,  f = i, k = 0

		dev_inteuler(dmod->spin,
				ddat->set[s].desc.lghtcrv.rend[i].t_integrate,
				ddat->set[s].desc.lghtcrv.rend[i].impulse,
				ddat->set[s].desc.lghtcrv.rend[i].n_integrate,
				ddat->set[s].desc.lghtcrv.rend[i].intspin,
				ddat->set[s].desc.lghtcrv.rend[i].ae,
				dmod->spin.pa, dpar->int_method, dpar->int_abstol);

		for (j=0; j<=2; j++)
			ddat->set[s].desc.lghtcrv.rend[i].intspin[j] += ddat->set[s].omegaoff[j].val;

		dev_cotrans2(ddat->set[s].desc.lghtcrv.rend[i].intspin,
				ddat->set[s].desc.lghtcrv.rend[i].ae,
				ddat->set[s].desc.lghtcrv.rend[i].intspin, -1);

		for (j=0; j<=2; j++)
			ddat->set[s].desc.lghtcrv.rend[i].spin[j] = ddat->set[s].desc.lghtcrv.rend[i].orbspin[j] +
			ddat->set[s].desc.lghtcrv.rend[i].intspin[j];
	}
}
/* Determine which spin impulses will be encountered in evolving the spin state
 * from initial spin epoch t0 to epoch t of a particular frame or lightcurve
 * point; then create lists of epochs and impulses, starting at t0 and ending
 * at t, with the impulses negated if we're evolving backwards in time.     */

__device__ void dev_realize_impulse(struct spin_t spin, double t,
		double t_integrate[], double impulse[][3], int *n_integrate, int s, int f, int k)
{
	int j, n;
	k = 0;
	t_integrate[k] = spin.t0;
	for (j=0; j<=2; j++)
		impulse[k][j] = 0.0;
	if (t >= spin.t0) {

		/* Integrating forward in time, so add the spin impulses  */
		for (n=0; n<spin.n_impulse; n++) {
			if (spin.t_impulse[n] > spin.t0 && spin.t_impulse[n] <= t) {
				k++;
				t_integrate[k] = spin.t_impulse[n];
				for (j=0; j<=2; j++)
					impulse[k][j] = spin.impulse[n][j].val;
			}
		}
		if (t_integrate[k] < t) {
			k++;
			t_integrate[k] = t;
			for (j=0; j<=2; j++)
				impulse[k][j] = 0.0;
		}
	} else {

		/* Integrating backwards in time, so subtract the spin impulses  */
		for (n=spin.n_impulse-1; n>=0; n--) {
			if (spin.t_impulse[n] < spin.t0 && spin.t_impulse[n] >= t) {
				k++;
				t_integrate[k] = spin.t_impulse[n];
				for (j=0; j<=2; j++)
					impulse[k][j] = -spin.impulse[n][j].val;
			}
		}
		if (t_integrate[k] > t) {
			k++;
			t_integrate[k] = t;
			for (j=0; j<=2; j++)
				impulse[k][j] = 0.0;
		}
	}
	*n_integrate = k + 1;
}
__global__ void update_spin_angle_krnl(struct mod_t *dmod,
		double3 *angle_omega_save)
{
	/* Single-threaded kernel */
	/*	angle_omega_save[0].x,y,z = original anglesave[3]
	 * 	angle_omega_save[1].x,y,z = original omegasave[3]
	 * 		 */
	if(threadIdx.x == 0) {
		dmod->spin.angle[0].val = angle_omega_save[0].x;
		dmod->spin.angle[1].val = angle_omega_save[0].y;
		dmod->spin.angle[2].val = angle_omega_save[0].z;
		dmod->spin.omega[0].val = angle_omega_save[1].x;
		dmod->spin.omega[1].val = angle_omega_save[1].y;
		dmod->spin.omega[2].val = angle_omega_save[1].z;
	}
}
__global__ void update_spin_angle_f_krnl(struct mod_t *dmod,
		float3 *angle_omega_save)
{
	/* Single-threaded kernel */
	/*	angle_omega_save[0].x,y,z = original anglesave[3]
	 * 	angle_omega_save[1].x,y,z = original omegasave[3]
	 * 		 */
	if(threadIdx.x == 0) {
		dmod->spin.angle[0].val = (double)angle_omega_save[0].x;
		dmod->spin.angle[1].val = (double)angle_omega_save[0].y;
		dmod->spin.angle[2].val = (double)angle_omega_save[0].z;
		dmod->spin.omega[0].val = (double)angle_omega_save[1].x;
		dmod->spin.omega[1].val = (double)angle_omega_save[1].y;
		dmod->spin.omega[2].val = (double)angle_omega_save[1].z;
	}
}
__global__ void realize_angleoff_krnl(struct dat_t *ddat)
{
	/* Single-threaded kernel - # of datasets nsets */
	/* Kernel implements the '=' state for each component of the angle offse */
	//int s = blockIdx.x * blockDim.x + threadIdx.x;

	int j, s_angleoff, s;

	if (threadIdx.x == 0) {
		for (j=0; j<=2; j++) {

			/* If a dataset has state '=' for component j of the angle offset, go back-
			 * wards in datafile until we reach a dataset for which component j of the
			 * angle offset has state 'f' or 'c' rather than '='.
			 *         s_angleoff is the number of the dataset we find.   */

			s_angleoff = -1;

			for (s=0; s<ddat->nsets; s++) {
				if (ddat->set[s].angleoff[j].state != '=')
					s_angleoff = s;
				else if (s_angleoff < 0)
					printf("can't use \"=\" state for the first dataset's angle offsets\n");
				else
					ddat->set[s].angleoff[j].val = ddat->set[s_angleoff].angleoff[j].val;
			}
		}
	}
}
__global__ void realize_omegaoff_krnl(struct dat_t *ddat)
{
	/* Multi-threaded kernel - # of datasets nsets */
	/* Implements the '=' state for each component of the spin offset   */
	int s = blockIdx.x * blockDim.x + threadIdx.x;
	int j, s_omegaoff;

	if (s < ddat->nsets) {
		for (j=0; j<=2; j++) {

			/* If a dataset has state = '=' for component j of the spin offset, go
			 * backwards in the datafile until we reach a dataset for which
			 * component j of the spin offset has state 'f' or 'c' rather than '='.
			 *         s_omegaoff is the number of the dataset we find.		 */

			s_omegaoff = -1;

			if (ddat->set[s].omegaoff[j].state != '=')
				s_omegaoff = s;
			else if (s_omegaoff < 0)
				printf("can't use \"=\" state for the first dataset's spin offsets\n");
			else
				ddat->set[s].omegaoff[j].val = ddat->set[s_omegaoff].omegaoff[j].val;
		}
	}
}

__host__ void realize_spin_gpu(
		struct par_t *dpar,
		struct mod_t *dmod,
		struct dat_t *ddat,
		unsigned char *htype,
		int *nframes,
		int *nviews,
		int nsets,
		cudaStream_t *rs_stream)
{
	int s, f;
	dim3 nsetsBLK, nsetsTHD, BLK, THD;
	double3 *angle_omega_save;
	THD.x = maxThreadsPerBlock;

	gpuErrchk(cudaMalloc((void**)&angle_omega_save, sizeof(double3)*2));

	/* Calculate launch parameters for all kernels going over all vertices */
	nsetsBLK.x = floor((THD.x - 1 + nsets) / THD.x);

	/* Get the three components of the angle and spin offsets for all datasets,
	 * with any "=" states taken into account  */
	realize_angleoff_krnl<<<1,1>>>(ddat);
	checkErrorAfterKernelLaunch("realize_angleoff_krnl (realize_spin_cuda_streams2.cu)");

	realize_omegaoff_krnl<<<nsetsBLK,THD>>>(ddat);
	checkErrorAfterKernelLaunch("realize_omegaoff_krnl, (realize_spin_cuda_streams2.cu");

	/* Note: Maybe turn the dataset loop into cudaStreams later */
	/* Determine the model spin state for each dataset in turn */

	for (s=0; s<nsets; s++) {

		/* Add this dataset's angle offsets to the model Euler angles. Later
		 * we'll add the spin offsets for each frame separately, after updating
		 * the intrinsic spin vector to each epoch. Save the original Euler
		 * angles to be restored later.          */
		/* Launch kernel do add angle offsets to Euler angles.  Three threads total */
		add_offsets_to_euler_krnl<<<1,1>>>(dmod,ddat,angle_omega_save,s);
		checkErrorAfterKernelLaunch("add_offsets_to_euler_streams2_krnl");

		BLK.x = floor((THD.x - 1 + nviews[s]) / THD.x);

		switch (htype[s]) {
		case DOPPLER:
			/* Create lists of epochs and impulses, starting at initial spin
			 * epoch t0 and ending at this view's epoch t, that will be
			 * "encountered" in evolving the spin state from t0 to t, with the
			 * impulses negated if we're evolving backwards in time. These
			 * lists will be used by the inteuler routine to break up evolution
			 * of the spin state) into integrations over several smaller time
			 * intervals, punctuated by spin impulses.                        */
			/* Integrate Euler's equations to get models intrinsic spin vector
			 * at the (light-time corrected) epoch of each view.
			 * dpar->int_method tells inteuler which integration method to use.
			 * If dmod->spin.pa == 1, Euler's equations aren't used (principal-
			 * axis rotator).
			 * Input dmod->spin is initial spin specification given in mod file.
			 * Output is frame[f].view[k].ae, the transformation matrix from
			 * ecliptic to body coordinates at epoch frame[f].view[k].t, and
			 * frame[f].view[k].intspin, the intrinsic spin vector (in body-
			 * fixed coordinates) at this epoch.         */
			/* Loop through every frame and launch a stream kernel with nview
			 * threads  */
			for (f=0; f<nframes[s]; f++)
				realize_spin_dop_krnl<<<BLK,THD,0,rs_stream[f]>>>(dmod,
						ddat, dpar, nviews[s], s, f);
			checkErrorAfterKernelLaunch("realize_spin_dop_streams2_krnl");

			break;
		case DELAY:
			/* See "case DOPPLER" above for more extensive comments, since the
			 * Doppler and delay-Doppler procedures are identical.  */
			/* Deal with spin impulses  */
			/* Get the model's intrinsic spin vector (in body coordinates)
			 * at the (light-time corrected) epoch of each view.            */
			/* Apply dataset's spin offsets (also in body coordinates)
			 * to the intrinsic spin vector of this view.                    */

			for (f=0; f<nframes[s]; f++)
				realize_spin_deldop_krnl<<<BLK,THD,0,rs_stream[f]>>>(
						dmod, ddat, dpar, nviews[s], s, f);
			checkErrorAfterKernelLaunch("realize_spin_deldop_streams2_krnl");

			break;
		case POS:
			/* See "case DOPPLER" above for more extensive comments, since the
			 * Doppler and POS procedures are identical. */
			/* Deal with spin impulses */
			/* Get model's intrinsic spin vector (in body coordinates)
			 * at the (light-time corrected) epoch of each view. */
			/* Apply dataset's spin offsets (also in body coordinates)
			 * to the intrinsic spin vector of this view. */

			for (f=0; f<nframes[s]; f++)
				realize_spin_poset_krnl<<<BLK,THD,0,rs_stream[f]>>>(
						dmod, ddat, dpar, nviews[s], s, f);
			checkErrorAfterKernelLaunch("realize_spin_poset_streams2_krnl");

			break;
		case LGHTCRV:
			/* See "case DOPPLER" above for more extensive comments, since the
			 * procedure for each Doppler frame is identical to the procedure
			 * for each calculated lightcurve point (except that calculated
			 * lightcurve points don't have multiple "views").	 */

			int i, ncalc;
			ncalc = nframes[s];

			/* Deal with spin impulses */
			/* Get model's intrinsic spin vector (in body coordinates)
			 * at (light-time corrected) epoch of lightcurve point.*/
			/* Apply this dataset's spin offsets (also in body coordinates)
				to the intrinsic spin vector of this point. */

			for (i=0; i<ncalc; i++)
				realize_spin_lghtcrv_krnl<<<1,1,0,rs_stream[i]>>>(
						dmod, ddat, dpar, s, i); // f = i, k = 0
			checkErrorAfterKernelLaunch("realize_spin_lghtcrv_streams_krnl");

			break;
		default:
			bailout("realize_spin_cuda_streams: can't handle this type yet\n");
		}

		/* Final kernel launch in realize_spin_cuda */
		update_spin_angle_krnl<<<1,1>>>(dmod, angle_omega_save);
		checkErrorAfterKernelLaunch("update_spin_angle_streams_krnl");

	}
	cudaFree(angle_omega_save);
}

__host__ void realize_spin_gpu_f(
		struct par_t *dpar,
		struct mod_t *dmod,
		struct dat_t *ddat,
		unsigned char *htype,
		int *nframes,
		int *nviews,
		int nsets,
		cudaStream_t *rs_stream)
{
	/* This version eliminates doubles as much as possible and it changes the
	 * grid sizing. Now there are nframes-threads and nviews-streams.  */
	int s, f, k;
	dim3 nsetsBLK, nsetsTHD, BLK, THD;
	float3 *angle_omega_save;
	THD.x = maxThreadsPerBlock;

	gpuErrchk(cudaMalloc((void**)&angle_omega_save, sizeof(float3)*2));

	/* Calculate launch parameters for all kernels going over all vertices */
	nsetsBLK.x = floor((THD.x - 1 + nsets) / THD.x);

	/* Get the three components of the angle and spin offsets for all datasets,
	 * with any "=" states taken into account  */
	realize_angleoff_krnl<<<1,1>>>(ddat);
	checkErrorAfterKernelLaunch("realize_angleoff_krnl (realize_spin_cuda_streams2.cu)");

	realize_omegaoff_krnl<<<nsetsBLK,THD>>>(ddat);
	checkErrorAfterKernelLaunch("realize_omegaoff_krnl, (realize_spin_cuda_streams2.cu");

	/* Note: Maybe turn the dataset loop into cudaStreams later */
	/* Determine the model spin state for each dataset in turn */

	for (s=0; s<nsets; s++) {

		/* Add this dataset's angle offsets to the model Euler angles. Later
		 * we'll add the spin offsets for each frame separately, after updating
		 * the intrinsic spin vector to each epoch. Save the original Euler
		 * angles to be restored later.          */
		/* Launch kernel do add angle offsets to Euler angles.  Three threads total */
		add_offsets_to_euler_f_krnl<<<1,1>>>(dmod,ddat,angle_omega_save,s);
		checkErrorAfterKernelLaunch("add_offsets_to_euler_streams2_F_krnl");

		BLK.x = floor((THD.x - 1 + nframes[s]) / THD.x);

		switch (htype[s]) {
		case DOPPLER:
			/* Create lists of epochs and impulses, starting at initial spin
			 * epoch t0 and ending at this view's epoch t, that will be
			 * "encountered" in evolving the spin state from t0 to t, with the
			 * impulses negated if we're evolving backwards in time. These
			 * lists will be used by the inteuler routine to break up evolution
			 * of the spin state) into integrations over several smaller time
			 * intervals, punctuated by spin impulses.                        */
			/* Integrate Euler's equations to get models intrinsic spin vector
			 * at the (light-time corrected) epoch of each view.
			 * dpar->int_method tells inteuler which integration method to use.
			 * If dmod->spin.pa == 1, Euler's equations aren't used (principal-
			 * axis rotator).
			 * Input dmod->spin is initial spin specification given in mod file.
			 * Output is frame[f].view[k].ae, the transformation matrix from
			 * ecliptic to body coordinates at epoch frame[f].view[k].t, and
			 * frame[f].view[k].intspin, the intrinsic spin vector (in body-
			 * fixed coordinates) at this epoch.         */
			/* Loop through every frame and launch a stream kernel with nview
			 * threads  */
			for (k=0; k<nviews[s]; k++)
				realize_spin_dop_f_krnl<<<BLK,THD,0,rs_stream[k]>>>(dmod,
						ddat, dpar, nframes[s], s, k);
			checkErrorAfterKernelLaunch("realize_spin_dop_streams_f_krnl");

			break;
		case DELAY:
			/* See "case DOPPLER" above for more extensive comments, since the
			 * Doppler and delay-Doppler procedures are identical.  */
			/* Deal with spin impulses  */
			/* Get the model's intrinsic spin vector (in body coordinates)
			 * at the (light-time corrected) epoch of each view.            */
			/* Apply dataset's spin offsets (also in body coordinates)
			 * to the intrinsic spin vector of this view.                    */

			for (k=0; k<nviews[s]; k++)
				realize_spin_deldop_f_krnl<<<BLK,THD,0,rs_stream[k]>>>(
						dmod, ddat, dpar, nframes[s], s, k);
			checkErrorAfterKernelLaunch("realize_spin_deldop_streams_f_krnl");

			break;
		case POS:
			/* See "case DOPPLER" above for more extensive comments, since the
			 * Doppler and POS procedures are identical. */
			/* Deal with spin impulses */
			/* Get model's intrinsic spin vector (in body coordinates)
			 * at the (light-time corrected) epoch of each view. */
			/* Apply dataset's spin offsets (also in body coordinates)
			 * to the intrinsic spin vector of this view. */

			for (k=0; k<nviews[s]; k++)
				realize_spin_poset_f_krnl<<<BLK,THD,0,rs_stream[k]>>>(
						dmod, ddat, dpar, nframes[s], s, k);
			checkErrorAfterKernelLaunch("realize_spin_poset_streams2_krnl");

			break;
		case LGHTCRV:
			/* See "case DOPPLER" above for more extensive comments, since the
			 * procedure for each Doppler frame is identical to the procedure
			 * for each calculated lightcurve point (except that calculated
			 * lightcurve points don't have multiple "views").	 */

			int ncalc;
			ncalc = nframes[s];

			/* Deal with spin impulses */
			/* Get model's intrinsic spin vector (in body coordinates)
			 * at (light-time corrected) epoch of lightcurve point.*/
			/* Apply this dataset's spin offsets (also in body coordinates)
				to the intrinsic spin vector of this point. */
			THD.x = maxThreadsPerBlock;
			BLK.x = floor((THD.x + ncalc)/THD.x);
			realize_spin_lghtcrv_f_krnl<<<BLK,THD>>>(
					dmod, ddat, dpar, s, ncalc); // f = i, k = 0
			checkErrorAfterKernelLaunch("realize_spin_lghtcrv_streams2_f_krnl");

			break;
		default:
			bailout("realize_spin_cuda_streams: can't handle this type yet\n");
		}

		/* Final kernel launch in realize_spin_cuda */
		update_spin_angle_f_krnl<<<1,1>>>(dmod, angle_omega_save);
		checkErrorAfterKernelLaunch("update_spin_angle_streams_f_krnl");

	}
	cudaFree(angle_omega_save);
}





