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
#include "head.h"
}

__device__ int rsaf_nframes, rsaf_nviews;
__device__ double rsaf_anglesave[3], rsaf_omegasave[3];

__device__ void dev_realize_impulse(struct spin_t spin, double t,
		double t_integrate[], double impulse[][3], int *n_integrate, int s, int f, int k);

__global__ void add_offsets_to_euler_af_krnl(struct mod_t *dmod, struct dat_t *ddat, int s)
{
	/* Three thread-kernel */
	int j = blockIdx.x * blockDim.x + threadIdx.x;

	if (j <= 2) {
		rsaf_anglesave[j] = dmod->spin.angle[j].val;
		rsaf_omegasave[j] = dmod->spin.omega[j].val;
		dmod->spin.angle[j].val += ddat->set[s].angleoff[j].val;
	}
}
__global__ void get_nframes_af_krnl(struct dat_t *ddat, int s)
{
	/* Single-thread kernel */
	if (threadIdx.x == 0) {
		switch(ddat->set[s].type) {
		case DELAY:
			rsaf_nframes = ddat->set[s].desc.deldop.nframes;
			rsaf_nviews  = ddat->set[s].desc.deldop.nviews;
			break;
		case DOPPLER:
			rsaf_nframes = ddat->set[s].desc.doppler.nframes;
			rsaf_nviews  = ddat->set[s].desc.doppler.nviews;
			break;
		case POS:
			rsaf_nframes = ddat->set[s].desc.poset.nframes;
			rsaf_nviews  = ddat->set[s].desc.poset.nviews;
			break;
		case LGHTCRV:
			rsaf_nframes = ddat->set[s].desc.lghtcrv.ncalc;
			break;
		}
	}
}
__global__ void realize_spin_dop_af_krnl(
		struct mod_t *dmod,
		struct dat_t *ddat,
		struct par_t *dpar,
		int s,
		int nframes,
		int nviews)
{
	/* nframes-threaded kernel */
	int j, k, frm = threadIdx.x;
	if (frm < nframes) {
		/* Loop through all available views */
		for (k=0; k<nviews; k++) {
			dev_realize_impulse(dmod->spin,
					ddat->set[s].desc.doppler.frame[frm].view[k].t,
					ddat->set[s].desc.doppler.frame[frm].t_integrate,
					ddat->set[s].desc.doppler.frame[frm].impulse,
					&ddat->set[s].desc.doppler.frame[frm].n_integrate,s,frm,k);

			dev_inteuler(dmod->spin,
					ddat->set[s].desc.doppler.frame[frm].t_integrate,
					ddat->set[s].desc.doppler.frame[frm].impulse,
					ddat->set[s].desc.doppler.frame[frm].n_integrate,
					ddat->set[s].desc.doppler.frame[frm].view[k].intspin,
					ddat->set[s].desc.doppler.frame[frm].view[k].ae,
					dmod->spin.pa, dpar->int_method, dpar->int_abstol);

			for (j=0; j<=2; j++)
				ddat->set[s].desc.doppler.frame[frm].view[k].intspin[j] += ddat->set[s].omegaoff[j].val;

			dev_cotrans2(ddat->set[s].desc.doppler.frame[frm].view[k].intspin,
					ddat->set[s].desc.doppler.frame[frm].view[k].ae,
					ddat->set[s].desc.doppler.frame[frm].view[k].intspin, -1);

			for (j=0; j<=2; j++)
				ddat->set[s].desc.doppler.frame[frm].view[k].spin[j] = ddat->set[s].desc.doppler.frame[frm].view[k].orbspin[j] +
				ddat->set[s].desc.doppler.frame[frm].view[k].intspin[j];
		}
	}
}
__global__ void realize_spin_deldop_af_krnl(
		struct mod_t *dmod,
		struct dat_t *ddat,
		struct par_t *dpar,
		int s,
		int nframes,
		int nviews)
{
	/* nframes-threaded kernel */
	int j, k, frm = threadIdx.x;
	if (frm < nframes) {

		/* Loop through views */
		for (k=0; k<nviews; k++) {
			dev_realize_impulse(dmod->spin,
					ddat->set[s].desc.deldop.frame[frm].view[k].t,
					ddat->set[s].desc.deldop.frame[frm].t_integrate,
					ddat->set[s].desc.deldop.frame[frm].impulse,
					&ddat->set[s].desc.deldop.frame[frm].n_integrate,
					s, frm, k);

			dev_inteuler(dmod->spin,
					ddat->set[s].desc.deldop.frame[frm].t_integrate,
					ddat->set[s].desc.deldop.frame[frm].impulse,
					ddat->set[s].desc.deldop.frame[frm].n_integrate,
					ddat->set[s].desc.deldop.frame[frm].view[k].intspin,
					ddat->set[s].desc.deldop.frame[frm].view[k].ae,
					dmod->spin.pa, dpar->int_method, dpar->int_abstol);

			for (j=0; j<=2; j++)
				ddat->set[s].desc.deldop.frame[frm].view[k].intspin[j] += ddat->set[s].omegaoff[j].val;

			dev_cotrans2(ddat->set[s].desc.deldop.frame[frm].view[k].intspin,
					ddat->set[s].desc.deldop.frame[frm].view[k].ae,
					ddat->set[s].desc.deldop.frame[frm].view[k].intspin, -1);

			for (j=0; j<=2; j++)
				ddat->set[s].desc.deldop.frame[frm].view[k].spin[j] = ddat->set[s].desc.deldop.frame[frm].view[k].orbspin[j] +
				ddat->set[s].desc.deldop.frame[frm].view[k].intspin[j];
		}
	}
}
__global__ void realize_spin_poset_af_krnl(
		struct mod_t *dmod,
		struct dat_t *ddat,
		struct par_t *dpar,
		int s,
		int nframes,
		int nviews)
{
	/* Single-threaded kernel */
	int j, k, frm = threadIdx.x;
	if (frm < nframes)	{
		/* Loop through views */
		for (k=0; k<nviews; k++) {
			dev_realize_impulse(dmod->spin,
					ddat->set[s].desc.poset.frame[frm].view[k].t,
					ddat->set[s].desc.poset.frame[frm].t_integrate,
					ddat->set[s].desc.poset.frame[frm].impulse,
					&ddat->set[s].desc.poset.frame[frm].n_integrate,s,frm,k);

			dev_inteuler(dmod->spin,
					ddat->set[s].desc.poset.frame[frm].t_integrate,
					ddat->set[s].desc.poset.frame[frm].impulse,
					ddat->set[s].desc.poset.frame[frm].n_integrate,
					ddat->set[s].desc.poset.frame[frm].view[k].intspin,
					ddat->set[s].desc.poset.frame[frm].view[k].ae,
					dmod->spin.pa, dpar->int_method, dpar->int_abstol);

			for (j=0; j<=2; j++)
				ddat->set[s].desc.poset.frame[frm].view[k].intspin[j] +=
						ddat->set[s].omegaoff[j].val;

			dev_cotrans2(ddat->set[s].desc.poset.frame[frm].view[k].intspin,
					ddat->set[s].desc.poset.frame[frm].view[k].ae,
					ddat->set[s].desc.poset.frame[frm].view[k].intspin, -1);

			for (j=0; j<=2; j++)
				ddat->set[s].desc.poset.frame[frm].view[k].spin[j] =
						ddat->set[s].desc.poset.frame[frm].view[k].orbspin[j] +
						ddat->set[s].desc.poset.frame[frm].view[k].intspin[j];
		}
	}
}
__global__ void realize_spin_lghtcrv_af_krnl(
		struct mod_t *dmod,
		struct dat_t *ddat,
		struct par_t *dpar,
		int s,
		int nframes)
{
	/* nframes (or ncalc)-threaded kernel */
	int j, i = threadIdx.x;
	if (i < nframes)
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
__global__ void update_spin_angle_af_krnl(struct mod_t *dmod)
{
	int j = blockIdx.x * blockDim.x + threadIdx.x;

	if(j < 3) {
		dmod->spin.angle[j].val = rsaf_anglesave[j];
		dmod->spin.omega[j].val = rsaf_omegasave[j];
	}
}
__host__ void realize_spin_cuda_af( struct par_t *dpar, struct mod_t *dmod, struct dat_t *ddat, int nsets)
{

	int s, nframes, nviews;
	unsigned char *dtype;
	dim3 nsetsBLK, nsetsTHD, BLK, THD;

	cudaCalloc((void**)&dtype, sizeof(unsigned char), nsets);

	/* Calculate launch parameters for all kernels going over all vertices */
	nsetsBLK.x = floor((maxThreadsPerBlock-1+nsets)/maxThreadsPerBlock);
	nsetsTHD.x = maxThreadsPerBlock;

	/* Get and copy the type of each dataset for use on host side switch statement below*/
	get_types_krnl<<<nsetsBLK,nsetsTHD>>>(ddat, dtype);
	checkErrorAfterKernelLaunch("get_types_krnl, called from realize_spin_cuda_af.cu");

	/* Get the three components of the angle and spin offsets for all datasets,
	 * with any "=" states taken into account  */
	realize_angleoff_krnl<<<nsetsBLK,nsetsTHD>>>(ddat);
	checkErrorAfterKernelLaunch("realize_angleoff_krnl, called from realize_spin_cuda_af.cu");

	realize_omegaoff_krnl<<<nsetsBLK,nsetsTHD>>>(ddat);
	checkErrorAfterKernelLaunch("realize_omegaoff_krnl, called from realize_spin_cuda_af.cu");

	/* Determine the model spin state for each dataset in turn */
	for (s=0; s<nsets; s++) {

		/* Add this dataset's angle offsets to the model Euler angles. Later
		 * we'll add the spin offsets for each frame separately, after updating
		 * the intrinsic spin vector to each epoch. Save the original Euler
		 * angles to be restored later.          */
		/* Launch kernel do add angle offsets to Euler angles.  Three threads total */
		/* The deviceSync call is necessary to access cudaCalloc'd memory from the host */
		THD.x = 3;
		add_offsets_to_euler_af_krnl<<<BLK,THD>>>(dmod, ddat, s);
		checkErrorAfterKernelLaunch("add_offsets_to_euler_af_krnl");
		deviceSyncAfterKernelLaunch("add_offsets_to_euler_af_krnl");

		switch (dtype[s]) {
		case DOPPLER:
			/* Get nframes and nviews for this set */
			get_nframes_af_krnl<<<1,1>>>(ddat, s);
			checkErrorAfterKernelLaunch("get_nframes_krnl, line 375");
			gpuErrchk(cudaMemcpyFromSymbol(&nframes, rsaf_nframes, sizeof(int),
					0, cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpyFromSymbol(&nviews, rsaf_nviews, sizeof(int),
					0, cudaMemcpyDeviceToHost));

			/* Loop through every view for every (smeared) frame  */

			/* Create lists of epochs and impulses, starting at initial
			 * spin epoch t0 and ending at this view's epoch t, that
			 * will be "encountered" in evolving the spin state from t0
			 * to t, with the impulses negated if we're evolving back-
			 * wards in time. These lists will be used by the inteuler
			 * routine to break up evolution of the spin state) into
			 * integrations over several smaller time intervals,
			 * punctuated by spin impulses.                        */
			/* Integrate Euler's equations to get models intrinsic spin
			 * vector at the (light-time corrected) epoch of each view.
			 * dpar->int_method tells inteuler which integration
			 * method to use. If dmod->spin.pa == 1, Euler's
			 * equations aren't used (principal-axis rotator).
			 * Input dmod->spin is the initial spin specification
			 * given in the mod file. Output is frame[f].view[k].ae,
			 * the transformation matrix from ecliptic to body
			 * coordinates at epoch frame[f].view[k].t, and
			 * frame[f].view[k].intspin, the intrinsic spin vector (in
			 * body-fixed coordinates) at this epoch.         */
			THD.x = nframes;
			realize_spin_dop_af_krnl<<<1,THD>>>(dmod, ddat, dpar, s,
					nframes, nviews);
			checkErrorAfterKernelLaunch("realize_impulse_dop_af_krnl");

			break;
		case DELAY:
			/* See "case DOPPLER" above for more extensive comments, since the
			 * Doppler and delay-Doppler procedures are identical.  */
			/* Get nframes and nviews for this set */
			/* Get nframes and nviews for this set */
			get_nframes_af_krnl<<<1,1>>>(ddat, s);
			checkErrorAfterKernelLaunch("get_nframes_krnl, line 375");
			gpuErrchk(cudaMemcpyFromSymbol(&nframes, rsaf_nframes, sizeof(int),
					0, cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpyFromSymbol(&nviews, rsaf_nviews, sizeof(int),
					0, cudaMemcpyDeviceToHost));

			/* Deal with spin impulses  */
			/* Get the model's intrinsic spin vector (in body coordinates)
			 * at the (light-time corrected) epoch of each view.            */
			/* Apply dataset's spin offsets (also in body coordinates)
			 * to the intrinsic spin vector of this view.                    */
			THD.x = nframes;
			realize_spin_deldop_af_krnl<<<1,THD>>>(dmod, ddat, dpar, s,
					nframes, nviews);
			checkErrorAfterKernelLaunch("realize_impulse_deldop_af_krnl");

			break;
		case POS:
			/* See "case DOPPLER" above for more extensive comments, since the
			 * Doppler and POS procedures are identical. */
			/* Get nframes and nviews for this set */
			get_nframes_af_krnl<<<1,1>>>(ddat, s);
			checkErrorAfterKernelLaunch("get_nframes_krnl, line 375");
			gpuErrchk(cudaMemcpyFromSymbol(&nframes, rsaf_nframes, sizeof(int),
					0, cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpyFromSymbol(&nviews, rsaf_nviews, sizeof(int),
					0, cudaMemcpyDeviceToHost));

			/* Deal with spin impulses */
			/* Get model's intrinsic spin vector (in body coordinates)
			 * at the (light-time corrected) epoch of each view. */
			/* Apply dataset's spin offsets (also in body coordinates)
			 * to the intrinsic spin vector of this view. */
			THD.x = nframes;
			realize_spin_poset_af_krnl<<<1,THD>>>(dmod, ddat, dpar, s,
					nframes, nviews);
			checkErrorAfterKernelLaunch("realize_impulse_af_krnl");

			break;

		case LGHTCRV:
			/* See "case DOPPLER" above for more extensive comments, since the
			 * procedure for each Doppler frame is identical to the procedure
			 * for each calculated lightcurve point (except that calculated
			 * lightcurve points don't have multiple "views").	 */
			get_nframes_af_krnl<<<1,1>>>(ddat, s);
			checkErrorAfterKernelLaunch("get_nframes_krnl, line 375");
			gpuErrchk(cudaMemcpyFromSymbol(&nframes, rsaf_nframes, sizeof(int),
					0, cudaMemcpyDeviceToHost));

			int ncalc;
			ncalc = nframes;

			/* Deal with spin impulses */
			/* Get model's intrinsic spin vector (in body coordinates)
			 * at (light-time corrected) epoch of lightcurve point.*/
			/* Apply this dataset's spin offsets (also in body coordinates)
				to the intrinsic spin vector of this point. */
			THD.x = ncalc;
			realize_spin_lghtcrv_af_krnl<<<1,THD>>>(dmod, ddat, dpar,
					s, nframes); // f = i, k = 0
			checkErrorAfterKernelLaunch("realize_impulse_af_krnl");

			break;
		default:
			bailout("realize_spin_cuda_af: can't handle this type yet\n");
		}

		/* Final kernel launch in realize_spin_cuda */
		THD.x = 3;
		update_spin_angle_af_krnl<<<BLK,THD>>>(dmod);
		checkErrorAfterKernelLaunch("update_spin_angle_krnl, line 476");
	}
	cudaFree(dtype);
}







