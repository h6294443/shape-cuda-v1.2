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

__device__ int dnframes, dnviews;
__device__ double anglesave[3], omegasave[3];

__global__ void get_types_krnl(struct dat_t *ddat, unsigned char *dtype) {
	/* nsets-threaded kernel */
	int s = blockIdx.x * blockDim.x + threadIdx.x;

	if (s < ddat->nsets) {
		dtype[s] = ddat->set[s].type;
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
__global__ void add_offsets_to_euler_krnl(struct mod_t *dmod, struct dat_t *ddat, int s)
{
	/* Three thread-kernel */
	int j = blockIdx.x * blockDim.x + threadIdx.x;

	if (j <= 2) {
		anglesave[j] = dmod->spin.angle[j].val;
		omegasave[j] = dmod->spin.omega[j].val;
		dmod->spin.angle[j].val += ddat->set[s].angleoff[j].val;
	}
}
__global__ void get_nframes_krnl(struct dat_t *ddat, int s)
{
	/* Single-thread kernel */
	if (threadIdx.x == 0) {
		dnframes = ddat->set[s].desc.doppler.nframes;
		dnviews  = ddat->set[s].desc.doppler.nviews; }
}
__global__ void get_deldop_nframes_krnl(struct dat_t *ddat, int s)
{
	/* Single-thread kernel */
	if (threadIdx.x == 0) {
		dnframes = ddat->set[s].desc.deldop.nframes;
		dnviews  = ddat->set[s].desc.deldop.nviews;}
}
__global__ void get_pos_nframes_krnl(struct dat_t *ddat, int s)
{
	/* Single-thread kernel */
	if (threadIdx.x == 0) {
		dnframes = ddat->set[s].desc.poset.nframes;
		dnviews  = ddat->set[s].desc.poset.nviews; }
}
__global__ void get_lghtcrv_nframes_krnl(struct dat_t *ddat, int s)
{
	/* Single-thread kernel */
	if (threadIdx.x == 0) dnframes = ddat->set[s].desc.lghtcrv.ncalc;
}
__global__ void realize_spin_dop_krnl(struct mod_t *dmod, struct dat_t *ddat,
		struct par_t *dpar, int s, int f, int k)
{
	/* Single threaded kernel */
	int j;
	if (threadIdx.x == 0) {

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
		struct par_t *dpar, int s, int f, int k)
{
	/* Single-threaded kernel */
	int j;
	if (threadIdx.x == 0) {
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
		struct par_t *dpar, int s, int f, int k)
{
	/* Single-threaded kernel */
	int j;
	if (threadIdx.x == 0)
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
		struct par_t *dpar, int s, int i, int k)
{
	/* Single-threaded kernel */
	int j;
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
__global__ void update_spin_angle_krnl(struct mod_t *dmod)
{
	int j = blockIdx.x * blockDim.x + threadIdx.x;

	if(j < 3) {
		dmod->spin.angle[j].val = anglesave[j];
		dmod->spin.omega[j].val = omegasave[j];
	}
}
__host__ void realize_spin_cuda( struct par_t *dpar, struct mod_t *dmod, struct dat_t *ddat, int nsets)
{

	int nframes, nviews;
	int s, f, k;
	unsigned char type[nsets], *dtype;
	dim3 nsetsBLK, nsetsTHD, BLK, THD;

	gpuErrchk(cudaMalloc((void**)&dtype, sizeof(unsigned char) * nsets));

	/* Calculate launch parameters for all kernels going over all vertices */
	nsetsBLK.x = floor((maxThreadsPerBlock - 1 + nsets) / maxThreadsPerBlock);
	nsetsTHD.x = maxThreadsPerBlock; // Thread block dimensions

	/* Get and copy the type of each dataset for use on host side switch statement below*/
	get_types_krnl<<<nsetsBLK,nsetsTHD>>>(ddat, dtype);
	checkErrorAfterKernelLaunch("get_types_krnl, line 346");
	gpuErrchk(cudaMemcpy(&type, dtype, (sizeof(unsigned char) * nsets), cudaMemcpyDeviceToHost));

	/* Get the three components of the angle and spin offsets for all datasets,
	 * with any "=" states taken into account  */
	realize_angleoff_krnl<<<1,1>>>(ddat);
	checkErrorAfterKernelLaunch("realize_angleoff_krnl, line 352");

	realize_omegaoff_krnl<<<nsetsBLK,nsetsTHD>>>(ddat);
	checkErrorAfterKernelLaunch("realize_omegaoff_krnl, line 355");

	/* Note: Maybe turn the dataset loop into cudaStreams later */
	/* Determine the model spin state for each dataset in turn */

	for (s=0; s<nsets; s++) {

		/* Add this dataset's angle offsets to the model Euler angles. Later
		 * we'll add the spin offsets for each frame separately, after updating
		 * the intrinsic spin vector to each epoch. Save the original Euler
		 * angles to be restored later.          */
		/* Launch kernel do add angle offsets to Euler angles.  Three threads total */
		THD.x = 3;
		add_offsets_to_euler_krnl<<<BLK,THD>>>(dmod, ddat, s);
		checkErrorAfterKernelLaunch("add_offsets_to_euler_krnl, line 369");

		switch (type[s]) {
		case DOPPLER:
			/* Get nframes and nviews for this set */
			get_nframes_krnl<<<1,1>>>(ddat, s);
			checkErrorAfterKernelLaunch("get_nframes_krnl, line 375");
			gpuErrchk(cudaMemcpyFromSymbol(&nframes, dnframes, sizeof(nframes), 0, cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpyFromSymbol(&nviews, dnviews, sizeof(nviews), 0, cudaMemcpyDeviceToHost));

			/* Loop through every view for every (smeared) frame  */
			for (f=0; f<nframes; f++)
				for (k=0; k<nviews; k++) {
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

					realize_spin_dop_krnl<<<1,1>>>(dmod, ddat, dpar, s, f, k);
					checkErrorAfterKernelLaunch("realize_impulse_dop_krnl, line 403");
				}
			break;

		case DELAY:
			/* See "case DOPPLER" above for more extensive comments, since the
			 * Doppler and delay-Doppler procedures are identical.  */
			/* Get nframes and nviews for this set */
			get_deldop_nframes_krnl<<<1,1>>>(ddat, s);
			checkErrorAfterKernelLaunch("get_deldop_nframes_krnl, line 412");
			gpuErrchk(cudaMemcpyFromSymbol(&nframes, dnframes, sizeof(nframes), 0, cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpyFromSymbol(&nviews, dnviews, sizeof(nviews), 0, cudaMemcpyDeviceToHost));

			for (f=0; f<nframes; f++)
				for (k=0; k<nviews; k++) {
					/* Deal with spin impulses  */
					/* Get the model's intrinsic spin vector (in body coordinates)
					 * at the (light-time corrected) epoch of each view.            */
					/* Apply dataset's spin offsets (also in body coordinates)
					 * to the intrinsic spin vector of this view.                    */
					realize_spin_deldop_krnl<<<1,1>>>(dmod, ddat, dpar, s, f, k);
					checkErrorAfterKernelLaunch("realize_impulse_deldop_krnl, line 424");
				}
			break;
		case POS:
			/* See "case DOPPLER" above for more extensive comments, since the
			 * Doppler and POS procedures are identical. */
			/* Get nframes and nviews for this set */
			get_pos_nframes_krnl<<<1,1>>>(ddat, s);
			checkErrorAfterKernelLaunch("get_pos_nframes_krnl, line 432");
			gpuErrchk(cudaMemcpyFromSymbol(&nframes, dnframes, sizeof(nframes), 0, cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpyFromSymbol(&nviews, dnviews, sizeof(nviews), 0, cudaMemcpyDeviceToHost));

			for (f=0; f<nframes; f++)
				for (k=0; k<nviews; k++) {
					/* Deal with spin impulses */
					/* Get model's intrinsic spin vector (in body coordinates)
					 * at the (light-time corrected) epoch of each view. */
					/* Apply dataset's spin offsets (also in body coordinates)
					 * to the intrinsic spin vector of this view. */
					realize_spin_poset_krnl<<<1,1>>>(dmod, ddat, dpar, s, f, k);
					checkErrorAfterKernelLaunch("realize_impulse_krnl, line 444");
				}
			break;

		case LGHTCRV:
			/* See "case DOPPLER" above for more extensive comments, since the
			 * procedure for each Doppler frame is identical to the procedure
			 * for each calculated lightcurve point (except that calculated
			 * lightcurve points don't have multiple "views").	 */
			get_lghtcrv_nframes_krnl<<<1,1>>>(ddat, s);
			checkErrorAfterKernelLaunch("get_lghtcrv_nframes_krnl, line 454");
			gpuErrchk(cudaMemcpyFromSymbol(&nframes, dnframes, sizeof(nframes),
					0, cudaMemcpyDeviceToHost));
			int i, ncalc;
			ncalc = nframes;

			for (i=1; i<=ncalc; i++) {
				/* Deal with spin impulses */
				/* Get model's intrinsic spin vector (in body coordinates)
				 * at (light-time corrected) epoch of lightcurve point.*/
				/* Apply this dataset's spin offsets (also in body coordinates)
				to the intrinsic spin vector of this point. */
				realize_spin_lghtcrv_krnl<<<1,1>>>(dmod, ddat, dpar, s, i, 0); // f = i, k = 0
				checkErrorAfterKernelLaunch("realize_impulse_krnl, line 466");
			}
			break;
		default:
			bailout("realize_spin: can't handle this type yet\n");
		}

		/* Final kernel launch in realize_spin_cuda */
		THD.x = 3;
		update_spin_angle_krnl<<<BLK,THD>>>(dmod);
		checkErrorAfterKernelLaunch("update_spin_angle_krnl, line 476");
	}
	cudaFree(dtype);
}







