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


__global__ void realize_spin_dop_mgpu_krnl(struct mod_t *dmod, struct dat_t *ddat,
		struct par_t *dpar, int nviews, int s, int size, int oddflg)
{
	/* nfrm_half0/nfrm_half1-threaded kernel for dual-GPU operation. */
	int k, hf = blockIdx.x * blockDim.x + threadIdx.x;
	int f = 2 * hf + oddflg;
	int j;
	if (hf < size) {

		for (k=0; k<nviews; k++) {
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
}

__global__ void realize_spin_deldop_mgpu_krnl(struct mod_t *dmod, struct dat_t *ddat,
		struct par_t *dpar, int nviews, int s, int size, int oddflg)
{
	/* nfrm_half0/nfrm_half1-threaded kernel for dual-GPU operation. */
	int k, hf = blockIdx.x * blockDim.x + threadIdx.x;
	int f = 2 * hf + oddflg;
	int j;
	if (hf < size) {
		for (k=0; k<nviews; k++) {
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
}

__global__ void realize_spin_poset_mgpu_krnl(struct mod_t *dmod, struct dat_t *ddat,
		struct par_t *dpar, int nviews, int s, int size, int oddflg)
{
	/* nfrm_half0/nfrm_half1-threaded kernel for dual-GPU operation. */
	int k, hf = blockIdx.x * blockDim.x + threadIdx.x;
	int f = 2 * hf + oddflg;
	int j;
	if (hf < size)
	{
		for (k=0; k<nviews; k++) {
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
}

__global__ void realize_spin_lghtcrv_krnl(struct mod_t *dmod, struct dat_t *ddat,
		struct par_t *dpar, int s, int size, int oddflg)
{
	/* nfrm_half0/nfrm_half1-threaded kernel for multi-GPU operation. */
	int hf = blockIdx.x * blockDim.x + threadIdx.x;
	int i = 2 * hf + oddflg + 1;
	int j;

	if (hf < size)
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

__host__ void realize_spin_mgpu(
		struct par_t *dpar,
		struct mod_t *dmod,
		struct dat_t *ddat,
		unsigned char *htype,
		int *nframes,
		int *nviews,
		int nsets,
		cudaStream_t *gpu0_stream,
		cudaStream_t *gpu1_stream)
{
	int s, f, hf, nfrm_alloc, nfrm_half0, nfrm_half1;
	dim3 nsetsBLK, nsetsTHD, BLK, BLK_half0, BLK_half1, THD, THD64;
	double3 *angle_omega_save;
	THD.x = maxThreadsPerBlock;
	THD64.x = 64;

	gpuErrchk(cudaMalloc((void**)&angle_omega_save, sizeof(double3)*2));

	/* Calculate launch parameters for all kernels going over all vertices */
	nsetsBLK.x = floor((THD.x - 1 + nsets) / THD.x);

	/* The following three kernels happen on GPU0 */
	gpuErrchk(cudaSetDevice(GPU0));

	/* Get the three components of the angle and spin offsets for all datasets,
	 * with any "=" states taken into account  */
	realize_angleoff_krnl<<<1,1>>>(ddat);
	checkErrorAfterKernelLaunch("realize_angleoff_krnl (realize_spin_mgpu.cu)");

	realize_omegaoff_krnl<<<nsetsBLK,THD>>>(ddat);
	checkErrorAfterKernelLaunch("realize_omegaoff_krnl, (realize_spin_mgpu.cu");

	/* Note: Maybe turn the dataset loop into cudaStreams later */
	/* Determine the model spin state for each dataset in turn */

	for (s=0; s<nsets; s++) {
		gpuErrchk(cudaSetDevice(GPU0));
		if (htype[s]==LGHTCRV)
			nfrm_alloc = nframes[s]+1;
		else
			nfrm_alloc = nframes[s];

		nfrm_half0 = nfrm_alloc/2 + nfrm_alloc%2;
		nfrm_half1 = nfrm_alloc/2;
		BLK_half0.x = floor((THD64.x - 1 + nfrm_half0) / THD64.x);
		BLK_half1.x = floor((THD64.x - 1 + nfrm_half1) / THD64.x);
		BLK.x = floor((THD.x - 1 + nviews[s]) / THD.x);

		/* Add this dataset's angle offsets to the model Euler angles. Later
		 * we'll add the spin offsets for each frame separately, after updating
		 * the intrinsic spin vector to each epoch. Save the original Euler
		 * angles to be restored later.          */
		/* Launch kernel do add angle offsets to Euler angles.  Three threads total */
		add_offsets_to_euler_krnl<<<1,1>>>(dmod,ddat,angle_omega_save,s);
		checkErrorAfterKernelLaunch("add_offsets_to_euler_krnl (realize_spin_mgpu");

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

			gpuErrchk(cudaSetDevice(GPU0));
			realize_spin_dop_mgpu_krnl<<<BLK_half0,THD64,0,gpu0_stream[0]>>>
					(dmod, ddat, dpar, nviews[s], s, nfrm_half0, 0);

			gpuErrchk(cudaSetDevice(GPU1));
			realize_spin_dop_mgpu_krnl<<<BLK_half1,THD64,0,gpu1_stream[0]>>>
					(dmod, ddat, dpar, nviews[s], s, nfrm_half1, 1);

			checkErrorAfterKernelLaunch("realize_spin_dop_mgpu_krnl");

			break;
		case DELAY:
			/* See "case DOPPLER" above for more extensive comments, since the
			 * Doppler and delay-Doppler procedures are identical.  */
			/* Deal with spin impulses  */
			/* Get the model's intrinsic spin vector (in body coordinates)
			 * at the (light-time corrected) epoch of each view.            */
			/* Apply dataset's spin offsets (also in body coordinates)
			 * to the intrinsic spin vector of this view.                    */

			gpuErrchk(cudaSetDevice(GPU0));
			realize_spin_deldop_mgpu_krnl<<<BLK_half0,THD64,0,gpu0_stream[0]>>>
					(dmod, ddat, dpar, nviews[s], s, nfrm_half0, 0);

			gpuErrchk(cudaSetDevice(GPU1));
			realize_spin_deldop_mgpu_krnl<<<BLK_half1,THD64,0,gpu1_stream[0]>>>
					(dmod, ddat, dpar, nviews[s], s, nfrm_half1, 1);

			checkErrorAfterKernelLaunch("realize_spin_deldop_mgpu_krnl");

			break;
		case POS:
			/* See "case DOPPLER" above for more extensive comments, since the
			 * Doppler and POS procedures are identical. */
			/* Deal with spin impulses */
			/* Get model's intrinsic spin vector (in body coordinates)
			 * at the (light-time corrected) epoch of each view. */
			/* Apply dataset's spin offsets (also in body coordinates)
			 * to the intrinsic spin vector of this view. */

			gpuErrchk(cudaSetDevice(GPU0));
			realize_spin_poset_mgpu_krnl<<<BLK_half0,THD64,0,gpu0_stream[0]>>>(
					dmod, ddat, dpar, nviews[s], s, nfrm_half0, 0);

			gpuErrchk(cudaSetDevice(GPU1));
			realize_spin_poset_mgpu_krnl<<<BLK_half1,THD64,0,gpu1_stream[0]>>>(
					dmod, ddat, dpar, nviews[s], s, nfrm_half1, 0);

			checkErrorAfterKernelLaunch("realize_spin_poset_streams2_krnl");

			break;
		case LGHTCRV:
			/* See "case DOPPLER" above for more extensive comments, since the
			 * procedure for each Doppler frame is identical to the procedure
			 * for each calculated lightcurve point (except that calculated
			 * lightcurve points don't have multiple "views").	 */

			int i, ncalc;
			ncalc = nfrm_alloc;

			/* Deal with spin impulses */
			/* Get model's intrinsic spin vector (in body coordinates)
			 * at (light-time corrected) epoch of lightcurve point.*/
			/* Apply this dataset's spin offsets (also in body coordinates)
				to the intrinsic spin vector of this point. */
			gpuErrchk(cudaSetDevice(GPU0));
			realize_spin_lghtcrv_krnl<<<BLK_half0,THD64,0,gpu0_stream[0]>>>
					(dmod, ddat, dpar, s, nfrm_half0, 0);
			gpuErrchk(cudaSetDevice(GPU1));
			realize_spin_lghtcrv_krnl<<<BLK_half1,THD64,0,gpu1_stream[0]>>>
					(dmod, ddat, dpar, s, nfrm_half1, 1);
			checkErrorAfterKernelLaunch("realize_spin_lghtcrv_mgpu_krnl");

			break;
		default:
			bailout("realize_spin_cuda_mgpu: can't handle this type yet\n");
		}

		/* Final kernel launch in realize_spin_cuda */
		update_spin_angle_krnl<<<1,1>>>(dmod, angle_omega_save);
		checkErrorAfterKernelLaunch("update_spin_angle_streams_krnl");

	}
	cudaFree(angle_omega_save);
}
