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

__device__ int _nframes, _nviews;

__global__ void add_offsets_to_euler_streams_krnl(struct mod_t *dmod,
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

		switch(ddat->set[s].type) {
		case DELAY:
			_nframes = ddat->set[s].desc.deldop.nframes;
			_nviews  = ddat->set[s].desc.deldop.nviews;
			break;
		case DOPPLER:
			_nframes = ddat->set[s].desc.doppler.nframes;
			_nviews  = ddat->set[s].desc.doppler.nviews;
			break;
		case POS:
			_nframes = ddat->set[s].desc.poset.nframes;
			_nviews  = ddat->set[s].desc.poset.nviews;
			break;
		case LGHTCRV:
			_nframes = ddat->set[s].desc.lghtcrv.ncalc;
			break;
		default:
			printf("Unrecognized data type in add_offsets_to_euler_streams_krnl");
		}
	}
}
__global__ void realize_spin_dop_streams_krnl(struct mod_t *dmod, struct dat_t *ddat,
		struct par_t *dpar, int s, int f)
{
	/* nview-threaded kernel */
	int k = blockIdx.x * blockDim.x + threadIdx.x;
	int j;
	if (k < _nviews) {

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
__global__ void realize_spin_deldop_streams_krnl(struct mod_t *dmod, struct dat_t *ddat,
		struct par_t *dpar, int s, int f)
{
	/* nview-threaded kernel */
	int k = blockIdx.x * blockDim.x + threadIdx.x;
	int j;
	if (k < _nviews) {
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
__global__ void realize_spin_poset_streams_krnl(struct mod_t *dmod, struct dat_t *ddat,
		struct par_t *dpar, int s, int f)
{
	/* nview-threaded kernel */
	int k = blockIdx.x * blockDim.x + threadIdx.x;
	int j;
	if (k < _nviews)
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
__global__ void realize_spin_lghtcrv_streams_krnl(struct mod_t *dmod, struct dat_t *ddat,
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
/* Determine which spin impulses will be encountered in evolving the spin state
 * from initial spin epoch t0 to epoch t of a particular frame or lightcurve
 * point; then create lists of epochs and impulses, starting at t0 and ending
 * at t, with the impulses negated if we're evolving backwards in time.     */
//
//__device__ void dev_realize_impulse(struct spin_t spin, double t,
//		double t_integrate[], double impulse[][3], int *n_integrate, int s, int f, int k)
//{
//	int j, n;
//	k = 0;
//	t_integrate[k] = spin.t0;
//	for (j=0; j<=2; j++)
//		impulse[k][j] = 0.0;
//	if (t >= spin.t0) {
//
//		/* Integrating forward in time, so add the spin impulses  */
//		for (n=0; n<spin.n_impulse; n++) {
//			if (spin.t_impulse[n] > spin.t0 && spin.t_impulse[n] <= t) {
//				k++;
//				t_integrate[k] = spin.t_impulse[n];
//				for (j=0; j<=2; j++)
//					impulse[k][j] = spin.impulse[n][j].val;
//			}
//		}
//		if (t_integrate[k] < t) {
//			k++;
//			t_integrate[k] = t;
//			for (j=0; j<=2; j++)
//				impulse[k][j] = 0.0;
//		}
//	} else {
//
//		/* Integrating backwards in time, so subtract the spin impulses  */
//		for (n=spin.n_impulse-1; n>=0; n--) {
//			if (spin.t_impulse[n] < spin.t0 && spin.t_impulse[n] >= t) {
//				k++;
//				t_integrate[k] = spin.t_impulse[n];
//				for (j=0; j<=2; j++)
//					impulse[k][j] = -spin.impulse[n][j].val;
//			}
//		}
//		if (t_integrate[k] > t) {
//			k++;
//			t_integrate[k] = t;
//			for (j=0; j<=2; j++)
//				impulse[k][j] = 0.0;
//		}
//	}
//	*n_integrate = k + 1;
//}
__global__ void update_spin_angle_streams_krnl(struct mod_t *dmod,
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
__host__ void realize_spin_cuda_streams( struct par_t *dpar, struct mod_t *dmod, struct dat_t *ddat, int nsets)
{
	int nframes, nviews;
	int s, f;
	unsigned char type[nsets], *dtype;
	dim3 nsetsBLK, nsetsTHD, BLK, THD;
	double3 *angle_omega_save;
	THD.x = maxThreadsPerBlock;

	gpuErrchk(cudaMalloc((void**)&dtype, sizeof(unsigned char) * nsets));
	gpuErrchk(cudaMalloc((void**)&angle_omega_save, sizeof(double3)*2));

	/* Calculate launch parameters for all kernels going over all vertices */
	nsetsBLK.x = floor((maxThreadsPerBlock - 1 + nsets) / maxThreadsPerBlock);
	nsetsTHD.x = maxThreadsPerBlock; // Thread block dimensions

	/* Get and copy the type of each dataset for use on host side switch statement below*/
	get_types_krnl<<<nsetsBLK,nsetsTHD>>>(ddat, dtype);
	checkErrorAfterKernelLaunch("get_types_krnl, line 346");
	gpuErrchk(cudaMemcpy(&type, dtype, (sizeof(unsigned char) * nsets), cudaMemcpyDeviceToHost));
	cudaFree(dtype);

	/* Get the three components of the angle and spin offsets for all datasets,
	 * with any "=" states taken into account  */
	realize_angleoff_krnl<<<1,1>>>(ddat);
	checkErrorAfterKernelLaunch("realize_angleoff_krnl (realize_spin_cuda_streams.cu)");

	realize_omegaoff_krnl<<<nsetsBLK,nsetsTHD>>>(ddat);
	checkErrorAfterKernelLaunch("realize_omegaoff_krnl, (realize_spin_cuda_streams.cu");

	/* Note: Maybe turn the dataset loop into cudaStreams later */
	/* Determine the model spin state for each dataset in turn */

	for (s=0; s<nsets; s++) {

		/* Add this dataset's angle offsets to the model Euler angles. Later
		 * we'll add the spin offsets for each frame separately, after updating
		 * the intrinsic spin vector to each epoch. Save the original Euler
		 * angles to be restored later.          */
		/* Launch kernel do add angle offsets to Euler angles.  Three threads total */
		add_offsets_to_euler_streams_krnl<<<1,1>>>(dmod,ddat,angle_omega_save,s);
		checkErrorAfterKernelLaunch("add_offsets_to_euler_streams_krnl");
		gpuErrchk(cudaMemcpyFromSymbol(&nframes, _nframes, sizeof(int),
				0, cudaMemcpyDeviceToHost));
		gpuErrchk(cudaMemcpyFromSymbol(&nviews, _nviews, sizeof(int), 0,
				cudaMemcpyDeviceToHost));

		/* Create streams and calculate the launch parameters for all nframes-
		 * threaded kernels to follow 		 */
		cudaStream_t rs_stream[nframes];
		for (f=0; f<nframes; f++)
			cudaStreamCreate(&rs_stream[f]);
		BLK.x = floor((THD.x - 1 + nviews) / THD.x);

		switch (type[s]) {
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
			for (f=0; f<nframes; f++)
				realize_spin_dop_streams_krnl<<<BLK,THD,0,rs_stream[f]>>>(dmod,
						ddat, dpar, s, f);
			checkErrorAfterKernelLaunch("realize_spin_dop_streams_krnl");

			break;
		case DELAY:
			/* See "case DOPPLER" above for more extensive comments, since the
			 * Doppler and delay-Doppler procedures are identical.  */
			/* Deal with spin impulses  */
			/* Get the model's intrinsic spin vector (in body coordinates)
			 * at the (light-time corrected) epoch of each view.            */
			/* Apply dataset's spin offsets (also in body coordinates)
			 * to the intrinsic spin vector of this view.                    */

			for (f=0; f<nframes; f++)
				realize_spin_deldop_streams_krnl<<<BLK,THD,0,rs_stream[f]>>>(
						dmod, ddat, dpar, s, f);
			checkErrorAfterKernelLaunch("realize_spin_deldop_streams_krnl");

			break;
		case POS:
			/* See "case DOPPLER" above for more extensive comments, since the
			 * Doppler and POS procedures are identical. */
			/* Deal with spin impulses */
			/* Get model's intrinsic spin vector (in body coordinates)
			 * at the (light-time corrected) epoch of each view. */
			/* Apply dataset's spin offsets (also in body coordinates)
			 * to the intrinsic spin vector of this view. */

			for (f=0; f<nframes; f++)
				realize_spin_poset_streams_krnl<<<BLK,THD,0,rs_stream[f]>>>(
						dmod, ddat, dpar, s, f);
			checkErrorAfterKernelLaunch("realize_spin_poset_streams_krnl");

			break;
		case LGHTCRV:
			/* See "case DOPPLER" above for more extensive comments, since the
			 * procedure for each Doppler frame is identical to the procedure
			 * for each calculated lightcurve point (except that calculated
			 * lightcurve points don't have multiple "views").	 */

			int i, ncalc;
			ncalc = nframes;

			/* Deal with spin impulses */
			/* Get model's intrinsic spin vector (in body coordinates)
			 * at (light-time corrected) epoch of lightcurve point.*/
			/* Apply this dataset's spin offsets (also in body coordinates)
				to the intrinsic spin vector of this point. */

			for (i=0; i<ncalc; i++)
				realize_spin_lghtcrv_streams_krnl<<<1,1,0,rs_stream[i]>>>(
						dmod, ddat, dpar, s, i); // f = i, k = 0
			checkErrorAfterKernelLaunch("realize_spin_lghtcrv_streams_krnl");

			break;
		default:
			bailout("realize_spin_cuda_streams: can't handle this type yet\n");
		}

		/* Final kernel launch in realize_spin_cuda */
		update_spin_angle_streams_krnl<<<1,1>>>(dmod, angle_omega_save);
		checkErrorAfterKernelLaunch("update_spin_angle_streams_krnl");

		/* Destroy this set's streams */
		for (f=0; f<nframes; f++)
			cudaStreamDestroy(rs_stream[f]);
	}
	cudaFree(angle_omega_save);
}







