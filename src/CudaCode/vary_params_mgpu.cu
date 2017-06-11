/*****************************************************************************************
                                                                            vary_params.c

This routine is called by every processing node for every trial value of every floating
parameter during a fit, in order to implement the "vary_radalb" "vary_optalb"
"vary_delcor0" and "vary_dopscale" parameters.  The code, which is essentially lifted from
calc_fits.c, computes up to four means:

a) mean distance towards Earth of the subradar point relative to the COM,
   for delay-Doppler frames whose 0th-order delay correction polynomial coefficient is not
   held constant; this is used to adjust the 0th-order delay correction polynomial
   coefficient if the "vary_delcor0" parameter is turned on.

b) mean "radar" projected area for (delay-)Doppler frames that are treated as absolute
   photometry; this is used to adjust the radar albedo (R) if the "vary_radalb" parameter
   is turned on.

c) mean "optical" unshadowed projected area for calculated lightcurve points that are
   treated as absolute photometry; this is used to adjust the optical albedo (R or w) if
   the "vary_optalb" parameter is turned on.  Note that plane-of-sky datasets are not used
   here, since these frames are always treated as relative photometry.

d) mean cos(subradar latitude) for (delay-)Doppler frames in datasets whose Doppler
   scaling parameter is allowed to float; this is used to adjust those parameters if the
   "vary_dopscale" parameter is turned on.

When a branch node calls this routine, it returns its datasets' summed contributions (NOT
mean contributions) to the four output parameters, deldop_zmax, rad_xsec, opt_brightness,
and cos_subradarlat.

When the root node calls this routine, it first computes its datasets' summed
contributions to these four parameters; then it receives and adds in the contributions
from the branch nodes; and finally it returns the mean (NOT summed) parameters.

Before calling vary_params, the model's size/shape and spin states must be realized
(realize_mod and realize_spin); if albedos are being varied jointly with other parameters,
the photometric state must also be realized (realize_photo); and in either case the
0th-order delay correction polynomial coefficients and the Doppler scaling factors must be
reset to their saved values via the appropriate calls to realize_delcor and
realize_dopscale, respectively.

Modified 2017 March 27 by ME:
	Split off again from previous CUDA code to create a cudaStreams version.
	cudaStreams provide another level of parallelism by executing functions
	inside a stream in that specific order, but other streams are independent
	and can thus perform their own tasks at the same time.  The CUDA runtime
	driver will keep loading the GPU with parallel streamed tasks until
	capacity is reached or the tasks run out.
	A special note on the code structure in this version:  There are a lot of
	if blocks and for loops over frames that may seem weirdly placed or
	inefficient, i.e. multiple for-loops through frames right after one another
	with just one or two lines of code inside each loop.  This is done
	deliberately to launch streamed kernels in parallel.  Alteration to the
	code could break this, resulting in a loss of parallelism and therefore:
	speed.

Modified 2016 November 6 by ME:
	Split off from vary_params to create a version that performs almost exclusively
	on the GPU

Modified 2015 June 10 by CM:
    Implement smearing

Modified 2014 February 12 by CM:
    Add "ilaw" argument to the apply_photo routine

Modified 2012 March 23 by CM:
    Implement Doppler scaling

Modified 2011 September 10 by CM:
    Two small aesthetic changes in the lightcurve section of the code

Modified 2010 June 15 by CM:
    Revise arguments to pos2deldop and pos2doppler routines

Modified 2010 April 12 by CM:
    Include overflow region when computing cross sections
    Added comment about calling realize_delcor before calling vary_params

Modified 2009 March 29 by CM:
    For MPI_Recv calls, mpi_par[0] is no longer equal to the MPI action,
        since the message tag argument already serves that purpose (as of
        2008 April 10) -- so the other mpi_par elements are renumbered
    Add "warn_badradar" argument to pos2deldop and pos2doppler routines

Modified 2008 April 10 by CM:
    Use message tag argument to MPI_Recv to identify the MPI action

Modified 2007 August 18 by CM:
    Rename MPI_TAG to MPI_TAG_1 to avoid name conflict with mpich headers

Modified 2007 August 4 by CM:
    Add orbit_offset and body arguments to posvis routine and remove
        facet argument
    Add orbit_xoff, orbit_yoff, orbit_dopoff, and body arguments to
        pos2deldop and pos2doppler routines
    Add body argument to apply_photo routine

Written 2006 October 1 by CM
 *****************************************************************************************/

extern "C" {
#include "../shape/head.h"
}
__device__ float mgpu_deldop_cross_section=0.0, mgpu_doppler_cross_section=0.0,
		mgpu_sum_rad_xsec=0.0, mgpu_sum_cos_subradarlat=0.0, mgpu_sum_deldop_zmax,
		mgpu_sum_opt_brightness;
__device__ double mgpu_deldop_zmax=0.0, mgpu_rad_xsec=0.0, mgpu_opt_brightness=0.0,
		mgpu_cos_subradarlat=0.0;


__global__ void mgpu_init_krnl(
		struct par_t *dpar,
		struct dat_t *ddat,
		int *compute_zmax,
		int *compute_cosdelta,
		int *compute_brightness,
		unsigned char *dtype,
		int nsets) {
	/* Single-threaded kernel, to be performed by GPU0 */
	int s;
	if (threadIdx.x == 0) {
		/* Initialize __device__ (file scope) variables to zero */
		mgpu_deldop_cross_section = 0.0;
		mgpu_doppler_cross_section = 0.0;
		mgpu_sum_rad_xsec = 0.0;
		mgpu_sum_cos_subradarlat = 0.0;
		mgpu_sum_deldop_zmax = 0.0;
		mgpu_sum_opt_brightness  = 0.0;
		mgpu_deldop_zmax = 0.0;
		mgpu_rad_xsec = 0.0;
		mgpu_opt_brightness = 0.0;
		mgpu_cos_subradarlat = 0.0;

		for (s=0; s<nsets; s++) {
			switch(dtype[s]) {
			case DELAY:
				compute_cosdelta[s] = (dpar->vary_dopscale != VARY_NONE
						&& ddat->set[s].desc.deldop.dopscale.state != 'c');
				compute_zmax[s] = (dpar->vary_delcor0 != VARY_NONE
						&& ddat->set[s].desc.deldop.delcor.a[0].state != 'c');
				break;
			case DOPPLER:
				compute_cosdelta[s] = (dpar->vary_dopscale != VARY_NONE &&
						ddat->set[s].desc.doppler.dopscale.state != 'c');
				break;
			case POS:
				break;
			case LGHTCRV:
				compute_brightness[s] = (dpar->vary_optalb != VARY_NONE
						&& ddat->set[s].desc.lghtcrv.cal.state == 'c');
			}
		}
	}
}

__global__ void zmax_finalize_mgpu_krnl(float value) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		mgpu_sum_deldop_zmax += value;
	}
}

__global__ void compute_xsec_doppler_mgpu_krnl(struct dat_t *ddat, float frm_xsec,
		int s, int f) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		switch (ddat->set[s].type) {
		case DELAY:
			mgpu_deldop_cross_section = __double2float_rd(ddat->set[s].desc.deldop.frame[f].overflow_xsec);
			mgpu_deldop_cross_section += frm_xsec; // fit is the end result of parallel reduction
			mgpu_deldop_cross_section *= __double2float_rd(ddat->set[s].desc.deldop.frame[f].cal.val);
			mgpu_sum_rad_xsec += mgpu_deldop_cross_section *
					__double2float_rd(ddat->set[s].desc.deldop.frame[f].weight);
			break;
		case DOPPLER:
			mgpu_doppler_cross_section = __double2float_rd(ddat->set[s].desc.doppler.frame[f].overflow_xsec);
			mgpu_doppler_cross_section += frm_xsec;
			mgpu_doppler_cross_section *= __double2float_rd(ddat->set[s].desc.doppler.frame[f].cal.val);
			mgpu_sum_rad_xsec += mgpu_doppler_cross_section *
					__double2float_rd(ddat->set[s].desc.doppler.frame[f].weight);

			break;
		}
	}
}

__global__ void xsec_deldop_mgpu_finish_krnl(float value) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0)
		mgpu_sum_rad_xsec += value;
}

__global__ void compute_cosdelta_mgpu_krnl(struct dat_t *ddat, int s, int size,
		int oddflg) {

	/* nfrm_half0/nfrm_half1-threaded kernel for dual-GPU operation */

	int hf = blockIdx.x * blockDim.x + threadIdx.x;
	int f = 2 * hf + oddflg;

	if (hf < size) {
		int j, view;
		double weight, cos_delta, oa[3][3], to_earth[3];

		switch(ddat->set[s].type){
		case DELAY:
			view = ddat->set[s].desc.deldop.v0;
			/* oa = matrix to transform body-fixed to observer coordinates  */
			/* to_earth = normalized target-to-Earth vector in body-fixed coords  */
			dev_mtrnsps(oa, ddat->set[s].desc.deldop.frame[f].view[view].ae);
			dev_mmmul(oa, ddat->set[s].desc.deldop.frame[f].view[view].oe, oa);
			for (j=0; j<=2; j++)
				to_earth[j] = oa[2][j];
			cos_delta = sqrt(to_earth[0]*to_earth[0] + to_earth[1]*to_earth[1]);
			weight = ddat->set[s].desc.deldop.frame[f].weight;
			mgpu_sum_cos_subradarlat += cos_delta*weight;
			break;
		case DOPPLER:
			view = ddat->set[s].desc.doppler.v0;
			/* oa = matrix to transform body-fixed to observer coordinates  */
			/* to_earth = normalized target-to-Earth vector in body-fixed coords  */
			dev_mtrnsps(oa, ddat->set[s].desc.doppler.frame[f].view[view].ae);
			dev_mmmul(oa, ddat->set[s].desc.doppler.frame[f].view[view].oe, oa);
			for (j=0; j<=2; j++)
				to_earth[j] = oa[2][j];
			cos_delta = sqrt(to_earth[0]*to_earth[0] + to_earth[1]*to_earth[1]);
			weight = ddat->set[s].desc.doppler.frame[f].weight;
			mgpu_sum_cos_subradarlat += cos_delta*weight;
		}
	}
}

__global__ void posclr_mgpu_krnl(struct pos_t **pos, int *posn, int f, int hf,
		int bdflag)
{
	/* Multi-threaded kernel (2*pos->n + 1)^2 threads) */
	/* The bdflag indicates whether to use the double array for brightness values */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int n = posn[hf];
	int nx = 2 * n + 1;
	int i = (offset % nx) - n;
	int j = (offset / nx) - n;

	if (offset < (nx*nx)) {
		/* For each POS pixel, zero out the optical brightness (b) and
		 * cos(scattering angle), reset the z coordinate (distance from COM towards
		 * Earth) to a dummy value, and reset the body, component, and facet onto
		 * which the pixel center projects to  dummy values                  */
		pos[hf]->body[i][j] = pos[hf]->comp[i][j] = pos[hf]->f[i][j] = -1;

		pos[hf]->b_s[offset] = pos[hf]->cose_s[offset] = pos[hf]->cosi_s[offset] = 0.0;
		pos[hf]->z_s[offset] = -HUGENUMBER;
		if (bdflag)
			pos[hf]->b_d[offset] = 0.0;

		/* In the x direction, reset the model's leftmost and rightmost
		 * pixel number to dummy values, and similarly for the y direction   */
		pos[hf]->xlim[0] = pos[hf]->ylim[0] =  n;
		pos[hf]->xlim[1] = pos[hf]->ylim[1] = -n;

		/* For a bistatic situation (lightcurve or plane-of-sky dataset), zero out
		 * cos(incidence angle) and reset the distance towards the sun, the body,
		 * component, and facet numbers as viewed from the sun, and the model's
		 * maximum projected extent as viewed from the sun to dummy values    */
		if (pos[hf]->bistatic) {
			pos[hf]->bodyill[i][j] = pos[hf]->compill[i][j] = pos[hf]->fill[i][j] = -1;

			pos[hf]->cosill_s[offset] = 0.0;
			pos[hf]->zill_s[offset] = -HUGENUMBER;

			pos[hf]->xlim2[0] = pos[hf]->ylim2[0] =  n;
			pos[hf]->xlim2[1] = pos[hf]->ylim2[1] = -n;
		}
	}
}

__global__ void set_four_parameters_krnl(struct dat_t *ddat) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		if (ddat->sum_deldop_zmax_weights > 0.0)
			mgpu_deldop_zmax = mgpu_sum_deldop_zmax / ddat->sum_deldop_zmax_weights;
		else
			mgpu_deldop_zmax = 0.0;
		if (ddat->sum_rad_xsec_weights > 0.0) {
			mgpu_rad_xsec = mgpu_sum_rad_xsec / ddat->sum_rad_xsec_weights;			}
		else
			mgpu_rad_xsec = 0.0;
		if (ddat->sum_opt_brightness_weights > 0.0)
			mgpu_opt_brightness = mgpu_sum_opt_brightness / ddat->sum_opt_brightness_weights;
		else
			mgpu_opt_brightness = 0.0;
		if (ddat->sum_cos_subradarlat_weights > 0.0)
			mgpu_cos_subradarlat = mgpu_sum_cos_subradarlat / ddat->sum_cos_subradarlat_weights;
		else
			mgpu_cos_subradarlat = 0.0;
	}
}

__global__ void delay_params_mgpu_krnl(struct par_t *dpar, struct dat_t *ddat,
		struct pos_t **pos,
		struct deldopfrm_t **frame,
		int *compute_xsec,
		int *posn,
		int *ndel,
		int *ndop,
		int s,
		int size,
		int oddflg) {
	/* nfrm-half0/nfrm_half1-threaded kernel for multi-gpu  */

	int hf = blockIdx.x * blockDim.x + threadIdx.x;
	int f = 2*hf + oddflg;

	if (hf < size) {
		compute_xsec[hf] = (dpar->vary_radalb != VARY_NONE
				&& ddat->set[s].desc.deldop.frame[f].cal.state == 'c');

		pos[hf] = &ddat->set[s].desc.deldop.frame[f].pos;
		posn[hf] = pos[hf]->n;
		ndel[hf] = ddat->set[s].desc.deldop.frame[f].ndel;
		ndop[hf] = ddat->set[s].desc.deldop.frame[f].ndop;
		frame[hf] = &ddat->set[s].desc.deldop.frame[f];
	}
}

__global__ void dop_params_mgpu_krnl(struct par_t *dpar, struct dat_t *ddat,
		struct pos_t **pos,
		struct dopfrm_t **frame,
		int *compute_xsec,
		int *posn,
		int *ndop,
		int s,
		int size,
		int oddflg) {

	/* nfrm_half0/nfrm_half1-threaded multi-GPU kernel */
	int hf = blockIdx.x * blockDim.x + threadIdx.x;
	int f = 2*hf + oddflg;

	if (hf < size) {
		compute_xsec[hf] = (dpar->vary_radalb != VARY_NONE &&
				ddat->set[s].desc.doppler.frame[f].cal.state == 'c');
		pos[hf] = &ddat->set[s].desc.doppler.frame[f].pos;
		posn[hf] = pos[hf]->n;
		ndop[hf] = ddat->set[s].desc.doppler.frame[f].ndop;
		frame[hf] = &ddat->set[s].desc.doppler.frame[f];
	}
}

__global__ void lghtcrv_params_mgpu_krnl(struct dat_t *ddat, struct pos_t **pos,
		int *posn, int *bistatic, int s, int size, int oddflg) {
	/* nframes-threaded kernel */
	int hf = blockIdx.x * blockDim.x + threadIdx.x;
	int f = 2 * hf + oddflg + 1;	/* +1 is lightcurve offset present everywhere) */
	/* This means GPU0 now handles odd frames (but still first frame in set!) */
	if (hf < size) {
		pos[hf] = &ddat->set[s].desc.lghtcrv.rend[f].pos;
		posn[hf] = pos[hf]->n;
		bistatic[hf] = pos[hf]->bistatic;
	}
}

__global__ void set_posmtrx_mgpu_krnl(struct dat_t *ddat,
		struct pos_t **pos,
		unsigned char type,
		int s,
		int f,
		int hf) {
	/* 9-threaded kernel, streamed, multi-gpu */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % 3;
	int j = offset / 3;
	int indx;

	if (offset < 9) {
		switch (type) {
		case DELAY:
			indx = ddat->set[s].desc.deldop.v0;
			pos[hf]->ae[i][j] = ddat->set[s].desc.deldop.frame[f].view[indx].ae[i][j];
			pos[hf]->oe[i][j] = ddat->set[s].desc.deldop.frame[f].view[indx].oe[i][j];
			if (i==1 && j==1)
				pos[hf]->bistatic = 0;
			break;
		case DOPPLER:
			indx = ddat->set[s].desc.doppler.v0;
			pos[hf]->ae[i][j] =	ddat->set[s].desc.doppler.frame[f].view[indx].ae[i][j];
			pos[hf]->oe[i][j] =	ddat->set[s].desc.doppler.frame[f].view[indx].oe[i][j];
			if (i==1 && j==1)
				pos[hf]->bistatic = 0;
			break;
		case LGHTCRV:
//			pos[hf]->ae[i][j] = ddat->set[s].desc.lghtcrv.rend[f].ae[i][j];
			pos[hf]->oe[i][j] = ddat->set[s].desc.lghtcrv.rend[f].oe[i][j];
			pos[hf]->se[i][j] = ddat->set[s].desc.lghtcrv.rend[f].se[i][j];
			if (i==1 && j==1)
				pos[hf]->bistatic = 1;
			break;
		}
	}
}

__global__ void get_xylim_mgpu_krnl(struct par_t *dpar, struct pos_t **pos,
		int4 *xylim, int nframes, int oddflg) {
	/* nfrm_half0/nfrm_half1-threaded kernel for dual-GPU operation. */

	int hf = blockIdx.x * blockDim.x + threadIdx.x;
	if (hf < nframes) {
		xylim[hf].w = pos[hf]->xlim[0];
		xylim[hf].x = pos[hf]->xlim[1];
		xylim[hf].y = pos[hf]->ylim[0];
		xylim[hf].z = pos[hf]->ylim[1];
	}
}


__host__ void vary_params_mgpu(
		struct par_t *dpar,
		struct mod_t *dmod,
		struct dat_t *ddat,
		int action,
		double *deldop_zmax,
		double *rad_xsec,
		double *opt_brightness,
		double *cos_subradarlat,
		int *hnframes,
		int *hlc_n,
		int *nviews,
		unsigned char *htype,
		unsigned char *dtype,
		int nf,
		int nsets,
		cudaStream_t *gpu0_stream,
		cudaStream_t *gpu1_stream,
		int max_frames)
{
	/* This iteration of vary_params combines previous new features like streams,
	 * re-used allocations, and streamed reductions with dual-GPU mode. For this
	 * to work properly, each GPU has write access only to memory on its own
	 * side. So GPU0 cannot write to memory allocated to GPU1 (but can read it).
	 * read_dat allocated even frames to GPU0 and odd frames to GPU1.
	 * Consequently, every operation that writes to even frame memory structures
	 * MUST be performed by GPU0. Conversely, ever write operation to odd frame
	 * structures MUST be performed by GPU1.
	 * Inputs:
	 * int action 				- dpar->action
	 * double *deldop_zmax 		- &deldop_zmax_save
	 * double *rad_xsec 		- &rad_xsec_save
	 * double *optbrightness	- &opt_brightness_save
	 * double *cos_subradarlat	- &cos_subradarlat_save
	 * int nsets 				- ddat->nsets
	 */

	int c=0, f, hf, s, *compute_brightness, *compute_zmax, *bistatic0, ncalc,
			*bistatic1, bistatic_all, *compute_cosdelta, n, nx,
			lghtcrv_bistatic, nfrm_alloc, nfrm_alloc_max, nThreads,
			compute_zmax_flag;
	nfrm_alloc_max = max_frames + 1;
	int nfrm_half0_max = nfrm_alloc_max/2 + nfrm_alloc_max%2;
	int nfrm_half1_max = nfrm_alloc_max/2;
	int hcomp_xsec0[nfrm_half0_max], hcomp_xsec1[nfrm_half1_max],
		npxls[nfrm_alloc_max], ddsize[nfrm_alloc_max],
		hndop0[nfrm_half0_max], hndop1[nfrm_half1_max],
		hndel0[nfrm_half0_max], hndel1[nfrm_half1_max],
		hposn0[nfrm_half0_max], hposn1[nfrm_half1_max],
		lc_xspan[nfrm_alloc_max],
		hbistatic0[nfrm_half0_max], hbistatic1[nfrm_half1_max],
		nThreadspx1[nfrm_alloc_max], hcomp_cosdelta[nsets], hcomp_zmax[nsets+1],
		hcomp_brightness[nsets+1];
	int2 span[nfrm_alloc_max];
	int4 *xylim0, *xylim1, hxylim[nfrm_alloc_max], hxylim0[nfrm_half0_max], hxylim1[nfrm_half1_max];
	float zmax, *pixels_per_km0, *pixels_per_km1, xsec[nfrm_alloc_max];
	float3 orbit_offset;
	double *u;
	double3 *so0, *so1;
	struct pos_t **pos0, **pos1;
	struct dopfrm_t **dframe0, **dframe1;
	struct deldopfrm_t **ddframe0, **ddframe1;
	int *posn0, *posn1, *ndel0, *ndel1, *ndop0, *ndop1, *compute_xsec0, *compute_xsec1,
		*outbndarr0, *outbndarr1;
	int nfrm_half0, nfrm_half1;

	cudaEvent_t start1, stop1;
	float milliseconds;

	dim3 pxBLK,THD,BLKncalc,THD9,BLKpx1, BLK_half0, BLK_half1,
		 BLK[nfrm_alloc_max], ddBLK[nfrm_alloc_max], THD64;
	THD.x = maxThreadsPerBlock;	THD9.x = 9;	THD64.x = 64;

	/* Initialize */
	orbit_offset.x = orbit_offset.y = orbit_offset.z = 0.0;

	/* GPU0 allocations */
	gpuErrchk(cudaSetDevice(GPU0));
	gpuErrchk(cudaMalloc((void**)&pos0, sizeof(struct pos_t*) * nfrm_half0_max));
	gpuErrchk(cudaMalloc((void**)&dframe0, sizeof(struct dopfrm_t*)*nfrm_half0_max));
	gpuErrchk(cudaMalloc((void**)&ddframe0, sizeof(struct deldopfrm_t*)*nfrm_half0_max));
	gpuErrchk(cudaMalloc((void**)&posn0, 	   sizeof(int) * nfrm_half0_max));
	gpuErrchk(cudaMalloc((void**)&ndel0, 	   sizeof(int) * nfrm_half0_max));
	gpuErrchk(cudaMalloc((void**)&ndop0, 	   sizeof(int) * nfrm_half0_max));
	gpuErrchk(cudaMalloc((void**)&outbndarr0,   sizeof(int) * nfrm_half0_max));
	gpuErrchk(cudaMalloc((void**)&compute_xsec0,sizeof(int) * nfrm_half0_max));
	gpuErrchk(cudaMalloc((void**)&xylim0, 	   sizeof(int4)* nfrm_half0_max));
	gpuErrchk(cudaMalloc((void**)&bistatic0, 	sizeof(int)* nfrm_half0_max));
	gpuErrchk(cudaMalloc((void**)&so0, sizeof(double3)*(nfrm_half0_max*3)));
	gpuErrchk(cudaMalloc((void**)&pixels_per_km0, sizeof(int)*nfrm_half0_max));

	/* GPU1 allocations */
	gpuErrchk(cudaSetDevice(GPU1));
	gpuErrchk(cudaMalloc((void**)&pos1, sizeof(struct pos_t*) * nfrm_half1_max));
	gpuErrchk(cudaMalloc((void**)&dframe1, sizeof(struct dopfrm_t*)*nfrm_half1_max));
	gpuErrchk(cudaMalloc((void**)&ddframe1, sizeof(struct deldopfrm_t*)*nfrm_half1_max));
	gpuErrchk(cudaMalloc((void**)&posn1, sizeof(int) * nfrm_half1_max));
	gpuErrchk(cudaMalloc((void**)&ndel1, sizeof(int) * nfrm_half1_max));
	gpuErrchk(cudaMalloc((void**)&ndop1, sizeof(int) * nfrm_half1_max));
	gpuErrchk(cudaMalloc((void**)&outbndarr1, sizeof(int) * nfrm_half1_max));
	gpuErrchk(cudaMalloc((void**)&compute_xsec1, sizeof(int) * nfrm_half1_max));
	gpuErrchk(cudaMalloc((void**)&xylim1, sizeof(int4)* nfrm_half1_max));
	gpuErrchk(cudaMalloc((void**)&bistatic1, sizeof(int)* nfrm_half1_max));
	gpuErrchk(cudaMalloc((void**)&so1, sizeof(double3)*(nfrm_half1_max*3)));
	gpuErrchk(cudaMalloc((void**)&pixels_per_km1, sizeof(int)*nfrm_half1_max));

	/* Now back to GPU0 for standard allocations */
	gpuErrchk(cudaSetDevice(GPU0));

	/* Some arrays are allocated for the maximum number of frames in
	 * any one set.  That way they are allocated once and deallocated once.
	 * They are re-used for each loop.	 */
	gpuErrchk(cudaMalloc((void**)&compute_brightness, sizeof(int) * (nsets+1)));
	gpuErrchk(cudaMalloc((void**)&compute_zmax, sizeof(int) * (nsets+1)));
	gpuErrchk(cudaMalloc((void**)&compute_cosdelta, sizeof(int) * (nsets+1)));

	gpuErrchk(cudaMalloc((void**)&u, sizeof(double) * nfrm_alloc_max));

	/* Initialize the device file-scope variables (GPU0) */
	mgpu_init_krnl<<<1,1>>>(dpar, ddat, compute_zmax, compute_cosdelta,
			compute_brightness, dtype,nsets);
	checkErrorAfterKernelLaunch("mgpu_init_krnl3");

	/* Process each dataset in turn */
	for (s=0; s<nsets; s++) {
		/* Get the allocation right as the indices for lghtcrv start with 1
		 * instead of 0 everywhere else. To not run into problems at the end
		 * or start of the array, we allocate one more than strictly needed */
		if (htype[s]==LGHTCRV)			nfrm_alloc = hnframes[s] + 1;
		else							nfrm_alloc = hnframes[s];

		nfrm_half0 = nfrm_alloc/2 + nfrm_alloc%2;
		nfrm_half1 = nfrm_alloc/2;

		/* Set up initial kernel launch parameter */
		BLK_half0 = floor ((THD64.x - 1 + nfrm_half0) / THD64.x);
		BLK_half1 = floor ((THD64.x - 1 + nfrm_half1) / THD64.x);

		switch (htype[s]) {
		case DELAY:
			/* Get computation flags */
			gpuErrchk(cudaMemcpy(&hcomp_cosdelta, compute_cosdelta,
					sizeof(int)*(nsets+1), cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpy(&hcomp_zmax, compute_zmax,
					sizeof(int)*(nsets+1), cudaMemcpyDeviceToHost));

			/* Launch kernel to get all relevant parameters. In dual-GPU mode,
			 * GPU0 gets all even frames and GPU1 gets all odd frames. Each
			 * GPU uses its own streams and memory constructs. */
			f = 0;	/* Virtual frame index for the device copies of half size */
			gpuErrchk(cudaSetDevice(GPU0));
			delay_params_mgpu_krnl<<<BLK_half0,THD64,0,gpu0_stream[0]>>>(dpar, ddat,
					pos0, ddframe0, compute_xsec0, posn0, ndel0, ndop0, s,
					nfrm_half0, 0);
			gpuErrchk(cudaSetDevice(GPU1));
			delay_params_mgpu_krnl<<<BLK_half1,THD64,0,gpu1_stream[0]>>>(dpar, ddat,
					pos1, ddframe1, compute_xsec1, posn1, ndel1, ndop1, s,
					nfrm_half1, 1);
			checkErrorAfterKernelLaunch("delay_params_mgpu_krnl");

			/* Copy GPU0 stuff to host arrays */
			gpuErrchk(cudaSetDevice(GPU0));
			gpuErrchk(cudaMemcpy(&hposn0, posn0, sizeof(int)*nfrm_half0,
					cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpy(&hcomp_xsec0, compute_xsec0, sizeof(int)*
					nfrm_half0, cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpy(&hndop0, ndop0, sizeof(int)*nfrm_half0,
					cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpy(&hndel0, ndel0, sizeof(int)*nfrm_half0,
					cudaMemcpyDeviceToHost));

			/* Copy GPU1 stuff to host arrays, then back to GPU0 */
			gpuErrchk(cudaSetDevice(GPU1));
			gpuErrchk(cudaMemcpy(&hposn1, posn1, sizeof(int)*nfrm_half1,
					cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpy(&hcomp_xsec1, compute_xsec1, sizeof(int)*
					nfrm_half1, cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpy(&hndop1, ndop1, sizeof(int)*nfrm_half1,
					cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpy(&hndel1, ndel1, sizeof(int)*nfrm_half1,
					cudaMemcpyDeviceToHost));
			gpuErrchk(cudaSetDevice(GPU0));

			/* Calculate launch parameters */
			hf = 0;
			for (f=0; f<nfrm_alloc; f++) {
				/* Even frame */
				npxls[f] = (2*hposn0[hf] + 1)*(2*hposn0[hf] + 1);
				BLK[f].x = floor((THD.x - 1 + npxls[f]) / THD.x);
				ddsize[f]= hndel0[hf] * hndop0[hf];
				ddBLK[f] = floor((THD.x -1 + ddsize[f]) / THD.x);

				/* Increase f, check for bounds, calculate odd frame */
				f++; if (f==nfrm_alloc)	break;
				npxls[f] = (2*hposn1[hf] + 1)*(2*hposn1[hf] + 1);
				BLK[f].x = floor((THD.x - 1 + npxls[f]) / THD.x);
				ddsize[f]= hndel1[hf] * hndop1[hf];
				ddBLK[f] = floor((THD.x -1 + ddsize[f]) / THD.x);
				hf++;
			}

			/* Launch posclr_streams_krnl to initialize POS view */
			f = 0;		/* Even frames first */
			for (hf=0; hf<nfrm_half0; hf++) {
				if (hcomp_zmax[s] || hcomp_xsec0[hf]) {
					set_posmtrx_mgpu_krnl<<<1,THD9,0,gpu0_stream[hf]>>>
							(ddat, pos0, htype[s], s, f, hf);
					posclr_mgpu_krnl<<<BLK[f],THD, 0, gpu0_stream[hf]>>>
							(pos0, posn0, f, hf, 0);
				}
				f += 2;	if (f >= nfrm_alloc)	break;
			}

			/* Now the odd frames */
			f = 1;
			gpuErrchk(cudaSetDevice(GPU1));
			//for (int f=1; f<nfrm_alloc; f+=2) {
			for (hf=0; hf<nfrm_half1; hf++) {
				if (hcomp_zmax[s] || hcomp_xsec1[hf]) {
					set_posmtrx_mgpu_krnl<<<1,THD9,0,gpu1_stream[hf]>>>
							(ddat, pos1, htype[s], s, f, hf);
					posclr_mgpu_krnl<<<BLK[f],THD, 0, gpu1_stream[hf]>>>
							(pos1, posn1, f, hf, 0);
				}
				f += 2;	if (f >= nfrm_alloc)	break;
			} checkErrorAfterKernelLaunch("set_posmtrx_mgpu_krnl and "
					"posclr_mgpu_krnl (Delay-Doppler)");

			/* Determine which POS pixels cover the target, and get distance
			 * toward Earth of each POS pixel. Pass the frame streams, too. */
			posvis_mgpu(dpar, dmod, ddat, pos0, pos1, orbit_offset,
					hposn0, hposn1, outbndarr0, outbndarr1, s, hnframes[s],
					0, nf, 0, c, htype[s], gpu0_stream, gpu1_stream);

			/* Clear out fit arrays, even frames first */
			hf = 0;
			gpuErrchk(cudaSetDevice(GPU0));
			for (f=0; f<nfrm_alloc; f++) {
				if (hcomp_xsec0[hf]) {
					clrvect_krnl1<<<ddBLK[f],THD, 0, gpu0_stream[hf]>>>(ddat,
							ddsize[f], s, f);
				}
				hf++;	f++;	if (f>=nfrm_alloc)	break;
			}
			/* Now the odd frames */
			hf = 0;
			gpuErrchk(cudaSetDevice(GPU1));
			for (f=1; f<nfrm_alloc; f++) {
				/* Odd frame */
				if (hcomp_xsec1[hf]) {
					clrvect_krnl2<<<ddBLK[f],THD, 0, gpu1_stream[hf]>>>(ddat,
							ddsize[f], s, f);
				}
				hf++;	f++;	if (f>=nfrm_alloc)	break;
			} checkErrorAfterKernelLaunch("clrvect_mgpu_krnl");
			/* Call the CUDA pos2deldop function */
			pos2deldop_mgpu(dpar, dmod, ddat, pos0, pos1, ddframe0, ddframe1,
					xylim0, xylim1, ndel0, ndel1, ndop0, ndop1, 0.0, 0.0, 0.0,
					0, s, hnframes[s], 0,	outbndarr0, outbndarr1,
					gpu0_stream, gpu1_stream);

			/* Debug start */
//			dbg_print_deldop_fit(ddat, s, 0, "GPU_deldop_fit_f0.csv");
//			dbg_print_deldop_fit(ddat, s, 1, "GPU_deldop_fit_f1.csv");
//			dbg_print_deldop_fit(ddat, s, 2, "GPU_deldop_fit_f2.csv");
//			dbg_print_deldop_fit(ddat, s, 3, "GPU_deldop_fit_f3.csv");
			/* Debug end */

			/* Calculate zmax for all frames (assumption: all pos in this set
			 * have the same pixel dimensions).  The dual-GPU version of this
			 * routine splits everything up into odd and even frames, but the
			 * specifics are contained within the compute_zmax_mgpu routine.  */
			gpuErrchk(cudaSetDevice(GPU0));
			if (hcomp_zmax[s]) {
				zmax = compute_zmax_mgpu(ddat, pos0, pos1, hnframes[s], nfrm_half0,
						nfrm_half1, npxls[0], s, gpu0_stream, gpu1_stream);
				zmax_finalize_mgpu_krnl<<<1,1>>>(zmax);
				checkErrorAfterKernelLaunch("zmax_finalize_krnl");
			}

			if (TIMING) {
				/* Create the timer events */
				cudaEventCreate(&start1);
				cudaEventCreate(&stop1);
				cudaEventRecord(start1);
			}
			/* Calculate radar cross section for each frame in set */
			xsec[0] = compute_deldop_xsec_mgpu(ddat, hnframes[s], nfrm_half0,
					nfrm_half1, ddsize[0], s, gpu0_stream, gpu1_stream);
			xsec_deldop_mgpu_finish_krnl<<<1,1>>>(xsec[0]);
			checkErrorAfterKernelLaunch("compute_xsec_final_streams2_krnl");

			if (TIMING) {
				cudaEventRecord(stop1);
				cudaEventSynchronize(stop1);
				milliseconds = 0;
				cudaEventElapsedTime(&milliseconds, start1, stop1);
				printf("Deldop xsec_streams: %3.3f ms with %i frames.\n", milliseconds, hnframes[s]);
				cudaEventDestroy(start1);
				cudaEventDestroy(stop1);
			}

			if (hcomp_cosdelta[s]) {
				for (hf=0; hf<nfrm_half0; hf++)
					compute_cosdelta_mgpu_krnl<<<BLK_half0, THD64, 0, gpu0_stream[hf]>>>
						(ddat, s, nfrm_half0, 0);
				gpuErrchk(cudaSetDevice(GPU1));
				for (hf=0; hf<nfrm_half1; hf++)
					compute_cosdelta_mgpu_krnl<<<BLK_half1, THD64, 0, gpu1_stream[hf]>>>
						(ddat, s, nfrm_half1, 1);
			}
			checkErrorAfterKernelLaunch("compute_cosdelta_streams_krnl");
			break;
		case DOPPLER:
			/* Get computation flags */
			gpuErrchk(cudaSetDevice(GPU0));
			gpuErrchk(cudaMemcpy(&hcomp_cosdelta, compute_cosdelta,
					sizeof(int)*(nsets+1), cudaMemcpyDeviceToHost));

			/* Launch kernel to get all relevant parameters. In dual-GPU mode,
			 * GPU0 gets all even frames and GPU1 gets all odd frames. Each
			 * GPU uses its own streams and memory constructs. */
			/* Even frames first */
			dop_params_mgpu_krnl<<<BLK_half0,THD64,0,gpu0_stream[0]>>>(dpar, ddat, pos0,
					dframe0, compute_xsec0, posn0, ndop0, s,nfrm_half0, 0);
			gpuErrchk(cudaSetDevice(GPU1));
			dop_params_mgpu_krnl<<<BLK_half1,THD64,0,gpu1_stream[0]>>>(dpar, ddat, pos1,
					dframe1, compute_xsec1, posn1, ndop1, s,nfrm_half1, 1);
			checkErrorAfterKernelLaunch("dop_params_mgpu_krnl");

			/* Copy GPU0 stuff to host arrays */
			gpuErrchk(cudaSetDevice(GPU0));
			gpuErrchk(cudaMemcpy(&hposn0, posn0, sizeof(int)*nfrm_half0,
					cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpyAsync(&hcomp_xsec0, compute_xsec0, sizeof(int)*
					nfrm_half0, cudaMemcpyDeviceToHost, gpu0_stream[1]));
			gpuErrchk(cudaMemcpyAsync(&hndop0, ndop0, sizeof(int)*nfrm_half0,
					cudaMemcpyDeviceToHost, gpu0_stream[2]));

			/* Copy GPU1 stuff to host arrays, then back to GPU0 */
			gpuErrchk(cudaSetDevice(GPU1));
			gpuErrchk(cudaMemcpyAsync(&hposn1, posn1, sizeof(int)*nfrm_half1,
					cudaMemcpyDeviceToHost, gpu0_stream[0]));
			gpuErrchk(cudaMemcpyAsync(&hcomp_xsec1, compute_xsec1, sizeof(int)*
					nfrm_half1, cudaMemcpyDeviceToHost, gpu0_stream[1]));
			gpuErrchk(cudaMemcpyAsync(&hndop1, ndop1, sizeof(int)*nfrm_half1,
					cudaMemcpyDeviceToHost, gpu0_stream[2]));

			/* Calculate launch parameters */
			hf = 0;
			for (f=0; f<nfrm_alloc; f++) {
				/* Even frame */
				npxls[f] = (2*hposn0[hf] + 1)*(2*hposn0[hf] + 1);
				BLK[f].x = floor((THD.x - 1 + npxls[f]) / THD.x);
				ddBLK[f] = floor((THD.x -1 + hndop0[hf]) / THD.x);

				/* Increase f, check for bounds, calculate odd frame */
				f++; if (f==nfrm_alloc)	break;
				npxls[f] = (2*hposn1[hf] + 1)*(2*hposn1[hf] + 1);
				BLK[f].x = floor((THD.x - 1 + npxls[f]) / THD.x);
				ddBLK[f] = floor((THD.x -1 + hndop1[hf]) / THD.x);
				hf++;
			}
			/* Launch posclr_streams_krnl to initialize POS view,
			 * even frames first */
			f = 0;
			gpuErrchk(cudaSetDevice(GPU0));
			for (hf=0; hf<nfrm_half0; hf++) {
				if (hcomp_xsec0[hf]) {
					set_posmtrx_mgpu_krnl<<<1,THD9,0,gpu0_stream[hf]>>>
							(ddat, pos0, htype[s], s, f, hf);
					posclr_mgpu_krnl<<<BLK[f],THD, 0, gpu0_stream[hf]>>>
							(pos0, posn0, f, hf, 0);
				}
				f += 2;	if (f >= nfrm_alloc)	break;
			}
			/* Now odd frames */
			f = 1;
			gpuErrchk(cudaSetDevice(GPU1));
			for (hf=0; hf<nfrm_half1; hf++) {
				if (hcomp_xsec1[hf]) {
					set_posmtrx_mgpu_krnl<<<1,THD9,0,gpu1_stream[hf]>>>
							(ddat, pos1, htype[s], s, f, hf);
					posclr_mgpu_krnl<<<BLK[f],THD, 0, gpu1_stream[hf]>>>
							(pos1, posn1, f, hf, 0);
				}
				f += 2;	if (f >= nfrm_alloc)	break;
			} checkErrorAfterKernelLaunch("set_posmtrx_mgpu_krnl and "
					"posclr_mgpu_krnl (Doppler)");

			/* Determine which POS pixels cover the target, and get distance
			 * toward Earth of each POS pixel. Pass the frame streams, too. */
			posvis_mgpu(dpar, dmod, ddat, pos0, pos1, orbit_offset,
					hposn0, hposn1, outbndarr0, outbndarr1, s, hnframes[s],
					0, nf, 0, c, htype[s], gpu0_stream, gpu1_stream);

			/* Clear out fit arrays, even frames first */
			hf = 0;
			gpuErrchk(cudaSetDevice(GPU0));
			for (f=0; f<nfrm_alloc; f++) {
				if (hcomp_xsec0[hf]) {
					clrvect_krnl3<<<ddBLK[f],THD, 0, gpu0_stream[hf]>>>(ddat,
							hndop0[hf], s, f);
				}
				hf++;	f++;	if (f==nfrm_alloc)	break;
			}
			/* Now the odd frames */
			hf = 0;
			gpuErrchk(cudaSetDevice(GPU1));
			for (f=1; f<nfrm_alloc; f++) {
				/* Odd frame */
				if (hcomp_xsec1[hf]) {
					clrvect_krnl4<<<ddBLK[f],THD, 0, gpu1_stream[hf]>>>(ddat,
							hndop1[hf], s, f);
				}
				hf++;	f++;	if (f==nfrm_alloc)	break;
			} checkErrorAfterKernelLaunch("clrvect_mgpu_krnl");

			pos2doppler_mgpu(dpar, dmod, ddat, pos0, pos1, dframe0, dframe1,
					xylim0, xylim1, 0.0, 0.0, 0.0, ndop0, ndop1, 0, s,
					hnframes[s], 0, outbndarr0, outbndarr1,
					gpu0_stream, gpu1_stream);

			/* Calculate the Doppler cross-section if applicable */
			gpuErrchk(cudaSetDevice(GPU0));
			for (hf=0; hf<nfrm_half0; hf++) {
				int f = 2*hf;
				if (hcomp_xsec0[hf]) {
					/* Compute cross section */
					xsec[f]=0.0;
					xsec[f] = compute_doppler_xsec(ddat, hndop0[hf], s, f);
				}
			}
			gpuErrchk(cudaSetDevice(GPU1));
			for (hf=0; hf<nfrm_half1; hf++) {
				int f = 2*hf + 1;
				if (hcomp_xsec1[hf]) {
					/* Compute cross section */
					xsec[f]=0.0;
					xsec[f] = compute_doppler_xsec(ddat, hndop1[hf], s, f);
				}
			}

			/* Finalize the xsec calculations and calculate cosdelta if specified */
			gpuErrchk(cudaSetDevice(GPU0));
			for (f=0; f<nfrm_alloc; f++)
				if (hcomp_xsec0[f/2+f%2]||hcomp_xsec1[f/2])
					compute_xsec_doppler_mgpu_krnl<<<1,1,0,gpu0_stream[f]>>>(ddat, xsec[f], s, f);

			if (hcomp_cosdelta[s]) {
				for (hf=0; hf<nfrm_half0; hf++)
					compute_cosdelta_mgpu_krnl<<<BLK_half0, THD64, 0, gpu0_stream[hf]>>>
					(ddat, s, nfrm_half0, 0);
				gpuErrchk(cudaSetDevice(GPU1));
				for (hf=0; hf<nfrm_half1; hf++)
					compute_cosdelta_mgpu_krnl<<<BLK_half1, THD64, 0, gpu1_stream[hf]>>>
					(ddat, s, nfrm_half1, 1);
			}
			gpuErrchk(cudaSetDevice(GPU0));
			break;
		case POS:
			break;
		case LGHTCRV:
			/* Lightcurves start with frame/pos index 1 instead of 0 (like
			 * radar!). For dual-GPU operation, GPU0 will handle odd frames
			 * and GPU1 will handle even frames. This is a noted departure
			 * from the opposite in radar frames.
			 * This route allows the same frame calculations as before, simply
			 * by adding +1 to all f figures.
			 * Another way to think of it:  GPU0 handles the first frame in any
			 * lightcurve set and GPU1 handles the second.
			 *
			 * 			nfrm_half0 >= nfrm_half1			 */

			gpuErrchk(cudaSetDevice(GPU0));
			/* Figure out the compute_brightness flag first */
			gpuErrchk(cudaMemcpy(&hcomp_brightness, compute_brightness,
					sizeof(int)*(nsets+1), cudaMemcpyDeviceToHost));

			if (hcomp_brightness[s]) {

				/* Launch nframes-threaded kernel to get all relevant
				 * parameters. GPU0 goes first, then GPU1 */

				lghtcrv_params_mgpu_krnl<<<BLK_half0,THD64, 0, gpu0_stream[0]>>>
						(ddat, pos0, posn0, bistatic0, s, nfrm_half0, 0);
				gpuErrchk(cudaSetDevice(GPU1));
				lghtcrv_params_mgpu_krnl<<<BLK_half1,THD64, 0, gpu1_stream[0]>>>
						(ddat, pos1, posn1, bistatic1, s, nfrm_half1, 1);
				checkErrorAfterKernelLaunch("lghtcrv_params_mgpu_krnl");

				/* Copy information back from both GPUs so I can combine it */
				gpuErrchk(cudaMemcpy(&hposn1, posn1, sizeof(int)*nfrm_half1,
						cudaMemcpyDeviceToHost));
				gpuErrchk(cudaMemcpy(&hbistatic1, bistatic1, sizeof(int)*
						nfrm_half1, cudaMemcpyDeviceToHost));
				gpuErrchk(cudaSetDevice(GPU0));
				gpuErrchk(cudaMemcpy(&hposn0, posn0, sizeof(int)*nfrm_half0,
						cudaMemcpyDeviceToHost));
				gpuErrchk(cudaMemcpy(&hbistatic0, bistatic0, sizeof(int)*
						nfrm_half0, cudaMemcpyDeviceToHost));

				/* Calculate launch parameters from the half-arrays just copied */
				hf = 0;
				for (f=1; f<nfrm_alloc; f++) {
					lc_xspan[f] = 2*hposn0[hf] + 1;
					npxls[f] = (2*hposn0[hf]+1)*(2*hposn0[hf]+1);
					BLK[f].x = floor((THD.x - 1 + npxls[f]) / THD.x);

					/* Increase f, check for bounds, calculate even frames */
					f++; if (f>=nfrm_alloc)	break;
					lc_xspan[f] = 2*hposn1[hf] + 1;
					npxls[f] = (2*hposn1[hf]+1)*(2*hposn1[hf]+1);
					BLK[f].x = floor((THD.x - 1 + npxls[f]) / THD.x);
					hf++;
				}

				/* Clear POS view, odd (GPU0) frames first */
				f = 1;
				gpuErrchk(cudaSetDevice(GPU0));
				for (hf=0; hf<nfrm_half0; hf++) {
					if (hcomp_xsec0[hf]) {
						set_posmtrx_mgpu_krnl<<<1,THD9,0,gpu0_stream[hf]>>>
								(ddat, pos0, htype[s], s, f, hf);
						posclr_mgpu_krnl<<<BLK[f],THD, 0, gpu0_stream[hf]>>>
								(pos0, posn0, f, hf, 0);
					}
					f += 2;	if (f >= nfrm_alloc)	break;
				}
				/* Now odd frames */
				f = 1;
				gpuErrchk(cudaSetDevice(GPU1));
				for (hf=0; hf<nfrm_half1; hf++) {
					if (hcomp_xsec1[hf]) {
						set_posmtrx_mgpu_krnl<<<1,THD9,0,gpu1_stream[hf]>>>
								(ddat, pos1, htype[s], s, f, hf);
						posclr_mgpu_krnl<<<BLK[f],THD, 0, gpu1_stream[hf]>>>
								(pos1, posn1, f, hf, 0);
					}
					f += 2;	if (f >= nfrm_alloc)	break;
				} checkErrorAfterKernelLaunch("set_posmtrx_mgpu_krnl and "
						"posclr_mgpu_krnl (Lightcurve)");

				/* Determine which POS pixels cover the target */
				gpuErrchk(cudaSetDevice(GPU0));
				posvis_mgpu(dpar, dmod, ddat, pos0, pos1, orbit_offset,
						hposn0, hposn1, outbndarr0, outbndarr1, s, nfrm_alloc,
						0, nf, 0, c, htype[s], gpu0_stream, gpu1_stream);

				/* Now view the model from the source (sun) and get the facet
				 * number and distance toward the source of each pixel in this
				 * projected view; use this information to determine which POS
				 * pixels are shadowed */
				f = 1;
				for (hf=0; hf<nfrm_half1; hf++)
					if (hbistatic0[hf]||hbistatic1[hf])
						bistatic_all = 1;
				if ((nfrm_alloc/2%2)==0)
					if (hbistatic0[nfrm_half0-1])
						bistatic_all = 1;

				if (bistatic_all)
					posvis_mgpu(dpar, dmod, ddat, pos0, pos1, orbit_offset,
							hposn0, hposn1, outbndarr0, outbndarr1, s, nfrm_alloc,
							1, nf, 0, c, htype[s], gpu0_stream, gpu1_stream);

				if (bistatic_all) {
					/* Initialize this stream for the posmask kernel to follow */
					posmask_init_mgpu_krnl<<<BLK_half0,THD64,0,gpu0_stream[0]>>>
							(pos0, so0, pixels_per_km0, nfrm_half0, 0);
					gpuErrchk(cudaSetDevice(GPU1));
					posmask_init_mgpu_krnl<<<BLK_half1,THD64,0,gpu1_stream[0]>>>
							(pos1, so1, pixels_per_km1, nfrm_half1, 1);
					checkErrorAfterKernelLaunch("posmask_init_mgpu_krnl");
					gpuErrchk(cudaSetDevice(GPU0));

					/* Now call posmask kernel for this stream, then loop to
					 * next stream and repeat 	 */
					for (hf=0; hf<nfrm_half0; hf++) {
						f = 2*hf+1;	/* Allows for lghtcrv offset */
						posmask_krnl<<<BLK[f],THD,0,gpu0_stream[hf]>>>(
								dpar, pos0, so0, pixels_per_km0, posn0,
								npxls[f], lc_xspan[f], hf);
					} checkErrorAfterKernelLaunch("posmask_krnl for GPU0");

					gpuErrchk(cudaSetDevice(GPU1));
					for (hf=0; hf<nfrm_half1; hf++) {
						f = 2*hf+2;	/* Allows for lghtcrv offset */
						posmask_krnl<<<BLK[f],THD,0,gpu1_stream[hf]>>>(
								dpar, pos1, so1, pixels_per_km1, posn1,
								npxls[f], lc_xspan[f], hf);
					} checkErrorAfterKernelLaunch("posmask_krnl for GPU1");
				}

				//BLKpx1.x = floor((THD.x - 1 + hnframes[s])/THD.x);
				gpuErrchk(cudaSetDevice(GPU0));
				get_xylim_mgpu_krnl<<<BLK_half0,THD64,0,gpu0_stream[0]>>>
						(dpar, pos0, xylim0, nfrm_half0, 0);
				gpuErrchk(cudaMemcpyAsync(&hxylim0, xylim0, sizeof(int4)*
						(nfrm_half0), cudaMemcpyDeviceToHost, gpu0_stream[0]));
				gpuErrchk(cudaSetDevice(GPU1));
				get_xylim_mgpu_krnl<<<BLK_half1,THD64,0,gpu1_stream[0]>>>
						(dpar, pos1, xylim1, nfrm_half1, 1);
				gpuErrchk(cudaMemcpyAsync(&hxylim1, xylim1, sizeof(int4)*
						(nfrm_half1), cudaMemcpyDeviceToHost, gpu0_stream[1]));
				checkErrorAfterKernelLaunch("get_xylim_mgpu_krnl");
				gpuErrchk(cudaSetDevice(GPU0));

				/* Calculate launch parameters for all frames */
				hf = 0;
				for (f=1; f<=hnframes[s]; f++) {
					span[f].x = hxylim0[hf].x - hxylim0[hf].w + 1;
					span[f].y = hxylim0[hf].z - hxylim0[hf].y + 1;
					nThreadspx1[f] = span[f].x * span[f].y;
					BLK[f].x = floor ((THD.x -1 + nThreadspx1[f]) / THD.x);

					/* Next up even frames */
					f++;	if (f>hnframes[s])	break;
					span[f].x = hxylim1[hf].x - hxylim1[hf].w + 1;
					span[f].y = hxylim1[hf].z - hxylim1[hf].y + 1;
					nThreadspx1[f] = span[f].x * span[f].y;
					BLK[f].x = floor ((THD.x -1 + nThreadspx1[f]) / THD.x);
				}

				/* Compute model brightness for this lightcurve point */
				/* lghtcrv->y[ncalc]: calculated points for interpolation,
				 * ncalc-points total 					 */
				apply_photo_mgpu_f(dmod, ddat, pos0, pos1, xylim0, xylim1, span,
						BLK, nThreadspx1, 0, s, nfrm_alloc, nfrm_half0,
						nfrm_half1, npxls, gpu0_stream, gpu1_stream);
//				apply_photo_cuda_streams_f(dmod, ddat, pos, xylim, span, BLK, nThreadspx1,
//							0, s, nfrm_alloc, npxls, gpu0_stream);

				/* Now that we have calculated the model lightcurve brightnesses
				 * y at each of the epochs x, we use cubic spline interpolation
				 * (Numerical Recipes routines spline and splint) to get model
				 * lightcurve brightness fit[i] at each OBSERVATION epoch t[i],
				 * with i=1,2,...,n.  This will allow us (in routine chi2) to
				 * compare model to data (fit[i] to obs[i]) to get chi squared.
				 * Note that vector y2 contains the second derivatives of the
				 * interpolating function at the calculation epochs x. */
				gpuErrchk(cudaSetDevice(GPU0));
				/* First make a pointer for u and cudaMalloc device memory for it */
				gpuErrchk(cudaMemset(u, 0, nfrm_alloc_max*sizeof(double)));

				lghtcrv_spline_streams_test_krnl<<<1,1>>>(ddat, s, 2.0e30,
						2.0e30, u, hnframes[s]);
				checkErrorAfterKernelLaunch("lghtcrv_spline_streams_krnl");

				/* Change launch parameters from ncalc threads to n threads */
				lghtcrv_splint_streams3_test_krnl<<<1,1>>>(ddat, s, hlc_n[s], hnframes[s]);
				checkErrorAfterKernelLaunch("lghtcrv_splint_streams_krnl");
			}
			break;
		default:
			bailout("vary_params.c: can't handle this dataset type yet\n");
		}
	}

	/* Calculate the zmax, radar cross-section, optical brightness, and cosine
	 * subradar latitude */
	double dd_zmax, rd_xsec, opt_brtns, cs_sb_rdr_lat;
	set_four_parameters_krnl<<<1,1>>>(ddat);
	checkErrorAfterKernelLaunch("vps_set_four_parameters, line ");
	gpuErrchk(cudaMemcpyFromSymbol(&dd_zmax, mgpu_deldop_zmax,
			sizeof(double), 0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&rd_xsec, mgpu_rad_xsec,
			sizeof(double), 0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&opt_brtns, mgpu_opt_brightness,
			sizeof(double), 0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&cs_sb_rdr_lat, mgpu_cos_subradarlat,
			sizeof(double), 0, cudaMemcpyDeviceToHost));

	*deldop_zmax = dd_zmax;
	*rad_xsec = rd_xsec;
	*opt_brightness = opt_brtns;
	*cos_subradarlat = cs_sb_rdr_lat;

	/* GPU0 de-allocations */
	gpuErrchk(cudaSetDevice(GPU0));
	gpuErrchk(cudaFree(pos0));
	gpuErrchk(cudaFree(dframe0));
	gpuErrchk(cudaFree(ddframe0));
	gpuErrchk(cudaFree(posn0));
	gpuErrchk(cudaFree(ndel0));
	gpuErrchk(cudaFree(ndop0));
	gpuErrchk(cudaFree(outbndarr0));
	gpuErrchk(cudaFree(compute_xsec0));
	gpuErrchk(cudaFree(xylim0));
	gpuErrchk(cudaFree(bistatic0));
	gpuErrchk(cudaFree(so0));
	gpuErrchk(cudaFree(pixels_per_km0));
	gpuErrchk(cudaFree(compute_brightness));
	gpuErrchk(cudaFree(compute_zmax));
	gpuErrchk(cudaFree(compute_cosdelta));

	/* GPU1 allocations */
	gpuErrchk(cudaSetDevice(GPU1));
	gpuErrchk(cudaFree(pos1));
	gpuErrchk(cudaFree(dframe1));
	gpuErrchk(cudaFree(ddframe1));
	gpuErrchk(cudaFree(posn1));
	gpuErrchk(cudaFree(ndel1));
	gpuErrchk(cudaFree(ndop1));
	gpuErrchk(cudaFree(outbndarr1));
	gpuErrchk(cudaFree(compute_xsec1));
	gpuErrchk(cudaFree(xylim1));
	gpuErrchk(cudaFree(bistatic1));
	gpuErrchk(cudaFree(so1));
	gpuErrchk(cudaFree(pixels_per_km1));

	/* Now back to GPU0 */
	gpuErrchk(cudaSetDevice(GPU0));

}

