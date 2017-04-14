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
#include "head.h"
}
__device__ float st_deldop_cross_section=0.0, st_doppler_cross_section=0.0,
		st_sum_rad_xsec=0.0, st_sum_cos_subradarlat=0.0, st_sum_deldop_zmax,
		st_sum_opt_brightness;
__device__ double st_deldop_zmax=0.0, st_rad_xsec=0.0, st_opt_brightness=0.0,
		st_cos_subradarlat=0.0;
__device__ int vpst_n, vpst_nf, vps_weight, vpst_dd_compute_zmax, vps_ncalc,
		vpst_dd_compute_cosdelta;

__global__ void vps_get_lghtcrv_params_krnl(struct dat_t *ddat, int s) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		vpst_n = ddat->set[s].desc.lghtcrv.n;
		vps_ncalc = ddat->set[s].desc.lghtcrv.ncalc;
		vps_weight = ddat->set[s].desc.lghtcrv.weight;
	}
}
__global__ void vpst_init_krnl1() {
	/* Single-threaded kernel */
	int s;
	if (threadIdx.x == 0) {
		/* Initialize __device__ (file scope) variables to zero */
		st_deldop_cross_section = 0.0;
		st_doppler_cross_section = 0.0;
		st_sum_rad_xsec = 0.0;
		st_sum_cos_subradarlat = 0.0;
		st_sum_deldop_zmax = 0.0;
		st_sum_opt_brightness  = 0.0;
		st_deldop_zmax = 0.0;
		st_rad_xsec = 0.0;
		st_opt_brightness = 0.0;
		st_cos_subradarlat = 0.0;
	}
}
__global__ void vpst_init_krnl2(
		struct par_t *dpar,
		struct mod_t *dmod,
		struct dat_t *ddat,
		struct vertices_t **verts,
		unsigned char *type,
		int *nframes,
		int *lc_n,
		int *compute_zmax,
		int *compute_cosdelta,
		int *compute_brightness,
		int nsets,
		int c) {
	/* Single-threaded kernel */
	int s;
	if (threadIdx.x == 0) {
		/* Initialize __device__ (file scope) variables to zero */
		st_deldop_cross_section = 0.0;
		st_doppler_cross_section = 0.0;
		st_sum_rad_xsec = 0.0;
		st_sum_cos_subradarlat = 0.0;
		st_sum_deldop_zmax = 0.0;
		st_sum_opt_brightness  = 0.0;
		st_deldop_zmax = 0.0;
		st_rad_xsec = 0.0;
		st_opt_brightness = 0.0;
		st_cos_subradarlat = 0.0;
		verts[0] = &dmod->shape.comp[c].real;
		vpst_nf = verts[0]->nf;

		for (s=0; s<nsets; s++) {
			type[s] = ddat->set[s].type;

			switch(type[s]) {
			case DELAY:
				nframes[s] = ddat->set[s].desc.deldop.nframes;
				lc_n[s] = 0;
				compute_cosdelta[s] = (dpar->vary_dopscale != VARY_NONE
						&& ddat->set[s].desc.deldop.dopscale.state != 'c');
				compute_zmax[s] = (dpar->vary_delcor0 != VARY_NONE
						&& ddat->set[s].desc.deldop.delcor.a[0].state != 'c');
				break;
			case DOPPLER:
				nframes[s] = ddat->set[s].desc.doppler.nframes;
				lc_n[s] = 0;
				compute_cosdelta[s] = (dpar->vary_dopscale != VARY_NONE &&
						ddat->set[s].desc.doppler.dopscale.state != 'c');
				break;
			case POS:
				nframes[s] = ddat->set[s].desc.poset.nframes;
				lc_n[s] = 0;
				break;
			case LGHTCRV:
				nframes[s] = ddat->set[s].desc.lghtcrv.ncalc;
				lc_n[s] = ddat->set[s].desc.lghtcrv.n;
				compute_brightness[s] = (dpar->vary_optalb != VARY_NONE
						&& ddat->set[s].desc.lghtcrv.cal.state == 'c');

			}
		}
	}
}
__global__ void vpst_init_krnl3(
		struct par_t *dpar,
		struct dat_t *ddat,
		int *compute_zmax,
		int *compute_cosdelta,
		int *compute_brightness,
		unsigned char *dtype,
		int nsets) {
	/* Single-threaded kernel */
	int s;
	if (threadIdx.x == 0) {
		/* Initialize __device__ (file scope) variables to zero */
		st_deldop_cross_section = 0.0;
		st_doppler_cross_section = 0.0;
		st_sum_rad_xsec = 0.0;
		st_sum_cos_subradarlat = 0.0;
		st_sum_deldop_zmax = 0.0;
		st_sum_opt_brightness  = 0.0;
		st_deldop_zmax = 0.0;
		st_rad_xsec = 0.0;
		st_opt_brightness = 0.0;
		st_cos_subradarlat = 0.0;

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

__global__ void set_ae_oe_streams_krnl(struct pos_t *pos, double ae[3][3],
		double oe[3][3], double se[3][3], unsigned char type) {
	/* 9-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % 3;
	int j = offset / 3;

	if (offset < 9) {
		pos->ae[i][j] =	ae[i][j];
		pos->oe[i][j] =	oe[i][j];

		if (type == LGHTCRV)
			pos->se[i][j] = se[i][j];
		//			pos->se[i][j] = ddat->set[s].desc.lghtcrv.rend[f].se[i][j];

		//
		//		switch(type) {
		//
		//		case DOPPLER:
		//			vp_pos->ae[i][j] =	ddat->set[s].desc.doppler.frame[f].view[ddat->set[s].desc.doppler.v0].ae[i][j];
		//			vp_pos->oe[i][j] =	ddat->set[s].desc.doppler.frame[f].view[ddat->set[s].desc.doppler.v0].oe[i][j];
		//			break;
		//		case LGHTCRV:
		//			vp_pos->ae[i][j] = ddat->set[s].desc.lghtcrv.rend[f].ae[i][j];
		//			vp_pos->oe[i][j] = ddat->set[s].desc.lghtcrv.rend[f].oe[i][j];
		//
		//		}
		//		/* The following is a single-thread task */
		//		if (threadIdx.x == 0) {
		//			if (ddat->set[s].type == LGHTCRV)
		//				vp_pos->bistatic = 1;
		//			else
		//				vp_pos->bistatic = 0;
		//		}
	}
}
__global__ void clrvect_streams_krnl(struct dat_t *ddat, int size, int s, int f) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;

	if (offset < size) {
		if (ddat->set[s].type == DELAY)
			ddat->set[s].desc.deldop.frame[f].fit_s[offset] = 0.0;
		if (ddat->set[s].type == DOPPLER)
			ddat->set[s].desc.doppler.frame[f].fit_s[offset] = 0.0;
	}
}
__global__ void deviceMemset(struct dat_t *ddat, int s, int f, const float val, int *size)
{
	volatile int tidx = threadIdx.x + blockIdx.x * blockDim.x;
	volatile int stride = gridDim.x * blockDim.x;

	for (int i = tidx; i < size[f]; i+=stride) {
		ddat->set[s].desc.deldop.frame[f].fit_s[i] = val; }
}
__global__ void zmax_finalize_streams_krnl(struct dat_t *ddat, float value, int s, int f) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		st_sum_deldop_zmax += value * __double2float_rn(ddat->set[s].desc.deldop.frame[f].weight);
	}
}
__global__ void compute_xsec_final_streams_krnl(struct dat_t *ddat, float frm_xsec,
		int s, int f) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		switch (ddat->set[s].type) {
		case DELAY:
			st_deldop_cross_section = __double2float_rd(ddat->set[s].desc.deldop.frame[f].overflow_xsec);
			st_deldop_cross_section += frm_xsec; // fit is the end result of parallel reduction
			st_deldop_cross_section *= __double2float_rd(ddat->set[s].desc.deldop.frame[f].cal.val);
			st_sum_rad_xsec += st_deldop_cross_section *
					__double2float_rd(ddat->set[s].desc.deldop.frame[f].weight);
			break;
		case DOPPLER:
			st_doppler_cross_section = __double2float_rd(ddat->set[s].desc.doppler.frame[f].overflow_xsec);
			st_doppler_cross_section += frm_xsec;
			st_doppler_cross_section *= __double2float_rd(ddat->set[s].desc.doppler.frame[f].cal.val);
			st_sum_rad_xsec += st_doppler_cross_section *
					__double2float_rd(ddat->set[s].desc.doppler.frame[f].weight);

			break;
		}
	}
}
__global__ void compute_cosdelta_streams_krnl(struct dat_t *ddat, int s, int f) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
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
			st_sum_cos_subradarlat += cos_delta*weight;
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
			st_sum_cos_subradarlat += cos_delta*weight;
		}
	}
}
__global__ void compute_cosdelta_streams_f_krnl(struct dat_t *ddat, int s, int f) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		int j, view;
		float weight, cos_delta, oa[3][3], to_earth[3];

		switch(ddat->set[s].type){
		case DELAY:
			view = ddat->set[s].desc.deldop.v0;
			/* oa = matrix to transform body-fixed to observer coordinates  */
			/* to_earth = normalized target-to-Earth vector in body-fixed coords  */
			dev_mtrnsps4(oa, ddat->set[s].desc.deldop.frame[f].view[view].ae);
			dev_mmmul4(oa, ddat->set[s].desc.deldop.frame[f].view[view].oe, oa);
			for (j=0; j<=2; j++)
				to_earth[j] = oa[2][j];
			cos_delta = sqrt(to_earth[0]*to_earth[0] + to_earth[1]*to_earth[1]);
			weight = __double2float_rn(ddat->set[s].desc.deldop.frame[f].weight);
			st_sum_cos_subradarlat += cos_delta*weight;
			break;
		case DOPPLER:
			view = ddat->set[s].desc.doppler.v0;
			/* oa = matrix to transform body-fixed to observer coordinates  */
			/* to_earth = normalized target-to-Earth vector in body-fixed coords  */
			dev_mtrnsps4(oa, ddat->set[s].desc.doppler.frame[f].view[view].ae);
			dev_mmmul4(oa, ddat->set[s].desc.doppler.frame[f].view[view].oe, oa);
			for (j=0; j<=2; j++)
				to_earth[j] = oa[2][j];
			cos_delta = sqrt(to_earth[0]*to_earth[0] + to_earth[1]*to_earth[1]);
			weight = __double2float_rn(ddat->set[s].desc.doppler.frame[f].weight);
			st_sum_cos_subradarlat += cos_delta*weight;
		}
	}
}
__global__ void posclr_streams_krnl(struct pos_t **pos, int *posn, int f)
{
	/* Multi-threaded kernel (2*pos->n + 1)^2 threads) */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int n = posn[f];
	int nx = 2 * n + 1;
	int i = (offset % nx) - n;
	int j = (offset / nx) - n;

	if (offset < (nx*nx)) {
		/* For each POS pixel, zero out the optical brightness (b) and
		 * cos(scattering angle), reset the z coordinate (distance from COM towards
		 * Earth) to a dummy value, and reset the body, component, and facet onto
		 * which the pixel center projects to  dummy values                  */
		pos[f]->body[i][j] = pos[f]->comp[i][j] = pos[f]->f[i][j] = -1;

		pos[f]->b_s[offset] = pos[f]->cose_s[offset] = 0.0;
		pos[f]->z_s[offset] = -HUGENUMBER;

		/* In the x direction, reset the model's leftmost and rightmost
		 * pixel number to dummy values, and similarly for the y direction   */
		pos[f]->xlim[0] = pos[f]->ylim[0] =  n;
		pos[f]->xlim[1] = pos[f]->ylim[1] = -n;

		/* For a bistatic situation (lightcurve or plane-of-sky dataset), zero out
		 * cos(incidence angle) and reset the distance towards the sun, the body,
		 * component, and facet numbers as viewed from the sun, and the model's
		 * maximum projected extent as viewed from the sun to dummy values    */
		if (pos[f]->bistatic) {
			pos[f]->bodyill[i][j] = pos[f]->compill[i][j] = pos[f]->fill[i][j] = -1;

			pos[f]->cosill_s[offset] = 0.0;
			pos[f]->zill_s[offset] = -HUGENUMBER;

			pos[f]->xlim2[0] = pos[f]->ylim2[0] =  n;
			pos[f]->xlim2[1] = pos[f]->ylim2[1] = -n;
		}
	}
}
__global__ void posmask_init_streams_krnl(struct pos_t **pos, double3 *so,
		double *pixels_per_km, int f) {
	/* This single-threaded kernel performs the first few tasks (outside the
	 * pixel loop) of routine posmask.	 */
	if (threadIdx.x == 0) {
		dev_mtrnsps2(so, pos[f]->oe, f);
		dev_mmmul2(so, pos[f]->se, so, f);
		pixels_per_km[f] = 1/pos[f]->km_per_pixel;
	}
}
__global__ void posmask_init_streams_f_krnl(struct pos_t **pos, float3 *so,
		float *pixels_per_km, int f) {
	/* This single-threaded kernel performs the first few tasks (outside the
	 * pixel loop) of routine posmask.	 */
	if (threadIdx.x == 0) {
		dev_mtrnsps3(so, pos[f]->oe, f);
		dev_mmmul3(so, pos[f]->se, so, f);
		pixels_per_km[f] = 1/(__double2float_rn(pos[f]->km_per_pixel));
	}
}
__global__ void posmask_init_streams2_krnl(struct pos_t **pos, double3 *so,
		float *pixels_per_km, int f) {
	/* This single-threaded kernel performs the first few tasks (outside the
	 * pixel loop) of routine posmask.	 */
	if (threadIdx.x == 0) {
		dev_mtrnsps2(so, pos[f]->oe, f);
		dev_mmmul2(so, pos[f]->se, so, f);
		pixels_per_km[f] = (float) 1/pos[f]->km_per_pixel;
	}
}
__global__ void posmask_streams_krnl(
		struct par_t *dpar,
		struct pos_t **pos,
		double3 *so,
		double *pixels_per_km,
		int *posn,
		int nThreads,
		int xspan,
		int f)
{
	/* multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int n = posn[f];
	int i = offset % xspan - n;
	int j = offset / xspan - n;
	double tol = dpar->mask_tol;
	int im, jm, i1, j1, i2, j2, i_sign, j_sign, pxa, pxa1, pxa2;
	double xk[3], i0_dbl, j0_dbl, zill, t, u, bignum;
	bignum = 0.99*HUGENUMBER;  /* z = -HUGENUMBER for blank-sky pixels */

	/*  Loop through all POS pixels  */
	if (offset < nThreads) {

		if (pos[f]->cose_s[offset] != 0.0) {     /* if there's something there */
			xk[0] = i*pos[f]->km_per_pixel;     /* calculate 3D position */
			xk[1] = j*pos[f]->km_per_pixel;
			xk[2] = pos[f]->z_s[offset];

			/* Given the observer coordinates x of of POS pixel (i,j), find
			 * which pixel (im,jm) this corresponds to in the projected view as
			 * seen from the source (sun)             */
			dev_cotrans6(xk, so, xk, 1, f);	/* go into source coordinates */
			i0_dbl = xk[0]*pixels_per_km[f];     /* unrounded (double precision) */
			j0_dbl = xk[1]*pixels_per_km[f];
			im = dev_vp_iround( i0_dbl);            /* center of nearest pixel in mask */
			jm = dev_vp_iround( j0_dbl);
			pxa = (jm+posn[f])*xspan + (im+posn[f]); /* The float single pointer address */

			/* If center of projected pixel "seen" from source (as determined
			 * by routine posvis) lies within the boundaries of the mask,
			 * projects onto model rather than onto blank space, and represents
			 * a body, component, and facet different from those seen in the
			 * POS, calculate distance from mask pixel to source and compare to
			 * distance from POS pixel to source.                             */
			if (fabs(i0_dbl) < n && fabs(j0_dbl) < n
					&& pos[f]->zill_s[pxa] > -bignum
					&&(pos[f]->f[i][j]    != pos[f]->fill[im][jm]    ||
							pos[f]->comp[i][j] != pos[f]->compill[im][jm] ||
							pos[f]->body[i][j] != pos[f]->bodyill[im][jm]    )) {

				/* Rather than using distance towards source of CENTER of mask
				 * pixel, use bilinear interpolation to get distance towards
				 * source where the line between source and POS pixel's center
				 * intersects the mask pixel.                                */
				i1 = (int) floor( i0_dbl);
				j1 = (int) floor( j0_dbl);
				pxa = (j1+posn[f])*xspan + (i1+posn[f]);
				pxa1 = (j1+posn[f]+1)*xspan + (i1+posn[f]);

				if (pos[f]->zill_s[pxa]     > -bignum &&
						pos[f]->zill_s[pxa+1]   > -bignum &&
						pos[f]->zill_s[pxa1]   > -bignum &&
						pos[f]->zill_s[pxa1+1] > -bignum    ) {

					/* Do standard bilinear interpolation: None of the four
					 * surrounding "grid square" pixels in the mask is
					 * blank sky                           */
					t = i0_dbl - i1;
					u = j0_dbl - j1;
					zill = (1 - t)*(1 - u)*pos[f]->zill_s[pxa]
					      + t*(1 - u)*pos[f]->zill_s[pxa+1]
					      + t*u*pos[f]->zill_s[pxa1+1]
					      + (1 - t)*u*pos[f]->zill_s[pxa1];
				} else {
					/* The following code block is a kludge: One or more of the
					 * four surrounding "grid square" pixels in mask is blank
					 * sky, so standard bilinear interpolation won't work.  */
					pxa = (jm+posn[f])*xspan + (im+posn[f]);
					zill = pos[f]->zill_s[pxa];

					i_sign = (i0_dbl >= im) ? 1 : -1;
					i2 = im + i_sign;
					pxa1 = (jm+posn[f])*xspan + (i2+posn[f]);

					if (abs(i2) <= n && pos[f]->zill_s[pxa1] > -bignum) {
						zill += fabs(i0_dbl - im)
           				  * (pos[f]->zill_s[pxa1] - pos[f]->zill_s[pxa]);
					} else {
						i2 = im - i_sign;
						pxa1 = (jm+posn[f])*xspan + (i2+posn[f]);
						if (abs(i2) <= n && pos[f]->zill_s[pxa1] > -bignum)
							zill -= fabs(i0_dbl - im)
							* (pos[f]->zill_s[pxa1] - pos[f]->zill_s[pxa]);
					}

					j_sign = (j0_dbl >= jm) ? 1 : -1;
					j2 = jm + j_sign;
					pxa1 = (j2+posn[f])*xspan + (im+posn[f]);

					if (abs(j2) <= n && pos[f]->zill_s[pxa1] > -bignum) {
						zill += fabs(j0_dbl - jm)
                        	  * (pos[f]->zill_s[pxa1] - pos[f]->zill_s[pxa]);
					} else {
						j2 = jm - j_sign;
						pxa1 = (j2+posn[f])*xspan + (im+posn[f]);

						if (abs(j2) <= n && pos[f]->zill_s[pxa1] > -bignum)
							zill -= fabs(j0_dbl - jm)
							* (pos[f]->zill_s[pxa1] - pos[f]->zill_s[pxa]);
					}
				}
				/* If interpolated point within mask pixel is at least tol km
				 * closer to source than is the center of POS pixel, the facet
				 * represented by the mask pixel is shadowing the POS pixel:
				 * represent this by setting
				 * 		cos(scattering angle) = 0.0 for the POS pixel.      */
				if (zill - xk[2] > tol)
					pos[f]->cose_s[offset] = 0.0;
			}
		}
	}
}
__global__ void posmask_streams_f_krnl(
		struct par_t *dpar,
		struct pos_t **pos,
		float3 *so,
		float *pixels_per_km,
		int nThreads,
		int xspan,
		int f)
{
	/* multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	__shared__ int n, kmpxl;
	int i = offset % xspan - n;
	int j = offset / xspan - n;
	float tol = dpar->mask_tol;
	int im, jm, i1, j1, i2, j2, i_sign, j_sign, pxa, pxa1, pxa2;
	float xk[3], i0_dbl, j0_dbl, zill, t, u, bignum;
	bignum = 0.99*HUGENUMBER;  /* z = -HUGENUMBER for blank-sky pixels */

	/* Set the shared variables once per block */
	if (threadIdx.x == 0) {
		n = pos[f]->n;
		kmpxl = pos[f]->km_per_pixel;
	}

	/*  Loop through all POS pixels  */
	if (offset < nThreads) {

		if (pos[f]->cose_s[offset] != 0.0) {     /* if there's something there */
			xk[0] = i*kmpxl;     /* calculate 3D position */
			xk[1] = j*kmpxl;
			xk[2] = pos[f]->z_s[offset];

			/* Given the observer coordinates x of of POS pixel (i,j), find
			 * which pixel (im,jm) this corresponds to in the projected view as
			 * seen from the source (sun)             */
			dev_cotrans3(xk, so, xk, 1, f);	/* go into source coordinates */
			i0_dbl = xk[0]*pixels_per_km[f];     /* unrounded (double precision) */
			j0_dbl = xk[1]*pixels_per_km[f];
			im = dev_vp_iroundf( i0_dbl);            /* center of nearest pixel in mask */
			jm = dev_vp_iroundf( j0_dbl);
			pxa = (jm+n)*xspan + (im+n); /* The float single pointer address */

			/* If center of projected pixel "seen" from source (as determined
			 * by routine posvis) lies within the boundaries of the mask,
			 * projects onto model rather than onto blank space, and represents
			 * a body, component, and facet different from those seen in the
			 * POS, calculate distance from mask pixel to source and compare to
			 * distance from POS pixel to source.                             */
			if (fabs(i0_dbl) < n && fabs(j0_dbl) < n
					&& pos[f]->zill_s[pxa] > -bignum
					&&(pos[f]->f[i][j]    != pos[f]->fill[im][jm]    ||
							pos[f]->comp[i][j] != pos[f]->compill[im][jm] ||
							pos[f]->body[i][j] != pos[f]->bodyill[im][jm]    )) {

				/* Rather than using distance towards source of CENTER of mask
				 * pixel, use bilinear interpolation to get distance towards
				 * source where the line between source and POS pixel's center
				 * intersects the mask pixel.                                */
				i1 = (int) floor( i0_dbl);
				j1 = (int) floor( j0_dbl);
				pxa = (j1+n)*xspan + (i1+n);
				pxa1 = (j1+n+1)*xspan + (i1+n);

				if (pos[f]->zill_s[pxa]     > -bignum &&
						pos[f]->zill_s[pxa+1]   > -bignum &&
						pos[f]->zill_s[pxa1]   > -bignum &&
						pos[f]->zill_s[pxa1+1] > -bignum    ) {

					/* Do standard bilinear interpolation: None of the four
					 * surrounding "grid square" pixels in the mask is
					 * blank sky                           */
					t = i0_dbl - i1;
					u = j0_dbl - j1;
					zill = (1 - t)*(1 - u)*pos[f]->zill_s[pxa]
					      + t*(1 - u)*pos[f]->zill_s[pxa+1]
					      + t*u*pos[f]->zill_s[pxa1+1]
					      + (1 - t)*u*pos[f]->zill_s[pxa1];
				} else {
					/* The following code block is a kludge: One or more of the
					 * four surrounding "grid square" pixels in mask is blank
					 * sky, so standard bilinear interpolation won't work.  */
					pxa = (jm+n)*xspan + (im+n);
					zill = pos[f]->zill_s[pxa];

					i_sign = (i0_dbl >= im) ? 1 : -1;
					i2 = im + i_sign;
					pxa1 = (jm+n)*xspan + (i2+n);

					if (abs(i2) <= n && pos[f]->zill_s[pxa1] > -bignum) {
						zill += fabs(i0_dbl - im)
           				  * (pos[f]->zill_s[pxa1] - pos[f]->zill_s[pxa]);
					} else {
						i2 = im - i_sign;
						pxa1 = (jm+n)*xspan + (i2+n);
						if (abs(i2) <= n && pos[f]->zill_s[pxa1] > -bignum)
							zill -= fabs(i0_dbl - im)
							* (pos[f]->zill_s[pxa1] - pos[f]->zill_s[pxa]);
					}

					j_sign = (j0_dbl >= jm) ? 1 : -1;
					j2 = jm + j_sign;
					pxa1 = (j2+n)*xspan + (im+n);

					if (abs(j2) <= n && pos[f]->zill_s[pxa1] > -bignum) {
						zill += fabs(j0_dbl - jm)
                        	  * (pos[f]->zill_s[pxa1] - pos[f]->zill_s[pxa]);
					} else {
						j2 = jm - j_sign;
						pxa1 = (j2+n)*xspan + (im+n);

						if (abs(j2) <= n && pos[f]->zill_s[pxa1] > -bignum)
							zill -= fabs(j0_dbl - jm)
							* (pos[f]->zill_s[pxa1] - pos[f]->zill_s[pxa]);
					}
				}
				/* If interpolated point within mask pixel is at least tol km
				 * closer to source than is the center of POS pixel, the facet
				 * represented by the mask pixel is shadowing the POS pixel:
				 * represent this by setting
				 * 		cos(scattering angle) = 0.0 for the POS pixel.      */
				if (zill - xk[2] > tol)
					pos[f]->cose_s[offset] = 0.0;
			}
		}
	}
}
__global__ void posmask_streams2_krnl(
		struct par_t *dpar,
		struct pos_t **pos,
		double3 *so,
		float *pixels_per_km,
		int *posn,
		int nThreads,
		int xspan,
		int f)
{
	/* multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int n = posn[f];
	int i = offset % xspan - n;
	int j = offset / xspan - n;
	float tol = dpar->mask_tol;
	float kmpxl = (float)pos[f]->km_per_pixel;
	int im, jm, i1, j1, i2, j2, i_sign, j_sign, pxa, pxa1, pxa2;
	//double xk[3], i0_dbl, j0_dbl, zill, t, u, bignum;
	float3 xk;
	float i0_f, j0_f, zill, t, u, bignum;
	bignum = 0.99*HUGENUMBER;  /* z = -HUGENUMBER for blank-sky pixels */

	/*  Loop through all POS pixels  */
	if (offset < nThreads) {

		if (pos[f]->cose_s[offset] != 0.0) {     /* if there's something there */
			xk.x = i*kmpxl;     /* calculate 3D position */
			xk.y = j*kmpxl;
			xk.z = pos[f]->z_s[offset];

			/* Given the observer coordinates x of of POS pixel (i,j), find
			 * which pixel (im,jm) this corresponds to in the projected view as
			 * seen from the source (sun)             */
			dev_cotrans7(&xk, so, xk, 1, f);	/* go into source coordinates */
			i0_f = xk.x*pixels_per_km[f];     /* unrounded (double precision) */
			j0_f = xk.y*pixels_per_km[f];
			im = dev_vp_iroundf(i0_f);            /* center of nearest pixel in mask */
			jm = dev_vp_iroundf(j0_f);
			pxa = (jm+posn[f])*xspan + (im+posn[f]); /* The float single pointer address */

			/* If center of projected pixel "seen" from source (as determined
			 * by routine posvis) lies within the boundaries of the mask,
			 * projects onto model rather than onto blank space, and represents
			 * a body, component, and facet different from those seen in the
			 * POS, calculate distance from mask pixel to source and compare to
			 * distance from POS pixel to source.                             */
			if (fabs(i0_f) < n && fabs(j0_f) < n
					&& pos[f]->zill_s[pxa] > -bignum
					&&(pos[f]->f[i][j]    != pos[f]->fill[im][jm]    ||
							pos[f]->comp[i][j] != pos[f]->compill[im][jm] ||
							pos[f]->body[i][j] != pos[f]->bodyill[im][jm]    )) {

				/* Rather than using distance towards source of CENTER of mask
				 * pixel, use bilinear interpolation to get distance towards
				 * source where the line between source and POS pixel's center
				 * intersects the mask pixel.                                */
				i1 = (int) floor( i0_f);
				j1 = (int) floor( j0_f);
				pxa = (j1+posn[f])*xspan + (i1+posn[f]);
				pxa1 = (j1+posn[f]+1)*xspan + (i1+posn[f]);

				if (pos[f]->zill_s[pxa]     > -bignum &&
						pos[f]->zill_s[pxa+1]   > -bignum &&
						pos[f]->zill_s[pxa1]   > -bignum &&
						pos[f]->zill_s[pxa1+1] > -bignum    ) {

					/* Do standard bilinear interpolation: None of the four
					 * surrounding "grid square" pixels in the mask is
					 * blank sky                           */
					t = i0_f - i1;
					u = j0_f - j1;
					zill = (1 - t)*(1 - u)*pos[f]->zill_s[pxa]
					      + t*(1 - u)*pos[f]->zill_s[pxa+1]
					      + t*u*pos[f]->zill_s[pxa1+1]
					      + (1 - t)*u*pos[f]->zill_s[pxa1];
				} else {
					/* The following code block is a kludge: One or more of the
					 * four surrounding "grid square" pixels in mask is blank
					 * sky, so standard bilinear interpolation won't work.  */
					pxa = (jm+posn[f])*xspan + (im+posn[f]);
					zill = pos[f]->zill_s[pxa];

					i_sign = (i0_f >= im) ? 1 : -1;
					i2 = im + i_sign;
					pxa1 = (jm+posn[f])*xspan + (i2+posn[f]);

					if (abs(i2) <= n && pos[f]->zill_s[pxa1] > -bignum) {
						zill += fabs(i0_f - im)
           				  * (pos[f]->zill_s[pxa1] - pos[f]->zill_s[pxa]);
					} else {
						i2 = im - i_sign;
						pxa1 = (jm+posn[f])*xspan + (i2+posn[f]);
						if (abs(i2) <= n && pos[f]->zill_s[pxa1] > -bignum)
							zill -= fabs(i0_f - im)
							* (pos[f]->zill_s[pxa1] - pos[f]->zill_s[pxa]);
					}

					j_sign = (j0_f >= jm) ? 1 : -1;
					j2 = jm + j_sign;
					pxa1 = (j2+posn[f])*xspan + (im+posn[f]);

					if (abs(j2) <= n && pos[f]->zill_s[pxa1] > -bignum) {
						zill += fabs(j0_f - jm)
                        	  * (pos[f]->zill_s[pxa1] - pos[f]->zill_s[pxa]);
					} else {
						j2 = jm - j_sign;
						pxa1 = (j2+posn[f])*xspan + (im+posn[f]);

						if (abs(j2) <= n && pos[f]->zill_s[pxa1] > -bignum)
							zill -= fabs(j0_f - jm)
							* (pos[f]->zill_s[pxa1] - pos[f]->zill_s[pxa]);
					}
				}
				/* If interpolated point within mask pixel is at least tol km
				 * closer to source than is the center of POS pixel, the facet
				 * represented by the mask pixel is shadowing the POS pixel:
				 * represent this by setting
				 * 		cos(scattering angle) = 0.0 for the POS pixel.      */
				if (zill - xk.z > tol)
					pos[f]->cose_s[offset] = 0.0;
			}
		}
	}
}
__global__ void lghtcrv_copy_y_streams_krnl(struct dat_t *ddat, double *host_value,
		int set, int size) {
	/* Multi-threaded kernel, not streamed. */
	int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
	if (i <= size) {
		ddat->set[set].desc.lghtcrv.y[i] = host_value[i];
	}
}
__global__ void lghtcrv_copy_y_streams_f_krnl(struct dat_t *ddat, float *host_value,
		int set, int size) {
	/* Multi-threaded kernel, not streamed. */
	int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
	if (i <= size) {
		ddat->set[set].desc.lghtcrv.y[i] = (double)host_value[i];
	}
}
__global__ void lghtcrv_spline_streams_krnl(struct dat_t *ddat, int set, double
		yp1, double ypn, double *u) {
	/*(double *x  - lghtcrv->x
	 * double *y  - lghtcrv->y
	 * int n      - calc
	 * double yp1 - 2.0e30
	 * double ypn - 2.0e30
	 * double *y2 - lghtcrv->y2)*/

	/* Multi-threaded kernel */
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int n = vps_ncalc;
	int k = n - 1 - i;
	double *x = ddat->set[set].desc.lghtcrv.x;
	double *y = ddat->set[set].desc.lghtcrv.y;
	double *y2 = ddat->set[set].desc.lghtcrv.y2;
	double p,qn,sig,un;

	/* Perform single-thread task */
	if (i == 0) {
		if (yp1 > 0.99e30)
			y2[1]=u[1]=0.0;
		else {
			y2[1] = -0.5;
			u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
		}
	}
	__syncthreads();

	if (i > 1 && i < n-1) {

		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	__syncthreads();

	/* Perform another single-thread task */
	if (i == 1) {
		if (ypn > 0.99e30)
			qn=un=0.0;
		else {
			qn=0.5;
			un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
		}
		y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
	}
	__syncthreads();

	if (k <= (n-1) && k >= 1)
		y2[k]=y2[k]*y2[k+1]+u[k];

	__syncthreads();
}
__global__ void lghtcrv_spline_streams_f_krnl(struct dat_t *ddat, int set, float
		yp1, float ypn, float *u) {
	/*(double *x  - lghtcrv->x
	 * double *y  - lghtcrv->y
	 * int n      - calc
	 * double yp1 - 2.0e30
	 * double ypn - 2.0e30
	 * double *y2 - lghtcrv->y2)*/

	/* Multi-threaded kernel */
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int n = vps_ncalc;
	int k = n - 1 - i;

	/* The following will be turned into floats */
	float *x = ddat->set[set].desc.lghtcrv.x_s;
	float *y = ddat->set[set].desc.lghtcrv.y_s;
	float *y2 = ddat->set[set].desc.lghtcrv.y2_s;
	float p,qn,sig,un;

	/* Perform single-thread task */
	if (i == 0) {
		if (yp1 > 0.99e30)
			y2[1]=u[1]=0.0;
		else {
			y2[1] = -0.5;
			u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
		}
	}
	__syncthreads();

	if (i > 1 && i < n-1) {

		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	__syncthreads();

	/* Perform another single-thread task */
	if (i == 1) {
		if (ypn > 0.99e30)
			qn=un=0.0;
		else {
			qn=0.5;
			un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
		}
		y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
	}
	__syncthreads();

	if (k <= (n-1) && k >= 1)
		y2[k]=y2[k]*y2[k+1]+u[k];

	__syncthreads();
}
__global__ void lghtcrv_splint_streams2_krnl(struct dat_t *ddat, int set)
{
	/* This is an n-threaded kernel where n = lghtcrv-> */
	/* Parameters:
	 * double *xa  - lghtcrv->x
	 * double *ya  - lghtcrv->y
	 * double *y2a - lghtcrv->y2
	 * int n       - ncalc
	 * double x    - lghtcrv->t[i][lghtcrv->v0]
	 * double *y   - lghtcrv->fit[i]	 *
	 */

	int i = blockIdx.x * blockDim.x + threadIdx.x;
	double *xa	= ddat->set[set].desc.lghtcrv.x;
	double *ya	= ddat->set[set].desc.lghtcrv.y;
	double *y2a	= ddat->set[set].desc.lghtcrv.y2;
	double x 	= ddat->set[set].desc.lghtcrv.t[i][ddat->set[set].desc.lghtcrv.v0];
	double *y 	= &ddat->set[set].desc.lghtcrv.t[i][ddat->set[set].desc.lghtcrv.v0];
	int n = vps_ncalc;

	int klo,khi,k;
	double h,b,a;

	if (i < vpst_n) {

		klo=1;
		khi=n;
		while (khi-klo > 1) {
			k=(khi + klo) >> 1;
			if (xa[k] > x) khi=k;
			else klo=k;
		}
		h = xa[khi] - xa[klo];
		if (h == 0.0) printf("Bad XA input to routine SPLINT");
		a=(xa[khi] - x) / h;
		b=(x - xa[klo]) / h;
		*y=a*ya[klo] + b*ya[khi] + ((a*a*a-a) * y2a[klo] + (b*b*b-b) * y2a[khi]) * (h*h) /6.0;
	}
}
__global__ void lghtcrv_splint_streams_f_krnl(struct dat_t *ddat, int set)
{
	/* This is an n-threaded kernel where n = lghtcrv->ncalc */
	/* This version uses floats instead of doubles */

	int i = blockIdx.x * blockDim.x + threadIdx.x;
	float *xa	= ddat->set[set].desc.lghtcrv.x_s;
	float *ya	= ddat->set[set].desc.lghtcrv.y_s;
	float *y2a	= ddat->set[set].desc.lghtcrv.y2_s;
	float x 	= __double2float_rn(ddat->set[set].desc.lghtcrv.t[i][ddat->set[set].desc.lghtcrv.v0]);
	double *y 	= &ddat->set[set].desc.lghtcrv.t[i][ddat->set[set].desc.lghtcrv.v0];
	int n = vps_ncalc;

	int klo,khi,k;
	float h,b,a;

	if (i < vpst_n) {

		klo=1;
		khi=n;
		while (khi-klo > 1) {
			k=(khi + klo) >> 1;
			if (xa[k] > x) khi=k;
			else klo=k;
		}
		h = xa[khi] - xa[klo];
		if (h == 0.0) printf("Bad XA input to routine SPLINT");
		a=(xa[khi] - x) / h;
		b=(x - xa[klo]) / h;
		*y=a*ya[klo] + b*ya[khi] + ((a*a*a-a) * y2a[klo] + (b*b*b-b) * y2a[khi]) * (h*h) /6.0;
	}
}
__global__ void vps_set_four_parameters_krnl(struct dat_t *ddat) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		if (ddat->sum_deldop_zmax_weights > 0.0)
			st_deldop_zmax = st_sum_deldop_zmax / ddat->sum_deldop_zmax_weights;
		else
			st_deldop_zmax = 0.0;
		if (ddat->sum_rad_xsec_weights > 0.0) {
			st_rad_xsec = st_sum_rad_xsec / ddat->sum_rad_xsec_weights;			}
		else
			st_rad_xsec = 0.0;
		if (ddat->sum_opt_brightness_weights > 0.0)
			st_opt_brightness = st_sum_opt_brightness / ddat->sum_opt_brightness_weights;
		else
			st_opt_brightness = 0.0;
		if (ddat->sum_cos_subradarlat_weights > 0.0)
			st_cos_subradarlat = st_sum_cos_subradarlat / ddat->sum_cos_subradarlat_weights;
		else
			st_cos_subradarlat = 0.0;
	}
}

__host__ void vary_params_cuda_streams( struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, int action, double *deldop_zmax, double
		*rad_xsec, double *opt_brightness, double *cos_subradarlat, int nsets)
{
	/* Inputs:
	 * int action 				- dpar->action
	 * double *deldop_zmax 		- &deldop_zmax_save
	 * double *rad_xsec 		- &rad_xsec_save
	 * double *optbrightness	- &opt_brightness_save
	 * double *cos_subradarlat	- &cos_subradarlat_save
	 * int nsets 				- ddat->nsets
	 */

	double orbit_offset[3] = {0.0, 0.0, 0.0};
	int c=0, f, s, compute_brightness, compute_zmax, bistatic_all,
			compute_cosdelta, n, ncalc, nx, lghtcrv_bistatic, nframes;
	dim3 pxBLK,THD,BLKncalc;
	THD.x = maxThreadsPerBlock;
	unsigned char type;
	int nThreads, *posn, *ndel, *ndop;

	float zmax;
	double weight, *pixels_per_km, *lghtcrv_y;
	double3 *so;
	struct pos_t **pos;
	struct vertices_t **verts;

	/* Initialize the device file-scope variables */
	vpst_init_krnl1<<<1,1>>>();
	checkErrorAfterKernelLaunch("vpst_init_krnl");
	deviceSyncAfterKernelLaunch("vpst_init_krnl");

	/* Process each dataset in turn */
	for (s=0; s<nsets; s++) {
		deviceSyncAfterKernelLaunch("top of posvis_cuda_streams");
		type = ddat->set[s].type;
		if (type==DELAY)	nframes = ddat->set[s].desc.deldop.nframes;
		if (type==DOPPLER)	nframes = ddat->set[s].desc.doppler.nframes;
		if (type==POS)		nframes = ddat->set[s].desc.poset.nframes;
		if (type==LGHTCRV)	{
			nframes = ddat->set[s].desc.lghtcrv.ncalc;
			lghtcrv_y = (double *) malloc(nframes*sizeof(double));
		}
		int compute_xsec[nframes], npxls[nframes], ddsize[nframes], hndop[nframes],
			hndel[nframes], lc_xspan[nframes], *outbndarr;
		float xsec[nframes];
		dim3 BLK[nframes], ddBLK[nframes];
		cudaStream_t vp_stream[nframes];
		cudaCalloc1((void**)&pos, 		sizeof(struct pos_t*), nframes);
		cudaCalloc1((void**)&posn, 		sizeof(int), 		   nframes);
		cudaCalloc1((void**)&ndel, 		sizeof(int), 		   nframes);
		cudaCalloc1((void**)&ndop, 		sizeof(int), 		   nframes);
		cudaCalloc1((void**)&outbndarr,  sizeof(int), 		   nframes);

		switch (type) {
		case DELAY:
			/* Get computation flags */
			compute_zmax = (dpar->vary_delcor0 != VARY_NONE
					&& ddat->set[s].desc.deldop.delcor.a[0].state != 'c');
			compute_cosdelta = (dpar->vary_dopscale != VARY_NONE
					&& ddat->set[s].desc.deldop.dopscale.state != 'c');

			for (f=0; f<nframes; f++) {
				compute_xsec[f] = (dpar->vary_radalb != VARY_NONE
						&& ddat->set[s].desc.deldop.frame[f].cal.state == 'c');
				pos[f] = &ddat->set[s].desc.deldop.frame[f].pos;
				posn[f] = pos[f]->n;
				//fit_s[f] = &ddat->set[s].desc.deldop.frame[f].fit_s;
				cudaStreamCreate(&vp_stream[f]);

				/* Calculate launch parameters that will be needed later */
				nx = 2*posn[f] + 1;
				nThreads = nx*nx;
				npxls[f] = nThreads;
				BLK[f].x = floor((THD.x - 1 + nThreads) / THD.x);

				ndel[f] = hndel[f] = ddat->set[s].desc.deldop.ndel;
				ndop[f] = hndop[f] = ddat->set[s].desc.deldop.ndop;
				ddsize[f]= ndel[f] * ndop[f];
				ddBLK[f] = floor((THD.x -1 + ddsize[f]) / THD.x);

				/* Assign the ae and oe matrices */
				for (int i=0; i<=2; i++)
					for (int j=0; j<=2; j++) {
						ddat->set[s].desc.deldop.frame[f].pos.ae[i][j] =
								ddat->set[s].desc.deldop.frame[f].view[ddat->set[s].desc.deldop.v0].ae[i][j];
						ddat->set[s].desc.deldop.frame[f].pos.oe[i][j] =
								ddat->set[s].desc.deldop.frame[f].view[ddat->set[s].desc.deldop.v0].oe[i][j];
					}
				pos[f]->bistatic = 0;
			}

			/* Launch posclr_streams_krnl to initialize POS view */
			for (f=0; f<nframes; f++) {
				/* Start the if block for computing zmax and/or cross-section */
				if (compute_zmax || compute_xsec[f]) {
					posclr_streams_krnl<<<BLK[f],THD, 0, vp_stream[f]>>>(pos,posn,f);
				}
			} checkErrorAfterKernelLaunch("posclr_streams_krnl (Delay-Doppler)");

			/* Determine which POS pixels cover the target, and get distance
			 * toward Earth of each POS pixel. Pass the frame streams, too. */
			posvis_cuda_streams(dpar, dmod, ddat, orbit_offset,s,f,	0, 0, c,
					outbndarr, vp_stream);

			for (f=0; f<nframes; f++) {
				if (compute_zmax || compute_xsec[f]) {
					/* Zero out the fit delay-Doppler image and call pos2deldop
					 * to create the fit image by mapping power from the plane
					 * of sky to delay-Doppler space.    				  */

					clrvect_streams_krnl<<<ddBLK[f],THD, 0, vp_stream[f]>>>(ddat,
							ddsize[f], s, f);
				}/* End frames loop again to call pos2deldop streams version */
			} checkErrorAfterKernelLaunch("clrvect_streams_krnl");

			/* Call the CUDA pos2deldop function */
			pos2deldop_cuda_streams(dpar, dmod, ddat, pos, ndel, ndop,
					0.0, 0.0, 0.0, 0, s, nframes, 0, outbndarr, vp_stream);

			for (f=0; f<nframes; f++) {
				if (compute_zmax || compute_xsec[f]) {
					/* Compute distance toward Earth of the subradar point  */
					if (compute_zmax) {
						zmax = 0.0;
						/* Call parallel reduction function to quickly find the
						 * distance toward Earth of the subradar point */
						zmax = compute_pos_zmax(ddat, npxls[f], s, f);
						zmax_finalize_streams_krnl<<<1,1,0,vp_stream[f]>>>(ddat, zmax, s, f);
					}
				}
			}
			for (f=0; f<nframes; f++) {
				/*  Compute cross section  */
				if (compute_xsec) {
					xsec[f] = 0.0;
					xsec[f] = compute_deldop_xsec_snglkrnl(ddat, hndel[f], hndop[f], s, f);
					if (compute_xsec[f])
						compute_xsec_final_streams_krnl<<<1,1,0,vp_stream[f]>>>(ddat, xsec[f], s, f);
				}
			}
			checkErrorAfterKernelLaunch("compute_xsec_final_streams_krnl");

			if (compute_cosdelta)
				for (f=0; f<nframes; f++)
					compute_cosdelta_streams_krnl<<<1,1>>>(ddat, s, f);
			checkErrorAfterKernelLaunch("compute_cosdelta_streams_krnl");

			for (f=0; f<nframes; f++)
				cudaStreamDestroy(vp_stream[f]);

			break;
		case DOPPLER:
			/* Get computation flags */
			compute_cosdelta = (dpar->vary_dopscale != VARY_NONE &&
					ddat->set[s].desc.doppler.dopscale.state != 'c');

			for (f=0; f<nframes; f++) {
				compute_xsec[f] = (dpar->vary_radalb != VARY_NONE &&
						ddat->set[s].desc.doppler.frame[f].cal.state == 'c');
				pos[f] = &ddat->set[s].desc.doppler.frame[f].pos;
				posn[f] = pos[f]->n;
				cudaStreamCreate(&vp_stream[f]);

				/* Calculate launch parameters that will be needed later */
				nx = 2*posn[f] + 1;
				nThreads = nx*nx;
				npxls[f] = nThreads;
				BLK[f].x = floor((THD.x - 1 + nThreads) / THD.x);

				hndop[f] = ndop[f] = ddat->set[s].desc.doppler.ndop;
				ddBLK[f] = floor((THD.x -1 + ndop[f]) / THD.x);

				/* Assign the ae and oe matrices */
				for (int i=0; i<=2; i++)
					for (int j=0; j<=2; j++) {
						ddat->set[s].desc.doppler.frame[f].pos.ae[i][j] =
								ddat->set[s].desc.doppler.frame[f].view[ddat->set[s].desc.doppler.v0].ae[i][j];
						ddat->set[s].desc.doppler.frame[f].pos.oe[i][j] =
								ddat->set[s].desc.doppler.frame[f].view[ddat->set[s].desc.doppler.v0].oe[i][j];
					}
				pos[f]->bistatic = 0;
			}
			/* Clear out the plane-of-sky first */
			for (f=0; f<nframes; f++)
				if (compute_xsec[f])
					posclr_streams_krnl<<<BLK[f],THD, 0, vp_stream[f]>>>(pos,posn,f);

			checkErrorAfterKernelLaunch("posclr_streams_krnl (Doppler) ");

			/* Determine which POS pixels cover the target, and get distance
			 * toward Earth of each POS pixel. Pass the frame streams, too. */
			posvis_cuda_streams(dpar, dmod, ddat, orbit_offset,s,f,	0, 0, c,
					outbndarr, vp_stream);

			for (f=0; f<nframes; f++) {
				if (compute_xsec[f]) {
					/* Zero out the fit delay-Doppler image and call pos2deldop
					 * to create the fit image by mapping power from the plane
					 * of sky to delay-Doppler space.    				      */

					clrvect_streams_krnl<<<ddBLK[f],THD, 0, vp_stream[f]>>>(ddat,
							hndop[f], s, f);

					/* End frames loop again to call pos2deldop streams version */
				}
			} checkErrorAfterKernelLaunch("clrvect_streams_krnl");

			pos2doppler_cuda_streams(dpar, dmod, ddat, pos, 0.0, 0.0, 0.0,
					ndop, 0, s, nframes, 0, outbndarr, vp_stream);

			/* Calculate the Doppler cross-section if applicable */
			for (f=0; f<nframes; f++) {
				if (compute_xsec[f]) {
					/* Compute cross section */
					xsec[f]=0.0;
					xsec[f] = compute_doppler_xsec(ddat, ndop[f], s, f);
				}
			}

			/* Finalize the xsec calculations and calculate cosdelta if specified */
			for (f=0; f<nframes; f++) {
				if (compute_xsec[f])
					compute_xsec_final_streams_krnl<<<1,1,0,vp_stream[f]>>>(ddat, xsec[f], s, f);
				if (compute_cosdelta)
					compute_cosdelta_streams_krnl<<<1,1,0,vp_stream[f]>>>(ddat, s, f);
			}
			break;
		case POS:
			break;
		case LGHTCRV:
			/* Figure out the compute_brightness flag first */
			compute_brightness = (dpar->vary_optalb != VARY_NONE
							&& ddat->set[s].desc.lghtcrv.cal.state == 'c');
//			get_lghtcrv_cb_krnl<<<1,1>>>(dpar, ddat, s);
//			checkErrorAfterKernelLaunch("get_lghtcrv_cb_krnl");
//			gpuErrchk(cudaMemcpyFromSymbol(&compute_brightness, dcompute_brightness,
//					sizeof(int), 0, cudaMemcpyDeviceToHost));

			if (compute_brightness) {
				cudaCalloc1((void**)&so, sizeof(double3), (nframes*3));
				cudaCalloc1((void**)&pixels_per_km, sizeof(int), nframes);
				/* Launch single-thread kernel to get lghtcrv parameters */
				vps_get_lghtcrv_params_krnl<<<1,1>>>(ddat, s);
				checkErrorAfterKernelLaunch("vps_get_lghtcrv_params_krnl");
				gpuErrchk(cudaMemcpyFromSymbol(&n, vpst_n, sizeof(int),
						0, cudaMemcpyDeviceToHost));
				gpuErrchk(cudaMemcpyFromSymbol(&ncalc, vps_ncalc, sizeof(int),
						0, cudaMemcpyDeviceToHost));
				gpuErrchk(cudaMemcpyFromSymbol(&weight, vps_weight, sizeof(double),
						0, cudaMemcpyDeviceToHost));

				for (f=0; f<ncalc; f++) {

					pos[f] = &ddat->set[s].desc.lghtcrv.rend[f].pos;
					posn[f] = pos[f]->n;
					cudaStreamCreate(&vp_stream[f]);

					/* Calculate launch parameters that will be needed later */
					lc_xspan[f] = nx = 2*posn[f] + 1;
					nThreads = nx*nx;
					npxls[f] = nThreads;
					BLK[f].x = floor((THD.x - 1 + nThreads) / THD.x);

					/* Assign the ae and oe matrices */
					for (int i=0; i<=2; i++)
						for (int j=0; j<=2; j++) {
							pos[f]->ae[i][j] =	ddat->set[s].desc.lghtcrv.rend[f].ae[i][j];
							pos[f]->oe[i][j] =	ddat->set[s].desc.lghtcrv.rend[f].oe[i][j];
						}
					pos[f]->bistatic = 1;
				}

				/* Clear out the plane-of-sky first */
				for (f=0; f<ncalc; f++)
					posclr_streams_krnl<<<BLK[f],THD, 0, vp_stream[f]>>>(pos,posn,f);
				checkErrorAfterKernelLaunch("posclr_streams_krnl (Doppler) ");

				/* Determine which POS pixels cover the target */
				posvis_cuda_streams(dpar, dmod, ddat, orbit_offset,s,f,	0, 0, c,
						outbndarr, vp_stream);

				/* Now view the model from the source (sun) and get the facet
				 * number and distance toward the source of each pixel in this
				 * projected view; use this information to determine which POS
				 * pixels are shadowed */
				for (f=0; f<nframes; f++)
					if (pos[f]->bistatic)
						bistatic_all = 1;
				if (bistatic_all)
					posvis_cuda_streams(dpar, dmod, ddat, orbit_offset,s,f,	1,
							0, c, outbndarr, vp_stream);

				if (bistatic_all) {
					for (f=0; f<ncalc; f++) {
						/* Initialize this stream for the posmask kernel to follow */
						posmask_init_streams_krnl<<<1,1,0,vp_stream[f]>>>(pos,
								so, pixels_per_km, f);

						/* Now call posmask kernel for this stream, then loop
						 * to next stream and repeat 	 */
						posmask_streams_krnl<<<BLK[f],THD,0,vp_stream[f]>>>(
								dpar, pos, so, pixels_per_km, posn, npxls[f],
								lc_xspan[f], f);

					} checkErrorAfterKernelLaunch("posmask_streams_ krnl");
				}

				/* Compute model brightness for this lightcurve point */
				/* lghtcrv->y[ncalc]: calculated points for interpolation,
				 * ncalc-points total 					 */
				for (f=0; f<nframes; f++)
					lghtcrv_y[f] = apply_photo_cuda(dmod, ddat, 0, s, f);

				/* Now launch a kernel to copy it over to the actual lghtcrv */
				BLKncalc.x = floor((THD.x - 1 + ncalc) / THD.x);
				lghtcrv_copy_y_streams_krnl<<<BLKncalc,THD>>>(ddat, lghtcrv_y, s, nframes);
				checkErrorAfterKernelLaunch("lghtcrv_copy_y_streams_krnl");

				/* Now that we have calculated the model lightcurve brightnesses
				 * y at each of the epochs x, we use cubic spline interpolation
				 * (Numerical Recipes routines spline and splint) to get model
				 * lightcurve brightness fit[i] at each OBSERVATION epoch t[i],
				 * with i=1,2,...,n.  This will allow us (in routine chi2) to
				 * compare model to data (fit[i] to obs[i]) to get chi squared.
				 * Note that vector y2 contains the second derivatives of the
				 * interpolating function at the calculation epochs x. */

				/* First make a pointer for u and cudaMalloc device memory for it */
				double *u;
				gpuErrchk(cudaMalloc((void**)&u, sizeof(double) * n));

				lghtcrv_spline_streams_krnl<<<BLKncalc,THD>>>(ddat, s, 2.0e30, 2.0e30, u);
				checkErrorAfterKernelLaunch("lghtcrv_spline_streams_krnl");

				/* Change launch parameters from ncalc threads to n threads */
				BLKncalc.x = floor((THD.x - 1 + n) / THD.x);
				lghtcrv_splint_streams2_krnl<<<BLKncalc,THD>>>(ddat, s);
				checkErrorAfterKernelLaunch("lghtcrv_splint_streams_krnl");
				/* Cleanup */
				cudaFree(u);
				cudaFree(so);
				cudaFree(pixels_per_km);
				/* Destroy streams and free memory */
				for (f=0; f<nframes; f++)
					cudaStreamDestroy(vp_stream[f]);
			}
			break;
		default:
			bailout("vary_params.c: can't handle this dataset type yet\n");
		}



		cudaFree(pos);
		cudaFree(posn);
		cudaFree(ndel);
		cudaFree(ndop);
		cudaFree(outbndarr);
	}

	/* Calculate the zmax, radar cross-section, optical brightness, and cosine
	 * subradar latitude */
	double dd_zmax, rd_xsec, opt_brtns, cs_sb_rdr_lat;
	vps_set_four_parameters_krnl<<<1,1>>>(ddat);
	checkErrorAfterKernelLaunch("vps_set_four_parameters, line ");
	gpuErrchk(cudaMemcpyFromSymbol(&dd_zmax, st_deldop_zmax,
			sizeof(double), 0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&rd_xsec, st_rad_xsec,
			sizeof(double), 0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&opt_brtns, st_opt_brightness,
			sizeof(double), 0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&cs_sb_rdr_lat, st_cos_subradarlat,
			sizeof(double), 0, cudaMemcpyDeviceToHost));

	*deldop_zmax = dd_zmax;
	*rad_xsec = rd_xsec;
	*opt_brightness = opt_brtns;
	*cos_subradarlat = cs_sb_rdr_lat;
}

__global__ void vpst_delay_params_krnl(struct par_t *dpar, struct dat_t *ddat,
		struct pos_t **pos, int *compute_xsec, int *posn, int *ndel, int *ndop,
		int s, int size) {
	/* nframes-threaded kernel */
	int f = blockIdx.x * blockDim.x + threadIdx.x;

	if (f < size) {
		compute_xsec[f] = (dpar->vary_radalb != VARY_NONE
				&& ddat->set[s].desc.deldop.frame[f].cal.state == 'c');
		pos[f] = &ddat->set[s].desc.deldop.frame[f].pos;
		posn[f] = pos[f]->n;
		ndel[f] = ddat->set[s].desc.deldop.frame[f].ndel;
		ndop[f] = ddat->set[s].desc.deldop.frame[f].ndop;
	}
}
__global__ void vpst_dop_params_krnl(struct par_t *dpar, struct dat_t *ddat,
		struct pos_t **pos, int *compute_xsec, int *posn, int *ndop, int s,
		int size) {
	/* nframes-threaded kernel */
	int f = blockIdx.x * blockDim.x + threadIdx.x;

	if (f < size) {
		compute_xsec[f] = (dpar->vary_radalb != VARY_NONE &&
				ddat->set[s].desc.doppler.frame[f].cal.state == 'c');
		pos[f] = &ddat->set[s].desc.doppler.frame[f].pos;
		posn[f] = pos[f]->n;
		ndop[f] = ddat->set[s].desc.doppler.frame[f].ndop;
	}
}
__global__ void vpst_lghtcrv_params_krnl(struct par_t *dpar, struct dat_t *ddat,
		struct pos_t **pos, int *posn, int *bistatic, int s, int size) {
	/* nframes-threaded kernel */
	int f = blockIdx.x * blockDim.x + threadIdx.x;

	if (f < size) {
		pos[f] = &ddat->set[s].desc.lghtcrv.rend[f].pos;
		posn[f] = pos[f]->n;
		bistatic[f] = pos[f]->bistatic;
	}
}
__global__ void vpst_set_posmtrx_streams_krnl(struct dat_t *ddat,
		struct pos_t **pos,
		unsigned char type,
		int s,
		int f) {
	/* 9-threaded kernel, streamed */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % 3;
	int j = offset / 3;
	int indx;

	if (offset < 9) {
		switch (type) {
		case DELAY:
			indx = ddat->set[s].desc.deldop.v0;
			pos[f]->ae[i][j] = ddat->set[s].desc.deldop.frame[f].view[indx].ae[i][j];
			pos[f]->oe[i][j] = ddat->set[s].desc.deldop.frame[f].view[indx].oe[i][j];
			if (i==1 && j==1)
				pos[f]->bistatic = 0;
			break;
		case DOPPLER:
			indx = ddat->set[s].desc.doppler.v0;
			pos[f]->ae[i][j] =	ddat->set[s].desc.doppler.frame[f].view[indx].ae[i][j];
			pos[f]->oe[i][j] =	ddat->set[s].desc.doppler.frame[f].view[indx].oe[i][j];
			if (i==1 && j==1)
				pos[f]->bistatic = 0;
			break;
		case LGHTCRV:
			pos[f]->ae[i][j] = ddat->set[s].desc.lghtcrv.rend[f].ae[i][j];
			pos[f]->oe[i][j] = ddat->set[s].desc.lghtcrv.rend[f].oe[i][j];
			pos[f]->se[i][j] = ddat->set[s].desc.lghtcrv.rend[f].se[i][j];
			if (i==1 && j==1)
				pos[f]->bistatic = 1;
			break;
		}
	}
}

__host__ void vary_params_cuda_streams2(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, int action, double *deldop_zmax, double
		*rad_xsec, double *opt_brightness, double *cos_subradarlat, int nsets)
{
	/* This version creates its own streams and is less efficient. It also
	 * gets all parameters needed separately.
	 * Inputs:
	 * int action 				- dpar->action
	 * double *deldop_zmax 		- &deldop_zmax_save
	 * double *rad_xsec 		- &rad_xsec_save
	 * double *optbrightness	- &opt_brightness_save
	 * double *cos_subradarlat	- &cos_subradarlat_save
	 * int nsets 				- ddat->nsets
	 */

	float3 orbit_offset;
	int c=0, f, s, *compute_brightness, *compute_zmax, *bistatic, bistatic_all,
			*compute_cosdelta, *compute_xsec, n, ncalc, nx, lghtcrv_bistatic,
			nfrm_alloc, nf;
	dim3 pxBLK,THD,BLKncalc,THD9;
	THD.x = maxThreadsPerBlock;	THD9.x = 9;
	unsigned char *type;
	unsigned char htype[nsets+1];	/* the +1 is safety margin */
	int nThreads, *posn, *ndel, *ndop, *lc_n, *nframes;
	int hnframes[nsets+1], hlc_n[nsets+1], hcomp_cosdelta[nsets],
	hcomp_zmax[nsets+1], hcomp_brightness[nsets+1];

	float zmax;
	double weight, *pixels_per_km, *lghtcrv_y;
	double3 *so;
	struct pos_t **pos;
	struct vertices_t **verts;

	double orbt_off[3] = {0.0, 0.0, 0.0};

	/* Initialize */
	orbit_offset.x = orbit_offset.y = orbit_offset.z = 0.0;

	/* Allocate memory for the set-wide parameters and verts */
	gpuErrchk(cudaMalloc((void**)&type, sizeof(unsigned char) * (nsets+1)));
	gpuErrchk(cudaMalloc((void**)&nframes, sizeof(int) * (nsets+1)));
	gpuErrchk(cudaMalloc((void**)&lc_n, sizeof(int) * (nsets+1)));
	gpuErrchk(cudaMalloc((void**)&verts, sizeof(struct vertices_t*) * 2));
	gpuErrchk(cudaMalloc((void**)&compute_brightness, sizeof(int) * (nsets+1)));
	gpuErrchk(cudaMalloc((void**)&compute_zmax, sizeof(int) * (nsets+1)));
	gpuErrchk(cudaMalloc((void**)&compute_cosdelta, sizeof(int) * (nsets+1)));

	/* Initialize the device file-scope variables */
	vpst_init_krnl2<<<1,1>>>(dpar, dmod, ddat, verts, type, nframes, lc_n,
			compute_zmax, compute_cosdelta, compute_brightness, nsets, c);
	checkErrorAfterKernelLaunch("vpst_init_krnl");

	/* Now create host copies of type, nframes, verts->nf, and lc_n */
	gpuErrchk(cudaMemcpy(&htype, type, sizeof(unsigned char)*(nsets+1),
			cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(&hnframes, nframes, sizeof(int)*(nsets+1),
			cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(&hlc_n, lc_n, sizeof(int)*(nsets+1),
			cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&nf, vpst_nf, sizeof(int), 0,
			cudaMemcpyDeviceToHost));

	/* Process each dataset in turn */
	for (s=0; s<nsets; s++) {

		/* Get the allocation right as the indices for lghtcrv start with 1
		 * instead of 0 everywhere else. To not run into problems at the end
		 * or start of the array, we allocate one more than strictly needed */
		if (htype[s]==LGHTCRV) {
			nfrm_alloc = hnframes[s] + 1;
			lghtcrv_y = (double *) malloc(nfrm_alloc*sizeof(double));
		} else
			nfrm_alloc = hnframes[s];

		int hcomp_xsec[nfrm_alloc], npxls[nfrm_alloc], ddsize[nfrm_alloc],
			hndop[nfrm_alloc], hndel[nfrm_alloc], lc_xspan[nfrm_alloc], *outbndarr,
			hposn[nfrm_alloc], hbistatic[nfrm_alloc];
		float xsec[nfrm_alloc];	/* This should really only ever by nframes */
		dim3 BLK[nfrm_alloc], ddBLK[nfrm_alloc];
		cudaStream_t vp_stream[nfrm_alloc];
		gpuErrchk(cudaMalloc((void**)&pos, 		   sizeof(struct pos_t*) * nfrm_alloc));
		gpuErrchk(cudaMalloc((void**)&posn, 	   sizeof(int) * nfrm_alloc));
		gpuErrchk(cudaMalloc((void**)&ndel, 	   sizeof(int) * nfrm_alloc));
		gpuErrchk(cudaMalloc((void**)&ndop, 	   sizeof(int) * nfrm_alloc));
		gpuErrchk(cudaMalloc((void**)&outbndarr,   sizeof(int) * nfrm_alloc));
		gpuErrchk(cudaMalloc((void**)&compute_xsec,sizeof(int) * nfrm_alloc));

		/* Set up initial kernel launch parameter */
		BLK[0].x = floor((THD.x - 1 + nfrm_alloc) / THD.x);

		switch (htype[s]) {
		case DELAY:
			/* Get computation flags */
			gpuErrchk(cudaMemcpy(&hcomp_cosdelta, compute_cosdelta,
					sizeof(int)*(nsets+1), cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpy(&hcomp_zmax, compute_zmax,
					sizeof(int)*(nsets+1), cudaMemcpyDeviceToHost));

			/* Launch nframes-threaded kernel to get all relevant parameters */
			vpst_delay_params_krnl<<<BLK[0],THD>>>(dpar, ddat, pos,
					compute_xsec, posn, ndel, ndop, s, nfrm_alloc);
			checkErrorAfterKernelLaunch("vpst_delay_params_krnl");
			gpuErrchk(cudaMemcpy(&hposn, posn, sizeof(int)*nfrm_alloc,
					cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpy(&hcomp_xsec, compute_xsec,
					sizeof(int)*nfrm_alloc, cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpy(&hndop, ndop, sizeof(int)*nfrm_alloc,
					cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpy(&hndel, ndel, sizeof(int)*nfrm_alloc,
					cudaMemcpyDeviceToHost));

			/* Create streams and calculate parameters*/
			for (f=0; f<nfrm_alloc; f++) {
				cudaStreamCreate(&vp_stream[f]);
				npxls[f] = (2*hposn[f] + 1)*(2*hposn[f] + 1);
				BLK[f].x = floor((THD.x - 1 + npxls[f]) / THD.x);
				ddsize[f]= hndel[f] * hndop[f];
				ddBLK[f] = floor((THD.x -1 + ddsize[f]) / THD.x);
			}

			/* Launch posclr_streams_krnl to initialize POS view */
			for (f=0; f<nfrm_alloc; f++) {
				/* Start the if block for computing zmax and/or cross-section */
				if (hcomp_zmax[s] || hcomp_xsec[f]) {
					vpst_set_posmtrx_streams_krnl<<<1,THD9,0,vp_stream[f]>>>(ddat, pos, htype[s],
							s, f);
					posclr_streams_krnl<<<BLK[f],THD, 0, vp_stream[f]>>>(pos,posn,f);
				}
			} checkErrorAfterKernelLaunch("vpst_set_posmtrx_streams_krnl and "
					"posclr_streams_krnl (Delay-Doppler)");

			/* Determine which POS pixels cover the target, and get distance
			 * toward Earth of each POS pixel. Pass the frame streams, too. */
			posvis_cuda_streams2(dpar, dmod, ddat, pos, verts, orbit_offset,
					posn, outbndarr, s,	hnframes[s], 0, nf, 0, c, htype[s],
					vp_stream);

			for (f=0; f<nfrm_alloc; f++) {
				if (hcomp_zmax[s] || hcomp_xsec[f]) {
					/* Zero out the fit delay-Doppler image and call pos2deldop
					 * to create the fit image by mapping power from the plane
					 * of sky to delay-Doppler space.    				  */
					clrvect_streams_krnl<<<ddBLK[f],THD, 0, vp_stream[f]>>>(ddat,
							ddsize[f], s, f);
				}/* End frames loop again to call pos2deldop streams version */
			} checkErrorAfterKernelLaunch("clrvect_streams_krnl");

			/* Call the CUDA pos2deldop function */
			pos2deldop_cuda_streams(dpar, dmod, ddat, pos, ndel, ndop,
					0.0, 0.0, 0.0, 0, s, hnframes[s], 0, outbndarr, vp_stream);

			for (f=0; f<nfrm_alloc; f++) {
				if (hcomp_zmax[s] || hcomp_xsec[f]) {
					/* Compute distance toward Earth of the subradar point  */
					if (hcomp_zmax[s]) {
						zmax = 0.0;
						/* Call parallel reduction function to quickly find the
						 * distance toward Earth of the subradar point */
						zmax = compute_pos_zmax(ddat, npxls[f], s, f);
						zmax_finalize_streams_krnl<<<1,1,0,vp_stream[f]>>>(ddat, zmax, s, f);
					}
				}
			}
			for (f=0; f<nfrm_alloc; f++) {
				/*  Compute cross section  */
				if (hcomp_xsec[f]) {
					xsec[f] = 0.0;
					xsec[f] = compute_deldop_xsec_snglkrnl(ddat, hndel[f], hndop[f], s, f);
					compute_xsec_final_streams_krnl<<<1,1,0,vp_stream[f]>>>(ddat, xsec[f], s, f);
				}
			}
			checkErrorAfterKernelLaunch("compute_xsec_final_streams_krnl");

			if (hcomp_cosdelta[s])
				for (f=0; f<nfrm_alloc; f++)
					compute_cosdelta_streams_krnl<<<1,1>>>(ddat, s, f);
			checkErrorAfterKernelLaunch("compute_cosdelta_streams_krnl");

			for (f=0; f<nfrm_alloc; f++)
				cudaStreamDestroy(vp_stream[f]);

			break;
		case DOPPLER:
			/* Get computation flags */
			gpuErrchk(cudaMemcpy(&hcomp_cosdelta, compute_cosdelta,
					sizeof(int)*(nsets+1), cudaMemcpyDeviceToHost));

			/* Launch nframes-threaded kernel to get all relevant parameters */
			vpst_dop_params_krnl<<<BLK[0],THD>>>(dpar, ddat, pos,
					compute_xsec, posn, ndop, s, nfrm_alloc);
			checkErrorAfterKernelLaunch("vpst_dop_params_krnl");
			gpuErrchk(cudaMemcpy(&hposn, posn, sizeof(int)*nfrm_alloc,
					cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpy(&hcomp_xsec, compute_xsec,
					sizeof(int)*nfrm_alloc, cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpy(&hndop, ndop, sizeof(int)*nfrm_alloc,
					cudaMemcpyDeviceToHost));

			/* Calculate launch parameters and create streams */
			for (f=0; f<nfrm_alloc; f++) {
				cudaStreamCreate(&vp_stream[f]);
				npxls[f] = (2*hposn[f] + 1)*(2*hposn[f] + 1);
				BLK[f].x = floor((THD.x - 1 + npxls[f]) / THD.x);
				ddBLK[f] = floor((THD.x -1 + hndop[f]) / THD.x);
			}

			/* Launch posclr_streams_krnl to initialize POS view */
			for (f=0; f<nfrm_alloc; f++) {
				/* Start the if block for computing zmax and/or cross-section */
				if (hcomp_xsec[f]) {
					vpst_set_posmtrx_streams_krnl<<<1,THD9,0,vp_stream[f]>>>(ddat, pos, htype[s],
							s, f);
					posclr_streams_krnl<<<BLK[f],THD,0,vp_stream[f]>>>(pos,posn,f);
				}
			} checkErrorAfterKernelLaunch("vpst_set_posmtrx_streams_krnl and "
					"posclr_streams_krnl (Doppler)");

			/* Determine which POS pixels cover the target, and get distance
			 * toward Earth of each POS pixel. Pass the frame streams, too. */
			posvis_cuda_streams2(dpar, dmod, ddat, pos, verts, orbit_offset,
					posn, outbndarr, s,	hnframes[s], 0, nf, 0, c, htype[s],
					vp_stream);

			for (f=0; f<nfrm_alloc; f++) {
				if (hcomp_xsec[f]) {
					/* Zero out the fit delay-Doppler image and call pos2deldop
					 * to create the fit image by mapping power from the plane
					 * of sky to delay-Doppler space.    				      */

					clrvect_streams_krnl<<<ddBLK[f],THD, 0, vp_stream[f]>>>(ddat,
							hndop[f], s, f);
					/* End frames loop again to call pos2deldop streams version */
				}
			} checkErrorAfterKernelLaunch("clrvect_streams_krnl");

			pos2doppler_cuda_streams(dpar, dmod, ddat, pos, 0.0, 0.0, 0.0,
					ndop, 0, s, hnframes[s], 0, outbndarr, vp_stream);

			/* Calculate the Doppler cross-section if applicable */
			for (f=0; f<nfrm_alloc; f++) {
				if (hcomp_xsec[f]) {
					/* Compute cross section */
					xsec[f]=0.0;
					xsec[f] = compute_doppler_xsec(ddat, hndop[f], s, f);
				}
			}
			/* Finalize the xsec calculations and calculate cosdelta if specified */
			for (f=0; f<nfrm_alloc; f++) {
				if (hcomp_xsec[f])
					compute_xsec_final_streams_krnl<<<1,1,0,vp_stream[f]>>>(ddat, xsec[f], s, f);
				if (compute_cosdelta)
					compute_cosdelta_streams_krnl<<<1,1,0,vp_stream[f]>>>(ddat, s, f);
			}
			break;
		case POS:
			break;
		case LGHTCRV:
			/* Figure out the compute_brightness flag first */
			gpuErrchk(cudaMemcpy(&hcomp_brightness, compute_brightness,
					sizeof(int)*(nsets+1), cudaMemcpyDeviceToHost));

			if (hcomp_brightness[s]) {
				cudaCalloc1((void**)&so, sizeof(double3), (nfrm_alloc*3));
				cudaCalloc1((void**)&pixels_per_km, sizeof(int), nfrm_alloc);

				/* Launch nframes-threaded kernel to get all relevant parameters */
				vpst_lghtcrv_params_krnl<<<BLK[0],THD>>>(dpar, ddat, pos,
						posn, bistatic, s, nfrm_alloc);
				checkErrorAfterKernelLaunch("vpst_lghtcrv_params_krnl");
				gpuErrchk(cudaMemcpy(&hposn, posn, sizeof(int)*nfrm_alloc,
						cudaMemcpyDeviceToHost));
				gpuErrchk(cudaMemcpy(&hbistatic, bistatic, sizeof(int)*nfrm_alloc,
						cudaMemcpyDeviceToHost));

				/* Calculate launch parameters and create streams */
				for (f=1; f<nfrm_alloc; f++) {
					cudaStreamCreate(&vp_stream[f]);
					lc_xspan[f] = 2*posn[f] + 1;
					npxls[f] = (2*posn[f]+1)*(2*posn[f]+1);
					BLK[f].x = floor((THD.x - 1 + npxls[f]) / THD.x);
				}

				/* Clear out the plane-of-sky first */
				for (f=1; f<nfrm_alloc; f++)
					vpst_set_posmtrx_streams_krnl<<<1,THD9,0,vp_stream[f-1]>>>(
							ddat, pos, htype[s], s, f);

					posclr_streams_krnl<<<BLK[f],THD, 0, vp_stream[f-1]>>>(pos,posn,f);
				checkErrorAfterKernelLaunch("vpst_set_posmtrx_streams_krnl and"
						"posclr_streams_krnl (Lightcurve) ");

				/* Determine which POS pixels cover the target */
				posvis_cuda_streams2(dpar, dmod, ddat, pos, verts, orbit_offset,
						posn, outbndarr, s,	hnframes[s], 0, nf, 0, c, htype[s],
						vp_stream);

				/* Now view the model from the source (sun) and get the facet
				 * number and distance toward the source of each pixel in this
				 * projected view; use this information to determine which POS
				 * pixels are shadowed */
				for (f=1; f<nfrm_alloc; f++)
					if (hbistatic[f])
						bistatic_all = 1;
				if (bistatic_all)
					posvis_cuda_streams2(dpar, dmod, ddat, pos, verts,
							orbit_offset, posn, outbndarr, s, hnframes[s], 1,
							nf, 0, c, htype[s],	vp_stream);

				if (bistatic_all) {
					for (f=1; f<nfrm_alloc; f++) {
						/* Initialize this stream for the posmask kernel to follow */
						posmask_init_streams_krnl<<<1,1,0,vp_stream[f-1]>>>(pos,
								so, pixels_per_km, f);

						/* Now call posmask kernel for this stream, then loop
						 * to next stream and repeat 	 */
						posmask_streams_krnl<<<BLK[f],THD,0,vp_stream[f-1]>>>(
								dpar, pos, so, pixels_per_km, posn, npxls[f],
								lc_xspan[f], f);

					} checkErrorAfterKernelLaunch("posmask_streams_ krnl");
				}

				/* Compute model brightness for this lightcurve point */
				/* lghtcrv->y[ncalc]: calculated points for interpolation,
				 * ncalc-points total 					 */
				for (f=1; f<nfrm_alloc; f++)
					lghtcrv_y[f] = apply_photo_cuda(dmod, ddat, 0, s, f);

				/* Now launch a kernel to copy it over to the actual lghtcrv */
				BLKncalc.x = floor((THD.x - 1 + hnframes[s]) / THD.x);
				lghtcrv_copy_y_streams_krnl<<<BLKncalc,THD>>>(ddat, lghtcrv_y, s, nfrm_alloc);
				checkErrorAfterKernelLaunch("lghtcrv_copy_y_streams_krnl");

				/* Now that we have calculated the model lightcurve brightnesses
				 * y at each of the epochs x, we use cubic spline interpolation
				 * (Numerical Recipes routines spline and splint) to get model
				 * lightcurve brightness fit[i] at each OBSERVATION epoch t[i],
				 * with i=1,2,...,n.  This will allow us (in routine chi2) to
				 * compare model to data (fit[i] to obs[i]) to get chi squared.
				 * Note that vector y2 contains the second derivatives of the
				 * interpolating function at the calculation epochs x. */

				/* First make a pointer for u and cudaMalloc device memory for it */
				double *u;
				gpuErrchk(cudaMalloc((void**)&u, sizeof(double) * n));

				lghtcrv_spline_streams_krnl<<<BLKncalc,THD>>>(ddat, s, 2.0e30, 2.0e30, u);
				checkErrorAfterKernelLaunch("lghtcrv_spline_streams_krnl");

				/* Change launch parameters from ncalc threads to n threads */
				BLKncalc.x = floor((THD.x - 1 + n) / THD.x);
				lghtcrv_splint_streams2_krnl<<<BLKncalc,THD>>>(ddat, s);
				checkErrorAfterKernelLaunch("lghtcrv_splint_streams_krnl");
				/* Cleanup */
				cudaFree(u);
				cudaFree(so);
				cudaFree(pixels_per_km);
				free(lghtcrv_y);
				/* Destroy streams and free memory */
				for (f=0; f<hnframes[s]; f++)
					cudaStreamDestroy(vp_stream[f]);
			}
			break;
		default:
			bailout("vary_params.c: can't handle this dataset type yet\n");
		}



		cudaFree(pos);
		cudaFree(posn);
		cudaFree(ndel);
		cudaFree(ndop);
		cudaFree(outbndarr);
		cudaFree(compute_xsec);
	}

	/* Calculate the zmax, radar cross-section, optical brightness, and cosine
	 * subradar latitude */
	double dd_zmax, rd_xsec, opt_brtns, cs_sb_rdr_lat;
	vps_set_four_parameters_krnl<<<1,1>>>(ddat);
	checkErrorAfterKernelLaunch("vps_set_four_parameters, line ");
	gpuErrchk(cudaMemcpyFromSymbol(&dd_zmax, st_deldop_zmax,
			sizeof(double), 0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&rd_xsec, st_rad_xsec,
			sizeof(double), 0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&opt_brtns, st_opt_brightness,
			sizeof(double), 0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&cs_sb_rdr_lat, st_cos_subradarlat,
			sizeof(double), 0, cudaMemcpyDeviceToHost));

	*deldop_zmax = dd_zmax;
	*rad_xsec = rd_xsec;
	*opt_brightness = opt_brtns;
	*cos_subradarlat = cs_sb_rdr_lat;

	cudaFree(type);
	cudaFree(nframes);
	cudaFree(lc_n);
	cudaFree(verts);
	cudaFree(compute_brightness);
	cudaFree(compute_zmax);
	cudaFree(compute_cosdelta);
}

__host__ void vary_params_cuda_streams3(
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
		struct vertices_t **verts,
		unsigned char *htype,
		unsigned char *dtype,
		int nf,
		int nsets,
		cudaStream_t *vp_stream)
{
	/* This third iteration uses streams that are passed as argument.
	 * It also does not calculate/copy the various parameters but accepts
	 * them as arguments. Many doubles are floats or CUDA internal types like
	 * float3 or int4.
	 * Inputs:
	 * int action 				- dpar->action
	 * double *deldop_zmax 		- &deldop_zmax_save
	 * double *rad_xsec 		- &rad_xsec_save
	 * double *optbrightness	- &opt_brightness_save
	 * double *cos_subradarlat	- &cos_subradarlat_save
	 * int nsets 				- ddat->nsets
	 */

	float3 orbit_offset;
	int c=0, f, s, *compute_brightness, *compute_zmax, *bistatic, bistatic_all,
			*compute_cosdelta, *compute_xsec, n, ncalc, nx, lghtcrv_bistatic,
			nfrm_alloc, nThreads, *posn, *ndel, *ndop;
	dim3 pxBLK,THD,BLKncalc,THD9;
	THD.x = maxThreadsPerBlock;	THD9.x = 9;
	int hcomp_cosdelta[nsets], hcomp_zmax[nsets+1], hcomp_brightness[nsets+1];

	float zmax;
	double weight, *pixels_per_km, *lghtcrv_y;
	double3 *so;
	struct pos_t **pos;
	double orbt_off[3] = {0.0, 0.0, 0.0};

	/* Initialize */
	orbit_offset.x = orbit_offset.y = orbit_offset.z = 0.0;

	/* Allocate memory for the set-wide parameters and verts */
	gpuErrchk(cudaMalloc((void**)&compute_brightness, sizeof(int) * (nsets+1)));
	gpuErrchk(cudaMalloc((void**)&compute_zmax, sizeof(int) * (nsets+1)));
	gpuErrchk(cudaMalloc((void**)&compute_cosdelta, sizeof(int) * (nsets+1)));

	/* Initialize the device file-scope variables */
	vpst_init_krnl3<<<1,1>>>(dpar, ddat, compute_zmax, compute_cosdelta,
			compute_brightness, dtype,nsets);
	checkErrorAfterKernelLaunch("vpst_init_krnl3");


	/* Process each dataset in turn */
	for (s=0; s<nsets; s++) {

		/* Get the allocation right as the indices for lghtcrv start with 1
		 * instead of 0 everywhere else. To not run into problems at the end
		 * or start of the array, we allocate one more than strictly needed */
		if (htype[s]==LGHTCRV) {
			nfrm_alloc = hnframes[s] + 1;
			lghtcrv_y = (double *) malloc(nfrm_alloc*sizeof(double));
		} else
			nfrm_alloc = hnframes[s];

		int hcomp_xsec[nfrm_alloc], npxls[nfrm_alloc], ddsize[nfrm_alloc],
			hndop[nfrm_alloc], hndel[nfrm_alloc], lc_xspan[nfrm_alloc], *outbndarr,
			hposn[nfrm_alloc], hbistatic[nfrm_alloc];
		float xsec[nfrm_alloc];	/* This should really only ever by nframes */
		dim3 BLK[nfrm_alloc], ddBLK[nfrm_alloc];
		gpuErrchk(cudaMalloc((void**)&pos, 		   sizeof(struct pos_t*) * nfrm_alloc));
		gpuErrchk(cudaMalloc((void**)&posn, 	   sizeof(int) * nfrm_alloc));
		gpuErrchk(cudaMalloc((void**)&ndel, 	   sizeof(int) * nfrm_alloc));
		gpuErrchk(cudaMalloc((void**)&ndop, 	   sizeof(int) * nfrm_alloc));
		gpuErrchk(cudaMalloc((void**)&outbndarr,   sizeof(int) * nfrm_alloc));
		gpuErrchk(cudaMalloc((void**)&compute_xsec,sizeof(int) * nfrm_alloc));

		/* Set up initial kernel launch parameter */
		BLK[0].x = floor((THD.x - 1 + nfrm_alloc) / THD.x);

		switch (htype[s]) {
		case DELAY:
			/* Get computation flags */
			gpuErrchk(cudaMemcpy(&hcomp_cosdelta, compute_cosdelta,
					sizeof(int)*(nsets+1), cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpy(&hcomp_zmax, compute_zmax,
					sizeof(int)*(nsets+1), cudaMemcpyDeviceToHost));

			/* Launch nframes-threaded kernel to get all relevant parameters */
			vpst_delay_params_krnl<<<BLK[0],THD>>>(dpar, ddat, pos,
					compute_xsec, posn, ndel, ndop, s, nfrm_alloc);
			checkErrorAfterKernelLaunch("vpst_delay_params_krnl");
			gpuErrchk(cudaMemcpy(&hposn, posn, sizeof(int)*nfrm_alloc,
					cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpy(&hcomp_xsec, compute_xsec,
					sizeof(int)*nfrm_alloc, cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpy(&hndop, ndop, sizeof(int)*nfrm_alloc,
					cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpy(&hndel, ndel, sizeof(int)*nfrm_alloc,
					cudaMemcpyDeviceToHost));

			/* Create streams and calculate parameters*/
			for (f=0; f<nfrm_alloc; f++) {
				npxls[f] = (2*hposn[f] + 1)*(2*hposn[f] + 1);
				BLK[f].x = floor((THD.x - 1 + npxls[f]) / THD.x);
				ddsize[f]= hndel[f] * hndop[f];
				ddBLK[f] = floor((THD.x -1 + ddsize[f]) / THD.x);
			}

			/* Launch posclr_streams_krnl to initialize POS view */
			for (f=0; f<nfrm_alloc; f++) {
				/* Start the if block for computing zmax and/or cross-section */
				if (hcomp_zmax[s] || hcomp_xsec[f]) {
					vpst_set_posmtrx_streams_krnl<<<1,THD9,0,vp_stream[f]>>>(ddat, pos, htype[s],
							s, f);
					posclr_streams_krnl<<<BLK[f],THD, 0, vp_stream[f]>>>(pos,posn,f);
				}
			} checkErrorAfterKernelLaunch("vpst_set_posmtrx_streams_krnl and "
					"posclr_streams_krnl (Delay-Doppler)");

			/* Determine which POS pixels cover the target, and get distance
			 * toward Earth of each POS pixel. Pass the frame streams, too. */
			posvis_cuda_streams2(dpar, dmod, ddat, pos, verts, orbit_offset,
					hposn, outbndarr, s, hnframes[s], 0, nf, 0, c, htype[s],
					vp_stream);

			for (f=0; f<nfrm_alloc; f++) {
				if (hcomp_zmax[s] || hcomp_xsec[f]) {
					/* Zero out the fit delay-Doppler image and call pos2deldop
					 * to create the fit image by mapping power from the plane
					 * of sky to delay-Doppler space.    				  */
					clrvect_streams_krnl<<<ddBLK[f],THD, 0, vp_stream[f]>>>(ddat,
							ddsize[f], s, f);
				}/* End frames loop again to call pos2deldop streams version */
			} checkErrorAfterKernelLaunch("clrvect_streams_krnl");

			/* Call the CUDA pos2deldop function */
			pos2deldop_cuda_streams(dpar, dmod, ddat, pos, ndel, ndop,
					0.0, 0.0, 0.0, 0, s, hnframes[s], 0, outbndarr, vp_stream);

			for (f=0; f<nfrm_alloc; f++) {
				if (hcomp_zmax[s] || hcomp_xsec[f]) {
					/* Compute distance toward Earth of the subradar point  */
					if (hcomp_zmax[s]) {
						zmax = 0.0;
						/* Call parallel reduction function to quickly find the
						 * distance toward Earth of the subradar point */
						zmax = compute_pos_zmax(ddat, npxls[f], s, f);
						zmax_finalize_streams_krnl<<<1,1,0,vp_stream[f]>>>(ddat, zmax, s, f);
					}
				}
			}
			for (f=0; f<nfrm_alloc; f++) {
				/*  Compute cross section  */
				if (hcomp_xsec[f]) {
					xsec[f] = 0.0;
					xsec[f] = compute_deldop_xsec_snglkrnl(ddat, hndel[f], hndop[f], s, f);
					compute_xsec_final_streams_krnl<<<1,1,0,vp_stream[f]>>>(ddat, xsec[f], s, f);
				}
			}
			checkErrorAfterKernelLaunch("compute_xsec_final_streams_krnl");

			if (hcomp_cosdelta[s])
				for (f=0; f<nfrm_alloc; f++)
					compute_cosdelta_streams_krnl<<<1,1>>>(ddat, s, f);
			checkErrorAfterKernelLaunch("compute_cosdelta_streams_krnl");

			break;
		case DOPPLER:
			/* Get computation flags */
			gpuErrchk(cudaMemcpy(&hcomp_cosdelta, compute_cosdelta,
					sizeof(int)*(nsets+1), cudaMemcpyDeviceToHost));

			/* Launch nframes-threaded kernel to get all relevant parameters */
			vpst_dop_params_krnl<<<BLK[0],THD>>>(dpar, ddat, pos,
					compute_xsec, posn, ndop, s, nfrm_alloc);
			checkErrorAfterKernelLaunch("vpst_dop_params_krnl");
			gpuErrchk(cudaMemcpy(&hposn, posn, sizeof(int)*nfrm_alloc,
					cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpy(&hcomp_xsec, compute_xsec,
					sizeof(int)*nfrm_alloc, cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpy(&hndop, ndop, sizeof(int)*nfrm_alloc,
					cudaMemcpyDeviceToHost));

			/* Calculate launch parameters and create streams */
			for (f=0; f<nfrm_alloc; f++) {
				npxls[f] = (2*hposn[f] + 1)*(2*hposn[f] + 1);
				BLK[f].x = floor((THD.x - 1 + npxls[f]) / THD.x);
				ddBLK[f] = floor((THD.x -1 + hndop[f]) / THD.x);
			}

			/* Launch posclr_streams_krnl to initialize POS view */
			for (f=0; f<nfrm_alloc; f++) {
				/* Start the if block for computing zmax and/or cross-section */
				if (hcomp_xsec[f]) {
					vpst_set_posmtrx_streams_krnl<<<1,THD9,0,vp_stream[f]>>>(ddat, pos, htype[s],
							s, f);
					posclr_streams_krnl<<<BLK[f],THD,0,vp_stream[f]>>>(pos,posn,f);
				}
			} checkErrorAfterKernelLaunch("vpst_set_posmtrx_streams_krnl and "
					"posclr_streams_krnl (Doppler)");

			/* Determine which POS pixels cover the target, and get distance
			 * toward Earth of each POS pixel. Pass the frame streams, too. */
			posvis_cuda_streams2(dpar, dmod, ddat, pos, verts, orbit_offset,
					hposn, outbndarr, s,	hnframes[s], 0, nf, 0, c, htype[s],
					vp_stream);

			for (f=0; f<nfrm_alloc; f++) {
				if (hcomp_xsec[f]) {
					/* Zero out the fit delay-Doppler image and call pos2deldop
					 * to create the fit image by mapping power from the plane
					 * of sky to delay-Doppler space.    				      */

					clrvect_streams_krnl<<<ddBLK[f],THD, 0, vp_stream[f]>>>(ddat,
							hndop[f], s, f);
					/* End frames loop again to call pos2deldop streams version */
				}
			} checkErrorAfterKernelLaunch("clrvect_streams_krnl");

			pos2doppler_cuda_streams(dpar, dmod, ddat, pos, 0.0, 0.0, 0.0,
					ndop, 0, s, hnframes[s], 0, outbndarr, vp_stream);

			/* Calculate the Doppler cross-section if applicable */
			for (f=0; f<nfrm_alloc; f++) {
				if (hcomp_xsec[f]) {
					/* Compute cross section */
					xsec[f]=0.0;
					xsec[f] = compute_doppler_xsec(ddat, hndop[f], s, f);
				}
			}
			/* Finalize the xsec calculations and calculate cosdelta if specified */
			for (f=0; f<nfrm_alloc; f++) {
				if (hcomp_xsec[f])
					compute_xsec_final_streams_krnl<<<1,1,0,vp_stream[f]>>>(ddat, xsec[f], s, f);
				if (compute_cosdelta)
					compute_cosdelta_streams_krnl<<<1,1,0,vp_stream[f]>>>(ddat, s, f);
			}
			break;
		case POS:
			break;
		case LGHTCRV:
			/* Figure out the compute_brightness flag first */
			gpuErrchk(cudaMemcpy(&hcomp_brightness, compute_brightness,
					sizeof(int)*(nsets+1), cudaMemcpyDeviceToHost));

			if (hcomp_brightness[s]) {
				cudaCalloc1((void**)&so, sizeof(double3), (nfrm_alloc*3));
				cudaCalloc1((void**)&pixels_per_km, sizeof(int), nfrm_alloc);

				/* Launch nframes-threaded kernel to get all relevant parameters */
				vpst_lghtcrv_params_krnl<<<BLK[0],THD>>>(dpar, ddat, pos,
						posn, bistatic, s, nfrm_alloc);
				checkErrorAfterKernelLaunch("vpst_lghtcrv_params_krnl");
				gpuErrchk(cudaMemcpy(&hposn, posn, sizeof(int)*nfrm_alloc,
						cudaMemcpyDeviceToHost));
				gpuErrchk(cudaMemcpy(&hbistatic, bistatic, sizeof(int)*nfrm_alloc,
						cudaMemcpyDeviceToHost));

				/* Calculate launch parameters and create streams */
				for (f=1; f<nfrm_alloc; f++) {
					lc_xspan[f] = 2*posn[f] + 1;
					npxls[f] = (2*posn[f]+1)*(2*posn[f]+1);
					BLK[f].x = floor((THD.x - 1 + npxls[f]) / THD.x);
				}

				/* Clear out the plane-of-sky first */
				for (f=1; f<nfrm_alloc; f++) {
					vpst_set_posmtrx_streams_krnl<<<1,THD9,0,vp_stream[f-1]>>>(
							ddat, pos, htype[s], s, f);

					posclr_streams_krnl<<<BLK[f],THD, 0, vp_stream[f-1]>>>(pos,posn,f);
				}
				checkErrorAfterKernelLaunch("vpst_set_posmtrx_streams_krnl and"
						"posclr_streams_krnl (Lightcurve) ");

				/* Determine which POS pixels cover the target */
				posvis_cuda_streams2(dpar, dmod, ddat, pos, verts, orbit_offset,
						hposn, outbndarr, s,	hnframes[s], 0, nf, 0, c, htype[s],
						vp_stream);

				/* Now view the model from the source (sun) and get the facet
				 * number and distance toward the source of each pixel in this
				 * projected view; use this information to determine which POS
				 * pixels are shadowed */
				for (f=1; f<nfrm_alloc; f++)
					if (hbistatic[f])
						bistatic_all = 1;
				if (bistatic_all)
					posvis_cuda_streams2(dpar, dmod, ddat, pos, verts,
							orbit_offset, hposn, outbndarr, s, hnframes[s], 1,
							nf, 0, c, htype[s],	vp_stream);

				if (bistatic_all) {
					for (f=1; f<nfrm_alloc; f++) {
						/* Initialize this stream for the posmask kernel to follow */
						posmask_init_streams_krnl<<<1,1,0,vp_stream[f-1]>>>(pos,
								so, pixels_per_km, f);

						/* Now call posmask kernel for this stream, then loop
						 * to next stream and repeat 	 */
						posmask_streams_krnl<<<BLK[f],THD,0,vp_stream[f-1]>>>(
								dpar, pos, so, pixels_per_km, posn, npxls[f],
								lc_xspan[f], f);

					} checkErrorAfterKernelLaunch("posmask_streams_ krnl");
				}

				/* Compute model brightness for this lightcurve point */
				/* lghtcrv->y[ncalc]: calculated points for interpolation,
				 * ncalc-points total 					 */
				for (f=1; f<nfrm_alloc; f++)
					lghtcrv_y[f] = apply_photo_cuda(dmod, ddat, 0, s, f);

				/* Now launch a kernel to copy it over to the actual lghtcrv */
				BLKncalc.x = floor((THD.x - 1 + hnframes[s]) / THD.x);
				lghtcrv_copy_y_streams_krnl<<<BLKncalc,THD>>>(ddat, lghtcrv_y, s, hnframes[s]);
				checkErrorAfterKernelLaunch("lghtcrv_copy_y_streams_krnl");

				/* Now that we have calculated the model lightcurve brightnesses
				 * y at each of the epochs x, we use cubic spline interpolation
				 * (Numerical Recipes routines spline and splint) to get model
				 * lightcurve brightness fit[i] at each OBSERVATION epoch t[i],
				 * with i=1,2,...,n.  This will allow us (in routine chi2) to
				 * compare model to data (fit[i] to obs[i]) to get chi squared.
				 * Note that vector y2 contains the second derivatives of the
				 * interpolating function at the calculation epochs x. */

				/* First make a pointer for u and cudaMalloc device memory for it */
				double *u;
				gpuErrchk(cudaMalloc((void**)&u, sizeof(double) * n));

				lghtcrv_spline_streams_krnl<<<BLKncalc,THD>>>(ddat, s, 2.0e30, 2.0e30, u);
				checkErrorAfterKernelLaunch("lghtcrv_spline_streams_krnl");

				/* Change launch parameters from ncalc threads to n threads */
				BLKncalc.x = floor((THD.x - 1 + n) / THD.x);
				lghtcrv_splint_streams2_krnl<<<BLKncalc,THD>>>(ddat, s);
				checkErrorAfterKernelLaunch("lghtcrv_splint_streams_krnl");
				/* Cleanup */
				cudaFree(u);
				cudaFree(so);
				cudaFree(pixels_per_km);
				free(lghtcrv_y);
				/* Destroy streams and free memory */
				for (f=0; f<hnframes[s]; f++)
					cudaStreamDestroy(vp_stream[f]);
			}
			break;
		default:
			bailout("vary_params.c: can't handle this dataset type yet\n");
		}



		cudaFree(pos);
		cudaFree(posn);
		cudaFree(ndel);
		cudaFree(ndop);
		cudaFree(outbndarr);
		cudaFree(compute_xsec);
	}

//	dbg_print_pos_z(ddat, 0, 0, 75, "Streams_pos_z_s0f0.csv");
//	dbg_print_pos_cose_s(ddat, 0, 0, 75, "Streams_pos_cose_s_s0f0.csv");

	/* Calculate the zmax, radar cross-section, optical brightness, and cosine
	 * subradar latitude */
	double dd_zmax, rd_xsec, opt_brtns, cs_sb_rdr_lat;
	vps_set_four_parameters_krnl<<<1,1>>>(ddat);
	checkErrorAfterKernelLaunch("vps_set_four_parameters, line ");
	gpuErrchk(cudaMemcpyFromSymbol(&dd_zmax, st_deldop_zmax,
			sizeof(double), 0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&rd_xsec, st_rad_xsec,
			sizeof(double), 0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&opt_brtns, st_opt_brightness,
			sizeof(double), 0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&cs_sb_rdr_lat, st_cos_subradarlat,
			sizeof(double), 0, cudaMemcpyDeviceToHost));

	*deldop_zmax = dd_zmax;
	*rad_xsec = rd_xsec;
	*opt_brightness = opt_brtns;
	*cos_subradarlat = cs_sb_rdr_lat;



	cudaFree(compute_brightness);
	cudaFree(compute_zmax);
	cudaFree(compute_cosdelta);
}

__host__ void vary_params_cuda_streams3f(
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
		struct vertices_t **verts,
		unsigned char *htype,
		unsigned char *dtype,
		int nf,
		int nsets,
		cudaStream_t *vp_stream)
{
	/* This third iteration uses streams that are passed as argument.
	 * It also does not calculate/copy the various parameters but accepts
	 * them as arguments. This "f" sub-version eliminates all remaining doubles
	 * Inputs:
	 * int action 				- dpar->action
	 * double *deldop_zmax 		- &deldop_zmax_save
	 * double *rad_xsec 		- &rad_xsec_save
	 * double *optbrightness	- &opt_brightness_save
	 * double *cos_subradarlat	- &cos_subradarlat_save
	 * int nsets 				- ddat->nsets
	 */

	float3 orbit_offset;
	int c=0, f, s, *compute_brightness, *compute_zmax, *bistatic, bistatic_all,
			*compute_cosdelta, *compute_xsec, n, ncalc, nx, lghtcrv_bistatic,
			nfrm_alloc, nThreads, *posn, *ndel, *ndop;
	dim3 pxBLK,THD,BLKncalc,THD9;
	THD.x = maxThreadsPerBlock;	THD9.x = 9;
	int hcomp_cosdelta[nsets], hcomp_zmax[nsets+1], hcomp_brightness[nsets+1];

	float zmax;
	float weight, *pixels_per_km, *lghtcrv_y;
	float3 *so, orbit_xydopoff;
	struct pos_t **pos;
	//double orbt_off[3] = {0.0, 0.0, 0.0};

	/* Initialize */
	orbit_offset.x = orbit_offset.y = orbit_offset.z = 0.0;
	orbit_xydopoff.x = orbit_xydopoff.y = orbit_xydopoff.z = 0.0;

	/* Allocate memory for the set-wide parameters and verts */
	gpuErrchk(cudaMalloc((void**)&compute_brightness, sizeof(int) * (nsets+1)));
	gpuErrchk(cudaMalloc((void**)&compute_zmax, sizeof(int) * (nsets+1)));
	gpuErrchk(cudaMalloc((void**)&compute_cosdelta, sizeof(int) * (nsets+1)));

	/* Initialize the device file-scope variables */
	vpst_init_krnl3<<<1,1>>>(dpar, ddat, compute_zmax, compute_cosdelta,
			compute_brightness, dtype,nsets);
	checkErrorAfterKernelLaunch("vpst_init_krnl3");


	/* Process each dataset in turn */
	for (s=0; s<nsets; s++) {

		/* Get the allocation right as the indices for lghtcrv start with 1
		 * instead of 0 everywhere else. To not run into problems at the end
		 * or start of the array, we allocate one more than strictly needed */
		if (htype[s]==LGHTCRV) {
			nfrm_alloc = hnframes[s] + 1;
			lghtcrv_y = (float *) malloc(nfrm_alloc*sizeof(float));
		} else
			nfrm_alloc = hnframes[s];

		int hcomp_xsec[nfrm_alloc], npxls[nfrm_alloc], ddsize[nfrm_alloc],
			hndop[nfrm_alloc], hndel[nfrm_alloc], lc_xspan[nfrm_alloc], *outbndarr,
			hposn[nfrm_alloc], hbistatic[nfrm_alloc];
		float xsec[nfrm_alloc];	/* This should really only ever by nframes */
		dim3 BLK[nfrm_alloc], ddBLK[nfrm_alloc];
		gpuErrchk(cudaMalloc((void**)&pos, 		   sizeof(struct pos_t*) * nfrm_alloc));
		gpuErrchk(cudaMalloc((void**)&posn, 	   sizeof(int) * nfrm_alloc));
		gpuErrchk(cudaMalloc((void**)&ndel, 	   sizeof(int) * nfrm_alloc));
		gpuErrchk(cudaMalloc((void**)&ndop, 	   sizeof(int) * nfrm_alloc));
		gpuErrchk(cudaMalloc((void**)&outbndarr,   sizeof(int) * nfrm_alloc));
		gpuErrchk(cudaMalloc((void**)&compute_xsec,sizeof(int) * nfrm_alloc));

		/* Set up initial kernel launch parameter */
		BLK[0].x = floor((THD.x - 1 + nfrm_alloc) / THD.x);

		switch (htype[s]) {
		case DELAY:
			/* Get computation flags */
			gpuErrchk(cudaMemcpy(&hcomp_cosdelta, compute_cosdelta,
					sizeof(int)*(nsets+1), cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpy(&hcomp_zmax, compute_zmax,
					sizeof(int)*(nsets+1), cudaMemcpyDeviceToHost));

			/* Launch nframes-threaded kernel to get all relevant parameters */
			vpst_delay_params_krnl<<<BLK[0],THD>>>(dpar, ddat, pos,
					compute_xsec, posn, ndel, ndop, s, nfrm_alloc);
			checkErrorAfterKernelLaunch("vpst_delay_params_krnl");
			gpuErrchk(cudaMemcpy(&hposn, posn, sizeof(int)*nfrm_alloc,
					cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpy(&hcomp_xsec, compute_xsec,
					sizeof(int)*nfrm_alloc, cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpy(&hndop, ndop, sizeof(int)*nfrm_alloc,
					cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpy(&hndel, ndel, sizeof(int)*nfrm_alloc,
					cudaMemcpyDeviceToHost));

			/* Create streams and calculate parameters*/
			for (f=0; f<nfrm_alloc; f++) {
				npxls[f] = (2*hposn[f] + 1)*(2*hposn[f] + 1);
				BLK[f].x = floor((THD.x - 1 + npxls[f]) / THD.x);
				ddsize[f]= hndel[f] * hndop[f];
				ddBLK[f] = floor((THD.x -1 + ddsize[f]) / THD.x);
			}

			/* Launch posclr_streams_krnl to initialize POS view */
			for (f=0; f<nfrm_alloc; f++) {
				/* Start the if block for computing zmax and/or cross-section */
				if (hcomp_zmax[s] || hcomp_xsec[f]) {
					vpst_set_posmtrx_streams_krnl<<<1,THD9,0,vp_stream[f]>>>(ddat, pos, htype[s],
							s, f);
					posclr_streams_krnl<<<BLK[f],THD, 0, vp_stream[f]>>>(pos,posn,f);
				}
			} checkErrorAfterKernelLaunch("vpst_set_posmtrx_streams_krnl and "
					"posclr_streams_krnl (Delay-Doppler)");

			/* Determine which POS pixels cover the target, and get distance
			 * toward Earth of each POS pixel. Pass the frame streams, too. */
			posvis_cuda_streams2(dpar, dmod, ddat, pos, verts, orbit_offset,
					hposn, outbndarr, s, hnframes[s], 0, nf, 0, c, htype[s],
					vp_stream);

			for (f=0; f<nfrm_alloc; f++) {
				if (hcomp_zmax[s] || hcomp_xsec[f]) {
					/* Zero out the fit delay-Doppler image and call pos2deldop
					 * to create the fit image by mapping power from the plane
					 * of sky to delay-Doppler space.    				  */
					clrvect_streams_krnl<<<ddBLK[f],THD, 0, vp_stream[f]>>>(ddat,
							ddsize[f], s, f);
				}/* End frames loop again to call pos2deldop streams version */
			} checkErrorAfterKernelLaunch("clrvect_streams_krnl");

			/* Call the CUDA pos2deldop function */
			pos2deldop_cuda_streams_f(dpar, dmod, ddat, pos, ndel, ndop,
					orbit_xydopoff,
					0,
					s,
					hnframes[s], 0, outbndarr, vp_stream);

			for (f=0; f<nfrm_alloc; f++) {
				if (hcomp_zmax[s] || hcomp_xsec[f]) {
					/* Compute distance toward Earth of the subradar point  */
					if (hcomp_zmax[s]) {
						zmax = 0.0;
						/* Call parallel reduction function to quickly find the
						 * distance toward Earth of the subradar point */
						zmax = compute_pos_zmax(ddat, npxls[f], s, f);
						zmax_finalize_streams_krnl<<<1,1,0,vp_stream[f]>>>(ddat, zmax, s, f);
					}
				}
			}
			for (f=0; f<nfrm_alloc; f++) {
				/*  Compute cross section  */
				if (hcomp_xsec[f]) {
					xsec[f] = 0.0;
					xsec[f] = compute_deldop_xsec_snglkrnl(ddat, hndel[f], hndop[f], s, f);
					compute_xsec_final_streams_krnl<<<1,1,0,vp_stream[f]>>>(ddat, xsec[f], s, f);
				}
			}
			checkErrorAfterKernelLaunch("compute_xsec_final_streams_krnl");

			if (hcomp_cosdelta[s])
				for (f=0; f<nfrm_alloc; f++)
					compute_cosdelta_streams_f_krnl<<<1,1>>>(ddat, s, f);
			checkErrorAfterKernelLaunch("compute_cosdelta_streams_krnl");

			break;
		case DOPPLER:
			/* Get computation flags */
			gpuErrchk(cudaMemcpy(&hcomp_cosdelta, compute_cosdelta,
					sizeof(int)*(nsets+1), cudaMemcpyDeviceToHost));

			/* Launch nframes-threaded kernel to get all relevant parameters */
			vpst_dop_params_krnl<<<BLK[0],THD>>>(dpar, ddat, pos,
					compute_xsec, posn, ndop, s, nfrm_alloc);
			checkErrorAfterKernelLaunch("vpst_dop_params_krnl");
			gpuErrchk(cudaMemcpy(&hposn, posn, sizeof(int)*nfrm_alloc,
					cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpy(&hcomp_xsec, compute_xsec,
					sizeof(int)*nfrm_alloc, cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpy(&hndop, ndop, sizeof(int)*nfrm_alloc,
					cudaMemcpyDeviceToHost));

			/* Calculate launch parameters and create streams */
			for (f=0; f<nfrm_alloc; f++) {
				npxls[f] = (2*hposn[f] + 1)*(2*hposn[f] + 1);
				BLK[f].x = floor((THD.x - 1 + npxls[f]) / THD.x);
				ddBLK[f] = floor((THD.x -1 + hndop[f]) / THD.x);
			}

			/* Launch posclr_streams_krnl to initialize POS view */
			for (f=0; f<nfrm_alloc; f++) {
				/* Start the if block for computing zmax and/or cross-section */
				if (hcomp_xsec[f]) {
					vpst_set_posmtrx_streams_krnl<<<1,THD9,0,vp_stream[f]>>>(ddat, pos, htype[s],
							s, f);
					posclr_streams_krnl<<<BLK[f],THD,0,vp_stream[f]>>>(pos,posn,f);
				}
			} checkErrorAfterKernelLaunch("vpst_set_posmtrx_streams_krnl and "
					"posclr_streams_krnl (Doppler)");

			/* Determine which POS pixels cover the target, and get distance
			 * toward Earth of each POS pixel. Pass the frame streams, too. */
			posvis_cuda_streams2(dpar, dmod, ddat, pos, verts, orbit_offset,
					hposn, outbndarr, s,	hnframes[s], 0, nf, 0, c, htype[s],
					vp_stream);

			for (f=0; f<nfrm_alloc; f++) {
				if (hcomp_xsec[f]) {
					/* Zero out the fit delay-Doppler image and call pos2deldop
					 * to create the fit image by mapping power from the plane
					 * of sky to delay-Doppler space.    				      */
					clrvect_streams_krnl<<<ddBLK[f],THD, 0, vp_stream[f]>>>(ddat,
							hndop[f], s, f);
					/* End frames loop again to call pos2deldop streams version */
				}
			} checkErrorAfterKernelLaunch("clrvect_streams_krnl");

			pos2doppler_cuda_streams_f(dpar, dmod, ddat, pos, orbit_xydopoff,
					ndop, 0, s, hnframes[s], 0, outbndarr, vp_stream);

			/* Calculate the Doppler cross-section if applicable */
			for (f=0; f<nfrm_alloc; f++) {
				if (hcomp_xsec[f]) {
					/* Compute cross section */
					xsec[f]=0.0;
					xsec[f] = compute_doppler_xsec(ddat, hndop[f], s, f);
				}
			}
			/* Finalize the xsec calculations and calculate cosdelta if specified */
			for (f=0; f<nfrm_alloc; f++) {
				if (hcomp_xsec[f])
					compute_xsec_final_streams_krnl<<<1,1,0,vp_stream[f]>>>(ddat, xsec[f], s, f);
				if (compute_cosdelta)
					compute_cosdelta_streams_f_krnl<<<1,1,0,vp_stream[f]>>>(ddat, s, f);
			}
			break;
		case POS:
			break;
		case LGHTCRV:
			/* Figure out the compute_brightness flag first */
			gpuErrchk(cudaMemcpy(&hcomp_brightness, compute_brightness,
					sizeof(int)*(nsets+1), cudaMemcpyDeviceToHost));

			if (hcomp_brightness[s]) {
				cudaCalloc1((void**)&so, sizeof(float3), (nfrm_alloc*3));
				cudaCalloc1((void**)&pixels_per_km, sizeof(float), nfrm_alloc);

				/* Launch nframes-threaded kernel to get all relevant parameters */
				vpst_lghtcrv_params_krnl<<<BLK[0],THD>>>(dpar, ddat, pos,
						posn, bistatic, s, nfrm_alloc);
				checkErrorAfterKernelLaunch("vpst_lghtcrv_params_krnl");
				gpuErrchk(cudaMemcpy(&hposn, posn, sizeof(int)*nfrm_alloc,
						cudaMemcpyDeviceToHost));
				gpuErrchk(cudaMemcpy(&hbistatic, bistatic, sizeof(int)*nfrm_alloc,
						cudaMemcpyDeviceToHost));

				/* Calculate launch parameters and create streams */
				for (f=1; f<nfrm_alloc; f++) {
					lc_xspan[f] = 2*posn[f] + 1;
					npxls[f] = (2*posn[f]+1)*(2*posn[f]+1);
					BLK[f].x = floor((THD.x - 1 + npxls[f]) / THD.x);
				}

				/* Clear out the plane-of-sky first */
				for (f=1; f<nfrm_alloc; f++) {
					vpst_set_posmtrx_streams_krnl<<<1,THD9,0,vp_stream[f-1]>>>(
							ddat, pos, htype[s], s, f);

					posclr_streams_krnl<<<BLK[f],THD, 0, vp_stream[f-1]>>>(pos,posn,f);
				}
				checkErrorAfterKernelLaunch("vpst_set_posmtrx_streams_krnl and"
						"posclr_streams_krnl (Lightcurve) ");

				/* Determine which POS pixels cover the target */
				posvis_cuda_streams2(dpar, dmod, ddat, pos, verts, orbit_offset,
						hposn, outbndarr, s,	hnframes[s], 0, nf, 0, c, htype[s],
						vp_stream);

				/* Now view the model from the source (sun) and get the facet
				 * number and distance toward the source of each pixel in this
				 * projected view; use this information to determine which POS
				 * pixels are shadowed */
				for (f=1; f<nfrm_alloc; f++)
					if (hbistatic[f])
						bistatic_all = 1;
				if (bistatic_all)
					posvis_cuda_streams2(dpar, dmod, ddat, pos, verts,
							orbit_offset, hposn, outbndarr, s, hnframes[s], 1,
							nf, 0, c, htype[s],	vp_stream);

				if (bistatic_all) {
					for (f=1; f<nfrm_alloc; f++) {
						/* Initialize this stream for the posmask kernel to follow */
						posmask_init_streams_f_krnl<<<1,1,0,vp_stream[f-1]>>>(pos,
								so, pixels_per_km, f);

						/* Now call posmask kernel for this stream, then loop
						 * to next stream and repeat 	 */
						posmask_streams_f_krnl<<<BLK[f],THD,0,vp_stream[f-1]>>>(
								dpar, pos, so, pixels_per_km, npxls[f],
								lc_xspan[f], f);

					} checkErrorAfterKernelLaunch("posmask_streams_ krnl");
				}

				/* Compute model brightness for this lightcurve point */
				/* lghtcrv->y[ncalc]: calculated points for interpolation,
				 * ncalc-points total 					 */
				for (f=1; f<nfrm_alloc; f++)
					lghtcrv_y[f] = apply_photo_cuda(dmod, ddat, 0, s, f);

				/* Now launch a kernel to copy it over to the actual lghtcrv */
				BLKncalc.x = floor((THD.x - 1 + hnframes[s]) / THD.x);
				lghtcrv_copy_y_streams_f_krnl<<<BLKncalc,THD>>>(ddat, lghtcrv_y, s, hnframes[s]);
				checkErrorAfterKernelLaunch("lghtcrv_copy_y_streams_krnl");

				/* Now that we have calculated the model lightcurve brightnesses
				 * y at each of the epochs x, we use cubic spline interpolation
				 * (Numerical Recipes routines spline and splint) to get model
				 * lightcurve brightness fit[i] at each OBSERVATION epoch t[i],
				 * with i=1,2,...,n.  This will allow us (in routine chi2) to
				 * compare model to data (fit[i] to obs[i]) to get chi squared.
				 * Note that vector y2 contains the second derivatives of the
				 * interpolating function at the calculation epochs x. */

				/* First make a pointer for u and cudaMalloc device memory for it */
				float *u;
				gpuErrchk(cudaMalloc((void**)&u, sizeof(float) * n));

				lghtcrv_spline_streams_f_krnl<<<BLKncalc,THD>>>(ddat, s, 2.0e30, 2.0e30, u);
				checkErrorAfterKernelLaunch("lghtcrv_spline_streams_f_krnl");

				/* Change launch parameters from ncalc threads to n threads */
				BLKncalc.x = floor((THD.x - 1 + n) / THD.x);
				lghtcrv_splint_streams_f_krnl<<<BLKncalc,THD>>>(ddat, s);
				checkErrorAfterKernelLaunch("lghtcrv_splint_streams_krnl");
				/* Cleanup */
				cudaFree(u);
				cudaFree(so);
				cudaFree(pixels_per_km);
				free(lghtcrv_y);
			}
			break;
		default:
			bailout("vary_params.c: can't handle this dataset type yet\n");
		}

		cudaFree(pos);
		cudaFree(posn);
		cudaFree(ndel);
		cudaFree(ndop);
		cudaFree(outbndarr);
		cudaFree(compute_xsec);
	}

//	dbg_print_pos_z(ddat, 0, 0, 75, "Streams_pos_z_s0f0.csv");
//	dbg_print_pos_cose_s(ddat, 0, 0, 75, "Streams_pos_cose_s_s0f0.csv");

	/* Calculate the zmax, radar cross-section, optical brightness, and cosine
	 * subradar latitude */
	double dd_zmax, rd_xsec, opt_brtns, cs_sb_rdr_lat;
	vps_set_four_parameters_krnl<<<1,1>>>(ddat);
	checkErrorAfterKernelLaunch("vps_set_four_parameters, line ");
	gpuErrchk(cudaMemcpyFromSymbol(&dd_zmax, st_deldop_zmax,
			sizeof(double), 0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&rd_xsec, st_rad_xsec,
			sizeof(double), 0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&opt_brtns, st_opt_brightness,
			sizeof(double), 0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&cs_sb_rdr_lat, st_cos_subradarlat,
			sizeof(double), 0, cudaMemcpyDeviceToHost));

	*deldop_zmax 		= dd_zmax;
	*rad_xsec 			= rd_xsec;
	*opt_brightness		= opt_brtns;
	*cos_subradarlat	= cs_sb_rdr_lat;

	cudaFree(compute_brightness);
	cudaFree(compute_zmax);
	cudaFree(compute_cosdelta);
}

