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

__device__ unsigned char dtype;
__device__ int dcompute_xsec, dcompute_zmax, dcompute_brightness,
			   dcompute_cosdelta, vary_params_dnframes, dpos_n,
			   dndop, dndel, vp_xlim0, vp_xlim1, vp_ylim0, vp_ylim1, vp_n,
			   dncalc, dlghtcrv_bistatic, dlghtcrv_n;
__device__ double dweight, vp_oa[3][3], vp_to_earth[3], vp_deldop_zmax,
			vp_rad_xsec, vp_opt_brightness, vp_cos_subradarlat;
__device__ float vpzmax, sum_deldop_zmax, sum_rad_xsec,	sum_opt_brightness,
				sum_cos_subradarlat, deldop_cross_section, doppler_cross_section;
__device__ struct pos_t *vp_pos;

//__host__ int NearestPowerOf2(int n);

__global__ void vp_init_vars() {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		sum_deldop_zmax = 0.0;
		sum_rad_xsec = 0.0;
		sum_opt_brightness = 0.0;
		sum_cos_subradarlat = 0.0;
	}
}
__global__ void get_data_type_krnl(struct dat_t *ddat, int s) {
	/* nset-threaded kernel */
	if (threadIdx.x == 0) {
		dtype = ddat->set[s].type;
	}
}
__global__ void get_compute_flags_krnl(struct par_t *dpar, struct dat_t *ddat,
		int s, int f) {
	/* Single-threaded kernel */
	if (threadIdx.x ==0) {
		switch (ddat->set[s].type) {
		case DELAY:
			dcompute_zmax = (dpar->vary_delcor0 != VARY_NONE
					&& ddat->set[s].desc.deldop.delcor.a[0].state != 'c');
			dcompute_xsec = (dpar->vary_radalb != VARY_NONE
					&& ddat->set[s].desc.deldop.frame[f].cal.state == 'c');
			dcompute_cosdelta = (dpar->vary_dopscale != VARY_NONE
					&& ddat->set[s].desc.deldop.dopscale.state != 'c');
			vp_pos = &ddat->set[s].desc.deldop.frame[f].pos;
			dweight = ddat->set[s].desc.deldop.frame[f].weight;
			dndel = ddat->set[s].desc.deldop.frame[f].ndel;
			dndop = ddat->set[s].desc.deldop.frame[f].ndop;
			break;
		case DOPPLER:
			dcompute_xsec = (dpar->vary_radalb != VARY_NONE &&
					ddat->set[s].desc.doppler.frame[f].cal.state == 'c');
			dcompute_cosdelta = (dpar->vary_dopscale != VARY_NONE &&
					ddat->set[s].desc.doppler.dopscale.state != 'c');
			vp_pos = &ddat->set[s].desc.doppler.frame[f].pos;
			dweight = ddat->set[s].desc.doppler.frame[f].weight;
			dndop = ddat->set[s].desc.doppler.frame[f].ndop;
			break;
//		case LGHTCRV:
//			vp_pos = &ddat->set[s].desc.lghtcrv.rend[f].pos;
//			dlghtcrv_bistatic = vp_pos->bistatic;
//			dlghtcrv_n = ddat->set[s].desc.lghtcrv.n;
//			break;
		}
	}
}
__global__ void lghtcrv_set_pos_krnl(struct dat_t *ddat, int s, int f) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		vp_pos = &ddat->set[s].desc.lghtcrv.rend[f].pos; /* Backup - delete this later */
		dlghtcrv_bistatic = vp_pos->bistatic;
		dlghtcrv_n = ddat->set[s].desc.lghtcrv.n;

	}
}
__global__ void get_lghtcrv_cb_krnl(struct par_t *dpar, struct dat_t *ddat,
		int s) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		dcompute_brightness = (dpar->vary_optalb != VARY_NONE
				&& ddat->set[s].desc.lghtcrv.cal.state == 'c');
	}
}
__global__ void get_vary_params_nframes_krnl(struct dat_t *ddat, int s) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		switch (ddat->set[s].type) {
		case DELAY:
			vary_params_dnframes = ddat->set[s].desc.deldop.nframes;
			break;
		case DOPPLER:
			vary_params_dnframes = ddat->set[s].desc.doppler.nframes;
			break;
		case LGHTCRV:
			vp_n = ddat->set[s].desc.lghtcrv.n;
			dncalc = ddat->set[s].desc.lghtcrv.ncalc;
			dweight = ddat->set[s].desc.lghtcrv.weight;
		}
	}
}
__global__ void set_ae_oe_bistatic_krnl(struct dat_t *ddat, int s,
		int f) {
	/* 9-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % 3;
	int j = offset / 3;

	if (offset < 9) {

		switch(ddat->set[s].type) {
		case DELAY:
			vp_pos->ae[i][j] =	ddat->set[s].desc.deldop.frame[f].view[ddat->set[s].desc.deldop.v0].ae[i][j];
			vp_pos->oe[i][j] =	ddat->set[s].desc.deldop.frame[f].view[ddat->set[s].desc.deldop.v0].oe[i][j];
			break;
		case DOPPLER:
			vp_pos->ae[i][j] =	ddat->set[s].desc.doppler.frame[f].view[ddat->set[s].desc.doppler.v0].ae[i][j];
			vp_pos->oe[i][j] =	ddat->set[s].desc.doppler.frame[f].view[ddat->set[s].desc.doppler.v0].oe[i][j];
			break;
		case LGHTCRV:
			vp_pos->ae[i][j] = ddat->set[s].desc.lghtcrv.rend[f].ae[i][j];
			vp_pos->oe[i][j] = ddat->set[s].desc.lghtcrv.rend[f].oe[i][j];
			vp_pos->se[i][j] = ddat->set[s].desc.lghtcrv.rend[f].se[i][j];
		}
		/* The following is a single-thread task */
		if (threadIdx.x == 0) {
			if (ddat->set[s].type == LGHTCRV)
				vp_pos->bistatic = 1;
			else
				vp_pos->bistatic = 0;
		}
	}
}
__global__ void get_pos_n_krnl()
{
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		dpos_n = vp_pos->n;
	}
}
__global__ void vp_get_xylims_krnl(struct mod_t *dmod, struct dat_t
		*ddat, int set, int frm) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		switch(ddat->set[set].type) {
		case DELAY:
			vp_xlim0 = ddat->set[set].desc.deldop.frame[frm].pos.xlim[0];
			vp_xlim1 = ddat->set[set].desc.deldop.frame[frm].pos.xlim[1];
			vp_ylim0 = ddat->set[set].desc.deldop.frame[frm].pos.ylim[0];
			vp_ylim1 = ddat->set[set].desc.deldop.frame[frm].pos.ylim[1];
			break;
		case DOPPLER:
			vp_xlim0 = ddat->set[set].desc.doppler.frame[frm].pos.xlim[0];
			vp_xlim1 = ddat->set[set].desc.doppler.frame[frm].pos.xlim[1];
			vp_ylim0 = ddat->set[set].desc.doppler.frame[frm].pos.ylim[0];
			vp_ylim1 = ddat->set[set].desc.doppler.frame[frm].pos.ylim[1];
			break;
		case LGHTCRV:
			vp_xlim0 = ddat->set[set].desc.lghtcrv.rend[frm].pos.xlim[0];
			vp_xlim1 = ddat->set[set].desc.lghtcrv.rend[frm].pos.xlim[1];
			vp_ylim0 = ddat->set[set].desc.lghtcrv.rend[frm].pos.ylim[0];
			vp_ylim1 = ddat->set[set].desc.lghtcrv.rend[frm].pos.ylim[1];
			break;
		}
		vpzmax = -HUGENUMBER;
	}
}
__global__ void deldop_clrvect_krnl(struct dat_t *ddat, int s, int f, int size) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	if (offset < size) {
		ddat->set[s].desc.deldop.frame[f].fit_s[offset] = 0.0;
	}
}
__global__ void zmax_finalize_krnl(float value) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		sum_deldop_zmax += value*dweight;
	}
}
__global__ void compute_xsec_final_krnl(struct dat_t *ddat, float frm_xsec,
		int s, int f) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		switch (ddat->set[s].type) {
		case DELAY:
			deldop_cross_section = __double2float_rd(ddat->set[s].desc.deldop.frame[f].overflow_xsec);
			deldop_cross_section += frm_xsec; // fit is the end result of parallel reduction
			deldop_cross_section *= ddat->set[s].desc.deldop.frame[f].cal.val;
			sum_rad_xsec += deldop_cross_section*ddat->set[s].desc.deldop.frame[f].weight;
			break;
		case DOPPLER:
			doppler_cross_section = __double2float_rd(ddat->set[s].desc.doppler.frame[f].overflow_xsec);
			doppler_cross_section += frm_xsec;
			doppler_cross_section *= ddat->set[s].desc.doppler.frame[f].cal.val;
			sum_rad_xsec += doppler_cross_section*ddat->set[s].desc.doppler.frame[f].weight;

			break;
		}
	}
}
__global__ void compute_cosdelta_krnl(struct dat_t *ddat, int s, int f) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		int j, view;
		double cos_delta;

		switch(ddat->set[s].type){
		case DELAY:
			view = ddat->set[s].desc.deldop.v0;
			/* oa = matrix to transform body-fixed to observer coordinates  */
			/* to_earth = normalized target-to-Earth vector in body-fixed coords  */
			dev_mtrnsps( vp_oa, ddat->set[s].desc.deldop.frame[f].view[view].ae);
			dev_mmmul( vp_oa, ddat->set[s].desc.deldop.frame[f].view[view].oe, vp_oa);
			for (j=0; j<=2; j++)
				vp_to_earth[j] = vp_oa[2][j];
			cos_delta = sqrt(vp_to_earth[0]*vp_to_earth[0] + vp_to_earth[1]*vp_to_earth[1]);
			dweight = ddat->set[s].desc.deldop.frame[f].weight;
			sum_cos_subradarlat += cos_delta*dweight;
			break;
		case DOPPLER:
			view = ddat->set[s].desc.doppler.v0;
			/* oa = matrix to transform body-fixed to observer coordinates  */
			/* to_earth = normalized target-to-Earth vector in body-fixed coords  */
			dev_mtrnsps(vp_oa, ddat->set[s].desc.doppler.frame[f].view[view].ae);
			dev_mmmul(vp_oa, ddat->set[s].desc.doppler.frame[f].view[view].oe, vp_oa);
			for (j=0; j<=2; j++)
				vp_to_earth[j] = vp_oa[2][j];
			cos_delta = sqrt(vp_to_earth[0]*vp_to_earth[0] + vp_to_earth[1]*vp_to_earth[1]);
			dweight = ddat->set[s].desc.doppler.frame[f].weight;
			sum_cos_subradarlat += cos_delta*dweight;
		}
	}
}
__global__ void posclr_krnl(int n, int nx)
{
	/* Multi-threaded kernel (2*pos->n + 1)^2 threads) */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = (offset % nx) - n;
	int j = (offset / nx) - n;

	if (offset < (nx*nx)) {
		/* For each POS pixel, zero out the optical brightness (b) and
		 * cos(scattering angle), reset the z coordinate (distance from COM towards
		 * Earth) to a dummy value, and reset the body, component, and facet onto
		 * which the pixel center projects to  dummy values                  */
//		vp_pos->b[i][j] = vp_pos->cose[i][j] = 0.0;
//		vp_pos->z[i][j] = -HUGENUMBER;
		vp_pos->body[i][j] = vp_pos->comp[i][j] = vp_pos->f[i][j] = -1;

		vp_pos->b_s[offset] = vp_pos->cose_s[offset] = 0.0;
		vp_pos->z_s[offset] = -HUGENUMBER;

		/* In the x direction, reset the model's leftmost and rightmost
		 * pixel number to dummy values, and similarly for the y direction   */
		vp_pos->xlim[0] = vp_pos->ylim[0] =  n;
		vp_pos->xlim[1] = vp_pos->ylim[1] = -n;

		/* For a bistatic situation (lightcurve or plane-of-sky dataset), zero out
		 * cos(incidence angle) and reset the distance towards the sun, the body,
		 * component, and facet numbers as viewed from the sun, and the model's
		 * maximum projected extent as viewed from the sun to dummy values    */
		if (vp_pos->bistatic) {
			vp_pos->cosill[i][j] = 0.0;
			vp_pos->zill[i][j] = -HUGENUMBER;
			vp_pos->bodyill[i][j] = vp_pos->compill[i][j] = vp_pos->fill[i][j] = -1;

			vp_pos->cosill_s[offset] = 0.0;
			vp_pos->zill_s[offset] = -HUGENUMBER;

			vp_pos->xlim2[0] = vp_pos->ylim2[0] =  n;
			vp_pos->xlim2[1] = vp_pos->ylim2[1] = -n;
		}
	}
	__syncthreads();
}
__global__ void doppler_clrvect_krnl(struct dat_t *ddat, int s, int f) {
	/* ndop-threaded kernel */
	int i = 1 + blockIdx.x * blockDim.x + threadIdx.x;
	if (i <= dndop){
		ddat->set[s].desc.doppler.frame[f].fit_s[i] = 0.0;
	}
}
__global__ void doppler_compute_xsec_krnl(struct dat_t *ddat, int s,
		int f) {
	/* ndop-threaded kernel */
	/* Deprecated and for debug use only */
	int dop = blockIdx.x * blockDim.x + threadIdx.x;
	//float fit;	// the fit value for this thread/doppler bin only

	if (dop == 0)
		doppler_cross_section = __double2float_rd(ddat->set[s].desc.doppler.frame[f].overflow_xsec);
	__syncthreads();

	if (dop < dndop) {
		atomicAdd(&doppler_cross_section, ddat->set[s].desc.doppler.frame[f].fit_s[dop+1]);
	}
	__syncthreads();

	if (dop == 1) {
		doppler_cross_section *= ddat->set[s].desc.doppler.frame[f].cal.val;
		sum_rad_xsec += doppler_cross_section*dweight;
	}
 	__syncthreads();
}
__global__ void posmask_krnl(struct par_t *dpar, int nThreads, int xspan)
{
	/* multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int n = dpos_n;
	int i = offset % xspan - n;
	int j = offset / xspan - n;
	double tol = dpar->mask_tol;
	int im, jm, i1, j1, i2, j2, i_sign, j_sign;
	double xk[3], so[3][3], pixels_per_km, i0_dbl, j0_dbl, zill, t, u, bignum;

	if (offset == 0){
		bignum = 0.99*HUGENUMBER;  /* z = -HUGENUMBER for blank-sky pixels */

		dev_mtrnsps( so, vp_pos->oe);
		dev_mmmul( so, vp_pos->se, so);    /* so takes obs into src coords */
		pixels_per_km = 1/vp_pos->km_per_pixel;
	}
	__syncthreads();

	/*  Loop through all POS pixels  */
	if (offset < nThreads) {
		//	n = vp_pos->n;
		//	for (i=(-n); i<=n; i++) {               /* for each pixel in the */
		//		for (j=(-n); j<=n; j++) {             /* observer's view */
		if (vp_pos->cose_s[offset] != 0.0) {     /* if there's something there */
			xk[0] = i*vp_pos->km_per_pixel;     /* calculate 3D position */
			xk[1] = j*vp_pos->km_per_pixel;
			xk[2] = vp_pos->z_s[offset];

			/* Given the observer coordinates x of of POS pixel (i,j), find
			 * which pixel (im,jm) this corresponds to in the projected view as
			 * seen from the source (sun)             */

			dev_cotrans2( xk, so, xk, 1);           /* go into source coordinates */
			i0_dbl = xk[0]*pixels_per_km;     /* unrounded (double precision) */
			j0_dbl = xk[1]*pixels_per_km;
			im = dev_vp_iround( i0_dbl);            /* center of nearest pixel in mask */
			jm = dev_vp_iround( j0_dbl);

			/* If center of projected pixel "seen" from source (as determined
			 * by routine posvis) lies within the boundaries of the mask,
			 * projects onto model rather than onto blank space, and represents
			 * a body, component, and facet different from those seen in the
			 * POS, calculate distance from mask pixel to source and compare to
			 * distance from POS pixel to source.                             */

			if (fabs(i0_dbl) < n && fabs(j0_dbl) < n
					&& vp_pos->zill[im][jm] > -bignum
					&& (vp_pos->f[i][j]    != vp_pos->fill[im][jm]    ||
							vp_pos->comp[i][j] != vp_pos->compill[im][jm] ||
							vp_pos->body[i][j] != vp_pos->bodyill[im][jm]    )) {

				/* Rather than using distance towards source of CENTER of mask
				 * pixel, use bilinear interpolation to get distance towards
				 * source where the line between source and POS pixel's center
				 * intersects the mask pixel.                                */
				i1 = (int) floor( i0_dbl);
				j1 = (int) floor( j0_dbl);

				if (vp_pos->zill[i1][j1]     > -bignum &&
						vp_pos->zill[i1+1][j1]   > -bignum &&
						vp_pos->zill[i1][j1+1]   > -bignum &&
						vp_pos->zill[i1+1][j1+1] > -bignum    ) {

					/* Do standard bilinear interpolation: None of the four
					 * surrounding "grid square" pixels in the mask is
					 * blank sky                           */
					t = i0_dbl - i1;
					u = j0_dbl - j1;
					zill = (1 - t)*(1 - u)*vp_pos->zill[i1][j1]
                               + t*(1 - u)*vp_pos->zill[i1+1][j1]
                                     + t*u*vp_pos->zill[i1+1][j1+1]
	                           + (1 - t)*u*vp_pos->zill[i1][j1+1];
				} else {

					/* The following code block is a kludge: One or more of the
					 * four surrounding "grid square" pixels in mask is blank
					 * sky, so standard bilinear interpolation won't work.  */
					zill = vp_pos->zill[im][jm];

					i_sign = (i0_dbl >= im) ? 1 : -1;
					i2 = im + i_sign;
					if (abs(i2) <= n && vp_pos->zill[i2][jm] > -bignum) {
						zill += fabs(i0_dbl - im)
           				  * (vp_pos->zill[i2][jm] - vp_pos->zill[im][jm]);
					} else {
						i2 = im - i_sign;
						if (abs(i2) <= n && vp_pos->zill[i2][jm] > -bignum)
							zill -= fabs(i0_dbl - im)
							* (vp_pos->zill[i2][jm] - vp_pos->zill[im][jm]);
					}

					j_sign = (j0_dbl >= jm) ? 1 : -1;
					j2 = jm + j_sign;
					if (abs(j2) <= n && vp_pos->zill[im][j2] > -bignum) {
						zill += fabs(j0_dbl - jm)
                          * (vp_pos->zill[im][j2] - vp_pos->zill[im][jm]);
					} else {
						j2 = jm - j_sign;
						if (abs(j2) <= n && vp_pos->zill[im][j2] > -bignum)
							zill -= fabs(j0_dbl - jm)
							* (vp_pos->zill[im][j2] - vp_pos->zill[im][jm]);
					}
				}

				/* If interpolated point within mask pixel is at least tol km
				 * closer to source than is the center of POS pixel, the facet
				 * represented by the mask pixel is shadowing the POS pixel:
				 * represent this by setting
				 * 		cos(scattering angle) = 0.0 for the POS pixel.      */
				if (zill - xk[2] > tol)
					vp_pos->cose_s[offset] = 0.0;
			}
		}
	}
}
__global__ void posmask_universal_krnl(struct par_t *dpar, struct pos_t *pos,
		int nThreads, int xspan)
{
	/* multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int n = pos->n;
	int i = offset % xspan - n;
	int j = offset / xspan - n;
	double tol = dpar->mask_tol;
	int im, jm, i1, j1, i2, j2, i_sign, j_sign;
	double xk[3], so[3][3], pixels_per_km, i0_dbl, j0_dbl, zill, t, u, bignum;

	if (offset == 0){
		bignum = 0.99*HUGENUMBER;  /* z = -HUGENUMBER for blank-sky pixels */

		dev_mtrnsps( so, pos->oe);
		dev_mmmul( so, pos->se, so);    /* so takes obs into src coords */
		pixels_per_km = 1/pos->km_per_pixel;
	}
	__syncthreads();

	/*  Loop through all POS pixels  */
	if (offset < nThreads) {
		if (pos->cose_s[offset] != 0.0) {     /* if there's something there */
			xk[0] = i*pos->km_per_pixel;     /* calculate 3D position */
			xk[1] = j*pos->km_per_pixel;
			xk[2] = pos->z_s[offset];

			/* Given the observer coordinates x of of POS pixel (i,j), find
			 * which pixel (im,jm) this corresponds to in the projected view as
			 * seen from the source (sun)             */

			dev_cotrans2( xk, so, xk, 1);           /* go into source coordinates */
			i0_dbl = xk[0]*pixels_per_km;     /* unrounded (double precision) */
			j0_dbl = xk[1]*pixels_per_km;
			im = dev_vp_iround( i0_dbl);            /* center of nearest pixel in mask */
			jm = dev_vp_iround( j0_dbl);

			/* If center of projected pixel "seen" from source (as determined
			 * by routine posvis) lies within the boundaries of the mask,
			 * projects onto model rather than onto blank space, and represents
			 * a body, component, and facet different from those seen in the
			 * POS, calculate distance from mask pixel to source and compare to
			 * distance from POS pixel to source.                             */

			if (fabs(i0_dbl) < n && fabs(j0_dbl) < n
					&& pos->zill[im][jm] > -bignum
					&&(pos->f[i][j]    != pos->fill[im][jm]    ||
					   pos->comp[i][j] != pos->compill[im][jm] ||
					   pos->body[i][j] != pos->bodyill[im][jm]    )) {

				/* Rather than using distance towards source of CENTER of mask
				 * pixel, use bilinear interpolation to get distance towards
				 * source where the line between source and POS pixel's center
				 * intersects the mask pixel.                                */
				i1 = (int) floor( i0_dbl);
				j1 = (int) floor( j0_dbl);

				if (pos->zill[i1][j1] > -bignum && pos->zill[i1+1][j1] > -bignum &&
					pos->zill[i1][j1+1] > -bignum && pos->zill[i1+1][j1+1] > -bignum) {

					/* Do standard bilinear interpolation: None of the four
					 * surrounding "grid square" pixels in the mask is
					 * blank sky                           */
					t = i0_dbl - i1;
					u = j0_dbl - j1;
					zill = (1 - t)*(1 - u)*pos->zill[i1][j1]
                               + t*(1 - u)*pos->zill[i1+1][j1]
                                     + t*u*pos->zill[i1+1][j1+1]
	                           + (1 - t)*u*pos->zill[i1][j1+1];
				} else {

					/* The following code block is a kludge: One or more of the
					 * four surrounding "grid square" pixels in mask is blank
					 * sky, so standard bilinear interpolation won't work.  */
					zill = pos->zill[im][jm];

					i_sign = (i0_dbl >= im) ? 1 : -1;
					i2 = im + i_sign;
					if (abs(i2) <= n && pos->zill[i2][jm] > -bignum) {
						zill += fabs(i0_dbl - im)  *
								(pos->zill[i2][jm] - pos->zill[im][jm]);
					} else {
						i2 = im - i_sign;
						if (abs(i2) <= n && pos->zill[i2][jm] > -bignum)
							zill -= fabs(i0_dbl - im)
							* (pos->zill[i2][jm] - pos->zill[im][jm]);
					}

					j_sign = (j0_dbl >= jm) ? 1 : -1;
					j2 = jm + j_sign;
					if (abs(j2) <= n && pos->zill[im][j2] > -bignum) {
						zill += fabs(j0_dbl - jm)
                          * (pos->zill[im][j2] - pos->zill[im][jm]);
					} else {
						j2 = jm - j_sign;
						if (abs(j2) <= n && pos->zill[im][j2] > -bignum)
							zill -= fabs(j0_dbl - jm)
							* (pos->zill[im][j2] - pos->zill[im][jm]);
					}
				}

				/* If interpolated point within mask pixel is at least tol km
				 * closer to source than is the center of POS pixel, the facet
				 * represented by the mask pixel is shadowing the POS pixel:
				 * represent this by setting
				 * 		cos(scattering angle) = 0.0 for the POS pixel.      */
				if (zill - xk[2] > tol)
					pos->cose_s[offset] = 0.0;
			}
		}
	}
}
__global__ void lghtcrv_copy_y_krnl(struct dat_t *ddat, double host_value,
		int set, int i) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		ddat->set[set].desc.lghtcrv.y[i] = host_value;
	}
}
__global__ void lghtcrv_spline_krnl(struct dat_t *ddat, int set, double
		yp1, double ypn, double *u) {
	/*(double *x  - lghtcrv->x
	 * double *y  - lghtcrv->y
	 * int n      - calc
	 * double yp1 - 2.0e30
	 * double ypn - 2.0e30
	 * double *y2 - lghtcrv->y2)*/

	/* Multi-threaded kernel */
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int n = dncalc;
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
__global__ void lghtcrv_splint_krnl(struct dat_t *ddat, int set)
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
	int n = dncalc;

	int klo,khi,k;
	double h,b,a;

	if (i < dlghtcrv_n) {

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
__global__ void vp_set_four_parameters_krnl(struct dat_t *ddat) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		if (ddat->sum_deldop_zmax_weights > 0.0)
			vp_deldop_zmax = sum_deldop_zmax / ddat->sum_deldop_zmax_weights;
		else
			vp_deldop_zmax = 0.0;
		if (ddat->sum_rad_xsec_weights > 0.0) {
			vp_rad_xsec = sum_rad_xsec / ddat->sum_rad_xsec_weights;			}
		else
			vp_rad_xsec = 0.0;
		if (ddat->sum_opt_brightness_weights > 0.0)
			vp_opt_brightness = sum_opt_brightness / ddat->sum_opt_brightness_weights;
		else
			vp_opt_brightness = 0.0;
		if (ddat->sum_cos_subradarlat_weights > 0.0)
			vp_cos_subradarlat = sum_cos_subradarlat / ddat->sum_cos_subradarlat_weights;
		else
			vp_cos_subradarlat = 0.0;
	}
}

__host__ void vary_params_cuda( struct par_t *dpar, struct mod_t *dmod,
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
	int c=0, f, s, i, compute_xsec, compute_brightness, compute_zmax,
			compute_cosdelta, n, ncalc, pos_n, nx, lghtcrv_bistatic, nframes,
			xlim[2], ylim[2], xspan, lghtcrv_n;
	double weight;
	dim3 BLK,THD;
	unsigned char type;
	int ndel, ndop, nThreads;
	struct pos_t *pos;	/* Experimental for lghtcrv use and reuse */

	/*  Initialize variables  */
	vp_init_vars<<<1,1>>>();
	checkErrorAfterKernelLaunch("vp_init_krnl, line ");

	/* Process each dataset in turn */
	for (s=0; s<nsets; s++) {
		/* Get set's data type from GPU to host so we can SWITCH it */
		get_data_type_krnl<<<1,1>>>(ddat, s);
		checkErrorAfterKernelLaunch ("get_data_type_krnl, line ");
		gpuErrchk(cudaMemcpyFromSymbol(&type, dtype, sizeof(dtype), 0,
				cudaMemcpyDeviceToHost));

		switch (type) {
		case DELAY:
			/* Get nframes */
			get_vary_params_nframes_krnl<<<1,1>>>(ddat, s);
			checkErrorAfterKernelLaunch ("get_data_type_krnl, line ");
			gpuErrchk(cudaMemcpyFromSymbol(&nframes, vary_params_dnframes,
				sizeof(vary_params_dnframes), 0, cudaMemcpyDeviceToHost));

			for (f=0; f<nframes; f++) {
				/* Get the compute_zmax and compute x_sec flags in a single-
				 * threaded kernel */
				get_compute_flags_krnl<<<1,1>>>(dpar, ddat, s, f);
				checkErrorAfterKernelLaunch("get_compute_flags_krnl, line ");
				gpuErrchk(cudaMemcpyFromSymbol(&compute_zmax, dcompute_zmax,
						sizeof(dcompute_zmax), 0, cudaMemcpyDeviceToHost));
				gpuErrchk(cudaMemcpyFromSymbol(&compute_xsec, dcompute_xsec,
						sizeof(dcompute_xsec), 0, cudaMemcpyDeviceToHost));
				gpuErrchk(cudaMemcpyFromSymbol(&compute_cosdelta, dcompute_cosdelta,
						sizeof(dcompute_cosdelta), 0, cudaMemcpyDeviceToHost));
				gpuErrchk(cudaMemcpyFromSymbol(&weight, dweight,
						sizeof(dweight), 0, cudaMemcpyDeviceToHost));
				gpuErrchk(cudaMemcpyFromSymbol(&ndel, dndel,
						sizeof(dndel), 0, cudaMemcpyDeviceToHost));
				gpuErrchk(cudaMemcpyFromSymbol(&ndop, dndop,
						sizeof(dndop), 0, cudaMemcpyDeviceToHost));

				/* Start the if block for computing zmax and/or cross-section */
				if (compute_zmax || compute_xsec) {
					/* Launch 9-threaded kernel to set up ae[3][3] and oe[3][3]
					 * and also set the bistatic flag 		 */
					THD.x = 9;
					set_ae_oe_bistatic_krnl<<<BLK,THD>>>(ddat, s, f);
					checkErrorAfterKernelLaunch("deldop_set_ae_oe_bistatic_krnl, line ");

					/* Need to get pos->n for kernel launch first */
					get_pos_n_krnl<<<1,1>>>();
					checkErrorAfterKernelLaunch("get_deldop_pos_n_krnl, line ");
					gpuErrchk(cudaMemcpyFromSymbol(&pos_n, dpos_n, sizeof(dpos_n),
							0, cudaMemcpyDeviceToHost));

					/* Configure & launch posclr_krnl to initialize POS view */
					nThreads = (2*pos_n+1)*(2*pos_n+1);
					BLK.x = floor((maxThreadsPerBlock - 1 + nThreads) /
							maxThreadsPerBlock);
					THD.x = maxThreadsPerBlock; // Thread block dimensions
					nx = 2*pos_n + 1;
					posclr_krnl<<<BLK,THD>>>(pos_n, nx);
					checkErrorAfterKernelLaunch("posclr_krnl, line ");

					/*  Determine which POS pixels cover the target, and
	                    get the distance toward Earth of each POS pixel   */
					posvis_cuda_2(dpar, dmod, ddat, orbit_offset,s,f,
							0, 0, c);

					/* Zero out the fit delay-Doppler image and call pos2deldop
					 * to create the fit image by mapping power from the plane
					 * of sky to delay-Doppler space.    				      *
					 * (Call to clrmat has been removed. It simply zeroed out *
					 * frame->fit. I am using a separate array (fit) which is *
					 * automatically zeroed at creation time via cudaCalloc   */
					/* Configure and launch the clrmat kernel for deldop data */
					/* This will be a (ndel*ndop)-threaded kernel */

					/* Call the CUDA pos2deldop function */
					/* par, mod->photo,orbit_xoff,orbit_yoff,orbit_dopoff=0.0
					 * deldop, body = 0,s,f,v=0
					 * fit is the single pointer replacing the frame->fit array
					 * which did not copy well (?).  We use this separate fit
					 * array instead		 */

					/* Clear out the fit_s array first */
					BLK.x = floor((maxThreadsPerBlock - 1 + (ndel*ndop)) /
							maxThreadsPerBlock);
					THD.x = maxThreadsPerBlock; // Thread block dimensions
					deldop_clrvect_krnl<<<BLK,THD>>>(ddat, s, f, (ndel*ndop));
					checkErrorAfterKernelLaunch("deldop_clrvect_krnl in vary_params_cuda");

					pos2deldop_cuda_2(dpar,dmod,ddat,0.0,0.0,0.0,0,s,f,0);

					/* Compute distance toward Earth of the subradar point  */
					if (compute_zmax) {
						float zmax = 0.0;
						int size = 0;
						size = (2*pos_n+1)*(2*pos_n+1);

						/* Configure and launch the parallel reduction kernel to
						 * find distance toward Earth of the subradar point */
						zmax = compute_pos_zmax(ddat, size, s, f);
						zmax_finalize_krnl<<<1,1>>>(zmax);
						checkErrorAfterKernelLaunch("zmax_finalize_krnl in vary_params_cuda");
					}

					/*  Compute cross section  */
					if (compute_xsec) {
						float xsec=0.0;
						xsec = compute_deldop_xsec_snglkrnl(ddat, ndel, ndop, s, f);
						compute_xsec_final_krnl<<<1,1>>>(ddat, xsec, s, f);
						checkErrorAfterKernelLaunch("compute_xsec_final_krnl (deldop)");
					}
				}
				if (compute_cosdelta) {
					/* Launch single-thread kernel to compute sum_cos_subradarlat */
					compute_cosdelta_krnl<<<1,1>>>(ddat, s, f);
					checkErrorAfterKernelLaunch("deldop_compute_cosdelta_krnl, line ");
				}
			}
			break;
		case DOPPLER:
			/* Get nframes */
			get_vary_params_nframes_krnl<<<1,1>>>(ddat, s);
			checkErrorAfterKernelLaunch ("get_data_type_krnl, line ");
			gpuErrchk(cudaMemcpyFromSymbol(&nframes, vary_params_dnframes,
					sizeof(vary_params_dnframes), 0, cudaMemcpyDeviceToHost));

			for (f=0; f<nframes; f++) {
				/* Launch single-thread kernel to get our compute flags first */
				get_compute_flags_krnl<<<1,1>>>(dpar, ddat, s, f);
				checkErrorAfterKernelLaunch("doppler_get_compute_flags, line ");
				gpuErrchk(cudaMemcpyFromSymbol(&compute_xsec, dcompute_xsec,
						sizeof(dcompute_xsec), 0, cudaMemcpyDeviceToHost));
				gpuErrchk(cudaMemcpyFromSymbol(&compute_cosdelta, dcompute_cosdelta,
						sizeof(dcompute_cosdelta), 0, cudaMemcpyDeviceToHost));
				gpuErrchk(cudaMemcpyFromSymbol(&weight, dweight,
						sizeof(dweight), 0, cudaMemcpyDeviceToHost));
				gpuErrchk(cudaMemcpyFromSymbol(&ndop, dndop,
						sizeof(dndop), 0, cudaMemcpyDeviceToHost));

				if (compute_xsec) {
					/* Launch 9-threaded kernel to set up ae[3][3] and oe[3][3]
					 * and also set the bistatic flag 		 */
					THD.x = 9;
					set_ae_oe_bistatic_krnl<<<BLK,THD>>>(ddat, s, f);
					checkErrorAfterKernelLaunch("deldop_set_ae_oe_bistatic_krnl, line ");

					/* Need to get pos->n for kernel launch first */
					get_pos_n_krnl<<<1,1>>>();
					checkErrorAfterKernelLaunch("get_deldop_pos_n_krnl, line ");
					gpuErrchk(cudaMemcpyFromSymbol(&pos_n, dpos_n, sizeof(dpos_n),
							0, cudaMemcpyDeviceToHost));

					/* Configure & launch posclr_krnl to initialize POS view */
					BLK.x = floor((maxThreadsPerBlock - 1 + (2*pos_n+1)*(2*pos_n+1)) /
							maxThreadsPerBlock);
					THD.x = maxThreadsPerBlock; // Thread block dimensions
					nx = 2*pos_n + 1;
					posclr_krnl<<<BLK,THD>>>(pos_n, nx);
					checkErrorAfterKernelLaunch("posclr_krnl, line ");

					/* Determine which POS pixels cover the target  */
					posvis_cuda_2(dpar,dmod,ddat,orbit_offset,s,f,0,0,c);

					/* Zero out the fit Doppler spectrum, then call pos2doppler to create the fit
					 * spectrum by mapping power from the plane of the sky to Doppler space.      */
					BLK.x = floor((maxThreadsPerBlock - 1 + ndop) / maxThreadsPerBlock);
					doppler_clrvect_krnl<<<BLK,THD>>>(ddat, s, f);
					checkErrorAfterKernelLaunch("doppler_clrvect_krnl, line ");

					pos2doppler_cuda_2(dpar,dmod,ddat,0.0,0.0,0.0,0,s,f,0);

					/* Compute cross section */
					float xsec=0.0;
					xsec = compute_doppler_xsec(ddat, ndop, s, f);
					compute_xsec_final_krnl<<<1,1>>>(ddat, xsec, s, f);
					checkErrorAfterKernelLaunch("compute_xsec_final_krnl (Doppler)");
				}
				if (compute_cosdelta) {
					/* Launch single-thread kernel to compute sum_cos_subradarlat */
					compute_cosdelta_krnl<<<1,1>>>(ddat, s, f);
					checkErrorAfterKernelLaunch("doppler_compute_cosdelta_krnl, line ");
				}
			}
			break;
		case POS:
			break;
		case LGHTCRV:
			/* Figure out the compute_brightness flag first */
			get_lghtcrv_cb_krnl<<<1,1>>>(dpar, ddat, s);
			checkErrorAfterKernelLaunch("get_lghtcrv_cb_krnl");
			gpuErrchk(cudaMemcpyFromSymbol(&compute_brightness, dcompute_brightness,
					sizeof(int), 0, cudaMemcpyDeviceToHost));

			if (compute_brightness) {
				/* Launch single-thread kernel to get lghtcrv parameters */
				get_vary_params_nframes_krnl<<<1,1>>>(ddat, s);
				checkErrorAfterKernelLaunch("lghtcrv_get_params, line ");
				gpuErrchk(cudaMemcpyFromSymbol(&n, vp_n, sizeof(vp_n),
							0, cudaMemcpyDeviceToHost));
				gpuErrchk(cudaMemcpyFromSymbol(&ncalc, dncalc, sizeof(ncalc),
							0, cudaMemcpyDeviceToHost));
				gpuErrchk(cudaMemcpyFromSymbol(&weight, dweight, sizeof(weight),
							0, cudaMemcpyDeviceToHost));
				double lghtcrv_y;

				for (i=1; i<=ncalc; i++) {

					/* Launch kernel to get compute flags and set pos */
					lghtcrv_set_pos_krnl<<<1,1>>>(ddat, s, i);
					checkErrorAfterKernelLaunch("lghtcrv_set_pos_krnl");
					gpuErrchk(cudaMemcpyFromSymbol(&lghtcrv_bistatic, dlghtcrv_bistatic,
									sizeof(lghtcrv_bistatic), 0, cudaMemcpyDeviceToHost));
					gpuErrchk(cudaMemcpyFromSymbol(&lghtcrv_n, dlghtcrv_n,
							sizeof(lghtcrv_n), 0, cudaMemcpyDeviceToHost));

					/* Launch 9-threaded kernel to set up ae[3][3] and oe[3][3]
					 * and also set the bistatic flag 		 */
					THD.x = 9;
					set_ae_oe_bistatic_krnl<<<BLK,THD>>>(ddat, s, f);
					checkErrorAfterKernelLaunch("lghtcrv_set_ae_oe_bistatic_krnl, line ");

					/* Need to get pos->n for kernel launch first */
					get_pos_n_krnl<<<1,1>>>();
					checkErrorAfterKernelLaunch("get_lghtcrv_pos_n_krnl, line ");
					gpuErrchk(cudaMemcpyFromSymbol(&pos_n, dpos_n, sizeof(dpos_n),
							0, cudaMemcpyDeviceToHost));

					/* Configure & launch posclr_krnl to initialize POS view */
					xspan = 2*pos_n+1;
					nThreads = xspan*xspan;
					BLK.x = floor((maxThreadsPerBlock-1 + nThreads)/maxThreadsPerBlock);
					THD.x = maxThreadsPerBlock;
					posclr_krnl<<<BLK,THD>>>(pos_n, xspan);
					checkErrorAfterKernelLaunch("posclr_krnl, line ");

					/* Determine which POS pixels cover the target */
//					for (c=0; c<mod->shape.ncomp; c++)
					if (STREAMS)
						printf("in vary params_cuda, fix this (lightcurve posvis call");
					//	posvis_cuda_streams(dpar, dmod, ddat, orbit_offset,s,i, 0, 0, c);
					else
						posvis_cuda_2(dpar, dmod, ddat, orbit_offset,s,i, 0, 0, c);

				 /* Now view the model from the source (sun) and get the facet
				  * number and distance toward the source of each pixel in this
				  * projected view; use this information to determine which POS
				  * pixels are shadowed */
					if (lghtcrv_bistatic) {
//						for (c=0; c<mod->shape.ncomp; c++)
						posvis_cuda_2(dpar,dmod,ddat,orbit_offset,s,i,1,0,c);

						/* Launch parameters still same as before for posclr */
						//posmask_universal_krnl<<<BLK,THD>>>(dpar, nThreads, xspan);
						posmask_krnl<<<BLK,THD>>>(dpar, nThreads, xspan);
						checkErrorAfterKernelLaunch("posmask_krnl (vary_params_cuda.cu)");
					}
					/* Compute model brightness for this lightcurve point */
					/* lghtcrv->y[ncalc]: calculated points for interpolation,
					 * ncalc-points total 					 */
					lghtcrv_y = apply_photo_cuda(dmod, ddat, 0, s, i);

					/* Now launch a kernel to copy it over to the actual lghtcrv */
					lghtcrv_copy_y_krnl<<<1,1>>>(ddat, lghtcrv_y, s, i);
					checkErrorAfterKernelLaunch("vp_copy_lghtcrv_y_krnl, line ");
				}

				/* Now that we have calculated the model lightcurve brightnesses
				 * y at each of the epochs x, we use cubic spline interpolation
				 * (Numerical Recipes routines spline and splint) to get model
				 * lightcurve brightness fit[i] at each OBSERVATION epoch t[i],
				 * with i=1,2,...,n.  This will allow us (in routine chi2) to
				 * compare model to data (fit[i] to obs[i]) to get chi squared.
				 * Note that vector y2 contains the second derivatives of the
				 * interpolating function at the calculation epochs x. */

				/* To-Do:   The splint kernel can be sped up by implementing a
				 * 			proper parallel reduction.				 */

				/* First make a pointer for u and cudaMalloc device memory for it */
				double *u;
				gpuErrchk(cudaMalloc((void**)&u, sizeof(double) * lghtcrv_n));

				BLK.x = floor((maxThreadsPerBlock - 1 + ncalc) /
						maxThreadsPerBlock);
				THD.x = maxThreadsPerBlock; // Thread block dimensions
				lghtcrv_spline_krnl<<<BLK,THD>>>(ddat, s, 2.0e30, 2.0e30, u);
				checkErrorAfterKernelLaunch("lghtcrv_spline_krnl, line ");

				BLK.x = floor((maxThreadsPerBlock - 1 + lghtcrv_n) /
						maxThreadsPerBlock);
				lghtcrv_splint_krnl<<<BLK,THD>>>(ddat, s);
				checkErrorAfterKernelLaunch("lghtcrv_splint_krnl, line ");
				/* Cleanup */
				cudaFree(u);
			}
			break;
		default:
			bailout("vary_params.c: can't handle this dataset type yet\n");
		}
	}

	/* Launch single-threaded kernel to wrap things up in vary_params_cuda,
	 * includes manipulation of the hostside variables deldop_zmax, rad_xsec,
	 * opt_brightness, and cos_subradarlat		 */
	double dd_zmax, rd_xsec, opt_brtns, cs_sb_rdr_lat;
	vp_set_four_parameters_krnl<<<1,1>>>(ddat);
	checkErrorAfterKernelLaunch("vp_set_four_parameters, line ");
	gpuErrchk(cudaMemcpyFromSymbol(&dd_zmax, vp_deldop_zmax,
			sizeof(double), 0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&rd_xsec, vp_rad_xsec,
			sizeof(double), 0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&opt_brtns, vp_opt_brightness,
			sizeof(double), 0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&cs_sb_rdr_lat, vp_cos_subradarlat,
			sizeof(double), 0, cudaMemcpyDeviceToHost));

	*deldop_zmax = dd_zmax;
	*rad_xsec = rd_xsec;
	*opt_brightness = opt_brtns;
	*cos_subradarlat = cs_sb_rdr_lat;
	}

