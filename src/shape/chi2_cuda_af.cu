/*****************************************************************************************
                                                                                   chi2.c

Compute chi-square for a model that was previously created.

Modified 2016 December 8 by ME:
	Converted to "FIT" action-only, CUDA code

Modified 2015 June 10 by CM:
    Implement smearing for the "fit" and "write" actions

Modified 2013 April 24 by CM:
    Adjust names of output files so they are in alphanumeric order if > 100 per dataset

Modified 2012 December 6 by CM:
    Fix bug introduced on December 6: take "calval" into account when determining whether
        or not model signal is greater than 'chi2fit0_thresh' sigmas

Modified 2012 December 5 by CM:
    When the "write_chi2fit0" parameter is turned on, display the number of degrees of
        freedom with zero model signal in addition to the chi2 contribution, and list
        chi2 and dof for all data combined
    Implement the 'chi2fit0_thresh' parameter

Modified 2012 March 24 by CM:
    When the root node calls chi2 with list_breakdown = 1 (e.g., for the "write" action),
        print a warning if the value of any (delay-)Doppler dataset's Doppler scaling
        factor is out of the allowed range

Modified 2010 August 25 by CM:
    Move TINYCALFACT definition to const.h

Modified 2010 August 10 by CM:
    Implement the "radfitmin" and "radobsmin" parameters: these are the
        pixel values that map to black for all fit and obs pgm images
        that are output for delay-Doppler frames

Modified 2010 March 20 by CM:
    For the "write" action for lightcurve datasets, include magnitude
        uncertainties as a new column in output files fit_MM.dat

Modified 2009 August 9 by CM:
    For the "write" action with the "listfit" parameter turned on, replace
        a ".rdf" or ".fits" or ".fit" suffix with ".fitdat"

Modified 2009 April 3 by CM:
    When the root node calls chi2 with list_breakdown = 1 (e.g., for the
        "write" action), print a warning if any plane-of-sky fit image is
        too small to "contain" all nonzero pixels in the POS sky rendering
        or if the model is too wide in delay-Doppler space for any
        (delay-)Doppler fit frame to be correctly constructed
    For MPI_Recv calls, mpi_par[0] is no longer equal to the MPI action,
        since the message tag argument already serves that purpose (as of
        2008 April 10) -- so the other mpi_par elements are renumbered

Modified 2008 June 29 by CM:
    For the "write" and "orbit" actions, zero out fit pixels/bins in
        delay-Doppler, Doppler, and plane-of-sky frames for which those
        pixels/bins have been zeroed out in the pixel mask

Modified 2008 April 10 by CM:
    Use message tag argument to MPI_Recv to identify the MPI action

Modified 2007 September 13 by CM:
    Implement "write_chi2fit0" parameter: for the "write" and "orbit"
        actions, output chi2 for delay-Doppler, Doppler, and plane-of-sky
        pixels/bins with zero model power; do this both for each individual
        frame, for all delay-Doppler data taken together, for all Doppler
        data taken together, and for all plane-of-sky data taken together

Modified 2007 August 18 by CM:
    Rename MPI_TAG to MPI_TAG_1 to avoid name conflict with mpich headers

Modified 2007 August 10 by CM:
    Create chi2_deldop, chi2_doppler, chi2_poset, and chi2_lghtcrv
        routines so that the code isn't one giant switch statement
    When a calfact value is negative, reset it to a tiny positive value
        rather than to zero in order to protect against division by zero
    Implement the "radfitmax" and "radobsmax" parameters for applying the
        same brightness scale to all fit and obs pgm images output for
        delay-Doppler frames

Modified 2007 August 4 by CM:
    Carry out "write" action steps for the "orbit" action as well

Modified 2007 January 6 by CM:
    For the "write" action for lightcurve datasets, output rotation phases
        to files fit_MM.dat and calc_MM.dat.

Modified 2006 October 1 by CM:
    For lightcurve datasets, chi-square is now computed in intensity space
        rather than in magnitude space.  (File output is still in
        magnitudes.)

Modified 2006 June 21 by CM:
    In delay-Doppler section, changed delres to del_per_pixel and dopres to
        dop_per_pixel
    For POS renderings and plane-of-sky fit frames, changed res to
        km_per_pixel
    When the root node calls chi2 with list_breakdown = 1, print a warning
        if the model extends beyond the boundaries of the POS frames, if
        any photometric parameter has an illegal value, or if any ellipsoid
        diameter is tiny or negative.  (This change will cause such
        warnings to be displayed for the "write" action.)

Modified 2006 June 18 by CM:
    Allow each delay-Doppler frame within a dataset to have different
        dimensions after vignetting
    Allow each Doppler frame within a dataset to have different
        dimensions after vignetting
    Allow plane-of-sky frames to be rectangular rather than square
    For the "write" action, adjust output for delay-Doppler, Doppler, and
        plane-of-sky frames to allow for masked-out pixels

Modified 2005 June 27 by CM:
    Renamed "round" function to "iround" to avoid conflict

Modified 2005 April 25 by CM:
    For the "write" action, compute the one-sigma percentage uncertainty
        on chi2 correctly; the expression used until now, 100*sqrt(2/dof),
        is only valid when all data are weighted equally.

Modified 2005 March 17 by CM:
    For the "fit" action with parallel processing, check that root receives
        the responses to the correct broadcast
    Compute degrees of freedom in routine read_dat rather than here, so
        that it only needs to be done once per fit rather than repeatedly
        (dof for individual data types is still computed here)
    Allow weights and degrees of freedom to be floating-point rather than
        integer; but if they are integer after all, print integer-rounded
        values to the screen rather than floating-point
    For the "write" action with "speckle" turned on, write Doppler file
        fit_MM_NN.dat with spectral values normalized to the input sdev
        value, not to the sdev value increased for self-noise.

Modified 2005 March 6 by CM:
    For the "write" action, write calibration factors for plane-of-sky
        datasets to disk if the "list_posetcal" parameter is turned on

Modified 2005 March 2 by CM:
    Rename some "sdev" and "var" values to be "oneovervar" (1/variance)
    For the "write" action, adjust arguments of revised "resampim" routine
    For the "write" action, rotate plane-of-sky fit/obs/res frames so that
        north is upward, unless poset_scaling = NONE

Modified 2005 February 22 by CM:
    For the "write" action, fix bug (incorrect array dimensions) in
        scaling fit vs. obs pgm image brightness for plane-of-sky datasets
    For the "write" action, add the new "image_rebin" argument to function
        resampim to handle output images which have much coarser resolution
        than the model POS frames from which they are constructed
        (i.e., which are greatly undersampled).  This situation probably
        won't arise often: The real importance of "image_rebin" is for
        dealing with plane-of-sky fit frames in routine calc_fits.

Modified 2005 January 25 by CM:
    Take care of uninitialized variables

Modified 2005 January 20 by CM:
    For the "write" action, implement the bilinear and bicubic
        interpolation options for the "dd_scaling" parameter
    For the "write" action, write model, data, and residual pgm images
        for POS datasets
    Add pgm file output (obs/fit/res) for POS datasets
    Correct the expression for chi-square for POS datasets
        (the calibration factor was being ignored)
    If the calibration factor for a POS frame is a floating parameter,
        don't allow it to be negative: instead set it to zero.  Also
        display a warning, unless this is a fit and we're not at the end
        of an iteration.

Modified 2005 January 12 by CM:
    For the "fit" action with parallel processing, revise the code so
        that "list_breakdown" will work: For each dataset which is handled
        by a branch node rather than by root, root broadcasts a request for
        that branch node to run chi2 for just that one dataset and to send
        the results to root.  This is necessary because otherwise each
        branch node calls chi2 just once to process *all* of its datasets
        in one shot, without keeping results for different data types
        (Doppler, delay-Doppler, etc.) separate.

Modified 2004 December 1 by CM:
    For the "write" action, adjust the "listfit" option so that
        the "listfit_path" directory can be created on the fly
        if necessary

Modified 2004 November 30 by CM:
    For the "write" action, implement the "listfit" option to write
        out the model "data" files for all delay-Doppler frames

Modified 2004 July 26 by CM:
    If the calibration factor for a Doppler or delay-Doppler frame
        is a floating parameter, don't allow it to be negative: instead
        set it to zero.  Also display a warning, unless this is a fit
        and we're not at the end of an iteration.

Modified 2004 June 27 by CM:
    Fixed bug: For 'dd_scaling' = 'block' the number of pixels per side
        for the resampled delay-Doppler obs/fit/residual images wasn't
        being defined

Modified 2004 April 3 by CM:
    Add the "list_breakdown" argument so that we can list chi2 by
        data type (Doppler, delay-Doppler, POS, lightcurves) for
        both the "fit" and the "write" action as desired

Modified 2004 March 20 by CM:
    For the "write" action, move final summary screen display from
        the main program (shape.c) to here
    For the "write" action, display summaries for each form of data
        (Doppler, delay-Doppler, POS, lightcurves) taken separately,
        in addition to the overall summary

Modified 2004 February 29 by CM:
    For lightcurve output, replace JD244 variable with jdoffset parameter
    For the "write" action, move all screen and file output which
        relies on updated calibration factors to this routine

Modified 2003 May 10 by CM:
    Account for contributions to chi-square from model (delay-)Doppler
        power which lies outside the data frame

Modified 2003 April 24 by CM:
    Display chi-square for Doppler datasets when action = "write"
    Implement the new "weight" parameter
 *****************************************************************************************/
extern "C" {
#include "head.h"
}
__device__ int c2_print_breakdown, dof;
__device__ unsigned char c2_write_chi2fit0, c2_badradar,
		c2_badphoto, c2_baddopscale, c2_badposet, c2_posbnd, c2_baddiam;
__device__ double c2_dof_deldop, c2_dof_doppler, c2_dof_poset, c2_dof_lghtcrv,
					c2_dof;

__host__ double chi2_deldop_cuda_af(struct par_t *dpar, struct dat_t *ddat, int s,
		int list_breakdown, double *chi2_all_deldop, double *chi2_fit0_deldop,
		double *dof_fit0_deldop);
__host__ double chi2_doppler_cuda_af(struct par_t *dpar, struct dat_t *ddat, int s,
		int list_breakdown, double *chi2_all_doppler, double *chi2_fit0_doppler,
		double *dof_fit0_doppler);

__device__ double c2_chi2, c2_chi2_frame, c2_chi2_set, c2_weight,
	c2_chi2_all_doppler, chi2_fit0_doppler, dof_fit0_doppler, c2_chi2_all_deldop,
	chi2_fit0_deldop, dof_fit0_deldop;
__device__ int c2_nsets, c2_nframes, c2_ndel, c2_ndop;
__device__ unsigned char c2_type;
__device__ float o2, m2, om;
__device__ struct deldop_t *c2_deldop;
__device__ struct doppler_t *c2_doppler;
__device__ float dbg_sum_fit=0.0, dbg_sum_obs=0.0, dbg_sum_oov=0.0;
__device__ double sum_oov=0.0, sum_oovs=0.0, dom=0.0, dm2=0.0, do2=0.0, dopo2=0.0,
		dopm2=0.0, dopom=0.0, sum_fit=0.0, sum_obs=0.0, sum_one=0.0;

__global__ void c2_init_krnl(struct dat_t *ddat) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		ddat->chi2 =  0.0;
		c2_nsets = ddat->nsets;
		c2_chi2_set = 0.0;
		c2_weight = 0.0;
		c2_chi2_all_doppler = 0.0;
	}
}
__global__ void c2_get_type_krnl(struct dat_t *ddat, int s) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0)
		c2_type = ddat->set[s].type;
}
__global__ void c2_add_chi2_krnl(struct dat_t *ddat, int s) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0)
		ddat->chi2 += ddat->set[s].chi2;
}
__global__ void c2_retrieve_chi2_krnl(struct dat_t *ddat) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0)
		c2_chi2 = ddat->chi2;
}
__global__ void c2_get_frames_krnl(struct dat_t *ddat, int s) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		switch (c2_type) {
		case DELAY:
			c2_deldop = &ddat->set[s].desc.deldop;
			c2_nframes = c2_deldop->nframes;
			c2_ndel = c2_deldop->ndel;
			c2_ndop = c2_deldop->ndop;

			break;
		case DOPPLER:
			c2_doppler = &ddat->set[s].desc.doppler;
			c2_nframes = c2_doppler->nframes;
			c2_ndel = NULL;
			c2_ndop = c2_doppler->ndop;
			break;
		}
	}
}
__global__ void c2_get_ndop_krnl(int f) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		switch (c2_type) {
		case DELAY:
			c2_ndel = c2_deldop->frame[f].ndel;
			c2_ndop = c2_deldop->frame[f].ndop;
			/* Initialize contributions to chi2 to values that account for overflow
			 * beyond limits of data frame. Contributions were computed by
			 * pos2deldop_cuda_2.               */
			o2 = c2_deldop->frame[f].overflow_o2;
			m2 = c2_deldop->frame[f].overflow_m2;
			om = 0.0;
			c2_weight = c2_deldop->frame[f].weight;
			break;
		case DOPPLER:
			c2_ndel = NULL;
			c2_ndop = c2_doppler->frame[f].ndop;
			o2 = c2_doppler->frame[f].overflow_o2;
			m2 = c2_doppler->frame[f].overflow_m2;
			om = 0.0;
			c2_weight = c2_doppler->frame[f].weight;
			break;
		}

	}
}
__global__ void c2_deldop_add_o2_krnl(struct par_t *dpar, struct
		dat_t *ddat, int s, int f) {
	/* ndel*ndop-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % c2_ndel + 1;
	int j = offset / c2_ndel + 1;
	float temp;

	if (offset < (c2_ndel*c2_ndop)) {
		/* The following two lines implement this:
		 * o2 += obs[i][j]*obs[i][j]*oneovervar[i][j];	 */
		temp = c2_deldop->frame[f].obs[i][j] * c2_deldop->frame[f].obs[i][j] *
				c2_deldop->frame[f].oneovervar[i][j];
		atomicAdd(&o2, temp);

		/* The following two lines implement this:
		 * m2 += fit[i][j]*fit[i][j]*oneovervar[i][j];		 */
		temp = c2_deldop->frame[f].fit_s[offset] * c2_deldop->frame[f].fit_s[offset] *
				c2_deldop->frame[f].oneovervar[i][j];
		atomicAdd(&m2, temp);

		/* The following two lines implement this:
		 * om += fit[i][j]*obs[i][j]*oneovervar[i][j];		 */
		temp = c2_deldop->frame[f].fit_s[offset] * c2_deldop->frame[f].obs[i][j] *
				c2_deldop->frame[f].oneovervar[i][j];
		atomicAdd(&om, temp);
	}
	__syncthreads();
}
__global__ void c2_deldop_dbg_krnl(struct par_t *dpar, struct dat_t *ddat,
		int s, int f) {
	/* ndel*ndop-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;

	if (offset < (c2_ndel*c2_ndop)) {
		atomicAdd(&dbg_sum_fit, (float)c2_deldop->frame[f].fit_s[offset]);
//		atomicAdd(&dbg_sum_obs, (float)c2_deldop->frame[f].obs_s[offset]);
//		atomicAdd(&dbg_sum_oov, (float)c2_deldop->frame[f].oneovervar_s[offset]);
		atomicAdd(&c2_occ1, 1);
	}
}

__global__ void c2_deldop_add_o2_serial_krnl(struct par_t *dpar,
		struct dat_t *ddat, int s, int f) {
	/* Single threaded kernel for debugging */
	if (threadIdx.x == 0) {
		double temp;
		int ndel = c2_ndel;
		int ndop = c2_ndop;
		int i, j, off;
		for (i=1; i<=ndel; i++) {
			for (j=1; j<=ndop; j++){
				off = (j-1)*ndel + (i-1);
				temp = c2_deldop->frame[f].obs[i][j] * c2_deldop->frame[f].obs[i][j] *
						c2_deldop->frame[f].oneovervar[i][j];
				do2 += temp;
				temp = c2_deldop->frame[f].fit_s[off] * c2_deldop->frame[f].fit_s[off] *
						c2_deldop->frame[f].oneovervar[i][j];
				dm2 += temp;
				temp = c2_deldop->frame[f].fit_s[off] * c2_deldop->frame[f].obs[i][j] *
						c2_deldop->frame[f].oneovervar[i][j];
				dom += temp;
			}
		}
	}
	__syncthreads();
}

__global__ void c2_add_deldop_contributions_krnl(struct par_t *dpar, struct
		dat_t *ddat, int s, int f) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {

		double calval, err, o2_fit0, dof_fit0, err_fit0, thresh_fit0;
		c2_chi2_frame = 0.0;
		int off, i, j;

		/* If this frame's calibration factor is allowed to float, set it to
		 * minimize chi-square, the sum over all pixels of
		 * 		        { (obs - calfact*fit)^2 / variance }              */

		if (c2_deldop->frame[f].cal.state == 'f') {
			if (om > 0.0)	c2_deldop->frame[f].cal.val = om/m2;
			else			c2_deldop->frame[f].cal.val = TINYCALFACT;
		}

		/*  Compute chi-square for this frame  */
		calval = c2_deldop->frame[f].cal.val;
		err = c2_weight*(o2 - 2*calval*om + calval*calval*m2);
		c2_deldop->frame[f].chi2 = err;
		c2_chi2_frame += err;

		/* Compute chi-square contributions and deg. of freedom due to pixels
		 * whose model signal is less than or equal to 'chi2fit0_thresh'
		 * standard deviations of the noise in the data frame   */
		o2_fit0 = 0.0;
		dof_fit0 = 0.0;
		err_fit0 = 0.0;
		thresh_fit0 = dpar->chi2fit0_thresh * c2_deldop->frame[f].sdev;
		if (dpar->write_chi2fit0) {
			for (i=0; i<c2_ndel; i++)
				for (j=0; j<c2_ndop; j++)
					off = j*c2_ndel + i; // For the unrolled fit, obs pointers
					if (calval*c2_deldop->frame[f].fit_s[off] <= thresh_fit0) {
						o2_fit0 += c2_deldop->frame[f].obs[i][j]*
								c2_deldop->frame[f].obs[i][j]*
								c2_deldop->frame[f].oneovervar[i][j];
						if (c2_deldop->frame[f].oneovervar[i][j] > 0.0)
							dof_fit0 += c2_weight;
					}
			err_fit0 = c2_weight*o2_fit0;
			chi2_fit0_deldop += err_fit0;
			dof_fit0_deldop += dof_fit0;
		}
	}
}

__global__ void c2_add_doppler_contributions_serial_krnl(struct par_t *dpar,
		struct dat_t *ddat, int s, int f) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		int j;
		double calval, err, o2_fit0, dof_fit0, err_fit0, thresh_fit0;
		c2_chi2_frame = 0.0;

		/* If this frame's calibration factor is allowed to float, set it to
		 * minimize chi-square, the sum over all bins of
		        { (obs - calfact*fit)^2 / variance }.                     */

		if (c2_doppler->frame[f].cal.state == 'f') {
			if (om > 0.0) 	c2_doppler->frame[f].cal.val = om/m2;
			else {
				c2_doppler->frame[f].cal.val = TINYCALFACT;
				if (dpar->action != FIT )
					printf("WARNING: set %2d frame %2d had negative calfact reset to %10.4e\n",
							s, f, c2_doppler->frame[f].cal.val);	}
		}

		/* Compute chi-square for this frame  */
		calval = c2_doppler->frame[f].cal.val;
		err = c2_weight*(o2 - 2*calval*om + calval*calval*m2);
		c2_doppler->frame[f].chi2 = err;
		c2_chi2_frame += err;
		//if (list_breakdown)
		c2_chi2_all_doppler += err;

		/* Compute chi-square contributions and dof due to bins whose model
		 * signal is =< 'chi2fit0_thresh' standard deviations of the noise in
		 * the data frame   */
		o2_fit0 = 0.0;
		dof_fit0 = 0.0;
		err_fit0 = 0.0;
		thresh_fit0 = dpar->chi2fit0_thresh * c2_doppler->frame[f].sdev;
		if (dpar->write_chi2fit0) {
			for (j=0; j<c2_ndop; j++)
				if (calval*c2_doppler->frame[f].fit[j] <= thresh_fit0) {
					o2_fit0 += c2_doppler->frame[f].obs[j]*c2_doppler->frame[f].obs[j]
					  *c2_doppler->frame[f].oneovervar[j];
					if (c2_doppler->frame[f].oneovervar[j] > 0.0)
						dof_fit0 += c2_weight;
				}
			err_fit0 = c2_weight*o2_fit0;
			chi2_fit0_doppler += err_fit0;
			dof_fit0_doppler += dof_fit0;
		}
	}
}
__global__ void c2_doppler_add_o2_krnl(struct par_t *dpar, struct
		dat_t *ddat, int s, int f) {//, double *atomo2, double
	//	*atomm2,  double *atomom) {
	/* ndop-threaded kernel */
	int j = blockIdx.x * blockDim.x + threadIdx.x;
	float temp;

	if (j < c2_ndop) {
		/* Add contributions from power within frame limits */

		/* Next 2 lines implement: o2 += obs[j]*obs[j]*oneovervar[j];	 */
		temp = c2_doppler->frame[f].obs[j] * c2_doppler->frame[f].obs[j] *
				c2_doppler->frame[f].oneovervar[j];
		atomicAdd(&o2, temp);
		//atomo2[j+1] = temp;


		/* Next 2 lines implement: m2 += fit[j]*fit[j]*oneovervar[j];		 */
		temp = c2_doppler->frame[f].fit_s[j] * c2_doppler->frame[f].fit_s[j] *
				c2_doppler->frame[f].oneovervar[j];
		atomicAdd(&m2, temp);
		//atomm2[j+1] = temp;

		/* Next 2 lines implement: om += fit[j]*obs[j]*oneovervar[j];		 */
		temp = c2_doppler->frame[f].fit_s[j] * c2_doppler->frame[f].obs[j] *
				c2_doppler->frame[f].oneovervar[j];
		atomicAdd(&om, temp);
		//atomom[j+1] = temp;
	}
	__syncthreads();
}
__global__ void c2_doppler_add_o2_serial_krnl(int f) {
	/* Single-threaded kernel */
	int j;
	if (threadIdx.x == 0) {
		for (j=0; j<=c2_ndop; j++) {
			dopo2 += c2_doppler->frame[f].obs[j] * c2_doppler->frame[f].obs[j] *
					c2_doppler->frame[f].oneovervar[j];
			dopm2 += c2_doppler->frame[f].fit_s[j] * c2_doppler->frame[f].fit_s[j] *
					c2_doppler->frame[f].oneovervar[j];
			dopom += c2_doppler->frame[f].fit_s[j] * c2_doppler->frame[f].obs[j] *
					c2_doppler->frame[f].oneovervar[j];
			sum_fit += c2_doppler->frame[f].fit_s[j];
			sum_obs += c2_doppler->frame[f].obs[j];
			sum_one += c2_doppler->frame[f].oneovervar[j];
		}
	}
	__syncthreads();
}
__global__ void c2_set_chi2_krnl(struct dat_t *ddat, double chi2, int s) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0)
		ddat->set[s].chi2 = chi2;
}
__global__ void c2_dbg3_krnl(struct dat_t *ddat, double **obs, double **oov,
		float *fit, int s, int f) {
	/* ndel*ndop-threaded kernel */
	int off = blockIdx.x * blockDim.x + threadIdx.x;
	int idel = off % c2_ndel + 1;
	int idop = off / c2_ndel + 1;
	if (off < (c2_ndel*c2_ndop)) {
		fit[off] = ddat->set[s].desc.deldop.frame[f].fit_s[off];
		obs[idel][idop] = ddat->set[s].desc.deldop.frame[f].obs[idel][idop];
		oov[idel][idop] = ddat->set[s].desc.deldop.frame[f].oneovervar[idel][idop];
	}
	__syncthreads();
}
__global__ void c2_dbg4_krnl(double *obs, double *oov,	double *fit, int f) {
	/* ndel*ndop-threaded kernel */
	int off = blockIdx.x * blockDim.x + threadIdx.x;
	int idop = off + 1;
	if (off < c2_ndop) {
		fit[idop] = c2_doppler->frame[f].fit_s[idop];
		obs[idop] = c2_doppler->frame[f].obs[idop];
		oov[idop] = c2_doppler->frame[f].oneovervar[idop];
	}
	__syncthreads();
}
__global__ void c2_get_prntflgs_krnl(struct par_t *dpar, struct dat_t *ddat) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		c2_print_breakdown = (ddat->dof_deldop > SMALLVAL || ddat->dof_doppler  > SMALLVAL
				|| ddat->dof_poset    > SMALLVAL	|| ddat->dof_lghtcrv  > SMALLVAL);
		c2_dof_deldop = ddat->dof_deldop;
		c2_dof_doppler = ddat->dof_doppler;
		c2_dof_poset = ddat->dof_poset;
		c2_dof_lghtcrv = ddat->dof_lghtcrv;
		c2_dof = ddat->dof;
		c2_write_chi2fit0 = dpar->write_chi2fit0;
		c2_badradar = dpar->badradar;
		c2_badphoto = dpar->badphoto;
		c2_baddopscale = dpar->baddopscale;
		c2_badposet = dpar->badposet;
		c2_posbnd = dpar->posbnd;
		c2_baddiam = dpar->baddiam;
	}
}

__host__ double chi2_cuda_af(struct par_t *dpar, struct dat_t *ddat, int list_breakdown)
{
	int s, nsets, print_breakdown;
	unsigned char write_chi2fit0, baddiam, badphoto, badposet, baddopscale,
			posbnd, badradar;
	unsigned char type;
	dim3 BLK,THD;

	double chi2_all_doppler, chi2_all_deldop, /*chi2_all_poset, chi2_all_lghtcrv,*/
		   chi2_fit0_doppler, chi2_fit0_deldop, /*chi2_fit0_poset, chi2_branch, */
		   dof_fit0_doppler, dof_fit0_deldop, /*dof_fit0_poset, */chi2, dof;
	double dof_deldop, dof_doppler/*, dof_poset, dof_lghtcrv*/;
	char dofstring[MAXLEN], dof0string[MAXLEN];

	/*  Initialize variables that accumulate chi-square values  */
	chi2_all_deldop = chi2_all_doppler = /*chi2_all_poset = chi2_all_lghtcrv = */
			chi2_fit0_deldop = chi2_fit0_doppler /*= chi2_fit0_poset*/ = 0.0;
	dof_fit0_deldop = dof_fit0_doppler/* = dof_fit0_poset*/ = 0.0;

	/* Initialize variables that accumulate chi-square values  */
	c2_init_krnl<<<1,1>>>(ddat);
	checkErrorAfterKernelLaunch("c2_init_krnl, chi2_cuda");
	gpuErrchk(cudaMemcpyFromSymbol(&nsets, c2_nsets, sizeof(int),
			0, cudaMemcpyDeviceToHost));

	/* Any references to MPI stuff has been removed
	 * for this "FIT" only cuda version.   */
	/* Loop through all datasets, carry out chi-square computations, and
	 * provide screen and image output                            */
	for (s=0; s<nsets; s++) {
		/* Launch single-threaded kernel to get type */
		c2_get_type_krnl<<<1,1>>>(ddat, s);
		checkErrorAfterKernelLaunch("c2_get_type_krnl, chi2_cuda");
		gpuErrchk(cudaMemcpyFromSymbol(&type, c2_type, sizeof(unsigned char),
				0, cudaMemcpyDeviceToHost));

		switch (type) {
		case DELAY:
			chi2 = chi2_deldop_cuda_af(dpar, ddat, s, list_breakdown,
				   &chi2_all_deldop, &chi2_fit0_deldop, &dof_fit0_deldop);
			c2_set_chi2_krnl<<<1,1>>>(ddat, chi2, s);
			checkErrorAfterKernelLaunch("c2_set_chi2_krnl, chi2_cuda");
			break;
		case DOPPLER:
			chi2 = chi2_doppler_cuda_af(dpar, ddat, s, list_breakdown,
				   &chi2_all_doppler, &chi2_fit0_deldop, &dof_fit0_doppler);
			c2_set_chi2_krnl<<<1,1>>>(ddat, chi2, s);
			checkErrorAfterKernelLaunch("c2_set_chi2_krnl, chi2_cuda");
			break;
		case POS:
			printf("\nWrite chi2_poset_cuda!\n");
			//			dat->set[s].chi2 = chi2_poset(dpar, s);
			break;
		case LGHTCRV:
			printf("\nWrite chi2_lghtcrv_cuda!\n");
			//			dat->set[s].chi2 = chi2_lghtcrv(dpar, s);
			break;
		default:
			printf("chi2_cuda.cu: can't handle this type yet\n");
		}

		/* Single-thread kernel adds ddat->set[s].chi2 to ddat->chi2 */
		c2_add_chi2_krnl<<<1,1>>>(ddat, s);
		checkErrorAfterKernelLaunch("c2_add_chi2_krnl, chi2_cuda");

	}  /* end for loop over datasets */
	/* Launch single-threaded kernel to retrieve ddat->chi2 to return it */
	c2_retrieve_chi2_krnl<<<1,1>>>(ddat);
	checkErrorAfterKernelLaunch("c2_retrieve_chi2_krnl, chi2_cuda");
	gpuErrchk(cudaMemcpyFromSymbol(&chi2, c2_chi2, sizeof(double),
			0, cudaMemcpyDeviceToHost));
	/*.......................................................................*/


	/* Call kernel to get flags from ddat */
	c2_get_prntflgs_krnl<<<1,1>>>(dpar, ddat);
	checkErrorAfterKernelLaunch("c2_get_prntflgs_krnl");
	gpuErrchk(cudaMemcpyFromSymbol(&print_breakdown, c2_print_breakdown,
			sizeof(int), 0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&dof_deldop, c2_dof_deldop, sizeof(double),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&dof_doppler, c2_dof_doppler, sizeof(double),
			0, cudaMemcpyDeviceToHost));
//	gpuErrchk(cudaMemcpyFromSymbol(&dof_poset, c2_dof_poset, sizeof(double),
//			0, cudaMemcpyDeviceToHost));
//	gpuErrchk(cudaMemcpyFromSymbol(&dof_lgthcrv, c2_dof_lgthcrvr, sizeof(double),
//			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&write_chi2fit0, c2_write_chi2fit0,
			sizeof(unsigned char), 0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&dof, c2_dof, sizeof(double),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&baddiam, c2_baddiam, sizeof(unsigned char),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&badphoto, c2_badphoto, sizeof(unsigned char),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&posbnd, c2_posbnd, sizeof(unsigned char),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&badposet, c2_badposet, sizeof(unsigned char),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&badradar, c2_badradar, sizeof(unsigned char),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&baddopscale, c2_baddopscale, sizeof(unsigned char),
			0, cudaMemcpyDeviceToHost));

	if (list_breakdown) {
		if (print_breakdown) {
			printf("#\n");
			if (dof_deldop > SMALLVAL) {
				intifpossible( dofstring, MAXLEN, dof_deldop, SMALLVAL, "%f");
				printf("delay   chi2 = %e for %s dof (reduced chi2 = %f)\n",
						chi2_all_deldop, dofstring, chi2_all_deldop/dof_deldop);
				if (write_chi2fit0) {
					intifpossible( dof0string, MAXLEN, dof_fit0_deldop, SMALLVAL, "%f");
					printf("              (%e outside model for %s dof)\n",
							chi2_fit0_deldop, dof0string);
				}
			}
			if (dof_doppler > SMALLVAL) {
				intifpossible( dofstring, MAXLEN, dof_doppler, SMALLVAL, "%f");
				printf("Doppler chi2 = %e for %s dof (reduced chi2 = %f)\n",
						chi2_all_doppler, dofstring, chi2_all_doppler/dof_doppler);
				if (write_chi2fit0) {
					intifpossible( dof0string, MAXLEN, dof_fit0_doppler, SMALLVAL, "%f");
					printf("              (%e outside model for %s dof)\n",
							chi2_fit0_doppler, dof0string);
				}
			}
			//		if (dof_poset > SMALLVAL) {
			//			intifpossible( dofstring, MAXLEN, dof_poset, SMALLVAL, "%f");
			//			printf("POS     chi2 = %e for %s dof (reduced chi2 = %f)\n",
			//					chi2_all_poset, dofstring, chi2_all_poset/dof_poset);
			//			if (write_chi2fit0) {
			//				intifpossible( dof0string, MAXLEN, dof_fit0_poset, SMALLVAL, "%f");
			//				printf("              (%e outside model for %s dof)\n",
			//						chi2_fit0_poset, dof0string);
			//			}
			//		}
			//		if (dof_lghtcrv > SMALLVAL) {
			//			intifpossible( dofstring, MAXLEN, dof_lghtcrv, SMALLVAL, "%f");
			//			printf("lghtcrv chi2 = %e for %s dof (reduced chi2 = %f)\n",
			//					chi2_all_lghtcrv, dofstring, chi2_all_lghtcrv/dof_lghtcrv);
			//		}
			intifpossible( dofstring, MAXLEN, dof, SMALLVAL, "%f");
			printf("ALLDATA chi2 = %e for %s dof (reduced chi2 = %f)",
					chi2, dofstring, chi2/dof);
		} else {
			intifpossible( dofstring, MAXLEN, ddat->dof, SMALLVAL, "%f");
			printf("        chi2 = %e for %s dof (reduced chi2 = %f)",
					chi2, dofstring, chi2/dof);
		}

	if (baddiam)		printf("  (BAD DIAMS)");
	if (badphoto)		printf("  (BAD PHOTO) (chi2_cuda)");
	if (posbnd)			printf("  (BAD POS)");
	if (badposet)		printf("  (BAD POSET)");
	if (badradar)		printf("  (BAD RADAR)");
	if (baddopscale)	printf("  (BAD DOPSCALE)");
	printf("\n");
	if (print_breakdown &&  write_chi2fit0) {
		intifpossible( dof0string, MAXLEN,
				dof_fit0_deldop + dof_fit0_doppler/* + dof_fit0_poset*/,
				SMALLVAL, "%f");
		printf("              (%e outside model for %s dof)\n",
				chi2_fit0_deldop + chi2_fit0_doppler/* + chi2_fit0_poset*/, dof0string);
	}
	printf("#\n");
	fflush(stdout);
	}
	/*.......................................................................*/

	return chi2;
}

__host__ double chi2_deldop_cuda_af(struct par_t *dpar, struct dat_t *ddat, int s,
		int list_breakdown, double *chi2_all_deldop, double *chi2_fit0_deldop,
		double *dof_fit0_deldop)
{
	int f, ndel, ndop, nframes;
	double chi2_set, chi2_frame;
	dim3 BLK,THD;
	chi2_set = 0.0;

	/* Launch single-threaded kernel to get # of frames for this set */
	c2_get_frames_krnl<<<1,1>>>(ddat, s);
	checkErrorAfterKernelLaunch("c2_get_frames_krnl, chi2_cuda");
	gpuErrchk(cudaMemcpyFromSymbol(&nframes, c2_nframes, sizeof(int),
			0, cudaMemcpyDeviceToHost));

	/*  Loop through all frames for this dataset  */
	for (f=0; f<nframes; f++) {
		/* Launch single-threaded kernel to get ndel/ndop for this frame */
		c2_get_ndop_krnl<<<1,1>>>(s);
		checkErrorAfterKernelLaunch("c2_get_frames_krnl, chi2_cuda");
		gpuErrchk(cudaMemcpyFromSymbol(&ndel, c2_ndel, sizeof(int),
				0, cudaMemcpyDeviceToHost));
		gpuErrchk(cudaMemcpyFromSymbol(&ndop, c2_ndop, sizeof(int),
				0, cudaMemcpyDeviceToHost));

		/* Add contributions from power within limits of data frame.
		 * This kernel also takes care of the frame's calibration factor and
		 * computes chi2 for this frame */
		BLK.x = floor((maxThreadsPerBlock - 1 + (ndel*ndop)) / maxThreadsPerBlock);
		THD.x = maxThreadsPerBlock; // Thread block dimensions
		c2_deldop_add_o2_krnl<<<BLK,THD>>>(dpar, ddat, s, f);
		checkErrorAfterKernelLaunch("c2_add_deldop_contributions_krnl, line ");

		c2_add_deldop_contributions_krnl<<<1,1>>>(dpar, ddat, s, f);
		checkErrorAfterKernelLaunch("c2_add_deldop_contributions_krnl, line ");
		gpuErrchk(cudaMemcpyFromSymbol(&chi2_frame, c2_chi2_frame, sizeof(double),
				0, cudaMemcpyDeviceToHost));
		chi2_set += chi2_frame;
		if (list_breakdown)
			*chi2_all_deldop += chi2_frame;
	}  /* end for loop over frames */
	return chi2_set;
}

__host__ double chi2_doppler_cuda_af(struct par_t *dpar, struct dat_t *ddat, int s,
		int list_breakdown, double *chi2_all_doppler, double *chi2_fit0_doppler,
		double *dof_fit0_doppler)
{
	int f, ndop, nframes;
	double chi2_set, chi2_frame;
	dim3 BLK,THD;
	/* Initialize chi-square for dataset  */
	chi2_set = 0.0;

	/* Launch single-threaded kernel to get # of frames for this set */
	c2_get_frames_krnl<<<1,1>>>(ddat, s);
	checkErrorAfterKernelLaunch("c2_get_frames_krnl, chi2_cuda");
	gpuErrchk(cudaMemcpyFromSymbol(&nframes, c2_nframes, sizeof(int),
			0, cudaMemcpyDeviceToHost));

	/*  Loop through all frames for this dataset  */
	for (f=0; f<nframes; f++) {

		/* Launch single-threaded kernel to get ndel/ndop for this frame */
		c2_get_ndop_krnl<<<1,1>>>(s);
		checkErrorAfterKernelLaunch("c2_get_ndop_krnl, chi2_cuda");
		gpuErrchk(cudaMemcpyFromSymbol(&ndop, c2_ndop, sizeof(int),
				0, cudaMemcpyDeviceToHost));
		chi2_frame = 0.0;
		/* Add contributions from power within limits of data frame.
		 * This kernel also takes care of the frame's calibration factor and
		 * computes chi2 for this frame */
		BLK.x = floor((maxThreadsPerBlock - 1 + ndop) / maxThreadsPerBlock);
		THD.x = maxThreadsPerBlock; // Thread block dimensions
		c2_doppler_add_o2_krnl<<<BLK,THD>>>(dpar, ddat, s, f);//, atomo2, atomm2, atomom);
		checkErrorAfterKernelLaunch("c2_add_o2_krnl, line ");

		c2_add_doppler_contributions_serial_krnl<<<1,1>>>(dpar, ddat, s, f);
		checkErrorAfterKernelLaunch("c2_add_doppler_contributions_serial_krnl, line ");
		gpuErrchk(cudaMemcpyFromSymbol(&chi2_frame, c2_chi2_frame, sizeof(double),
				0, cudaMemcpyDeviceToHost));
		if (list_breakdown)
					*chi2_all_doppler += chi2_frame;

		chi2_set += chi2_frame;
	}
	return chi2_set;
}


//__host__ double chi2_poset_cuda(struct par_t *dpar, struct poset_t *poset,
//		int s, double *chi2_all_poset, double *chi2_fit0_poset,
//		double *dof_fit0_poset)
//{
//	int f, i, j, n_pix, n_pos, nrow, ncol;
//	double chi2_set, err, err_fit0, o2, m2, om, o2_fit0, calval, fit255, obs255,
//	xoff, yoff, resamp_fact, resamp_x0, resamp_y0, resamp_width,
//	resamp_angle, weight, dof, dof_fit0, thresh_fit0;
//	double **obs, **fit, **res, **oneovervar, **resamp_fit, **resamp_obs,
//	**resamp_res;
//
//	/*  Initialize chi-square for dataset  */
//
//	chi2_set = 0.0;
//
//	/*  Loop through all frames for this dataset  */
//
//	for (f=0; f<poset->nframes; f++) {
//		ncol = poset->frame[f].ncol;
//		nrow = poset->frame[f].nrow;
//		obs = poset->frame[f].obs.b;
//		fit = poset->frame[f].fit.b;
//		oneovervar = poset->frame[f].oneovervar;  /* 1/variance */
//		weight = poset->frame[f].weight;
//		dof = poset->frame[f].dof;
//
//		o2 = m2 = om = 0.0;
//		for (i=1; i<=ncol; i++)
//			for (j=1; j<=nrow; j++) {
//				o2 += obs[i][j]*obs[i][j]*oneovervar[i][j];
//				m2 += fit[i][j]*fit[i][j]*oneovervar[i][j];
//				om += fit[i][j]*obs[i][j]*oneovervar[i][j];
//			}
//
//		/*  The calibration factor always floats for plane-of-sky frames:
//        set it to minimize chi-square, the sum over all pixels of
//        { (obs - calfact*fit)^2 / variance }.                          */
//
//		if (om > 0.0)
//			poset->frame[f].cal.val = om/m2;
//		else
//			poset->frame[f].cal.val = TINYCALFACT;
//
//		/*  Compute chi-square for this frame  */
//
//		calval = poset->frame[f].cal.val;
//		err = weight*(o2 - 2*calval*om + calval*calval*m2);
//		poset->frame[f].chi2 = err;
//		chi2_set += err;
//
//		/*  Compute the chi-square contributions and number of degrees of freedom
//        due to pixels whose model signal is less than or equal to
//        'chi2fit0_thresh' standard deviations of the noise in the data frame   */
//
//		o2_fit0 = 0.0;
//		dof_fit0 = 0.0;
//		err_fit0 = 0.0;
//		thresh_fit0 = par->chi2fit0_thresh;  /* "sdev" = 1.0 for plane-of-sky data */
//
//	}  /* end for loop over frames */
//
//	return chi2_set;
//}
//
//
//__host__ double chi2_lghtcrv_cuda(struct par_t *dpar, struct lghtcrv_t *lghtcrv,
//		int s, double *chi2_all_lghtcrv)
//{
//	int i, n, ncalc;
//	double chi2_set, err, o2, m2, om, calval, weight, dof, obsmag, fitmag, obsmagerr;
//
//	n = lghtcrv->n;
//	ncalc = lghtcrv->ncalc;
//
//	/*  Compute contributions to chi-square  */
//
//	o2 = m2 = om = 0.0;
//	for (i=1; i<=n; i++) {
//		o2 += lghtcrv->obs[i] * lghtcrv->obs[i] * lghtcrv->oneovervar[i];
//		m2 += lghtcrv->fit[i] * lghtcrv->fit[i] * lghtcrv->oneovervar[i];
//		om += lghtcrv->fit[i] * lghtcrv->obs[i] * lghtcrv->oneovervar[i];
//	}
//
//	/*  If this lightcurve's calibration factor is allowed to float,
//      set it to minimize chi-square, the sum over all points of
//      { (obs - calfact*fit)^2 / variance }.                         */
//dpar
//	if (lghtcrv->cal.state == 'f') {
//		if (om > 0.0)
//			lghtcrv->cal.val = om/m2;
//		else
//			lghtcrv->cal.val = TINYCALFACT;
//	}
//
//	/*  Compute chi-square for dataset  */
//
//	calval = lghtcrv->cal.val;
//	weight = lghtcrv->weight;
//	dof = lghtcrv->dof;
//	err = weight*(o2 - 2*calval*om + calval*calval*m2);
//	chi2_set = err;
//
//	return chi2_set;
//}
