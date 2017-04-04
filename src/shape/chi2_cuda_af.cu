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
__device__ int c2af_print_breakdown, c2af_nframes;

__device__ unsigned char c2af_write_chi2fit0, c2af_badradar,
		c2af_badphoto, c2af_baddopscale, c2af_badposet, c2af_posbnd,
		c2af_baddiam;

__device__ double c2af_dof_deldop, c2af_dof_doppler, c2af_dof_poset,
	c2af_dof_lghtcrv, c2af_dof, c2af_chi2_all_doppler, c2af_chi2_fit0_doppler,
	c2af_dof_fit0_doppler, c2af_chi2_fit0_deldop, c2af_dof_fit0_deldop, c2af_chi2;
__device__ unsigned char c2af_type;



__host__ double chi2_deldop_cuda_af(struct par_t *dpar, struct dat_t *ddat, int s,
		int list_breakdown, double *chi2_all_deldop, double *chi2_fit0_deldop,
		double *dof_fit0_deldop);
__host__ double chi2_doppler_cuda_af(struct par_t *dpar, struct dat_t *ddat, int s,
		int list_breakdown, double *chi2_all_doppler, double *chi2_fit0_doppler,
		double *dof_fit0_doppler);

__global__ void c2af_init_krnl(struct dat_t *ddat) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		ddat->chi2 =  0.0;
		//c2af_chi2_set = 0.0;
		//c2_weight = 0.0;
		c2af_chi2_all_doppler = 0.0;
	}
}
__global__ void c2af_get_type_krnl(struct dat_t *ddat, int s) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0)
		c2af_type = ddat->set[s].type;
}
__global__ void c2af_add_chi2_krnl(struct dat_t *ddat, int s) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0)
		ddat->chi2 += ddat->set[s].chi2;
}
__global__ void c2af_retrieve_chi2_krnl(struct dat_t *ddat) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0)
		c2af_chi2 = ddat->chi2;
}
__global__ void c2af_get_frames_krnl(
		struct dat_t *ddat,
		int s) {
	/* Single-threaded kernel */
	if (threadIdx.x ==0) {
		switch (c2af_type) {
		case DELAY:
			c2af_nframes = ddat->set[s].desc.deldop.nframes;
			break;
		case DOPPLER:
			c2af_nframes = ddat->set[s].desc.doppler.nframes;
			break;
		}
	}
}
__global__ void c2af_deldop_set_shortcuts_krnl(
		struct dat_t *ddat,
		struct deldopfrm_t **frame,
		int *ndel,
		int *ndop,
		float4 *o2_m2_om_wt,
		int s,
		int nframes) {
	/* nframes-threaded kernel */
	int frm = threadIdx.x;

	if (frm < nframes) {

		frame[frm] = &ddat->set[s].desc.deldop.frame[frm];
		ndel[frm] = frame[frm]->ndel;
		ndop[frm] = frame[frm]->ndop;
		/* Initialize contributions to chi2 to values that account for overflow
		 * beyond limits of data frame. Contributions were computed by
		 * pos2deldop_cu_af               */
		o2_m2_om_wt[frm].w = frame[frm]->overflow_o2;
		o2_m2_om_wt[frm].x = frame[frm]->overflow_m2;
		o2_m2_om_wt[frm].y = 0.0;
		o2_m2_om_wt[frm].z = frame[frm]->weight;
	}
}
__global__ void c2af_doppler_set_shortcuts_krnl(
		struct dat_t *ddat,
		struct dopfrm_t **frame,
		int *ndop,
		float4 *o2_m2_om_wt,
		int s,
		int nframes) {
	/* nframes-threaded kernel */
	int frm = threadIdx.x;

	if (frm < nframes) {

		frame[frm] = &ddat->set[s].desc.doppler.frame[frm];
		ndop[frm] = frame[frm]->ndop;
		o2_m2_om_wt[frm].w = frame[frm]->overflow_o2;
		o2_m2_om_wt[frm].x = frame[frm]->overflow_m2;
		o2_m2_om_wt[frm].y = 0.0;
		o2_m2_om_wt[frm].z = frame[frm]->weight;
	}
}
__global__ void c2af_deldop_prep_o2_m2_krnl(
		struct dat_t *ddat,
		struct deldopfrm_t **frame,
		float **temp_o2,
		float **temp_m2,
		float **temp_om,
		int *ndel,
		int *ndop,
		int s,
		int nframes,
		int frame_size) {

	/* ndel*ndop-threaded kernel */
	int total_offset = blockIdx.x * blockDim.x + threadIdx.x;
	int frm = total_offset / frame_size;
	int offset = total_offset % frame_size;
	int i = offset % ndel[0] + 1;
	int j = offset / ndel[0] + 1;

	if ((offset < frame_size) && (frm < nframes)) {
		/* The following two lines implement this:
		 * o2 += obs[i][j]*obs[i][j]*oneovervar[i][j];	 */
		temp_o2[frm][offset] = frame[frm]->obs[i][j] * frame[frm]->obs[i][j] *
						  frame[frm]->oneovervar[i][j];
		//atomicAdd(&o2, temp);

		/* The following two lines implement this:
		 * m2 += fit[i][j]*fit[i][j]*oneovervar[i][j];		 */
		temp_m2[frm][offset] = frame[frm]->fit_s[offset] * frame[frm]->fit_s[offset] *
						 frame[frm]->oneovervar[i][j];
		//atomicAdd(&m2, temp);

		/* The following two lines implement this:
		 * om += fit[i][j]*obs[i][j]*oneovervar[i][j];		 */
		temp_om[frm][offset] = frame[frm]->fit_s[offset] * frame[frm]->obs[i][j] *
						  frame[frm]->oneovervar[i][j];
		//atomicAdd(&om, temp);
	}
}
__global__ void c2af_add_deldop_contributions_krnl(
		struct par_t *dpar,
		struct deldopfrm_t **frame,
		float **temp_o2,
		float **temp_m2,
		float **temp_om,
		double *chi2_frame,
		float4 *o2_m2_om_wt,
		int *ndel,
		int *ndop,
		int s,
		int nframes) {

	/* nframes-threaded kernel */
	int frm = threadIdx.x;
	float calval, err, o2_fit0, dof_fit0, err_fit0, thresh_fit0;
	float o2, m2, om, wt;
	int off, i, j;

	if (frm < nframes) {
		o2 = o2_m2_om_wt[frm].w + temp_o2[frm][0];
		m2 = o2_m2_om_wt[frm].x + temp_m2[frm][0];
		om = o2_m2_om_wt[frm].y + temp_om[frm][0];
		wt = o2_m2_om_wt[frm].z;

		chi2_frame[frm] = 0.0;

		/* If frame's calibration factor is allowed to float, set it to
		 * minimize chi-square, the sum over all pixels of
		 * 		        { (obs - calfact*fit)^2 / variance }              */

		if (frame[frm]->cal.state == 'f') {
			if (om > 0.0)	frame[frm]->cal.val = om/m2;
			else			frame[frm]->cal.val = TINYCALFACT;
		}

		/*  Compute chi-square for this frame  */
		calval = frame[frm]->cal.val;
		err = wt*(o2 - 2*calval*om + calval*calval*m2);
		frame[frm]->chi2 = err;
		chi2_frame[frm] += err;

		/* Compute chi-square contributions and deg. of freedom due to pixels
		 * whose model signal is less than or equal to 'chi2fit0_thresh'
		 * standard deviations of the noise in the data frame   */
		o2_fit0 = dof_fit0 = err_fit0 = 0.0;
		thresh_fit0 = dpar->chi2fit0_thresh * frame[frm]->sdev;
		if (dpar->write_chi2fit0) {
			for (i=0; i<ndel[0]; i++)
				for (j=0; j<ndop[0]; j++)
					off = j*ndel[0] + i; // For the unrolled fit, obs pointers
					if (calval*frame[frm]->fit_s[off] <= thresh_fit0) {
						o2_fit0 += frame[frm]->obs[i][j]*frame[frm]->obs[i][j]*
								frame[frm]->oneovervar[i][j];
						if (frame[frm]->oneovervar[i][j] > 0.0)
							dof_fit0 += wt;
					}
			err_fit0 = wt*o2_fit0;
			c2af_chi2_fit0_deldop += err_fit0;
			c2af_dof_fit0_deldop += dof_fit0;
		}
	}
}
__global__ void c2af_add_doppler_contributions_krnl(
		struct par_t *dpar,
		struct dopfrm_t **frame,
		double *chi2_frame,
		float *o2_arr,
		float *m2_arr,
		float *om_arr,
		float4 *o2_m2_om_wt,
		int *ndop,
		int s,
		int nframes) {
	/* nframes-threaded kernel */
	int j, frm = threadIdx.x;
	double calval, err, o2_fit0, dof_fit0, err_fit0, thresh_fit0;
	float o2, m2, om, wt;

	if (frm < nframes) {
		o2 = o2_m2_om_wt[frm].w + o2_arr[frm];
		m2 = o2_m2_om_wt[frm].x + m2_arr[frm];
		om = o2_m2_om_wt[frm].y + om_arr[frm];
		wt = o2_m2_om_wt[frm].z;
		chi2_frame[frm] = 0.0;

		/* If this frame's calibration factor is allowed to float, set it to
		 * minimize chi-square, the sum over all bins of
		        { (obs - calfact*fit)^2 / variance }.                     */
		if (frame[frm]->cal.state == 'f') {
			if (om > 0.0) 	frame[frm]->cal.val = om/m2;
			else {
				frame[frm]->cal.val = TINYCALFACT;
				if (dpar->action != FIT )
					printf("WARNING: set %2d frame %2d had negative calfact reset to %10.4e\n",
							s, frm, frame[frm]->cal.val);	}
		}

		/* Compute chi-square for this frame  */
		calval = frame[frm]->cal.val;
		err = wt*(o2 - 2*calval*om + calval*calval*m2);
		frame[frm]->chi2 = err;
		chi2_frame[frm] += err;

		/* Compute chi-square contributions and dof due to bins whose model
		 * signal is =< 'chi2fit0_thresh' standard deviations of the noise in
		 * the data frame   */
		o2_fit0 = 0.0;
		dof_fit0 = 0.0;
		err_fit0 = 0.0;
		thresh_fit0 = dpar->chi2fit0_thresh * frame[frm]->sdev;
		if (dpar->write_chi2fit0) {
			for (j=0; j<ndop[0]; j++)
				if (calval*frame[frm]->fit[j] <= thresh_fit0) {
					o2_fit0 += frame[frm]->obs[j]*frame[frm]->obs[j]
					  *frame[frm]->oneovervar[j];
					if (frame[frm]->oneovervar[j] > 0.0)
						dof_fit0 += wt;
				}
			err_fit0 = wt*o2_fit0;
			c2af_chi2_fit0_doppler += err_fit0;
			c2af_dof_fit0_doppler += dof_fit0;
		}
	}
}
__global__ void c2af_doppler_add_o2_krnl(
		struct par_t *dpar,
		struct dopfrm_t **frame,
		int *ndop,
		float *o2,
		float *m2,
		float *om,
		int s,
		int frmsz,
		int nframes) {

	/* ndop-threaded kernel */
	int total_offset = blockIdx.x * blockDim.x + threadIdx.x;
	int frm = total_offset / frmsz;
	int j = total_offset % frmsz;
	float temp1, temp2, temp3;

	if ((j < ndop[0]) && (frm<nframes)) {
		/* Add contributions from power within frame limits */

		/* Next 2 lines implement: o2 += obs[j]*obs[j]*oneovervar[j];	 */
		temp1 = frame[frm]->obs[j] * frame[frm]->obs[j] *
				frame[frm]->oneovervar[j];
		atomicAdd(&o2[frm], temp1);

		/* Next 2 lines implement: m2 += fit[j]*fit[j]*oneovervar[j];		 */
		temp2 = frame[frm]->fit_s[j] * frame[frm]->fit_s[j] *
				frame[frm]->oneovervar[j];
		atomicAdd(&m2[frm], temp2);

		/* Next 2 lines implement: om += fit[j]*obs[j]*oneovervar[j];		 */
		temp3 = frame[frm]->fit_s[j] * frame[frm]->obs[j] *
				frame[frm]->oneovervar[j];
		atomicAdd(&om[frm], temp3);
	}
}
__global__ void c2af_set_chi2_krnl(struct dat_t *ddat, double chi2, int s) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0)
		ddat->set[s].chi2 = chi2;
}
__global__ void c2af_get_prntflgs_krnl(struct par_t *dpar, struct dat_t *ddat) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		c2af_print_breakdown = (ddat->dof_deldop > SMALLVAL || ddat->dof_doppler  > SMALLVAL
				|| ddat->dof_poset    > SMALLVAL	|| ddat->dof_lghtcrv  > SMALLVAL);
		c2af_dof_deldop = ddat->dof_deldop;
		c2af_dof_doppler = ddat->dof_doppler;
		c2af_dof_poset = ddat->dof_poset;
		c2af_dof_lghtcrv = ddat->dof_lghtcrv;
		c2af_dof = ddat->dof;
		c2af_write_chi2fit0 = dpar->write_chi2fit0;
		c2af_badradar = dpar->badradar;
		c2af_badphoto = dpar->badphoto;
		c2af_baddopscale = dpar->baddopscale;
		c2af_badposet = dpar->badposet;
		c2af_posbnd = dpar->posbnd;
		c2af_baddiam = dpar->baddiam;
	}
}

__host__ double chi2_cuda_af(struct par_t *dpar, struct dat_t *ddat,
		int list_breakdown, int nsets)
{
	int s, print_breakdown;
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
	chi2 = 0.0;
	/* Initialize variables that accumulate chi-square values  */
	c2af_init_krnl<<<1,1>>>(ddat);
	checkErrorAfterKernelLaunch("c2_init_krnl, chi2_cuda");

	/* Loop through all datasets, carry out chi-square computations, and
	 * provide screen and image output                            */
	for (s=0; s<nsets; s++) {
		/* Launch single-threaded kernel to get type */
		c2af_get_type_krnl<<<1,1>>>(ddat, s);
		checkErrorAfterKernelLaunch("c2af_get_type_krnl");
		gpuErrchk(cudaMemcpyFromSymbol(&type, c2af_type, sizeof(unsigned char),
				0, cudaMemcpyDeviceToHost));

		switch (type) {
		case DELAY:
			chi2 = chi2_deldop_cuda_af(dpar, ddat, s, list_breakdown,
				   &chi2_all_deldop, &chi2_fit0_deldop, &dof_fit0_deldop);
			c2af_set_chi2_krnl<<<1,1>>>(ddat, chi2, s);
			checkErrorAfterKernelLaunch("c2af_set_chi2_krnl");
			break;
		case DOPPLER:
			chi2 = chi2_doppler_cuda_af(dpar, ddat, s, list_breakdown,
				   &chi2_all_doppler, &chi2_fit0_deldop, &dof_fit0_doppler);
			c2af_set_chi2_krnl<<<1,1>>>(ddat, chi2, s);
			checkErrorAfterKernelLaunch("c2af_set_chi2_krnl");
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
		c2af_add_chi2_krnl<<<1,1>>>(ddat, s);
		checkErrorAfterKernelLaunch("c2_add_chi2_krnl, chi2_cuda");

	}  /* end for loop over datasets */
	/* Launch single-threaded kernel to retrieve ddat->chi2 to return it */
	c2af_retrieve_chi2_krnl<<<1,1>>>(ddat);
	checkErrorAfterKernelLaunch("c2af_retrieve_chi2_krnl");
	gpuErrchk(cudaMemcpyFromSymbol(&chi2, c2af_chi2, sizeof(double),
			0, cudaMemcpyDeviceToHost));
	/*.......................................................................*/


	/* Call kernel to get flags from ddat */
	c2af_get_prntflgs_krnl<<<1,1>>>(dpar, ddat);
	checkErrorAfterKernelLaunch("c2_get_prntflgs_krnl");
	gpuErrchk(cudaMemcpyFromSymbol(&print_breakdown, c2af_print_breakdown,
			sizeof(int), 0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&dof_deldop, c2af_dof_deldop, sizeof(double),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&dof_doppler, c2af_dof_doppler, sizeof(double),
			0, cudaMemcpyDeviceToHost));
//	gpuErrchk(cudaMemcpyFromSymbol(&dof_poset, c2_dof_poset, sizeof(double),
//			0, cudaMemcpyDeviceToHost));
//	gpuErrchk(cudaMemcpyFromSymbol(&dof_lgthcrv, c2_dof_lgthcrvr, sizeof(double),
//			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&write_chi2fit0, c2af_write_chi2fit0,
			sizeof(unsigned char), 0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&dof, c2af_dof, sizeof(double),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&baddiam, c2af_baddiam, sizeof(unsigned char),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&badphoto, c2af_badphoto, sizeof(unsigned char),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&posbnd, c2af_posbnd, sizeof(unsigned char),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&badposet, c2af_badposet, sizeof(unsigned char),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&badradar, c2af_badradar, sizeof(unsigned char),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&baddopscale, c2af_baddopscale, sizeof(unsigned char),
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
	int *ndel, *ndop, nframes, nThreads, frmsz;
	double chi2_set, *chi2_frame;
	float4 *o2_m2_om_wt;
	float **temp_o2, **temp_m2, **temp_om;
	struct deldopfrm_t **frame;
	dim3 BLK,THD;
	chi2_set = 0.0;

	/* Launch nframes-threaded kernel to get # of frames for this set */
	c2af_get_frames_krnl<<<1,1>>>(ddat, s);
	checkErrorAfterKernelLaunch("c2af_deldop_get_frames_krnl");
	gpuErrchk(cudaMemcpyFromSymbol(&nframes, c2af_nframes, sizeof(int),
			0, cudaMemcpyDeviceToHost));

	/* Allocate memory */
	cudaCalloc1((void**)&ndel, 		 sizeof(int),    nframes);
	cudaCalloc1((void**)&ndop, 		 sizeof(int), 	 nframes);
	cudaCalloc1((void**)&o2_m2_om_wt, sizeof(float4), nframes);
	cudaCalloc1((void**)&frame, 		 sizeof(struct deldopfrm_t), nframes);
	cudaCalloc1((void**)&temp_o2, 	 sizeof(float*), nframes);
	cudaCalloc1((void**)&temp_m2,   	 sizeof(float*), nframes);
	cudaCalloc1((void**)&temp_om, 	 sizeof(float*), nframes);
	cudaCalloc1((void**)&chi2_frame,  sizeof(double), nframes);

	/* Launch nframes-threaded kernel to get ndel/ndop for this frame and
	 * set up other radar parameter shortcuts. See kernel for details. */
	THD.x = nframes;
	c2af_deldop_set_shortcuts_krnl<<<1,THD>>>(ddat, frame, ndel, ndop, \
			o2_m2_om_wt, s, nframes);
	checkErrorAfterKernelLaunch("c2af_deldop_set_shortcuts_krnl");
	deviceSyncAfterKernelLaunch("c2af_deldop_set_shortcuts_krnl");

	/* Allocate the inner loops for temp arrays */
	frmsz = ndel[0] * ndop[0];
	nThreads = nframes * frmsz;
	for (int frm=0; frm<nframes; frm++) {
		cudaCalloc1((void**)&temp_o2[frm], sizeof(float), frmsz);
		cudaCalloc1((void**)&temp_m2[frm], sizeof(float), frmsz);
		cudaCalloc1((void**)&temp_om[frm], sizeof(float), frmsz);
	}

	/* Configure and launch kernel that prepares three arrays (o2, m2, om)
	 * for parallel reduction 	 */
	BLK.x = floor((maxThreadsPerBlock-1+nThreads)/maxThreadsPerBlock);
	THD.x = maxThreadsPerBlock;

	c2af_deldop_prep_o2_m2_krnl<<<BLK,THD>>>(ddat,frame,temp_o2,temp_m2,temp_om,
			ndel,ndop,s,nframes,frmsz);
	checkErrorAfterKernelLaunch("c2af_deldop_prep_o2_m2_krnl");

	/* Now launch the parallel reduction function to get the summations for
	 * all three double-arrays (all frames). Results for each frame and each
	 * array is in temp_o2[frm][0].	 */
	c2af_deldop_add_o2_m2(temp_o2, temp_m2, temp_om, frmsz, nframes);

	THD.x = nframes;
	c2af_add_deldop_contributions_krnl<<<1,THD>>>(dpar, frame, temp_o2, temp_m2,
			temp_om, chi2_frame, o2_m2_om_wt, ndel, ndop, s, nframes);
	checkErrorAfterKernelLaunch("c2af_add_deldop_contributions_krnl");
	deviceSyncAfterKernelLaunch("c2af_add_deldop_contributions_krnl");
	gpuErrchk(cudaMemcpyFromSymbol(&chi2_fit0_deldop, c2af_chi2_fit0_deldop,
			sizeof(double),	0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&dof_fit0_deldop, c2af_dof_fit0_deldop,
				sizeof(double),	0, cudaMemcpyDeviceToHost));

	for (int frm=0; frm<nframes; frm++) {
		chi2_set += chi2_frame[frm];
		if (list_breakdown)
			*chi2_all_deldop += chi2_frame[frm];
	}

	/* De-allocate memory */
	cudaFree(ndel);
	cudaFree(ndop);
	cudaFree(o2_m2_om_wt);
	cudaFree(frame);
	cudaFree(temp_o2);
	cudaFree(temp_m2);
	cudaFree(temp_om);
	cudaFree(chi2_frame);
	return chi2_set;
}

__host__ double chi2_doppler_cuda_af(struct par_t *dpar, struct dat_t *ddat, int s,
		int list_breakdown, double *chi2_all_doppler, double *chi2_fit0_doppler,
		double *dof_fit0_doppler)
{
	int *ndop, nframes, frmsz, nThreads;
	double chi2_set, *chi2_frame;
	struct dopfrm_t **frame;
	float *temp_o2, *temp_m2, *temp_om;
	float4 *o2_m2_om_wt;
	dim3 BLK,THD;

	/* Initialize chi-square for dataset  */
	chi2_set = 0.0;

	/* Launch nframes-threaded kernel to get # of frames for this set */
	c2af_get_frames_krnl<<<1,1>>>(ddat, s);
	checkErrorAfterKernelLaunch("c2af_get_frames_krnl");
	gpuErrchk(cudaMemcpyFromSymbol(&nframes, c2af_nframes, sizeof(int),
			0, cudaMemcpyDeviceToHost));

	/* Allocate memory */
	cudaCalloc1((void**)&ndop, 		 sizeof(int), 	 nframes);
	cudaCalloc1((void**)&o2_m2_om_wt, sizeof(float4), nframes);
	cudaCalloc1((void**)&frame, 		 sizeof(struct dopfrm_t), nframes);
	cudaCalloc1((void**)&temp_o2, 	 sizeof(float), nframes);
	cudaCalloc1((void**)&temp_m2,   	 sizeof(float), nframes);
	cudaCalloc1((void**)&temp_om, 	 sizeof(float), nframes);
	cudaCalloc1((void**)&chi2_frame,  sizeof(double), nframes);

	/* Launch nframes-threaded kernel to get ndel/ndop for this frame and
	 * set up other radar parameter shortcuts. See kernel for details. */
	THD.x = nframes;
	c2af_doppler_set_shortcuts_krnl<<<1,THD>>>(ddat, frame, ndop, o2_m2_om_wt,
			s, nframes);
	checkErrorAfterKernelLaunch("c2af_doppler_set_shortcuts_krnl");
	deviceSyncAfterKernelLaunch("c2af_doppler_set_shortcuts_krnl");

	frmsz = ndop[0];
	nThreads = nframes * frmsz;

	/* Add contributions from power within limits of data frame.
	 * This kernel also takes care of the frame's calibration factor and
	 * computes chi2 for this frame */
	BLK.x = floor((maxThreadsPerBlock-1+nThreads)/maxThreadsPerBlock);
	THD.x = maxThreadsPerBlock;

	c2af_doppler_add_o2_krnl<<<BLK,THD>>>(dpar, frame, ndop, temp_o2,
			temp_m2, temp_om, s, frmsz, nframes);
	checkErrorAfterKernelLaunch("c2af_add_o2_krnl");

	THD.x = nframes;
	c2af_add_doppler_contributions_krnl<<<1,THD>>>(dpar, frame, chi2_frame,
			temp_o2, temp_m2, temp_om, o2_m2_om_wt, ndop, s, nframes);
	checkErrorAfterKernelLaunch("c2af_add_doppler_contributions_krnl, line ");
	deviceSyncAfterKernelLaunch("c2af_add_doppler_contributions_krnl");
	gpuErrchk(cudaMemcpyFromSymbol(&chi2_fit0_doppler, c2af_chi2_fit0_deldop,
			sizeof(double),	0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&dof_fit0_doppler, c2af_dof_fit0_deldop,
			sizeof(double),	0, cudaMemcpyDeviceToHost));

	for (int frm=0; frm<nframes; frm++) {
		chi2_set += chi2_frame[frm];
		if (list_breakdown)
			*chi2_all_doppler += chi2_frame[frm];
	}

	/* De-allocate memory */
	cudaFree(ndop);
	cudaFree(o2_m2_om_wt);
	cudaFree(frame);
	cudaFree(temp_o2);
	cudaFree(temp_m2);
	cudaFree(temp_om);
	cudaFree(chi2_frame);
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
