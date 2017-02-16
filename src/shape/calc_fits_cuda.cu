/*****************************************************************************************
                                                                              calc_fits.c

As the name implies, this routine calculates the fits to each data frame for the current
set of model parameters.  For example, for each delay-Doppler frame it calls routine
posvis to create the model plane-of-sky image and then routine pos2deldop to create the
model delay-Doppler image from this POS image.

calc_fits also performs some of the screen and file output required by the "write" action;
in particular, it carries out tasks that require information associated with plane-of-sky
renderings, since such information is quickly overwritten if the "pos_scope" parameter is
set to "global" (i.e., if all frames and lightcurve points share the same memory for their
"pos" structures).

Modified 2015 June 10 by CM:
    Implement smearing for the "fit" and "write" actions

Modified 2014 February 14 by CM:
    Add "ilaw" argument to the apply_photo routine

Modified 2013 July 28 by CM:
    For the "write" action, output ppm POS images when the "write_highlight" parameter is
        turned on

Modified 2013 July 7 by CM:
    For the "write" action for lightcurve points and plane-of-sky frames, display the
        body-fixed longitude and latitude of the phase-angle bisector

Modified 2013 June 25 by CM:
    Allow POS images written for optical data to be annotated with principal-axis shafts
        and the angular momentum vector
    For POS images (sky renderings), display the name of the image file and the maximum
        pixel value in the plot_surface routine (called by the write_pos routine) rather
        than here

Modified 2013 April 24 by CM:
    Implement the "listpos_deldop" "listpos_opt" and "listpos_path" parameters
    Adjust names of output images so they are in alphanumeric order if > 100 per dataset

Modified 2012 April 2 by CM:
    Correct instantaneous maximum breadth calculation for Doppler scaling factor

Modified 2011 August 14 by CM:
    Display sidereal spin vector at each epoch, even for a PA rotator, if
        any spin impulses are used

Modified 2010 September 1 by CM:
    Initialize variables to avoid compilation warnings

Modified 2010 July 29 by CM:
    Fix bug introduced in calc_lghtcrv: rotation phases weren't being
        displayed for the "write" action
    For the "write" action for lightcurve datasets, include shadowed
        regions in projected area (and geometric albedo calculation)
        and display percentage of projected area that's shadowed

Modified 2010 June 15 by CM:
    Revise arguments to pos2deldop and pos2doppler routines

Modified 2010 May 28 by CM:
    Fix bug introduced with preceding change: in calc_lghtcrv, only
        deallocate memory for the "write" action (since it wasn't
        allocated in the first place for other actions)

Modified 2010 May 24 by CM:
    For the "write" action for lightcurves, output the projected area and
        (for absolute photometry) geometric albedo

Modified 2010 April 12 by CM:
    For the "write" action, include overflow region when computing
        cross sections

Modified 2009 July 29 by CM:
    For the "write" action, fix bug: output ppm images rather than pgm
        images if the "plot_angmom" parameter is turned on
    For the "write" action, pass an argument to the "write_pos" routine
        explicitly telling it whether or not to produce a colored image

Modified 2009 April 3 by CM:
    Initialize the "posbnd_logfactor" parameter and later set it for
        models that extend beyond the POS frame
    Add "badposet" and "badposet_logfactor" parameters: initialize them
        here and then use the new "checkposet" routine to adjust them for
        plane-of-sky fit images that are too small to "contain" the
        target
    Add "badradar" and "badradar_logfactor" parameters: initialize them
        here and then use the "pos2deldop" and "pos2doppler" routines
        (which are now int rather than void) to adjust them for models that
        are too wide in delay-Doppler space for the routines to handle
    Add "warn_badradar" argument to pos2deldop and pos2doppler routines
    For the "write" action, display each plane-of-sky fit frame's linear
        dimensions, the linear dimensions of the rectangular subset that
        contains the target, and the linear COM offsets

Modified 2008 December 12 by CM:
    For the "write" action for NPA rotators, list Euler angles (giving
        the body-fixed axes' orientations in ecliptic coordinates) and
        spin vector components (in body-fixed coordinates) for each
        observation epoch
    For the "write" action for NPA rotators, ensure that maximum breadth
        is nonnegative

Modified 2007 August 10 by CM:
    Eliminated unused variables and cleaned up a printf format
    For POS model frames (sky renderings) associated with lightcurve points
        and with plane-of-sky data frames, don't display the maximum pixel
        value unless the "optposmax" parameter is nonzero

Modified 2007 August 4 by CM:
    Add comp matrix for POS frames
    Add orbit_offset and body arguments to posvis routine and remove
        facet argument
    Add orbit_xoff, orbit_yoff, orbit_dopoff and body parameters to
        pos2deldop and pos2doppler routines
    Add body argument to apply_photo routine

Modified 2007 January 17 by CM:
    For the "write" action, display instantaneous folded zero-crossing
        bandwidth for Doppler and delay-Doppler frames

Modified 2007 January 11 by CM:
    In calc_lghtcrv for the "write" action, count lightcurve points
        from 0 rather than 1, as is already done for lightcurve POS images
        (and for Doppler, delay-Doppler, and plane-of-sky frames)

Modified 2007 January 6 by CM:
    In calc_lghtcrv for the "write" action, save rotation phase for each
        calculated lightcurve point so they can be output by routine chi2,
        and use cubic spline interpolation to obtain rotation phase at
        each observation epoch.  Also display range of rotation phases
        if only one calculated point per lightcurve is displayed in full

Modified 2006 October 1 by CM:
    In calc_lghtcrv, model lightcurve points are now intensities
        (relative to the solar intensity) rather than magnitudes
    In calc_lghtcrv and calc_poset, apply_photo routine has been revised
        to account for the POS pixel area and the 1 AU Sun-target distance

Modified 2006 September 1 by CM and MCN:
    When "exclude_seen" parameter is used, add check that facet number
        pos->f[i][j] is nonnegative
    For the "write" action, don't display cross sections and albedos
        for uncalibrated (delay-)Doppler frames

Modified 2006 June 21 by CM:
    In calc_deldop, changed delres to del_per_pixel and dopres to
        dop_per_pixel
    In calc_doppler, changed dopres to dop_per_bin
    For POS renderings and plane-of-sky fit frames, changed res to
        km_per_pixel

Modified 2006 June 18 by CM:
    Allow each delay-Doppler frame within a dataset to have different
        dimensions after vignetting
    Allow each Doppler frame within a dataset to have different
        dimensions after vignetting
    Allow plane-of-sky frames to be rectangular rather than square,
        and no longer require an odd number of pixels per side
    Eliminate range datasets

Modified 2006 March 10 by CM:
    Add "speckle" argument to pos2deldop and pos2doppler routines

Modified 2005 October 6 by CM:
    For lightcurve datasets, replace SUNMAG constant by "sun_appmag"
        parameter, so that absolute photometry with filters other than
        V band can be used

Modified 2005 July 25 by CM:
    For "write" action, display the model radar cross section and albedo
        for each delay-Doppler and Doppler frame

Modified 2005 July 22 by CM:
    Created five separate routines for writing POS frames as images
        so that they can be called separately if the "mark_unseen"
        parameter is turned on for the "write" action (since in this
        case we must first process all datasets to see which model
        facets were "seen" and only then can write the POS images)

Modified 2005 July 14 by CM:
    Fix bug in computing LE-to-COM delay and distance, LE-to-TE
        delay and distance, and instantantaneous bandwidth and breadth

Modified 2005 July 13 by CM:
    For "write" action for lightcurve points and plane-of-sky frames,
        display the body-fixed longitude and latitude of the
        Sun-to-asteroid line

Modified 2005 July 5 by CM:
    Remove the "dir" argument from pos2deldop and pos2doppler and add
        the "set" argument

Modified 2005 July 3 by CM:
    For "write" action for lightcurve datasets, implement the
        "lcrv_writeall" parameter, which produces screen display for
        every model lightcurve point rather than just the one point
        which falls closest to the midpoint of the observations.

Modified 2005 June 25 by CM:
    For "write" action for delay-Doppler frames, display the delay and
        distance between the leading edge and the center of mass and
        between the leading edge and the trailing edge;
        for delay-Doppler and Doppler frames, display the instantaneous
        zero-crossing bandwidth and maximum breadth.  All of the above
        are obtained from the model's delay-Doppler limits as
        determined PRIOR to convolution with the delay and Doppler
        response functions.

Modified 2005 June 22 by CM:
    Keep track of which model facets have been "seen" (i.e., are visible
        from Earth, are unshadowed, and have sufficiently low scattering
        and incidence angles) in at least one data frame or lightcurve
        point

Modified 2005 April 23 by CM:
    For the "write" action, list whether or not epochs have been corrected
        for one-way light travel time

Modified 2005 March 1 by CM:
    Adjust arguments to the revised "resampim" routine to permit rotation
        of resampled plane-of-sky frames
    Initialize the "posbnd" parameter (flag indicating that the model
        extends beyond the model POS frame) to 0 here rather than in
        bestfit.c so that it can used for actions other than "fit"
    Fix bug in calc_poset which was incorrectly flagging the model as
        being too small for the model POS frame

Modified 2005 February 21 by CM:
    Use the new "poset_resample" parameter to allow interpolation methods
        other than bilinear for constructing plane-of-sky fit images for
        plane-of-sky data frames
    Add the new "image_rebin" argument to function resampim to handle
        plane-of-sky fit frames which have much coarser resolution
        than the model POS frames from which they are constructed
        (i.e., which are greatly undersampled)
    For "write" action, display maximum pixel value for model POS images
        for plane-of-sky frames and calculated lightcurve images
        (in case someone wants to use the "optposmax" parameter to
        truncate the image brightness)

Modified 2005 February 6 by CM:
    For "write" action, display rotation phase
    For "write" action, fix bug in computing the angular body-fixed
        coordinates of the line of sight for lightcurve datasets

Modified 2005 January 25 by CM:
    Take care of unused and uninitialized variables

Modified 2005 January 24 by CM:
    Add "calc_poset" routine to handle POS datasets
    For "write" action, display the angular body-fixed coordinates of
        the line of sight
    For "write" action, display calendar dates in addition to Julian dates
    For "write" action, display the date for range datasets

Modified 2004 December 19 by CM:
    For "write" action, display the projected area for each Doppler and
        delay-Doppler frame

Modified 2004 May 3 by CM:
    For "write" action, display the (delay-)Doppler corrections for each
        frame

Modified 2004 April 9 by CM:
    For "write" action, display the solar azimuth angles (N->E in the POS)

Modified 2004 March 27 by CM:
    Eliminate output of range (rng) plane-of-sky images for
        delay-Doppler frames
    For "write" action, display the epoch, solar phase angle and
        apparent spin vector direction at the midpoint of lightcurve
        datasets
    For "write" action, if "plot_spinvec" parameter is turned on, 
        POS pgm images include an arrow indicating the target's
        intrinsic spin vector.
    For "write" action, if "plot_subradar" parameter is turned on, 
        POS pgm images for (delay-)Doppler datasets include an X
        indicating the target's subradar point.
    For "write" action, if "plot_com" parameter is turned on, 
        POS pgm images for (delay-)Doppler datasets include a cross
        indicating the target's projected COM.
    For "write" action, if "plot_pa" parameter vector has any
        component(s) turned on, POS ppm images for (delay-)Doppler
        datasets include colored cylindrical shaft(s) indicating the
        positive end of the corresponding principal axis/axes.

Modified 2004 Feb 29 by CM:
    Add comments for lightcurves
    Remove "sdev" argument to routine gamma_trans
    Compute lightcurve magnitudes rather than negative magnitudes
    Eliminate the "curve_mm" lightcurve output file, since it nearly
        duplicates the "fit.mm" file (except that the cal factor
        isn't included)
    Move the lightcurve calculations to the new "calc_lghtcrv" routine
    Eliminate the unused dat argument to calc_deldop, calc_doppler,
        and calc_range
    Eliminate "type" argument to the "apply_photo" routine, and
        add the "phase" (solar phase angle) argument
    Label lightcurve POS images as 0 through (ncalc-1) rather than
        1 through ncalc, similar to (delay-)Doppler pgm images

Modified 2003 July 30 by CM:
    Add three parameters for rotating/flipping output pgm files
        for delay-Doppler images (fit, data, residuals)

Modified 2003 May 16 by CM:
    Add listres parameter for producing output files containing
        residual matrices

Modified 2003 May 13 by CM:
    Don't resample and recenter residual pgm images if dd_scaling = none
    Correct a bug in normalizing file output for Doppler fits

Modified 2003 May 10 by CM:
    Add scalefitobs parameter so that user can choose whether to scale
        the data and fit pgm images separately (default), to the maximum
        value of the two taken together, to the maximum fit value, or to
        the maximum data value

Modified 2003 May 7 by CM:
    Add sinc2width argument to pos2deldop and pos2doppler

Modified 2003 April 29 by CM:
    Don't truncate residuals to integer values before making pgm images
    Add nsinc2 argument to pos2deldop and pos2doppler

Modified 2003 April 28 by CM:
    Display two angles for the spin vector, not just one

Modified 2003 April 24 by CM:
    Move "delcom" from delay-Doppler datasets to individual frames

Modified 2003 April 23 by CM:
    Removed "deldopoffs" call from calc_deldop and "dopoffs" call from
        calc_deldop, since these calls are now included in realize_delcor
 *****************************************************************************************/
extern "C" {
#include "head.h"
}

__host__ void calc_deldop_cuda(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, int s);
__host__ void calc_doppler_cuda(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, int s);
//__host__ void calc_poset_cuda( struct par_t *par, struct mod_t *mod, int s);
//__host__ void calc_lghtcrv_cuda(struct par_t *par, struct mod_t *mod, struct
//		lghtcrv_t *lghtcrv, int s);

__device__ int cf_nf, cf_nsets, cf_nframes, cf_nviews, cf_ndel, cf_ndop,
			   cf_v0_index, cf_pos_n, cf_exclude_seen, cf_xlim0, cf_xlim1,
			   cf_ylim0, cf_ylim1;
__device__ unsigned char cf_type;
__device__ float cf_overflow_o2_store, cf_overflow_m2_store, cf_overflow_xsec_store,
		cf_overflow_delmean_store, cf_overflow_dopmean_store;
__device__ struct deldop_t *cf_deldop;
__device__ struct deldopfrm_t *cf_del_frame;
__device__ struct deldopview_t *cf_del_view0;
__device__ struct pos_t *cf_pos;
__device__ struct doppler_t *cf_doppler;
__device__ struct dopfrm_t *cf_dop_frame;
__device__ struct dopview_t *cf_dop_view0;

__global__ void cf_init_devpar_krnl(struct par_t *dpar, struct mod_t
		*dmod, struct dat_t *ddat) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		dpar->posbnd = 0;
		dpar->badposet = 0;
		dpar->badradar = 0;
		dpar->posbnd_logfactor = 0.0;
		dpar->badposet_logfactor = 0.0;
		dpar->badradar_logfactor = 0.0;
		cf_nf = dmod->shape.comp[0].real.nf;
		cf_nsets = ddat->nsets;
	}
}
__global__ void cf_init_seen_flags_krnl(struct mod_t *dmod) {
	/* nf-threaded kernel */
	int f = blockIdx.x * blockDim.x + threadIdx.x;

	if (f < cf_nf)
		dmod->shape.comp[0].real.f[f].seen = 0;
}
__global__ void cf_get_set_type_krnl(struct dat_t *ddat, int s) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0)
		cf_type = ddat->set[s].type;
}
__global__ void cf_set_final_pars_krnl(struct par_t *dpar, struct
		dat_t *ddat) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		dpar->posbnd_logfactor /= ddat->dof;
		dpar->badposet_logfactor /= ddat->dof_poset;
		dpar->badradar_logfactor /= (ddat->dof_deldop + ddat->dof_doppler);
	}
}

__host__ void calc_fits_cuda(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat)
{
	int s, nf, nsets;
	unsigned char type;
	dim3 BLK,THD;

	/* Initialize flags that indicate the model extends beyond POS frame, that
	 * plane-of-sky fit images are too small to "contain" the target, and that
	 * model is too wide in (delay-)Doppler space to create (delay-)Doppler fit
	 * frames.  Note that this also gets mod->shape.nf and nsets            */

	cf_init_devpar_krnl<<<1,1>>>(dpar, dmod, ddat);
	checkErrorAfterKernelLaunch("cf_init_devpar_krnl (calc_fits_cuda)");
	gpuErrchk(cudaMemcpyFromSymbol(&nf, cf_nf, sizeof(int),
				0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&nsets, cf_nsets, sizeof(int),
					0, cudaMemcpyDeviceToHost));

	/* Initialize the flags that indicate whether or not each facet of each
	 * model component is ever visible and unshadowed from Earth
	 * Note:  Single component only for now.  */
	//for (c=0; c<mod->shape.ncomp; c++)
	BLK.x = floor((maxThreadsPerBlock - 1 + nf)/maxThreadsPerBlock);
	THD.x = maxThreadsPerBlock;
	cf_init_seen_flags_krnl<<<BLK,THD>>>(dmod);
	checkErrorAfterKernelLaunch("cf_init_seen_flags_krnl (calc_fits_cuda)");

	/* Calculate the fits for each dataset in turn  */
	for (s=0; s<nsets; s++) {

		/* Get data type */
		cf_get_set_type_krnl<<<1,1>>>(ddat, s);
		checkErrorAfterKernelLaunch("cf_init_seen_flags_krnl (calc_fits_cuda)");
		gpuErrchk(cudaMemcpyFromSymbol(&type, cf_type, sizeof(unsigned char),
				0, cudaMemcpyDeviceToHost));

		switch (type) {
		case DELAY:
			calc_deldop_cuda(dpar, dmod, ddat, s);
			break;
		case DOPPLER:
			calc_doppler_cuda(dpar, dmod, ddat, s);
			break;
		case POS:
			printf("Write calc_poset_cuda!");
//			calc_poset_cuda(dpar, dmod, s);
			break;
		case LGHTCRV:
			printf("Write calc_lghtcrv_cuda!");
//			calc_lghtcrv_cuda(dpar, dmod, s);
			break;
		default:
			printf("calc_fits_cuda.c: can't handle this type yet\n");
		}
	}

	/* Complete calculations of values that will be used during a fit to
	 * increase the objective function for models with bad properties   */
	cf_set_final_pars_krnl<<<1,1>>>(dpar, ddat);
	checkErrorAfterKernelLaunch("cf_set_final_pars_krnl (calc_fits_cuda)");

}

__global__ void cf_get_frames_krnl(struct dat_t *ddat, int s) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		switch(ddat->set[s].type) {
		case DELAY:
			cf_nframes = ddat->set[s].desc.deldop.nframes;
			break;
		case DOPPLER:
			cf_nframes = ddat->set[s].desc.doppler.nframes;
			break;
		case POS:
			cf_nframes = ddat->set[s].desc.poset.nframes;
			break;
		case LGHTCRV:
			cf_nframes = ddat->set[s].desc.lghtcrv.ncalc;
			break;
		}

	}
}
__global__ void cf_set_shortcuts_krnl(struct dat_t *ddat, int s, int f) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		switch(cf_type) {
		case DELAY:
			cf_deldop 	 = &ddat->set[s].desc.deldop;
			cf_del_frame = &cf_deldop->frame[f];
			cf_ndop 	 = cf_del_frame->ndop;
			cf_ndel 	 = cf_del_frame->ndel;
			cf_del_view0 = &cf_del_frame->view[cf_deldop->v0];
			cf_nviews	 = cf_deldop->nviews;
			cf_v0_index  = cf_deldop->v0;
			cf_pos		 = &cf_del_frame->pos;
			break;
		case DOPPLER:
			cf_doppler 	 = &ddat->set[s].desc.doppler;
			cf_dop_frame = &cf_doppler->frame[f];
			cf_dop_view0 = &cf_dop_frame->view[cf_doppler->v0];
			cf_nviews	 = cf_doppler->nviews;
			cf_v0_index  = cf_doppler->v0;
			cf_ndop		 = cf_dop_frame->ndop;
			cf_pos		 = &cf_dop_frame->pos;
			break;
		}
		cf_overflow_o2_store = 0.0;
		cf_overflow_m2_store = 0.0;
		cf_overflow_xsec_store = 0.0;
		cf_overflow_delmean_store = 0.0;
		cf_overflow_dopmean_store = 0.0;
	}
}
__global__ void cf_set_pos_ae_krnl(int v) {
	/* 9-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % 3;
	int j = offset / 3;

	if (offset < 9) {
		switch (cf_type) {
		case DELAY:
			cf_pos->ae[i][j] = cf_del_frame->view[v].ae[i][j];
			cf_pos->oe[i][j] = cf_del_frame->view[v].oe[i][j];
			break;
		case DOPPLER:
			cf_pos->ae[i][j] = cf_dop_frame->view[v].ae[i][j];
			cf_pos->oe[i][j] = cf_dop_frame->view[v].oe[i][j];
			break;
		}

		/* Single-thread task */
		if (offset == 0) {
			cf_pos->bistatic = 0;
			cf_pos_n = cf_pos->n;
		}
	}
}
__global__ void cf_posclr_krnl(int n, int nx)
{
	/* Multi-threaded kernel (2*pos->n + 1)^2 threads) */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = (offset % nx) - n;
	int j = (offset / nx) - n;

	if (offset < (2*cf_pos->n+1)*(2*cf_pos->n+1)) {
		/* For each POS pixel, zero out the optical brightness (b) and
		 * cos(scattering angle), reset the z coordinate (distance from COM towards
		 * Earth) to a dummy value, and reset the body, component, and facet onto
		 * which the pixel center projects to  dummy values                  */
		cf_pos->b[i][j] = cf_pos->cose[i][j] = 0.0;
		cf_pos->z[i][j] = -HUGENUMBER;
		cf_pos->body[i][j] = cf_pos->comp[i][j] = cf_pos->f[i][j] = -1;

		cf_pos->b_s[offset] = cf_pos->cose_s[offset] = 0.0;
		cf_pos->z_s[offset] = -HUGENUMBER;

		/* In the x direction, reset the model's leftmost and rightmost
		 * pixel number to dummy values, and similarly for the y direction   */
		cf_pos->xlim[0] = cf_pos->ylim[0] =  n;
		cf_pos->xlim[1] = cf_pos->ylim[1] = -n;

		/* For a bistatic situation (lightcurve or plane-of-sky dataset), zero out
		 * cos(incidence angle) and reset the distance towards the sun, the body,
		 * component, and facet numbers as viewed from the sun, and the model's
		 * maximum projected extent as viewed from the sun to dummy values    */
		if (cf_pos->bistatic) {
			cf_pos->cosill[i][j] = 0.0;
			cf_pos->zill[i][j] = -HUGENUMBER;
			cf_pos->bodyill[i][j] = cf_pos->compill[i][j] = cf_pos->fill[i][j] = -1;

			cf_pos->cosill_s[offset] = 0.0;
			cf_pos->zill_s[offset] = 0.0;

			cf_pos->xlim2[0] = cf_pos->ylim2[0] =  n;
			cf_pos->xlim2[1] = cf_pos->ylim2[1] = -n;
		}
	}
}
__global__ void cf_set_posbnd_krnl(struct par_t *dpar) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		dpar->posbnd = 1;
		switch (cf_type) {
		case DELAY:
			dpar->posbnd_logfactor += cf_del_frame->dof * cf_pos->posbnd_logfactor;
			break;
		case DOPPLER:
			dpar->posbnd_logfactor += cf_dop_frame->dof * cf_pos->posbnd_logfactor;
			break;
		}

	}
}
__global__ void cf_get_exclude_seen_krnl(struct par_t *dpar, struct
		mod_t *dmod) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		cf_exclude_seen = dpar->exclude_seen;
		cf_xlim0 = cf_pos->xlim[0];
		cf_xlim1 = cf_pos->xlim[1];
		cf_ylim0 = cf_pos->ylim[0];
		cf_ylim1 = cf_pos->ylim[1];
	}
}
__global__ void cf_mark_pixels_seen_krnl(struct par_t *dpar,
		struct mod_t *dmod, int npixels, int xspan) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int k = (offset % xspan) + cf_xlim0;
	int l = (offset / xspan) + cf_ylim0;
	int facetnum;

	if (offset < npixels) {
		if ((cf_pos->cose_s[offset] > dpar->mincosine_seen)
				&& (cf_pos->f[k][l] >= 0)) {
			facetnum = cf_pos->f[k][l];
			//c = cf_pos->comp[k][l];
			dmod->shape.comp[0].real.f[facetnum].seen = 1;
		}
	}
}
__global__ void cf_set_badradar_krnl(struct par_t *dpar) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		switch (cf_type) {
		case DELAY:
			dpar->badradar = 1;
			dpar->badradar_logfactor += cf_del_frame->dof *
				cf_del_frame->badradar_logfactor / cf_deldop->nviews;
			break;
		case DOPPLER:
			dpar->badradar = 1;
			dpar->badradar_logfactor += cf_dop_frame->dof *
				cf_dop_frame->badradar_logfactor / cf_doppler->nviews;
			break;
		}
	}
}
__global__ void cf_add_fit_store_krnl1(struct dat_t *ddat, float *fit_store,
		int nThreads, int s, int f) {
	/* ndel*ndop-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;

	if (offset < (nThreads)) {
		switch (cf_type) {
		case DELAY:
			fit_store[offset] += ddat->set[s].desc.deldop.frame[f].fit_s[offset];
			break;
		case DOPPLER:
			fit_store[offset] += ddat->set[s].desc.doppler.frame[f].fit_s[offset];
			break;
		}
	}
}
__global__ void cf_add_fit_store_krnl2() {
	/* ndel*ndop-threaded kernel */
	if (threadIdx.x == 0) {
		switch (cf_type) {
		case DELAY:
			cf_overflow_o2_store += cf_del_frame->overflow_o2;
			cf_overflow_m2_store += cf_del_frame->overflow_m2;
			cf_overflow_xsec_store += cf_del_frame->overflow_xsec;
			cf_overflow_delmean_store += cf_del_frame->overflow_delmean;
			cf_overflow_dopmean_store += cf_del_frame->overflow_dopmean;
			break;
		case DOPPLER:
			cf_overflow_o2_store += cf_dop_frame->overflow_o2;
			cf_overflow_m2_store += cf_dop_frame->overflow_m2;
			cf_overflow_xsec_store += cf_dop_frame->overflow_xsec;
			cf_overflow_dopmean_store += cf_dop_frame->overflow_dopmean;
		}

	}
	__syncthreads();
}
__global__ void cf_finish_fit_store_krnl(struct dat_t *ddat, float *fit_store,
		int s, int f, int nThreads) {
	/* multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;

	if (offset < nThreads)
		switch (cf_type) {
		case DELAY:
			ddat->set[s].desc.deldop.frame[f].fit_s[offset] = fit_store[offset];
			break;
		case DOPPLER:
			ddat->set[s].desc.doppler.frame[f].fit_s[offset] = fit_store[offset];
			break;
		}
}
__global__ void cf_finish_fit_krnl2() {
	/* Single-threaded Kernel */
	int nviews;
	if (threadIdx.x == 0) {
		switch (cf_type) {
		case DELAY:
			nviews = cf_deldop->nviews;
			cf_del_frame->overflow_o2 = cf_overflow_o2_store / nviews;
			cf_del_frame->overflow_m2 = cf_overflow_m2_store / nviews;
			cf_del_frame->overflow_xsec = cf_overflow_xsec_store / nviews;
			cf_del_frame->overflow_delmean = cf_overflow_delmean_store / nviews;
			cf_del_frame->overflow_dopmean = cf_overflow_dopmean_store / nviews;
			break;
		case DOPPLER:
			nviews = cf_doppler->nviews;
			cf_dop_frame->overflow_o2 = cf_overflow_o2_store / nviews;
			cf_dop_frame->overflow_m2 = cf_overflow_m2_store / nviews;
			cf_dop_frame->overflow_xsec = cf_overflow_xsec_store / nviews;
			cf_dop_frame->overflow_dopmean = cf_overflow_dopmean_store / nviews;
		}

	}
}
__global__ void cf_gamma_trans_krnl(struct par_t *dpar, struct dat_t *ddat,
		int s, int f, int nThreads) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	if (offset < nThreads) {
		/*  Carry out a gamma transformation on the fit image if requested  */
		if (dpar->dd_gamma != 1.0) {
			switch (cf_type) {
			case DELAY:
				dev_gamma_trans(&ddat->set[s].desc.deldop.frame[f].fit_s[offset],
						dpar->dd_gamma);
				break;
			case DOPPLER:
				//cf_dop_frame->fit[offset] = fit[offset];
				break;
			}

		}
	}
}
__global__ void cf_copy_fit_back_krnl(struct dat_t *ddat, float *fit, int s,
		int f, int nThreads) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int idop;

	if (offset < nThreads) {
		switch (ddat->set[s].type) {
		case DELAY:
			int idel = offset % cf_ndel - 1;
			idop = offset / cf_ndel - 1;
			if (idel>=0 && idel <= cf_ndel && idop>=0 && idop<=cf_ndop) {
				//			ddat->set[s].desc.deldop.frame[f].fit[idel][idop] = fit[offset];
				ddat->set[s].desc.deldop.frame[f].fit_s[idop*cf_ndel+idel] = fit[offset];
				//fit_copy[offset] = fit[offset];
			}
			if (idop == cf_ndop-2)
				ddat->set[s].desc.deldop.frame[f].fit_s[(idop+1)*cf_ndel+idel] = 0.0;
			if (idel == cf_ndel-2)
				ddat->set[s].desc.deldop.frame[f].fit_s[idop*cf_ndel+idel+1] = 0.0;
			break;
		case DOPPLER:
			idop = offset;
			if (idop>=0 && idop<=cf_ndop)
				ddat->set[s].desc.doppler.frame[f].fit[idop] = fit[offset];
			break;
		case POS:
			printf("Write code for POS in vp_copy_fit_back_krnl in "
					"vary_params_cuda");
			break;
		case LGHTCRV:
			printf("Write code for LGTHCRV in vp_copy_fit_back_krnl in "
					"vary_params_cuda");
			break;
		}
	}
}

__host__ void calc_deldop_cuda(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, int s)
{
	double orbit_offset[3] = {0.0, 0.0, 0.0};
	int nframes, nviews, ndel, ndop, v0_index, pos_n, nx, exclude_seen, xlim0,
		xlim1, ylim0, ylim1, f, v, v2;
	float *fit_store;
	dim3 BLK,THD;

	/* Get # of frames for this deldop */
	cf_get_frames_krnl<<<1,1>>>(ddat, s);
	checkErrorAfterKernelLaunch("cf_get_nframes_krnl (calc_deldop_cuda)");
	gpuErrchk(cudaMemcpyFromSymbol(&nframes, cf_nframes, sizeof(int),
			0, cudaMemcpyDeviceToHost));

	for (f=0; f<nframes; f++) {
		/* Set deldop, frame, view0, and pos */
		cf_set_shortcuts_krnl<<<1,1>>>(ddat, s, f);
		checkErrorAfterKernelLaunch("cf_deldop_set_shortcuts_krnl (calc_deldop_cuda)");
		gpuErrchk(cudaMemcpyFromSymbol(&nviews, cf_nviews, sizeof(int),
				0, cudaMemcpyDeviceToHost));
		gpuErrchk(cudaMemcpyFromSymbol(&ndel, cf_ndel, sizeof(int),
				0, cudaMemcpyDeviceToHost));
		gpuErrchk(cudaMemcpyFromSymbol(&ndop, cf_ndop, sizeof(int),
				0, cudaMemcpyDeviceToHost));
		gpuErrchk(cudaMemcpyFromSymbol(&v0_index, cf_v0_index, sizeof(int),
				0, cudaMemcpyDeviceToHost));

		/* If smearing is being modeled, initialize variables that
		 * will be used to sum results calculated for individual views.  */
		if (nviews > 1) {
			/* Allocate fit_store as a single pointer, originally a double
			 * pointer. This also initializes the entire array to zero. */
			cudaCalloc((void**)&fit_store, sizeof(float), ndel*ndop);
		}

		/*  Loop over all views for this (smeared) frame, going in an order that
        ends with the view corresponding to the epoch listed for this frame
        in the obs file; this way we can use the calculated information for
        that view in the "write" action screen and disk output that follows   */

		for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
			v = v2 % nviews;

			/* Launch 9-threaded kernel to set pos->ae,pos->oe,pos->bistatic.*/
			THD.x = 9;
			cf_set_pos_ae_krnl<<<BLK,THD>>>(v);
			checkErrorAfterKernelLaunch("cf_deldop_set_pos_ae_krnl "
					"(calc_deldop_cuda)");
			gpuErrchk(cudaMemcpyFromSymbol(&pos_n, cf_pos_n, sizeof(int),
					0, cudaMemcpyDeviceToHost));

			/* Configure & launch posclr_krnl to initialize POS view */
			BLK.x = floor((maxThreadsPerBlock - 1 + (2*pos_n+1)*(2*pos_n+1)) /
					maxThreadsPerBlock);
			THD.x = maxThreadsPerBlock; // Thread block dimensions
			nx = 2*pos_n + 1;
			cf_posclr_krnl<<<BLK,THD>>>(pos_n, nx);
			checkErrorAfterKernelLaunch("cf_posclr_krnl (calc_fits_cuda)");

			/* Call posvis_cuda_2 to get facet number, scattering angle,
			 * distance toward Earth at center of each POS pixel; set flag
			 * posbnd if any model portion extends beyond POS frame limits.*/
			/* NOTE: Limited to single component for now */

			if (posvis_cuda_2(dpar, dmod, ddat, orbit_offset, s, f, 0, 0, 0) &&
					v == v0_index) {
				/* Call single-threaded kernel to set dpar->posbnd and
				 * dpar->posbnd_logfactor */
				cf_set_posbnd_krnl<<<1,1>>>(dpar);
				checkErrorAfterKernelLaunch("cf_deldop_set_posbnd_krnl (calc_fits_cuda)");
			}

			/* Launch single-threaded kernel to get dpar->exclude_seen */
			cf_get_exclude_seen_krnl<<<1,1>>>(dpar, dmod);
			checkErrorAfterKernelLaunch("cf_get_exclude_seen_krnl (calc_fits_cuda)");
			gpuErrchk(cudaMemcpyFromSymbol(&exclude_seen, cf_exclude_seen, sizeof(int),
					0, cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpyFromSymbol(&xlim0, cf_xlim0, sizeof(int),
					0, cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpyFromSymbol(&xlim1, cf_xlim1, sizeof(int),
					0, cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpyFromSymbol(&ylim0, cf_ylim0, sizeof(int),
					0, cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpyFromSymbol(&ylim1, cf_ylim1, sizeof(int),
					0, cudaMemcpyDeviceToHost));

			/* Go through all POS pixels which are visible with low enough
			 * scattering angle and mark the facets which project onto their
			 * centers as having been "seen" at least once                   */
			/* I'll launch a multi-threaded kernel here:
			 * (xlim1 - xlim0 + 1)^2 threads			 */

			if (s != exclude_seen && v == v0_index) {

				int xspan = xlim1 - xlim0 + 1;
				int yspan = ylim1 - ylim0 + 1;
				int nThreads = xspan * yspan;

				/* Configure & launch posclr_krnl to initialize POS view */
				BLK.x = floor((maxThreadsPerBlock - 1 + nThreads) /
						maxThreadsPerBlock);
				THD.x = maxThreadsPerBlock; // Thread block dimensions
				cf_mark_pixels_seen_krnl<<<BLK,THD>>>(dpar, dmod,
						nThreads, xspan);
				checkErrorAfterKernelLaunch("cf_posclr_krnl (calc_fits_cuda)");
			}

			/* Zero out the fit delay-Doppler image, then call pos2deldop to
			 * create the fit image by mapping power from the plane of the sky
			 * to delay-Doppler space.                             */
			BLK.x = floor((maxThreadsPerBlock - 1 + ndel*ndop) /
					maxThreadsPerBlock);
			THD.x = maxThreadsPerBlock; // Thread block dimensions
			clrvect_krnl<<<BLK,THD>>>(ddat, s, f, ndel*ndop);
			checkErrorAfterKernelLaunch("clrvect_krnl, calc_fits_cuda:882");

			if (pos2deldop_cuda_2(dpar,dmod,ddat,0.0,0.0,0.0,0, s,f,v)) {
				/* Call single-threaded kernel to set badradar flag and
				 * associated badradar_logfactor			 */
				cf_set_badradar_krnl<<<1,1>>>(dpar);
				checkErrorAfterKernelLaunch("cf_deldop_set_badradar_krnl (calc_fits_cuda)");
			}

			/* If smearing is being modeled, include delay-Doppler calculations
			 * from this view in the summed results for this frame  */
			if (nviews > 1) {
				/* Launch ndel*ndop-threaded kernel to add fit[i][j] to
				 * fit_store[i][j]*/
				BLK.x = floor((maxThreadsPerBlock - 1 + (ndel)*(ndop)) /
						maxThreadsPerBlock);
				THD.x = maxThreadsPerBlock; // Thread block dimensions
				cf_add_fit_store_krnl1<<<BLK,THD>>>(ddat,fit_store,(ndel*ndop),s,f);
				checkErrorAfterKernelLaunch("cf_add_fit_store_krnl1 (calc_fits_cuda)");
				cf_add_fit_store_krnl2<<<1,1>>>();
				checkErrorAfterKernelLaunch("cf_add_fit_store_krnl2 (calc_fits_cuda)");
			}
		}

		/* If smearing is being modeled, compute mean values over all views for
		 * this frame and store them in the standard frame structure     */
		/* This kernel also carries out the gamma transformation on the fit
		 * image if the par->dd_gamma flag is not set  */
		if (nviews > 1) {
			/* Launch ndel*ndop-threaded kernel to add fit[i][j] to
			 * fit_store[i][j]*/
			BLK.x = floor((maxThreadsPerBlock - 1 + (ndel)*(ndop)) /
					maxThreadsPerBlock);
			THD.x = maxThreadsPerBlock; // Thread block dimensions
			cf_finish_fit_store_krnl<<<BLK,THD>>>(ddat, fit_store,s, f, (ndel*ndop));
			checkErrorAfterKernelLaunch("cf_deldop_finish_fit_store (calc_fits_cuda)");
			cf_finish_fit_krnl2<<<1,1>>>();
			checkErrorAfterKernelLaunch("cf_finish_fit_krnl2");
			cf_gamma_trans_krnl<<<BLK,THD>>>(dpar, ddat, s, f, (ndel*ndop));
			checkErrorAfterKernelLaunch("cf_gamma_trans_krnl");
			cudaFree(fit_store);
		}
	}  /* end loop over frames */
}


__host__ void calc_doppler_cuda(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, int s)
{
	double orbit_offset[3] = {0.0, 0.0, 0.0};
	float *fit, *fit_store;
	int ndop, v0_index, pos_n, xlim0, xlim1, ylim0, ylim1, exclude_seen,
		nviews, nframes, nx, f, v, v2;
	dim3 BLK,THD;

	/* Get # of frames for this deldop */
	cf_get_frames_krnl<<<1,1>>>(ddat, s);
	checkErrorAfterKernelLaunch("cf_get_nframes_krnl (calc_deldop_cuda)");
	gpuErrchk(cudaMemcpyFromSymbol(&nframes, cf_nframes, sizeof(int),
			0, cudaMemcpyDeviceToHost));

	for (f=0; f<nframes; f++) {
		/* Set deldop, frame, view0, and pos */
		cf_set_shortcuts_krnl<<<1,1>>>(ddat, s, f);
		checkErrorAfterKernelLaunch("cf_deldop_1st_krnl (calc_deldop_cuda)");
		gpuErrchk(cudaMemcpyFromSymbol(&nviews, cf_nviews, sizeof(int),
				0, cudaMemcpyDeviceToHost));
		gpuErrchk(cudaMemcpyFromSymbol(&ndop, cf_ndop, sizeof(int),
				0, cudaMemcpyDeviceToHost));
		gpuErrchk(cudaMemcpyFromSymbol(&v0_index, cf_v0_index, sizeof(int),
				0, cudaMemcpyDeviceToHost));

		/* First allocate the fit array */
		cudaCalloc((void**)&fit, sizeof(float), ndop);

		/* If smearing is being modeled, initialize variables that
		 * will be used to sum results calculated for individual views.  */
		if (nviews > 1) {
			/* Allocate fit_store as a single pointer, originally a double
			 * pointer. This also initializes the entire array to zero. */
			cudaCalloc((void**)&fit_store, sizeof(float), ndop);
		}

		/* Loop over all views for this (smeared) frame, going in an order that
		 * ends with the view corresponding to the epoch listed for this frame
		 * in the obs file; this way we can use the calculated information for
		 * that view in the "write" action screen and disk output that follows*/

		for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
			v = v2 % nviews;

			/* Launch 9-threaded kernel to set pos->ae,pos->oe,pos->bistatic.*/
			THD.x = 9;
			cf_set_pos_ae_krnl<<<BLK,THD>>>(v);
			checkErrorAfterKernelLaunch("cf_deldop_set_pos_ae_krnl "
					"(calc_doppler_cuda)");
			gpuErrchk(cudaMemcpyFromSymbol(&pos_n, cf_pos_n, sizeof(int),
					0, cudaMemcpyDeviceToHost));

			/* Configure & launch posclr_krnl to initialize POS view */
			BLK.x = floor((maxThreadsPerBlock - 1 + (2*pos_n+1)*(2*pos_n+1)) /
					maxThreadsPerBlock);
			THD.x = maxThreadsPerBlock; // Thread block dimensions
			nx = 2*pos_n + 1;
			cf_posclr_krnl<<<BLK,THD>>>(pos_n, nx);
			checkErrorAfterKernelLaunch("cf_posclr_krnl (calc_fits_doppler)");

			/* Call posvis_cuda_2 to get facet number, scattering angle,
			 * distance toward Earth at center of each POS pixel; set flag
			 * posbnd if any model portion extends beyond POS frame limits.*/
			/* NOTE: Limited to single component for now */

			if (posvis_cuda_2(dpar, dmod, ddat, orbit_offset, s, f, 0, 0, 0) &&
					v == v0_index) {
				/* Call single-threaded kernel to set dpar->posbnd and
				 * dpar->posbnd_logfactor */
				cf_set_posbnd_krnl<<<1,1>>>(dpar);
				checkErrorAfterKernelLaunch("cf_deldop_set_posbnd_krnl (calc_fits_cuda)");
			}

			/* Launch single-threaded kernel to get dpar->exclude_seen */
			cf_get_exclude_seen_krnl<<<1,1>>>(dpar, dmod);
			checkErrorAfterKernelLaunch("cf_get_exclude_seen_krnl (calc_fits_cuda)");
			gpuErrchk(cudaMemcpyFromSymbol(&exclude_seen, cf_exclude_seen, sizeof(int),
					0, cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpyFromSymbol(&xlim0, cf_xlim0, sizeof(int),
					0, cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpyFromSymbol(&xlim1, cf_xlim1, sizeof(int),
					0, cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpyFromSymbol(&ylim0, cf_ylim0, sizeof(int),
					0, cudaMemcpyDeviceToHost));
			gpuErrchk(cudaMemcpyFromSymbol(&ylim1, cf_ylim1, sizeof(int),
					0, cudaMemcpyDeviceToHost));

			/* Go through all POS pixels visible with low enough scattering
			 * angle and mark the facets which project onto their centers as
			 * having been "seen" at least once                        */
			/* I'll launch a multi-threaded kernel here:
			 * (xlim1 - xlim0 + 1)^2 threads			 */
			if (s != exclude_seen && v == v0_index) {

				int xspan = xlim1 - xlim0 + 1;
				int yspan = ylim1 - ylim0 + 1;
				int nThreads = xspan * yspan;

				/* Configure & launch posclr_krnl to initialize POS view */
				BLK.x = floor((maxThreadsPerBlock - 1 + nThreads) /
						maxThreadsPerBlock);
				THD.x = maxThreadsPerBlock; // Thread block dimensions
				cf_mark_pixels_seen_krnl<<<BLK,THD>>>(dpar, dmod,
						nThreads, xspan);
				checkErrorAfterKernelLaunch("cf_mark_pixels_krnl (calc_fits_cuda)");
			}

			/* Zero out fit Doppler spectrum, then call pos2doppler to create
			 * the fit image by mapping power from the plane of the sky to
			 * Doppler space.                             */
			/* Zero out the fit delay-Doppler image, then call pos2deldop to
			 * create the fit image by mapping power from the plane of the sky
			 * to delay-Doppler space.                             */
			BLK.x = floor((maxThreadsPerBlock - 1 + ndop) /
					maxThreadsPerBlock);
			THD.x = maxThreadsPerBlock; // Thread block dimensions
			clrvect_krnl<<<BLK,THD>>>(ddat, s, f, ndop);
			checkErrorAfterKernelLaunch("clrvect_krnl, calc_fits_cuda:1060");

			if (pos2doppler_cuda_2(dpar,dmod,ddat,0.0,0.0,0.0,0, s,f,v)) {
				/* Call single-threaded kernel to set badradar flag and
				 * associated badradar_logfactor			 */
				cf_set_badradar_krnl<<<1,1>>>(dpar);
				checkErrorAfterKernelLaunch("cf_set_badradar_krnl (calc_fits_cuda)");
			}

			/*  If smearing is being modeled, include the Doppler
          calculations from this view in the summed results for this frame  */

			if (nviews > 1) {
				/* Launch ndel*ndop-threaded kernel to add fit[i][j] to
				 * fit_store[i][j]*/
				BLK.x = floor((maxThreadsPerBlock - 1 + ndop) /
						maxThreadsPerBlock);
				THD.x = maxThreadsPerBlock; // Thread block dimensions
				cf_add_fit_store_krnl1<<<BLK,THD>>>(ddat,fit_store,ndop,s,f);
				checkErrorAfterKernelLaunch("cf_add_fit_store_krnl1 (calc_fits_cuda)");
				cf_add_fit_store_krnl2<<<1,1>>>();
				checkErrorAfterKernelLaunch("cf_add_fit_store_krnl2 (calc_fits_cuda)");
			}
		}

		/* If smearing is being modeled, compute mean values over all views for
		 * this frame and store them in the standard frame structure     */
		/* This kernel also carries out the gamma transformation on the fit
		 * image if the par->dd_gamma flag is not set  */
		if (nviews > 1) {
			/* Launch ndop-threaded kernel to add fit[i] to
			 * fit_store[i]*/
			BLK.x = floor((maxThreadsPerBlock - 1 + ndop) /
					maxThreadsPerBlock);
			THD.x = maxThreadsPerBlock; // Thread block dimensions
			cf_finish_fit_store_krnl<<<BLK,THD>>>(ddat, fit_store,s,f,ndop);
			checkErrorAfterKernelLaunch("cf_deldop_finish_fit_store (calc_fits_cuda)");
			cf_finish_fit_krnl2<<<1,1>>>();
			checkErrorAfterKernelLaunch("cf_finish_fit_krnl2");
			cf_gamma_trans_krnl<<<BLK,THD>>>(dpar, ddat, s, f, ndop);
			checkErrorAfterKernelLaunch("cf_gamma_trans_krnl");
			cudaFree(fit_store);
		}
	}  /* end loop over frames */
}


//void calc_poset( struct par_t *par, struct mod_t *mod, struct poset_t *poset,
//		int s)
//{
//	const char *monthName[12] = {"Jan", "Feb", "Mar", "Apr", "May", "Jun",
//			"Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};
//	double orbit_offset[3] = {0.0, 0.0, 0.0};
//
//	FILE *fpopt;
//	char tempstring[MAXLEN], name[MAXLEN];
//	int year, mon, day, hour, min, sec, f, c, i, j, k, l, nrow_fit, ncol_fit, n_pos,
//	facetnum, x, y, v, v2;
//	double w[3], spin_colat, spin_azim, xoff, yoff, resamp_fact, resamp_x0, resamp_y0,
//	xcom_fit, ycom_fit, resamp_xwidth, resamp_ywidth, resamp_angle, oa[3][3],
//	to_earth[3], to_earth_lat, to_earth_long, rotphase, sa[3][3], to_sun[3],
//	to_sun_lat, to_sun_long, pab[3], pab_lat, pab_long, intensityfactor,
//	phi, theta, psi, intspin_body[3], badposet_logfactor_view;
//	double **fit_store;
//	struct posetfrm_t *frame;
//	struct posetview_t *view0;
//	struct pos_t *pos;
//
//	for (f=0; f<poset->nframes; f++) {
//
//		frame = &poset->frame[f];
//		view0 = &frame->view[poset->v0];
//		pos = &frame->pos;
//
//		ncol_fit = frame->ncol;
//		nrow_fit = frame->nrow;
//
//		/*  If smearing is being modeled, initialize variables that
//        will be used to sum results calculated for individual views  */
//
//		if (poset->nviews > 1) {
//			fit_store = matrix( 1, ncol_fit, 1, nrow_fit);
//			for (i=1; i<=ncol_fit; i++)
//				for (j=1; j<=nrow_fit; j++)
//					fit_store[i][j] = 0.0;
//		}
//
//		/*  Loop over all views for this (smeared) frame, going in an order that
//        ends with the view corresponding to the epoch listed for this frame
//        in the obs file; this way we can use the calculated information for
//        that view in the "write" action screen and disk output that follows   */
//
//		for (v2=poset->v0+1; v2<=poset->v0+poset->nviews; v2++) {
//			v = v2 % poset->nviews;
//
//			for (i=0; i<=2; i++)
//				for (j=0; j<=2; j++) {
//					pos->ae[i][j] = frame->view[v].ae[i][j];
//					pos->oe[i][j] = frame->view[v].oe[i][j];
//					pos->se[i][j] = frame->view[v].se[i][j];
//				}
//			pos->bistatic = 1;
//
//			/*  Initialize the plane-of-sky view  */
//
//			posclr( pos);
//
//			/*  Call routine posvis to get the facet number, scattering angle,
//          incidence angle, and distance toward Earth at the center of
//          each POS pixel; set the posbnd parameter to 1 if any portion
//          of the model extends beyond the POS frame limits.              */
//
//			for (c=0; c<mod->shape.ncomp; c++)
//				if (posvis( &mod->shape.comp[c].real, orbit_offset, pos,
//						(int) par->pos_smooth, 0, 0, c) && v == poset->v0) {
//					par->posbnd = 1;
//					if (pos->bistatic)
//						par->posbnd_logfactor += 0.5 * frame->dof * pos->posbnd_logfactor;
//					else
//						par->posbnd_logfactor += frame->dof * pos->posbnd_logfactor;
//				}
//
//			/*  Now view the model from the source (sun) and get the facet number
//          and distance toward the source of each pixel in this projected view;
//          use this information to determine which POS pixels are shadowed       */
//
//			if (pos->bistatic) {
//				for (c=0; c<mod->shape.ncomp; c++)
//					if (posvis( &mod->shape.comp[c].real, orbit_offset, pos,
//							0, 1, 0, c)) {
//						par->posbnd = 1;
//						par->posbnd_logfactor += 0.5 * frame->dof * pos->posbnd_logfactor;
//					}
//
//				/*  Identify and mask out shadowed POS pixels  */
//
//				posmask( pos, par->mask_tol);
//			}
//
//			/*  Go through all POS pixels which are visible and unshadowed with
//          sufficiently low scattering and incidence angles, and mark the facets
//          which project onto their centers as having been "seen" at least once   */
//
//			if (s != par->exclude_seen && v == poset->v0) {
//				for (k=pos->xlim[0]; k<=pos->xlim[1]; k++)
//					for (l=pos->ylim[0]; l<=pos->ylim[1]; l++) {
//						if ((pos->cose[k][l] > par->mincosine_seen)
//								&& (pos->cosi[k][l] > par->mincosine_seen)
//								&& (pos->f[k][l] >= 0)) {
//							facetnum = pos->f[k][l];
//							c = pos->comp[k][l];
//							mod->shape.comp[c].real.f[facetnum].seen = 1;
//						}
//					}
//			}
//
//			/*  Compute the sky rendering  */
//
//			intensityfactor = pow( pos->km_per_pixel/AU, 2.0);
//			apply_photo( mod, poset->ioptlaw, frame->view[v].solar_phase,
//					intensityfactor, pos, 0);
//
//			/*  Resample the sky rendering to get the model plane-of-sky image    */
//			/*  (if using bicubic interpolation or cubic convolution, force       */
//			/*  all model pixel values to be nonnegative)                         */
//			/*                                                                    */
//			/*  Implement the x and y COM offsets, xoff and yoff, by first        */
//			/*  using them to compute xcom_fit and ycom_fit -- the COM position   */
//			/*  in the fit image, relative to the center of the fit image -- and  */
//			/*  then shifting the resampled region in the *opposite* direction    */
//			/*  by the appropriate proportional amount.  Then implement the       */
//			/*  "northangle" setting (clockwise heading of north) by rotating     */
//			/*  the resampling grid *counterclockwise* by northangle.             */
//
//			n_pos = pos->n;
//			xoff = frame->off[0].val;
//			yoff = frame->off[1].val;
//			xcom_fit = (frame->colcom_vig - (ncol_fit + 1)/2.0) + xoff;
//			ycom_fit = (frame->rowcom_vig - (nrow_fit + 1)/2.0) + yoff;
//			resamp_fact = frame->fit.km_per_pixel / pos->km_per_pixel;
//			resamp_x0 = -xcom_fit*resamp_fact;
//			resamp_y0 = -ycom_fit*resamp_fact;
//			resamp_xwidth = resamp_fact*(ncol_fit - 1);
//			resamp_ywidth = resamp_fact*(nrow_fit - 1);
//			resamp_angle = -frame->northangle;
//			resampim( frame->pos.b, -n_pos, n_pos, -n_pos, n_pos,
//					frame->fit.b, 1, ncol_fit, 1, nrow_fit,
//					resamp_x0, resamp_xwidth, resamp_y0, resamp_ywidth, resamp_angle,
//					(int) par->poset_resample, (int) par->image_rebin);
//			if (par->poset_resample == BICUBIC || par->poset_resample == CUBICCONV) {
//				for (k=1; k<=ncol_fit; k++)
//					for (l=1; l<=nrow_fit; l++)
//						frame->fit.b[k][l] = MAX( 0.0, frame->fit.b[k][l]);
//			}
//
//			/*  Set the badposet flag and increase badposet_logfactor if the model   */
//			/*  plane-of-sky image is too small to "contain" all of the sky          */
//			/*  rendering's nonzero pixels.                                          */
//
//			if (checkposet( pos->b, -n_pos, n_pos, -n_pos, n_pos,
//					resamp_x0, resamp_xwidth, resamp_y0, resamp_ywidth, resamp_angle,
//					&badposet_logfactor_view)) {
//				par->badposet = 1;
//				par->badposet_logfactor += frame->dof * badposet_logfactor_view
//						/ poset->nviews;
//			}
//
//			/*  If smearing is being modeled, include the plane-of-sky
//          calculations from this view in the summed results for this frame  */
//
//			if (poset->nviews > 1)
//				for (i=1; i<=ncol_fit; i++)
//					for (j=1; j<=nrow_fit; j++)
//						fit_store[i][j] += frame->fit.b[i][j];
//
//		}
//
//		/*  If smearing is being modeled, compute mean values over all views
//        for this frame and store them in the standard frame structure     */
//
//		if (poset->nviews > 1) {
//			for (i=1; i<=ncol_fit; i++)
//				for (j=1; j<=nrow_fit; j++)
//					frame->fit.b[i][j] = fit_store[i][j] / poset->nviews;
//			free_matrix( fit_store, 1, ncol_fit, 1, nrow_fit);
//		}
//
//
//	}  /* end loop over frames */
//}
//
//
//void calc_lghtcrv( struct par_t *par, struct mod_t *mod, struct lghtcrv_t *lghtcrv,
//		int s)
//{
//	const char *monthName[12] = {"Jan", "Feb", "Mar", "Apr", "May", "Jun",
//			"Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};
//	double orbit_offset[3] = {0.0, 0.0, 0.0};
//
//	FILE *fpopt;
//	char tempstring[MAXLEN], name[MAXLEN];
//	int year, mon, day, hour, min, sec, n, ncalc, c, i, i_mid, j, k, l, facetnum,
//	n_cross360, n_projectedpixels, n_shadowedpixels, x, y, v;
//	double epoch_mid, epoch_diff_min, epoch_diff, w[3], spin_colat, spin_azim, oa[3][3],
//	rotphase, sa[3][3], to_sun[3], to_sun_lat, to_sun_long, pab[3], pab_lat,
//	pab_long, intensityfactor, phi, theta, psi, intspin_body[3], posbnd_logfactor,
//	projected_area, lambertdisk_intensity, interp;
//	double **to_earth, *to_earth_lat, *to_earth_long, *rotphase_unwrapped;
//	struct crvrend_t *rend;
//	struct pos_t *pos;
//
//	/*  Initialize variables to avoid compilation warning  */
//
//	i_mid = 0;
//	epoch_mid = epoch_diff = epoch_diff_min = 0.0;
//	n_cross360 = 0;
//	to_earth = NULL;
//	to_earth_lat = to_earth_long = rotphase_unwrapped = NULL;
//
//	/*  Initialize variables dealing with bad models  */
//
//	posbnd_logfactor = 0.0;
//
//	/*  Get n, the number of observed points for this lightcurve,
//      and ncalc, the number of epochs at which model lightcurve
//      brightnesses are to be computed                            */
//
//	n = lghtcrv->n;
//	ncalc = lghtcrv->ncalc;
//
//	/*  Calculate the model lightcurve values at each of the user-specified
//      epochs x[i], with i=1,2,...,ncalc; these may or may not be the same as the
//      epochs t[i] (i=1,2,...,n) at which actual lightcurve observations were made.  */
//
//	for (i=1; i<=ncalc; i++) {
//
//		rend = &lghtcrv->rend[i];
//		pos = &rend->pos;
//
//		for (j=0; j<=2; j++)
//			for (k=0; k<=2; k++) {
//				pos->ae[j][k] = rend->ae[j][k];
//				pos->oe[j][k] = rend->oe[j][k];
//				pos->se[j][k] = rend->se[j][k];
//			}
//		pos->bistatic = 1;
//
//		/*  Initialize the plane-of-sky view  */
//
//		posclr( pos);
//
//		/*  Call routine posvis to get the facet number, scattering angle,
//        incidence angle, and distance toward Earth at the center of
//        each POS pixel; set the posbnd parameter to 1 if any portion
//        of the model extends beyond the POS frame limits.              */
//
//		for (c=0; c<mod->shape.ncomp; c++)
//			if (posvis( &mod->shape.comp[c].real, orbit_offset, pos,
//					(int) par->pos_smooth, 0, 0, c)) {
//				par->posbnd = 1;
//				if (pos->bistatic)
//					posbnd_logfactor += 0.5 * pos->posbnd_logfactor;
//				else
//					posbnd_logfactor += pos->posbnd_logfactor;
//			}
//
//		/*  Now view the model from the source (sun) and get the facet number
//        and distance toward the source of each pixel in this projected view;
//        use this information to determine which POS pixels are shadowed       */
//
//		if (pos->bistatic) {
//			for (c=0; c<mod->shape.ncomp; c++)
//				if (posvis( &mod->shape.comp[c].real, orbit_offset, pos,
//						0, 1, 0, c)) {
//					par->posbnd = 1;
//					posbnd_logfactor += 0.5 * pos->posbnd_logfactor;
//				}
//
//			/*  Identify and mask out shadowed POS pixels  */
//
//			posmask( pos, par->mask_tol);
//		}
//
//		/*  Go through all POS pixels which are visible and unshadowed with
//        sufficiently low scattering and incidence angles, and mark the facets
//        which project onto their centers as having been "seen" at least once   */
//
//		if (s != par->exclude_seen) {
//			for (k=pos->xlim[0]; k<=pos->xlim[1]; k++)
//				for (l=pos->ylim[0]; l<=pos->ylim[1]; l++) {
//					if ((pos->cose[k][l] > par->mincosine_seen)
//							&& (pos->cosi[k][l] > par->mincosine_seen)
//							&& (pos->f[k][l] >= 0)) {
//						facetnum = pos->f[k][l];
//						c = pos->comp[k][l];
//						mod->shape.comp[c].real.f[facetnum].seen = 1;
//					}
//				}
//		}
//
//		/*  Compute the model brightness for this model lightcurve point  */
//
//		intensityfactor = pow( pos->km_per_pixel/AU, 2.0);
//		lghtcrv->y[i] = apply_photo( mod, lghtcrv->ioptlaw, lghtcrv->solar_phase[i],
//				intensityfactor, pos, 0);
//
//		/*  Finished with this calculated lightcurve point  */
//
//	}
//
//	/*  Now that we have calculated the model lightcurve brightnesses y at each
//      of the epochs x, we use cubic spline interpolation (Numerical Recipes
//      routines spline and splint) to get model lightcurve brightness fit[i]
//      at each OBSERVATION epoch t[i], with i=1,2,...,n.  This will allow us
//      (in routine chi2) to compare model to data (fit[i] to obs[i]) to get
//      chi-square.  Note that vector y2 contains the second derivatives of
//      the interpolating function at the calculation epochs x.
//
//      Smearing is handled by interpolating the brightness at the time t of
//      each individual view and then taking the mean of all views that
//      correspond to a given observed lightcurve point.                         */
//
//	spline( lghtcrv->x, lghtcrv->y, ncalc, 2.0e30, 2.0e30, lghtcrv->y2);
//	for (i=1; i<=n; i++) {
//		lghtcrv->fit[i] = 0.0;
//		for (v=0; v<lghtcrv->nviews; v++) {
//			splint( lghtcrv->x, lghtcrv->y, lghtcrv->y2, ncalc,
//					lghtcrv->t[i][v], &interp);
//			lghtcrv->fit[i] += interp;
//		}
//		lghtcrv->fit[i] /= lghtcrv->nviews;
//	}
//
//	/*  Deal with flags for model that extends beyond the POS frame  */
//
//	par->posbnd_logfactor += lghtcrv->dof * (posbnd_logfactor/ncalc);
//
//}
