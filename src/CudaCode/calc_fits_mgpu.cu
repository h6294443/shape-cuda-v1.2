/*****************************************************************************************
                                                                              calc_fits.c
This version of calc_fits runs on a dual-GPU setup.  Each GPU has its own
memory space and streams.  GPU0 will handle the even frames (for radar) and odd
frames for lightcurves (essentially, GPU0 always handles the first frame in a
set, then jumps over the next frame to do the one after that, and so on.  GPU1
handles the rest of the frames.  As such, both GPUs' calculations are inter-
woven to work on alternating frames.  Each GPU gets its own half-array copies
of variables.  This carries through to host code as well.
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
#include "../shape/head.h"
}
/* Sub-function prototypes for the main function */
__host__ void calc_deldop_mgpu(struct par_t *dpar, struct mod_t *dmod, struct
		dat_t *ddat, int s, int nfrm_alloc, int nfrm_half0, int nfrm_half1,
		dim3 BLK_half0, dim3 BLK_half1,	int nviews,	unsigned char type,	int nf,
		struct pos_t **pos0, struct pos_t **pos1, struct deldopfrm_t **frame0,
		struct deldopfrm_t **frame1, struct deldopview_t **view0_0,	struct
		deldopview_t **view0_1,	float **fit_store, int *ndel0, int *ndel1, int
		*ndop0, int *ndop1, int *posn0,	int *posn1,	int4 *xylim0, int4 *xylim1,
		float *overflow,int *outbndarr0,
		int *outbndarr1, cudaStream_t *gpu0_stream, cudaStream_t *gpu1_stream);

__host__ void calc_doppler_mgpu(struct par_t *dpar,	struct mod_t *dmod, struct
		dat_t *ddat, int s,	int nfrm_alloc, int nfrm_half0, int nfrm_half1,
		dim3 BLK_half0, dim3 BLK_half1,	int nviews,	unsigned char type,	int nf,
		struct pos_t **pos0, struct pos_t **pos1, struct dopfrm_t **frame0,
		struct dopfrm_t **frame1, struct dopview_t **view0_0, struct dopview_t
		**view0_1, float **fit_store, int *ndop0, int *ndop1, int *posn0, int
		*posn1, int4 *xylim0, int4 *xylim1,	float *overflow, int *outbndarr0,
		int *outbndarr1, cudaStream_t *gpu0_stream,	cudaStream_t *gpu1_stream);

__host__ void calc_lghtcrv_mgpu(struct par_t *dpar, struct mod_t *dmod, struct
		dat_t *ddat, int s, int nfrm_alloc, int nfrm_half0, int nfrm_half1, dim3
		BLK_half0, dim3 BLK_half1, int nviews, unsigned char type, int lc_n, int
		nf, struct pos_t **pos0, struct pos_t **pos1, struct crvrend_t **rend0,
		struct crvrend_t **rend1, int *posn0, int *posn1, float *pixels_per_km0,
		float *pixels_per_km1, float *overflow, int4 *xylim0, int4 *xylim1, int
		*outbndarr0, int *outbndarr1, double3 *so0, double3 *so1, double *u,
		cudaStream_t *gpu0_stream, cudaStream_t *gpu1_stream);

//__host__ void calc_deldop_mgpu(struct par_t *dpar, struct mod_t *dmod,
//		struct dat_t *ddat, struct vertices_t **verts, int s, int nframes,
//		int nviews, unsigned char type, int nf, struct pos_t **pos, struct
//		deldopfrm_t **frame, struct deldopview_t **view0, float **fit_store,
//		int *ndel, int *ndop, int *posn, int *bistatic, int4 *xylim, float
//		*overflow, int *outbndarr, cudaStream_t *cf_stream);
//__host__ void calc_doppler_mgpu(struct par_t *dpar, struct mod_t *dmod,
//		struct dat_t *ddat, struct vertices_t **verts, int s, int nframes, int
//		nviews, unsigned char type, int nf, struct pos_t **pos, struct dopfrm_t
//		**frame, struct dopview_t **view0, float **fit_store, int *ndop, int
//		*posn, int *bistatic, int4 *xylim, float *overflow, int	*outbndarr,
//		cudaStream_t *cf_stream);
//__host__ void calc_lghtcrv_mgpu(struct par_t *dpar, struct mod_t *dmod,
//		struct dat_t *ddat, struct vertices_t **verts, int s, int nframes, int
//		nviews, unsigned char type, int lc_n, int nf, struct pos_t **pos,
//		struct crvrend_t **rend, int *posn, int *bistatic, float *pxlpkm, float
//		*overflow, int4 *xylim, int *outbndarr, double3 *so, double *u,
//		cudaStream_t *cf_stream);
//
__device__ int mgpu_v0_index, mgpu_exclude_seen, lc_bistatic_all=0;
__device__ double mgpu_lghtcrv_posbnd_logfactor;

__host__ void calc_fits_mgpu(
		struct par_t *dpar,
		struct mod_t *dmod,
		struct dat_t *ddat,
		int *nviews,
		int *hnframes,
		int *lc_n,
		unsigned char *type,
		int nsets,
		int nf,
		cudaStream_t *gpu0_stream,
		cudaStream_t *gpu1_stream,
		int max_frames)
{
	int s, nfrm_alloc, *ndel0, *ndel1, *ndop0, *ndop1, *posn0, *posn1,
		*outbndarr0, *outbndarr1, nfrm_alloc_max;
	int4 *xylim0, *xylim1;
	float **fit_store, *overflow, *pixels_per_km0, *pixels_per_km1;
	double *u;
	double3 *so0, *so1;
	nfrm_alloc_max = max_frames+1;
	int nfrm_half0, nfrm_half0_max = nfrm_alloc_max/2 + nfrm_alloc_max%2;
	int nfrm_half1, nfrm_half1_max = nfrm_alloc_max/2;

	struct pos_t **pos0, **pos1;
	struct deldopfrm_t **ddframe0, **ddframe1;
	struct dopfrm_t **dframe0, **dframe1;
	struct crvrend_t **rend0, **rend1;
	struct deldopview_t **ddview0_0, **ddview0_1;
	struct dopview_t **dview0_0, **dview0_1;
	dim3 BLK,BLK_half0, BLK_half1, THD, THD64;
	THD.x = maxThreadsPerBlock;
	THD64.x = 64;

	/* GPU0 allocations */
	gpuErrchk(cudaSetDevice(GPU0));
	gpuErrchk(cudaMalloc((void**)&pos0, sizeof(pos_t*) * nfrm_half0_max));
	gpuErrchk(cudaMalloc((void**)&ddframe0, sizeof(deldopfrm_t*) * nfrm_half0_max));
	gpuErrchk(cudaMalloc((void**)&ddview0_0, sizeof(deldopview_t*) * nfrm_half0_max));
	gpuErrchk(cudaMalloc((void**)&dframe0, sizeof(dopfrm_t*) * nfrm_half0_max));
	gpuErrchk(cudaMalloc((void**)&dview0_0, sizeof(dopview_t*) * nfrm_half0_max));
	gpuErrchk(cudaMalloc((void**)&rend0, sizeof(crvrend_t*) * nfrm_half0_max));
	gpuErrchk(cudaMalloc((void**)&ndel0, sizeof(int) * nfrm_half0_max));
	gpuErrchk(cudaMalloc((void**)&ndop0, sizeof(int) * nfrm_half0_max));
	gpuErrchk(cudaMalloc((void**)&posn0, sizeof(int) * nfrm_half0_max));
	gpuErrchk(cudaMalloc((void**)&xylim0, sizeof(int4) * nfrm_half0_max));
	gpuErrchk(cudaMalloc((void**)&outbndarr0, sizeof(int) * nfrm_half0_max));
	gpuErrchk(cudaMalloc((void**)&so0, 	 	sizeof(double3) *((nfrm_half0_max*3)+1)));
	gpuErrchk(cudaMalloc((void**)&overflow, sizeof(float) * 6));
	gpuErrchk(cudaMalloc((void**)&u, sizeof(double) * nfrm_alloc_max));
	gpuErrchk(cudaMalloc((void**)&pixels_per_km0,   sizeof(float) * nfrm_half0_max));
	/* The following is allocated as Cuda managed memory because it's a double-pointer */
	cudaCalloc((void**)&fit_store, sizeof(float*), nfrm_alloc_max);

	/* GPU1 allocations */
	gpuErrchk(cudaSetDevice(GPU1));
	gpuErrchk(cudaMalloc((void**)&pos1, sizeof(pos_t*) * nfrm_half1_max));
	gpuErrchk(cudaMalloc((void**)&ddframe1, sizeof(deldopfrm_t*) * nfrm_half1_max));
	gpuErrchk(cudaMalloc((void**)&ddview0_1, sizeof(deldopview_t*) * nfrm_half1_max));
	gpuErrchk(cudaMalloc((void**)&dframe1, sizeof(dopfrm_t*) * nfrm_half1_max));
	gpuErrchk(cudaMalloc((void**)&dview0_1, sizeof(dopview_t*) * nfrm_half1_max));
	gpuErrchk(cudaMalloc((void**)&rend1, sizeof(crvrend_t*) * nfrm_half1_max));
	gpuErrchk(cudaMalloc((void**)&ndel1, sizeof(int) * nfrm_half1_max));
	gpuErrchk(cudaMalloc((void**)&ndop1, sizeof(int) * nfrm_half1_max));
	gpuErrchk(cudaMalloc((void**)&posn1, sizeof(int) * nfrm_half1_max));
	gpuErrchk(cudaMalloc((void**)&xylim1, sizeof(int4) * nfrm_half1_max));
	gpuErrchk(cudaMalloc((void**)&outbndarr1, sizeof(int) * nfrm_half1_max));
	gpuErrchk(cudaMalloc((void**)&pixels_per_km1,   sizeof(float) * nfrm_half1_max));
	gpuErrchk(cudaSetDevice(GPU0));

	/* Initialize flags that indicate the model extends beyond POS frame, that
	 * plane-of-sky fit images are too small to "contain" the target, and that
	 * model is too wide in (delay-)Doppler space to create (delay-)Doppler fit
	 * frames.  Note that this also gets mod->shape.nf and nsets            */

	cfs_init_devpar_krnl<<<1,1>>>(dpar);
	checkErrorAfterKernelLaunch("cfs_init_devpar_krnl2");

	/* Initialize the flags that indicate whether or not each facet of each
	 * model component is ever visible and unshadowed from Earth
	 * Note:  Single component only for now.  */
	//for (c=0; c<mod->shape.ncomp; c++)
	BLK.x = floor((THD.x - 1 + nf)/THD.x);
	cf_init_seen_flags_krnl<<<BLK,THD>>>(dmod,nf);
	checkErrorAfterKernelLaunch("cf_init_seen_flags_krnl (calc_fits_cuda_streams)");

	/* Calculate the fits for each dataset in turn  */
	for (s=0; s<nsets; s++) {
		/* Figure out this set's launch parameters and frame splits */
		if (type[s]==LGHTCRV)			nfrm_alloc = hnframes[s] + 1;
		else							nfrm_alloc = hnframes[s];

		nfrm_half0 = nfrm_alloc/2 + nfrm_alloc%2;
		nfrm_half1 = nfrm_alloc/2;
		BLK_half0 = floor ((THD64.x - 1 + nfrm_half0) / THD64.x);
		BLK_half1 = floor ((THD64.x - 1 + nfrm_half1) / THD64.x);

		switch (type[s]) {
		case DELAY:
			calc_deldop_mgpu(dpar, dmod, ddat, s, nfrm_alloc, nfrm_half0,
					nfrm_half1, BLK_half0, BLK_half1, nviews[s], type[s], nf,
					pos0, pos1, ddframe0, ddframe1, ddview0_0, ddview0_1,
					fit_store, ndel0, ndel1, ndop0, ndop1, posn0, posn1, xylim0,
					xylim1,	overflow, outbndarr0, outbndarr1, gpu0_stream,
					gpu1_stream);
//			calc_deldop_gpu(dpar, dmod, ddat, verts, s, nframes[s],
//					nviews[s], type[s], nf, pos, ddframe, ddview0, fit_store,
//					ndel, ndop,	posn, bistatic, xylim, overflow, outbndarr,
//					cf_stream);
			break;
		case DOPPLER:
			calc_doppler_mgpu(dpar,	dmod, ddat, s, nfrm_alloc, nfrm_half0,
					nfrm_half1, BLK_half0, BLK_half1, nviews[s], type[s], nf,
					pos0, pos1, dframe0, dframe1, dview0_0, dview0_1, fit_store,
					ndop0, ndop1,  posn0, posn1, xylim0, xylim1, overflow,
					outbndarr0, outbndarr1, gpu0_stream, gpu1_stream);
//			calc_doppler_gpu(dpar, dmod, ddat, verts, s, nframes[s],
//					nviews[s], type[s], nf, pos, dframe, dview0, fit_store,
//					ndop, posn, bistatic, xylim, overflow, outbndarr,
//					cf_stream );
			break;
//		case POS:
//			printf("Write calc_poset_cuda!");
////			calc_poset_cuda(dpar, dmod, s);
//			break;
		case LGHTCRV:
			calc_lghtcrv_mgpu(dpar, dmod, ddat, s, nfrm_alloc, nfrm_half0,
					nfrm_half1, BLK_half0, BLK_half1, nviews[s], type[s],
					lc_n[s], nf, pos0, pos1, rend0, rend1, posn0, posn1,
					pixels_per_km0, pixels_per_km1, overflow, xylim0, xylim1,
					outbndarr0, outbndarr1, so0, so1, u, gpu0_stream, gpu1_stream);
//			calc_lghtcrv_gpu(dpar, dmod, ddat, verts, s, nframes[s],
//					nviews[s], type[s], lc_n[s],nf, pos, rend, posn, bistatic,
//					pxlpkm, overflow, xylim, outbndarr, so, u, cf_stream);
			break;
		default:
			printf("calc_fits_mgpu.c: can't handle this type yet\n");
		}
	}

	/* Complete calculations of values that will be used during a fit to
	 * increase the objective function for models with bad properties   */
	cf_set_final_pars_krnl<<<1,1>>>(dpar, ddat);
	checkErrorAfterKernelLaunch("cf_set_final_pars_krnl (calc_fits_cuda)");

	/* GPU0 de-allocations */				/* GPU1 de-allocations */
	gpuErrchk(cudaSetDevice(GPU0));			gpuErrchk(cudaSetDevice(GPU1));
	gpuErrchk(cudaFree(pos0));				gpuErrchk(cudaFree(pos1));
	gpuErrchk(cudaFree(ddframe0));			gpuErrchk(cudaFree(ddframe1));
	gpuErrchk(cudaFree(ddview0_0));			gpuErrchk(cudaFree(ddview0_1));
	gpuErrchk(cudaFree(dframe0));			gpuErrchk(cudaFree(dframe1));
	gpuErrchk(cudaFree(dview0_0));			gpuErrchk(cudaFree(dview0_1));
	gpuErrchk(cudaFree(rend0));				gpuErrchk(cudaFree(rend1));
	gpuErrchk(cudaFree(ndel0));				gpuErrchk(cudaFree(ndel1));
	gpuErrchk(cudaFree(ndop0));				gpuErrchk(cudaFree(ndop1));
	gpuErrchk(cudaFree(posn0));				gpuErrchk(cudaFree(posn1));
	gpuErrchk(cudaFree(xylim0));			gpuErrchk(cudaFree(xylim1));
	gpuErrchk(cudaFree(outbndarr0));		gpuErrchk(cudaFree(outbndarr1));
	gpuErrchk(cudaFree(pixels_per_km0));	gpuErrchk(cudaFree(pixels_per_km1));
	gpuErrchk(cudaFree(so0));				gpuErrchk(cudaSetDevice(GPU0));
	gpuErrchk(cudaFree(overflow));
	gpuErrchk(cudaFree(u));
	gpuErrchk(cudaFree(fit_store));


}

__global__ void init_deldop_calcfits_krnl(
		struct dat_t *ddat,
		struct deldopfrm_t **frame,
		struct pos_t **pos,
		struct deldopview_t **view0,
		int *ndop,
		int *ndel,
		int *posn,
		int s, int size, int oddflg) {

	/* nfrm_half0/nfrm_half1-threaded kernel for dual-GPU operation */
	/* 	overflow[0] - overflow_o2_store
	 * 	overflow[1] - overflow_m2_store
	 * 	overflow[2] - overflow_xsec_store
	 * 	overflow[3] - overflow_dopmean_store
	 * 	overflow[4] - overflow_delmean_store
	 */
	int hf = blockIdx.x * blockDim.x + threadIdx.x;
	int f = 2 * hf + oddflg;

	if (hf < size) {
		frame[hf] 	= &ddat->set[s].desc.deldop.frame[f];
		ndop[hf]  	= frame[hf]->ndop;
		ndel[hf]	= frame[hf]->ndel;
		view0[hf] 	= &frame[hf]->view[ddat->set[s].desc.deldop.v0];
		pos[hf]		= &frame[hf]->pos;
		posn[hf]    = pos[hf]->n;
	}
	__syncthreads();
	if (hf == 0)
		mgpu_v0_index= ddat->set[s].desc.deldop.v0;
}

__global__ void init_doppler_calcfits_krnl(
		struct dat_t *ddat,
		struct dopfrm_t **frame,
		struct pos_t **pos,
		struct dopview_t **view0,
		int *ndop,
		int *posn,
		int s, int size, int oddflg) {

	/* nfrm_half0/nfrm_half1-threaded kernel for dual-GPU operation */
	int hf = blockIdx.x * blockDim.x + threadIdx.x;
	int f = 2 * hf + oddflg;

	if (hf < size) {
		frame[hf] = &ddat->set[s].desc.doppler.frame[f];
		view0[hf] = &frame[hf]->view[ddat->set[s].desc.doppler.v0];
		ndop[hf]	 = frame[hf]->ndop;
		pos[hf]	 = &frame[hf]->pos;
		posn[hf]	= pos[hf]->n;
	}
	__syncthreads();
	if (hf == 0)
		mgpu_v0_index  = ddat->set[s].desc.doppler.v0;
}

__global__ void init_lghtcrv_calcfits_krnl(
		struct dat_t *ddat,
		struct crvrend_t **rend,
		struct pos_t **pos,
		int *posn,
		int s,
		int size,
		int oddflg) {

	/* nfrm_half0/nfrm_half1-threaded kernel for dual-GPU operation */
		int hf = blockIdx.x * blockDim.x + threadIdx.x;
		int f = 2 * hf + oddflg + 1;

		if (hf < size) {

			rend[hf] = &ddat->set[s].desc.lghtcrv.rend[f];//&cfs_lghtcrv->rend[f];
			pos[hf] = &rend[hf]->pos;
			posn[hf] = pos[hf]->n;
		}
}

__global__ void set_aeoe_calcfits_krnl(
		struct dat_t *ddat,
		struct pos_t **pos,
		int s,
		unsigned char type,
		int v,
		int size,
		int oddflg,
		int lghtcrvflg) {

	/* nfrm_half0/nfrm_half1-threaded kernel for dual-GPU operation */
	int i, j, hf = blockIdx.x * blockDim.x + threadIdx.x;
	int f = 2 * hf + oddflg + lghtcrvflg;

	if (hf < size) {
		switch (type) {
		case DELAY:
			for (i=0; i<3; i++)
				for (j=0; j<3; j++) {
					pos[hf]->ae[i][j] = ddat->set[s].desc.deldop.frame[f].view[v].ae[i][j];
					pos[hf]->oe[i][j] = ddat->set[s].desc.deldop.frame[f].view[v].oe[i][j];
				}
			break;
		case DOPPLER:
			for (i=0; i<3; i++)
				for (j=0; j<3; j++) {
					pos[hf]->ae[i][j] = ddat->set[s].desc.doppler.frame[f].view[v].ae[i][j];
					pos[hf]->oe[i][j] = ddat->set[s].desc.doppler.frame[f].view[v].oe[i][j];
				}
			break;
		case LGHTCRV:
			for (i=0; i<3; i++)
				for (j=0; j<3; j++) {
					pos[hf]->ae[i][j] = ddat->set[s].desc.lghtcrv.rend[f].ae[i][j];
					pos[hf]->oe[i][j] = ddat->set[s].desc.lghtcrv.rend[f].oe[i][j];
					pos[hf]->se[i][j] = ddat->set[s].desc.lghtcrv.rend[f].se[i][j];
				}
		}/* end of switch */
		/* Single-thread task */
		if (hf == 0) {
			if ((type == LGHTCRV) || (type == POS))
				pos[hf]->bistatic = 1;
			else if ((type == DELAY) || (type == DOPPLER))
				pos[hf]->bistatic = 0;
		}
	}
}

__global__ void set_posbnd_calcfits_krnl(
		struct par_t *dpar,
		struct dat_t *ddat,
		struct pos_t **pos,
		int4 *xylim,
		int *outbndarr,
		int s,
		unsigned char type,
		int size,
		int oddflg,
		int lghtcrvflg) {

	/* nfrm_half0/nfrm_half1-threaded kernel for dual GPU operation. */
	int hf = blockIdx.x * blockDim.x + threadIdx.x;
	int f = 2 * hf + oddflg + lghtcrvflg;

	if (hf < size){
		if (outbndarr[hf]==1) {
			dpar->posbnd = 1;
			switch (type) {
			case DELAY:
				dpar->posbnd_logfactor += ddat->set[s].desc.deldop.frame[f].dof *
				pos[hf]->posbnd_logfactor;
				break;
			case DOPPLER:
				dpar->posbnd_logfactor += ddat->set[s].desc.doppler.frame[f].dof *
				pos[hf]->posbnd_logfactor;
				break;
			case LGHTCRV:
				if (pos[hf]->bistatic) {
					/* Note that mgpu_lghtcrv_posbnd_logfactor is a __device__
					 * variable and therefore implemented independently on each
					 * GPU, under the same name. Use cudaSetDevice(gpuid) to
					 * select which copy is desired 					 */
					mgpu_lghtcrv_posbnd_logfactor += 0.5 * pos[hf]->posbnd_logfactor;
					atomicExch(&lc_bistatic_all, 1);
				}
				else
					mgpu_lghtcrv_posbnd_logfactor += pos[hf]->posbnd_logfactor;
				break;
			}
		}
		xylim[hf].w = pos[hf]->xlim[0];
		xylim[hf].x = pos[hf]->xlim[1];
		xylim[hf].y = pos[hf]->ylim[0];
		xylim[hf].z = pos[hf]->ylim[1];
	}
	__syncthreads();

	if (threadIdx.x == 0)
		mgpu_exclude_seen = dpar->exclude_seen;

}

__global__ void set_badradar_calcfits_krnl(
		struct par_t *dpar,
		struct dat_t *ddat,
		int s,
		int nfrm_alloc,
		int nfrm_half0,
		int nfrm_half1,
		int *outbndarr0,
		int *outbndarr1,
		unsigned char type) {
	/* Single-threaded kernel - using a loop to ensure exclusive access and to
	 * prevent a race condition */
	int hf = 0;

	if (threadIdx.x == 0) {

		switch (type) {
		case DELAY:
			for (int f=0; f<nfrm_alloc; f++) {
				if ((outbndarr0[hf]))	{
					dpar->badradar = 1;
					dpar->badradar_logfactor += ddat->set[s].desc.deldop.frame[f].dof *
							ddat->set[s].desc.deldop.frame[f].badradar_logfactor /
							ddat->set[s].desc.deldop.nviews;
				}
				if (hf<nfrm_half1) {
					f++;	if (f>=nfrm_alloc)	break;
					if ((outbndarr1[hf])) {
						dpar->badradar = 1;
						dpar->badradar_logfactor += ddat->set[s].desc.deldop.frame[f].dof *
								ddat->set[s].desc.deldop.frame[f].badradar_logfactor /
								ddat->set[s].desc.deldop.nviews;
					}
				hf++;
				}
			}
			break;
		case DOPPLER:
			for (int f=0; f<nfrm_alloc; f++) {
				if ((outbndarr0[hf]))	{
					dpar->badradar = 1;
					dpar->badradar_logfactor += ddat->set[s].desc.doppler.frame[f].dof *
							ddat->set[s].desc.doppler.frame[f].badradar_logfactor /
							ddat->set[s].desc.doppler.nviews;
				}
				if (hf<nfrm_half1) {
					f++;	if (f>=nfrm_alloc)	break;
					if ((outbndarr1[hf])) {
						dpar->badradar = 1;
						dpar->badradar_logfactor += ddat->set[s].desc.deldop.frame[f].dof *
								ddat->set[s].desc.deldop.frame[f].badradar_logfactor /
								ddat->set[s].desc.deldop.nviews;
					}
				}
			}
			break;
		}
	}
}

__global__ void finish_overflow_krnl(struct dat_t *ddat,
		float *overflow, int s, unsigned char type, int nfrm_alloc) {
	/* nfrm_alloc-threaded Kernel */
	int f = blockIdx.x * blockDim.x + threadIdx.x;
	int nviews;
	if (f < nfrm_alloc) {
		switch (type) {
		case DELAY:
			nviews = ddat->set[s].desc.deldop.nviews;
			ddat->set[s].desc.deldop.frame[f].overflow_o2 = overflow[0] / nviews;
			ddat->set[s].desc.deldop.frame[f].overflow_m2 = overflow[1] / nviews;
			ddat->set[s].desc.deldop.frame[f].overflow_xsec = overflow[2] / nviews;
			ddat->set[s].desc.deldop.frame[f].overflow_delmean = overflow[3] / nviews;
			ddat->set[s].desc.deldop.frame[f].overflow_dopmean = overflow[4] / nviews;
			break;
		case DOPPLER:
			nviews = ddat->set[s].desc.doppler.nviews;
			ddat->set[s].desc.doppler.frame[f].overflow_o2 = overflow[0] / nviews;
			ddat->set[s].desc.doppler.frame[f].overflow_m2 = overflow[1] / nviews;
			ddat->set[s].desc.doppler.frame[f].overflow_xsec = overflow[2] / nviews;
			ddat->set[s].desc.doppler.frame[f].overflow_dopmean = overflow[3] / nviews;
		}
	}
}

__global__ void gammatrans_calcfits_krnl(struct par_t *dpar, struct dat_t *ddat,
		int s, int f, int nThreads, unsigned char type) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	if (dpar->dd_gamma != 1.0) {
		if (offset < nThreads) {
			/*  Carry out a gamma transformation on the fit image if requested  */
			dev_gamma_trans_f(&ddat->set[s].desc.deldop.frame[f].fit_s[offset],
					__double2float_rn(dpar->dd_gamma));
		}
	}
}

__host__ void calc_deldop_mgpu(
		struct par_t *dpar,
		struct mod_t *dmod,
		struct dat_t *ddat,
		int s,
		int nfrm_alloc, int nfrm_half0, int nfrm_half1,
		dim3 BLK_half0, dim3 BLK_half1,
		int nviews,
		unsigned char type,
		int nf,
		struct pos_t **pos0,			struct pos_t **pos1,
		struct deldopfrm_t **frame0, 	struct deldopfrm_t **frame1,
		struct deldopview_t **view0_0,	struct deldopview_t **view0_1,
		float **fit_store,
		int *ndel0,						int *ndel1,
		int *ndop0,						int *ndop1,
		int *posn0, 					int *posn1,
		int4 *xylim0,					int4 *xylim1,
		float *overflow,
		int *outbndarr0, 				int *outbndarr1,
		cudaStream_t *gpu0_stream, 		cudaStream_t *gpu1_stream)
{
	float3 orbit_off3, orb_xydopoff;
	orbit_off3.x = orbit_off3.y = orbit_off3.z = 0.0;
	orb_xydopoff.x = orb_xydopoff.y = orb_xydopoff.z = 0.0;
	int hf, v0_index, exclude_seen, f, v2, c=0, yspan;
	dim3 BLKdd[nfrm_alloc], BLKpx[nfrm_alloc], THD, THD64, THD9;
	THD.x = maxThreadsPerBlock; THD64.x = 64; THD9.x = 9;
	int4 hxylim0[nfrm_half0], hxylim1[nfrm_half1];
	int hndop0[nfrm_half0], hndop1[nfrm_half1],
		hndel0[nfrm_half0], hndel1[nfrm_half1],
		hposn0[nfrm_half0], hposn1[nfrm_half1],
		/*houtbndarr0[nfrm_half0], houtbndarr1[nfrm_half1],*/
		xspan[nfrm_alloc], nThreadspx[nfrm_alloc], nThreadsdd[nfrm_alloc],
		nThreadspx1[nfrm_alloc], v[nviews+1];

	/* Launch the initialization kernel on each GPU and copy results back to
	 * split host arrays */
	gpuErrchk(cudaSetDevice(GPU0));
	init_deldop_calcfits_krnl<<<BLK_half0,THD64,0,gpu0_stream[0]>>>(ddat,frame0,
			pos0, view0_0, ndop0, ndel0, posn0, s, nfrm_half0, 0);
	gpuErrchk(cudaMemcpyAsync(hndel0, ndel0, sizeof(int)*nfrm_half0,
			cudaMemcpyDeviceToHost, gpu0_stream[0]));
	gpuErrchk(cudaMemcpyAsync(hndop0, ndop0, sizeof(int)*nfrm_half0,
			cudaMemcpyDeviceToHost, gpu0_stream[0]));
	gpuErrchk(cudaMemcpyAsync(hposn0, posn0, sizeof(int)*nfrm_half0,
			cudaMemcpyDeviceToHost, gpu0_stream[0]));

	gpuErrchk(cudaSetDevice(GPU1));
	init_deldop_calcfits_krnl<<<BLK_half1,THD64,0,gpu1_stream[0]>>>(ddat,frame1,
			pos1, view0_1, ndop1, ndel1, posn1, s, nfrm_half1, 1);
	gpuErrchk(cudaMemcpyAsync(hndel1, ndel1, sizeof(int)*nfrm_half1,
			cudaMemcpyDeviceToHost, gpu1_stream[0]));
	gpuErrchk(cudaMemcpyAsync(hndop1, ndop1, sizeof(int)*nfrm_half1,
			cudaMemcpyDeviceToHost, gpu1_stream[0]));
	gpuErrchk(cudaMemcpyAsync(hposn1, posn1, sizeof(int)*nfrm_half1,
			cudaMemcpyDeviceToHost, gpu1_stream[0]));

	gpuErrchk(cudaSetDevice(GPU0));
	gpuErrchk(cudaMemcpyFromSymbol(&v0_index, mgpu_v0_index, sizeof(int),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemsetAsync(overflow, 0, sizeof(float)*6, gpu0_stream[1]));
	checkErrorAfterKernelLaunch("init_deldop_calcfits_krnl");

	hf = 0;
	for (f=0; f<nfrm_alloc; f++) {
		/* Even frame first */
		nThreadsdd[f] = hndel0[hf]*hndop0[hf];
		BLKdd[f].x = floor((THD.x - 1 + nThreadsdd[f]) /	THD.x);
		nThreadspx[f] = (2 * hposn0[hf] + 1) * (2 * hposn0[hf] + 1);
		BLKpx[f].x = floor((THD.x -1 + nThreadspx[f]) / THD.x);
		/* Increase f, check for bounds, calculate odd frame */
		f++; if (f==nfrm_alloc)	break;
		nThreadsdd[f] = hndel1[hf]*hndop1[hf];
		BLKdd[f].x = floor((THD.x - 1 + nThreadsdd[f]) /	THD.x);
		nThreadspx[f] = (2 * hposn1[hf] + 1) * (2 * hposn1[hf] + 1);
		BLKpx[f].x = floor((THD.x -1 + nThreadspx[f]) / THD.x);
		hf++;
		/* If smearing is being modeled, initialize variables that
		 * will be used to sum results calculated for individual views.  */
		if (nviews > 1) {
			/* Allocate fit_store as a single pointer, originally a double
			 * pointer. This also initializes the entire array to zero. */
			gpuErrchk(cudaMalloc((void**)&fit_store[f], sizeof(float) * nThreadsdd[f]));
		}
	}
	/*  Loop over all views for this (smeared) frame, going in an order that
        ends with the view corresponding to the epoch listed for this frame
        in the obs file; this way we can use the calculated information for
        that view in the "write" action screen and disk output that follows   */
	for (v2=v0_index+1; v2<=v0_index+nviews; v2++)
		v[v2] = v2 % nviews;

	/* Set the ae and oe arrays for all frames on both GPUs */
	for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
		gpuErrchk(cudaSetDevice(GPU0));
		/* Launch nfrm_half0/nfrm_half1-threaded kernel to set pos->ae,pos->oe,pos->bistatic.*/
		set_aeoe_calcfits_krnl<<<BLK_half0,THD64,0,gpu0_stream[0]>>>(ddat,
				pos0, s,type, v[v2], nfrm_half0, 0, 0);
		gpuErrchk(cudaSetDevice(GPU1));
		set_aeoe_calcfits_krnl<<<BLK_half1,THD64,0,gpu1_stream[0]>>>(ddat,
				pos1, s, type, v[v2], nfrm_half1, 1, 0);
	}
	/* Initialize POS view for all frames on both GPUs by clearing out pos arrays */
	hf = 0;
	for (f=0; f<nfrm_alloc; f++) {
		for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
			gpuErrchk(cudaSetDevice(GPU0));
			posclr_mgpu_krnl<<<BLKpx[f],THD,0,gpu0_stream[hf]>>>(pos0, posn0, f, hf, 0);
			gpuErrchk(cudaSetDevice(GPU1));
			f++;	if (f>=nfrm_alloc)	break;
			posclr_mgpu_krnl<<<BLKpx[f],THD,0,gpu1_stream[hf]>>>(pos1, posn1, f, hf, 0);
			hf++;
		}
	} checkErrorAfterKernelLaunch("posclr_mgpu_krnl (calc_deldop_mgpu)");

	/* Call posvis_cuda_2 to get facet number, scattering angle, distance
	 * toward Earth at center of each POS pixel; set flag posbnd if any model
	 * portion extends beyond POS frame limits.*/
	/* NOTE: Limited to single component for now */
	for (v2=v0_index+1; v2<=v0_index+nviews; v2++)
		posvis_mgpu(dpar, dmod, ddat, pos0, pos1, orbit_off3,
				hposn0, hposn1, outbndarr0, outbndarr1, s, nfrm_alloc,
				0, nf, 0, c, type, gpu0_stream, gpu1_stream);
//
//	/* Now copy the out of bounds arrays back to the host */
//	gpuErrchk(cudaSetDevice(GPU0));
//	gpuErrchk(cudaMemcpy(&houtbndarr0, outbndarr0, sizeof(int)*nfrm_half0,
//			cudaMemcpyDeviceToHost));
//	gpuErrchk(cudaSetDevice(GPU1));
//	gpuErrchk(cudaMemcpy(&houtbndarr1, outbndarr1, sizeof(int)*nfrm_half1,
//			cudaMemcpyDeviceToHost));
//	gpuErrchk(cudaSetDevice(GPU0));

	/* The following code block also retries xylim values from the device and
	 * the dpar->exclude_seen flag */
	for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
		if (v[v2] == v0_index) {
			/* Set posbnd factors for frames if flag set by posvis */
			gpuErrchk(cudaSetDevice(GPU0));
			set_posbnd_calcfits_krnl<<<BLK_half0, THD64,0,gpu0_stream[0]>>>(dpar,
					ddat, pos0, xylim0, outbndarr0, s, type, nfrm_half0, 0, 0);
			gpuErrchk(cudaSetDevice(GPU1));
			set_posbnd_calcfits_krnl<<<BLK_half1, THD64,0,gpu1_stream[0]>>>(dpar,
					ddat, pos1, xylim1, outbndarr1, s, type, nfrm_half1, 1, 0);
		}
	} checkErrorAfterKernelLaunch("set_posbnd_krnl (calc_deldop_mgpu");

	/* Now copy xylim and exclude_seen back to host */
	gpuErrchk(cudaSetDevice(GPU0));
	gpuErrchk(cudaMemcpyFromSymbol(&exclude_seen, mgpu_exclude_seen,
			sizeof(int), 0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(&hxylim0, xylim0, sizeof(int4)*nfrm_half0,
			cudaMemcpyDeviceToHost));
	gpuErrchk(cudaSetDevice(GPU1));
	gpuErrchk(cudaMemcpy(&hxylim1, xylim1, sizeof(int4)*nfrm_half1,
			cudaMemcpyDeviceToHost));
	gpuErrchk(cudaSetDevice(GPU0));

	/* Calculate launch parameters for all frames */
	hf = 0;
	for (f=0; f<nfrm_alloc; f++) {
		xspan[f] = hxylim0[hf].x - hxylim0[hf].w + 1;
		yspan = hxylim0[hf].z - hxylim0[hf].y + 1;
		nThreadspx1[f] = xspan[f] * yspan;
		BLKpx[f].x = (THD.x -1 + nThreadspx1[f]) / THD.x;

		f++;	if (f>=nfrm_alloc)	break;
		xspan[f] = hxylim1[hf].x - hxylim1[hf].w + 1;
		yspan = hxylim1[hf].z - hxylim1[hf].y + 1;
		nThreadspx1[f] = xspan[f] * yspan;
		BLKpx[f].x = (THD.x -1 + nThreadspx1[f]) / THD.x;
		hf++;
	}

	hf = 0;
	for (f=0; f<nfrm_alloc; f++) {
		for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
			gpuErrchk(cudaSetDevice(GPU0));
			/* Go through all POS pixels which are visible with low enough
			 * scattering angle and mark the facets which project onto their
			 * centers as having been "seen" at least once                   */
			if (s != exclude_seen && v[v2] == v0_index)
				cf_mark_pixels_seen_krnl<<<BLKpx[f],THD,0,gpu0_stream[hf]>>>(
						dpar, dmod, pos0, xylim0, nThreadspx1[f], xspan[f], hf);

			/* Clear out fit arrays in preparation for pos2deldop calculation */
			clrvect_krnl<<<BLKdd[f],THD,0,gpu0_stream[hf+1]>>>(ddat, nThreadsdd[f], s, f);

			/* Now switch to GPU1 and repeat */
			f++;	if (f>=nfrm_alloc)	break;
			gpuErrchk(cudaSetDevice(GPU1));
			if (s != exclude_seen && v[v2] == v0_index)
				cf_mark_pixels_seen_krnl<<<BLKpx[f],THD,0,gpu1_stream[hf]>>>(
						dpar, dmod, pos1, xylim1, nThreadspx1[f], xspan[f], hf);
			clrvect_krnl<<<BLKdd[f],THD,0,gpu1_stream[hf+1]>>>(ddat, nThreadsdd[f], s, f);
			hf++;
		} checkErrorAfterKernelLaunch("cf_mark_pixels_seen and clrvct_krnl in "
				"calc_deldop_mgpu");
	}

	/* Call the delay-Doppler radar calculation for all frames, loop through views */
	for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
		pos2deldop_mgpu(dpar, dmod, ddat, pos0, pos1, frame0, frame1,
				xylim0, xylim1, ndel0, ndel1, ndop0, ndop1, 0.0, 0.0, 0.0,
				0, s, nfrm_alloc, v[v2], outbndarr0, outbndarr1,
				gpu0_stream, gpu1_stream);
	}
//	/* Now copy the out of bounds arrays back to the host */
//	gpuErrchk(cudaSetDevice(GPU0));
//	gpuErrchk(cudaMemcpy(&houtbndarr0, outbndarr0, sizeof(int)*nfrm_half0,
//			cudaMemcpyDeviceToHost));
//	gpuErrchk(cudaSetDevice(GPU1));
//	gpuErrchk(cudaMemcpy(&houtbndarr1, outbndarr1, sizeof(int)*nfrm_half1,
//			cudaMemcpyDeviceToHost));
//	gpuErrchk(cudaSetDevice(GPU0));

	/* The following checks the bad radar flags for each frame and sets logfactors
	 * accordingly */
	for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
		set_badradar_calcfits_krnl<<<1,1>>>(dpar, ddat, s, nfrm_alloc,
				nfrm_half0, nfrm_half1, outbndarr0, outbndarr1, type);
		checkErrorAfterKernelLaunch("set_badradar_calcfits_krnl (calc_fits_mgpu)");
	}

	for (f=0; f<nfrm_alloc; f++) {
		/* If smearing is being modeled, include delay-Doppler calculations
		 * from this view in the summed results for this frame  */
		if (nviews > 1) {
			/* add fit_s[offset] to fit_store[i][offset]*/

			cf_add_fit_store_krnl1<<<BLKdd[f],THD,0,gpu0_stream[f]>>>(
					ddat,fit_store,nThreadsdd[f],s,f, type);
			cf_add_fit_store_krnl2<<<1,1>>>(ddat, s, f, overflow, type);
		}
	} checkErrorAfterKernelLaunch("cf_add_fit_store_krnl1 and 2 (calc_deldop_mgpu");

	/* If smearing is being modeled, compute mean values over all views for
	 * this frame and store them in the standard frame structure     */
	/* This kernel also carries out the gamma transformation on the fit
	 * image if the par->dd_gamma flag is not set  */
	if (nviews > 1) {
		BLKpx[0] = floor((THD64.x -1 + nfrm_alloc)/THD64.x);
		hf = 0;
		for (f=0; f<nfrm_alloc; f++) {
			gpuErrchk(cudaSetDevice(GPU0));
			cf_finish_fit_store_krnl<<<BLKdd[f],THD,0,gpu0_stream[hf]>>>(
					ddat, fit_store, s, f, nThreadsdd[f], type);
			gammatrans_calcfits_krnl<<<BLKdd[f],THD,0,gpu0_stream[hf]>>>(dpar,
					ddat, s, f, nThreadsdd[f], type);
			f++;	if (f>=nfrm_alloc)	break;
			gpuErrchk(cudaSetDevice(GPU1));
			cf_finish_fit_store_krnl<<<BLKdd[f],THD,0,gpu1_stream[hf]>>>(
					ddat, fit_store, s, f, nThreadsdd[f], type);
			gammatrans_calcfits_krnl<<<BLKdd[f],THD,0,gpu1_stream[hf]>>>(dpar,
								ddat, s, f, nThreadsdd[f], type);
			hf++;
		}checkErrorAfterKernelLaunch("smearing kernels (calc_deldop_mgpu");

		/* This kernel handles all frames, all on GPU0 */
		gpuErrchk(cudaSetDevice(GPU0));
		finish_overflow_krnl<<<BLKpx[0],THD64,0,gpu0_stream[0]>>>(ddat, overflow,
				s, type, nfrm_alloc);
	} checkErrorAfterKernelLaunch("finish_overflow_krnl (calc_deldop_mgpu");
}

__host__ void calc_doppler_mgpu(
		struct par_t *dpar,
		struct mod_t *dmod,
		struct dat_t *ddat,
		int s,
		int nfrm_alloc, int nfrm_half0, int nfrm_half1,
		dim3 BLK_half0, dim3 BLK_half1,
		int nviews,
		unsigned char type,
		int nf,
		struct pos_t **pos0, 		struct pos_t **pos1,
		struct dopfrm_t **frame0, 	struct dopfrm_t **frame1,
		struct dopview_t **view0_0, struct dopview_t **view0_1,
		float **fit_store,
		int	*ndop0, 				int *ndop1,
		int *posn0, 				int *posn1,
		int4 *xylim0, 				int4 *xylim1,
		float *overflow,
		int	*outbndarr0,			int *outbndarr1,
		cudaStream_t *gpu0_stream,	cudaStream_t *gpu1_stream)
{
	float3 orbit_off3;
	orbit_off3.x = orbit_off3.y = orbit_off3.z = 0.0;
	int v0_index, exclude_seen, f, hf, v2, c=0, yspan;
	dim3 BLKd[nfrm_alloc], BLKpx[nfrm_alloc],THD,THD64;
	THD.x = maxThreadsPerBlock; THD64.x = 64;
	int4 hxylim0[nfrm_half0], hxylim1[nfrm_half1];
	int hndop0[nfrm_half0], hndop1[nfrm_half1], hndop[nfrm_alloc],
		hposn0[nfrm_half0],	hposn1[nfrm_half1], xspan[nfrm_alloc],
		nThreadspx[nfrm_alloc],	nThreadspx1[nfrm_alloc], v[nviews+1];

	/* Launch the initialization kernel on each GPU and copy results back to
	 * split host arrays */
	gpuErrchk(cudaSetDevice(GPU0));
	init_doppler_calcfits_krnl<<<BLK_half0,THD64,0,gpu0_stream[0]>>>(ddat,frame0,
			pos0, view0_0, ndop0, posn0, s, nfrm_half0, 0);
	gpuErrchk(cudaMemcpyAsync(hndop0, ndop0, sizeof(int)*nfrm_half0,
			cudaMemcpyDeviceToHost, gpu0_stream[0]));
	gpuErrchk(cudaMemcpyAsync(hposn0, posn0, sizeof(int)*nfrm_half0,
			cudaMemcpyDeviceToHost, gpu0_stream[0]));

	gpuErrchk(cudaSetDevice(GPU1));
	init_doppler_calcfits_krnl<<<BLK_half1,THD64,0,gpu1_stream[0]>>>(ddat,frame1,
			pos1, view0_1, ndop1, posn1, s, nfrm_half1, 1);
	gpuErrchk(cudaMemcpyAsync(hndop1, ndop1, sizeof(int)*nfrm_half1,
			cudaMemcpyDeviceToHost, gpu1_stream[0]));
	gpuErrchk(cudaMemcpyAsync(hposn1, posn1, sizeof(int)*nfrm_half1,
			cudaMemcpyDeviceToHost, gpu1_stream[0]));
	checkErrorAfterKernelLaunch("init_doppler_calcfits_krnl");
	gpuErrchk(cudaSetDevice(GPU0));
	gpuErrchk(cudaMemcpyFromSymbol(&v0_index, mgpu_v0_index, sizeof(int),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemsetAsync(overflow, 0, sizeof(float)*6, gpu0_stream[1]));

	/* Determine launch parameters of kernels to come */
	hf = 0;
	for (f=0; f<nfrm_alloc; f++) {
		/* Calculate launch parameters needed later */
		BLKd[f] = floor((THD.x - 1 + hndop0[hf]) / THD.x);
		nThreadspx[f] = (2 * hposn0[hf] + 1) * (2 * hposn0[hf] + 1);
		BLKpx[f] = floor((THD.x - 1 + nThreadspx[f]) / THD.x);
		hndop[f] = hndop0[hf];
		f++;	if (f>=nfrm_alloc)	break;
		BLKd[f] = floor((THD.x - 1 + hndop1[hf]) / THD.x);
		nThreadspx[f] = (2 * hposn1[hf] + 1) * (2 * hposn1[hf] + 1);
		BLKpx[f] = floor((THD.x - 1 + nThreadspx[f]) / THD.x);
		hndop[f] = hndop1[hf];
		hf++;
		/* If smearing is being modeled, initialize variables that
		 * will be used to sum results calculated for individual views.  */
		if (nviews > 1)
			/* Allocate fit_store as a single pointer, originally a double
			 * pointer. This also initializes the entire array to zero. */
			gpuErrchk(cudaMalloc((void**)&fit_store[f], sizeof(float)*hndop[f]));
	}

	/* Loop over all views for this (smeared) frame, going in an order that
	 * ends with the view corresponding to the epoch listed for this frame
	 * in the obs file; this way we can use the calculated information for
	 * that view in the "write" action screen and disk output that follows*/
	for (v2=v0_index+1; v2<=v0_index+nviews; v2++)
		v[v2] = v2 % nviews;

	/* Set the ae and oe arrays for all frames on both GPUs */
	for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
		gpuErrchk(cudaSetDevice(GPU0));
		/* Launch nfrm_half0/nfrm_half1-threaded kernel to set pos->ae,pos->oe,pos->bistatic.*/
		set_aeoe_calcfits_krnl<<<BLK_half0,THD64,0,gpu0_stream[0]>>>(ddat,
				pos0, s, type, v[v2], nfrm_half0, 0, 0);
		gpuErrchk(cudaSetDevice(GPU1));
		set_aeoe_calcfits_krnl<<<BLK_half1,THD64,0,gpu1_stream[0]>>>(ddat,
				pos1, s, type, v[v2], nfrm_half1, 1, 0);
	}
	/* Initialize POS view for all frames on both GPUs by clearing out pos arrays */
	hf = 0;
	for (f=0; f<nfrm_alloc; f++) {
		for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
			gpuErrchk(cudaSetDevice(GPU0));
			posclr_mgpu_krnl<<<BLKpx[f],THD,0,gpu0_stream[hf]>>>(pos0, posn0, f, hf, 0);
			gpuErrchk(cudaSetDevice(GPU1));
			f++;	if (f>=nfrm_alloc)	break;
			posclr_mgpu_krnl<<<BLKpx[f],THD,0,gpu1_stream[hf]>>>(pos1, posn1, f, hf, 0);
			hf++;
		}
	} checkErrorAfterKernelLaunch("posclr_mgpu_krnl (calc_doppler_mgpu)");

	/* Call posvis_cuda_2 to get facet number, scattering angle, distance
	 * toward Earth at center of each POS pixel; set flag posbnd if any model
	 * portion extends beyond POS frame limits.*/
	/* NOTE: Limited to single component for now */
	for (v2=v0_index+1; v2<=v0_index+nviews; v2++)
		posvis_mgpu(dpar, dmod, ddat, pos0, pos1, orbit_off3,
				hposn0, hposn1, outbndarr0, outbndarr1, s, nfrm_alloc,
				0, nf, 0, c, type, gpu0_stream, gpu1_stream);

	/* Set posbnd factors (both GPUs) and also retrieve xylim values  */
	for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
		if (v[v2] == v0_index) {
			gpuErrchk(cudaSetDevice(GPU0));
			set_posbnd_calcfits_krnl<<<BLK_half0, THD64,0,gpu0_stream[0]>>>(dpar,
					ddat, pos0, xylim0, outbndarr0, s, type, nfrm_half0, 0, 0);
			gpuErrchk(cudaSetDevice(GPU1));
			set_posbnd_calcfits_krnl<<<BLK_half1, THD64,0,gpu1_stream[0]>>>(dpar,
					ddat, pos1, xylim1, outbndarr1, s, type, nfrm_half1, 1, 0);
		}
	} checkErrorAfterKernelLaunch("set_posbnd_krnl (calc_doppler_mgpu");

	/* Now copy xylim and exclude_seen back to host */
	gpuErrchk(cudaSetDevice(GPU0));
	gpuErrchk(cudaMemcpyFromSymbol(&exclude_seen, mgpu_exclude_seen,
			sizeof(int), 0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(&hxylim0, xylim0, sizeof(int4)*nfrm_half0,
			cudaMemcpyDeviceToHost));
	gpuErrchk(cudaSetDevice(GPU1));
	gpuErrchk(cudaMemcpy(&hxylim1, xylim1, sizeof(int4)*nfrm_half1,
			cudaMemcpyDeviceToHost));
	gpuErrchk(cudaSetDevice(GPU0));

	/* Calculate launch parameters for all frames */
	hf=0;
	for (f=0; f<nfrm_alloc; f++) {
		xspan[f] = hxylim0[hf].x - hxylim0[hf].w + 1;
		yspan = hxylim0[hf].z - hxylim0[hf].y + 1;
		nThreadspx1[f] = xspan[f] * yspan;
		BLKpx[f].x = (THD.x -1 + nThreadspx1[f]) / THD.x;

		f++;	if (f>=nfrm_alloc)	break;
		xspan[f] = hxylim1[hf].x - hxylim1[hf].w + 1;
		yspan = hxylim1[hf].z - hxylim1[hf].y + 1;
		nThreadspx1[f] = xspan[f] * yspan;
		BLKpx[f].x = (THD.x -1 + nThreadspx1[f]) / THD.x;
		hf++;
	}

	hf = 0;
	for (f=0; f<nfrm_alloc; f++) {
		for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
			gpuErrchk(cudaSetDevice(GPU0));
			/* Go through all POS pixels which are visible with low enough
			 * scattering angle and mark the facets which project onto their
			 * centers as having been "seen" at least once                   */
			if (s != exclude_seen && v[v2] == v0_index)
				cf_mark_pixels_seen_krnl<<<BLKpx[f],THD,0,gpu0_stream[hf]>>>(
						dpar, dmod, pos0, xylim0, nThreadspx1[f], xspan[f], hf);

			/* Clear out fit arrays in preparation for pos2deldop calculation */
			clrvect_krnl<<<BLKd[f],THD,0,gpu0_stream[hf+1]>>>(ddat, hndop0[hf], s, f);

			/* Now switch to GPU1 and repeat */
			f++;	if (f>=nfrm_alloc)	break;
			gpuErrchk(cudaSetDevice(GPU1));
			if (s != exclude_seen && v[v2] == v0_index)
				cf_mark_pixels_seen_krnl<<<BLKpx[f],THD,0,gpu1_stream[hf]>>>(
						dpar, dmod, pos1, xylim1, nThreadspx1[f], xspan[f], hf);
			clrvect_krnl<<<BLKd[f],THD,0,gpu1_stream[hf+1]>>>(ddat, hndop1[hf], s, f);
			hf++;
		} checkErrorAfterKernelLaunch("cf_mark_pixels_seen and clrvct_krnl in "
				"calc_doppler_mgpu");
	}

	/* Call the delay-Doppler radar calculation for all frames, loop through views */
	for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
		pos2doppler_mgpu(dpar, dmod, ddat, pos0, pos1, frame0, frame1, xylim0,
				xylim1, 0.0,0.0,0.0, ndop0, ndop1, 0, s, nfrm_alloc, v[v2],
				outbndarr0, outbndarr1,	gpu0_stream, gpu1_stream);
	}
	/* The following checks the bad radar flags for each frame and sets logfactors
	 * accordingly */
	for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
		set_badradar_calcfits_krnl<<<1,1>>>(dpar, ddat, s, nfrm_alloc,
				nfrm_half0, nfrm_half1, outbndarr0, outbndarr1, type);
		checkErrorAfterKernelLaunch("set_badradar_calcfits_krnl (calc_doppler_mgpu)");
	}
	for (f=0; f<nfrm_alloc; f++) {
		/* If smearing is being modeled, include delay-Doppler calculations
		 * from this view in the summed results for this frame  */
		if (nviews > 1) {
			/* add fit_s[offset] to fit_store[i][offset]*/

			cf_add_fit_store_krnl1<<<BLKd[f],THD,0,gpu0_stream[f]>>>(
					ddat,fit_store,hndop[f],s,f, type);
			cf_add_fit_store_krnl2<<<1,1>>>(ddat, s, f, overflow, type);
		}
	} checkErrorAfterKernelLaunch("cf_add_fit_store_krnl1 and 2 (calc_doppler_mgpu");

	/* If smearing is being modeled, compute mean values over all views for
	 * this frame and store them in the standard frame structure     */
	/* This kernel also carries out the gamma transformation on the fit
	 * image if the par->dd_gamma flag is not set  */
	if (nviews > 1) {
		BLKpx[0] = floor((THD64.x -1 + nfrm_alloc)/THD64.x);
		hf = 0;
		for (f=0; f<nfrm_alloc; f++) {
			gpuErrchk(cudaSetDevice(GPU0));
			cf_finish_fit_store_krnl<<<BLKd[f],THD,0,gpu0_stream[hf]>>>(ddat,
					fit_store, s, f, hndop[f], type);
			gammatrans_calcfits_krnl<<<BLKd[f],THD,0,gpu0_stream[hf]>>>(dpar,
					ddat, s, f, hndop[f], type);
			f++;	if (f>=nfrm_alloc)	break;
			gpuErrchk(cudaSetDevice(GPU1));
			cf_finish_fit_store_krnl<<<BLKd[f],THD,0,gpu1_stream[hf]>>>(ddat,
					fit_store, s, f, hndop[f], type);
			gammatrans_calcfits_krnl<<<BLKd[f],THD,0,gpu1_stream[hf]>>>(dpar,
					ddat, s, f, hndop[f], type);
			hf++;
		}checkErrorAfterKernelLaunch("smearing kernels (calc_deldop_mgpu");

		/* This kernel handles all frames, all on GPU0 */
		gpuErrchk(cudaSetDevice(GPU0));
		finish_overflow_krnl<<<BLKpx[0],THD64,0,gpu0_stream[0]>>>(ddat, overflow,
				s, type, nfrm_alloc);
	} checkErrorAfterKernelLaunch("finish_overflow_krnl (calc_doppler_mgpu");
}

__host__ void calc_lghtcrv_mgpu(
		struct par_t *dpar,
		struct mod_t *dmod,
		struct dat_t *ddat,
		int s,
		int nfrm_alloc, int nfrm_half0, int nfrm_half1,
		dim3 BLK_half0, dim3 BLK_half1,
		int nviews,
		unsigned char type,
		int lc_n,
		int nf,
		struct pos_t **pos0,		struct pos_t **pos1,
		struct crvrend_t **rend0,	struct crvrend_t **rend1,
		int *posn0,					int *posn1,
		/*int *bistatic,*/
		float *pixels_per_km0,		float *pixels_per_km1,
		float *overflow,
		int4 *xylim0,				int4 *xylim1,
		int *outbndarr0,			int *outbndarr1,
		double3 *so0,				double3 *so1,
		double *u,
		cudaStream_t *gpu0_stream,
		cudaStream_t *gpu1_stream)
{
	int ncalc, c=0, n, nThreads, exclude_seen, f, hf, bistatic_all0,
			bistatic_all1;
	float3 orbit_off3;
	orbit_off3.x = orbit_off3.y = orbit_off3.z = 0.0;
	dim3 BLKpx[nfrm_alloc],THD,THD64;
	THD.x = maxThreadsPerBlock; THD64.x = 64;
	ncalc = nfrm_alloc;
	n = lc_n;
	int nThreadspx[nfrm_alloc], nThreadspx1[nfrm_alloc], hposn0[nfrm_half0],
		hposn1[nfrm_half1];
	int4 hxylim0[nfrm_half0], hxylim1[nfrm_half1];
	int2 span[nfrm_alloc];

	printf("calc_lghtcrv_mgpu in calc_fits_mgpu\n");

	/* Launch the initialization kernel on each GPU and copy results back to
	 * split host arrays */
	gpuErrchk(cudaSetDevice(GPU0));
	init_lghtcrv_calcfits_krnl<<<BLK_half0,THD64,0,gpu0_stream[0]>>>(ddat,rend0,
			pos0, posn0, s, nfrm_half0, 0);
	gpuErrchk(cudaMemcpyAsync(hposn0, posn0, sizeof(int)*nfrm_half0,
			cudaMemcpyDeviceToHost, gpu0_stream[0]));
	gpuErrchk(cudaSetDevice(GPU1));
	init_lghtcrv_calcfits_krnl<<<BLK_half1,THD64,0,gpu1_stream[0]>>>(ddat, rend1,
			pos1, posn1, s, nfrm_half1, 1);
	gpuErrchk(cudaMemcpyAsync(hposn1, posn1, sizeof(int)*nfrm_half1,
			cudaMemcpyDeviceToHost, gpu1_stream[0]));
	checkErrorAfterKernelLaunch("init_doppler_calcfits_krnl");
	gpuErrchk(cudaSetDevice(GPU0));
	gpuErrchk(cudaMemsetAsync(overflow, 0, sizeof(float)*6, gpu0_stream[1]));

	/* Calculate launch parameters for later */
	hf = 0;
	for (f=1; f<=nfrm_alloc; f++) {
		span[f].x = (2 * hposn0[hf] + 1);
		nThreadspx[f] =  span[f].x * span[f].x;
		BLKpx[f] = floor((THD.x - 1 + nThreadspx[f]) / THD.x);
		f++;	if (f>=nfrm_alloc)	break;
		span[f].x = (2 * hposn1[hf] + 1);
		nThreadspx[f] =  span[f].x * span[f].x;
		BLKpx[f] = floor((THD.x - 1 + nThreadspx[f]) / THD.x);
		hf++;
	}

	/* Set the ae and oe arrays for all frames on both GPUs */
	gpuErrchk(cudaSetDevice(GPU0));
	/* Launch nfrm_half0/nfrm_half1-threaded kernel to set pos->ae,pos->oe,pos->bistatic.*/
	set_aeoe_calcfits_krnl<<<BLK_half0,THD64,0,gpu0_stream[0]>>>(ddat,
			pos0, s, type, 0, nfrm_half0, 0, 1);
	gpuErrchk(cudaSetDevice(GPU1));
	set_aeoe_calcfits_krnl<<<BLK_half1,THD64,0,gpu1_stream[0]>>>(ddat,
			pos1, s, type, 0, nfrm_half1, 1, 1);
	checkErrorAfterKernelLaunch("set_aeoe_calcfits_krnl (calc_lghtcrv_mgpu");

	/* Initialize POS view for all frames on both GPUs by clearing out pos arrays */
	hf = 0;
	for (f=1; f<nfrm_alloc; f++) {
		gpuErrchk(cudaSetDevice(GPU0));
		posclr_mgpu_krnl<<<BLKpx[f],THD,0,gpu0_stream[hf]>>>(pos0, posn0, f, hf, 0);
		gpuErrchk(cudaSetDevice(GPU1));
		f++;	if (f>=nfrm_alloc)	break;
		posclr_mgpu_krnl<<<BLKpx[f],THD,0,gpu1_stream[hf]>>>(pos1, posn1, f, hf, 0);
		hf++;
	} checkErrorAfterKernelLaunch("posclr_mgpu_krnl (calc_lghtcrv_mgpu)");

	/* Call routine posvis to get  facet number, scattering & incidence
	 * angle, distance toward Earth at center of each POS pixel; set posbnd
	 * parameter = 1 if any model portion extends beyond POS frame limits*/
	posvis_mgpu(dpar, dmod, ddat, pos0, pos1, orbit_off3,
			hposn0, hposn1, outbndarr0, outbndarr1, s, nfrm_alloc,
			0, nf, 0, c, type, gpu0_stream, gpu1_stream);
	/* Set posbnd factors (both GPUs) and also retrieve xylim values  */
	gpuErrchk(cudaSetDevice(GPU0));
	set_posbnd_calcfits_krnl<<<BLK_half0, THD64,0,gpu0_stream[0]>>>(dpar,
			ddat, pos0, xylim0, outbndarr0, s, type, nfrm_half0, 0, 1);
	gpuErrchk(cudaSetDevice(GPU1));
	set_posbnd_calcfits_krnl<<<BLK_half1, THD64,0,gpu1_stream[0]>>>(dpar,
			ddat, pos1, xylim1, outbndarr1, s, type, nfrm_half1, 1, 1);

	checkErrorAfterKernelLaunch("set_posbnd_krnl (calc_lghtcrv_mgpu");

	/* Now copy xylim and exclude_seen back to host */
	gpuErrchk(cudaSetDevice(GPU0));
	gpuErrchk(cudaMemcpyFromSymbol(&exclude_seen, mgpu_exclude_seen,
			sizeof(int), 0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&bistatic_all0, lc_bistatic_all, sizeof(int),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(&hxylim0, xylim0, sizeof(int4)*nfrm_half0,
			cudaMemcpyDeviceToHost));
	gpuErrchk(cudaSetDevice(GPU1));
	gpuErrchk(cudaMemcpy(&hxylim1, xylim1, sizeof(int4)*nfrm_half1,
			cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&bistatic_all1, lc_bistatic_all, sizeof(int),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaSetDevice(GPU0));

	/* Now view model from source (sun) and get facet number and distance
	 * toward source of each pixel in this projected view; use this
	 * information to determine which POS pixels are shadowed       */
	/* Because posvis_gpu processes all frames at the same time, if
	 * any of the frames are bistatic, all of them get calculated again  */

	if (bistatic_all0 || bistatic_all1) {
		posvis_mgpu(dpar, dmod, ddat, pos0, pos1, orbit_off3,
				hposn0, hposn1, outbndarr0, outbndarr1, s, nfrm_alloc,
				1, nf, 0, c, type, gpu0_stream, gpu1_stream);

		/* Set posbnd factors (both GPUs) and also retrieve xylim values  */
		gpuErrchk(cudaSetDevice(GPU0));
		set_posbnd_calcfits_krnl<<<BLK_half0, THD64,0,gpu0_stream[0]>>>(dpar,
				ddat, pos0, xylim0, outbndarr0, s, type, nfrm_half0, 0, 1);
		gpuErrchk(cudaSetDevice(GPU1));
		set_posbnd_calcfits_krnl<<<BLK_half1, THD64,0,gpu1_stream[0]>>>(dpar,
				ddat, pos1, xylim1, outbndarr1, s, type, nfrm_half1, 1, 1);
		checkErrorAfterKernelLaunch("set_posbnd_krnl (calc_lghtcrv_mgpu"
				"bistatic calculation");
		gpuErrchk(cudaSetDevice(GPU0));

		/* Initialize plane-of-sky for all frames for the posmask kernel */
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
			posmask_krnl<<<BLKpx[f],THD,0,gpu0_stream[hf]>>>(
					dpar, pos0, so0, pixels_per_km0, posn0,
					nThreadspx[f], span[f].x, hf);
		} checkErrorAfterKernelLaunch("posmask_krnl for GPU0");

		gpuErrchk(cudaSetDevice(GPU1));
		for (hf=0; hf<nfrm_half1; hf++) {
			f = 2*hf+2;	/* Allows for lghtcrv offset */
			posmask_krnl<<<BLKpx[f],THD,0,gpu1_stream[hf]>>>(
					dpar, pos1, so1, pixels_per_km1, posn1,
					nThreadspx[f], span[f].x, hf);
		} checkErrorAfterKernelLaunch("posmask_krnl for GPU1");
	}

	/* Calculate launch parameters for all frames */
	hf = 0;
	for (f=1; f<nfrm_alloc; f++) {
		span[f].x = hxylim0[hf].x - hxylim0[hf].w + 1;
		span[f].y = hxylim0[hf].z - hxylim0[hf].y + 1;
		nThreadspx1[f] = span[f].x * span[f].y;
		BLKpx[f].x = floor ((THD.x -1 + nThreadspx1[f]) / THD.x);

		f++;	if (f>=nfrm_alloc)	break;
		span[f].x = hxylim1[hf].x - hxylim1[hf].w + 1;
		span[f].y = hxylim1[hf].z - hxylim1[hf].y + 1;
		nThreadspx1[f] = span[f].x * span[f].y;
		BLKpx[f].x = floor ((THD.x -1 + nThreadspx1[f]) / THD.x);
		hf++;
	}

	if (s != exclude_seen) {
		hf = 0;
		for (f=1; f<nfrm_alloc; f++) {
			gpuErrchk(cudaSetDevice(GPU0));
			/* Go through all POS pixels which are visible with low enough
			 * scattering angle and mark the facets which project onto their
			 * centers as having been "seen" at least once                   */
			cf_mark_pixels_seen_krnl<<<BLKpx[f],THD,0,gpu0_stream[hf]>>>(
					dpar, dmod, pos0, xylim0, nThreadspx1[f], span[f].x, hf);
			/* Now switch to GPU1 and repeat */
			f++;	if (f>=nfrm_alloc)	break;
			gpuErrchk(cudaSetDevice(GPU1));
			cf_mark_pixels_seen_krnl<<<BLKpx[f],THD,0,gpu1_stream[hf]>>>(
					dpar, dmod, pos1, xylim1, nThreadspx1[f], span[f].x, hf);
			hf++;
			checkErrorAfterKernelLaunch("cf_mark_pixels_seen and clrvct_krnl in "
					"calc_lghtcrv_mgpu");
		}
	}

	/* Compute model brightness for this lightcurve point then copy to device  */
//	apply_photo_mgpu_f(dmod, ddat, pos0, pos1, xylim0, xylim1, span,
//			BLKpx, nThreadspx1, 0, s, nfrm_alloc, nfrm_half0,
//			nfrm_half1, nThreadspx, gpu0_stream, gpu1_stream);

	/* Now that we have calculated the model lightcurve brightnesses y at each
	 * of the epochs x, we use cubic spline interpolation (Numerical Recipes
	 * routines spline and splint) to get model lightcurve brightness fit[i] at
	 * each OBSERVATION epoch t[i], with i=1,2,...,n. This will allow us (in
	 * routine chi2) to compare model to data (fit[i] to obs[i]) to get chi-
	 * square. Note that vector y2 contains the second derivatives of the
	 * interpolating function at the calculation epochs x. Smearing is handled
	 * by interpolating the brightness at the time t of each individual view
	 * and then taking the mean of all views that correspond to a given
	 * observed lightcurve point.                         */
	/* Configure and launch an ncalc-threaded kernel that performs what NR
	 * function spline does.  Original call:
	 *
	 * spline( lghtcrv->x, lghtcrv->y, ncalc, 2.0e30, 2.0e30, lghtcrv->y2);
	 */
//	printf("nfrm_alloc: %i\n", nfrm_alloc);
//	gpuErrchk(cudaSetDevice(GPU0));
//	/* First make a pointer for u and cudaMalloc device memory for it */
//	gpuErrchk(cudaMemset(u, 0, nfrm_alloc*sizeof(double)));
//
//	lghtcrv_spline_streams_test_krnl<<<1,1>>>(ddat, s, 2.0e30,
//			2.0e30, u, (nfrm_alloc-1));
//	checkErrorAfterKernelLaunch("lghtcrv_spline_streams_krnl");
//
//	/* Change launch parameters from ncalc threads to n threads */
//	lghtcrv_splint_streams3_test_krnl<<<1,1>>>(ddat, s, n, (nfrm_alloc-1));
//	checkErrorAfterKernelLaunch("lghtcrv_splint_streams_krnl");


//	nThreads = nframes;
//	THD.x = maxThreadsPerBlock;
//	gpuErrchk(cudaMemset(u, 0, sizeof(double)*nfplus));
//	BLKpx[0].x = floor((THD.x - 1 + nThreads)/THD.x);
//	lghtcrv_spline_streams_test_krnl<<<1,1>>>(ddat, s, 2.0e30, 2.0e30, u, nframes);
//	checkErrorAfterKernelLaunch("lghtcrv_spline_streams_krnl");
//
//	/* Change launch parameters from ncalc threads to n threads */
//	BLKpx[0].x = floor((THD.x - 1 + n) / THD.x);
//	lghtcrv_splint_streams3_test_krnl<<<1,1>>>(ddat, s, n, nframes);
//	checkErrorAfterKernelLaunch("lghtcrv_splint_streams_krnl");
}
