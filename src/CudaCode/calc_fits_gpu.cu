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
#include "../shape/head.h"
}

__host__ void calc_deldop_gpu32(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, struct vertices_t **verts, int s, int nframes,
		int nviews, unsigned char type, int nf, struct pos_t **pos, struct
		deldopfrm_t **frame, struct deldopview_t **view0, float **fit_store,
		int *ndel, int *ndop, int *posn, int4 *xylim, float
		*overflow, int *outbndarr, cudaStream_t *cf_stream);
__host__ void calc_deldop_gpu64(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, struct vertices_t **verts, int s, int nframes,
		int nviews, unsigned char type, int nf, struct pos_t **pos, struct
		deldopfrm_t **frame, struct deldopview_t **view0, double **fit_store,
		int *ndel, int *ndop, int *posn, int4 *xylim, double
		*overflow, int *outbndarr, cudaStream_t *cf_stream);
__host__ void calc_doppler_gpu32(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, struct vertices_t **verts, int s, int nframes, int
		nviews, unsigned char type, int nf, struct pos_t **pos, struct dopfrm_t
		**frame, struct dopview_t **view0, float **fit_store, int *ndop, int
		*posn, int4 *xylim, float *overflow, int	*outbndarr,
		cudaStream_t *cf_stream);
__host__ void calc_doppler_gpu64(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, struct vertices_t **verts, int s, int nframes, int
		nviews, unsigned char type, int nf, struct pos_t **pos, struct dopfrm_t
		**frame, struct dopview_t **view0, double **fit_store, int *ndop, int
		*posn, int4 *xylim, double *overflow, int	*outbndarr,
		cudaStream_t *cf_stream);
__host__ void calc_lghtcrv_gpu32(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, struct vertices_t **verts, int s, int nframes, int
		nviews, unsigned char type, int lc_n, int nf, struct pos_t **pos,
		struct crvrend_t **rend, int *posn, int *bistatic, float *pxlpkm, float
		*overflow, int4 *xylim, int *outbndarr, double3 *so, double *u,
		cudaStream_t *cf_stream);
__host__ void calc_lghtcrv_gpu64(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, struct vertices_t **verts, int s, int nframes, int
		nviews, unsigned char type, int lc_n, int nf, struct pos_t **pos,
		struct crvrend_t **rend, int *posn, int *bistatic, double *pxlpkm, double
		*overflow, int4 *xylim, int *outbndarr, double3 *so, double *u,
		cudaStream_t *cf_stream);
void *calc_fits_pthread_sub(void *ptr);

__device__ int cfs_v0_index, cfs_exclude_seen;
__device__ double cfs_lghtcrv_posbnd_logfactor;

typedef struct calcfits_thread_t
{
    int thread_no;
	struct par_t *parameter;
    struct mod_t *model;
    struct dat_t *data;
    struct vertices_t **verts;
    unsigned char *type;
    int *nviews;
    int *nframes;
    int *hlc_n;
    int *GPUID;
    int gpuid;
    int nsets;
    int nf;
    int max_frames;
    cudaStream_t *gpu_stream;
} calcfits_data;

__device__ void dev_splint_cfs(double *xa,double *ya,double *y2a,int n,double x,double *y)
{
	int klo,khi,k;
	double h,b,a;

	klo = 1;
	khi = n;
	while (khi-klo > 1) {
		k = (khi+klo) >> 1;
		if (xa[k] > x) 	khi=k;
		else klo = k;
	}
	h = xa[khi] - xa[klo];
	if (h == 0.0) 	printf("Bad XA input to routine SPLINT");
	a = (xa[khi] - x) / h;
	b = (x - xa[klo]) / h;
	*y = a * ya[klo] + b * ya[khi] + ((a*a*a-a) * y2a[klo] + (b*b*b-b) *
			y2a[khi]) * (h*h)/6.0;
}
__global__ void cf_init_seen_flags_krnl(struct mod_t *dmod, int nf, int c) {
	/* nf-threaded kernel */
	int f = blockIdx.x * blockDim.x + threadIdx.x;

	if (f < nf)
		dmod->shape.comp[c].real.f[f].seen = 0;
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
__global__ void cfs_init_devpar_krnl(struct par_t *dpar) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		dpar->posbnd = 0;
		dpar->badposet = 0;
		dpar->badradar = 0;
		dpar->posbnd_logfactor = 0.0;
		dpar->badposet_logfactor = 0.0;
		dpar->badradar_logfactor = 0.0;
	}
}

__host__ void calc_fits_gpu(
		struct par_t *dpar,
		struct mod_t *dmod,
		struct dat_t *ddat,
		struct vertices_t **verts,
		int *nviews,
		int *nframes,
		int *lc_n,
		unsigned char *type,
		int nsets,
		int nf,
		cudaStream_t *cf_stream,
		int max_frames)
{
	int s, nfrm_alloc, *ndel, *ndop, *posn, *bistatic, *outbndarr;
	int4 *xylim;
	float **fit_store32, *overflow32, *pxlpkm32;
	double **fit_store64, *overflow64, *pxlpkm64;
	double *u;
	double3 *so;
	nfrm_alloc = max_frames+1;

	struct pos_t **pos;
	struct deldopfrm_t **ddframe;
	struct dopfrm_t **dframe;
	struct crvrend_t **rend;
	struct deldopview_t **ddview0;
	struct dopview_t **dview0;
	dim3 BLK,THD;
	THD.x = maxThreadsPerBlock;

//	gpuErrchk(cudaSetDevice(GPU0));

	/* Allocate memory for all arrays that are needed for any possible data set.
	 * This is done to avoid repeat allocations/deallocations	 */
	gpuErrchk(cudaMalloc((void**)&pos, sizeof(pos_t*) * nfrm_alloc));
	gpuErrchk(cudaMalloc((void**)&ddframe, sizeof(deldopfrm_t*) * nfrm_alloc));
	gpuErrchk(cudaMalloc((void**)&ddview0, sizeof(deldopview_t*) * nfrm_alloc));
	gpuErrchk(cudaMalloc((void**)&dframe, sizeof(dopfrm_t*) * nfrm_alloc));
	gpuErrchk(cudaMalloc((void**)&dview0, sizeof(dopview_t*) * nfrm_alloc));
	gpuErrchk(cudaMalloc((void**)&rend, sizeof(crvrend_t*) * nfrm_alloc));
	gpuErrchk(cudaMalloc((void**)&ndel, sizeof(int) * nfrm_alloc));
	gpuErrchk(cudaMalloc((void**)&ndop, sizeof(int) * nfrm_alloc));
	gpuErrchk(cudaMalloc((void**)&posn, sizeof(int) * nfrm_alloc));
	gpuErrchk(cudaMalloc((void**)&bistatic, sizeof(int) * nfrm_alloc));
	gpuErrchk(cudaMalloc((void**)&xylim, sizeof(int4) * nfrm_alloc));
	gpuErrchk(cudaMalloc((void**)&outbndarr, sizeof(int) * nfrm_alloc));
	gpuErrchk(cudaMalloc((void**)&u, sizeof(double) * nfrm_alloc));
	gpuErrchk(cudaMalloc((void**)&so, 	 	sizeof(double3) *((nfrm_alloc*3)+1)));

	if (FP64) {
		gpuErrchk(cudaMalloc((void**)&fit_store64, sizeof(double*) * nfrm_alloc));
		gpuErrchk(cudaMalloc((void**)&overflow64, sizeof(double) * 6));
		gpuErrchk(cudaMalloc((void**)&pxlpkm64,   sizeof(double) * nfrm_alloc));
	} else {
		gpuErrchk(cudaMalloc((void**)&fit_store32, sizeof(float*) * nfrm_alloc));
		gpuErrchk(cudaMalloc((void**)&overflow32, sizeof(float) * 6));
		gpuErrchk(cudaMalloc((void**)&pxlpkm32,   sizeof(float) * nfrm_alloc));
	}

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
	int c = 0;
	BLK.x = floor((THD.x - 1 + nf)/THD.x);
	cf_init_seen_flags_krnl<<<BLK,THD>>>(dmod,nf,c);
	checkErrorAfterKernelLaunch("cf_init_seen_flags_krnl (calc_fits_cuda_streams)");

	/* Calculate the fits for each dataset in turn  */
	for (s=0; s<nsets; s++) {
		if (FP64) {
			switch (type[s]) {
			case DELAY:
				calc_deldop_gpu64(dpar, dmod, ddat, verts, s, nframes[s], nviews[s],
						type[s], nf, pos, ddframe, ddview0, fit_store64, ndel, ndop,
						posn, xylim, overflow64, outbndarr, cf_stream);
				break;
			case DOPPLER:
				calc_doppler_gpu64(dpar, dmod, ddat, verts, s, nframes[s],
						nviews[s], type[s], nf, pos, dframe, dview0, fit_store64,
						ndop, posn, xylim, overflow64, outbndarr,
						cf_stream );
				break;
			case POS:
				printf("Write calc_poset_cuda!");
				//			calc_poset_cuda(dpar, dmod, s);
				break;
			case LGHTCRV:
				calc_lghtcrv_gpu64(dpar, dmod, ddat, verts, s, nframes[s],
						nviews[s], type[s], lc_n[s],nf, pos, rend, posn, bistatic,
						pxlpkm64, overflow64, xylim, outbndarr, so, u, cf_stream);
				break;
			default:
				printf("calc_fits_gpu.c: can't handle this type yet\n");
			}
		}
		else {
			switch (type[s]) {
			case DELAY:
				calc_deldop_gpu32(dpar, dmod, ddat, verts, s, nframes[s], nviews[s],
						type[s], nf, pos, ddframe, ddview0, fit_store32, ndel, ndop,
						posn, xylim, overflow32, outbndarr, cf_stream);
				break;
			case DOPPLER:
				calc_doppler_gpu32(dpar, dmod, ddat, verts, s, nframes[s],
						nviews[s], type[s], nf, pos, dframe, dview0, fit_store32,
						ndop, posn, xylim, overflow32, outbndarr,
						cf_stream );
				break;
			case POS:
				printf("Write calc_poset_cuda!");
				//			calc_poset_cuda(dpar, dmod, s);
				break;
			case LGHTCRV:
				calc_lghtcrv_gpu32(dpar, dmod, ddat, verts, s, nframes[s],
						nviews[s], type[s], lc_n[s],nf, pos, rend, posn, bistatic,
						pxlpkm32, overflow32, xylim, outbndarr, so, u, cf_stream);
				break;
			default:
				printf("calc_fits_cuda.c: can't handle this type yet\n");
			}
		}
	}

	/* Complete calculations of values that will be used during a fit to
	 * increase the objective function for models with bad properties   */
	cf_set_final_pars_krnl<<<1,1>>>(dpar, ddat);
	checkErrorAfterKernelLaunch("cf_set_final_pars_krnl (calc_fits_cuda)");

	if (FP64) {
		cudaFree(fit_store64);
		cudaFree(overflow64);
		cudaFree(pxlpkm64);
	} else {
		cudaFree(fit_store32);
		cudaFree(overflow32);
		cudaFree(pxlpkm32);
	}
	cudaFree(pos);
	cudaFree(ddframe);
	cudaFree(dframe);
	cudaFree(rend);
	cudaFree(ddview0);
	cudaFree(dview0);
	cudaFree(ndel);
	cudaFree(ndop);
	cudaFree(posn);
	cudaFree(bistatic);
	cudaFree(xylim);
	cudaFree(outbndarr);
	cudaFree(u);
	cudaFree(so);
}

__host__ void calc_fits_pthreads(
		struct par_t *dpar0,
		struct par_t *dpar1,
		struct mod_t *dmod0,
		struct mod_t *dmod1,
		struct dat_t *ddat0,
		struct dat_t *ddat1,
		struct vertices_t **verts0,
		struct vertices_t **verts1,
		int *nviews,
		int *nframes,
		int *lc_n,
		int *GPUID,
		unsigned char *type,
		int nsets,
		int nf,
		int max_frames,
		pthread_t thread1,
		pthread_t thread2,
		cudaStream_t *gpu0_stream,
		cudaStream_t *gpu1_stream)
{

	dim3 BLK,THD;
	/* Create and populate the structures used to pass data to the pthread
	 * sub functions	 */
	calcfits_data data1, data2;
	data1.thread_no = 1;
	data2.thread_no = 2;
	data1.gpuid = GPU0;
	data2.gpuid = GPU1;
	data1.gpu_stream = gpu0_stream;
	data2.gpu_stream = gpu1_stream;
	data1.GPUID = data2.GPUID = GPUID;
	data1.data = ddat0;
	data2.data = ddat1;
	data1.parameter = dpar0;
	data2.parameter = dpar1;
	data1.model = dmod0;
	data2.model = dmod1;
	data1.hlc_n = data2.hlc_n = lc_n;
	data1.type = data2.type = type;
	data1.nf = data2.nf = nf;
	data1.nframes = data2.nframes = nframes;
	data1.nsets = data2.nsets = nsets;
	data1.verts = verts0;
	data2.verts = verts1;
	data1.nviews = data2.nviews = nviews;
	data1.max_frames = data2.max_frames = max_frames;

	/* Initialize flags that indicate the model extends beyond POS frame, that
	 * plane-of-sky fit images are too small to "contain" the target, and that
	 * model is too wide in (delay-)Doppler space to create (delay-)Doppler fit
	 * frames.  Note that this also gets mod->shape.nf and nsets            */

	cfs_init_devpar_krnl<<<1,1>>>(dpar0);
	checkErrorAfterKernelLaunch("cfs_init_devpar_krnl");
	gpuErrchk(cudaSetDevice(GPU1));
	cfs_init_devpar_krnl<<<1,1>>>(dpar1);
	checkErrorAfterKernelLaunch("cfs_init_devpar_krnl");
	gpuErrchk(cudaSetDevice(GPU0));

	/* From here, launch the pthreaded subfunction */
	pthread_create(&thread1, NULL, calc_fits_pthread_sub,(void*)&data1);
	pthread_create(&thread2, NULL, calc_fits_pthread_sub,(void*)&data2);

	pthread_join(thread1, NULL);
	pthread_join(thread2, NULL);

	/* Complete calculations of values that will be used during a fit to
	 * increase the objective function for models with bad properties   */
	gpuErrchk(cudaSetDevice(GPU0));
	cf_set_final_pars_krnl<<<1,1>>>(dpar0, ddat0);
	checkErrorAfterKernelLaunch("cf_set_final_pars_krnl on GPU0");
	gpuErrchk(cudaSetDevice(GPU1));
	cf_set_final_pars_krnl<<<1,1>>>(dpar1, ddat1);
	checkErrorAfterKernelLaunch("cf_set_final_pars_krnl on GPU1");
	gpuErrchk(cudaSetDevice(GPU0));

}

void *calc_fits_pthread_sub(void *ptr) {

	int c=0, s, nfrm_alloc, *ndel, *ndop, *posn, *bistatic, *outbndarr;
	int4 *xylim;
	float **fit_store32, *overflow32, *pxlpkm32;
	double *u, **fit_store64, *overflow64, *pxlpkm64;
	double3 *so;
	dim3 BLK,THD;
	struct pos_t **pos;
	struct deldopfrm_t **ddframe;
	struct dopfrm_t **dframe;
	struct crvrend_t **rend;
	struct deldopview_t **ddview0;
	struct dopview_t **dview0;
	calcfits_data *data;
	data = (calcfits_data *) ptr;
	nfrm_alloc = data->max_frames + 1;
	THD.x = maxThreadsPerBlock;
	gpuErrchk(cudaSetDevice(data->gpuid));

	/* Initialize the flags that indicate whether or not each facet of each
	 * model component is ever visible and unshadowed from Earth
	 * Note:  Single component only for now.  */
	BLK.x = floor((THD.x - 1 + data->nf)/THD.x);
	//for (c=0; c<mod->shape.ncomp; c++)
	cf_init_seen_flags_krnl<<<BLK,THD>>>(data->model, data->nf, c);
	checkErrorAfterKernelLaunch("cf_init_seen_flags_krnl");

	/* Allocate memory for all arrays that are needed for any possible data set.
	 * This is done to avoid repeat allocations/deallocations	 */
	gpuErrchk(cudaMalloc((void**)&pos, sizeof(pos_t*) * nfrm_alloc));
	gpuErrchk(cudaMalloc((void**)&ddframe, sizeof(deldopfrm_t*) * nfrm_alloc));
	gpuErrchk(cudaMalloc((void**)&ddview0, sizeof(deldopview_t*) * nfrm_alloc));
	gpuErrchk(cudaMalloc((void**)&dframe, sizeof(dopfrm_t*) * nfrm_alloc));
	gpuErrchk(cudaMalloc((void**)&dview0, sizeof(dopview_t*) * nfrm_alloc));
	gpuErrchk(cudaMalloc((void**)&rend, sizeof(crvrend_t*) * nfrm_alloc));
	gpuErrchk(cudaMalloc((void**)&ndel, sizeof(int) * nfrm_alloc));
	gpuErrchk(cudaMalloc((void**)&ndop, sizeof(int) * nfrm_alloc));
	gpuErrchk(cudaMalloc((void**)&posn, sizeof(int) * nfrm_alloc));
	gpuErrchk(cudaMalloc((void**)&bistatic, sizeof(int) * nfrm_alloc));
	gpuErrchk(cudaMalloc((void**)&xylim, sizeof(int4) * nfrm_alloc));
	gpuErrchk(cudaMalloc((void**)&outbndarr, sizeof(int) * nfrm_alloc));
	gpuErrchk(cudaMalloc((void**)&u, sizeof(double) * nfrm_alloc));
	gpuErrchk(cudaMalloc((void**)&so, 	 	sizeof(double3) *((nfrm_alloc*3)+1)));


	if (FP64) {
		gpuErrchk(cudaMalloc((void**)&fit_store64, sizeof(double*) * nfrm_alloc));
		gpuErrchk(cudaMalloc((void**)&overflow64, sizeof(double) * 6));
		gpuErrchk(cudaMalloc((void**)&pxlpkm64,   sizeof(double) * nfrm_alloc));
	} else {
		gpuErrchk(cudaMalloc((void**)&fit_store32, sizeof(float*) * nfrm_alloc));
		gpuErrchk(cudaMalloc((void**)&overflow32, sizeof(float) * 6));
		gpuErrchk(cudaMalloc((void**)&pxlpkm32,   sizeof(float) * nfrm_alloc));
	}

	/* Calculate the fits for each dataset in turn  */
	for (s=0; s<data->nsets; s++) {
		if (data->GPUID[s]==data->gpuid) {
			if (FP64) {
				switch (data->type[s]) {
				case DELAY:
					calc_deldop_gpu64(data->parameter, data->model, data->data,
							data->verts, s, data->nframes[s], data->nviews[s],
							data->type[s], data->nf, pos, ddframe, ddview0,
							fit_store64, ndel, ndop, posn, xylim, overflow64, outbndarr,
							data->gpu_stream);
					break;
				case DOPPLER:
					calc_doppler_gpu64(data->parameter, data->model, data->data,
							data->verts, s, data->nframes[s], data->nviews[s],
							data->type[s], data->nf, pos, dframe, dview0, fit_store64,
							ndop,posn,xylim,overflow64,outbndarr,data->gpu_stream);
					break;
				case POS:
					printf("Write calc_poset_cuda!");
					//			calc_poset_cuda(dpar, dmod, s);
					break;
				case LGHTCRV:
					calc_lghtcrv_gpu64(data->parameter, data->model, data->data,
							data->verts, s, data->nframes[s], data->nviews[s],
							data->type[s], data->hlc_n[s], data->nf, pos, rend,
							posn, bistatic, pxlpkm64, overflow64, xylim, outbndarr, so,
							u, data->gpu_stream);
					break;
				default:
					printf("calc_fits_pthreads_sub.c: can't handle this type yet\n");
				}
			}
			else {
				switch (data->type[s]) {
				case DELAY:
					calc_deldop_gpu32(data->parameter, data->model, data->data,
							data->verts, s, data->nframes[s], data->nviews[s],
							data->type[s], data->nf, pos, ddframe, ddview0,
							fit_store32, ndel, ndop, posn, xylim, overflow32, outbndarr,
							data->gpu_stream);
					break;
				case DOPPLER:
					calc_doppler_gpu32(data->parameter, data->model, data->data,
							data->verts, s, data->nframes[s], data->nviews[s],
							data->type[s], data->nf, pos, dframe, dview0, fit_store32,
							ndop,posn,xylim,overflow32,outbndarr,data->gpu_stream);
					break;
				case POS:
					printf("Write calc_poset_cuda!");
					//			calc_poset_cuda(dpar, dmod, s);
					break;
				case LGHTCRV:
					calc_lghtcrv_gpu32(data->parameter, data->model, data->data,
							data->verts, s, data->nframes[s], data->nviews[s],
							data->type[s], data->hlc_n[s], data->nf, pos, rend,
							posn, bistatic, pxlpkm32, overflow32, xylim, outbndarr, so,
							u, data->gpu_stream);
					break;
				default:
					printf("calc_fits_pthreads_sub.c: can't handle this type yet\n");
				}
			}
		}
	}

	cudaFree(pos);
	cudaFree(ddframe);
	cudaFree(dframe);
	cudaFree(rend);
	cudaFree(ddview0);
	cudaFree(dview0);
	cudaFree(ndel);
	cudaFree(ndop);
	cudaFree(posn);
	cudaFree(bistatic);
	cudaFree(xylim);
	cudaFree(outbndarr);
	cudaFree(u);
	cudaFree(so);
	if (FP64) {
		cudaFree(fit_store64);
		cudaFree(overflow64);
		cudaFree(pxlpkm64);
	} else {
		cudaFree(fit_store32);
		cudaFree(overflow32);
		cudaFree(pxlpkm32);
	}

	pthread_exit(0);
}

__global__ void cfs_set_deldop_shortcuts_krnl32(struct dat_t *ddat,
		struct deldopfrm_t **frame, struct pos_t **pos,
		struct deldopview_t **view0, int *ndop, int *ndel, float *overflow,
		int *posn, int s, int size) {
	/* nfrm_alloc-threaded kernel */
	/* 	overflow[0] - overflow_o2_store
	 * 	overflow[1] - overflow_m2_store
	 * 	overflow[2] - overflow_xsec_store
	 * 	overflow[3] - overflow_dopmean_store
	 * 	overflow[4] - overflow_delmean_store
	 */
	int f = blockIdx.x * blockDim.x + threadIdx.x;

	if (f < size) {
		frame[f] 	  = &ddat->set[s].desc.deldop.frame[f];
		ndop[f]  	  = frame[f]->ndop;
		ndel[f]	 	  = frame[f]->ndel;
		view0[f] 	  = &frame[f]->view[ddat->set[s].desc.deldop.v0];
		pos[f]		  = &frame[f]->pos;
		posn[f]		  = pos[f]->n;
	}
	__syncthreads();
	if (f==0) {
		cfs_v0_index  = ddat->set[s].desc.deldop.v0;
		overflow[0] = 0.0;
		overflow[1] = 0.0;
		overflow[2] = 0.0;
		overflow[3] = 0.0;
		overflow[4] = 0.0;
	}
}
__global__ void cfs_set_deldop_shortcuts_krnl64(struct dat_t *ddat,
		struct deldopfrm_t **frame, struct pos_t **pos,
		struct deldopview_t **view0, int *ndop, int *ndel, double *overflow,
		int *posn, int s, int size) {
	/* nfrm_alloc-threaded kernel */
	/* 	overflow[0] - overflow_o2_store
	 * 	overflow[1] - overflow_m2_store
	 * 	overflow[2] - overflow_xsec_store
	 * 	overflow[3] - overflow_dopmean_store
	 * 	overflow[4] - overflow_delmean_store
	 */
	int f = blockIdx.x * blockDim.x + threadIdx.x;

	if (f < size) {
		frame[f] 	  = &ddat->set[s].desc.deldop.frame[f];
		ndop[f]  	  = frame[f]->ndop;
		ndel[f]	 	  = frame[f]->ndel;
		view0[f] 	  = &frame[f]->view[ddat->set[s].desc.deldop.v0];
		pos[f]		  = &frame[f]->pos;
		posn[f]		  = pos[f]->n;
	}
	__syncthreads();
	if (f==0) {
		cfs_v0_index  = ddat->set[s].desc.deldop.v0;
		overflow[0] = 0.0;
		overflow[1] = 0.0;
		overflow[2] = 0.0;
		overflow[3] = 0.0;
		overflow[4] = 0.0;
	}
}
__global__ void cfs_set_doppler_shortcuts_krnl32(struct dat_t *ddat,
		struct dopfrm_t **frame, struct pos_t **pos, struct dopview_t **view0,
		float *overflow, int *ndop, int *posn, int s, int size) {
	/* nframes-threaded kernel */
	int f = blockIdx.x * blockDim.x + threadIdx.x;

	if (f < size) {

		frame[f] = &ddat->set[s].desc.doppler.frame[f];
		view0[f] = &frame[f]->view[ddat->set[s].desc.doppler.v0];
		ndop[f]	 = frame[f]->ndop;
		pos[f]	 = &frame[f]->pos;
		posn[f]	= pos[f]->n;
	}
	__syncthreads();

	if (f==0) {
		cfs_v0_index  = ddat->set[s].desc.doppler.v0;
		overflow[0] = 0.0;
		overflow[1] = 0.0;
		overflow[2] = 0.0;
		overflow[3] = 0.0;
		overflow[4] = 0.0;
	}
}
__global__ void cfs_set_doppler_shortcuts_krnl64(struct dat_t *ddat,
		struct dopfrm_t **frame, struct pos_t **pos, struct dopview_t **view0,
		double *overflow, int *ndop, int *posn, int s, int size) {
	/* nframes-threaded kernel */
	int f = blockIdx.x * blockDim.x + threadIdx.x;

	if (f < size) {

		frame[f] = &ddat->set[s].desc.doppler.frame[f];
		view0[f] = &frame[f]->view[ddat->set[s].desc.doppler.v0];
		ndop[f]	 = frame[f]->ndop;
		pos[f]	 = &frame[f]->pos;
		posn[f]	= pos[f]->n;
	}
	__syncthreads();

	if (f==0) {
		cfs_v0_index  = ddat->set[s].desc.doppler.v0;
		overflow[0] = 0.0;
		overflow[1] = 0.0;
		overflow[2] = 0.0;
		overflow[3] = 0.0;
		overflow[4] = 0.0;
	}
}
__global__ void cfs_set_lghtcrv_shortcuts_krnl32(struct dat_t *ddat,
		struct crvrend_t **rend, struct pos_t **pos, float *overflow,
		int *posn, int s, int size) {
	/* nfrm_alloc-threaded kernel */
	int f = blockIdx.x * blockDim.x + threadIdx.x + 1;
	if (f <= size) {
		posn[0] = 0;
		rend[f] = &ddat->set[s].desc.lghtcrv.rend[f];//&cfs_lghtcrv->rend[f];
		pos[f] = &rend[f]->pos;
		posn[f] = pos[f]->n;
		if (f==1) {
			overflow[0] = 0.0;
			overflow[1] = 0.0;
			overflow[2] = 0.0;
			overflow[3] = 0.0;
			overflow[4] = 0.0;
		}
	}
}
__global__ void cfs_set_lghtcrv_shortcuts_krnl64(struct dat_t *ddat,
		struct crvrend_t **rend, struct pos_t **pos, double *overflow,
		int *posn, int s, int size, int *outbndarr) {
	/* nfrm_alloc-threaded kernel */
	int f = blockIdx.x * blockDim.x + threadIdx.x + 1;
	if (f <= size) {
		posn[0] = 0;
		rend[f] = &ddat->set[s].desc.lghtcrv.rend[f];//&cfs_lghtcrv->rend[f];
		pos[f] = &rend[f]->pos;
		posn[f] = pos[f]->n;
		outbndarr[0] = 0;
		if (f==1) {
			overflow[0] = 0.0;
			overflow[1] = 0.0;
			overflow[2] = 0.0;
			overflow[3] = 0.0;
			overflow[4] = 0.0;
		}
	}
}
__global__ void cfs_set_pos_ae_krnl(struct dat_t *ddat, struct pos_t **pos,
		int s, int size, unsigned char type, int v, int lghtcrv) {
	/* nfrm_alloc-threaded kernel */
	int i, j, f = blockIdx.x * blockDim.x + threadIdx.x + lghtcrv;

	if (f < (size+lghtcrv)) {
		switch (type) {
		case DELAY:
			for (i=0; i<3; i++) {
				for (j=0; j<3; j++) {
					pos[f]->ae[i][j] = ddat->set[s].desc.deldop.frame[f].view[v].ae[i][j];
					pos[f]->oe[i][j] = ddat->set[s].desc.deldop.frame[f].view[v].oe[i][j];
				}
			}
			break;
		case DOPPLER:
			for (i=0; i<3; i++) {
				for (j=0; j<3; j++) {
					pos[f]->ae[i][j] = ddat->set[s].desc.doppler.frame[f].view[v].ae[i][j];
					pos[f]->oe[i][j] = ddat->set[s].desc.doppler.frame[f].view[v].oe[i][j];
				}
			}
			break;
		case LGHTCRV:
//			for (i=0; i<3; i++) {
//				for (j=0; j<3; j++) {
//					pos[f]->ae[i][j] = ddat->set[s].desc.lghtcrv.rend[f].ae[i][j];
//					pos[f]->oe[i][j] = ddat->set[s].desc.lghtcrv.rend[f].oe[i][j];
//					pos[f]->se[i][j] = ddat->set[s].desc.lghtcrv.rend[f].se[i][j];
//				}
//			}
			pos[f]->ae[0][0] = ddat->set[s].desc.lghtcrv.rend[f].ae[0][0];
			pos[f]->ae[0][1] = ddat->set[s].desc.lghtcrv.rend[f].ae[0][1];
			pos[f]->ae[0][2] = ddat->set[s].desc.lghtcrv.rend[f].ae[0][2];
			pos[f]->ae[1][0] = ddat->set[s].desc.lghtcrv.rend[f].ae[1][0];
			pos[f]->ae[1][1] = ddat->set[s].desc.lghtcrv.rend[f].ae[1][1];
			pos[f]->ae[1][2] = ddat->set[s].desc.lghtcrv.rend[f].ae[1][2];
			pos[f]->ae[2][0] = ddat->set[s].desc.lghtcrv.rend[f].ae[2][0];
			pos[f]->ae[2][1] = ddat->set[s].desc.lghtcrv.rend[f].ae[2][1];
			pos[f]->ae[2][2] = ddat->set[s].desc.lghtcrv.rend[f].ae[2][2];

			pos[f]->oe[0][0] = ddat->set[s].desc.lghtcrv.rend[f].oe[0][0];
			pos[f]->oe[0][1] = ddat->set[s].desc.lghtcrv.rend[f].oe[0][1];
			pos[f]->oe[0][2] = ddat->set[s].desc.lghtcrv.rend[f].oe[0][2];
			pos[f]->oe[1][0] = ddat->set[s].desc.lghtcrv.rend[f].oe[1][0];
			pos[f]->oe[1][1] = ddat->set[s].desc.lghtcrv.rend[f].oe[1][1];
			pos[f]->oe[1][2] = ddat->set[s].desc.lghtcrv.rend[f].oe[1][2];
			pos[f]->oe[2][0] = ddat->set[s].desc.lghtcrv.rend[f].oe[2][0];
			pos[f]->oe[2][1] = ddat->set[s].desc.lghtcrv.rend[f].oe[2][1];
			pos[f]->oe[2][2] = ddat->set[s].desc.lghtcrv.rend[f].oe[2][2];

			pos[f]->se[0][0] = ddat->set[s].desc.lghtcrv.rend[f].se[0][0];
			pos[f]->se[0][1] = ddat->set[s].desc.lghtcrv.rend[f].se[0][1];
			pos[f]->se[0][2] = ddat->set[s].desc.lghtcrv.rend[f].se[0][2];
			pos[f]->se[1][0] = ddat->set[s].desc.lghtcrv.rend[f].se[1][0];
			pos[f]->se[1][1] = ddat->set[s].desc.lghtcrv.rend[f].se[1][1];
			pos[f]->se[1][2] = ddat->set[s].desc.lghtcrv.rend[f].se[1][2];
			pos[f]->se[2][0] = ddat->set[s].desc.lghtcrv.rend[f].se[2][0];
			pos[f]->se[2][1] = ddat->set[s].desc.lghtcrv.rend[f].se[2][1];
			pos[f]->se[2][2] = ddat->set[s].desc.lghtcrv.rend[f].se[2][2];
			break;
		}
		/* Single-thread task */
		if (f == (0+lghtcrv)) {
			if ((type == LGHTCRV) || (type == POS))
				pos[f]->bistatic = 1;
			else if ((type == DELAY) || (type == DOPPLER))
				pos[f]->bistatic = 0;

		}
	}
}
__global__ void cfs_set_posbnd_krnl(struct par_t *dpar, struct dat_t *ddat,
		struct pos_t **pos, int s, int f, unsigned char type) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		dpar->posbnd = 1;
		switch (type) {
		case DELAY:
			dpar->posbnd_logfactor += ddat->set[s].desc.deldop.frame[f].dof *
				pos[f]->posbnd_logfactor;
			break;
		case DOPPLER:
			dpar->posbnd_logfactor += ddat->set[s].desc.doppler.frame[f].dof *
				pos[f]->posbnd_logfactor;
			break;
		case LGHTCRV:
			if (pos[f]->bistatic)
				cfs_lghtcrv_posbnd_logfactor += 0.5 * pos[f]->posbnd_logfactor;
			else
				cfs_lghtcrv_posbnd_logfactor += pos[f]->posbnd_logfactor;
			break;
		}
	}
}
__global__ void cfs_get_exclude_seen_krnl(struct par_t *dpar, struct pos_t **pos,
		int4 *xylim, int size, int lghtcrv) {
	/* nframes-threaded kernel */
	int f = blockIdx.x * blockDim.x + threadIdx.x + lghtcrv;

	if (f<(size+lghtcrv)) {
		if (lghtcrv) {
			xylim[0].w = xylim[0].x = xylim[0].y = xylim[0].z = 0;
		}
		cfs_exclude_seen = dpar->exclude_seen;
		xylim[f].w = pos[f]->xlim[0];
		xylim[f].x = pos[f]->xlim[1];
		xylim[f].y = pos[f]->ylim[0];
		xylim[f].z = pos[f]->ylim[1];
		if (f==0)
		cfs_exclude_seen = dpar->exclude_seen;
	}
}
__global__ void cf_mark_pixels_seen_krnl32(struct par_t *dpar,
		struct mod_t *dmod, struct pos_t **pos, int4 *xylim, int npixels,
		int xspan, int f) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int k = (offset % xspan) + xylim[f].w;
	int l = (offset / xspan) + xylim[f].y;
	int facetnum;

	if (offset < npixels) {
		if ((pos[f]->cose_s[offset] > dpar->mincosine_seen)
				&& (pos[f]->f[k][l] >= 0)) {
			facetnum = pos[f]->f[k][l];
			//c = cf_pos->comp[k][l];
			dmod->shape.comp[0].real.f[facetnum].seen = 1;
		}
	}
}
__global__ void cf_mark_pixels_seen_krnl64(struct par_t *dpar,
		struct mod_t *dmod, struct pos_t **pos, int4 *xylim, int npixels,
		int xspan, int f) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int k = (offset % xspan) + xylim[f].w;
	int l = (offset / xspan) + xylim[f].y;
	int facetnum;

	if (offset < npixels) {
		if ((pos[f]->cose[k][l] > dpar->mincosine_seen)
				&& (pos[f]->f[k][l] >= 0)) {
			facetnum = pos[f]->f[k][l];
			//c = cf_pos->comp[k][l];
			dmod->shape.comp[0].real.f[facetnum].seen = 1;
		}
	}
}
__global__ void cf_set_badradar_krnl(struct par_t *dpar,
		struct dat_t *ddat, int s, int f, unsigned char type) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		switch (type) {
		case DELAY:
			dpar->badradar = 1;
			dpar->badradar_logfactor += ddat->set[s].desc.deldop.frame[f].dof *
					ddat->set[s].desc.deldop.frame[f].badradar_logfactor /
					ddat->set[s].desc.deldop.nviews;
			break;
		case DOPPLER:
			dpar->badradar = 1;
			dpar->badradar_logfactor += ddat->set[s].desc.doppler.frame[f].dof *
					ddat->set[s].desc.doppler.frame[f].badradar_logfactor /
					ddat->set[s].desc.doppler.nviews;
			break;
		}
	}
}
__global__ void cf_add_fit_store_krnl1_32(struct dat_t *ddat, float **fit_store,
		int nThreads, int s, int f, unsigned char type) {
	/* ndel*ndop-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;

	if (offset < (nThreads)) {
		switch (type) {
		case DELAY:
			fit_store[f][offset] += ddat->set[s].desc.deldop.frame[f].fit_s[offset];
			break;
		case DOPPLER:
			fit_store[f][offset] += ddat->set[s].desc.doppler.frame[f].fit_s[offset];
			break;
		}
	}
}
__global__ void cf_add_fit_store_krnl1_64(struct dat_t *ddat, double **fit_store,
		int nThreads, int s, int f, unsigned char type, int *ndel) {
	/* ndel*ndop-threaded kernel */
	int idel, idop, offset = blockIdx.x * blockDim.x + threadIdx.x;
	idel = offset % ndel[f] + 1;
	idop = offset / ndel[f] + 1;

	if (offset < (nThreads)) {
		switch (type) {
		case DELAY:
			fit_store[f][offset] += ddat->set[s].desc.deldop.frame[f].fit[idel][idop];
			break;
		case DOPPLER:
			fit_store[f][offset] += ddat->set[s].desc.doppler.frame[f].fit[idop];
			break;
		}
	}
}
__global__ void cf_add_fit_store_krnl2_32(struct dat_t *ddat, int s, int f,
		float *overflow, unsigned char type) {
	/* ndel*ndop-threaded kernel */
	if (threadIdx.x == 0) {
		switch (type) {
		case DELAY:
			atomicAdd(&overflow[0], (float)ddat->set[s].desc.deldop.frame[f].overflow_o2);
			atomicAdd(&overflow[1], (float)ddat->set[s].desc.deldop.frame[f].overflow_m2);
			atomicAdd(&overflow[2], (float)ddat->set[s].desc.deldop.frame[f].overflow_xsec);
			atomicAdd(&overflow[3], (float)ddat->set[s].desc.deldop.frame[f].overflow_delmean);
			atomicAdd(&overflow[4], (float)ddat->set[s].desc.deldop.frame[f].overflow_dopmean);
			break;
		case DOPPLER:
			atomicAdd(&overflow[0], (float)ddat->set[s].desc.doppler.frame[f].overflow_o2);
			atomicAdd(&overflow[1], (float)ddat->set[s].desc.doppler.frame[f].overflow_m2);
			atomicAdd(&overflow[2], (float)ddat->set[s].desc.doppler.frame[f].overflow_xsec);
			atomicAdd(&overflow[3], (float)ddat->set[s].desc.doppler.frame[f].overflow_dopmean);
		}
	}
}
__global__ void cf_add_fit_store_krnl2_64(struct dat_t *ddat, int s, int f,
		double *overflow, unsigned char type) {
	/* ndel*ndop-threaded kernel */
	if (threadIdx.x == 0) {
		switch (type) {
		case DELAY:
			atomicAdd(&overflow[0], ddat->set[s].desc.deldop.frame[f].overflow_o2);
			atomicAdd(&overflow[1], ddat->set[s].desc.deldop.frame[f].overflow_m2);
			atomicAdd(&overflow[2], ddat->set[s].desc.deldop.frame[f].overflow_xsec);
			atomicAdd(&overflow[3], ddat->set[s].desc.deldop.frame[f].overflow_delmean);
			atomicAdd(&overflow[4], ddat->set[s].desc.deldop.frame[f].overflow_dopmean);
			break;
		case DOPPLER:
			atomicAdd(&overflow[0], ddat->set[s].desc.doppler.frame[f].overflow_o2);
			atomicAdd(&overflow[1], ddat->set[s].desc.doppler.frame[f].overflow_m2);
			atomicAdd(&overflow[2], ddat->set[s].desc.doppler.frame[f].overflow_xsec);
			atomicAdd(&overflow[3], ddat->set[s].desc.doppler.frame[f].overflow_dopmean);
		}
	}
}
__global__ void cf_finish_fit_store_krnl32(struct dat_t *ddat, float **fit_store,
		int s, int f, int nThreads, unsigned char type) {
	/* multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;

	if (offset < nThreads) {
		switch (type) {
		case DELAY:
			ddat->set[s].desc.deldop.frame[f].fit_s[offset] = fit_store[f][offset];
			break;
		case DOPPLER:
			ddat->set[s].desc.doppler.frame[f].fit_s[offset] = fit_store[f][offset];
			break;
		}
	}
}
__global__ void cf_finish_fit_store_krnl64(struct dat_t *ddat, double **fit_store,
		int s, int f, int nThreads, unsigned char type, int *ndel) {
	/* multi-threaded kernel */
	int idel, idop, offset = blockIdx.x * blockDim.x + threadIdx.x;
	idel = offset % ndel[f] + 1;
	idop = offset / ndel[f] + 1;
	if (offset < nThreads) {
		switch (type) {
		case DELAY:
			ddat->set[s].desc.deldop.frame[f].fit[idel][idop] = fit_store[f][offset];
			break;
		case DOPPLER:
			ddat->set[s].desc.doppler.frame[f].fit[idop] = fit_store[f][offset];
			break;
		}
	}
}
__global__ void cf_finish_fit_krnl2_32(struct dat_t *ddat,
		float *overflow, int s, int f, unsigned char type) {
	/* Single-threaded Kernel */
	int nviews;
	if (threadIdx.x == 0) {
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
__global__ void cf_finish_fit_krnl2_64(struct dat_t *ddat,
		double *overflow, int s, int f, unsigned char type) {
	/* Single-threaded Kernel */
	int nviews;
	if (threadIdx.x == 0) {
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
__global__ void cf_gamma_trans_krnl(struct par_t *dpar, struct dat_t *ddat,
		int s, int f, int nThreads, unsigned char type, int *ndel, int dbl) {
	/* Multi-threaded kernel */
	int idel, idop, offset = blockIdx.x * blockDim.x + threadIdx.x;
	idel = offset % ndel[f] + 1;
	idop = offset / ndel[f] + 1;
	if (offset < nThreads) {
		/*  Carry out a gamma transformation on the fit image if requested  */
		if (dpar->dd_gamma != 1.0) {
			switch (type) {
			case DELAY:
				if (dbl)
					dev_gamma_trans64(&ddat->set[s].desc.deldop.frame[f].fit[idel][idop],
							dpar->dd_gamma);
				else
					dev_gamma_trans32(&ddat->set[s].desc.deldop.frame[f].fit_s[offset],
							dpar->dd_gamma);
				break;
			case DOPPLER:
				//cf_dop_frame->fit[offset] = fit[offset];
				break;
			}
		}
	}
}

__host__ void calc_deldop_gpu32(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, struct vertices_t **verts, int s, int nframes,
		int nviews, unsigned char type, int nf, struct pos_t **pos, struct
		deldopfrm_t **frame, struct deldopview_t **view0, float **fit_store,
		int *ndel, int *ndop, int *posn, int4 *xylim, float
		*overflow, int *outbndarr, cudaStream_t *cf_stream)
{
	float3 orbit_off3;
	orbit_off3.x = orbit_off3.y = orbit_off3.z = 0.0;
	int v0_index, exclude_seen, f, v2, c=0, yspan;
	dim3 BLKdd[nframes], BLKpx[nframes], THD, THD9, BLKfrm,THD64;
	THD.x = maxThreadsPerBlock; THD9.x = 9;	THD64.x = 64;
	int4 hxylim[nframes];
	int hndop[nframes], hndel[nframes], hposn[nframes],
		houtbndarr[nframes], xspan[nframes], nThreadspx[nframes], nThreadsdd[nframes],
		nThreadspx1[nframes], v[nviews+1];
	BLKfrm = floor((THD64.x - 1 + nframes)/THD64.x);

	/* Set deldop, frame, view0, and pos in nframes streamed kernels */
	cfs_set_deldop_shortcuts_krnl32<<<BLKfrm,THD64>>>(ddat, frame, pos,
			view0, ndop, ndel, overflow, posn, s, nframes);
	checkErrorAfterKernelLaunch("cfs_set_deldop_shortcuts_krnl");
	gpuErrchk(cudaMemcpy(&hndel, ndel, sizeof(int)*nframes, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(&hndop, ndop, sizeof(int)*nframes, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(&hposn, posn, sizeof(int)*nframes, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&v0_index, cfs_v0_index, sizeof(int),
				0, cudaMemcpyDeviceToHost));

	for (f=0; f<nframes; f++) {
		nThreadsdd[f] = hndel[f]*hndop[f];
		BLKdd[f].x = floor((THD.x - 1 + nThreadsdd[f]) /	THD.x);
		nThreadspx[f] = (2 * hposn[f] + 1) * (2 * hposn[f] + 1);
		BLKpx[f].x = floor((THD.x -1 + nThreadspx[f]) / THD.x);

		/* If smearing is being modeled, initialize variables that
		 * will be used to sum results calculated for individual views.  */
		if (nviews > 1)
			/* Allocate fit_store as a single pointer, originally a double
			 * pointer. This also initializes the entire array to zero. */
			gpuErrchk(cudaMalloc((void**)&fit_store[f], sizeof(float) * nThreadsdd[f]));
	}
	/*  Loop over all views for this (smeared) frame, going in an order that
        ends with the view corresponding to the epoch listed for this frame
        in the obs file; this way we can use the calculated information for
        that view in the "write" action screen and disk output that follows   */
	for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
		v[v2] = v2 % nviews;
		/* Launch 9-threaded kernel to set pos->ae,pos->oe,pos->bistatic.*/
		cfs_set_pos_ae_krnl<<<BLKfrm,THD64>>>(ddat, pos, s, nframes,
				type, v[v2], 0);
	}
	for (f=0; f<nframes; f++)
		/* Launch posclr_krnl to initialize POS view */
		posclr_krnl<<<BLKpx[f],THD,0,cf_stream[f]>>>(pos, posn, f,FP64, 0);
	checkErrorAfterKernelLaunch("posclr_krnl (calc_fits_gpu2)");

	/* Call posvis to get facet number, scattering angle, distance
	 * toward Earth at center of each POS pixel; set flag posbnd if any model
	 * portion extends beyond POS frame limits.
	 * NOTE: Limited to single component for now */
	for (v2=v0_index+1; v2<=v0_index+nviews; v2++)
		posvis_gpu32(dpar, dmod, ddat, pos, verts, orbit_off3, hposn,
				outbndarr, s, nframes, 0, nf, 0, c, type, cf_stream, 0);

	gpuErrchk(cudaMemcpy(&houtbndarr, outbndarr, sizeof(int)*nframes,
			cudaMemcpyDeviceToHost));
//
//	f=0;
//	dbg_print_pos_arrays_full32(pos, f, nThreadspx[f], hposn[f]);



	for (f=0; f<nframes; f++)
		for (v2=v0_index+1; v2<=v0_index+nviews; v2++)
			if ((houtbndarr[f]) && (v[v2] == v0_index)) {
			/* Call single-threaded kernel to set dpar->posbnd and dpar->posbnd_logfactor */
			cfs_set_posbnd_krnl<<<1,1,0,cf_stream[f]>>>(dpar, ddat, pos, s, f, type);
			houtbndarr[f]=0;
			} checkErrorAfterKernelLaunch("cfs_set_posbnd_krnl");

	/* Get xlim and ylim and exclude_seen flag */
	cfs_get_exclude_seen_krnl<<<BLKfrm,THD64>>>(dpar,pos,xylim,nframes, 0);
	checkErrorAfterKernelLaunch("cfs_get_exclude_krnl");

	/* Now copy the flag and all frame pos's xlim and ylim values back from GPU */
	gpuErrchk(cudaMemcpyFromSymbol(&exclude_seen, cfs_exclude_seen, sizeof(int),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(&hxylim, xylim, sizeof(int4)*nframes, cudaMemcpyDeviceToHost));

	/* Calculate launch parameters for all frames */
	for (f=0; f<nframes; f++) {
		xspan[f] = hxylim[f].x - hxylim[f].w + 1;
		yspan = hxylim[f].z - hxylim[f].y + 1;
		nThreadspx1[f] = xspan[f] * yspan;
		BLKpx[f].x = (THD.x -1 + nThreadspx1[f]) / THD.x;
	}
	for (f=0; f<nframes; f++) {
		for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
			/* Go through all POS pixels which are visible with low enough
			 * scattering angle and mark the facets which project onto their
			 * centers as having been "seen" at least once                   */
			if (s != exclude_seen && v[v2] == v0_index)
				cf_mark_pixels_seen_krnl32<<<BLKpx[f],THD,0,cf_stream[f]>>>(
						dpar, dmod, pos, xylim, nThreadspx1[f], xspan[f], f);

			/* Zero out the fit delay-Doppler image, then call pos2deldop to
			 * create the fit image by mapping power from the plane of the sky
			 * to delay-Doppler space.                             */
			clrvect_krnl<<<BLKdd[f],THD,0,cf_stream[f]>>>(ddat, nThreadsdd[f], s, f, FP64);
		} checkErrorAfterKernelLaunch("clrvect_krnl and cf_mark_pixels_seen_streams_krnl");
	}

	for (v2=v0_index+1; v2<=v0_index+nviews; v2++)
		pos2deldop_gpu32(dpar, dmod, ddat, pos, frame, xylim, ndel, ndop,
				0.0,0.0,0.0,0, s, nframes, v[v2], outbndarr, cf_stream);

	/* Copy the badradar flag returns for all frames to a host copy */
	gpuErrchk(cudaMemcpy(&houtbndarr, outbndarr, sizeof(int)*nframes,
			cudaMemcpyDeviceToHost));

	for (f=0; f<nframes; f++) {
		for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
			if (houtbndarr[f]) {
				/* Call single-threaded kernel to set badradar flag and
				 * associated badradar_logfactor			 */
				cf_set_badradar_krnl<<<1,1,0,cf_stream[f]>>>(dpar, ddat, s, f, type);
				checkErrorAfterKernelLaunch("cf_deldop_set_badradar_krnl (calc_fits_cuda)");
			}
		}
	}
	for (f=0; f<nframes; f++) {
		/* If smearing is being modeled, include delay-Doppler calculations
		 * from this view in the summed results for this frame  */
		if (nviews > 1) {
			/* Launch ndel*ndop-threaded kernel to add fit[i][j] to fit_store[i][j]*/
			cf_add_fit_store_krnl1_32<<<BLKdd[f],THD,0,cf_stream[f]>>>(
					ddat,fit_store,nThreadsdd[f],s,f, type);
			cf_add_fit_store_krnl2_32<<<1,1>>>(ddat, s, f, overflow, type);
		}
	} checkErrorAfterKernelLaunch("cf_add_fit_store_streams_krnl1 and 2");

	/* If smearing is being modeled, compute mean values over all views for
	 * this frame and store them in the standard frame structure.  This kernel
	 * also carries out the gamma transformation on the fit image if the
	 * par->dd_gamma flag is not set  */
	if (nviews > 1) {
		for (f=0; f<nframes; f++) {

			cf_finish_fit_store_krnl32<<<BLKdd[f],THD,0,cf_stream[f]>>>(
					ddat, fit_store, s, f, nThreadsdd[f], type);

			cf_finish_fit_krnl2_32<<<1,1,0,cf_stream[f]>>>(ddat, overflow, s, f, type);
		} checkErrorAfterKernelLaunch("cf_finish_fit_store_streams kernels");
		cudaFree(fit_store);
	}

	for (f=0; f<nframes; f++) {
		cf_gamma_trans_krnl<<<BLKdd[f],THD,0,cf_stream[f]>>>(dpar, ddat, s, f,
				nThreadsdd[f], type, ndel, FP64);
	} checkErrorAfterKernelLaunch("cf_gamma_trans_krnl");
}

__host__ void calc_deldop_gpu64(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, struct vertices_t **verts, int s, int nframes,
		int nviews, unsigned char type, int nf, struct pos_t **pos, struct
		deldopfrm_t **frame, struct deldopview_t **view0, double **fit_store,
		int *ndel, int *ndop, int *posn, int4 *xylim, double
		*overflow, int *outbndarr, cudaStream_t *cf_stream)
{
	double3 orbit_off3;
	orbit_off3.x = orbit_off3.y = orbit_off3.z = 0.0;
	int v0_index, exclude_seen, f, v2, c=0, yspan;
	dim3 BLKdd[nframes], BLKpx[nframes], THD, THD9, BLKfrm,THD64;
	THD.x = maxThreadsPerBlock; THD9.x = 9;	THD64.x = 64;
	int4 hxylim[nframes];
	int hndop[nframes], hndel[nframes], hposn[nframes],
		houtbndarr[nframes], xspan[nframes], nThreadspx[nframes], nThreadsdd[nframes],
		nThreadspx1[nframes], v[nviews+1];
	BLKfrm = floor((THD64.x - 1 + nframes)/THD64.x);

	/* Set deldop, frame, view0, and pos in nframes streamed kernels */
	cfs_set_deldop_shortcuts_krnl64<<<BLKfrm,THD64>>>(ddat, frame, pos,
			view0, ndop, ndel, overflow, posn, s, nframes);
	checkErrorAfterKernelLaunch("cfs_set_deldop_shortcuts_krnl");
	gpuErrchk(cudaMemcpy(&hndel, ndel, sizeof(int)*nframes, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(&hndop, ndop, sizeof(int)*nframes, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(&hposn, posn, sizeof(int)*nframes, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&v0_index, cfs_v0_index, sizeof(int),
				0, cudaMemcpyDeviceToHost));

	for (f=0; f<nframes; f++) {
		nThreadsdd[f] = hndel[f]*hndop[f];
		BLKdd[f].x = floor((THD.x - 1 + nThreadsdd[f]) /	THD.x);
		nThreadspx[f] = (2 * hposn[f] + 1) * (2 * hposn[f] + 1);
		BLKpx[f].x = floor((THD.x -1 + nThreadspx[f]) / THD.x);

		/* If smearing is being modeled, initialize variables that
		 * will be used to sum results calculated for individual views.  */
		if (nviews > 1)
			/* Allocate fit_store as a single pointer, originally a double
			 * pointer. This also initializes the entire array to zero. */
			gpuErrchk(cudaMalloc((void**)&fit_store[f], sizeof(double) * nThreadsdd[f]));
	}
	/*  Loop over all views for this (smeared) frame, going in an order that
        ends with the view corresponding to the epoch listed for this frame
        in the obs file; this way we can use the calculated information for
        that view in the "write" action screen and disk output that follows   */
	for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
		v[v2] = v2 % nviews;
		/* Launch 9-threaded kernel to set pos->ae,pos->oe,pos->bistatic.*/
		cfs_set_pos_ae_krnl<<<BLKfrm,THD64>>>(ddat, pos, s, nframes,
				type, v[v2], 0);
	}
	for (f=0; f<nframes; f++)
		/* Launch posclr_krnl to initialize POS view */
		posclr_krnl<<<BLKpx[f],THD,0,cf_stream[f]>>>(pos, posn, f, FP64, 0);
	checkErrorAfterKernelLaunch("posclr_krnl");

	/* Call posvis to get facet number, scattering angle, distance
	 * toward Earth at center of each POS pixel; set flag posbnd if any model
	 * portion extends beyond POS frame limits.
	 * NOTE: Limited to single component for now */
	for (v2=v0_index+1; v2<=v0_index+nviews; v2++)
		posvis_gpu64(dpar, dmod, ddat, pos, verts, orbit_off3, hposn,
				outbndarr, s, nframes, 0, nf, 0, c, type, cf_stream, 0);

	gpuErrchk(cudaMemcpy(&houtbndarr, outbndarr, sizeof(int)*nframes,
			cudaMemcpyDeviceToHost));
//
	f=0;
	dbg_print_pos_arrays_full64(pos, f, nThreadspx[f], hposn[f]);


	for (f=0; f<nframes; f++)
		for (v2=v0_index+1; v2<=v0_index+nviews; v2++)
			if ((houtbndarr[f]) && (v[v2] == v0_index)) {
			/* Call single-threaded kernel to set dpar->posbnd and dpar->posbnd_logfactor */
			cfs_set_posbnd_krnl<<<1,1,0,cf_stream[f]>>>(dpar, ddat, pos, s, f, type);
			houtbndarr[f]=0;
			} checkErrorAfterKernelLaunch("cfs_set_posbnd_krnl");

	/* Get xlim and ylim and exclude_seen flag */
	cfs_get_exclude_seen_krnl<<<BLKfrm,THD64>>>(dpar,pos,xylim,nframes, 0);
	checkErrorAfterKernelLaunch("cfs_get_exclude_krnl");

	/* Now copy the flag and all frame pos's xlim and ylim values back from GPU */
	gpuErrchk(cudaMemcpyFromSymbol(&exclude_seen, cfs_exclude_seen, sizeof(int),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(&hxylim, xylim, sizeof(int4)*nframes, cudaMemcpyDeviceToHost));

	/* Calculate launch parameters for all frames */
	for (f=0; f<nframes; f++) {
		xspan[f] = hxylim[f].x - hxylim[f].w + 1;
		yspan = hxylim[f].z - hxylim[f].y + 1;
		nThreadspx1[f] = xspan[f] * yspan;
		BLKpx[f].x = (THD.x -1 + nThreadspx1[f]) / THD.x;
	}
	for (f=0; f<nframes; f++) {
		for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
			/* Go through all POS pixels which are visible with low enough
			 * scattering angle and mark the facets which project onto their
			 * centers as having been "seen" at least once                   */
			if (s != exclude_seen && v[v2] == v0_index)
				cf_mark_pixels_seen_krnl64<<<BLKpx[f],THD,0,cf_stream[f]>>>(
						dpar, dmod, pos, xylim, nThreadspx1[f], xspan[f], f);

			/* Zero out the fit delay-Doppler image, then call pos2deldop to
			 * create the fit image by mapping power from the plane of the sky
			 * to delay-Doppler space.                             */
			clrvect_krnl<<<BLKdd[f],THD,0,cf_stream[f]>>>(ddat, nThreadsdd[f], s, f, FP64);
		} checkErrorAfterKernelLaunch("clrvect_krnl and cf_mark_pixels_seen_streams_krnl");
	}

	for (v2=v0_index+1; v2<=v0_index+nviews; v2++)
		pos2deldop_gpu64(dpar, dmod, ddat, pos, frame, xylim, ndel, ndop,
				0.0,0.0,0.0,0, s, nframes, v[v2], outbndarr, cf_stream);

	/* Copy the badradar flag returns for all frames to a host copy */
	gpuErrchk(cudaMemcpy(&houtbndarr, outbndarr, sizeof(int)*nframes,
			cudaMemcpyDeviceToHost));

//dbg_print_deldop_fit(ddat, 0, 0, "1080Ti_deldop_fit.csv");


	for (f=0; f<nframes; f++) {
		for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
			if (houtbndarr[f]) {
				/* Call single-threaded kernel to set badradar flag and
				 * associated badradar_logfactor			 */
				cf_set_badradar_krnl<<<1,1,0,cf_stream[f]>>>(dpar, ddat, s, f, type);
				checkErrorAfterKernelLaunch("cf_deldop_set_badradar_krnl (calc_fits_cuda)");
			}
		}
	}
	for (f=0; f<nframes; f++) {
		/* If smearing is being modeled, include delay-Doppler calculations
		 * from this view in the summed results for this frame  */
		if (nviews > 1) {
			/* Launch ndel*ndop-threaded kernel to add fit[i][j] to fit_store[i][j]*/
			cf_add_fit_store_krnl1_64<<<BLKdd[f],THD,0,cf_stream[f]>>>(
					ddat,fit_store,nThreadsdd[f],s,f, type, ndel);
			cf_add_fit_store_krnl2_64<<<1,1>>>(ddat, s, f, overflow, type);
		}
	} checkErrorAfterKernelLaunch("cf_add_fit_store_streams_krnl1 and 2");

	/* If smearing is being modeled, compute mean values over all views for
	 * this frame and store them in the standard frame structure.  This kernel
	 * also carries out the gamma transformation on the fit image if the
	 * par->dd_gamma flag is not set  */
	if (nviews > 1) {
		for (f=0; f<nframes; f++) {

			cf_finish_fit_store_krnl64<<<BLKdd[f],THD,0,cf_stream[f]>>>(
					ddat, fit_store, s, f, nThreadsdd[f], type, ndel);

			cf_finish_fit_krnl2_64<<<1,1,0,cf_stream[f]>>>(ddat, overflow, s, f, type);
		} checkErrorAfterKernelLaunch("cf_finish_fit_store_streams kernels");
		cudaFree(fit_store);
	}

	for (f=0; f<nframes; f++) {
		cf_gamma_trans_krnl<<<BLKdd[f],THD,0,cf_stream[f]>>>(dpar, ddat, s, f,
				nThreadsdd[f], type, ndel, FP64);
	} checkErrorAfterKernelLaunch("cf_gamma_trans_krnl");


}

__host__ void calc_doppler_gpu32(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, struct vertices_t **verts, int s, int nframes,
		int nviews, unsigned char type, int nf, struct pos_t **pos, struct
		dopfrm_t **frame, struct dopview_t **view0, float **fit_store, int
		*ndop, int *posn, int4 *xylim, float *overflow, int
		*outbndarr, cudaStream_t *cf_stream)
{
	float3 orbit_off3;
	orbit_off3.x = orbit_off3.y = orbit_off3.z = 0.0;
	int v0_index, exclude_seen, f, v2, c=0, yspan;
	dim3 BLKd[nframes], BLKpx[nframes],THD,THD9,THD64,BLKfrm;
	THD.x = maxThreadsPerBlock; THD9.x = 9; THD64.x = 64;
	int4 hxylim[nframes];
	int hndop[nframes], hposn[nframes], houtbndarr[nframes],
		xspan[nframes], nThreadspx[nframes], nThreadspx1[nframes], v[nviews+1];
	BLKfrm = floor((THD64.x - 1 + nframes)/THD64.x);

	/* Set doppler, frame, view0, and pos in nframes-threaded kernels */
	cfs_set_doppler_shortcuts_krnl32<<<BLKfrm,THD64>>>(ddat, frame,
			pos, view0, overflow, ndop, posn, s, nframes);
	checkErrorAfterKernelLaunch("cfs_set_doppler_shortcuts_krnl2");
	gpuErrchk(cudaMemcpy(&hndop, ndop, sizeof(int)*nframes, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(&hposn, posn, sizeof(int)*nframes, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&v0_index, cfs_v0_index, sizeof(int), 0,
			cudaMemcpyDeviceToHost));

	for (f=0; f<nframes; f++) {
		/* Calculate launch parameters needed later */
		BLKd[f] = floor((THD.x - 1 + hndop[f]) / THD.x);
		nThreadspx[f] = (2 * hposn[f] + 1) * (2 * hposn[f] + 1);
		BLKpx[f] = floor((THD.x - 1 + nThreadspx[f]) / THD.x);

		/* If smearing is being modeled, initialize variables that
		 * will be used to sum results calculated for individual views.  */
		if (nviews > 1)
			/* Allocate fit_store as a single pointer, originally a double
			 * pointer. This also initializes the entire array to zero. */
			gpuErrchk(cudaMalloc((void**)&fit_store[f], sizeof(float)*ndop[f]));
	}

	/* Loop over all views for this (smeared) frame, going in an order that
	 * ends with the view corresponding to the epoch listed for this frame
	 * in the obs file; this way we can use the calculated information for
	 * that view in the "write" action screen and disk output that follows*/
	for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
		v[v2] = v2 % nviews;
		/* Launch 9-threaded kernel to set pos->ae,pos->oe,pos->bistatic.*/
		cfs_set_pos_ae_krnl<<<BLKfrm,THD64>>>(ddat, pos, s, nframes,
				type, v[v2], 0);
	} checkErrorAfterKernelLaunch("cfs_set_pos_ae_krnl2 in calc_doppler_gpu2");

	for (f=0; f<nframes; f++)
		/* Launch posclr_krnl to initialize POS view */
		posclr_krnl<<<BLKpx[f],THD,0,cf_stream[f]>>>(pos, posn, f, FP64, 0);
	checkErrorAfterKernelLaunch("posclr_krnl (calc_doppler_gpu2)");

	/* Call posvis_cuda_2 to get facet number, scattering angle, distance
	 * toward Earth at center of each POS pixel; set flag posbnd if any model
	 * portion extends beyond POS frame limits.*/
	/* NOTE: Limited to single component for now */
	for (v2=v0_index+1; v2<=v0_index+nviews; v2++)
		posvis_gpu32(dpar, dmod, ddat, pos, verts, orbit_off3, hposn,
				outbndarr, s, nframes, 0, nf, 0, c, type, cf_stream, 0);

	/* Copy the posbnd flag returns for all frames to a host copy */
	gpuErrchk(cudaMemcpy(&houtbndarr, outbndarr, sizeof(int)*nframes,
			cudaMemcpyDeviceToHost));

	for (f=0; f<nframes; f++)
		for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
			if ((houtbndarr[f]) && (v[v2] == v0_index)) {
				/* Call single-threaded kernel to set dpar->posbnd and dpar->posbnd_logfactor */
				cfs_set_posbnd_krnl<<<1,1,0,cf_stream[f]>>>(dpar, ddat, pos, s, f, type);
				houtbndarr[f]=0;
			} checkErrorAfterKernelLaunch("cfs_set_posbnd_krnl");
		}

	/* Get xlim and ylim and exclude_seen flag */
	cfs_get_exclude_seen_krnl<<<BLKfrm,THD64>>>(dpar,pos,xylim,nframes, 0);
	checkErrorAfterKernelLaunch("cfs_get_exclude_krnl");

	/* Now copy the flag and all frame pos's xlim and ylim values back from GPU */
	gpuErrchk(cudaMemcpyFromSymbol(&exclude_seen, cfs_exclude_seen, sizeof(int),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(&hxylim, xylim, sizeof(int4)*nframes, cudaMemcpyDeviceToHost));

	/* Calculate launch parameters for all frames */
	for (f=0; f<nframes; f++) {
		xspan[f] = hxylim[f].x - hxylim[f].w + 1;
		yspan = hxylim[f].z - hxylim[f].y + 1;
		nThreadspx1[f] = xspan[f] * yspan;
		BLKpx[f].x = (THD.x -1 + nThreadspx1[f]) / THD.x;
	}

	for (f=0; f<nframes; f++) {
		for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
			/* Go through all POS pixels which are visible with low enough
			 * scattering angle and mark the facets which project onto
			 * their centers as having been "seen" at least once           */
			if (s != exclude_seen && v[v2] == v0_index)
				cf_mark_pixels_seen_krnl32<<<BLKpx[f],THD,0,cf_stream[f]>>>(
						dpar, dmod, pos, xylim, nThreadspx1[f], xspan[f], f);

			/* Zero out the fit delay-Doppler image, then call pos2deldop to
			 * create the fit image by mapping power from the plane of the sky
			 * to delay-Doppler space.                             */
			clrvect_krnl<<<BLKd[f],THD,0,cf_stream[f]>>>(ddat, hndop[f], s, f, FP64);
		} checkErrorAfterKernelLaunch("clrvect_krnl and cf_mark_pixels_seen_streams_krnl");
	}


	/* Call pos2deldop to calculate the Doppler radar fit image */
	for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
		pos2doppler_gpu32(dpar, dmod, ddat, pos, frame, xylim,0.0,0.0,0.0, ndop,
				0, s, nframes, v[v2], outbndarr, cf_stream);
	}
	/* Copy the badradar flag returns for all frames to a host copy */
	gpuErrchk(cudaMemcpy(&houtbndarr, outbndarr, sizeof(int)*nframes,
			cudaMemcpyDeviceToHost));

	for (f=0; f<nframes; f++) {
		for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
			if (houtbndarr[f]) {
				/* Call single-threaded kernel to set badradar flag and
				 * associated badradar_logfactor			 */
				cf_set_badradar_krnl<<<1,1,0,cf_stream[f]>>>(dpar, ddat, s, f, type);
				checkErrorAfterKernelLaunch("cf_set_badradar_streams_krnl (calc_doppler_cuda_streams)");
			}
		}
	}

	for (f=0; f<nframes; f++) {
		/* If smearing is being modeled, include delay-Doppler calculations
		 * from this view in the summed results for this frame  */
		if (nviews > 1) {
			/* Launch ndop-threaded kernel to add fit[i][j] to
			 * fit_store[i][j]*/
			cf_add_fit_store_krnl1_32<<<BLKd[f],THD,0,cf_stream[f]>>>(
					ddat, fit_store, ndop[f], s, f, type);

			cf_add_fit_store_krnl2_32<<<1,1>>>(ddat, s, f, overflow, type);
		}
	} checkErrorAfterKernelLaunch("cf_add_fit_store_streams_krnl1 and 2 (calc_doppler_cuda_streams");

	/* If smearing is being modeled, compute mean values over all views for
	 * this frame and store them in the standard frame structure     */
	/* This kernel also carries out the gamma transformation on the fit
	 * image if the par->dd_gamma flag is not set  */
	if (nviews > 1) {
		for (f=0; f<nframes; f++) {

			cf_finish_fit_store_krnl32<<<BLKd[f],THD,0,cf_stream[f]>>>(
					ddat, fit_store, s, f, ndop[f], type);

			cf_finish_fit_krnl2_32<<<1,1,0,cf_stream[f]>>>(ddat, overflow, s, f, type);

//			cf_gamma_trans_krnl<<<BLKd[f],THD,0,cf_stream[f]>>>(dpar,
//					ddat, s, f, ndop[f], type);
		} checkErrorAfterKernelLaunch("cf_finish_fit_store_streams kernels and "
				"cf_gamma_trans_krnl (calc_doppler_cuda_streams");
		cudaFree(fit_store);
	}
}

__host__ void calc_doppler_gpu64(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, struct vertices_t **verts, int s, int nframes,
		int nviews, unsigned char type, int nf, struct pos_t **pos, struct
		dopfrm_t **frame, struct dopview_t **view0, double **fit_store, int
		*ndop, int *posn, int4 *xylim, double *overflow, int
		*outbndarr, cudaStream_t *cf_stream)
{
	double3 orbit_off3;
	orbit_off3.x = orbit_off3.y = orbit_off3.z = 0.0;
	int v0_index, exclude_seen, f, v2, c=0, yspan;
	dim3 BLKd[nframes], BLKpx[nframes],THD,THD9,THD64,BLKfrm;
	THD.x = maxThreadsPerBlock; THD9.x = 9; THD64.x = 64;
	int4 hxylim[nframes];
	int hndop[nframes], hposn[nframes], houtbndarr[nframes],
		xspan[nframes], nThreadspx[nframes], nThreadspx1[nframes], v[nviews+1];
	BLKfrm = floor((THD64.x - 1 + nframes)/THD64.x);

	/* Set doppler, frame, view0, and pos in nframes-threaded kernels */
	cfs_set_doppler_shortcuts_krnl64<<<BLKfrm,THD64>>>(ddat, frame,
			pos, view0, overflow, ndop, posn, s, nframes);
	checkErrorAfterKernelLaunch("cfs_set_doppler_shortcuts_krnl64");
	gpuErrchk(cudaMemcpy(&hndop, ndop, sizeof(int)*nframes, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(&hposn, posn, sizeof(int)*nframes, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&v0_index, cfs_v0_index, sizeof(int), 0,
			cudaMemcpyDeviceToHost));

	for (f=0; f<nframes; f++) {
		/* Calculate launch parameters needed later */
		BLKd[f] = floor((THD.x - 1 + hndop[f]) / THD.x);
		nThreadspx[f] = (2 * hposn[f] + 1) * (2 * hposn[f] + 1);
		BLKpx[f] = floor((THD.x - 1 + nThreadspx[f]) / THD.x);

		/* If smearing is being modeled, initialize variables that
		 * will be used to sum results calculated for individual views.  */
		if (nviews > 1)
			/* Allocate fit_store as a single pointer, originally a double
			 * pointer. This also initializes the entire array to zero. */
			gpuErrchk(cudaMalloc((void**)&fit_store[f], sizeof(double)*ndop[f]));
	}

	/* Loop over all views for this (smeared) frame, going in an order that
	 * ends with the view corresponding to the epoch listed for this frame
	 * in the obs file; this way we can use the calculated information for
	 * that view in the "write" action screen and disk output that follows*/
	for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
		v[v2] = v2 % nviews;
		/* Launch 9-threaded kernel to set pos->ae,pos->oe,pos->bistatic.*/
		cfs_set_pos_ae_krnl<<<BLKfrm,THD64>>>(ddat, pos, s, nframes,
				type, v[v2], 0);
	} checkErrorAfterKernelLaunch("cfs_set_pos_ae_krnl");

	for (f=0; f<nframes; f++)
		/* Launch posclr_krnl to initialize POS view */
		posclr_krnl<<<BLKpx[f],THD,0,cf_stream[f]>>>(pos, posn, f, FP64, 0);
	checkErrorAfterKernelLaunch("posclr_krnl ");

	/* Call posvis_cuda_2 to get facet number, scattering angle, distance
	 * toward Earth at center of each POS pixel; set flag posbnd if any model
	 * portion extends beyond POS frame limits.*/
	/* NOTE: Limited to single component for now */
	for (v2=v0_index+1; v2<=v0_index+nviews; v2++)
		posvis_gpu64(dpar, dmod, ddat, pos, verts, orbit_off3, hposn,
				outbndarr, s, nframes, 0, nf, 0, c, type, cf_stream, 0);

	/* Copy the posbnd flag returns for all frames to a host copy */
	gpuErrchk(cudaMemcpy(&houtbndarr, outbndarr, sizeof(int)*nframes,
			cudaMemcpyDeviceToHost));

	for (f=0; f<nframes; f++)
		for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
			if ((houtbndarr[f]) && (v[v2] == v0_index)) {
				/* Call single-threaded kernel to set dpar->posbnd and dpar->posbnd_logfactor */
				cfs_set_posbnd_krnl<<<1,1,0,cf_stream[f]>>>(dpar, ddat, pos, s, f, type);
				houtbndarr[f]=0;
			} checkErrorAfterKernelLaunch("cfs_set_posbnd_krnl");
		}

	/* Get xlim and ylim and exclude_seen flag */
	cfs_get_exclude_seen_krnl<<<BLKfrm,THD64>>>(dpar,pos,xylim,nframes, 0);
	checkErrorAfterKernelLaunch("cfs_get_exclude_krnl");

	/* Now copy the flag and all frame pos's xlim and ylim values back from GPU */
	gpuErrchk(cudaMemcpyFromSymbol(&exclude_seen, cfs_exclude_seen, sizeof(int),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(&hxylim, xylim, sizeof(int4)*nframes, cudaMemcpyDeviceToHost));

	/* Calculate launch parameters for all frames */
	for (f=0; f<nframes; f++) {
		xspan[f] = hxylim[f].x - hxylim[f].w + 1;
		yspan = hxylim[f].z - hxylim[f].y + 1;
		nThreadspx1[f] = xspan[f] * yspan;
		BLKpx[f].x = (THD.x -1 + nThreadspx1[f]) / THD.x;
	}

	for (f=0; f<nframes; f++) {
		for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
			/* Go through all POS pixels which are visible with low enough
			 * scattering angle and mark the facets which project onto
			 * their centers as having been "seen" at least once           */
			if (s != exclude_seen && v[v2] == v0_index)
				cf_mark_pixels_seen_krnl64<<<BLKpx[f],THD,0,cf_stream[f]>>>(
						dpar, dmod, pos, xylim, nThreadspx1[f], xspan[f], f);

			/* Zero out the fit delay-Doppler image, then call pos2deldop to
			 * create the fit image by mapping power from the plane of the sky
			 * to delay-Doppler space.                             */
			clrvect_krnl<<<BLKd[f],THD,0,cf_stream[f]>>>(ddat, hndop[f], s, f, FP64);
		} checkErrorAfterKernelLaunch("clrvect_krnl and cf_mark_pixels_seen_krnl");
	}


	/* Call pos2deldop to calculate the Doppler radar fit image */
	for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
		pos2doppler_gpu64(dpar, dmod, ddat, pos, frame, xylim,0.0,0.0,0.0, ndop,
				0, s, nframes, v[v2], outbndarr, cf_stream);
	}
	/* Copy the badradar flag returns for all frames to a host copy */
	gpuErrchk(cudaMemcpy(&houtbndarr, outbndarr, sizeof(int)*nframes,
			cudaMemcpyDeviceToHost));

	for (f=0; f<nframes; f++) {
		for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
			if (houtbndarr[f]) {
				/* Call single-threaded kernel to set badradar flag and
				 * associated badradar_logfactor			 */
				cf_set_badradar_krnl<<<1,1,0,cf_stream[f]>>>(dpar, ddat, s, f, type);
				checkErrorAfterKernelLaunch("cf_set_badradar_krnl");
			}
		}
	}

	for (f=0; f<nframes; f++) {
		/* If smearing is being modeled, include delay-Doppler calculations
		 * from this view in the summed results for this frame  */
		if (nviews > 1) {
			/* Launch ndop-threaded kernel to add fit[i][j] to
			 * fit_store[i][j]*/
			cf_add_fit_store_krnl1_64<<<BLKd[f],THD,0,cf_stream[f]>>>(
					ddat, fit_store, ndop[f], s, f, type, ndop);

			cf_add_fit_store_krnl2_64<<<1,1>>>(ddat, s, f, overflow, type);
		}
	} checkErrorAfterKernelLaunch("cf_add_fit_store_krnl1_64 and 2_64");

	/* If smearing is being modeled, compute mean values over all views for
	 * this frame and store them in the standard frame structure     */
	/* This kernel also carries out the gamma transformation on the fit
	 * image if the par->dd_gamma flag is not set  */
	if (nviews > 1) {
		for (f=0; f<nframes; f++) {

			cf_finish_fit_store_krnl64<<<BLKd[f],THD,0,cf_stream[f]>>>(
					ddat, fit_store, s, f, ndop[f], type, ndop);

			cf_finish_fit_krnl2_64<<<1,1,0,cf_stream[f]>>>(ddat, overflow, s, f, type);

//			cf_gamma_trans_krnl<<<BLKd[f],THD,0,cf_stream[f]>>>(dpar,
//					ddat, s, f, ndop[f], type);
		} checkErrorAfterKernelLaunch("cf_finish_fit_store_kernels and "
				"cf_gamma_trans_krnl");
		cudaFree(fit_store);
	}
}

////void calc_poset( struct par_t *par, struct mod_t *mod, struct poset_t *poset,
////		int s)
////{
////	const char *monthName[12] = {"Jan", "Feb", "Mar", "Apr", "May", "Jun",
////			"Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};
////	double orbit_offset[3] = {0.0, 0.0, 0.0};
////
////	FILE *fpopt;
////	char tempstring[MAXLEN], name[MAXLEN];
////	int year, mon, day, hour, min, sec, f, c, i, j, k, l, nrow_fit, ncol_fit, n_pos,
////	facetnum, x, y, v, v2;
////	double w[3], spin_colat, spin_azim, xoff, yoff, resamp_fact, resamp_x0, resamp_y0,
////	xcom_fit, ycom_fit, resamp_xwidth, resamp_ywidth, resamp_angle, oa[3][3],
////	to_earth[3], to_earth_lat, to_earth_long, rotphase, sa[3][3], to_sun[3],
////	to_sun_lat, to_sun_long, pab[3], pab_lat, pab_long, intensityfactor,
////	phi, theta, psi, intspin_body[3], badposet_logfactor_view;
////	double **fit_store;
////	struct posetfrm_t *frame;
////	struct posetview_t *view0;
////	struct pos_t *pos;
////
////	for (f=0; f<poset->nframes; f++) {
////
////		frame = &poset->frame[f];
////		view0 = &frame->view[poset->v0];
////		pos = &frame->pos;
////
////		ncol_fit = frame->ncol;
////		nrow_fit = frame->nrow;
////
////		/*  If smearing is being modeled, initialize variables that
////        will be used to sum results calculated for individual views  */
////
////		if (poset->nviews > 1) {
////			fit_store = matrix( 1, ncol_fit, 1, nrow_fit);
////			for (i=1; i<=ncol_fit; i++)
////				for (j=1; j<=nrow_fit; j++)
////					fit_store[i][j] = 0.0;
////		}
////
////		/*  Loop over all views for this (smeared) frame, going in an order that
////        ends with the view corresponding to the epoch listed for this frame
////        in the obs file; this way we can use the calculated information for
////        that view in the "write" action screen and disk output that follows   */
////
////		for (v2=poset->v0+1; v2<=poset->v0+poset->nviews; v2++) {
////			v = v2 % poset->nviews;
////
////			for (i=0; i<=2; i++)
////				for (j=0; j<=2; j++) {
////					pos->ae[i][j] = frame->view[v].ae[i][j];
////					pos->oe[i][j] = frame->view[v].oe[i][j];
////					pos->se[i][j] = frame->view[v].se[i][j];
////				}
////			pos->bistatic = 1;
////
////			/*  Initialize the plane-of-sky view  */
////
////			posclr( pos);
////
////			/*  Call routine posvis to get the facet number, scattering angle,
////          incidence angle, and distance toward Earth at the center of
////          each POS pixel; set the posbnd parameter to 1 if any portion
////          of the model extends beyond the POS frame limits.              */
////
////			for (c=0; c<mod->shape.ncomp; c++)
////				if (posvis( &mod->shape.comp[c].real, orbit_offset, pos,
////						(int) par->pos_smooth, 0, 0, c) && v == poset->v0) {
////					par->posbnd = 1;
////					if (pos->bistatic)
////						par->posbnd_logfactor += 0.5 * frame->dof * pos->posbnd_logfactor;
////					else
////						par->posbnd_logfactor += frame->dof * pos->posbnd_logfactor;
////				}
////
////			/*  Now view the model from the source (sun) and get the facet number
////          and distance toward the source of each pixel in this projected view;
////          use this information to determine which POS pixels are shadowed       */
////
////			if (pos->bistatic) {
////				for (c=0; c<mod->shape.ncomp; c++)
////					if (posvis( &mod->shape.comp[c].real, orbit_offset, pos,
////							0, 1, 0, c)) {
////						par->posbnd = 1;
////						par->posbnd_logfactor += 0.5 * frame->dof * pos->posbnd_logfactor;
////					}
////
////				/*  Identify and mask out shadowed POS pixels  */
////
////				posmask( pos, par->mask_tol);
////			}
////
////			/*  Go through all POS pixels which are visible and unshadowed with
////          sufficiently low scattering and incidence angles, and mark the facets
////          which project onto their centers as having been "seen" at least once   */
////
////			if (s != par->exclude_seen && v == poset->v0) {
////				for (k=pos->xlim[0]; k<=pos->xlim[1]; k++)
////					for (l=pos->ylim[0]; l<=pos->ylim[1]; l++) {
////						if ((pos->cose[k][l] > par->mincosine_seen)
////								&& (pos->cosi[k][l] > par->mincosine_seen)
////								&& (pos->f[k][l] >= 0)) {
////							facetnum = pos->f[k][l];
////							c = pos->comp[k][l];
////							mod->shape.comp[c].real.f[facetnum].seen = 1;
////						}
////					}
////			}
////
////			/*  Compute the sky rendering  */
////
////			intensityfactor = pow( pos->km_per_pixel/AU, 2.0);
////			apply_photo( mod, poset->ioptlaw, frame->view[v].solar_phase,
////					intensityfactor, pos, 0);
////
////			/*  Resample the sky rendering to get the model plane-of-sky image    */
////			/*  (if using bicubic interpolation or cubic convolution, force       */
////			/*  all model pixel values to be nonnegative)                         */
////			/*                                                                    */
////			/*  Implement the x and y COM offsets, xoff and yoff, by first        */
////			/*  using them to compute xcom_fit and ycom_fit -- the COM position   */
////			/*  in the fit image, relative to the center of the fit image -- and  */
////			/*  then shifting the resampled region in the *opposite* direction    */
////			/*  by the appropriate proportional amount.  Then implement the       */
////			/*  "northangle" setting (clockwise heading of north) by rotating     */
////			/*  the resampling grid *counterclockwise* by northangle.             */
////
////			n_pos = pos->n;
////			xoff = frame->off[0].val;
////			yoff = frame->off[1].val;
////			xcom_fit = (frame->colcom_vig - (ncol_fit + 1)/2.0) + xoff;
////			ycom_fit = (frame->rowcom_vig - (nrow_fit + 1)/2.0) + yoff;
////			resamp_fact = frame->fit.km_per_pixel / pos->km_per_pixel;
////			resamp_x0 = -xcom_fit*resamp_fact;
////			resamp_y0 = -ycom_fit*resamp_fact;
////			resamp_xwidth = resamp_fact*(ncol_fit - 1);
////			resamp_ywidth = resamp_fact*(nrow_fit - 1);
////			resamp_angle = -frame->northangle;
////			resampim( frame->pos.b, -n_pos, n_pos, -n_pos, n_pos,
////					frame->fit.b, 1, ncol_fit, 1, nrow_fit,
////					resamp_x0, resamp_xwidth, resamp_y0, resamp_ywidth, resamp_angle,
////					(int) par->poset_resample, (int) par->image_rebin);
////			if (par->poset_resample == BICUBIC || par->poset_resample == CUBICCONV) {
////				for (k=1; k<=ncol_fit; k++)
////					for (l=1; l<=nrow_fit; l++)
////						frame->fit.b[k][l] = MAX( 0.0, frame->fit.b[k][l]);
////			}
////
////			/*  Set the badposet flag and increase badposet_logfactor if the model   */
////			/*  plane-of-sky image is too small to "contain" all of the sky          */
////			/*  rendering's nonzero pixels.                                          */
////
////			if (checkposet( pos->b, -n_pos, n_pos, -n_pos, n_pos,
////					resamp_x0, resamp_xwidth, resamp_y0, resamp_ywidth, resamp_angle,
////					&badposet_logfactor_view)) {
////				par->badposet = 1;
////				par->badposet_logfactor += frame->dof * badposet_logfactor_view
////						/ poset->nviews;
////			}
////
////			/*  If smearing is being modeled, include the plane-of-sky
////          calculations from this view in the summed results for this frame  */
////
////			if (poset->nviews > 1)
////				for (i=1; i<=ncol_fit; i++)
////					for (j=1; j<=nrow_fit; j++)
////						fit_store[i][j] += frame->fit.b[i][j];
////
////		}
////
////		/*  If smearing is being modeled, compute mean values over all views
////        for this frame and store them in the standard frame structure     */
////
////		if (poset->nviews > 1) {
////			for (i=1; i<=ncol_fit; i++)
////				for (j=1; j<=nrow_fit; j++)
////					frame->fit.b[i][j] = fit_store[i][j] / poset->nviews;
////			free_matrix( fit_store, 1, ncol_fit, 1, nrow_fit);
////		}
////
////
////	}  /* end loop over frames */
////}
////
////

__host__ void calc_lghtcrv_gpu32(
		struct par_t *dpar,
		struct mod_t *dmod,
		struct dat_t *ddat,
		struct vertices_t **verts,
		int s,
		int nframes,
		int nviews,
		unsigned char type,
		int lc_n,
		int nf,
		struct pos_t **pos,
		struct crvrend_t **rend,
		int *posn,
		int *bistatic,
		float *pxlpkm,
		float *overflow,
		int4 *xylim,
		int *outbndarr,
		double3 *so,
		double *u,
		cudaStream_t *cf_stream)
{
	int ncalc, c=0, n, nThreads, exclude_seen, f;//, bistatic_all=0;
	int nfplus = nframes+1; /* This is to accomodate the +1 start in lghtcrv */
	float3 orbit_off3;
	orbit_off3.x = orbit_off3.y = orbit_off3.z = 0.0;
	dim3 BLKpx[nfplus],THD,THD9,THD64,BLKfrm;
	THD.x = maxThreadsPerBlock; THD9.x = 9; THD64.x = 64;
	BLKfrm = floor((THD64.x - 1 + nframes)/THD64.x);
	ncalc = nframes;
	n = lc_n;
	int /*hbistatic[nfplus], */houtbndarr[nfplus],
		nThreadspx[nfplus], nThreadspx1[nfplus], hposn[nfplus];
	int4 hxylim[nfplus];
	int2 span[nfplus];

	/* Set shortcuts and copy pos->n back for all frames */
	cfs_set_lghtcrv_shortcuts_krnl32<<<BLKfrm,THD64>>>(ddat,	rend, pos, overflow,
			posn, s, nframes);
	checkErrorAfterKernelLaunch("cfs_set_lghtcrv_shortcuts_krnl2");
	gpuErrchk(cudaMemcpy(&hposn, posn, sizeof(int)* nfplus,
			cudaMemcpyDeviceToHost));

	/* Calculate launch parameters for later */
	for (f=1; f<=nframes; f++) {
		span[f].x = (2 * hposn[f] + 1);
		nThreadspx[f] =  span[f].x * span[f].x;
		BLKpx[f] = floor((THD.x - 1 + nThreadspx[f]) / THD.x);
	}
	/* Calculate model lightcurve values at each user-specified epochs x[i],
	 * with i=1,2,...,ncalc; these may/may not be the same as epochs t[i]
	 * (i=1,2,...,n) at which actual lightcurve observations were made.  */
	cfs_set_pos_ae_krnl<<<BLKfrm,THD64>>>(ddat, pos, s, nframes,type, 0, 1);

	for (f=1; f<=ncalc; f++) {
		/* Clear the POS-view to initialize */
		posclr_krnl<<<BLKpx[f],THD,0,cf_stream[f-1]>>>(pos, posn, f, FP64, 1);
	} checkErrorAfterKernelLaunch("cfs_set_pos_ae_streams and posclr_streams_krnl");

	/* Call routine posvis to get  facet number, scattering & incidence
	 * angle, distance toward Earth at center of each POS pixel; set posbnd
	 * parameter = 1 if any model portion extends beyond POS frame limits*/
	posvis_gpu32(dpar, dmod, ddat, pos, verts, orbit_off3, hposn, outbndarr,
			s, (nframes+1), 0, nf, 0, c, type, cf_stream, 0);

	/* Copy the posbnd flag returns for all frames to a host copy */
	gpuErrchk(cudaMemcpy(&houtbndarr, outbndarr, sizeof(int)*(nfplus),
			cudaMemcpyDeviceToHost));

	/* Now check the outbndarr for the posbnd flag for each frame */
	for (f=1; f<=nframes; f++) {
		if (houtbndarr[f]) {
			/* Call single-threaded kernel to set dpar->posbnd and dpar->posbnd_logfactor */
			cfs_set_posbnd_krnl<<<1,1,0,cf_stream[f-1]>>>(dpar, ddat, pos, s, f, type);
			houtbndarr[f]=0;
		}
	}

	posvis_gpu32(dpar, dmod, ddat, pos, verts, orbit_off3, hposn,
			outbndarr, s, (nframes+1), 1, nf, 0, c, type, cf_stream, 1);
	/* Copy the posbnd flag returns for all frames to a host copy */
	gpuErrchk(cudaMemcpy(&houtbndarr, outbndarr, sizeof(int)*nfplus,
			cudaMemcpyDeviceToHost));

	/* Now check the outbndarr for the posbnd flag for each frame */
	for (f=1; f<=nframes; f++) {
		if (houtbndarr[f]) {
			/* Call single-threaded kernel to set dpar->posbnd and dpar->posbnd_logfactor */
			cfs_set_posbnd_krnl<<<1,1,0,cf_stream[f-1]>>>(dpar, ddat, pos, s, f, type);
		}
	}
	/* Initialize this stream for the posmask kernel to follow */
	posmask_init_krnl32<<<BLKfrm,THD64>>>(pos, so, pxlpkm, nframes);

	for (f=1; f<=nframes; f++) {
			/* Now call posmask kernel for this stream, then loop
			 * to next stream and repeat 	 */
			posmask_krnl32<<<BLKpx[f],THD,0,cf_stream[f-1]>>>(
					dpar, pos, so, pxlpkm, posn, nThreadspx[f],	span[f].x, f);

	} checkErrorAfterKernelLaunch("posmask_krnl");

	/* Go through all visible and unshadowed POS pixels with low enough
	 * scattering and incidence angles, and mark facets which project onto
	 * their centers as having been "seen" at least once. Get xlim and ylim
	 * and exclude_seen flag */
	cfs_get_exclude_seen_krnl<<<BLKfrm,THD64>>>(dpar,pos,xylim,nframes, 1);
	checkErrorAfterKernelLaunch("cfs_get_exclude_seen_krnl2 (calc_lghtcrv_gpu2");

	/* Now copy the flag and all frame pos's xlim and ylim values back from GPU */
	gpuErrchk(cudaMemcpyFromSymbol(&exclude_seen, cfs_exclude_seen, sizeof(int),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(&hxylim, xylim, sizeof(int4)*nfplus, cudaMemcpyDeviceToHost));

	/* Calculate launch parameters for all frames */
	int xspan_max = 0, yspan_max = 0, maxthds = 0;
	int4 maxxylim;
	maxxylim.w = maxxylim.y = 1e3;
	maxxylim.x = maxxylim.z = -1e3;
	for (f=1; f<=nframes; f++) {
		span[f].x = hxylim[f].x - hxylim[f].w + 1;
		span[f].y = hxylim[f].z - hxylim[f].y + 1;
		nThreadspx1[f] = span[f].x * span[f].y;
		BLKpx[f].x = floor ((THD.x -1 + nThreadspx1[f]) / THD.x);
		maxxylim.w = min(maxxylim.w, hxylim[f].w);
		maxxylim.x = max(maxxylim.x, hxylim[f].x);
		maxxylim.y = min(maxxylim.y, hxylim[f].y);
		maxxylim.z = max(maxxylim.z, hxylim[f].z);

	}
	xspan_max = maxxylim.x - maxxylim.w + 1;
	yspan_max = maxxylim.z - maxxylim.y + 1;
	maxthds = xspan_max * yspan_max;


	for (f=1; f<=nframes; f++) {
		/* Go through all POS pixels which are visible with low enough
		 * scattering angle and mark the facets which project onto
		 * their centers as having been "seen" at least once           */
		if (s != exclude_seen)
			cf_mark_pixels_seen_krnl32<<<BLKpx[f],THD,0,cf_stream[f-1]>>>(
					dpar, dmod, pos, xylim, nThreadspx1[f], span[f].x, f);
	} checkErrorAfterKernelLaunch("cf_mark_pixels_krnl (calc_lghtcrv_cuda)");

//	/* Compute model brightness for this lightcurve point then copy to device  */
//	apply_photo_gpu32(dmod, ddat, pos, xylim, span, BLKpx, nThreadspx,
//			0, s, nframes, nThreadspx1, maxthds, maxxylim, cf_stream);

	if (s==6) {
		int debug = 0;
		int frm = 5;
		if (debug)
			dbg_print_pos_arrays_full32(pos, frm,	nThreadspx[frm], hposn[frm]);
	}

	apply_photo_gpu48(dmod, ddat, pos, xylim, span, BLKpx, 0, s, nframes, maxthds, maxxylim, cf_stream);


//	dbg_print_facet_normals(dmod, nf, "FP32_facet_normals.csv");

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
//	if (s==6)
//		printf("\n6\n");

	nThreads = nframes;
	THD.x = maxThreadsPerBlock;
	gpuErrchk(cudaMemset(u, 0, sizeof(double)*nfplus));
	BLKpx[0].x = floor((THD.x - 1 + nThreads)/THD.x);
	lghtcrv_spline_krnl<<<1,1>>>(ddat, s, 2.0e30, 2.0e30, u, nframes);
	checkErrorAfterKernelLaunch("lghtcrv_spline_krnl");

	/* Change launch parameters from ncalc threads to n threads */
	BLKpx[0].x = floor((THD.x - 1 + n) / THD.x);
	lghtcrv_splint_krnl<<<1,1>>>(ddat, s, n, nframes);
	checkErrorAfterKernelLaunch("lghtcrv_splint_streams_krnl");
}

__host__ void calc_lghtcrv_gpu64(
		struct par_t *dpar,
		struct mod_t *dmod,
		struct dat_t *ddat,
		struct vertices_t **verts,
		int s,
		int nframes,
		int nviews,
		unsigned char type,
		int lc_n,
		int nf,
		struct pos_t **pos,
		struct crvrend_t **rend,
		int *posn,
		int *bistatic,
		double *pxlpkm,
		double *overflow,
		int4 *xylim,
		int *outbndarr,
		double3 *so,
		double *u,
		cudaStream_t *cf_stream)
{
	int ncalc, c=0, n, nThreads, exclude_seen, f;//, bistatic_all=0;
	int nfplus = nframes+1; /* This is to accomodate the +1 start in lghtcrv */
	double3 orbit_off3;
	orbit_off3.x = orbit_off3.y = orbit_off3.z = 0.0;
	dim3 BLKpx[nfplus],THD,THD9,THD64,BLKfrm;
	THD.x = maxThreadsPerBlock; THD9.x = 9; THD64.x = 64;
	BLKfrm = floor((THD64.x - 1 + nframes)/THD64.x);
	ncalc = nframes;
	n = lc_n;
	int /*hbistatic[nfplus], */houtbndarr[nfplus],
		nThreadspx[nfplus], nThreadspx1[nfplus], hposn[nfplus];
	int4 hxylim[nfplus];
	int2 span[nfplus];

	/* Set shortcuts and copy pos->n back for all frames */
	cfs_set_lghtcrv_shortcuts_krnl64<<<BLKfrm,THD64>>>(ddat, rend, pos, overflow,
			posn, s, nframes, outbndarr);
	checkErrorAfterKernelLaunch("cfs_set_lghtcrv_shortcuts_krnl32");
	gpuErrchk(cudaMemcpy(&hposn, posn, sizeof(int)* nfplus,
			cudaMemcpyDeviceToHost));

	/* Calculate launch parameters for later */
	for (f=1; f<=nframes; f++) {
		span[f].x = (2 * hposn[f] + 1);
		nThreadspx[f] =  span[f].x * span[f].x;
		BLKpx[f] = floor((THD.x - 1 + nThreadspx[f]) / THD.x);
	}
	/* Calculate model lightcurve values at each user-specified epochs x[i],
	 * with i=1,2,...,ncalc; these may/may not be the same as epochs t[i]
	 * (i=1,2,...,n) at which actual lightcurve observations were made.  */
	cfs_set_pos_ae_krnl<<<BLKfrm,THD64>>>(ddat, pos, s, nframes,type, 0, 1);

	for (f=1; f<=ncalc; f++) {
		/* Clear the POS-view to initialize */
		posclr_krnl<<<BLKpx[f],THD,0,cf_stream[f-1]>>>(pos, posn, f, FP64, 1);
	} checkErrorAfterKernelLaunch("cfs_set_pos_ae and posclr_krnl");

	/* Call routine posvis to get  facet number, scattering & incidence
	 * angle, distance toward Earth at center of each POS pixel; set posbnd
	 * parameter = 1 if any model portion extends beyond POS frame limits*/
	posvis_gpu64(dpar, dmod, ddat, pos, verts, orbit_off3, hposn, outbndarr,
			s, (nframes+1), 0, nf, 0, c, type, cf_stream, 0);

	/* Copy the posbnd flag returns for all frames to a host copy */
	gpuErrchk(cudaMemcpy(&houtbndarr, outbndarr, sizeof(int)*(nfplus),
			cudaMemcpyDeviceToHost));

	/* Now check the outbndarr for the posbnd flag for each frame */
	for (f=1; f<=nframes; f++) {
		if (houtbndarr[f]) {
			/* Call single-threaded kernel to set dpar->posbnd and dpar->posbnd_logfactor */
			cfs_set_posbnd_krnl<<<1,1,0,cf_stream[f-1]>>>(dpar, ddat, pos, s, f, type);
			houtbndarr[f]=0;
		}
	}

	/* Now view model from source (sun) and get facet number and distance
	 * toward source of each pixel in this projected view; use this
	 * information to determine which POS pixels are shadowed       */
	/* Because posvis_gpu processes all frames at the same time, if
	 * any of the frames are bistatic, all of them get calculated again  */
	posvis_gpu64(dpar, dmod, ddat, pos, verts, orbit_off3, hposn,
			outbndarr, s, (nframes+1), 1, nf, 0, c, type, cf_stream, 1);

	/* Copy the posbnd flag returns for all frames to a host copy */
	gpuErrchk(cudaMemcpy(&houtbndarr, outbndarr, sizeof(int)*nfplus,
			cudaMemcpyDeviceToHost));

	/* Now check the outbndarr for the posbnd flag for each frame */
	for (f=1; f<=nframes; f++) {
		if (houtbndarr[f]) {
			/* Call single-threaded kernel to set dpar->posbnd and dpar->posbnd_logfactor */
			cfs_set_posbnd_krnl<<<1,1,0,cf_stream[f-1]>>>(dpar, ddat, pos, s, f, type);
		}
	}
	/* Initialize this stream for the posmask kernel to follow */
	posmask_init_krnl64<<<BLKfrm,THD64>>>(pos, so, pxlpkm, nframes);

	for (f=1; f<=nframes; f++) {
			/* Now call posmask kernel for this stream, then loop
			 * to next stream and repeat 	 */
			posmask_krnl64<<<BLKpx[f],THD,0,cf_stream[f-1]>>>(
					dpar, pos, so, pxlpkm, posn, nThreadspx[f],	span[f].x, f);
	} checkErrorAfterKernelLaunch("posmask_krnl64");

	/* Go through all visible and unshadowed POS pixels with low enough
	 * scattering and incidence angles, and mark facets which project onto
	 * their centers as having been "seen" at least once. Get xlim and ylim
	 * and exclude_seen flag */
	cfs_get_exclude_seen_krnl<<<BLKfrm,THD64>>>(dpar,pos,xylim,nframes, 1);
	checkErrorAfterKernelLaunch("cfs_get_exclude_seen_krnl");

	/* Now copy the flag and all frame pos's xlim and ylim values back from GPU */
	gpuErrchk(cudaMemcpyFromSymbol(&exclude_seen, cfs_exclude_seen, sizeof(int),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(&hxylim, xylim, sizeof(int4)*nfplus, cudaMemcpyDeviceToHost));

	/* Calculate launch parameters for all frames */
	int xspan_max = 0, yspan_max = 0, maxthds = 0;
	int4 maxxylim;
	maxxylim.w = maxxylim.y = 1e3;
	maxxylim.x = maxxylim.z = -1e3;
	for (f=1; f<=nframes; f++) {
		span[f].x = hxylim[f].x - hxylim[f].w + 1;
		span[f].y = hxylim[f].z - hxylim[f].y + 1;
		nThreadspx1[f] = span[f].x * span[f].y;
		BLKpx[f].x = floor ((THD.x -1 + nThreadspx1[f]) / THD.x);
		maxxylim.w = min(maxxylim.w, hxylim[f].w);
		maxxylim.x = max(maxxylim.x, hxylim[f].x);
		maxxylim.y = min(maxxylim.y, hxylim[f].y);
		maxxylim.z = max(maxxylim.z, hxylim[f].z);

	}
	xspan_max = maxxylim.x - maxxylim.w + 1;
	yspan_max = maxxylim.z - maxxylim.y + 1;
	maxthds = xspan_max * yspan_max;

	for (f=1; f<=nframes; f++) {
		/* Go through all POS pixels which are visible with low enough
		 * scattering angle and mark the facets which project onto
		 * their centers as having been "seen" at least once           */
		if (s != exclude_seen)
			cf_mark_pixels_seen_krnl64<<<BLKpx[f],THD,0,cf_stream[f-1]>>>(
					dpar, dmod, pos, xylim, nThreadspx1[f], span[f].x, f);
	} checkErrorAfterKernelLaunch("cf_mark_pixels_krnl");

	/* Compute model brightness for this lightcurve point then copy to device  */
	apply_photo_gpu64(dmod, ddat, pos, xylim, span, BLKpx, nThreadspx,
			0, s, nframes, nThreadspx1, maxthds, maxxylim, cf_stream);
//	dbg_print_pos_arrays_full64(pos, 1, nThreadspx[1], hposn[1]);
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
	nThreads = nframes;
	THD.x = maxThreadsPerBlock;
	gpuErrchk(cudaMemset(u, 0, sizeof(double)*nfplus));
	BLKpx[0].x = floor((THD.x - 1 + nThreads)/THD.x);
	lghtcrv_spline_krnl<<<1,1>>>(ddat, s, 2.0e30, 2.0e30, u, nframes);
	checkErrorAfterKernelLaunch("lghtcrv_spline_krnl");

	/* Change launch parameters from ncalc threads to n threads */
	BLKpx[0].x = floor((THD.x - 1 + n) / THD.x);
	lghtcrv_splint_krnl<<<1,1>>>(ddat, s, n, nframes);
	checkErrorAfterKernelLaunch("lghtcrv_splint_krnl");
}
