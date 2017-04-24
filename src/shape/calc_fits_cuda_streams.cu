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

__host__ void calc_deldop_cuda_streams(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, struct vertices_t **verts, int s, int nframes, int
		nviews, unsigned char type, int nf, cudaStream_t *cf_stream);
__host__ void calc_doppler_cuda_streams(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, struct vertices_t **verts, int s, int nframes, int
		nviews, unsigned char type, int nf, cudaStream_t *cf_stream);
////__host__ void calc_poset_cuda( struct par_t *par, struct mod_t *mod, int s);
__host__ void calc_lghtcrv_cuda_streams(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, struct vertices_t **verts, int s, int nframes, int
		nviews, unsigned char type, int lc_n, int nf, cudaStream_t *cf_stream);
__host__ void calc_lghtcrv_cuda_streams2(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, struct vertices_t **verts, int s, int nframes,
		int nviews, unsigned char type, int lc_n, int nf, cudaStream_t *cf_stream);
__host__ void calc_lghtcrv_cuda_streams2f(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, struct vertices_t **verts, int s, int nframes,
		int nviews, unsigned char type, int lc_n, int nf, cudaStream_t *cf_stream);

__device__ int cfs_nf, cfs_nsets, cfs_v0_index, cfs_exclude_seen, cfs_ncalc, cfs_n;
__device__ double cfs_lghtcrv_posbnd_logfactor;
__device__ struct deldop_t *cfs_deldop;
__device__ struct doppler_t *cfs_doppler;
__device__ struct lghtcrv_t *cfs_lghtcrv;

__device__ void dev_spline_cfs(double *x,double *y,int n,double yp1,double ypn,double *y2, double *u)
{
	int i,k;
	double p,qn,sig,un;

	if (yp1 > 0.99e30)
		y2[1]=u[1]=0.0;
	else {
		y2[1] = -0.5;
		u[1] = (3.0 / (x[2]-x[1])) * ((y[2]-y[1]) / (x[2]-x[1])-yp1);
	}

	for (i=2;i<=n-1;i++) {
		sig = (x[i]-x[i-1]) / (x[i+1]-x[i-1]);
		p = sig * y2[i-1] + 2.0;
		y2[i] = (sig-1.0) / p;
		u[i] = (y[i+1]-y[i]) / (x[i+1]-x[i]) - (y[i]-y[i-1]) / (x[i]-x[i-1]);
		u[i] = (6.0 * u[i]/(x[i+1]-x[i-1]) - sig * u[i-1]) / p;
	}

	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
	}
	y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);

	for (k=n-1;k>=1;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
}
__device__ void dev_spline_cfs2(double *x,double *y,int n,double yp1,double ypn,double *y2, float *u)
{
	int i,k;
	float p,qn,sig,un;

	if (yp1 > 0.99e30)
		y2[1]=u[1]=0.0;
	else {
		y2[1] = -0.5;
		u[1] = (3.0 / (x[2]-x[1])) * ((y[2]-y[1]) / (x[2]-x[1])-yp1);
	}

	for (i=2;i<=n-1;i++) {
		sig = (x[i]-x[i-1]) / (x[i+1]-x[i-1]);
		p = sig * y2[i-1] + 2.0;
		y2[i] = (sig-1.0) / p;
		u[i] = (y[i+1]-y[i]) / (x[i+1]-x[i]) - (y[i]-y[i-1]) / (x[i]-x[i-1]);
		u[i] = (6.0 * u[i]/(x[i+1]-x[i-1]) - sig * u[i-1]) / p;
	}

	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
	}
	y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);

	for (k=n-1;k>=1;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
}
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
__global__ void cfs_init_devpar_krnl(struct par_t *dpar, struct mod_t
		*dmod, struct dat_t *ddat, struct vertices_t **verts) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		dpar->posbnd = 0;
		dpar->badposet = 0;
		dpar->badradar = 0;
		dpar->posbnd_logfactor = 0.0;
		dpar->badposet_logfactor = 0.0;
		dpar->badradar_logfactor = 0.0;
		verts[0] = &dmod->shape.comp[0].real;
		cfs_nf = verts[0]->nf;
		cfs_nsets = ddat->nsets;
	}
}
__global__ void cfs_init_devpar_krnl2(struct par_t *dpar) {
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
__global__ void cfs_get_set_type_krnl(struct dat_t *ddat, int nsets,
		unsigned char *type, int *nframes, int *nviews, int *lc_n) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		for (int s=0; s<nsets; s++) {
			type[s] = ddat->set[s].type;
			switch (type[s]) {
			case DELAY:
				nframes[s] = ddat->set[s].desc.deldop.nframes;
				nviews[s] = ddat->set[s].desc.deldop.nviews;
				lc_n[s] = 0;
				break;
			case DOPPLER:
				nframes[s] = ddat->set[s].desc.doppler.nframes;
				nviews[s] = ddat->set[s].desc.doppler.nviews;
				lc_n[s] = 0;
				break;
			case POS:
				nframes[s] = ddat->set[s].desc.poset.nframes;
				nviews[s] = ddat->set[s].desc.poset.nviews;
				lc_n[s] = 0;
				break;
			case LGHTCRV:
				nframes[s] = ddat->set[s].desc.lghtcrv.ncalc;
				nviews[s] = ddat->set[s].desc.lghtcrv.nviews;
				lc_n[s] = ddat->set[s].desc.lghtcrv.n;
				break;
			}
		}
	}
}

__host__ void calc_fits_cuda_streams(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat)
{
	int s, nf, nsets, f, *nframes, *nviews, *lc_n;
	unsigned char *type;
	dim3 BLK,THD;
	THD.x = maxThreadsPerBlock;
	struct vertices_t **verts;
	gpuErrchk(cudaMalloc((void**)&verts, sizeof(struct vertices_t*)*2));

	/* Initialize flags that indicate the model extends beyond POS frame, that
	 * plane-of-sky fit images are too small to "contain" the target, and that
	 * model is too wide in (delay-)Doppler space to create (delay-)Doppler fit
	 * frames.  Note that this also gets mod->shape.nf and nsets            */

	cfs_init_devpar_krnl<<<1,1>>>(dpar, dmod, ddat, verts);
	checkErrorAfterKernelLaunch("cfs_init_devpar_krnl");
	gpuErrchk(cudaMemcpyFromSymbol(&nf, cfs_nf, sizeof(int),
				0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&nsets, cfs_nsets, sizeof(int),
					0, cudaMemcpyDeviceToHost));

	/* Allocate temporary constructs (host and device) */
	unsigned char htype[nsets];
	int hnframes[nsets], hnviews[nsets], n[nsets];
	gpuErrchk(cudaMalloc((void**)&type, sizeof(unsigned char)*(nsets+1)));
	gpuErrchk(cudaMalloc((void**)&nframes, sizeof(int)*(nsets+1)));
	gpuErrchk(cudaMalloc((void**)&lc_n, sizeof(int)*(nsets+1)));
	gpuErrchk(cudaMalloc((void**)&nviews, sizeof(int)*(nsets+1))); // +1 is a safety margin

	/* Initialize the flags that indicate whether or not each facet of each
	 * model component is ever visible and unshadowed from Earth
	 * Note:  Single component only for now.  */
	//for (c=0; c<mod->shape.ncomp; c++)
	BLK.x = floor((THD.x - 1 + nf)/THD.x);
	cf_init_seen_flags_krnl<<<BLK,THD>>>(dmod,nf);
	checkErrorAfterKernelLaunch("cf_init_seen_flags_krnl (calc_fits_cuda_streams)");

	/* Get type, nframes, nviews */
	cfs_get_set_type_krnl<<<1,1>>>(ddat, nsets, type, nframes, nviews, lc_n);
	checkErrorAfterKernelLaunch("cf_init_seen_flags_krnl (calc_fits_cuda)");
	gpuErrchk(cudaMemcpy(&htype, type, sizeof(unsigned char)*nsets,
			cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(&hnframes, nframes, sizeof(int)*nsets,
			cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(&hnviews, nviews, sizeof(int)*nsets,
			cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(&n, lc_n, sizeof(int)*nsets,
			cudaMemcpyDeviceToHost));

	/* Calculate the fits for each dataset in turn  */
	for (s=0; s<nsets; s++) {

		/* Now create the streams we need */
		cudaStream_t cf_stream[hnframes[s]];
		for (f=0; f<hnframes[s]; f++)
			cudaStreamCreate(&cf_stream[f]);

		switch (htype[s]) {
		case DELAY:
			calc_deldop_cuda_streams(dpar, dmod, ddat, verts, s, hnframes[s],
					hnviews[s], htype[s], nf, cf_stream );
			break;
		case DOPPLER:
			calc_doppler_cuda_streams(dpar, dmod, ddat, verts, s, hnframes[s],
					hnviews[s], htype[s], nf, cf_stream );
			break;
		case POS:
			printf("Write calc_poset_cuda!");
//			calc_poset_cuda(dpar, dmod, s);
			break;
		case LGHTCRV:
			calc_lghtcrv_cuda_streams(dpar, dmod, ddat, verts, s, hnframes[s],
					hnviews[s], htype[s], n[s],nf, cf_stream);
			break;
		default:
			printf("calc_fits_cuda.c: can't handle this type yet\n");
		}
		/* Lastly, destroy streams */
		for (f=0; f<hnframes[s]; f++)
			cudaStreamDestroy(cf_stream[f]);
	}

	/* Complete calculations of values that will be used during a fit to
	 * increase the objective function for models with bad properties   */
	cf_set_final_pars_krnl<<<1,1>>>(dpar, ddat);
	checkErrorAfterKernelLaunch("cf_set_final_pars_krnl (calc_fits_cuda)");

	/* Free temporary memory  */
	//free(htype);
	//free(hnframes);
	//free(hnviews);
	cudaFree(type);
	cudaFree(nframes);
	cudaFree(nviews);
	cudaFree(lc_n);
	cudaFree(verts);
}
__host__ void calc_fits_cuda_streams2(
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
		cudaStream_t *cf_stream)
{
	int s;
	dim3 BLK,THD;
	THD.x = maxThreadsPerBlock;

	/* Initialize flags that indicate the model extends beyond POS frame, that
	 * plane-of-sky fit images are too small to "contain" the target, and that
	 * model is too wide in (delay-)Doppler space to create (delay-)Doppler fit
	 * frames.  Note that this also gets mod->shape.nf and nsets            */

	cfs_init_devpar_krnl2<<<1,1>>>(dpar);
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
		switch (type[s]) {
		case DELAY:
			calc_deldop_cuda_streams(dpar, dmod, ddat, verts, s, nframes[s],
					nviews[s], type[s], nf, cf_stream );
			break;
		case DOPPLER:
			calc_doppler_cuda_streams(dpar, dmod, ddat, verts, s, nframes[s],
					nviews[s], type[s], nf, cf_stream );
			break;
		case POS:
			printf("Write calc_poset_cuda!");
//			calc_poset_cuda(dpar, dmod, s);
			break;
		case LGHTCRV:
			if (FLOAT)
				calc_lghtcrv_cuda_streams2f(dpar, dmod, ddat, verts, s, nframes[s],
						nviews[s], type[s], lc_n[s], nf, cf_stream);
			else
				calc_lghtcrv_cuda_streams2(dpar, dmod, ddat, verts, s, nframes[s],
					nviews[s], type[s], lc_n[s],nf, cf_stream);
			break;
		default:
			printf("calc_fits_cuda.c: can't handle this type yet\n");
		}
	}
//
//
//	int debug = 0;
//	int size;
//	if (debug)
//		size = (2*50+1)*(2*50+1);
//	dbg_print_lghtcrv_pos_arrays(ddat, 0, 22, size, 50);


	/* Complete calculations of values that will be used during a fit to
	 * increase the objective function for models with bad properties   */
	cf_set_final_pars_krnl<<<1,1>>>(dpar, ddat);
	checkErrorAfterKernelLaunch("cf_set_final_pars_krnl (calc_fits_cuda)");

}
__global__ void cfs_set_deldop_shortcuts_krnl(struct dat_t *ddat,
		struct deldopfrm_t **frame, struct pos_t **pos,
		struct deldopview_t **view0, int *ndop, int *ndel, float *overflow,
		int *posn, int s, int f) {
	/* Single-threaded kernel */
	/* 	overflow[0] - overflow_o2_store
	 * 	overflow[1] - overflow_m2_store
	 * 	overflow[2] - overflow_xsec_store
	 * 	overflow[3] - overflow_dopmean_store
	 * 	overflow[4] - overflow_delmean_store
	 */
	if (threadIdx.x == 0) {
		if (f == 0)
			cfs_deldop 	 = &ddat->set[s].desc.deldop;

		frame[f] 	  = &cfs_deldop->frame[f];
		ndop[f]  	  = frame[f]->ndop;
		ndel[f]	 	  = frame[f]->ndel;
		view0[f] 	  = &frame[f]->view[cfs_deldop->v0];
		cfs_v0_index  = cfs_deldop->v0;
		pos[f]		  = &frame[f]->pos;
		posn[f]		  = pos[f]->n;
		overflow[0] = 0.0;
		overflow[1] = 0.0;
		overflow[2] = 0.0;
		overflow[3] = 0.0;
		overflow[4] = 0.0;

	}
}
__global__ void cfs_set_doppler_shortcuts_krnl(struct dat_t *ddat,
		struct dopfrm_t **frame, struct pos_t **pos, struct dopview_t **view0,
		float *overflow, int *ndop, int *posn, int s, int f) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		if (f == 0)
			cfs_doppler = &ddat->set[s].desc.doppler;

		frame[f] = &cfs_doppler->frame[f];
		view0[f] = &frame[f]->view[cfs_doppler->v0];
		cfs_v0_index  = cfs_doppler->v0;
		ndop[f]	 = frame[f]->ndop;
		pos[f]		 = &frame[f]->pos;
		posn[f]	= pos[f]->n;
		overflow[0] = 0.0;
		overflow[1] = 0.0;
		overflow[2] = 0.0;
		overflow[3] = 0.0;
		overflow[4] = 0.0;
	}
}
__global__ void cfs_set_lghtcrv_shortcuts_streams_krnl(struct dat_t *ddat,
		struct crvrend_t **rend, struct pos_t **pos, float *overflow,
		int *posn, int s, int f) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		if (f == 1)
			cfs_lghtcrv = &ddat->set[s].desc.lghtcrv;

		rend[f] = &ddat->set[s].desc.lghtcrv.rend[f];//&cfs_lghtcrv->rend[f];
		pos[f] = &rend[f]->pos;
		posn[f] = pos[f]->n;
		overflow[0] = 0.0;
		overflow[1] = 0.0;
		overflow[2] = 0.0;
		overflow[3] = 0.0;
		overflow[4] = 0.0;
	}
}
__global__ void cfs_spline_lghtcrv_krnl(double yp1, double ypn, double *u) {
	/* ncalc-threaded kernel */
	int k, i = blockIdx.x * blockDim.x + threadIdx.x + 1;
	double p, qn, sig, un;
	int n = cfs_ncalc;

	/* Single-threaded task first */
	if (i == 1) {
		if (yp1 > 0.99e30)
			cfs_lghtcrv->y2[1] = u[1] = 0.0;
		else {
			cfs_lghtcrv->y2[1] = -0.5;
			u[1] = (3.0 / (cfs_lghtcrv->x[2] - cfs_lghtcrv->x[1])) *
				   ((cfs_lghtcrv->y[2] - cfs_lghtcrv->y[1]) /
				    (cfs_lghtcrv->x[2] - cfs_lghtcrv->x[1]) - yp1);
		}
	}
	__syncthreads();

	if ((i > 1) && (i <= (n-1))) {
		sig = (cfs_lghtcrv->x[i]   - cfs_lghtcrv->x[i-1]) /
			  (cfs_lghtcrv->x[i+1] - cfs_lghtcrv->x[i-1]);

		p = sig * cfs_lghtcrv->y2[i-1] + 2.0;

		cfs_lghtcrv->y2[i] = (sig - 1.0) / p;

		u[i] = (cfs_lghtcrv->y[i+1] - cfs_lghtcrv->y[i]) / (cfs_lghtcrv->x[i+1] -
				cfs_lghtcrv->x[i]) - (cfs_lghtcrv->y[i]  -  cfs_lghtcrv->y[i-1]) /
		   	   (cfs_lghtcrv->x[i]  -  cfs_lghtcrv->x[i-1]);

		u[i] = (6.0 *u[i] / (cfs_lghtcrv->x[i+1] - cfs_lghtcrv->x[i-1]) -
				sig * u[i-1]) / p;
	}
	__syncthreads();

	/* Another single-threaded task */
	if (i == 1) {
		if (ypn > 0.99e30)
			qn = un = 0.0;
		else {
			qn = 0.5;
			un = (3.0 / (cfs_lghtcrv->x[n] - cfs_lghtcrv->x[n-1])) * (ypn -
						(cfs_lghtcrv->y[n] - cfs_lghtcrv->y[n-1]) /
						(cfs_lghtcrv->x[n] - cfs_lghtcrv->x[n-1]));
		}
		cfs_lghtcrv->y2[n]=(un - qn * u[n-1]) /
				(qn * cfs_lghtcrv->y2[n-1] + 1.0);

		for (k=n-1; k>=1; k--)
			cfs_lghtcrv->y2[k] = cfs_lghtcrv->y2[k] * cfs_lghtcrv->y2[k+1] + u[k];
	}
	__syncthreads();
}
__global__ void cfs_spline_lghtcrv_serial_krnl(double *u) {
	/* single-threaded kernel */
	if (threadIdx.x == 0)
		dev_spline_cfs( cfs_lghtcrv->x, cfs_lghtcrv->y, cfs_ncalc, 2.0e30, 2.0e30, cfs_lghtcrv->y2, u);
}
__global__ void cfs_spline_lghtcrv_serial2_krnl(float *u) {
	/* single-threaded kernel */
	if (threadIdx.x == 0)
		dev_spline_cfs2( cfs_lghtcrv->x, cfs_lghtcrv->y, cfs_ncalc, 2.0e30, 2.0e30, cfs_lghtcrv->y2, u);
}
__global__ void cfs_splint_lghtcrv_krnl(struct par_t *dpar) {
	/* ncalc-threaded kernel */
	int v, i = blockIdx.x * blockDim.x + threadIdx.x + 1;
	double interp;

	if ((i >= 1) && (i <= cfs_lghtcrv->n)) {

		cfs_lghtcrv->fit[i] = 0.0;

		for (v=0; v<cfs_lghtcrv->nviews; v++) {
			dev_splint_cfs(cfs_lghtcrv->x, cfs_lghtcrv->y, cfs_lghtcrv->y2, cfs_ncalc,
					cfs_lghtcrv->t[i][v], &interp);
			cfs_lghtcrv->fit[i] += interp;
		}
		cfs_lghtcrv->fit[i] /= cfs_lghtcrv->nviews;
	}
	__syncthreads();

	/* Single-threaded task: */
	if (i == 1) {
		/* Deal with flags for model that extends beyond the POS frame  */
		dpar->posbnd_logfactor += cfs_lghtcrv->dof *
				(cfs_lghtcrv_posbnd_logfactor/cfs_ncalc);
	}
}
__global__ void cfs_set_pos_ae_streams_krnl(struct pos_t **pos, int f,
		int *bistatic, unsigned char type, int v) {
	/* 9-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % 3;
	int j = offset / 3;

	if (offset < 9) {
		switch (type) {
		case DELAY:
			pos[f]->ae[i][j] = cfs_deldop->frame[f].view[v].ae[i][j];
			pos[f]->oe[i][j] = cfs_deldop->frame[f].view[v].oe[i][j];
			break;
		case DOPPLER:
			pos[f]->ae[i][j] = cfs_doppler->frame[f].view[v].ae[i][j];
			pos[f]->oe[i][j] = cfs_doppler->frame[f].view[v].oe[i][j];
			break;
		case LGHTCRV:
			pos[f]->ae[i][j] = cfs_lghtcrv->rend[f].ae[i][j];
			pos[f]->oe[i][j] = cfs_lghtcrv->rend[f].oe[i][j];
			pos[f]->se[i][j] = cfs_lghtcrv->rend[f].se[i][j];
		}
		/* Single-thread task */
		if (offset == 0) {
			if ((type == LGHTCRV) || (type == POS))
				pos[f]->bistatic = 1;
			else if ((type == DELAY) || (type == DOPPLER))
				pos[f]->bistatic = 0;

			bistatic[f] = pos[f]->bistatic;
		}
	}
}
__global__ void cfs_set_posbnd_streams_krnl(struct par_t *dpar, struct pos_t **pos,
		int f, unsigned char type) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		dpar->posbnd = 1;
		switch (type) {
		case DELAY:
			dpar->posbnd_logfactor += cfs_deldop->frame[f].dof * pos[f]->posbnd_logfactor;
			break;
		case DOPPLER:
			dpar->posbnd_logfactor += cfs_doppler->frame[f].dof * pos[f]->posbnd_logfactor;
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
__global__ void cfs_get_exclude_seen_streams_krnl(struct par_t *dpar, struct pos_t **pos,
		int4 *xylim, int f) {
	/* single-threaded kernel */

	if (threadIdx.x == 0) {
		cfs_exclude_seen = dpar->exclude_seen;
		xylim[f].w = pos[f]->xlim[0];
		xylim[f].x = pos[f]->xlim[1];
		xylim[f].y = pos[f]->ylim[0];
		xylim[f].z = pos[f]->ylim[1];
	}
}
__global__ void cf_mark_pixels_seen_streams_krnl(struct par_t *dpar,
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
__global__ void cfs_set_lghtcrv_y_streams_krnl(double *brightness_temp,
		int i) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		cfs_lghtcrv->y[i] = brightness_temp[i];
	}
}
__global__ void cfs_set_lghtcrv_y_streams2_krnl(double *brightness_temp,
		int nfplus) {
	/* Single-threaded kernel */
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < nfplus) {
		cfs_lghtcrv->y[i] = brightness_temp[i];
	}
}
__global__ void cf_set_badradar_streams_krnl(struct par_t *dpar, int f,
		unsigned char type) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		switch (type) {
		case DELAY:
			dpar->badradar = 1;
			dpar->badradar_logfactor += cfs_deldop->frame[f].dof *
			cfs_deldop->frame[f].badradar_logfactor / cfs_deldop->nviews;
			break;
		case DOPPLER:
			dpar->badradar = 1;
			dpar->badradar_logfactor += cfs_doppler->frame[f].dof *
				cfs_doppler->frame[f].badradar_logfactor / cfs_doppler->nviews;
			break;
		}
	}
}
__global__ void cf_add_fit_store_streams_krnl1(struct dat_t *ddat, float **fit_store,
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
__global__ void cf_add_fit_store_streams_krnl2(float *overflow, int f,
		unsigned char type) {
	/* ndel*ndop-threaded kernel */
	if (threadIdx.x == 0) {
		switch (type) {
		case DELAY:
			atomicAdd(&overflow[0], (float)cfs_deldop->frame[f].overflow_o2);
			atomicAdd(&overflow[1], (float)cfs_deldop->frame[f].overflow_m2);
			atomicAdd(&overflow[2], (float)cfs_deldop->frame[f].overflow_xsec);
			atomicAdd(&overflow[3], (float)cfs_deldop->frame[f].overflow_delmean);
			atomicAdd(&overflow[4], (float)cfs_deldop->frame[f].overflow_dopmean);
			break;
		case DOPPLER:
			atomicAdd(&overflow[0], (float)cfs_doppler->frame[f].overflow_o2);
			atomicAdd(&overflow[1], (float)cfs_doppler->frame[f].overflow_m2);
			atomicAdd(&overflow[2], (float)cfs_doppler->frame[f].overflow_xsec);
			atomicAdd(&overflow[3], (float)cfs_doppler->frame[f].overflow_dopmean);
		}
	}
}
__global__ void cf_finish_fit_store_streams_krnl(float **fit_store,
		int s, int f, int nThreads, unsigned char type) {
	/* multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;

	if (offset < nThreads) {
		switch (type) {
		case DELAY:
			cfs_deldop->frame[f].fit_s[offset] = fit_store[f][offset];
			break;
		case DOPPLER:
			cfs_doppler->frame[f].fit_s[offset] = fit_store[f][offset];
			break;
		}
	}
}
__global__ void cf_finish_fit_streams_krnl2(float *overflow, int f, unsigned char type) {
	/* Single-threaded Kernel */
	int nviews;
	if (threadIdx.x == 0) {
		switch (type) {
		case DELAY:
			nviews = cfs_deldop->nviews;
			cfs_deldop->frame[f].overflow_o2 = overflow[0] / nviews;
			cfs_deldop->frame[f].overflow_m2 = overflow[1] / nviews;
			cfs_deldop->frame[f].overflow_xsec = overflow[2] / nviews;
			cfs_deldop->frame[f].overflow_delmean = overflow[3] / nviews;
			cfs_deldop->frame[f].overflow_dopmean = overflow[4] / nviews;
			break;
		case DOPPLER:
			nviews = cfs_doppler->nviews;
			cfs_doppler->frame[f].overflow_o2 = overflow[0] / nviews;
			cfs_doppler->frame[f].overflow_m2 = overflow[1] / nviews;
			cfs_doppler->frame[f].overflow_xsec = overflow[2] / nviews;
			cfs_doppler->frame[f].overflow_dopmean = overflow[3] / nviews;
		}

	}
}
__global__ void cf_gamma_trans_streams_krnl(struct par_t *dpar, struct dat_t *ddat,
		int s, int f, int nThreads, unsigned char type, int flt) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	if (offset < nThreads) {
		/*  Carry out a gamma transformation on the fit image if requested  */
		if (dpar->dd_gamma != 1.0) {
			switch (type) {
			case DELAY:
				if (flt)
					dev_gamma_trans_f(&ddat->set[s].desc.deldop.frame[f].fit_s[offset],
							__double2float_rn(dpar->dd_gamma));
				else
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

__host__ void calc_deldop_cuda_streams(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, struct vertices_t **verts, int s, int nframes,
		int nviews, unsigned char type, int nf, cudaStream_t *cf_stream)
{
	//double orbit_offset[3] = {0.0, 0.0, 0.0};
	float3 orbit_off3, orb_xydopoff;
	orbit_off3.x = orbit_off3.y = orbit_off3.z = 0.0;
	orb_xydopoff.x = orb_xydopoff.y = orb_xydopoff.z = 0.0;
	int *ndel, *ndop, *posn, *bistatic, v0_index, exclude_seen, f, v2, c=0, yspan,
			*outbndarr;
	float **fit_store, *overflow;
	dim3 BLKdd[nframes], BLKpx[nframes], THD, THD9;
	THD.x = maxThreadsPerBlock; THD9.x = 9;
	int4 *xylim, hxylim[nframes];
	struct pos_t **pos;
	struct deldopfrm_t **frame;
	struct deldopview_t **view0;
	int hndop[nframes], hndel[nframes], hposn[nframes], hbistatic[nframes],
		houtbndarr[nframes], xspan[nframes], nThreadspx[nframes], nThreadsdd[nframes],
		nThreadspx1[nframes], v[nviews+1];

	gpuErrchk(cudaMalloc((void**)&pos, sizeof(pos_t*) * nframes));
	gpuErrchk(cudaMalloc((void**)&frame, sizeof(deldopfrm_t*) * nframes));
	gpuErrchk(cudaMalloc((void**)&view0, sizeof(deldopview_t*) * nframes));
	gpuErrchk(cudaMalloc((void**)&fit_store, sizeof(float*) * nframes));
	gpuErrchk(cudaMalloc((void**)&ndel, sizeof(int) * nframes));
	gpuErrchk(cudaMalloc((void**)&ndop, sizeof(int) * nframes));
	gpuErrchk(cudaMalloc((void**)&posn, sizeof(int) * nframes));
	gpuErrchk(cudaMalloc((void**)&bistatic, sizeof(int) * nframes));
	gpuErrchk(cudaMalloc((void**)&xylim, sizeof(int4) * nframes));
	gpuErrchk(cudaMalloc((void**)&overflow, sizeof(float) * 6));
	gpuErrchk(cudaMalloc((void**)&outbndarr, sizeof(int) * nframes));

	for (f=0; f<nframes; f++)
		/* Set deldop, frame, view0, and pos in nframes streamed kernels */
		cfs_set_deldop_shortcuts_krnl<<<1,1,0,cf_stream[f]>>>(ddat, frame, pos,
				view0, ndop, ndel, overflow, posn, s, f);
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

	for (f=0; f<nframes; f++) {
		for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
			/* Launch 9-threaded kernel to set pos->ae,pos->oe,pos->bistatic.*/
			cfs_set_pos_ae_streams_krnl<<<1,THD9,0,cf_stream[f]>>>(pos, f,
					bistatic, type, v[v2]);

			/* Launch posclr_krnl to initialize POS view */
			posclr_streams_krnl<<<BLKpx[f],THD,0,cf_stream[f]>>>(pos, posn, f);

		}
	} checkErrorAfterKernelLaunch("posclr_streams_krnl (calc_fits_cuda_streams)");
	gpuErrchk(cudaMemcpy(&hbistatic, bistatic, sizeof(int)*nframes,
			cudaMemcpyDeviceToHost));

	/* Call posvis_cuda_2 to get facet number, scattering angle, distance
	 * toward Earth at center of each POS pixel; set flag posbnd if any model
	 * portion extends beyond POS frame limits.*/
	/* NOTE: Limited to single component for now */

//	for (v2=0; v2<no_views; v2++)
//	posvis_cuda_streams(dpar, dmod, ddat, orbit_offset, s, nframes,
//						0, 0, c, outbndarr, cf_stream);
	posvis_cuda_streams2(dpar, dmod, ddat, pos, verts, orbit_off3, hposn,
			outbndarr, s, nframes, 0, nf, 0, c, type, cf_stream);

	gpuErrchk(cudaMemcpy(&houtbndarr, outbndarr, sizeof(int)*nframes,
			cudaMemcpyDeviceToHost));

	for (f=0; f<nframes; f++) {
		for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
			if ((houtbndarr[f]) && (v[v2] == v0_index)) {
			/* Call single-threaded kernel to set dpar->posbnd and dpar->posbnd_logfactor */
			cfs_set_posbnd_streams_krnl<<<1,1,0,cf_stream[f]>>>(dpar, pos, f, type);
			houtbndarr[f]=0;
			}
		}
		/* Get xlim and ylim and exclude_seen flag */
		cfs_get_exclude_seen_streams_krnl<<<1,1,0,cf_stream[f]>>>(dpar,pos,xylim,f);
	} checkErrorAfterKernelLaunch("cfs_set_posbnd_streams_krnl and"
			"cfs_get_exclude_seen_streams_krnl");

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
				cf_mark_pixels_seen_streams_krnl<<<BLKpx[f],THD,0,cf_stream[f]>>>(
						dpar, dmod, pos, xylim, nThreadspx1[f], xspan[f], f);

			/* Zero out the fit delay-Doppler image, then call pos2deldop to
			 * create the fit image by mapping power from the plane of the sky
			 * to delay-Doppler space.                             */
			clrvect_krnl<<<BLKdd[f],THD,0,cf_stream[f]>>>(ddat, s, f, nThreadsdd[f]);
		} checkErrorAfterKernelLaunch("clrvect_krnl and cf_mark_pixels_seen_streams_krnl");
	}

	for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
		if (FLOAT)
			pos2deldop_cuda_streams_f(dpar,dmod,ddat, pos, ndel, ndop, orb_xydopoff,0,
					s, nframes, v[v2], outbndarr, cf_stream);
		else
			pos2deldop_cuda_streams(dpar,dmod,ddat, pos, ndel, ndop, 0.0,0.0,0.0,0,
					s, nframes, v[v2], outbndarr, cf_stream);
	}

	/* Copy the badradar flag returns for all frames to a host copy */
	gpuErrchk(cudaMemcpy(&houtbndarr, outbndarr, sizeof(int)*nframes,
			cudaMemcpyDeviceToHost));

	for (f=0; f<nframes; f++) {
		for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
			if (houtbndarr[f]) {
				/* Call single-threaded kernel to set badradar flag and
				 * associated badradar_logfactor			 */
				cf_set_badradar_streams_krnl<<<1,1,0,cf_stream[f]>>>(dpar, f, type);
				checkErrorAfterKernelLaunch("cf_deldop_set_badradar_krnl (calc_fits_cuda)");
			}
		}
	}
	for (f=0; f<nframes; f++) {
		/* If smearing is being modeled, include delay-Doppler calculations
		 * from this view in the summed results for this frame  */
		if (nviews > 1) {
			/* Launch ndel*ndop-threaded kernel to add fit[i][j] to
			 * fit_store[i][j]*/

			cf_add_fit_store_streams_krnl1<<<BLKdd[f],THD,0,cf_stream[f]>>>(
					ddat,fit_store,nThreadsdd[f],s,f, type);

			cf_add_fit_store_streams_krnl2<<<1,1>>>(overflow, f, type);
		}
	} checkErrorAfterKernelLaunch("cf_add_fit_store_streams_krnl1 and 2");

	/* If smearing is being modeled, compute mean values over all views for
	 * this frame and store them in the standard frame structure     */
	/* This kernel also carries out the gamma transformation on the fit
	 * image if the par->dd_gamma flag is not set  */
	if (nviews > 1) {
		for (f=0; f<nframes; f++) {

			cf_finish_fit_store_streams_krnl<<<BLKdd[f],THD,0,cf_stream[f]>>>(
					fit_store, s, f, nThreadsdd[f], type);

			cf_finish_fit_streams_krnl2<<<1,1,0,cf_stream[f]>>>(overflow, f, type);

			cf_gamma_trans_streams_krnl<<<BLKdd[f],THD,0,cf_stream[f]>>>(dpar, ddat, s, f, nThreadsdd[f], type, FLOAT);
		} checkErrorAfterKernelLaunch("cf_finish_fit_store_streams kernels and "
				"cf_gamma_trans_krnl");
		cudaFree(fit_store);
	}
	cudaFree(ndel);
	cudaFree(ndop);
	cudaFree(posn);
	cudaFree(bistatic);
	cudaFree(xylim);
	cudaFree(overflow);
	cudaFree(outbndarr);
}

__host__ void calc_doppler_cuda_streams(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, struct vertices_t **verts, int s, int nframes,
		int nviews, unsigned char type, int nf, cudaStream_t *cf_stream)
{
	//double orbit_offset[3] = {0.0, 0.0, 0.0};
	float3 orbit_off3;
	orbit_off3.x = orbit_off3.y = orbit_off3.z = 0.0;
	float **fit_store, *overflow;
	int *ndop, *posn, *bistatic, v0_index, exclude_seen, f, v2, c=0, yspan,
			*outbndarr;
	dim3 BLKd[nframes], BLKpx[nframes],THD,THD9;
	THD.x = maxThreadsPerBlock; THD9.x = 9;
	int4 *xylim, hxylim[nframes];
	struct pos_t **pos;
	struct dopfrm_t **frame;
	struct dopview_t **view0;
	int hndop[nframes], hposn[nframes], hbistatic[nframes], houtbndarr[nframes],
		xspan[nframes], nThreadspx[nframes], nThreadspx1[nframes], v[nviews+1];

	/* Allocate the CUDA pointers */
	cudaCalloc1((void**)&pos, 		sizeof(struct pos_t*), 		nframes);
	cudaCalloc1((void**)&frame, 		sizeof(struct dopfrm_t*), 	nframes);
	cudaCalloc1((void**)&view0, 		sizeof(struct dopview_t*), 	nframes);
	cudaCalloc1((void**)&fit_store,	sizeof(float*), 			nframes);
	cudaCalloc1((void**)&ndop, 		sizeof(int), 				nframes);
	cudaCalloc1((void**)&posn, 		sizeof(int), 				nframes);
	cudaCalloc1((void**)&bistatic, 	sizeof(int), 				nframes);
	cudaCalloc1((void**)&xylim, 		sizeof(int4), 				nframes);
	cudaCalloc1((void**)&overflow, 	sizeof(float), 				nframes);
	cudaCalloc1((void**)&outbndarr, 	sizeof(int), 				nframes);

	for (f=0; f<nframes; f++)
		/* Set doppler, frame, view0, and pos in nframes streamed kernels */
		cfs_set_doppler_shortcuts_krnl<<<1,1,0,cf_stream[f]>>>(ddat, frame,
				pos, view0, overflow, ndop, posn, s, f);
	checkErrorAfterKernelLaunch("cfs_set_doppler_shortcuts_krnl");
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
			cudaCalloc1((void**)&fit_store[f], sizeof(float), ndop[f]);
	}

	/* Loop over all views for this (smeared) frame, going in an order that
	 * ends with the view corresponding to the epoch listed for this frame
	 * in the obs file; this way we can use the calculated information for
	 * that view in the "write" action screen and disk output that follows*/

	for (v2=v0_index+1; v2<=v0_index+nviews; v2++)
		v[v2] = v2 % nviews;

	for (f=0; f<nframes; f++) {
		for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
			/* Launch 9-threaded kernel to set pos->ae,pos->oe,pos->bistatic.*/
			cfs_set_pos_ae_streams_krnl<<<1,THD9,0,cf_stream[f]>>>(pos, f,
					bistatic, type, v[v2]);

			/* Clear the POS-view to initialize */
			posclr_streams_krnl<<<BLKpx[f],THD,0,cf_stream[f]>>>(pos, posn, f);
		}
	} checkErrorAfterKernelLaunch("cfs_set_pos_ae_streams and posclr_streams_krnl");

	/* Call posvis_cuda_2 to get facet number, scattering angle, distance
	 * toward Earth at center of each POS pixel; set flag posbnd if any model
	 * portion extends beyond POS frame limits.*/
	/* NOTE: Limited to single component for now */
	//for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
	posvis_cuda_streams2(dpar, dmod, ddat, pos, verts, orbit_off3, hposn,
			outbndarr, s, nframes, 0, nf, 0, c, type, cf_stream);

	/* Copy the posbnd flag returns for all frames to a host copy */
	gpuErrchk(cudaMemcpy(&houtbndarr, outbndarr, sizeof(int)*nframes,
			cudaMemcpyDeviceToHost));

	for (f=0; f<nframes; f++) {
		for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
			if ((houtbndarr[f]) && (v[v2] == v0_index)) {
				/* Call single-threaded kernel to set dpar->posbnd and dpar->posbnd_logfactor */
				cfs_set_posbnd_streams_krnl<<<1,1,0,cf_stream[f]>>>(dpar, pos, f, type);
				houtbndarr[f]=0;
			}
		}
		/* Get xlim and ylim and exclude_seen flag */
		cfs_get_exclude_seen_streams_krnl<<<1,1,0,cf_stream[f]>>>(dpar,pos,xylim,f);
	} checkErrorAfterKernelLaunch("cfs_set_posbnd_streams_krnl and"
			"cfs_get_exclude_seen_streams_krnl (calc_doppler_cuda_streams");

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
				cf_mark_pixels_seen_streams_krnl<<<BLKpx[f],THD,0,cf_stream[f]>>>(
						dpar, dmod, pos, xylim, nThreadspx1[f], xspan[f], f);

			/* Zero out the fit delay-Doppler image, then call pos2deldop to
			 * create the fit image by mapping power from the plane of the sky
			 * to delay-Doppler space.                             */
			clrvect_krnl<<<BLKd[f],THD,0,cf_stream[f]>>>(ddat, s, f, hndop[f]);
		} checkErrorAfterKernelLaunch("clrvect_krnl and cf_mark_pixels_seen_streams_krnl");
	}
	/* Call pos2deldop to calculate the Doppler radar fit image */
	for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
		if (FLOAT)
			pos2doppler_cuda_streams_f(dpar, dmod, ddat, pos, orbit_off3,
					ndop, 0, s, nframes, v[v2],	outbndarr, cf_stream);
		else
			pos2doppler_cuda_streams(dpar, dmod, ddat, pos, 0.0, 0.0, 0.0,
					ndop, 0, s, nframes, v[v2],	outbndarr, cf_stream);
	}
	/* Copy the badradar flag returns for all frames to a host copy */
	gpuErrchk(cudaMemcpy(&houtbndarr, outbndarr, sizeof(int)*nframes,
			cudaMemcpyDeviceToHost));

	for (f=0; f<nframes; f++) {
		for (v2=v0_index+1; v2<=v0_index+nviews; v2++) {
			if (houtbndarr[f]) {
				/* Call single-threaded kernel to set badradar flag and
				 * associated badradar_logfactor			 */
				cf_set_badradar_streams_krnl<<<1,1,0,cf_stream[f]>>>(dpar, f, type);
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
			cf_add_fit_store_streams_krnl1<<<BLKd[f],THD,0,cf_stream[f]>>>(
					ddat, fit_store, ndop[f], s, f, type);

			cf_add_fit_store_streams_krnl2<<<1,1>>>(overflow, f, type);
		}
	} checkErrorAfterKernelLaunch("cf_add_fit_store_streams_krnl1 and 2 (calc_doppler_cuda_streams");

	/* If smearing is being modeled, compute mean values over all views for
	 * this frame and store them in the standard frame structure     */
	/* This kernel also carries out the gamma transformation on the fit
	 * image if the par->dd_gamma flag is not set  */
	if (nviews > 1) {
		for (f=0; f<nframes; f++) {

			cf_finish_fit_store_streams_krnl<<<BLKd[f],THD,0,cf_stream[f]>>>(
					fit_store, s, f, ndop[f], type);

			cf_finish_fit_streams_krnl2<<<1,1,0,cf_stream[f]>>>(overflow, f, type);

			cf_gamma_trans_streams_krnl<<<BLKd[f],THD,0,cf_stream[f]>>>(dpar,
					ddat, s, f, ndop[f], type, FLOAT);
		} checkErrorAfterKernelLaunch("cf_finish_fit_store_streams kernels and "
				"cf_gamma_trans_krnl (calc_doppler_cuda_streams");
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
__host__ void calc_lghtcrv_cuda_streams(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, struct vertices_t **verts, int s, int nframes, int
		nviews, unsigned char type, int lc_n, int nf, cudaStream_t *cf_stream)
{
	int ncalc, c=0, i, *posn, *bistatic, n, nThreads, exclude_seen, f,
			bistatic_all, *outbndarr;
	int nfplus = nframes+1; /* This is to accomodate the +1 start in lghtcrv */
	double *pxlpkm, posbnd_logfactor, orbit_offset[3] = {0.0, 0.0, 0.0};
	float3 orbit_off3;
	orbit_off3.x = orbit_off3.y = orbit_off3.z = 0.0;
	double3 *so;
	int4 *xylim;
	dim3 BLKpx[nfplus],THD,THD9;
	THD.x = maxThreadsPerBlock; THD9.x = 9;
	struct pos_t **pos;
	struct crvrend_t **rend;
	ncalc = nframes;
	n = lc_n;
	float *overflow;
	int hbistatic[nfplus], houtbndarr[nfplus],
		nThreadspx[nfplus], nThreadspx1[nfplus], hposn[nfplus];
	int4 hxylim[nfplus];
	int2 span[nfplus];

	gpuErrchk(cudaMalloc((void**)&pos, 	  sizeof(struct pos_t*) * nfplus));
	gpuErrchk(cudaMalloc((void**)&rend, 	  sizeof(struct crvrend_t*) * nfplus));
	gpuErrchk(cudaMalloc((void**)&posn, 	  sizeof(int) * nfplus));
	gpuErrchk(cudaMalloc((void**)&bistatic, sizeof(int) * nfplus));
	gpuErrchk(cudaMalloc((void**)&pxlpkm,   sizeof(double) * nfplus));
	gpuErrchk(cudaMalloc((void**)&overflow, sizeof(float) * nfplus));
	gpuErrchk(cudaMalloc((void**)&xylim, 	  sizeof(int4) * nfplus));
	gpuErrchk(cudaMalloc((void**)&outbndarr,sizeof(int) * nfplus));
	gpuErrchk(cudaMalloc((void**)&so, 	  sizeof(double3) * ((nframes*3)+1)));

	/* Set shortcuts and copy pos->n back for all frames */
	for (f=1; f<=nframes; f++){
		cfs_set_lghtcrv_shortcuts_streams_krnl<<<1,1,0,cf_stream[f-1]>>>(ddat,
				rend, pos, overflow, posn, s, f);
	checkErrorAfterKernelLaunch("cfs_set_lghtcrv_shortcuts_streams_krnl");
	}

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

	for (f=1; f<=ncalc; f++) {
		/* Launch 9-threaded kernel to set pos->ae,pos->oe,pos->bistatic.*/
		cfs_set_pos_ae_streams_krnl<<<1,THD9,0,cf_stream[f-1]>>>(pos, f,
				bistatic, type, 0);

		/* Clear the POS-view to initialize */
		posclr_streams_krnl<<<BLKpx[f],THD,0,cf_stream[f-1]>>>(pos, posn, f);
	} checkErrorAfterKernelLaunch("cfs_set_pos_ae_streams and posclr_streams_krnl");

	/* Call routine posvis to get  facet number, scattering & incidence
	 * angle, distance toward Earth at center of each POS pixel; set posbnd
	 * parameter = 1 if any model portion extends beyond POS frame limits*/
	posvis_cuda_streams2(dpar, dmod, ddat, pos, verts, orbit_off3, hposn,
				outbndarr, s, nframes, 0, nf, 0, c, type, cf_stream);

	/* Copy the posbnd flag returns for all frames to a host copy */
	gpuErrchk(cudaMemcpy(&houtbndarr, outbndarr, sizeof(int)*(nfplus),
			cudaMemcpyDeviceToHost));

	/* Now check the outbndarr for the posbnd flag for each frame */
	for (f=1; f<=nframes; f++) {
		if (houtbndarr[f]) {
			/* Call single-threaded kernel to set dpar->posbnd and dpar->posbnd_logfactor */
			cfs_set_posbnd_streams_krnl<<<1,1,0,cf_stream[f-1]>>>(dpar, pos, f, type);
			houtbndarr[f]=0;
		}
	}
	gpuErrchk(cudaMemcpy(&hbistatic, bistatic, sizeof(int)*(nfplus), cudaMemcpyDeviceToHost));

	/* Now view model from source (sun) and get facet number and distance
	 * toward source of each pixel in this projected view; use this
	 * information to determine which POS pixels are shadowed       */
	/* Because posvis_cuda_streams processes all frames at the same time, if
	 * any of the frames are bistatic, all of them get calculated again  */
	for (f=1; f<=nframes; f++)
		if (hbistatic[f])
			bistatic_all = 1;

	if (bistatic_all) {
//		posvis_cuda_streams(dpar, dmod, ddat, orbit_offset, s, nframes, 1, 0, c,
//				outbndarr, cf_stream);
		posvis_cuda_streams2(dpar, dmod, ddat, pos, verts, orbit_off3, hposn,
				outbndarr, s, nframes, 1, nf, 0, c, type, cf_stream);

		/* Copy the posbnd flag returns for all frames to a host copy */
		gpuErrchk(cudaMemcpy(&houtbndarr, outbndarr, sizeof(int)*nfplus,
				cudaMemcpyDeviceToHost));

		/* Now check the outbndarr for the posbnd flag for each frame */
		for (f=1; f<=nframes; f++) {
			if (houtbndarr[f]) {
				/* Call single-threaded kernel to set dpar->posbnd and dpar->posbnd_logfactor */
				cfs_set_posbnd_streams_krnl<<<1,1,0,cf_stream[f-1]>>>(dpar, pos, f, type);
			}
			/* Initialize this stream for the posmask kernel to follow */
			posmask_init_streams_krnl<<<1,1,0,cf_stream[f-1]>>>(pos,so,	pxlpkm,	f);

			/* Now call posmask kernel for this stream, then loop
			 * to next stream and repeat 	 */
			posmask_streams_krnl<<<BLKpx[f],THD,0,cf_stream[f-1]>>>(
					dpar, pos, so, pxlpkm, posn, nThreadspx[f],	span[f].x, f);

		} checkErrorAfterKernelLaunch("posmask_streams_krnl");
	}
	/* Go through all visible and unshadowed POS pixels with low enough
	 * scattering and incidence angles, and mark facets which project onto
	 * their centers as having been "seen" at least once   */
	/* Get xlim and ylim and exclude_seen flag */
	for (f=1; f<=nframes; f++)
		cfs_get_exclude_seen_streams_krnl<<<1,1,0,cf_stream[f-1]>>>(dpar,pos,xylim,f);
	checkErrorAfterKernelLaunch("cfs_get_exclude_seen_streams_krnl (calc_lghtcrv_cuda_streams");

	/* Now copy the flag and all frame pos's xlim and ylim values back from GPU */
	gpuErrchk(cudaMemcpyFromSymbol(&exclude_seen, cfs_exclude_seen, sizeof(int),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(&hxylim, xylim, sizeof(int4)*nfplus, cudaMemcpyDeviceToHost));

	/* Calculate launch parameters for all frames */
	for (f=1; f<=nframes; f++) {
		span[f].x = hxylim[f].x - hxylim[f].w + 1;
		span[f].y = hxylim[f].z - hxylim[f].y + 1;
		nThreadspx1[f] = span[f].x * span[f].y;
		BLKpx[f].x = floor ((THD.x -1 + nThreadspx1[f]) / THD.x);
	}

	for (f=1; f<=nframes; f++) {
		/* Go through all POS pixels which are visible with low enough
		 * scattering angle and mark the facets which project onto
		 * their centers as having been "seen" at least once           */
		if (s != exclude_seen)
			cf_mark_pixels_seen_streams_krnl<<<BLKpx[f],THD,0,cf_stream[f-1]>>>(
					dpar, dmod, pos, xylim, nThreadspx1[f], span[f].x, f);
	} checkErrorAfterKernelLaunch("cf_mark_pixels_krnl (calc_lghtcrv_cuda)");

	/* Compute model brightness for this lightcurve point then copy to device  */
	apply_photo_cuda_streams(dmod, ddat, pos, xylim, span, BLKpx, nThreadspx1,
			0, s, nframes, nThreadspx, cf_stream);

//	dbg_print_lghtcrv_pos_arrays(ddat, s, 22, nThreads, pos_n);

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
	float *u;
	gpuErrchk(cudaMalloc((void**)&u, sizeof(float) * ncalc));
	int threads = 128;
	BLKpx[0].x = floor((threads-1+nThreads)/threads);
	THD.x = threads;
	//cf_spline_lghtcrv_krnl<<<BLK,THD>>>(2.0e30, 2.0e30, u);
	cfs_spline_lghtcrv_serial2_krnl<<<1,1>>>(u);
	checkErrorAfterKernelLaunch("cf_spline_lghtcrv_krnl");

	/* Start debug */
	/* Pull out lghtcrv->x, lghtcrv->y, lghtcrv->y2 (all of length ncalc) */
	//dbg_print_lghtcrv_xyy2(ddat, s, ncalc, "xyy2_arrays_CUDA.csv");


	/* Launch n-threaded kernel to do the following:
	 * 	- set each fit[i] = 0
	 * 	- loop through all views and splint
	 * 	- add interp to fit[i]
	 * 	- divide final fit[i] over nviews
	 */
	BLKpx[0].x = floor((threads-1+n)/threads);
	THD.x = threads;
	cfs_splint_lghtcrv_krnl<<<BLKpx[0],THD>>>(dpar);
	checkErrorAfterKernelLaunch("cf_splint_lghtcrv_krnl");

	cudaFree(pos);
	cudaFree(rend);
	cudaFree(posn);
	cudaFree(bistatic);
	cudaFree(pxlpkm);
	cudaFree(overflow);
	cudaFree(xylim);
	cudaFree(outbndarr);
	cudaFree(so);
	cudaFree(u);
}
__host__ void calc_lghtcrv_cuda_streams2(
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
		cudaStream_t *cf_stream)
{
	int ncalc, c=0, i, *posn, *bistatic, n, nThreads, exclude_seen, f,
			bistatic_all, *outbndarr;
	int nfplus = nframes+1; /* This is to accomodate the +1 start in lghtcrv */
	float *pxlpkm, posbnd_logfactor, orbit_offset[3] = {0.0, 0.0, 0.0};
	float3 orbit_off3;
	orbit_off3.x = orbit_off3.y = orbit_off3.z = 0.0;
	double3 *so;
	int4 *xylim;
	dim3 BLKpx[nfplus],THD,THD9;
	THD.x = maxThreadsPerBlock; THD9.x = 9;
	struct pos_t **pos;
	struct crvrend_t **rend;
	ncalc = nframes;
	n = lc_n;
	float *overflow;
	int hbistatic[nfplus], houtbndarr[nfplus],
		nThreadspx[nfplus], nThreadspx1[nfplus], hposn[nfplus];
	int4 hxylim[nfplus];
	int2 span[nfplus];

	gpuErrchk(cudaMalloc((void**)&pos, 	  sizeof(struct pos_t*) * nfplus));
	gpuErrchk(cudaMalloc((void**)&rend, 	  sizeof(struct crvrend_t*) * nfplus));
	gpuErrchk(cudaMalloc((void**)&posn, 	  sizeof(int) * nfplus));
	gpuErrchk(cudaMalloc((void**)&bistatic, sizeof(int) * nfplus));
	gpuErrchk(cudaMalloc((void**)&pxlpkm,   sizeof(float) * nfplus));
	gpuErrchk(cudaMalloc((void**)&overflow, sizeof(float) * nfplus));
	gpuErrchk(cudaMalloc((void**)&xylim, 	  sizeof(int4) * nfplus));
	gpuErrchk(cudaMalloc((void**)&outbndarr,sizeof(int) * nfplus));
	gpuErrchk(cudaMalloc((void**)&so, 	  sizeof(double3) * ((nframes*3)+1)));

	/* Set shortcuts and copy pos->n back for all frames */
	for (f=1; f<=nframes; f++){
		cfs_set_lghtcrv_shortcuts_streams_krnl<<<1,1,0,cf_stream[f-1]>>>(ddat,
				rend, pos, overflow, posn, s, f);
	checkErrorAfterKernelLaunch("cfs_set_lghtcrv_shortcuts_streams_krnl");
	}

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

	for (f=1; f<=ncalc; f++) {
		/* Launch 9-threaded kernel to set pos->ae,pos->oe,pos->bistatic.*/
		cfs_set_pos_ae_streams_krnl<<<1,THD9,0,cf_stream[f-1]>>>(pos, f,
				bistatic, type, 0);

		/* Clear the POS-view to initialize */
		posclr_streams_krnl<<<BLKpx[f],THD,0,cf_stream[f-1]>>>(pos, posn, f);
	} checkErrorAfterKernelLaunch("cfs_set_pos_ae_streams and posclr_streams_krnl");

	/* Call routine posvis to get  facet number, scattering & incidence
	 * angle, distance toward Earth at center of each POS pixel; set posbnd
	 * parameter = 1 if any model portion extends beyond POS frame limits*/

	//dbg_print_facet_normals(dmod, nf, "STR2_facet_normals.csv");


	posvis_cuda_streams2(dpar, dmod, ddat, pos, verts, orbit_off3, hposn,
				outbndarr, s, nframes, 0, nf, 0, c, type, cf_stream);
//	int npixels = (2*hposn[1]+1)*(2*hposn[1]+1);
//	dbg_print_pos_arrays_full(pos, 1, npixels, hposn[1]);
//	dbg_print_pos_arrays_full(pos, 2, npixels, hposn[2]);
//	dbg_print_pos_arrays_full(pos, 3, npixels, hposn[3]);

	/* Copy the posbnd flag returns for all frames to a host copy */
	gpuErrchk(cudaMemcpy(&houtbndarr, outbndarr, sizeof(int)*(nfplus),
			cudaMemcpyDeviceToHost));

	/* Now check the outbndarr for the posbnd flag for each frame */
	for (f=1; f<=nframes; f++) {
		if (houtbndarr[f]) {
			/* Call single-threaded kernel to set dpar->posbnd and dpar->posbnd_logfactor */
			cfs_set_posbnd_streams_krnl<<<1,1,0,cf_stream[f-1]>>>(dpar, pos, f, type);
			houtbndarr[f]=0;
		}
	}
	gpuErrchk(cudaMemcpy(&hbistatic, bistatic, sizeof(int)*(nfplus), cudaMemcpyDeviceToHost));

	/* Now view model from source (sun) and get facet number and distance
	 * toward source of each pixel in this projected view; use this
	 * information to determine which POS pixels are shadowed       */
	/* Because posvis_cuda_streams processes all frames at the same time, if
	 * any of the frames are bistatic, all of them get calculated again  */
	for (f=1; f<=nframes; f++)
		if (hbistatic[f])
			bistatic_all = 1;

	if (bistatic_all) {
		posvis_cuda_streams2(dpar, dmod, ddat, pos, verts, orbit_off3, hposn,
				outbndarr, s, nframes, 1, nf, 0, c, type, cf_stream);

		/* Copy the posbnd flag returns for all frames to a host copy */
		gpuErrchk(cudaMemcpy(&houtbndarr, outbndarr, sizeof(int)*nfplus,
				cudaMemcpyDeviceToHost));

		/* Now check the outbndarr for the posbnd flag for each frame */
		for (f=1; f<=nframes; f++) {
			if (houtbndarr[f]) {
				/* Call single-threaded kernel to set dpar->posbnd and dpar->posbnd_logfactor */
				cfs_set_posbnd_streams_krnl<<<1,1,0,cf_stream[f-1]>>>(dpar, pos, f, type);
			}
			/* Initialize this stream for the posmask kernel to follow */
			posmask_init_streams2_krnl<<<1,1,0,cf_stream[f-1]>>>(pos, so, pxlpkm, f);

			/* Now call posmask kernel for this stream, then loop
			 * to next stream and repeat 	 */
			posmask_streams2_krnl<<<BLKpx[f],THD,0,cf_stream[f-1]>>>(
					dpar, pos, so, pxlpkm, posn, nThreadspx[f],	span[f].x, f);

		} checkErrorAfterKernelLaunch("posmask_streams_krnl");
	}
//	for (f=21; f<=22; f++)
//		dbg_print_lghtcrv_pos_arrays(ddat, 0, f, nThreadspx[f], hposn[f]);
	/* Go through all visible and unshadowed POS pixels with low enough
	 * scattering and incidence angles, and mark facets which project onto
	 * their centers as having been "seen" at least once   */
	/* Get xlim and ylim and exclude_seen flag */
	for (f=1; f<=nframes; f++)
		cfs_get_exclude_seen_streams_krnl<<<1,1,0,cf_stream[f-1]>>>(dpar,pos,xylim,f);
	checkErrorAfterKernelLaunch("cfs_get_exclude_seen_streams_krnl (calc_lghtcrv_cuda_streams");

	/* Now copy the flag and all frame pos's xlim and ylim values back from GPU */
	gpuErrchk(cudaMemcpyFromSymbol(&exclude_seen, cfs_exclude_seen, sizeof(int),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(&hxylim, xylim, sizeof(int4)*nfplus, cudaMemcpyDeviceToHost));

	/* Calculate launch parameters for all frames */
	for (f=1; f<=nframes; f++) {
		span[f].x = hxylim[f].x - hxylim[f].w + 1;
		span[f].y = hxylim[f].z - hxylim[f].y + 1;
		nThreadspx1[f] = span[f].x * span[f].y;
		BLKpx[f].x = floor ((THD.x -1 + nThreadspx1[f]) / THD.x);
	}

	for (f=1; f<=nframes; f++) {
		/* Go through all POS pixels which are visible with low enough
		 * scattering angle and mark the facets which project onto
		 * their centers as having been "seen" at least once           */
		if (s != exclude_seen)
			cf_mark_pixels_seen_streams_krnl<<<BLKpx[f],THD,0,cf_stream[f-1]>>>(
					dpar, dmod, pos, xylim, nThreadspx1[f], span[f].x, f);
	} checkErrorAfterKernelLaunch("cf_mark_pixels_krnl (calc_lghtcrv_cuda)");

	/* Compute model brightness for this lightcurve point then copy to device  */
	apply_photo_cuda_streams(dmod, ddat, pos, xylim, span, BLKpx, nThreadspx1,
			0, s, nframes, nThreadspx, cf_stream);


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
	double *u;
	THD.x = maxThreadsPerBlock;
	gpuErrchk(cudaMalloc((void**)&u, sizeof(double)*(nframes+1)));
	BLKpx[0].x = floor((THD.x - 1 + nThreads)/THD.x);
//	lghtcrv_spline_streams_krnl<<<BLKpx[0],THD>>>(ddat, s, 2.0e30,
//			2.0e30, u, nframes);
	lghtcrv_spline_streams_test_krnl<<<1,1>>>(ddat, s, 2.0e30, 2.0e30, u, nframes);
	checkErrorAfterKernelLaunch("lghtcrv_spline_streams_krnl");

	/* Change launch parameters from ncalc threads to n threads */
	BLKpx[0].x = floor((THD.x - 1 + n) / THD.x);
	//lghtcrv_splint_streams3_krnl<<<BLKpx[0],THD>>>(ddat, s, n, nframes);
	lghtcrv_splint_streams3_test_krnl<<<1,1>>>(ddat, s, n, nframes);
	checkErrorAfterKernelLaunch("lghtcrv_splint_streams_krnl");


//	dbg_print_lc_fit(ddat, s, "STR2_lghtcrv_fit.csv", n);

//	f = 1;
//	int npixels = (2*hposn[f]+1)*(2*hposn[f]+1);
//		dbg_print_pos_arrays_full(pos, f, npixels, hposn[f]);
//		dbg_print_pos_bd(pos, f, npixels, hposn[f]);
//	dbg_print_lghtcrv_xyy2(ddat, s, nframes, "streams2_x_y_y2_arrays.csv");
	//	dbg_print_lghtcrv_pos_arrays(ddat, s, 1, nThreadspx[1], hposn[1]);
//	int threads = 128;
//	BLKpx[0].x = floor((threads-1+nThreads)/threads);
//	THD.x = threads;
//	//cf_spline_lghtcrv_krnl<<<BLK,THD>>>(2.0e30, 2.0e30, u);
//	cfs_spline_lghtcrv_serial_krnl<<<1,1>>>(u);
//	checkErrorAfterKernelLaunch("cf_spline_lghtcrv_krnl");

	/* Start debug */
	/* Pull out lghtcrv->x, lghtcrv->y, lghtcrv->y2 (all of length ncalc) */
	//dbg_print_lghtcrv_xyy2(ddat, s, ncalc, "xyy2_arrays_CUDA.csv");


	/* Launch n-threaded kernel to do the following:
	 * 	- set each fit[i] = 0
	 * 	- loop through all views and splint
	 * 	- add interp to fit[i]
	 * 	- divide final fit[i] over nviews
	 */
//	BLKpx[0].x = floor((threads-1+n)/threads);
//	THD.x = threads;
//	cfs_splint_lghtcrv_krnl<<<BLKpx[0],THD>>>(dpar);
	//	checkErrorAfterKernelLaunch("cf_splint_lghtcrv_krnl");
//
//	int npixels = (2*hposn[1]+1)*(2*hposn[1]+1);
//
//	dbg_print_pos_arrays_full(pos, 1, npixels, hposn[1]);
//	dbg_print_pos_arrays_full(pos, 2, npixels, hposn[2]);
//	dbg_print_pos_arrays_full(pos, 3, npixels, hposn[3]);

	cudaFree(pos);
	cudaFree(rend);
	cudaFree(posn);
	cudaFree(bistatic);
	cudaFree(pxlpkm);
	cudaFree(overflow);
	cudaFree(xylim);
	cudaFree(outbndarr);
	cudaFree(so);
	cudaFree(u);
}
__host__ void calc_lghtcrv_cuda_streams2f(
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
		cudaStream_t *cf_stream)
{
	int nfplus = nframes+1; /* This is to accomodate the +1 start in lghtcrv */
	int ncalc, c=0, i, *posn, *bistatic, n, nThreads, exclude_seen, f,
			bistatic_all, *outbndarr, hbistatic[nfplus], houtbndarr[nfplus],
		nThreadspx[nfplus], nThreadspx1[nfplus], hposn[nfplus] ;
	float *pxlpkm, posbnd_logfactor;
	float3 orbit_off3, *so;
	int2 span[nfplus];
	int4 *xylim, hxylim[nfplus];
	dim3 BLKpx[nfplus],THD,THD9;
	THD.x = maxThreadsPerBlock; THD9.x = 9;
	struct pos_t **pos;
	struct crvrend_t **rend;
	float *overflow;

	ncalc = nframes;
	n = lc_n;
	orbit_off3.x = orbit_off3.y = orbit_off3.z = 0.0;

	gpuErrchk(cudaMalloc((void**)&pos, 	  sizeof(struct pos_t*) * nfplus));
	gpuErrchk(cudaMalloc((void**)&rend, 	  sizeof(struct crvrend_t*) * nfplus));
	gpuErrchk(cudaMalloc((void**)&posn, 	  sizeof(int) * nfplus));
	gpuErrchk(cudaMalloc((void**)&bistatic, sizeof(int) * nfplus));
	gpuErrchk(cudaMalloc((void**)&pxlpkm,   sizeof(float) * nfplus));
	gpuErrchk(cudaMalloc((void**)&overflow, sizeof(float) * nfplus));
	gpuErrchk(cudaMalloc((void**)&xylim, 	  sizeof(int4) * nfplus));
	gpuErrchk(cudaMalloc((void**)&outbndarr,sizeof(int) * nfplus));
	gpuErrchk(cudaMalloc((void**)&so, 	  sizeof(double3) * ((nframes*3)+1)));

	/* Set shortcuts and copy pos->n back for all frames */
	for (f=1; f<=nframes; f++){
		cfs_set_lghtcrv_shortcuts_streams_krnl<<<1,1,0,cf_stream[f-1]>>>(ddat,
				rend, pos, overflow, posn, s, f);
	checkErrorAfterKernelLaunch("cfs_set_lghtcrv_shortcuts_streams_krnl");
	}
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

	for (f=1; f<=ncalc; f++) {
		/* Launch 9-threaded kernel to set pos->ae,pos->oe,pos->bistatic.*/
		cfs_set_pos_ae_streams_krnl<<<1,THD9,0,cf_stream[f-1]>>>(pos, f,
				bistatic, type, 0);

		/* Clear the POS-view to initialize */
		posclr_streams_krnl<<<BLKpx[f],THD,0,cf_stream[f-1]>>>(pos, posn, f);
	} checkErrorAfterKernelLaunch("cfs_set_pos_ae_streams and posclr_streams_krnl");

	/* Call routine posvis to get  facet number, scattering & incidence
	 * angle, distance toward Earth at center of each POS pixel; set posbnd
	 * parameter = 1 if any model portion extends beyond POS frame limits*/

	posvis_cuda_streams2(dpar, dmod, ddat, pos, verts, orbit_off3, hposn,
				outbndarr, s, nframes, 0, nf, 0, c, type, cf_stream);

	/* Copy the posbnd flag returns for all frames to a host copy */
	gpuErrchk(cudaMemcpy(&houtbndarr, outbndarr, sizeof(int)*(nfplus),
			cudaMemcpyDeviceToHost));

	/* Now check the outbndarr for the posbnd flag for each frame */
	for (f=1; f<=nframes; f++) {
		if (houtbndarr[f]) {
			/* Call single-threaded kernel to set dpar->posbnd and dpar->posbnd_logfactor */
			cfs_set_posbnd_streams_krnl<<<1,1,0,cf_stream[f-1]>>>(dpar, pos, f, type);
			houtbndarr[f]=0;
		}
	}
	gpuErrchk(cudaMemcpy(&hbistatic, bistatic, sizeof(int)*(nfplus), cudaMemcpyDeviceToHost));

	/* Now view model from source (sun) and get facet number and distance
	 * toward source of each pixel in this projected view; use this
	 * information to determine which POS pixels are shadowed       */
	/* Because posvis_cuda_streams processes all frames at the same time, if
	 * any of the frames are bistatic, all of them get calculated again  */
	for (f=1; f<=nframes; f++)
		if (hbistatic[f])	bistatic_all = 1;

	if (bistatic_all) {
//		posvis_cuda_streams(dpar, dmod, ddat, orbit_offset, s, nframes, 1, 0, c,
//				outbndarr, cf_stream);
		posvis_cuda_streams2(dpar, dmod, ddat, pos, verts, orbit_off3, hposn,
				outbndarr, s, nframes, 1, nf, 0, c, type, cf_stream);

		/* Copy the posbnd flag returns for all frames to a host copy */
		gpuErrchk(cudaMemcpy(&houtbndarr, outbndarr, sizeof(int)*nfplus,
				cudaMemcpyDeviceToHost));

		/* Now check the outbndarr for the posbnd flag for each frame */
		for (f=1; f<=nframes; f++) {
			if (houtbndarr[f]) {
				/* Call single-threaded kernel to set dpar->posbnd and dpar->posbnd_logfactor */
				cfs_set_posbnd_streams_krnl<<<1,1,0,cf_stream[f-1]>>>(dpar, pos, f, type);
			}
			/* Initialize this stream for the posmask kernel to follow */
			posmask_init_streams_f_krnl<<<1,1,0,cf_stream[f-1]>>>(pos, so, pxlpkm, f);

			/* Now call posmask kernel for this stream, then loop
			 * to next stream and repeat 	 */
			posmask_streams_f_krnl<<<BLKpx[f],THD,0,cf_stream[f-1]>>>(
					dpar, pos, so, pxlpkm, nThreadspx[f], span[f].x, f);
		} checkErrorAfterKernelLaunch("posmask_streams_krnl");
	}
	/* Go through all visible and unshadowed POS pixels with low enough
	 * scattering and incidence angles, and mark facets which project onto
	 * their centers as having been "seen" at least once   */
	/* Get xlim and ylim and exclude_seen flag */
	for (f=1; f<=nframes; f++)
		cfs_get_exclude_seen_streams_krnl<<<1,1,0,cf_stream[f-1]>>>(dpar,pos,xylim,f);
	checkErrorAfterKernelLaunch("cfs_get_exclude_seen_streams_krnl (calc_lghtcrv_cuda_streams");

	/* Now copy the flag and all frame pos's xlim and ylim values back from GPU */
	gpuErrchk(cudaMemcpyFromSymbol(&exclude_seen, cfs_exclude_seen, sizeof(int),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(&hxylim, xylim, sizeof(int4)*nfplus, cudaMemcpyDeviceToHost));

	/* Calculate launch parameters for all frames */
	for (f=1; f<=nframes; f++) {
		span[f].x = hxylim[f].x - hxylim[f].w + 1;
		span[f].y = hxylim[f].z - hxylim[f].y + 1;
		nThreadspx1[f] = span[f].x * span[f].y;
		BLKpx[f].x = floor ((THD.x -1 + nThreadspx1[f]) / THD.x);
	}

	for (f=1; f<=nframes; f++) {
		/* Go through all POS pixels which are visible with low enough
		 * scattering angle and mark the facets which project onto
		 * their centers as having been "seen" at least once           */
		if (s != exclude_seen)
			cf_mark_pixels_seen_streams_krnl<<<BLKpx[f],THD,0,cf_stream[f-1]>>>(
					dpar, dmod, pos, xylim, nThreadspx1[f], span[f].x, f);
	} checkErrorAfterKernelLaunch("cf_mark_pixels_krnl (calc_lghtcrv_cuda)");

	/* Compute model brightness for this lightcurve point then copy to device  */
	apply_photo_cuda_streams_f(dmod, ddat, pos, xylim, span, BLKpx, nThreadspx1,
			0, s, nframes, cf_stream);

//	dbg_print_lghtcrv_xyy2(ddat, s, nframes, "streams2_x_y_y2_arrays.csv");
//	dbg_print_lghtcrv_pos_arrays(ddat, s, 1, nThreadspx[1], hposn[1]);

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
	float *u;
	gpuErrchk(cudaMalloc((void**)&u, sizeof(float) * (nframes+1)));

	THD.x = maxThreadsPerBlock;
	BLKpx[0].x = floor((THD.x-1+nThreads)/THD.x);
	//cf_spline_lghtcrv_krnl<<<BLK,THD>>>(2.0e30, 2.0e30, u);
	lghtcrv_spline_streams_f_krnl<<<BLKpx[0],THD>>>(ddat, s, 2.0e30, 2.0e30, u, nframes);
	checkErrorAfterKernelLaunch("lghtcrv_spline_streams_f_krnl in calc_fits_cuda_streams");

	/* Change launch parameters from ncalc threads to n threads */
	BLKpx[0].x = floor((THD.x - 1 + lc_n) / THD.x);
	lghtcrv_splint_streams3f_krnl<<<BLKpx[0],THD>>>(ddat, s, n);
	checkErrorAfterKernelLaunch("lghtcrv_splint_streams_krnl");


	/* Start debug */
	/* Pull out lghtcrv->x, lghtcrv->y, lghtcrv->y2 (all of length ncalc) */
	//dbg_print_lghtcrv_xyy2(ddat, s, ncalc, "xyy2_arrays_CUDA.csv");

	cudaFree(pos);
	cudaFree(rend);
	cudaFree(posn);
	cudaFree(bistatic);
	cudaFree(pxlpkm);
	cudaFree(overflow);
	cudaFree(xylim);
	cudaFree(outbndarr);
	cudaFree(so);
	cudaFree(u);
}
