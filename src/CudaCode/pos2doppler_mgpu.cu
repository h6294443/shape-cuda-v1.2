/*****************************************************************************************
                                                                            pos2doppler.c

Takes a pos_t structure and a doppler_t structure and "fills in" a dopfrm_t structure
indexed by frm.  In other words, pos2doppler works from a model plane-of-sky image (with
an observer z-coordinate and a scattering angle assigned to each pixel) to produce a model
Doppler spectrum corresponding to data frame frm.

In the case of an orbiting binary system (the "orbit" action), pos2doppler only computes
power contributions from the orbiting body denoted by the "body" parameter: the routine is
called twice, once for each body.

pos2doppler takes contributions only from the rectangular plane-of-sky region defined by
pos.xlim and pos.ylim -- the smallest rectangle which completely "contains" the model in
the plane of the sky.  No power is contributed by parts of the model which extend beyond
the POS window; this is why such models are heavily penalized (the objective function is
doubled -- see function f in file bestfit.c).

idoplim is updated for frame frm to show the model Doppler region that contains nonzero
power.

Modified 2016 November 14 by ME:
	Implemented an all-GPU version of pos2doppler.

Modified 2015 June 3 by CM:
    Implement smearing for the "fit" and "write" actions by adding "v" (view) parameter
        and applying each run of pos2doppler to a single view rather than to an entire
        (smeared) frame

Modified 2014 February 10 by CM:
    Add "ilaw" argument to the radlaw routine to implement multiple radar scattering laws

Modified 2012 March 23 by CM:
    Implement Doppler scaling

Modified 2010 September 1 by CM:
    Add braces to an if-then-else statement to avoid compilation warning

Modified 2010 August 26 by CM:
    For the "map" action, change the "map_forward" parameter to "map_mode"
        and implement map_mode = 'facets'
    For the "map" action, implement the "map_verbose" parameter

Modified 2010 June 15 by CM:
    Pass the entire par_t parameter structure as an argument rather than
        just selected parameters
    Implement the map action

Modified 2009 November 15 by CM:
    Fix argument type in a printf statement

Modified 2009 April 3 by CM:
    Add the "warn_badradar" parameter (see below)
    If the sinc^2 Doppler response function extends to too many Doppler
        bins, set a flag and compute a factor by which the objective
        function should be increased (actually the logarithm of this
        factor).  If the "warn_badradar" parameter is turned on, print an
        explicit warning.
    If the model is too wide in Doppler space even for the overflow image,
        set a flag and compute a factor by which the objective function
        should be increased (actually the logarithm of this factor).  If
        the "warn_badradar" parameter is turned on, print an explicit
        warning.
    Make pos2doppler int rather than void in order to return the flag
        described above to the calling procedure

Modified 2007 August 4 by CM:
    Add orbit_xoff, orbit_yoff, and orbit_dopoff parameters, the x offset
        (POS image rows), y offset (POS image columns), and Doppler offset
        (spectral bins) of the center of mass due to orbital motion.
    Add body parameter to indicate (for the "orbit" action) which of the two
        orbiting bodies' power contributions should be computed
    Add c (component) argument to radlaw routine

Modified 2006 September 14 by CM:
    If the overflow region is too small, print a warning rather than
        halting the program

Modified 2006 June 21 by CM:
    Change dopres to dop_per_bin
    For POS renderings, change res to km_per_pixel

Modified 2006 June 18 by CM:
    Allow each Doppler frame in a dataset to have different dimensions
        after vignetting

Modified 2006 March 10 by CM:
    Pass the "speckle" parameter so that self-noise can be included when
        computing the chi squared contribution of the overflow region
    Compute overflow_xsec and overflow_dopmean so that these quantities
        can be used by the "delcorinit" action

Modified 2005 July 25 by CM:
    Fix bug in overall cross-section scale factor: return to Scott's scheme
        of normalizing the cross-section contributions from a given POS 
        pixel so that they sum to the cross section actually present on the 
        sky in that pixel

Modified 2005 July 20 by CM:
    Fix bug in computing floating-point Doppler limits in Hz
    Add "facet" argument to radlaw routine

Modified 2005 July 5 by CM:
    Eliminate "dir" argument (since we always add power to the model image
       and never subtract it)
    Add "set" (set number) argument in order to improve error messages

Modified 2005 June 27 by CM:
    Rename INFINITY constant to HUGENUMBER to avoid conflicts

Modified 2005 June 25 by CM:
    Rename old "doplim" to "idoplim"; this is the Doppler limits in
        (integer) bin numbers
    Add new "doplim" which is the floating-point Doppler limits in Hz,
        obtained PRIOR to convolution with the Doppler response function

Modified 2005 January 25 by CM:
    Take care of uninitialized variable

Modified 2003 May 11 by CM:
    Compute contributions to chi squared by model power which lies
    outside the limits of the data frame.

Modified 2003 May 5 by CM:
    For each POS pixel, compute the entire pixel's contribution to a
    given Doppler bin in the model spectrum so long as even one point
    at which we evaluate the sinc^2 response function is less than
    sinc2width/2.0 bins away from that Doppler bin.  In other words,
    err on the side of computing too many small contributions to each
    bin in the model spectrum, so as not to omit significant contributions
    just because a POS pixel's *center* isn't close enough in Doppler.

Modified 2003 April 29 by CM:
    Evaluate the sinc^2 Doppler response function at nsinc2 points
    per POS pixel dimension, not just at the pixel center.
    The sinc^2 function varies rapidly -- one null per Doppler bin
    away from the central peak -- so if the pixel width is more than
    about half the Doppler resolution, we want to take the mean of
    several points within the pixel.

Modified 2003 April 26 by CM:
    Zero out the sinc^2 Doppler response function beyond the
    nearest sinc2width bins rather than beyond the nearest 2 bins

Modified 2003 April 17 by CM:
    Now correctly scales the model Doppler spectrum to account for
    Doppler mismatching
 *****************************************************************************************/
extern "C" {
#include "../shape/head.h"
}

/* Declare __device__ vars and structs, which have file scope */
/* In multi-gpu mode, each GPU gets its own copy of a device variable. These
 * are independent from each other and must be read/written to separately
 * by switching GPUs. */
__device__ int mgpu_nsinc2_sq, mgpu_any_overflow, mgpu_in_bounds, mgpu_badradar_global;

/*	Note that both pos2deldop_cuda.cu and posvis_cuda.cu have the atomicMaxf
 * 	and atomicMinf device functions defined separately.  This is done due to
 * 	the way static device functions I handled I guess.  I tried putting them
 * 	into a separate file and a declaration in the shape-cuda.h header file,
 * 	but to no avail.  So here they are, duplicated in both files.			*/
__device__ static float atomicMaxf(float* address, float val)
{
	int* address_as_i = (int*) address;
	int old = *address_as_i, assumed;
	do {
		assumed = old;
		old = ::atomicCAS(address_as_i, assumed,
				__float_as_int(::fmaxf(val, __int_as_float(assumed))));
	} while (assumed != old);
	return __int_as_float(old);
}

__device__ static float atomicMinf(float* address, float val)
{
	int* address_as_i = (int*) address;
	int old = *address_as_i, assumed;
	do {
		assumed = old;
		old = ::atomicCAS(address_as_i, assumed,
				__float_as_int(::fminf(val, __int_as_float(assumed))));
	} while (assumed != old);
	return __int_as_float(old);
}

__global__ void pos2doppler_init_mgpu_krnl(
		struct dat_t *ddat,
		struct dopfrm_t **frame,
		struct pos_t **pos,
		int4 *xylim,
		float3 *w,
		float2 *doplim,
		int *ndop,
		int *idop0,
		int set,
		int v,
		int size,
		int oddflg,
		int *badradararr) {

	/* nfrm_half0 or nfrm_half1-threaded kernel; goes through all half-frames
	 * and assigns values.  Intended for multi-GPU operation. */

	int hf = blockIdx.x * blockDim.x + threadIdx.x;
	int f = 2*hf + oddflg;
	if (hf < size) {
		/* Initialize variables  */
		idop0[hf] = 0;
		mgpu_any_overflow = 0;
		frame[hf] = &ddat->set[set].desc.doppler.frame[f];
		frame[hf]->idoplim[0] = ndop[hf] + 999999;
		frame[hf]->idoplim[1] = -999999;
		frame[hf]->doplim[0] =  HUGENUMBER;
		frame[hf]->doplim[1] = -HUGENUMBER;
		badradararr[hf] = 0;
		frame[hf]->badradar_logfactor = 0.0;

		/*  Get w, the apparent spin vector in observer coordinates  */
		dev_cotrans4(w, frame[hf]->view[v].oe, frame[hf]->view[v].spin, 1, hf);

		/* Copy frame->doplim over to the float device variable */
		doplim[hf].x = frame[hf]->doplim[0];
		doplim[hf].y = frame[hf]->doplim[1];

		xylim[hf].w = pos[hf]->xlim[0];
		xylim[hf].x = pos[hf]->xlim[1];
		xylim[hf].y = pos[hf]->ylim[0];
		xylim[hf].z = pos[hf]->ylim[1];
	}
}

__global__ void pos2doppler_radar_parameters_mgpu_krnl(
		struct par_t *dpar,
		struct dat_t *ddat,
		struct dopfrm_t **frame,
		struct pos_t **pos,
		float4 *dop,
		float3 *w,
		float2 *axay,
		float2 *xyincr,
		float *dopshift,
		double orbit_dopoff,
		int set,
		int size,
		int v,
		int oddflg,
		int *badradararr) {

	/* nfrm_half0 or nfrm_half1-threaded kernel for multi-gpu operation. GPU0
	 * does the even frames while GPU1 does the odd frames. */
	int hf = blockIdx.x * blockDim.x + threadIdx.x;
	int f = 2*hf + oddflg;
	float dopfact;
	__shared__ int dpb, nsinc2, sinc2width;
	float kmppxl;


	/* Compute the Doppler bin increment per plane-of-sky pixel westward
	 * (ax) and northward (ay); these values are scaled by the "dopscale"
	 * parameter for this dataset. Then compute km2Hz, the Doppler increment
	 * (Hz)/ km perpendicular to projected spin axis in the plane of sky.
	 * This happens behind a check that this is the first thread in the thread
	 * block. The calculations do not depend on thread index (hf) */
	if (threadIdx.x == 0) {
		dpb = ddat->set[set].desc.doppler.dop_per_bin;
		nsinc2 = dpar->nsinc2;
		sinc2width = dpar->sinc2width;
		dopfact = ddat->set[set].desc.doppler.dopscale.val * KM2HZFACT * kmppxl
				* ddat->set[set].desc.doppler.Ftx / dpb;
	}
	__syncthreads();

	/* Now we switch to all threads (index hf) */
	if (hf < size) {
		kmppxl = pos[hf]->km_per_pixel;
		axay[hf].x = -w[hf].y*dopfact;
		axay[hf].y =  w[hf].x*dopfact;
		frame[hf]->view[v].km2Hz = sqrt(axay[hf].x*axay[hf].x +
				axay[hf].y*axay[hf].y) * dpb	/ kmppxl;

		/* Compute absolute value of the difference between maximum (or minimum)
		 * Doppler on any given POS pixel's edge and the Doppler at its center                               */
		/* dop.w - dopdiff_bl
		 * dop.x - dopdiff_max
		 * dop.y - dopDC_vig
		 * dop.z - dop_extra		 */
		if (w[hf].x != 0.0 || w[hf].y != 0.0)
			dop[hf].z = frame[hf]->view[v].km2Hz * 0.5 * kmppxl *
			sqrt(w[hf].x*w[hf].x + w[hf].y*w[hf].y) /
			MAX( fabs(w[hf].x), fabs(w[hf].y));
		else
			dop[hf].z = 0.0;

		/* We may be evaluating the sinc^2 Doppler response function at more than
		 * one point per POS pixel.
		 * 	xincr & yincr: Doppler bin increments between adjacent evaluation
		 * 			points in the x and y directions.
		 * 	dopdiff_bl: Doppler bin difference between bottom-leftmost (SE-most)
		 * 			evaluation point and the pixel center.
		 * 	dopdiff_max: max. positive Doppler bin difference between any
		 * 			evaluation point and the     pixel center.             */
		mgpu_nsinc2_sq = nsinc2 * nsinc2;
		xyincr[hf].x = axay[hf].x / nsinc2;
		xyincr[hf].y = axay[hf].y / nsinc2;
		dop[hf].w = -(nsinc2 - 1)*(xyincr[hf].x + xyincr[hf].y)/2;
		dop[hf].x = (nsinc2 - 1)*(fabs(xyincr[hf].x) + fabs(xyincr[hf].x))/2;
		if (2*dop[hf].x + sinc2width + 1 > MAXBINS) {
			badradararr[hf] = 1;
			frame[hf]->badradar_logfactor += log((2*dop[hf].x + sinc2width + 1) / MAXBINS);
			if (dpar->warn_badradar) {
				printf("\nWARNING in pos2doppler_mgpu.c for set %2d frame %2d:\n", set, hf);
				printf("        sinc^2 function evaluated at %d Doppler bins, must be <= %d\n",
						(int) ceil(2*dop[hf].x + sinc2width + 1), MAXBINS);
			}
		}
		/* Get the COM Doppler bin, corrected for ephemeris drift and adjusted for
		 * orbital motion                         */
		dopshift[hf] = frame[hf]->dopcom_vig + frame[hf]->view[v].dopoff + orbit_dopoff;
	}
}

__global__ void pos2doppler_pixel_mgpu_krnl(
		struct par_t *dpar,
		struct mod_t *dmod,
		struct dat_t *ddat,
		struct pos_t **pos,
		struct dopfrm_t **frame,
		float4 *dop,
		float2 *axay,
		float2 *doplim,
		float2 *xyincr,
		int *ndop,
		int *idop0,
		float *dopshift,
		int xspan,
		int set,
		int frame_size,
		int body,
		double orbit_xoff,
		double orbit_yoff,
		int hf) {

	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int x = offset % xspan + pos[hf]->xlim[0];
	int y = offset / xspan + pos[hf]->ylim[0];
	int n, nsinc2, sinc2width ;

	int idop, idop_min, idop_max, idop1, idop2, i, j, c, fac, k, zaddr;
	double tmp, amp, arg_left, sinc2arg, sinc2_mean, arg_bl, fit_contribution,
	sumweights, dop_contribution[MAXBINS], dopPOS;
	n = pos[hf]->n;
	nsinc2 = dpar->nsinc2;
	sinc2width = dpar->sinc2width;

	if (offset < frame_size) {

		zaddr = (y+n)*(2*n+1) + (x+n);

		/* Loop through all POS pixels within the rectangular plane-of-sky
		 * region spanned by the model; for each such pixel which isn't blank
		 * sky, compute the cross-section contributions to pixels in the model
		 * Doppler spectrum. Note that functions posclr and posvis flag
		 * blank-sky pixels by assigning "cose" = cos(scattering angle) = 0.
		 * Only compute contributions from POS pixels that project onto the
		 * right body, in case this is the "orbit" action (for which this
		 * routine is called twice, once for each of the 2 orbiting bodies).*/

		if (pos[hf]->cose_s[zaddr] > 0.0 && pos[hf]->body[x][y] == body) {

			/* Get the fp Doppler bin of POS pixel center: dopPOS. Also get the
			 * min and max int Doppler bins to which this pixel contributes
			 * power: idop_min and idop_max. Each POS pixel contributes power
			 * to *all* Doppler bins, but here we're zeroing out the sinc^2
			 * response function beyond the nearest sinc2width bins.
			 * Actually, if nsinc2 > 1, we'll distribute power to *at least*
			 * sinc2width Doppler bins: For pixels which span multiple bins
			 * we'll err on the side of computing more contributions rather
			 * than fewer. */
			/* dop.w - dopdiff_bl
			 * dop.x - dopdiff_max
			 * dop.y - dopDC_vig
			 * dop.z - dop_extra		 */
			dopPOS = axay[hf].x*(x - orbit_xoff) + axay[hf].y* (y-orbit_yoff) +
					dopshift[hf];
			idop_min = (int) floor(dopPOS - dop[hf].x + 1 - sinc2width/2.0);
			idop_max = (int) floor(dopPOS + dop[hf].x + sinc2width/2.0);

			/* Update the rectangular delay-Doppler region with nonzero power
			 * according to the model              */
			atomicMin(&frame[hf]->idoplim[0], idop_min);
			atomicMax(&frame[hf]->idoplim[1], idop_max);

			/* Update model's fp Doppler limits, as determined prior to
			 * convolution with the Doppler response function. At this point
			 * in the code, doplim is a pair of floating-point bin numbers
			 * which applies to POS pixel centers; when the loop over POS
			 * pixels is finished we will convert these limits to Hz and will
			 * widen the limits to account for nonzero POS pixel width. */
			/* Note that p2d_doplim[2] is a single-precision (float) copy of
			 * the original p2d_frame->doplim[2] (double-precision). This is
			 * necessary to get atomic operations to work.		 */
			atomicMinf(&doplim[hf].x, dopPOS);
			atomicMaxf(&doplim[hf].y, dopPOS);

			/* Check if all Doppler bins which will receive power from this POS
			 * pixel fall within the data frame; if not, initialize the
			 * "overflow" spectrum if necessary.  */
			if ( (idop_min >= 1) && (idop_max <= ndop[hf]) )
				mgpu_in_bounds = 1;
			else {
				mgpu_in_bounds = 0;
				if (!mgpu_any_overflow) {
					mgpu_any_overflow = 1;
					for (j=0; j<MAXOVERFLOW; j++)
						frame[hf]->fit_overflow[j] = 0.0;  // To-Do:  This might need attention.

					/* Center the COM in the overflow spectrum: bin [idop] in
					 * the fit frame corresponds to bin [idop+idop0] in the
					 * fit_overflow frame.  */
					idop0[hf] = MAXOVERFLOW/2 - (int) floor(dopshift[hf] + 0.5);
				}
			}

			/* Compute the sinc^2 factors for Doppler mismatching: Take the
			 * mean of nsinc2^2 points interior to the POS pixel. Do the two
			 * most common cases (nsinc2 = 1 or 2) without loops to gain speed.
			 * Note the SINC2 macro multiplies its argument by pi. Then add the
			 * cross-section contributions to the model spectrum.  */

			for (idop=idop_min; idop<=idop_max; idop++) {
				switch (nsinc2) {
				case 1:
					sinc2_mean = SINC2(dopPOS - idop);
					break;
				case 2:
					arg_bl = dopPOS + dop[hf].w - idop;   /* bl = bottom left */
					sinc2_mean = (SINC2(arg_bl) +
							SINC2(arg_bl+xyincr[hf].x) +
							SINC2(arg_bl+xyincr[hf].y) +
							SINC2(arg_bl+xyincr[hf].x+xyincr[hf].y)) / 4;
					break;
				default:
					arg_left = dopPOS + dop[hf].w - idop;
					sinc2_mean = 0.0;
					for (i=0; i<nsinc2; i++) {
						sinc2arg = arg_left;
						for (j=0; j<nsinc2; j++) {
							sinc2_mean += SINC2(sinc2arg);
							sinc2arg += xyincr[hf].x;
						}
						arg_left += xyincr[hf].y;
					}
					sinc2_mean /= mgpu_nsinc2_sq;
					break;
				}
				k = MIN( idop - idop_min, MAXBINS);
				dop_contribution[k] = sinc2_mean;
			}

			/* Compute the sum of Doppler weighting factors  */
			sumweights = 0.0;
			for (idop=idop_min; idop<=idop_max; idop++) {
				k = MIN( idop - idop_min, MAXBINS);
				sumweights += dop_contribution[k];
			}

			/* The radar cross section within this plane-of-sky pixel is
			 * [differential radar scattering law]*[POS pixel area in km^2].
			 * The differential radar scattering law (function radlaw
			 * = d[cross section]/d[area] ) includes a sec(theta) factor to
			 * account for the fact that the POS pixel area is projected area
			 * rather than physical area on the target surface.	 */
			amp = dev_radlaw(&dmod->photo, ddat->set[set].desc.doppler.iradlaw,
					pos[hf]->cose_s[zaddr], pos[hf]->comp[x][y], pos[hf]->f[x][y])
		    		 * pos[hf]->km_per_pixel * pos[hf]->km_per_pixel / sumweights;

			/* Only add POS pixel's power contributions to model Doppler spect-
			 * rum if NONE of those contributions fall outside spectrum limits*/
			if (mgpu_in_bounds) {
				/* Add the cross-section contributions to the model frame  */
				for (idop=idop_min; idop<=idop_max; idop++) {
					k = MIN( idop - idop_min, MAXBINS);
					fit_contribution = amp * dop_contribution[k];

					atomicAdd(&frame[hf]->fit_s[idop], fit_contribution);

					if (dpar->action == MAP) {
						if (dpar->map_mode == MAPMODE_DELDOP) {
							if (frame[hf]->map_fit[idop] > 0.0) {
								frame[hf]->map_pos[x][y] += fit_contribution;
								c = pos[hf]->comp[x][y];
								fac = pos[hf]->f[x][y];
								frame[hf]->map_facet_power[c][fac] += fit_contribution;
								if (dpar->map_verbose)
									printf("# POS (%3d, %3d) comp %d facet %4d contributes %e to Dop (%3d)\n",
											x+n, y+n, c, fac, fit_contribution, idop-1);
							}
						} else if (dpar->map_mode == MAPMODE_POS) {
							if (frame[hf]->map_pos[x][y] > 0.0) {
								frame[hf]->map_fit[idop] += fit_contribution;
								c = pos[hf]->comp[x][y];
								fac = pos[hf]->f[x][y];
								frame[hf]->map_facet_power[c][fac] += fit_contribution;
								if (dpar->map_verbose)
									printf("# POS (%3d, %3d) comp %d facet %4d contributes %e to Dop (%3d)\n",
											x+n, y+n, c, fac, fit_contribution, idop-1);
							}
						} else {
							if (frame[hf]->map_pos[x][y] > 0.0) {
								frame[hf]->map_fit[idop] += fit_contribution;
								if (dpar->map_verbose) {
									c = pos[hf]->comp[x][y];
									fac = pos[hf]->f[x][y];
									printf("# POS (%3d, %3d) comp %d facet %4d contributes %e to Dop (%3d)\n",
											x+n, y+n, c, fac, fit_contribution, idop-1);
								}
							}
						}
					}
				}
			} else {

				/* Add the cross-section contributions to the "overflow" spectrum  */
				if (dpar->action == MAP && dpar->map_mode != MAPMODE_DELDOP)
					if (frame[hf]->map_pos[x][y] > 0.0)
						dpar->map_overflow = 1;

				idop1 = MAX( idop_min, -idop0[hf]);
				idop2 = MIN( idop_max, -idop0[hf] + MAXOVERFLOW - 1);
				for (idop=idop1; idop<=idop2; idop++) {
					k = MIN(idop - idop_min, MAXBINS);
					fit_contribution = amp * dop_contribution[k];
					frame[hf]->fit_overflow[idop+idop0[hf]] += fit_contribution;  // might need atomics
					if (dpar->action == MAP && dpar->map_mode == MAPMODE_DELDOP)
						if (idop >= dpar->map_doplim[0] && idop <= dpar->map_doplim[1]) {
							frame[hf]->map_pos[x][y] += fit_contribution;
							c = pos[hf]->comp[x][y];
							fac = pos[hf]->f[x][y];
							frame[hf]->map_facet_power[c][fac] += fit_contribution;
							if (dpar->map_verbose)
								printf("# POS (%3d, %3d) comp %d facet %4d contributes %e to Dop (%3d)\n",
										x+n, y+n, c, fac, fit_contribution, idop-1);
						}
				}
			}
		}  /* if cos(scattering angle) > 0 and POS pixel projects onto the right body */
	}
}

__global__ void pos2doppler_finish_mgpu_krnl(
		struct par_t *dpar,
		struct dat_t *ddat,
		struct dopfrm_t **frame,
		float4 *dop,
		float2 *doplim,
		float *dopshift,
		int *idop0,
		int *ndop,
		int set,
		int size,
		int *badradararr) {

	/* nfrm_half0/nfrm_half1-threaded kernel for multi-GPU operation.
	 * GPU0 handles even frames while GPU1 handles odd frames */

	int hf = blockDim.x*blockIdx.x + threadIdx.x;
	int j, j1, j2;
	double lookfact, sdev_sq, variance, dopfactor;

	if (hf < size) {
		/* Copy float device variable over to the frame->doplim   */
		frame[hf]->doplim[0] = doplim[hf].x;
		frame[hf]->doplim[1] = doplim[hf].y;

		/* Convert model's Doppler limits from float bin numbers to Hz and
		 * widen the limits to account for nonzero POS pixel width    */
		frame[hf]->doplim[0] = (frame[hf]->doplim[0] - dopshift[hf])*
				ddat->set[set].desc.doppler.dop_per_bin - dop[hf].z;
		frame[hf]->doplim[1] = (frame[hf]->doplim[1] - dopshift[hf])*
				ddat->set[set].desc.doppler.dop_per_bin + dop[hf].z;

		/* Calculate overflow contributions to chi squared:
		 *   o2 = obs^2 contribution, m2 = model^2 contribution.
		 * Also compute summed cross section and mean Doppler bin for overflow
		 * region, for use with the "delcorinit" action   */

		frame[hf]->overflow_o2 = 0.0;
		frame[hf]->overflow_m2 = 0.0;
		frame[hf]->overflow_xsec = 0.0;
		frame[hf]->overflow_dopmean = 0.0;
		sdev_sq = frame[hf]->sdev*frame[hf]->sdev;
		variance = sdev_sq;
		lookfact = (frame[hf]->nlooks > 0.0) ? 1.0/frame[hf]->nlooks : 0.0;

		/* Note that mgpu_any_overflow is instantiated on both GPUs in a dual-
		 * GPU setup. The copies are independent of one another.		 */

		if (mgpu_any_overflow) {
			j1 = MAX(frame[hf]->idoplim[0] + idop0[hf], 0);
			j2 = MIN(frame[hf]->idoplim[1] + idop0[hf], MAXOVERFLOW - 1);
			for (j=j1; j<=j2; j++) {
				if (frame[hf]->fit_overflow[j] != 0.0) {
					if (dpar->speckle)
						variance = sdev_sq + lookfact*frame[hf]->fit_overflow[j]*
								frame[hf]->fit_overflow[j];
					frame[hf]->overflow_o2 += 1.0;
					frame[hf]->overflow_m2 += frame[hf]->fit_overflow[j]*
							frame[hf]->fit_overflow[j]/variance;
					frame[hf]->overflow_xsec += frame[hf]->fit_overflow[j];
					frame[hf]->overflow_dopmean += (j - idop0[hf])*
							frame[hf]->fit_overflow[j];
				}
			}
			if (frame[hf]->overflow_xsec != 0.0)
				frame[hf]->overflow_dopmean /= frame[hf]->overflow_xsec;

			/* Print a warning if the model extends even beyond the overflow spectrum  */
			if (((frame[hf]->idoplim[0] + idop0[hf]) < 0) ||
					((frame[hf]->idoplim[1] + idop0[hf]) >= MAXOVERFLOW) ) {
				badradararr[hf] = 1;
				dopfactor = (MAX(frame[hf]->idoplim[1] + idop0[hf], MAXOVERFLOW)
					- MIN(frame[hf]->idoplim[0]+idop0[hf],0))/(1.0*MAXOVERFLOW);
				frame[hf]->badradar_logfactor += log(dopfactor);
				if (dpar->warn_badradar) {
					printf("\nWARNING in pos2doppler_mgpu.cu for set %2d frame %2d:\n", set, hf);
					printf("        model Doppler spectrum extends too far beyond the data spectrum\n");
					printf("             data:  bins %2d to %2d\n", 1, ndop[hf]);
					printf("            model:  bins %2d to %2d\n",
							frame[hf]->idoplim[0], frame[hf]->idoplim[1]);
				}
			}
		}

		if (badradararr[hf]) mgpu_badradar_global = 1;
	}
}


__host__ int pos2doppler_mgpu(
		struct par_t *dpar,
		struct mod_t *dmod,
		struct dat_t *ddat,
		struct pos_t **pos0,
		struct pos_t **pos1,
		struct dopfrm_t **frame0,
		struct dopfrm_t **frame1,
		int4 *xylim0,
		int4 *xylim1,
		double orbit_xoff,
		double orbit_yoff,
		double orbit_dopoff,
		int *ndop0,
		int *ndop1,
		int body,
		int set,
		int nframes,
		int v,
		int *badradararr0,
		int *badradararr1,
		cudaStream_t *gpu0_stream,
		cudaStream_t *gpu1_stream)
{
	int nfrm_half0 = nframes/2 + nframes%2;
	int nfrm_half1 = nframes/2;
	int badradar0, badradar1, xspan[nframes], yspan, nThreads[nframes], f,
		*idop0_0, *idop0_1;
	dim3 BLK_half0, BLK_half1, BLK[nframes], THD;
	float *dopshift0, *dopshift1;
	float2 *doplim0, *doplim1, *axay0, *axay1, *xyincr0, *xyincr1;
	float3 *w0, *w1;
	float4 *dop0, *dop1;
	int4 hxylim0[nfrm_half0], hxylim1[nfrm_half1];

	THD.x = maxThreadsPerBlock;
	BLK_half0.x = floor((THD.x - 1 + nfrm_half0));
	BLK_half1.x = floor((THD.x - 1 + nfrm_half1));

	/* Allocate GPU0 memory */
	gpuErrchk(cudaSetDevice(GPU0));
	gpuErrchk(cudaMalloc((void**)&dopshift0,sizeof(float)*nfrm_half0));
	gpuErrchk(cudaMalloc((void**)&axay0, 	sizeof(float2)*nfrm_half0));
	gpuErrchk(cudaMalloc((void**)&xyincr0,  sizeof(float2)*nfrm_half0));
	gpuErrchk(cudaMalloc((void**)&doplim0,  sizeof(float2)*nfrm_half0));
	gpuErrchk(cudaMalloc((void**)&w0, 	   	sizeof(float3)*nfrm_half0));
	gpuErrchk(cudaMalloc((void**)&dop0, 	sizeof(float4)*nfrm_half0));
	gpuErrchk(cudaMalloc((void**)&idop0_0,	sizeof(int)*nfrm_half0));

	/* Allocate GPU1 memory */
	gpuErrchk(cudaSetDevice(GPU1));
	gpuErrchk(cudaMalloc((void**)&dopshift1,sizeof(float)*nfrm_half1));
	gpuErrchk(cudaMalloc((void**)&axay1, 	sizeof(float2)*nfrm_half1));
	gpuErrchk(cudaMalloc((void**)&xyincr1,  sizeof(float2)*nfrm_half1));
	gpuErrchk(cudaMalloc((void**)&doplim1,  sizeof(float2)*nfrm_half1));
	gpuErrchk(cudaMalloc((void**)&w1, 	   	sizeof(float3)*nfrm_half1));
	gpuErrchk(cudaMalloc((void**)&dop1, 	sizeof(float4)*nfrm_half1));
	gpuErrchk(cudaMalloc((void**)&idop0_1,	sizeof(int)*nfrm_half1));

	/* Launch nfrm_half threaded kernel to initialize. GPU1/ Odd frames first.*/
	pos2doppler_init_mgpu_krnl<<<BLK_half1,THD,0,gpu1_stream[0]>>>(ddat, frame1,
					pos1, xylim1, w1, doplim1, ndop1, idop0_1, set, v,
					nfrm_half1, 1, badradararr1);
	pos2doppler_radar_parameters_mgpu_krnl<<<BLK_half1,THD,0,gpu1_stream[0]>>>
			(dpar, ddat, frame1, pos1, dop1, w1, axay1, xyincr1, dopshift1,
			orbit_dopoff, set, nfrm_half1, v, 1, badradararr1);
	gpuErrchk(cudaMemcpyAsync(hxylim1, xylim1, nfrm_half1*sizeof(int4),
			cudaMemcpyDeviceToHost, gpu1_stream[0]));

	/* Now even frames / GPU0 */
	gpuErrchk(cudaSetDevice(GPU0));
	pos2doppler_init_mgpu_krnl<<<BLK_half0,THD,0,gpu0_stream[0]>>>(ddat,frame0,
			pos0, xylim0, w0, doplim0, ndop0, idop0_0, set, v, nfrm_half0, 0,
			badradararr0);
	pos2doppler_radar_parameters_mgpu_krnl<<<BLK_half0,THD,0,gpu0_stream[0]>>>
			(dpar, ddat, frame0, pos0, dop0, w0, axay0, xyincr0, dopshift0,
					orbit_dopoff, set, nfrm_half0, v, 0, badradararr0);
	gpuErrchk(cudaMemcpyAsync(hxylim0, xylim0, nfrm_half0*sizeof(int4),
			cudaMemcpyDeviceToHost, gpu0_stream[0]));

	checkErrorAfterKernelLaunch("pos2doppler_mgpu first two kernels");

	/* Figure out the kernel launch parameters for every frame, both GPUs */
	int hf = 0;
	for (int f=0; f<nframes; f++) {
		/* Even frames first */
		xspan[f] = hxylim0[hf].x - hxylim0[hf].w + 1;
		yspan = hxylim0[hf].z - hxylim0[hf].y + 1;
		nThreads[f] = xspan[f]*yspan;
		BLK[f].x = floor((THD.x -1 + nThreads[f]) / THD.x);

		/* Odd frames */
		f++;	if (f==nframes)	break;
		xspan[f] = hxylim1[hf].x - hxylim1[hf].w + 1;
		yspan = hxylim1[hf].z - hxylim1[hf].y + 1;
		nThreads[f] = xspan[f]*yspan;
		BLK[f].x = floor((THD.x -1 + nThreads[f]) / THD.x);

		hf++;
	}

	/* GPU0 first */
	f = 0;
	for (int hf=0; hf<nfrm_half0; hf++) {
		pos2doppler_pixel_mgpu_krnl<<<BLK[f],THD,0,gpu0_stream[hf]>>>(dpar,
				dmod, ddat, pos0, frame0, dop0, axay0, doplim0, xyincr0, ndop0,
				idop0_0, dopshift0, xspan[f], set, nThreads[f], body,
				orbit_xoff, orbit_yoff, hf);
		f += 2;
	}
	/* Now GPU1 */
	gpuErrchk(cudaSetDevice(GPU1));

	f = 1;
	for (int hf=0; hf<nfrm_half1; hf++) {
		pos2doppler_pixel_mgpu_krnl<<<BLK[f],THD,0,gpu1_stream[hf]>>>(dpar,
				dmod, ddat, pos1, frame1, dop1, axay1, doplim1, xyincr1, ndop1,
				idop0_1, dopshift1, xspan[f], set, nThreads[f], body,
				orbit_xoff, orbit_yoff, hf);
		f += 2;
	}

	/* Finish up the overflow, GPU1 with the odd frames first */
	pos2doppler_finish_mgpu_krnl<<<BLK_half1,THD,0,gpu1_stream[0]>>>(dpar,ddat,
			frame1, dop1, doplim1, dopshift1, idop0_1, ndop1, set, nfrm_half1,
			badradararr1);
	gpuErrchk(cudaMemcpyFromSymbol(&badradar1, mgpu_badradar_global, sizeof(int),
			0, cudaMemcpyDeviceToHost));

	/* Now GPU0 for the even frames */
	gpuErrchk(cudaSetDevice(GPU0));
	pos2doppler_finish_mgpu_krnl<<<BLK_half0,THD,0,gpu0_stream[0]>>>(dpar,ddat,
			frame0, dop0, doplim0, dopshift0, idop0_0, ndop0, set, nfrm_half0,
			badradararr0);
	gpuErrchk(cudaMemcpyFromSymbol(&badradar0, mgpu_badradar_global, sizeof(int),
				0, cudaMemcpyDeviceToHost));
	/* Check for errors in the kernel launches & copy the badradar flag back */
	checkErrorAfterKernelLaunch("pos2doppler_finish_mgpu_krnl");

	if (badradar1)	badradar0 = 1;


	/* Free GPU1 memory */
	gpuErrchk(cudaSetDevice(GPU1));
	cudaFree(dopshift1);
	cudaFree(axay1);
	cudaFree(xyincr1);
	cudaFree(doplim1);
	cudaFree(w1);
	cudaFree(dop1);
	cudaFree(idop0_1);

	/* Free GPU0 memory */
	gpuErrchk(cudaSetDevice(GPU0));
	cudaFree(dopshift0);
	cudaFree(axay0);
	cudaFree(xyincr0);
	cudaFree(doplim0);
	cudaFree(w0);
	cudaFree(dop0);
	cudaFree(idop0_0);

	return badradar0;

}
