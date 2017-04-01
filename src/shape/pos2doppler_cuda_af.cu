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
#include "head.h"
}

/* Declare __device__ vars and structs, which have file scope */
__device__ int afdop_nsinc2_sq, afdop_any_overflow, afdop_in_bounds,
			   afdop_badradar;

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
__global__ void pos2doppler_init_af_krnl(struct dat_t *ddat, int set, int v,
		struct dopfrm_t **frame, struct pos_t **pos, int nframes, int *ndop,
		int *idop0, float3 *w, float2 *doplim, int4 *xylim) {
	/* nframes-threaded kernel */
	int frm = threadIdx.x;
	if (frm < nframes) {
		/* Initialize variables  */
		idop0[frm] = 0;
		afdop_any_overflow = 0;
		frame[frm] = &ddat->set[set].desc.doppler.frame[frm];
		pos[frm] = &frame[frm]->pos;
		ndop[frm] = frame[frm]->ndop;
		frame[frm]->idoplim[0] = ndop[frm] + 999999;
		frame[frm]->idoplim[1] = -999999;
		frame[frm]->doplim[0] =  HUGENUMBER;
		frame[frm]->doplim[1] = -HUGENUMBER;
		afdop_badradar = 0;
		frame[frm]->badradar_logfactor = 0.0;

		/*  Get w, the apparent spin vector in observer coordinates  */
		dev_cotrans4(w, frame[frm]->view[v].oe, frame[frm]->view[v].spin, 1, frm);

		/* Copy frame->doplim over to the float device variable */
		doplim[frm].x = frame[frm]->doplim[0];
		doplim[frm].y = frame[frm]->doplim[1];

		/* Now get pos->xlim[0], pos->xlim[1], pos->ylim[0], pos->ylim[1] */
		xylim[frm].w = pos[frm]->xlim[0];
		xylim[frm].x = pos[frm]->xlim[1];
		xylim[frm].y = pos[frm]->ylim[0];
		xylim[frm].z = pos[frm]->ylim[1];
	}
}
__global__ void pos2doppler_radar_parameters_af_krnl(struct par_t *dpar,
		struct dat_t *ddat, struct dopfrm_t **frame, struct pos_t **pos,
		double orbit_dopoff, int set, int nframes, int v, float2 *axay,
		float2 *xyincr, float3 *w, float4 *dop, float *dopshift) {
	/* nframes-threaded kernel */
	int frm = threadIdx.x;
	double dopfact;
	if (frm < nframes) {
		/*  Compute the Doppler bin increment per plane-of-sky pixel westward (ax)
      and northward (ay); these values are scaled by the "dopscale" parameter
      for this dataset.  Then compute km2Hz, the Doppler increment (Hz) per
      km perpendicular to the projected spin axis in the plane of the sky.     */

		dopfact = ddat->set[set].desc.doppler.dopscale.val * KM2HZFACT * pos[0]->km_per_pixel
			* ddat->set[set].desc.doppler.Ftx / ddat->set[set].desc.doppler.dop_per_bin;
		axay[frm].x = -w[frm].y*dopfact;
		axay[frm].y =  w[frm].x*dopfact;
		frame[frm]->view[v].km2Hz = sqrt(axay[frm].x*axay[frm].x +
				axay[frm].y*axay[frm].y) * ddat->set[set].desc.doppler.dop_per_bin
			/ pos[frm]->km_per_pixel;

	/* Compute absolute value of the difference between maximum (or minimum)
	 * Doppler on any given POS pixel's edge and the Doppler at its center                               */
	/* dop.w - dopdiff_bl
	 * dop.x - dopdiff_max
	 * dop.y - dopDC_vig
	 * dop.z - dop_extra		 */

	if (w[frm].x != 0.0 || w[frm].y != 0.0)
		dop[frm].z = frame[frm]->view[v].km2Hz * 0.5 * pos[frm]->km_per_pixel
			* sqrt(w[frm].x*w[frm].x + w[frm].y*w[frm].y) /
			MAX( fabs(w[frm].x), fabs(w[frm].y));
	else
		dop[frm].z = 0.0;

	/*  We may be evaluating the sinc^2 Doppler response function at
      more than one point per POS pixel.  xincr and yincr are the
      Doppler bin increments between adjacent evaluation points in the
      x and y directions.  dopdiff_bl is the Doppler bin difference
      between the bottom-leftmost (southeasternmost) evaluation point
      and the pixel center.  dopdiff_max is the maximum positive
      Doppler bin difference between any evaluation point and the
      pixel center.                                                     */

	afdop_nsinc2_sq = dpar->nsinc2 * dpar->nsinc2;
	xyincr[frm].x = axay[frm].x / dpar->nsinc2;
	xyincr[frm].y = axay[frm].y / dpar->nsinc2;
	dop[frm].w = -(dpar->nsinc2 - 1)*(xyincr[frm].x + xyincr[frm].y)/2;
	dop[frm].x = (dpar->nsinc2 - 1)*(fabs(xyincr[frm].x) + fabs(xyincr[frm].x))/2;
	if (2*dop[frm].x + dpar->sinc2width + 1 > MAXBINS) {
		afdop_badradar = 1;
		frame[frm]->badradar_logfactor += log((2*dop[frm].x + dpar->sinc2width + 1) / MAXBINS);
		if (dpar->warn_badradar) {
			printf("\nWARNING in pos2doppler.c for set %2d frame %2d:\n", set, frm);
			printf("        sinc^2 function evaluated at %d Doppler bins, must be <= %d\n",
					(int) ceil(2*dop[frm].x + dpar->sinc2width + 1), MAXBINS);
		}
	}

	/*  Get the COM Doppler bin, corrected for ephemeris drift
      and adjusted for orbital motion                         */

	dopshift[frm] = frame[frm]->dopcom_vig + frame[frm]->view[v].dopoff + orbit_dopoff;
	}
}
__global__ void pos2doppler_get_global_frmsz_krnl(int *global_lim, int4 *xylim,
		int nframes) {
	/* nframes-threaded kernel */
	int f = threadIdx.x;
	if (f < nframes) {
		/* Initialize global_lim 	 */
		for (int i=0; i<4; i++)
			global_lim[i] = 0;

		/* Now calculate minimum for all frames */
		atomicMin(&global_lim[0], xylim[f].w);
		atomicMax(&global_lim[1], xylim[f].x);
		atomicMin(&global_lim[2], xylim[f].y);
		atomicMax(&global_lim[3], xylim[f].z);
	}
}
__global__ void pos2doppler_pixel_af_krnl(
		struct par_t *dpar,
		struct mod_t *dmod,
		struct dat_t *ddat,
		struct pos_t **pos,
		struct dopfrm_t **frame,
		int xspan, int set,	int nframes, int frame_size, int total_size, int body,
		double orbit_xoff, double orbit_yoff,
		float2 *axay, float2 *doplim, float2 *xyincr, float4 *dop,
		int *ndop, int *idop0, int *global_lim,
		float *dopshift) {
	/* Multi-threaded kernel */
	int total_offset = blockIdx.x * blockDim.x + threadIdx.x;
	int offset = total_offset % frame_size;
	int frm = total_offset / frame_size;
	int x = offset % xspan + global_lim[0]; // pos[frm]->xlim[0];
	int y = offset / xspan + global_lim[2]; // pos[frm]->ylim[0];
	int n ;

	int idop, idop_min, idop_max, idop1, idop2, i, j, c, f, k, zaddr;
	double tmp, amp, arg_left, sinc2arg, sinc2_mean, arg_bl, fit_contribution,
	sumweights, dop_contribution[MAXBINS], dopPOS;

	if ((offset < frame_size) && (frm < nframes)) {
		n = pos[frm]->n;

		zaddr = (y+n)*(2*n+1) + (x+n);

		/* Loop through all POS pixels within the rectangular plane-of-sky
		 * region spanned by the model; for each such pixel which isn't blank
		 * sky, compute the cross-section contributions to pixels in the model
		 * Doppler spectrum. Note that functions posclr and posvis flag
		 * blank-sky pixels by assigning "cose" = cos(scattering angle) = 0.
		 * Only compute contributions from POS pixels that project onto the
		 * right body, in case this is the "orbit" action (for which this
		 * routine is called twice, once for each of the 2 orbiting bodies).*/

		if (pos[frm]->cose_s[zaddr] > 0.0 && pos[frm]->body[x][y] == body) {

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
			dopPOS = axay[frm].x*(x - orbit_xoff) + axay[frm].y*(y - orbit_yoff) + dopshift[frm];
			idop_min = (int) floor(dopPOS - dop[frm].x + 1 - dpar->sinc2width/2.0);
			idop_max = (int) floor(dopPOS + dop[frm].x + dpar->sinc2width/2.0);

			/* Update the rectangular delay-Doppler region with nonzero power
			 * according to the model              */
			atomicMin(&frame[frm]->idoplim[0], idop_min);
			atomicMax(&frame[frm]->idoplim[1], idop_max);

			/* Update model's fp Doppler limits, as determined prior to
			 * convolution with the Doppler response function. At this point
			 * in the code, doplim is a pair of floating-point bin numbers
			 * which applies to POS pixel centers; when the loop over POS
			 * pixels is finished we will convert these limits to Hz and will
			 * widen the limits to account for nonzero POS pixel width. */
			/* Note that p2d_doplim[2] is a single-precision (float) copy of
			 * the original p2d_frame->doplim[2] (double-precision). This is
			 * necessary to get atomic operations to work.		 */
			atomicMinf(&doplim[frm].x, dopPOS);
			atomicMaxf(&doplim[frm].y, dopPOS);

			/* Check if all Doppler bins which will receive power from this POS
			 * pixel fall within the data frame; if not, initialize the
			 * "overflow" spectrum if necessary.  */
			if ( (idop_min >= 1) && (idop_max <= ndop[frm]) )
				afdop_in_bounds = 1;
			else {
				afdop_in_bounds = 0;
				if (!afdop_any_overflow) {
					afdop_any_overflow = 1;
					for (j=0; j<MAXOVERFLOW; j++)
						frame[frm]->fit_overflow[j] = 0.0;  // To-Do:  This might need attention.

					/* Center the COM in the overflow spectrum: bin [idop] in
					 * the fit frame corresponds to bin [idop+idop0] in the
					 * fit_overflow frame.  */
					idop0[frm] = MAXOVERFLOW/2 - (int) floor(dopshift[frm] + 0.5);
				}
			}

			/* Compute the sinc^2 factors for Doppler mismatching: Take the
			 * mean of nsinc2^2 points interior to the POS pixel. Do the two
			 * most common cases (nsinc2 = 1 or 2) without loops to gain speed.
			 * Note the SINC2 macro multiplies its argument by pi. Then add the
			 * cross-section contributions to the model spectrum.  */

			for (idop=idop_min; idop<=idop_max; idop++) {
				switch (dpar->nsinc2) {
				case 1:
					sinc2_mean = SINC2(dopPOS - idop);
					break;
				case 2:
					arg_bl = dopPOS + dop[frm].w - idop;   /* bl = bottom left */
					sinc2_mean = (SINC2(arg_bl) +
							SINC2(arg_bl+xyincr[frm].x) +
							SINC2(arg_bl+xyincr[frm].y) +
							SINC2(arg_bl+xyincr[frm].x+xyincr[frm].y)) / 4;
					break;
				default:
					arg_left = dopPOS + dop[frm].w - idop;
					sinc2_mean = 0.0;
					for (i=0; i<dpar->nsinc2; i++) {
						sinc2arg = arg_left;
						for (j=0; j<dpar->nsinc2; j++) {
							sinc2_mean += SINC2(sinc2arg);
							sinc2arg += xyincr[frm].x;
						}
						arg_left += xyincr[frm].y;
					}
					sinc2_mean /= afdop_nsinc2_sq;
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
					pos[frm]->cose_s[zaddr], pos[frm]->comp[x][y], pos[frm]->f[x][y])
		    		 * pos[frm]->km_per_pixel * pos[frm]->km_per_pixel  / sumweights;

			/* Only add POS pixel's power contributions to model Doppler spect-
			 * rum if NONE of those contributions fall outside spectrum limits*/
			if (afdop_in_bounds) {

				/*  Add the cross-section contributions to the model frame  */
				for (idop=idop_min; idop<=idop_max; idop++) {
					k = MIN( idop - idop_min, MAXBINS);
					fit_contribution = amp * dop_contribution[k];
					atomicAdd(&ddat->set[set].desc.doppler.frame[frm].fit_s[idop],
							fit_contribution);

					if (dpar->action == MAP) {
						if (dpar->map_mode == MAPMODE_DELDOP) {
							if (frame[frm]->map_fit[idop] > 0.0) {
								frame[frm]->map_pos[x][y] += fit_contribution;
								c = pos[frm]->comp[x][y];
								f = pos[frm]->f[x][y];
								frame[frm]->map_facet_power[c][f] += fit_contribution;
								if (dpar->map_verbose)
									printf("# POS (%3d, %3d) comp %d facet %4d contributes %e to Dop (%3d)\n",
											x+n, y+n, c, f, fit_contribution, idop-1);
							}
						} else if (dpar->map_mode == MAPMODE_POS) {
							if (frame[frm]->map_pos[x][y] > 0.0) {
								frame[frm]->map_fit[idop] += fit_contribution;
								c = pos[frm]->comp[x][y];
								f = pos[frm]->f[x][y];
								frame[frm]->map_facet_power[c][f] += fit_contribution;
								if (dpar->map_verbose)
									printf("# POS (%3d, %3d) comp %d facet %4d contributes %e to Dop (%3d)\n",
											x+n, y+n, c, f, fit_contribution, idop-1);
							}
						} else {
							if (frame[frm]->map_pos[x][y] > 0.0) {
								frame[frm]->map_fit[idop] += fit_contribution;
								if (dpar->map_verbose) {
									c = pos[frm]->comp[x][y];
									f = pos[frm]->f[x][y];
									printf("# POS (%3d, %3d) comp %d facet %4d contributes %e to Dop (%3d)\n",
											x+n, y+n, c, f, fit_contribution, idop-1);
								}
							}
						}
					}
				}
			} else {

				/*  Add the cross-section contributions to the "overflow" spectrum  */

				if (dpar->action == MAP && dpar->map_mode != MAPMODE_DELDOP)
					if (frame[frm]->map_pos[x][y] > 0.0)
						dpar->map_overflow = 1;

				idop1 = MAX( idop_min, -idop0[frm]);
				idop2 = MIN( idop_max, -idop0[frm] + MAXOVERFLOW - 1);
				for (idop=idop1; idop<=idop2; idop++) {
					k = MIN( idop - idop_min, MAXBINS);
					fit_contribution = amp * dop_contribution[k];
					frame[frm]->fit_overflow[idop+idop0[frm]] += fit_contribution;  // might need atomics
					if (dpar->action == MAP && dpar->map_mode == MAPMODE_DELDOP)
						if (idop >= dpar->map_doplim[0] && idop <= dpar->map_doplim[1]) {
							frame[frm]->map_pos[x][y] += fit_contribution;
							c = pos[frm]->comp[x][y];
							f = pos[frm]->f[x][y];
							frame[frm]->map_facet_power[c][f] += fit_contribution;
							if (dpar->map_verbose)
								printf("# POS (%3d, %3d) comp %d facet %4d contributes %e to Dop (%3d)\n",
										x+n, y+n, c, f, fit_contribution, idop-1);
						}
				}
			}
		}  /* if cos(scattering angle) > 0 and POS pixel projects onto the right body */
	}
}
__global__ void pos2doppler_finish_krnl(struct par_t *dpar, struct dat_t *ddat,
		struct dopfrm_t **frame, float4 *dop, float2 *doplim, float *dopshift,
		int *idop0, int *ndop, int set, int nframes) {
	/* Single-threaded kernel */
	int frm = threadIdx.x;
	int j, j1, j2;
	double lookfact, sdev_sq, variance, dopfactor;

	if (threadIdx.x <nframes) {
		/* Copy float device variable over to the frame->doplim   */
		frame[frm]->doplim[0] = doplim[frm].x;
		frame[frm]->doplim[1] = doplim[frm].y;

		/* Convert model's Doppler limits from float bin numbers to Hz and
		 * widen the limits to account for nonzero POS pixel width    */

		frame[frm]->doplim[0] = (frame[frm]->doplim[0] - dopshift[frm])*
				ddat->set[set].desc.doppler.dop_per_bin - dop[frm].z;
		frame[frm]->doplim[1] = (frame[frm]->doplim[1] - dopshift[frm])*
				ddat->set[set].desc.doppler.dop_per_bin + dop[frm].z;

		/* Calculate overflow contributions to chi squared:
		 *   o2 = obs^2 contribution, m2 = model^2 contribution.
		 * Also compute summed cross section and mean Doppler bin for overflow
		 * region, for use with the "delcorinit" action   */

		frame[frm]->overflow_o2 = 0.0;
		frame[frm]->overflow_m2 = 0.0;
		frame[frm]->overflow_xsec = 0.0;
		frame[frm]->overflow_dopmean = 0.0;
		sdev_sq = frame[frm]->sdev*frame[frm]->sdev;
		variance = sdev_sq;
		lookfact = (frame[frm]->nlooks > 0.0) ? 1.0/frame[frm]->nlooks : 0.0;
		if (afdop_any_overflow) {
			j1 = MAX(frame[frm]->idoplim[0] + idop0[frm], 0);
			j2 = MIN(frame[frm]->idoplim[1] + idop0[frm], MAXOVERFLOW - 1);
			for (j=j1; j<=j2; j++) {
				if (frame[frm]->fit_overflow[j] != 0.0) {
					if (dpar->speckle)
						variance = sdev_sq + lookfact*frame[frm]->fit_overflow[j]*
								frame[frm]->fit_overflow[j];
					frame[frm]->overflow_o2 += 1.0;
					frame[frm]->overflow_m2 += frame[frm]->fit_overflow[j]*
							frame[frm]->fit_overflow[j]/variance;
					frame[frm]->overflow_xsec += frame[frm]->fit_overflow[j];
					frame[frm]->overflow_dopmean += (j - idop0[frm])*frame[frm]->fit_overflow[j];
				}
			}
			if (frame[frm]->overflow_xsec != 0.0)
				frame[frm]->overflow_dopmean /= frame[frm]->overflow_xsec;

			/*  Print a warning if the model extends even beyond the overflow spectrum  */

			if ( ((frame[frm]->idoplim[0] + idop0[frm]) < 0) ||
					((frame[frm]->idoplim[1] + idop0[frm]) >= MAXOVERFLOW) ) {
				afdop_badradar = 1;
				dopfactor = (MAX(frame[frm]->idoplim[1] + idop0[frm], MAXOVERFLOW)
						- MIN(frame[frm]->idoplim[0] + idop0[frm], 0)         )
		                								  / (1.0*MAXOVERFLOW);
				frame[frm]->badradar_logfactor += log(dopfactor);
				if (dpar->warn_badradar) {
					printf("\nWARNING in pos2doppler.c for set %2d frame %2d:\n", set, frm);
					printf("        model Doppler spectrum extends too far beyond the data spectrum\n");
					printf("             data:  bins %2d to %2d\n", 1, ndop[frm]);
					printf("            model:  bins %2d to %2d\n",
							frame[frm]->idoplim[0], frame[frm]->idoplim[1]);
				}
			}
		}
	}
}
__host__ int pos2doppler_cuda_af( struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, double orbit_xoff, double orbit_yoff, double
		orbit_dopoff, int body, int set, int nframes, int v)
{
	int badradar, xspan, yspan, nThreads, frmsz, *global_lim, *idop0, *ndop;
	dim3 BLK, THD;
	struct dopfrm_t **frame;
	struct pos_t **pos;
	float *dopshift;
	float2 *doplim, *axay, *xyincr;
	float3 *w;
	float4 *dop;
	int4 *xylim;

	cudaCalloc((void**)&frame, 	  	 sizeof(struct dopfrm_t*),nframes);
	cudaCalloc((void**)&pos, 	  	 sizeof(struct pos_t*),   nframes);
	cudaCalloc((void**)&dopshift, 	 sizeof(float), 		  nframes);
	cudaCalloc((void**)&axay, 	  	 sizeof(float2), 		  nframes);
	cudaCalloc((void**)&xyincr,   	 sizeof(float2), 		  nframes);
	cudaCalloc((void**)&doplim,   	 sizeof(float2), 		  nframes);
	cudaCalloc((void**)&w, 		  	 sizeof(float3), 		  nframes);
	cudaCalloc((void**)&dop, 	  	 sizeof(float4), 		  nframes);
	cudaCalloc((void**)&xylim, 	  	 sizeof(int4), 			  nframes);
	cudaCalloc((void**)&global_lim,	 sizeof(int), 			        4);
	cudaCalloc((void**)&idop0,	   	 sizeof(int), 			  nframes);
	cudaCalloc((void**)&ndop,	   	 sizeof(int), 			  nframes);

	/* Launch nframes-threaded initialization kernel */
	THD.x = nframes;
	pos2doppler_init_af_krnl<<<1,THD>>>(ddat, set, v, frame, pos, nframes,
			ndop, idop0, w,	doplim,	xylim);
	checkErrorAfterKernelLaunch("pos2doppler_init_af_krnl");

	pos2doppler_radar_parameters_af_krnl<<<1,THD>>>(dpar, ddat, frame, pos,
			orbit_dopoff, set,nframes, v, axay, xyincr, w, dop,	dopshift);
	checkErrorAfterKernelLaunch("pos2doppler_radar_parameters_af_krnl");

	/* Figure out the largest pos->xlim/ylim window for the entire set */
	pos2doppler_get_global_frmsz_krnl<<<1,1>>>(global_lim, xylim, nframes);
	checkErrorAfterKernelLaunch("pos2doppler_get_global_frmsz_krnl");
	deviceSyncAfterKernelLaunch("pos2doppler_get_global_frmsz_krnl");

	/* Configure the pixel kernel  */
	xspan = global_lim[1] - global_lim[0] + 1; //xlim1 - xlim0 + 1;
	yspan = global_lim[3] - global_lim[2] + 1; //ylim1 - ylim0 + 1;
	frmsz = xspan * yspan;
	nThreads = frmsz * nframes;
	BLK.x = floor((maxThreadsPerBlock - 1 + nThreads) / maxThreadsPerBlock);
	THD.x = maxThreadsPerBlock; // Thread block dimensions

	pos2doppler_pixel_af_krnl<<<BLK,THD>>>(dpar, dmod, ddat, pos, frame, xspan,
			set, nframes, frmsz, nThreads, body, orbit_xoff, orbit_yoff,
			axay, doplim, xyincr, dop, ndop, idop0, global_lim, dopshift);
	checkErrorAfterKernelLaunch("pos2doppler_pixel_af_krnl");

	THD.x = nframes;
	/* Launch the single-thread kernel to finish up Doppler calculations */
	pos2doppler_finish_krnl<<<1,THD.x>>>(dpar, ddat, frame, dop, doplim,
			dopshift, idop0, ndop, set, nframes);
	checkErrorAfterKernelLaunch("pos2doppler_finish_krnl, line ");
	gpuErrchk(cudaMemcpyFromSymbol(&badradar, afdop_badradar, sizeof(badradar),
			0, cudaMemcpyDeviceToHost));

//	int debug = 0;
//	if (debug)
//		dbg_print_fit(ddat, set, 3);

	cudaFree(frame);
	cudaFree(pos);
	cudaFree(dopshift);
	cudaFree(axay);
	cudaFree(xyincr);
	cudaFree(doplim);
	cudaFree(w);
	cudaFree(dop);
	cudaFree(xylim);
	cudaFree(global_lim);
	cudaFree(idop0);

	return badradar;
}
