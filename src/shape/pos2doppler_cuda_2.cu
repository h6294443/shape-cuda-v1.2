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
__device__ struct dopfrm_t *p2d_frame;
__device__ struct pos_t *p2d_pos;
static __device__ double p2d_fit_overflow[MAXOVERFLOW];
__device__ int p2d_nsinc2_sq, p2d_any_overflow, p2d_in_bounds, p2d_ndop,
			   p2d_badradar, p2d_idop0, p2d_xlim0, p2d_xlim1, p2d_ylim0,
			   p2d_ylim1;
__device__ double p2d_ax, p2d_ay, p2d_xincr, p2d_yincr, p2d_dopshift,
				  p2d_dopdiff_bl, p2d_dopdiff_max, p2d_dop_extra,  p2d_w[3];
__device__ float p2d_doplim[2];

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
__global__ void pos2doppler_init_krnl(struct dat_t *ddat, int set,
		int frm, int v) {
	/* Single-threaded kernel */
	if (threadIdx.x ==0) {
		/* Initialize variables  */
		p2d_idop0 = 0;
		p2d_any_overflow = 0;
		p2d_frame = &ddat->set[set].desc.doppler.frame[frm];
		p2d_pos = &p2d_frame->pos;
		p2d_ndop = p2d_frame->ndop;
		p2d_frame->idoplim[0] = p2d_ndop + 999999;
		p2d_frame->idoplim[1] = -999999;
		p2d_frame->doplim[0] =  HUGENUMBER;
		p2d_frame->doplim[1] = -HUGENUMBER;
		p2d_badradar = 0;
		p2d_frame->badradar_logfactor = 0.0;

		/*  Get w, the apparent spin vector in observer coordinates  */
		dev_cotrans2(p2d_w, p2d_frame->view[v].oe, p2d_frame->view[v].spin, 1);

		/* Copy frame->doplim over to the float device variable */
		p2d_doplim[0] = p2d_frame->doplim[0];
		p2d_doplim[1] = p2d_frame->doplim[1];

		/* Now get pos->xlim[0], pos->xlim[1], pos->ylim[0], pos->ylim[1] */
		p2d_xlim0 = p2d_pos->xlim[0];
		p2d_xlim1 = p2d_pos->xlim[1];
		p2d_ylim0 = p2d_pos->ylim[0];
		p2d_ylim1 = p2d_pos->ylim[1];
	}
}
__global__ void pos2doppler_radar_parameters_krnl(struct par_t *dpar,
		struct dat_t *ddat, double orbit_dopoff, int set, int frm, int v) {
	/* Single-threaded kernel */
	double dopfact;
	if (threadIdx.x == 0) {
		/*  Compute the Doppler bin increment per plane-of-sky pixel westward (ax)
      and northward (ay); these values are scaled by the "dopscale" parameter
      for this dataset.  Then compute km2Hz, the Doppler increment (Hz) per
      km perpendicular to the projected spin axis in the plane of the sky.     */

		dopfact = ddat->set[set].desc.doppler.dopscale.val * KM2HZFACT * p2d_pos->km_per_pixel
			* ddat->set[set].desc.doppler.Ftx / ddat->set[set].desc.doppler.dop_per_bin;
		p2d_ax = -p2d_w[1]*dopfact;
		p2d_ay = p2d_w[0]*dopfact;
		p2d_frame->view[v].km2Hz = sqrt(p2d_ax*p2d_ax + p2d_ay*p2d_ay) * ddat->set[set].desc.doppler.dop_per_bin
			/ p2d_pos->km_per_pixel;

	/*  Compute the absolute value of the difference between the
      maximum (or minimum) Doppler on any given POS pixel's edge
      and the Doppler at its center                               */

	if (p2d_w[0] != 0.0 || p2d_w[1] != 0.0)
		p2d_dop_extra = p2d_frame->view[v].km2Hz * 0.5 * p2d_pos->km_per_pixel
		* sqrt(p2d_w[0]*p2d_w[0] + p2d_w[1]*p2d_w[1]) / MAX( fabs(p2d_w[0]), fabs(p2d_w[1]));
	else
		p2d_dop_extra = 0.0;

	/*  We may be evaluating the sinc^2 Doppler response function at
      more than one point per POS pixel.  xincr and yincr are the
      Doppler bin increments between adjacent evaluation points in the
      x and y directions.  dopdiff_bl is the Doppler bin difference
      between the bottom-leftmost (southeasternmost) evaluation point
      and the pixel center.  dopdiff_max is the maximum positive
      Doppler bin difference between any evaluation point and the
      pixel center.                                                     */

	p2d_nsinc2_sq = dpar->nsinc2 * dpar->nsinc2;
	p2d_xincr = p2d_ax / dpar->nsinc2;
	p2d_yincr = p2d_ay / dpar->nsinc2;
	p2d_dopdiff_bl = -(dpar->nsinc2 - 1)*(p2d_xincr + p2d_yincr)/2;
	p2d_dopdiff_max = (dpar->nsinc2 - 1)*(fabs(p2d_xincr) + fabs(p2d_yincr))/2;
	if (2*p2d_dopdiff_max + dpar->sinc2width + 1 > MAXBINS) {
		p2d_badradar = 1;
		p2d_frame->badradar_logfactor += log((2*p2d_dopdiff_max + dpar->sinc2width + 1) / MAXBINS);
		if (dpar->warn_badradar) {
			printf("\nWARNING in pos2doppler.c for set %2d frame %2d:\n", set, frm);
			printf("        sinc^2 function evaluated at %d Doppler bins, must be <= %d\n",
					(int) ceil(2*p2d_dopdiff_max + dpar->sinc2width + 1), MAXBINS);
		}
	}

	/*  Get the COM Doppler bin, corrected for ephemeris drift
      and adjusted for orbital motion                         */

	p2d_dopshift = p2d_frame->dopcom_vig + p2d_frame->view[v].dopoff + orbit_dopoff;
	}
}
__global__ void pos2doppler_pixel_krnl(struct par_t *dpar, struct mod_t
		*dmod, struct dat_t *ddat, int xspan, int set, int frame, int
		nThreads, int body, double orbit_xoff, double orbit_yoff) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int x = offset % xspan + p2d_pos->xlim[0];
	int y = offset / xspan + p2d_pos->ylim[0];

	int idop, idop_min, idop_max, idop1, idop2, i, j, c, f, k, zaddr;
	double tmp, amp, arg_left, sinc2arg, sinc2_mean, arg_bl, fit_contribution,
	sumweights, dop_contribution[MAXBINS], dopPOS;

	if (offset < nThreads) {

		zaddr = (y + p2d_pos->n)*(2*p2d_pos->n + 1) + (x + p2d_pos->n);

		/* Loop through all POS pixels within the rectangular plane-of-sky
		 * region spanned by the model; for each such pixel which isn't blank
		 * sky, compute the cross-section contributions to pixels in the model
		 * Doppler spectrum. Note that functions posclr and posvis flag
		 * blank-sky pixels by assigning "cose" = cos(scattering angle) = 0.
		 * Only compute contributions from POS pixels that project onto the
		 * right body, in case this is the "orbit" action (for which this
		 * routine is called twice, once for each of the 2 orbiting bodies).*/

		if (p2d_pos->cose_s[zaddr] > 0.0 && p2d_pos->body[x][y] == body) {

			/* Get the fp Doppler bin of POS pixel center: dopPOS. Also get the
			 * min and max int Doppler bins to which this pixel contributes
			 * power: idop_min and idop_max. Each POS pixel contributes power
			 * to *all* Doppler bins, but here we're zeroing out the sinc^2
			 * response function beyond the nearest sinc2width bins.
			 * Actually, if nsinc2 > 1, we'll distribute power to *at least*
			 * sinc2width Doppler bins: For pixels which span multiple bins
			 * we'll err on the side of computing more contributions rather
			 * than fewer. */

			dopPOS = p2d_ax*(x - orbit_xoff) + p2d_ay*(y - orbit_yoff) + p2d_dopshift;
			idop_min = (int) floor(dopPOS - p2d_dopdiff_max + 1 - dpar->sinc2width/2.0);
			idop_max = (int) floor(dopPOS + p2d_dopdiff_max + dpar->sinc2width/2.0);

			/* Update the rectangular delay-Doppler region with nonzero power
			 * according to the model              */
			atomicMin(&p2d_frame->idoplim[0], idop_min);
			atomicMax(&p2d_frame->idoplim[1], idop_max);

			/* Update model's fp Doppler limits, as determined prior to
			 * convolution with the Doppler response function. At this point
			 * in the code, doplim is a pair of floating-point bin numbers
			 * which applies to POS pixel centers; when the loop over POS
			 * pixels is finished we will convert these limits to Hz and will
			 * widen the limits to account for nonzero POS pixel width. */
			/* Note that p2d_doplim[2] is a single-precision (float) copy of
			 * the original p2d_frame->doplim[2] (double-precision). This is
			 * necessary to get atomic operations to work.		 */
			atomicMinf(&p2d_doplim[0], dopPOS);
			atomicMaxf(&p2d_doplim[1], dopPOS);

			/* Check if all Doppler bins which will receive power from this POS
			 * pixel fall within the data frame; if not, initialize the
			 * "overflow" spectrum if necessary.  */
			if ( (idop_min >= 1) && (idop_max <= p2d_ndop) )
				p2d_in_bounds = 1;
			else {
				p2d_in_bounds = 0;
				if (!p2d_any_overflow) {
					p2d_any_overflow = 1;
					for (j=0; j<MAXOVERFLOW; j++)
						p2d_fit_overflow[j] = 0.0;  // To-Do:  This might need attention.

					/* Center the COM in the overflow spectrum: bin [idop] in
					 * the fit frame corresponds to bin [idop+idop0] in the
					 * fit_overflow frame.  */
					p2d_idop0 = MAXOVERFLOW/2 - (int) floor(p2d_dopshift + 0.5);
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
					arg_bl = dopPOS + p2d_dopdiff_bl - idop;   /* bl = bottom left */
					sinc2_mean = (SINC2(arg_bl) +
							SINC2(arg_bl+p2d_xincr) +
							SINC2(arg_bl+p2d_yincr) +
							SINC2(arg_bl+p2d_xincr+p2d_yincr)) / 4;
					break;
				default:
					arg_left = dopPOS + p2d_dopdiff_bl - idop;
					sinc2_mean = 0.0;
					for (i=0; i<dpar->nsinc2; i++) {
						sinc2arg = arg_left;
						for (j=0; j<dpar->nsinc2; j++) {
							sinc2_mean += SINC2(sinc2arg);
							sinc2arg += p2d_xincr;
						}
						arg_left += p2d_yincr;
					}
					sinc2_mean /= p2d_nsinc2_sq;
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
					p2d_pos->cose_s[zaddr], p2d_pos->comp[x][y], p2d_pos->f[x][y])
		    		 * p2d_pos->km_per_pixel * p2d_pos->km_per_pixel  / sumweights;

			/* Only add POS pixel's power contributions to model Doppler spect-
			 * rum if NONE of those contributions fall outside spectrum limits*/
			if (p2d_in_bounds) {

				/*  Add the cross-section contributions to the model frame  */
				for (idop=idop_min; idop<=idop_max; idop++) {
					k = MIN( idop - idop_min, MAXBINS);
					fit_contribution = amp * dop_contribution[k];
					atomicAdd(&ddat->set[set].desc.doppler.frame[frame].fit_s[idop],
							fit_contribution);

					if (dpar->action == MAP) {
						if (dpar->map_mode == MAPMODE_DELDOP) {
							if (p2d_frame->map_fit[idop] > 0.0) {
								p2d_frame->map_pos[x][y] += fit_contribution;
								c = p2d_pos->comp[x][y];
								f = p2d_pos->f[x][y];
								p2d_frame->map_facet_power[c][f] += fit_contribution;
								if (dpar->map_verbose)
									printf("# POS (%3d, %3d) comp %d facet %4d contributes %e to Dop (%3d)\n",
											x+p2d_pos->n, y+p2d_pos->n, c, f, fit_contribution, idop-1);
							}
						} else if (dpar->map_mode == MAPMODE_POS) {
							if (p2d_frame->map_pos[x][y] > 0.0) {
								p2d_frame->map_fit[idop] += fit_contribution;
								c = p2d_pos->comp[x][y];
								f = p2d_pos->f[x][y];
								p2d_frame->map_facet_power[c][f] += fit_contribution;
								if (dpar->map_verbose)
									printf("# POS (%3d, %3d) comp %d facet %4d contributes %e to Dop (%3d)\n",
											x+p2d_pos->n, y+p2d_pos->n, c, f, fit_contribution, idop-1);
							}
						} else {
							if (p2d_frame->map_pos[x][y] > 0.0) {
								p2d_frame->map_fit[idop] += fit_contribution;
								if (dpar->map_verbose) {
									c = p2d_pos->comp[x][y];
									f = p2d_pos->f[x][y];
									printf("# POS (%3d, %3d) comp %d facet %4d contributes %e to Dop (%3d)\n",
											x+p2d_pos->n, y+p2d_pos->n, c, f, fit_contribution, idop-1);
								}
							}
						}
					}
				}
			} else {

				/*  Add the cross-section contributions to the "overflow" spectrum  */

				if (dpar->action == MAP && dpar->map_mode != MAPMODE_DELDOP)
					if (p2d_frame->map_pos[x][y] > 0.0)
						dpar->map_overflow = 1;

				idop1 = MAX( idop_min, -p2d_idop0);
				idop2 = MIN( idop_max, -p2d_idop0 + MAXOVERFLOW - 1);
				for (idop=idop1; idop<=idop2; idop++) {
					k = MIN( idop - idop_min, MAXBINS);
					fit_contribution = amp * dop_contribution[k];
					p2d_fit_overflow[idop+p2d_idop0] += fit_contribution;  // might need atomics
					if (dpar->action == MAP && dpar->map_mode == MAPMODE_DELDOP)
						if (idop >= dpar->map_doplim[0] && idop <= dpar->map_doplim[1]) {
							p2d_frame->map_pos[x][y] += fit_contribution;
							c = p2d_pos->comp[x][y];
							f = p2d_pos->f[x][y];
							p2d_frame->map_facet_power[c][f] += fit_contribution;
							if (dpar->map_verbose)
								printf("# POS (%3d, %3d) comp %d facet %4d contributes %e to Dop (%3d)\n",
										x+p2d_pos->n, y+p2d_pos->n, c, f, fit_contribution, idop-1);
						}
				}
			}
		}  /* if cos(scattering angle) > 0 and POS pixel projects onto the right body */
	}
}
__global__ void pos2doppler_finish_krnl(struct par_t *dpar, struct dat_t
		*ddat, int set, int frm) {
	/* Single-threaded kernel */
	int j, j1, j2;
	double lookfact, sdev_sq, variance, dopfactor;

	if (threadIdx.x ==0) {
		/* Copy float device variable over to the frame->doplim   */
		p2d_frame->doplim[0] = p2d_doplim[0];
		p2d_frame->doplim[1] = p2d_doplim[1];

		/* Convert model's Doppler limits from float bin numbers to Hz and
		 * widen the limits to account for nonzero POS pixel width    */

		p2d_frame->doplim[0] = (p2d_frame->doplim[0] - p2d_dopshift)*
				ddat->set[set].desc.doppler.dop_per_bin - p2d_dop_extra;
		p2d_frame->doplim[1] = (p2d_frame->doplim[1] - p2d_dopshift)*
				ddat->set[set].desc.doppler.dop_per_bin + p2d_dop_extra;

		/* Calculate overflow contributions to chi squared:
		 *   o2 = obs^2 contribution, m2 = model^2 contribution.
		 * Also compute summed cross section and mean Doppler bin for overflow
		 * region, for use with the "delcorinit" action   */

		p2d_frame->overflow_o2 = 0.0;
		p2d_frame->overflow_m2 = 0.0;
		p2d_frame->overflow_xsec = 0.0;
		p2d_frame->overflow_dopmean = 0.0;
		sdev_sq = p2d_frame->sdev*p2d_frame->sdev;
		variance = sdev_sq;
		lookfact = (p2d_frame->nlooks > 0.0) ? 1.0/p2d_frame->nlooks : 0.0;
		if (p2d_any_overflow) {
			j1 = MAX( p2d_frame->idoplim[0] + p2d_idop0, 0);
			j2 = MIN( p2d_frame->idoplim[1] + p2d_idop0, MAXOVERFLOW - 1);
			for (j=j1; j<=j2; j++) {
				if (p2d_fit_overflow[j] != 0.0) {
					if (dpar->speckle)
						variance = sdev_sq + lookfact*p2d_fit_overflow[j]*p2d_fit_overflow[j];
					p2d_frame->overflow_o2 += 1.0;
					p2d_frame->overflow_m2 += p2d_fit_overflow[j]*p2d_fit_overflow[j]/variance;
					p2d_frame->overflow_xsec += p2d_fit_overflow[j];
					p2d_frame->overflow_dopmean += (j - p2d_idop0)*p2d_fit_overflow[j];
				}
			}
			if (p2d_frame->overflow_xsec != 0.0)
				p2d_frame->overflow_dopmean /= p2d_frame->overflow_xsec;

			/*  Print a warning if the model extends even beyond the overflow spectrum  */

			if ( ((p2d_frame->idoplim[0] + p2d_idop0) < 0) ||
					((p2d_frame->idoplim[1] + p2d_idop0) >= MAXOVERFLOW) ) {
				p2d_badradar = 1;
				dopfactor = (MAX( p2d_frame->idoplim[1] + p2d_idop0, MAXOVERFLOW)
						- MIN( p2d_frame->idoplim[0] + p2d_idop0, 0)         )
		                								  / (1.0*MAXOVERFLOW);
				p2d_frame->badradar_logfactor += log(dopfactor);
				if (dpar->warn_badradar) {
					printf("\nWARNING in pos2doppler.c for set %2d frame %2d:\n", set, frm);
					printf("        model Doppler spectrum extends too far beyond the data spectrum\n");
					printf("             data:  bins %2d to %2d\n", 1, p2d_ndop);
					printf("            model:  bins %2d to %2d\n",
							p2d_frame->idoplim[0], p2d_frame->idoplim[1]);
				}
			}
		}
	}
}
__host__ int pos2doppler_cuda_2( struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, double orbit_xoff, double orbit_yoff, double
		orbit_dopoff, int body, int set, int frm, int v)
{
	int badradar, xlim0, xlim1, ylim0, ylim1, xspan, yspan, nThreads;
	dim3 BLK, THD;

	/* Launch single-threaded init kernel */
	pos2doppler_init_krnl<<<1,1>>>(ddat, set,frm, v);
	checkErrorAfterKernelLaunch("pos2doppler_init_krnl, line ");
	gpuErrchk(cudaMemcpyFromSymbol(&xlim0, p2d_xlim0, sizeof(xlim0),
				0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&xlim1, p2d_xlim1, sizeof(xlim1),
					0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&ylim0, p2d_ylim0, sizeof(ylim0),
					0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&ylim1, p2d_ylim1, sizeof(ylim1),
					0, cudaMemcpyDeviceToHost));

	/* Launch single-thread kernel for radar parameter calculation */
	pos2doppler_radar_parameters_krnl<<<1,1>>>(dpar, ddat, orbit_dopoff,
			set, frm, v);
	checkErrorAfterKernelLaunch("pos2doppler_radar_parameters_krnl, line ");

	xspan = xlim1 - xlim0 + 1;
	yspan = ylim1 - ylim0 + 1;
	nThreads = xspan * yspan;
	BLK.x = floor((maxThreadsPerBlock - 1 + nThreads)/maxThreadsPerBlock);
	THD.x = maxThreadsPerBlock; // Thread block dimensions

	pos2doppler_pixel_krnl<<<BLK,THD>>>(dpar, dmod, ddat, xspan, set,
			frm, nThreads, body, orbit_xoff, orbit_yoff);
	checkErrorAfterKernelLaunch("pos2doppler_pixel_krnl, line ");

	/* Launch the single-thread kernel to finish up Doppler calculations */
	pos2doppler_finish_krnl<<<1,1>>>(dpar, ddat, set, frm);
	checkErrorAfterKernelLaunch("pos2doppler_finish_krnl, line ");
	gpuErrchk(cudaMemcpyFromSymbol(&badradar, p2d_badradar, sizeof(badradar),
			0, cudaMemcpyDeviceToHost));

	return badradar;
}
