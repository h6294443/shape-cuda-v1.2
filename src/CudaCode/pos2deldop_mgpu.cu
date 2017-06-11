extern "C" {
#include "../shape/head.h"
}
/* TO DO: The static device float fit_overflow[MAXBINS][MAXBINS] is implemented
 * on each GPU. It must be recombined before the overflow kernel somehow. */
/* This multi-gpu version of pos2deldop cannot use file-scope __device__ variables.
 * (Actually it turns out that they can, but each gpu maintains its own copy.
 * I am leaving this code as it is - ME 6/1/17)
 * All previous __device__ variables are now contained in ddints, ddfloats,
 * dddoubles. Each GPU gets its own copy
 * 	ddints[0] - any_overflow
 * 	ddints[1] - codemethod
 * 	ddints[2] - spb
 * 	ddints[3] - stride
 * 	ddints[4] - spb_sq
 * 	ddints[5] - dopfftlen
 * 	ddints[6] - spb_over_stride
 * 	ddints[7] - nsinc2_sq
 *
 * 	ddfloats[0] - const1f
 * 	ddfloats[1] - const2f
 * 	ddfloats[2] - one_over_spbf
 * 	ddfloats[3] - delfactf
 * 	ddfloats[4] - dopfactf
 *
 * 	Doubles are similar.
 * 	Note that the previous static device float fit_overflow[MAX][MAX] where max
 * 	is MAXOVERFLOW (2000) is now two separate single-pointer arrays, unrolled
 * 	and split by device
 */

static __device__ float fit_overflow[MAXOVERFLOW][MAXOVERFLOW];
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

__global__ void pos2deldop_init_mgpu_krnl(
		struct dat_t *ddat,
		struct deldopfrm_t **frame,
		int *idel0,
		int *idop0,
		int *ndel,
		int *ndop,
		int set,
		int size,
		int oddflg,
		int *badradararr) {

	/* hf-threaded kernel. This kernel goes through all half-frames
	 * and assigns values.  Intended for multi-GPU operation. */
	int hf = blockIdx.x * blockDim.x + threadIdx.x;
	int f = 2 * hf + oddflg;	/* The actual frame */
	if (hf < size) {

		/*  Initialize variables to avoid compilation warnings  */
		idel0[hf] = idop0[hf] = 0;

		frame[hf] = &ddat->set[set].desc.deldop.frame[f];
		frame[hf]->idellim[0] = ndel[hf] + 999999;
		frame[hf]->idellim[1] = -999999;
		frame[hf]->idoplim[0] = ndop[hf] + 999999;
		frame[hf]->idoplim[1] = -999999;
		frame[hf]->dellim[0] =  HUGENUMBER;
		frame[hf]->dellim[1] = -HUGENUMBER;
		frame[hf]->doplim[0] =  HUGENUMBER;
		frame[hf]->doplim[1] = -HUGENUMBER;

		badradararr[hf] = 0;
		frame[hf]->badradar_logfactor = 0.0;
	}
}

__global__ void pos2deldop_data_sampling_mgpu_krnl(
		struct par_t *dpar,
		struct dat_t *ddat,
		struct deldopfrm_t **frame,
		struct pos_t **pos,
		float2 *axay,
		float2 *xyincr,
		float2 *deldopshift,
		float3 *w,
		float4 *dop,
		float4 *deldoplim,
		int4 *xylim,
		int *ddints,
		float *ddfloats,
		int set,
		int oddflg,
		int size,
		int v,
		double orbit_dopoff,
		int *badradararr) {

	/* multi-threaded kernel, goes through nfrm_half0 or nfrm_half1 threads */

	/* Get parameters related to data sampling and data reduction; then
	 * compute two more (both in units of delay bins = image rows):
	 *  	const1: half of the base width of the delay response function
	 *  	const2: half the delay difference between the first and last
	 * 				image rows within each baud  */

	int hf = blockDim.x * blockIdx.x + threadIdx.x;
	int f = 2 * hf + oddflg;

	/* First single-thread calculations */
	if (hf == 0) {
		ddints[0] = 0;										/* any_overflow		*/
		ddints[1] = ddat->set[set].desc.deldop.codemethod;	/* codemethod 		*/
		ddints[2] = ddat->set[set].desc.deldop.spb;			/* spb 				*/
		ddints[3] = ddat->set[set].desc.deldop.stride;		/* stride 			*/
		ddints[4] = ddints[2]*ddints[2];					/* spb_sq			*/
		ddints[5] = ddat->set[set].desc.deldop.dopfftlen;	/* dopfftlen		*/
		ddints[6] = ddints[2]/ddints[3];  					/* spb_over_stride	*/
		ddints[7] = dpar->nsinc2*dpar->nsinc2; 				/* nsinc2_sq		*/
		dop[hf].y = frame[hf]->dopDC_vig;

		if (ddints[1] != LONG_ORIG) {
			ddfloats[0] = (3*ddints[2] - 1)/(2.0*ddints[3]);
			ddfloats[1] = (ddints[2] - 1)/(2.0*ddints[3]);
		} else {
			ddfloats[0] = (double) (ddints[2]/ddints[3]);
			ddfloats[1] = 0.0;  /* not used with this code + reduction method */
		}

		/* Converts from km towards radar to delay bins. Note that ddfloats[4]
		 * assumes  that km_per_pixel remains the same for all frames */
		ddfloats[2] = 1.0/ddints[2];
		ddfloats[3] = -KM2US/ddat->set[set].desc.deldop.del_per_pixel;
		ddfloats[4] = ddat->set[set].desc.deldop.dopscale.val * KM2HZFACT *
				pos[0]->km_per_pixel * ddat->set[set].desc.deldop.Ftx /
				ddat->set[set].desc.deldop.dop_per_pixel;
	}
	__syncthreads();

	/* now onto multi-thread calculations */
	if (hf < size) {
		/*  Get w, the apparent spin vector in observer coordinates  */
		dev_cotrans4(w, frame[hf]->view[v].oe, frame[hf]->view[v].spin, 1, hf);

		/* Compute the Doppler bin increment per plane-of-sky pixel westward
		 * (ax) and northward (ay); these values are scaled by the "dopscale"
		 * parameter for this dataset. Then compute km2Hz, the Doppler
		 * increment (Hz) per km perpendicular to the projected spin axis in
		 * the plane of the sky.     */


		axay[hf].x = -w[hf].y * ddfloats[4];
		axay[hf].y =  w[hf].x * ddfloats[4];
		frame[hf]->view[v].km2Hz = sqrt(axay[hf].x*axay[hf].x+axay[hf].y*axay[hf].y) *
				ddat->set[set].desc.deldop.dop_per_pixel / pos[hf]->km_per_pixel;

		/* Compute the absolute value of the difference between the maximum (or minimum)
		 * Doppler on any given POS pixel's edge and the Doppler at its center             */
		if (w[hf].x != 0.0 || w[hf].y != 0.0)
			dop[hf].z = frame[hf]->view[v].km2Hz * 0.5 * pos[hf]->km_per_pixel
				* sqrt(w[hf].x * w[hf].x + w[hf].y * w[hf].y) /	MAX( fabs(w[hf].x),
						fabs(w[hf].y));
		else
			dop[hf].z = 0.0;

		/* We may be evaluating the sinc^2 Doppler response function at more than one point
		 * per POS pixel.  xincr and yincr are the Doppler bin increments between adjacent
		 * evaluation points in the x and y directions.  dopdiff_bl is the Doppler bin difference
		 * between the bottom-leftmost (southeasternmost) evaluation point and the pixel center.
		 * dopdiff_max is the maximum positive Doppler bin difference between any evaluation point
		 * and the pixel center.
		 * 	dop.w - dopdiff_bl
		 * 	dop.x - dopdiff_max
		 *	dop.y - dopDC_vig
		 * 	dop.z - dop_extra		 */

		xyincr[hf].x = axay[hf].x / dpar->nsinc2;
		xyincr[hf].y = axay[hf].y / dpar->nsinc2;
		dop[hf].w = -(dpar->nsinc2 - 1) * (xyincr[hf].x + xyincr[hf].y) / 2;
		dop[hf].x = (dpar->nsinc2 - 1) * (fabs(xyincr[hf].x) + fabs(xyincr[hf].y)) / 2;

		if (2 * dop[hf].x + dpar->sinc2width + 1 > MAXBINS) {
			badradararr[hf] = 1;
			frame[hf]->badradar_logfactor += log(
					(2 * dop[hf].x + dpar->sinc2width + 1) / MAXBINS);
			if (dpar->warn_badradar) {
				printf(
						"\nWARNING in pos2deldop_mgpu.c for set %2d frame %2d:\n",
						set, f);
				printf(
						"        sinc^2 function evaluated at %d Doppler bins, must be <= %d\n",
						(int) ceil(2 * dop[hf].x + dpar->sinc2width + 1),
						MAXBINS);
			}
		}

		/* Get the COM delay and Doppler bins, corrected for ephemeris drift and adjusted for
		 * orbital motion; the delay adjustment for orbital motion is done implicitly (the
		 * posvis routine has already adjusted the "z" values for all POS pixels), whereas the
		 * Doppler adjustment for orbital motion must be done here explicitly.                  */
		deldopshift[hf].x = frame[hf]->delcom_vig + frame[hf]->view[v].deloff;
		deldopshift[hf].y = frame[hf]->dopcom_vig + frame[hf]->view[v].dopoff
						+ orbit_dopoff;

		/* Now get pos->xlim[0], pos->xlim[1], pos->ylim[0], pos->ylim[1] */
		xylim[hf].w = pos[hf]->xlim[0];
		xylim[hf].x = pos[hf]->xlim[1];
		xylim[hf].y = pos[hf]->ylim[0];
		xylim[hf].z = pos[hf]->ylim[1];

		/* Copy frame[frm]->dellim[2] and frame[frm]->doplim[2] to the device variables */
		deldoplim[hf].w = __double2float_rd(frame[hf]->dellim[0]);
		deldoplim[hf].x = __double2float_rd(frame[hf]->dellim[1]);
		deldoplim[hf].y = __double2float_rd(frame[hf]->doplim[0]);
		deldoplim[hf].z = __double2float_rd(frame[hf]->doplim[1]);
	}
}

__global__ void pos2deldop_pixel_mgpu_krnl(
		struct par_t *dpar,
		struct mod_t *dmod,
		struct dat_t *ddat,
		struct pos_t **pos,
		struct deldopfrm_t **frame,
		float4 *deldoplim,
		float4 *dop,
		float2 *deldopshift,
		float2 *axay,
		float2 *xyincr,
		int *idel0,
		int *idop0,
		int *ndel,
		int *ndop,
		int *ddints,
		float *ddfloats,
		int xspan,
		int nThreads,
		int body,
		double orbit_xoff,
		double orbit_yoff,
		int set,
		int hf) {
	/* nThreads-threaded kernel */

	/*  Loop through all POS pixels within the rectangular plane-of-sky region spanned by the
	 *  model; for each such pixel which isn't blank sky, compute the cross-section contributions
	 *  to pixels in the model delay-Doppler frame. Note that functions posclr and posvis flag
	 *  blank-sky pixels by assigning "cose" = cos(scattering angle) = 0.
	 *  Only compute contributions from POS pixels that project onto the right body, in case this
	 *  is the "orbit" action (for which this routine is called twice, once for each of the two
	 *  orbiting bodies). */

	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int x = offset % xspan + pos[hf]->xlim[0];
	int y = offset / xspan + pos[hf]->ylim[0];
	int n;
	int fac, c, i, j, k, idel, idel_min, idel_max, idop_min, idop_max,
		idop, m, m_min, m_max, idel1, idel2, idop1, idop2, zaddr, in_bounds;

	double delPOS, dopPOS, codefactor, tmp, arg_sample, amp, arg_left,
		   sinc2arg, sinc2_mean, arg_bl, sumweights, fit_contribution;
	float del_contribution[MAXBINS], dop_contribution[MAXBINS];
	n = pos[hf]->n;
	if (offset < nThreads) {
		/* zaddr is the unrolled 1D pos->z_s[] array address  */
		zaddr = (y + n) * (2*n + 1) + (x + n);
		if (pos[hf]->cose_s[zaddr] > 0.0 && pos[hf]->body[x][y] == body) {

			/* Get the (floating-point) delay and Doppler bin of the POS pixel
			 * center: delPOS and dopPOS. Also get the minimum and maximum
			 * (integer) delay and Doppler bins to which this pixel contributes
			 * power: idel_min and idel_max, idop_min and idop_max. Each POS
			 * pixel contributes power to *all* Doppler columns, but here we're
			 * zeroing out the sinc^2 response function beyond the nearest
			 * sinc2width columns. Actually, if nsinc2 > 1, we'll distribute
			 * power to *at least* sinc2width Doppler bins: For pixels which
			 * span multiple bins we'll err on the side of computing more
			 * contributions rather          than fewer.       */

			delPOS = pos[hf]->z_s[zaddr] * ddfloats[3] + deldopshift[hf].x;
			idel_min = (int) floor(delPOS - ddfloats[0]) + 1;
			idel_max = (int) ceil(delPOS + ddfloats[0]) - 1;
			dopPOS = axay[hf].x*(x - orbit_xoff) + axay[hf].y*
					(y - orbit_yoff) + deldopshift[hf].y;
			idop_min = (int) floor(dopPOS - dop[hf].x + 1 - dpar->sinc2width/2.0);
			idop_max = (int) floor(dopPOS + dop[hf].x + dpar->sinc2width/2.0);

			/*  For the short code, sensitivity drops as we move away from DC. (This variation is slow,
			 *  so we can just evaluate the response at the center of the POS pixel.)
			 *  Note that the SINC2 macro multiplies its argument by pi.        */
			codefactor = (ddints[1] == SHORT) ? SINC2( (dopPOS-dop[hf].y)/
					ddints[5] ) : 1.0;

			/*  Update rectangular delay-Doppler region (row/column numbers) with !0 power according to model  */
			atomicMin(&frame[hf]->idellim[0], idel_min);
			atomicMax(&frame[hf]->idellim[1], idel_max);
			atomicMin(&frame[hf]->idoplim[0], idop_min);
			atomicMax(&frame[hf]->idoplim[1], idop_max);

			/*  Update the model's floating-point delay-Doppler limits, as determined prior to convolution
			 *  with the delay and Doppler response functions. At this point in the code, dellim and doplim
			 *  are pairs of floating-point row and column numbers which apply to POS pixel centers; when
			 *  the loop over POS pixels is finished we will convert them to usec and Hz, and will widen
			 *  the Doppler limits to account for nonzero POS pixel width.     */
			atomicMinf(&deldoplim[hf].w, (float) delPOS);
			atomicMaxf(&deldoplim[hf].x, (float) delPOS);
			atomicMinf(&deldoplim[hf].y, (float) dopPOS);
			atomicMaxf(&deldoplim[hf].z, (float) dopPOS);

			/*  Check whether or not all delay-Doppler pixels which will receive power from this POS pixel
			 *  fall within the data frame; if not, initialize the "overflow" image if necessary.         */
			if ((idel_min>=1) && (idel_max<=ndel[hf]) && (idop_min>=1) &&
					(idop_max<=ndop[hf]))
				in_bounds = 1;
			else {
				in_bounds = 0;
				if (!ddints[0]) {
					atomicExch(&ddints[0], 1);
					for (i=0; i<MAXOVERFLOW; i++)
						for (j=0; j<MAXOVERFLOW; j++)
							frame[hf]->fit_overflow[i][j] = 0.0;

					/* Center the COM in the overflow image:
					 * pixel [idel][idop] in the fit frame corresponds to
					 * pixel [idel+idel0][idop+idop0] in the fit_overflow frame*/

					idel0[hf] =MAXOVERFLOW/2-(int)floor(deldopshift[hf].x+0.5);
					idop0[hf] =MAXOVERFLOW/2-(int)floor(deldopshift[hf].y+0.5);
				}
			}

			/* Loop thru all delay bins this POS pixel contributes power (cross
			 * section), and compute delay response function for each bin       */
			for (idel=idel_min; idel<=idel_max; idel++) {
				if (ddints[1] != LONG_ORIG) {
					/* Get the delay response function for image row idel:
					 * sum the triangle-function contributions from each
					 * sample per baud, then divide the sum by spb and square.
					 * The triangle function for sample m  (0 <= m <= spb-1)
					 * has unit height and a half-width of spb/stride image
					 * rows,and is centered [-const2 + m/stride] rows later
					 * than the row center (idel).
					 * In the code block below, the arguments to macros TRI
					 * and TRI2 have been divided by half-width spb/stride,
					 * since those two macros are defined to give nonzero
					 * values for arguments between -1 and 1.  Each argument,
					 * then, is just (delPOS - [triangle-function center]) /
					 * half-width.
					 * Do the two most common cases (spb = 1 or 2) without
					 * loops in order to gain a bit of speed.  For the other
					 * cases, set m_min and m_max so as not to waste time
					 * computing contributions that are zero.             */

					switch (ddints[2]) {
					case 1:
						del_contribution[idel-idel_min] = TRI2( delPOS - idel );
						break;
					case 2:
						arg_sample = (delPOS - (idel - ddfloats[1])) /(ddints[6]);
						del_contribution[idel-idel_min] = TRI( arg_sample )
					                      				+ TRI( arg_sample - 0.5 );
						del_contribution[idel-idel_min] *= del_contribution[idel-idel_min]/4;
						break;
					default:
						del_contribution[idel-idel_min] = 0.0;
						m_min = MAX( (int) floor((delPOS - idel - ddfloats[1])
								* ddints[3]) , 0 );
						m_max = MIN( (int) ceil((delPOS - idel + ddfloats[0])
								* ddints[3]) , ddints[2] ) - 1;
						arg_sample = (delPOS - (idel - ddfloats[1])) /
								(ddints[2]/ddints[3]) -
								m_min*(ddfloats[2]);
						for (m=m_min; m<=m_max; m++) {
							del_contribution[idel-idel_min] += TRI( arg_sample );
							arg_sample -= (ddfloats[2]);
						}
						del_contribution[idel-idel_min] *=
								del_contribution[idel-idel_min]/(ddints[4]);
						break;
					}
				} else {

					/*  Long code with original (Harmon) reduction method: data for
					 *  each sample per baud are reduced separately,as if datataking
					 *  were just one sample per baud; then the image rows for spb/
					 *  stride samples are interleaved.  */
					del_contribution[idel-idel_min] = TRI2( (delPOS - idel) /
							(ddints[2]/ddints[3]) );
				}
			}

			/*  Next include the sinc^2 factor for Doppler mismatching: Take the
			 *  mean of nsinc2^2 points interior to the POS pixel. Do the two most
			 *  common cases (nsinc2 = 1 or 2) without loops in order to gain a bit
			 *  of speed. Note that the SINC2 macro multiplies its argument by pi*/

			for (idop=idop_min; idop<=idop_max; idop++) {
				switch (dpar->nsinc2) {
				case 1:
					sinc2_mean = SINC2( dopPOS - idop );
					break;
				case 2:
					arg_bl = dopPOS + dop[hf].w - idop;   /* bl = bottom left */
					sinc2_mean = ( SINC2( arg_bl ) +
							SINC2( arg_bl+xyincr[hf].x ) +
							SINC2( arg_bl+xyincr[hf].y ) +
							SINC2( arg_bl+xyincr[hf].x+xyincr[hf].y ) ) / 4;
					break;
				default:
					arg_left = dopPOS + dop[hf].w - idop;
					sinc2_mean = 0.0;
					for (i=0; i<dpar->nsinc2; i++) {
						sinc2arg = arg_left;
						for (j=0; j<dpar->nsinc2; j++) {
							sinc2_mean += SINC2( sinc2arg );
							sinc2arg += xyincr[hf].x;
						}
						arg_left += xyincr[hf].y;
					}
					sinc2_mean /= ddints[7];
					break;
				}
				k = MIN( idop - idop_min, MAXBINS);
				dop_contribution[k] = sinc2_mean;
			}

			/*  Compute the sum of delay-Doppler weighting factors  */
			sumweights = 0.0;
			for (idel=idel_min; idel<=idel_max; idel++)
				for (idop=idop_min; idop<=idop_max; idop++) {
					k = MIN( idop - idop_min, MAXBINS);
					sumweights += del_contribution[idel-idel_min]*dop_contribution[k];
				}

			/* The radar cross section within this plane-of-sky pixel is
			 * [differential radar scattering law]*[POS pixel area in km^2]
			 * The differential radar scattering law (function radlaw
			 * = d[cross section]/d[area] ) includes a sec(theta) factor to
			 * account for the fact that the POS pixel area is projected area
			 * rather than physical area on the target surface.      */

			amp = dev_radlaw( &dmod->photo, ddat->set[set].desc.deldop.iradlaw,
					pos[hf]->cose_s[zaddr], pos[hf]->comp[x][y], pos[hf]->f[x][y])
			       * pos[hf]->km_per_pixel * pos[hf]->km_per_pixel
			       * codefactor / sumweights;

			/* Only add this POS pixel's power contributions to model delay-
			 * Doppler frame if NONE of those contributions fall outside
			 * the frame limits.                                   */

			if (in_bounds) {

				/*  Add the cross-section contributions to the model frame  */
				for (idel=idel_min; idel<=idel_max; idel++)
					for (idop=idop_min; idop<=idop_max; idop++) {
						k = MIN( idop - idop_min, MAXBINS);
						fit_contribution = amp * del_contribution[idel-idel_min]
						                       * dop_contribution[k];

						atomicAdd(&frame[hf]->fit_s[(idop-1)*ndel[hf]+(idel-1)],
								(float)fit_contribution);
						if (dpar->action == MAP) {
							if (dpar->map_mode == MAPMODE_DELDOP) {
								if (frame[hf]->map_fit[idel][idop] > 0.0) {
									frame[hf]->map_pos[x][y] += fit_contribution;
									c = pos[hf]->comp[x][y];
									fac = pos[hf]->f[x][y];
									frame[hf]->map_facet_power[c][fac] += fit_contribution;
									if (dpar->map_verbose)
										printf("# POS (%3d, %3d) comp %d facet %4d contributes %e to d-D (%3d, %3d)\n",
												x+pos[hf]->n, y+pos[hf]->n, c, fac, fit_contribution, idel-1, idop-1);
								}
							} else if (dpar->map_mode == MAPMODE_POS) {
								if (frame[hf]->map_pos[x][y] > 0.0) {
									frame[hf]->map_fit[idel][idop] += fit_contribution;
									c = pos[hf]->comp[x][y];
									fac = pos[hf]->f[x][y];
									frame[hf]->map_facet_power[c][fac] += fit_contribution;
									if (dpar->map_verbose)
										printf("# POS (%3d, %3d) comp %d facet %4d contributes %e to d-D (%3d, %3d)\n",
												x+pos[hf]->n, y+pos[hf]->n, c, fac, fit_contribution, idel-1, idop-1);
								}
							} else {
								if (frame[hf]->map_pos[x][y] > 0.0) {
									frame[hf]->map_fit[idel][idop] += fit_contribution;
									if (dpar->map_verbose) {
										c = pos[hf]->comp[x][y];
										fac = pos[hf]->f[x][y];
										printf("# POS (%3d, %3d) comp %d facet %4d contributes %e to d-D (%3d, %3d)\n",
												x+pos[hf]->n, y+pos[hf]->n, c, fac, fit_contribution, idel-1, idop-1);
									}
								}
							}
						}
					}
			} else {

				/* Add the cross-section contributions to the "overflow" image */
				if (dpar->action == MAP && dpar->map_mode != MAPMODE_DELDOP)
					if (frame[hf]->map_pos[x][y] > 0.0)
						dpar->map_overflow = 1;
				idel1 = MAX( idel_min, -idel0[hf]);
				idel2 = MIN( idel_max, -idel0[hf] + MAXOVERFLOW - 1);
				idop1 = MAX( idop_min, -idop0[hf]);
				idop2 = MIN( idop_max, -idop0[hf] + MAXOVERFLOW - 1);
				for (idel=idel1; idel<=idel2; idel++)
					for (idop=idop1; idop<=idop2; idop++) {
						k = MIN( idop - idop_min, MAXBINS);
						fit_contribution = amp * del_contribution[idel-idel_min]
						                   * dop_contribution[k];
						fit_overflow[idel+idel0[hf]][idop+idop0[hf]] += fit_contribution;
						if (dpar->action == MAP && dpar->map_mode == MAPMODE_DELDOP)
							if (idel >= dpar->map_dellim[0] && idel <= dpar->map_dellim[1] &&
									idop >= dpar->map_doplim[0] && idop <= dpar->map_doplim[1]    ) {
								frame[hf]->map_pos[x][y] += fit_contribution;
								c = pos[hf]->comp[x][y];
								fac = pos[hf]->f[x][y];
								frame[hf]->map_facet_power[c][fac] += fit_contribution;
								if (dpar->map_verbose)
									printf("# POS (%3d, %3d) comp %d facet %4d contributes %e to d-D (%3d, %3d)\n",
											x+pos[hf]->n, y+pos[hf]->n, c, fac, fit_contribution, idel-1, idop-1);
							}
					}
			}
		}
	}
}

__global__ void pos2deldop_deldoplim_mgpu_krnl(
		struct dat_t *ddat,
		struct deldopfrm_t **frame,
		float4 *deldoplim,
		float4 *dop,
		float2 *deldopshift,
		int set,
		int size) {

	/* multi-threaded kernel. There are nfrm_half0 or nfrm_half1 threads to
	 * go through the half-frame arrays (multi-GPU application) */

	int hf = blockDim.x * blockIdx.x + threadIdx.x;
	float dlppxl, dpppxl;
	dlppxl = __double2float_rn(ddat->set[set].desc.deldop.del_per_pixel);
	dpppxl = __double2float_rn(ddat->set[set].desc.deldop.dop_per_pixel);

	if (hf < size) {
		frame[hf]->dellim[0] = deldoplim[hf].w;
		frame[hf]->dellim[1] = deldoplim[hf].x;
		frame[hf]->doplim[0] = deldoplim[hf].y;
		frame[hf]->doplim[1] = deldoplim[hf].z;

		/*  Convert the model's floating-point delay-Doppler limits from floating-
		 *  point row and column numbers to usec and Hz, and widen the Doppler limits
		 *  to account for nonzero POS pixel width  */
		frame[hf]->dellim[0] = (frame[hf]->dellim[0] - deldopshift[hf].x)*dlppxl;
		frame[hf]->dellim[1] = (frame[hf]->dellim[1] - deldopshift[hf].x)*dlppxl;
		frame[hf]->doplim[0] = (frame[hf]->doplim[0] - deldopshift[hf].y)*dpppxl
				- dop[hf].z;
		frame[hf]->doplim[1] = (frame[hf]->doplim[1] - deldopshift[hf].y)*dpppxl
				+ dop[hf].z;

	}
}

__global__ void pos2deldop_overflow_mgpu_krnl(
		struct par_t *dpar,
		struct deldopfrm_t **frame,
		int *idel0,
		int *idop0,
		int *ndel,
		int *ndop,
		int *ddints,
		int set,
		int oddflg,
		int size,
		int *badradararr) {
	/* nfrm_half0 or nfrm_half1-threaded kernel for for multi-GPU operation*/

	int hf = blockIdx.x * blockDim.x + threadIdx.x;
	int f = hf*2 + oddflg;
	int i, i1, i2, j, j1, j2;
	double lookfact, sdev_sq, variance,dopfactor, delfactor;

	if (hf < size) {
		/*  Calculate the overflow contributions to chi squared:
		 * 	 o2 = obs^2 contribution, m2 = model^2 contribution.
		 *
		 *  Also compute the summed cross section and the mean delay and Doppler
		 *  bins for the overflow region, for use with the "delcorinit" action    */

		frame[hf]->overflow_o2 = 0.0;
		frame[hf]->overflow_m2 = 0.0;
		frame[hf]->overflow_xsec = 0.0;
		frame[hf]->overflow_delmean = 0.0;
		frame[hf]->overflow_dopmean = 0.0;
		sdev_sq = frame[hf]->sdev*frame[hf]->sdev;
		variance = sdev_sq;
		lookfact = (frame[hf]->nlooks > 0.0) ? 1.0/frame[hf]->nlooks : 0.0;
		if (ddints[0]) {	/* If any overflow */
			i1 = MAX( frame[hf]->idellim[0] + idel0[hf], 0);
			i2 = MIN( frame[hf]->idellim[1] + idel0[hf], MAXOVERFLOW - 1);
			j1 = MAX( frame[hf]->idoplim[0] + idop0[hf], 0);
			j2 = MIN( frame[hf]->idoplim[1] + idop0[hf], MAXOVERFLOW - 1);
			for (i=i1; i<=i2; i++)
				for (j=j1; j<=j2; j++) {
					if (fit_overflow[i][j] != 0.0) {
						if (dpar->speckle)
							variance = sdev_sq + lookfact *
							frame[hf]->fit_overflow[i][j] * frame[hf]->fit_overflow[i][j];
						frame[hf]->overflow_o2 += 1.0;
						frame[hf]->overflow_m2 += frame[hf]->fit_overflow[i][j] *
								frame[hf]->fit_overflow[i][j]/variance;
						frame[hf]->overflow_xsec += frame[hf]->fit_overflow[i][j];
						frame[hf]->overflow_delmean += (i-idel0[hf]) *
								frame[hf]->fit_overflow[i][j];
						frame[hf]->overflow_dopmean += (j-idop0[hf]) *
								frame[hf]->fit_overflow[i][j];
					}
				}
			if (frame[hf]->overflow_xsec != 0.0) {
				frame[hf]->overflow_delmean /= frame[hf]->overflow_xsec;
				frame[hf]->overflow_dopmean /= frame[hf]->overflow_xsec;
			}

			/*  Print a warning if the model extends even beyond the overflow image  */

			if (((frame[hf]->idellim[0] + idel0[hf]) < 0)            ||
					((frame[hf]->idellim[1] + idel0[hf]) >= MAXOVERFLOW) ||
					((frame[hf]->idoplim[0] + idop0[hf]) < 0)            ||
					((frame[hf]->idoplim[1] + idop0[hf]) >= MAXOVERFLOW)    ) {

				badradararr[hf] = 1;
				delfactor = (MAX(frame[hf]->idellim[1] + idel0[hf], MAXOVERFLOW)
						- MIN(frame[hf]->idellim[0] + idel0[hf], 0)         )
		                				  / (1.0*MAXOVERFLOW);
				dopfactor = (MAX(frame[hf]->idoplim[1] + idop0[hf], MAXOVERFLOW)
						- MIN(frame[hf]->idoplim[0] + idop0[hf], 0)         )
		                				  / (1.0*MAXOVERFLOW);
				frame[hf]->badradar_logfactor += log(delfactor*dopfactor);
				if (dpar->warn_badradar) {
					printf("\nWARNING in pos2deldop_mgpu.cu for set %2d frame %2d:\n", set, f);
					printf("        model delay-Doppler image extends too far beyond the data image\n");
					printf("             data:  rows %2d to %2d , cols %2d to %2d\n", 1, ndel[hf], 1, ndop[hf]);
					printf("            model:  rows %2d to %2d , cols %2d to %2d\n",
							frame[hf]->idellim[0], frame[hf]->idellim[1],
							frame[hf]->idoplim[0], frame[hf]->idoplim[1]);
				}
			}
		}
	}
}

__host__ int pos2deldop_mgpu(
		struct par_t *dpar,
		struct mod_t *dmod,
		struct dat_t *ddat,
		struct pos_t **pos0,			/* Even frames on GPU0 */
		struct pos_t **pos1,			/* Odd frames on GPU1  */
		struct deldopfrm_t **frame0,	/* Even frames on GPU0 */
		struct deldopfrm_t **frame1,	/* Odd frames on GPU1  */
		int4 *xylim0,
		int4 *xylim1,
		int *ndel0,
		int *ndel1,
		int *ndop0,
		int *ndop1,
		double orbit_xoff,
		double orbit_yoff,
		double orbit_dopoff,
		int body,
		int set,
		int nframes,
		int v,
		int *badradararr0,
		int *badradararr1,
		cudaStream_t *gpu0_stream,
		cudaStream_t *gpu1_stream)
{
	/* Calculate the half-frame values real quick: */
	int nfrm_half0 = nframes/2 + nframes%2;
	int nfrm_half1 = nframes/2;
	int hf = 0;		/* Used as array index for the split arrays */

	int xspan[nframes], yspan, nThreads[nframes], *idop0_0, *idop0_1,
		*idel0_0, *idel0_1, *ddints0, *ddints1;
	float *ddfloats0, *ddfloats1;
	float2 *axay0, *axay1, *xyincr0, *xyincr1, *deldopshift0, *deldopshift1;
	float3 *w0, *w1;
	float4 *dop0, *dop1, *deldoplim0, *deldoplim1;
	int4 hxylim0[nfrm_half0], hxylim1[nfrm_half1];
	dim3 BLK[nframes], THD, BLKhalf0, BLKhalf1;
	THD.x = maxThreadsPerBlock;
	BLKhalf0.x = floor ((THD.x - 1 + nfrm_half0))/THD.x;
	BLKhalf1.x = floor ((THD.x - 1 + nfrm_half1))/THD.x;

	/* Allocate GPU0 arrays (even frames) */
	gpuErrchk(cudaSetDevice(GPU0));
	gpuErrchk(cudaMalloc((void**)&ddints0, sizeof(int)*8));
	gpuErrchk(cudaMalloc((void**)&ddfloats0, sizeof(float)*8));
	gpuErrchk(cudaMalloc((void**)&idop0_0, sizeof(int)*nfrm_half0));
	gpuErrchk(cudaMalloc((void**)&idel0_0, sizeof(int)*nfrm_half0));
	gpuErrchk(cudaMalloc((void**)&w0, sizeof(float3)*nfrm_half0));
	gpuErrchk(cudaMalloc((void**)&axay0, sizeof(float2)*nfrm_half0));
	gpuErrchk(cudaMalloc((void**)&xyincr0, sizeof(float2)*nfrm_half0));
	gpuErrchk(cudaMalloc((void**)&deldopshift0, sizeof(float2)*nfrm_half0));
	gpuErrchk(cudaMalloc((void**)&deldoplim0, sizeof(float4)*nfrm_half0));
	gpuErrchk(cudaMalloc((void**)&dop0, sizeof(float4)*nfrm_half0));

	/* Allocate GPU1 arrays (even frames), then switch back to GPU0 */
	gpuErrchk(cudaSetDevice(GPU1));
	gpuErrchk(cudaMalloc((void**)&ddints1, sizeof(int)*8));
	gpuErrchk(cudaMalloc((void**)&ddfloats1, sizeof(float)*8));
	gpuErrchk(cudaMalloc((void**)&idop0_1, sizeof(int)*nfrm_half1));
	gpuErrchk(cudaMalloc((void**)&idel0_1, sizeof(int)*nfrm_half1));
	gpuErrchk(cudaMalloc((void**)&w1, sizeof(float3)*nfrm_half1));
	gpuErrchk(cudaMalloc((void**)&axay1, sizeof(float2)*nfrm_half1));
	gpuErrchk(cudaMalloc((void**)&xyincr1, sizeof(float2)*nfrm_half1));
	gpuErrchk(cudaMalloc((void**)&deldopshift1, sizeof(float2)*nfrm_half1));
	gpuErrchk(cudaMalloc((void**)&deldoplim1, sizeof(float4)*nfrm_half1));
	gpuErrchk(cudaMalloc((void**)&dop1, sizeof(float4)*nfrm_half1));
	gpuErrchk(cudaSetDevice(GPU0));

	/* Launch single-thread initialization kernels for both GPUs. Even frames first.*/
	pos2deldop_init_mgpu_krnl<<<BLKhalf0,THD,0,gpu0_stream[0]>>>(ddat, frame0,
			idel0_0, idop0_0, ndel0, ndop0, set, nfrm_half0, 0, badradararr0);
	pos2deldop_data_sampling_mgpu_krnl<<<BLKhalf0,THD,0,gpu0_stream[0]>>>(dpar, ddat,
			frame0, pos0, axay0, xyincr0, deldopshift0, w0, dop0,
			deldoplim0, xylim0, ddints0, ddfloats0, set, 0, nfrm_half0,
			v, orbit_dopoff, badradararr0);
	gpuErrchk(cudaMemcpyAsync(hxylim0, xylim0, nfrm_half0*sizeof(int4),
			cudaMemcpyDeviceToHost, gpu0_stream[0]));

	/* Switch to GPU1 and do odd frames */
	gpuErrchk(cudaSetDevice(GPU1));
	pos2deldop_init_mgpu_krnl<<<BLKhalf1,THD,0,gpu1_stream[0]>>>(ddat, frame1,
			idel0_1, idop0_1, ndel1, ndop1, set, nfrm_half1, 1,	badradararr1);
	pos2deldop_data_sampling_mgpu_krnl<<<BLKhalf1,THD,0,gpu1_stream[0]>>>(dpar, ddat,
			frame1, pos1, axay1, xyincr1, deldopshift1, w1, dop1,
			deldoplim1, xylim1, ddints1, ddfloats1, set, 0, nfrm_half1,
			v, orbit_dopoff, badradararr1);
	gpuErrchk(cudaMemcpyAsync(hxylim1, xylim1, nfrm_half1*sizeof(int4),
			cudaMemcpyDeviceToHost, gpu1_stream[0]));

	checkErrorAfterKernelLaunch("pos2deldop_mgpu first two kernels");

	gpuErrchk(cudaSetDevice(GPU0));

	/* Figure out the kernel launch parameters for every frame, both GPUs */
	hf = 0;
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
	int f = 0;
	for (int hf=0; hf<nfrm_half0; hf++) {
		pos2deldop_pixel_mgpu_krnl<<<BLK[f],THD,0,gpu0_stream[hf]>>>(dpar,dmod,
				ddat, pos0, frame0, deldoplim0, dop0, deldopshift0, axay0,
				xyincr0, idel0_0, idop0_0, ndel0, ndop0, ddints0, ddfloats0,
				xspan[f], nThreads[f], body, orbit_xoff, orbit_yoff, set, hf);
		f += 2;
	}
	/* Now GPU1 */
	gpuErrchk(cudaSetDevice(GPU1));
	f = 1;
	for (int hf=0; hf<nfrm_half1; hf++) {
		pos2deldop_pixel_mgpu_krnl<<<BLK[f],THD,0,gpu1_stream[hf]>>>(dpar,dmod,
				ddat, pos1, frame1, deldoplim1, dop1, deldopshift1, axay1,
				xyincr1, idel0_1, idop0_1, ndel1, ndop1, ddints1, ddfloats1,
				xspan[f], nThreads[f], body, orbit_xoff, orbit_yoff, set, hf);
		f += 2;
	}

	/* Launch kernel to copy the deldop limits back to original doubles in
	 * the frame structures.	 */
	/* GPU1 first */
	pos2deldop_deldoplim_mgpu_krnl<<<BLKhalf1,THD,0,gpu1_stream[0]>>>(ddat, frame1,
			deldoplim1, dop1, deldopshift1, set, nfrm_half1);
	/* Launch kernel to take care of any bin overflow */
	pos2deldop_overflow_mgpu_krnl<<<BLKhalf1,THD,0,gpu1_stream[0]>>>(dpar, frame1,
			idel0_1, idop0_1, ndel1, ndop1, ddints1, set, 0, nfrm_half1,
			badradararr1);

	/* GPU0 for even frames next */
	gpuErrchk(cudaSetDevice(GPU0));
	pos2deldop_deldoplim_mgpu_krnl<<<BLKhalf0,THD,0,gpu0_stream[0]>>>(ddat, frame0,
				deldoplim0, dop0, deldopshift0, set, nfrm_half0);
	/* Launch kernel to take care of any bin overflow */
	pos2deldop_overflow_mgpu_krnl<<<BLKhalf0,THD,0,gpu0_stream[0]>>>(dpar, frame0,
			idel0_0, idop0_0, ndel0, ndop0, ddints0, set, 0, nfrm_half0,
			badradararr0);

	/* Free GPU0 memory */
	gpuErrchk(cudaFree(ddints0));
	gpuErrchk(cudaFree(ddfloats0));
	gpuErrchk(cudaFree(idop0_0));
	gpuErrchk(cudaFree(idel0_0));
	gpuErrchk(cudaFree(w0));
	gpuErrchk(cudaFree(axay0));
	gpuErrchk(cudaFree(xyincr0));
	gpuErrchk(cudaFree(deldopshift0));
	gpuErrchk(cudaFree(deldoplim0));
	gpuErrchk(cudaFree(dop0));

	/* Free GPU1 memory */
	gpuErrchk(cudaSetDevice(GPU1));
	gpuErrchk(cudaFree(ddints1));
	gpuErrchk(cudaFree(ddfloats1));
	gpuErrchk(cudaFree(idop0_1));
	gpuErrchk(cudaFree(idel0_1));
	gpuErrchk(cudaFree(w1));
	gpuErrchk(cudaFree(axay1));
	gpuErrchk(cudaFree(xyincr1));
	gpuErrchk(cudaFree(deldopshift1));
	gpuErrchk(cudaFree(deldoplim1));
	gpuErrchk(cudaFree(dop1));

	gpuErrchk(cudaSetDevice(GPU0));
}
