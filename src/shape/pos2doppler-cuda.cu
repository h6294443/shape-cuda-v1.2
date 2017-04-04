/*****************************************************************************************
                                                                                 posvis.c

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

This is the CUDA-enabled (and CUDA-only!) version of pos2doppler
 *****************************************************************************************/
extern "C" {
#include "../shape/head.h"
#include "../shape/shape-cuda.h"
}

/*	The following __device__ kernel is associated with the posvis_cuda() routine and the
 * 	__global__ posvis_facets kernel.  The kernel accomplishes roughly the same thing as
 * 	the two nested for-loops (i1 < i < i2, j1 < j < j2) in the original posvis single CPU
 * 	routine.
 * 	This kernel is called from another kernel - posvis_facets.  As such, it implements
 * 	dynamic parallelism and requires CUDA Compute Capability 3.5 to work.
 *
 * 	Input arguments:
 * 		*pos 		   = Pointer to address of pos substructure of *frame.
 *		*frame 		   = the current radar frame
 *		*fit_overflow  = array of floats for overflow calculation
 *		*radar		   = photo->radar
 *		*radtype	   = photo->radtype
 *		*counter	   = test variable for debugging
 * 		kernel_ints[0] = par->sinc2width;
 *		kernel_ints[1] = par->nsinc2;
 *		kernel_ints[2] = par->map_doplim[0];
 *		kernel_ints[3] = par->map_doplim[1];
 *		kernel_ints[4] = body;
 *		kernel_ints[5] = doppler->iradlaw;
 *		kernel_ints[6] = 0;		// init for the old idop0[0] // Note that idop0[0] = idop0
 *		kernel_ints[7] = 0;		// init for the old idop0[1] // and idop[1] = any_overflow
 *		kernel_ints[8] = frame->ndop;
 *
 *		kernel_uchars[0] = par->action;
 *		kernel_uchars[1] = par->map_mode;
 *		kernel_uchars[2] = par->map_verbose;
 *	 	kernel_uchars[3] = par->map_overflow;
 *
 *		kernel_doubles[0] = ax
 *		kernel_doubles[1] = ay
 *		kernel_doubles[2] = orbit_xoff;
 *		kernel_doubles[3] = orbit_yoff;
 *		kernel_doubles[4] = dopshift
 *		kernel_doubles[5] = dopdiff_max
 *		kernel_doubles[6] = dopdiff_bl
 *		kernel_doubles[7] = xincr
 *		kernel_doubles[8] = yincr
 *
 *		*fitf			  = float array of frame->fit values for atomic add  */

/*	Note that:  xlim[0] = positive limit ( 75)
 * 				xlim[1] = negative limit (-75)
 * 				ylim[0] = positive limit ( 75)
 * 				ylim[1] = negative limit (-75)	*/

__global__ void pos2doppler_pixel_kernel(struct pos_t *pos, struct dopfrm_t *frame,
		float *fit_overflow, union radscat_t *radar, unsigned char *radtype,
		int *kernel_ints, unsigned char *kernel_uchars,	double *kernel_doubles,
		float *fitf){

	/*	Find the facet index with thread and block indices.  First are global
	 * 	thread indices in x and y directions.  Then total offset.  Finally i
	 * 	and j which index i1 <= i <= i2 and j1 <= j <= j2	*/
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int x = pos->xlim[0] + (offset % (pos->xlim[1]-pos->xlim[0]));
	int y = pos->ylim[0] + (offset / (pos->xlim[1]-pos->xlim[0]));
	int idop, idop_min, idop_max, idop1, idop2;
	int i, j, k, c, f, in_bounds;
	double dopPOS, arg_bl, arg_left, sinc2arg, sinc2_mean;
	double dop_contribution[MAXBINS], tmp, sumweights, amp, fit_contribution;

	if (offset < (pos->xlim[0] - pos->xlim[1]) * (pos->ylim[0] - pos->ylim[1])){

		if (pos->cose[x][y] > 0.0 && pos->body[x][y] == kernel_ints[4]) {

			/*	Get the (floating-point) Doppler bin of the POS pixel center: dopPOS.
			 *  Also get the minimum and maximum (integer) Doppler bins to which this
			 *  pixel contributes power: idop_min and idop_max.  Strictly speaking, each
			 *  POS pixel contributes power to *all* Doppler bins, but here we're zeroing
			 *  out the sinc^2 response function beyond the nearest sinc2width bins.
			 *  Actually, if nsinc2 > 1, we'll distribute power to *at least* sinc2width
			 *  Doppler bins: For pixels which span multiple bins we'll err on the side
			 *  of computing more contributions rather than fewer.		 */

			dopPOS = kernel_doubles[0]*(x - kernel_doubles[2]) + kernel_doubles[1]*
					(y - kernel_doubles[3]) + kernel_doubles[4];
			idop_min = (int) floor(dopPOS - kernel_doubles[5] + 1 - kernel_ints[0]/2.0);
			idop_max = (int) floor(dopPOS + kernel_doubles[5] + kernel_ints[0]/2.0);

			/*  Update the rectangular delay-Doppler region with nonzero power according to the model              */
			if (idop_min < frame->idoplim[0])	frame->idoplim[0] = idop_min;
			if (idop_max > frame->idoplim[1])	frame->idoplim[1] = idop_max;

			/*	Update the model's floating-point Doppler limits, as determined prior to
			 * 	convolution with the Doppler response function. At this point in the code, doplim
			 * 	is a pair of floating-point bin numbers which applies to POS pixel centers; when the
			 * 	loop over POS pixels is finished we will convert these limits to Hz and will widen the
			 * 	limits to account for nonzero POS pixel width.	 */

			// To-Do:  This should also be atomic
			frame->doplim[0] = MIN( frame->doplim[0], dopPOS);
			frame->doplim[1] = MAX( frame->doplim[1], dopPOS);

			/*  Check whether or not all Doppler bins which will receive
		    power from this POS pixel fall within the data frame;
		    if not, initialize the "overflow" spectrum if necessary.  */

			if ( (idop_min >= 1) && (idop_max <= kernel_ints[8]) )
				in_bounds = 1;
			else {
				in_bounds = 0;
				if (!kernel_ints[7]) {
					kernel_ints[7] = 1;
					for (j=0; j<MAXOVERFLOW; j++)
						fit_overflow[j] = 0.0;

					/*	Center the COM in the overflow spectrum:
		            bin [idop] in the fit frame corresponds to
		            bin [idop+idop0] in the fit_overflow frame.  */
					kernel_ints[6] = MAXOVERFLOW/2 - (int) floor(kernel_doubles[4] + 0.5);
				}
			}

			/*  Compute the sinc^2 factors for Doppler mismatching:
            Take the mean of nsinc2^2 points interior to the POS pixel.
            Do the two most common cases (nsinc2 = 1 or 2) without
            loops in order to gain a bit of speed.  Note that the SINC2
           	macro multiplies its argument by pi.
            Then add the cross-section contributions to the model spectrum.  */
			sumweights = 0.0;	// sum of Doppler weighting factors
			for (idop=idop_min; idop<=idop_max; idop++) {
				switch (kernel_ints[1]) {
				case 1:
					sinc2_mean = SINC2(dopPOS - idop);
					break;
				case 2:
					arg_bl = dopPOS + kernel_doubles[6] - idop;   /* bl = bottom left */
					sinc2_mean = (SINC2(arg_bl) +
							SINC2(arg_bl+kernel_doubles[7]) +
							SINC2(arg_bl+kernel_doubles[8]) +
							SINC2(arg_bl+kernel_doubles[7]+kernel_doubles[8])) / 4;
					break;
				default:
					arg_left = dopPOS + kernel_doubles[6] - idop;
					sinc2_mean = 0.0;
					for (i=0; i<kernel_ints[1]; i++) {
						sinc2arg = arg_left;
						for (j=0; j<kernel_ints[1]; j++) {
							sinc2_mean += SINC2(sinc2arg);
							sinc2arg += kernel_doubles[7];
						}
						arg_left += kernel_doubles[8];
					}
					sinc2_mean /= kernel_ints[1] * kernel_ints[1];
					break;
				}
				k = MIN( idop - idop_min, MAXBINS);
				dop_contribution[k] = sinc2_mean;  //To-Do-needs atomic
				sumweights += dop_contribution[k];
			}

			/*	Radar cross section within this plane-of-sky pixel is [differential radar scattering
			 * 	law]*[POS pixel area in km^2]. Differential radar scattering law (function radlaw =
			 * 	d[cross section]/d[area] ) includes a sec(theta) factor to account for the fact that
			 * 	the POS pixel area is projected area rather than physical area on the target surface.	 */

			amp = radlaw_cuda(radar, radtype, kernel_ints[5], pos->cose[x][y],	pos->comp[x][y],
					pos->f[x][y]) * pos->km_per_pixel * pos->km_per_pixel / sumweights;

			/*  Only add this POS pixel's power contributions to the model Doppler spectrum if NONE of
			 * 	those contributions fall outside the spectrum limits.                                 */
			if (in_bounds) {
				/*  Add the cross-section contributions to the model frame  */
				for (idop=idop_min; idop<=idop_max; idop++) {
					k = MIN( idop - idop_min, MAXBINS);
					fit_contribution = amp * dop_contribution[k];
					atomicAdd(&fitf[idop], (float)fit_contribution);

					if (kernel_uchars[0] == MAP) {					// [0] = par->action
						if (kernel_uchars[1] == MAPMODE_DELDOP) { 	// [1] = par->map_mode
							if (frame->map_fit[idop] > 0.0) {
								frame->map_pos[x][y] += fit_contribution;
								c = pos->comp[x][y];
								f = pos->f[x][y];
								frame->map_facet_power[c][f] += fit_contribution; // needs proper atomic add
								if (kernel_uchars[2])				// [2] = par->map_verbose
									printf("# POS (%3d, %3d) comp %d facet %4d contributes %e to Dop (%3d)\n",
											x+pos->n, y+pos->n, c, f, fit_contribution, idop-1);
							}
						} else if (kernel_uchars[1] == MAPMODE_POS) { //[1] = par->map_mode
							if (frame->map_pos[x][y] > 0.0) {
								frame->map_fit[idop] += fit_contribution; // needs proper atomic add
								c = pos->comp[x][y];
								f = pos->f[x][y];
								frame->map_facet_power[c][f] += fit_contribution;  // needs proper atomic add
								if (kernel_uchars[2])					// [2] = par->map_verbose
									printf("# POS (%3d, %3d) comp %d facet %4d contributes %e to Dop (%3d)\n",
											x+pos->n, y+pos->n, c, f, fit_contribution, idop-1);
							}
						} else {
							if (frame->map_pos[x][y] > 0.0) {
								frame->map_fit[idop] += fit_contribution; // needs proper atomic add
								if (kernel_uchars[2]) {
									c = pos->comp[x][y];
									f = pos->f[x][y];
									printf("# POS (%3d, %3d) comp %d facet %4d contributes %e to Dop (%3d)\n",
											x+pos->n, y+pos->n, c, f, fit_contribution, idop-1);
								}
							}
						}
					}
				}
			} else {
				/*  Add the cross-section contributions to the "overflow" spectrum  */
				if (kernel_uchars[0] == MAP && kernel_uchars[1] != MAPMODE_DELDOP)	// [0] = par->action
					if (frame->map_pos[x][y] > 0.0)									// [1] = par->map_mode
						kernel_uchars[3] = 1;										// [3] = par->map_overflow

				idop1 = MAX( idop_min, -kernel_ints[6]);
				idop2 = MIN( idop_max, -kernel_ints[6] + MAXOVERFLOW - 1);
				for (idop=idop1; idop<=idop2; idop++) {
					k = MIN( idop - idop_min, MAXBINS);
					fit_contribution = amp * dop_contribution[k];
					atomicAdd(&fit_overflow[idop+kernel_ints[6]], (float)fit_contribution);
					//fit_overflow[idop+kernel_ints[6]] += fit_contribution;
					if (kernel_uchars[0] == MAP && kernel_uchars[1] == MAPMODE_DELDOP)	// [0] = par->action
						if (idop >= kernel_ints[2] && idop <= kernel_ints[3]) {
							frame->map_pos[x][y] += fit_contribution;
							c = pos->comp[x][y];
							f = pos->f[x][y];
							frame->map_facet_power[c][f] += fit_contribution;  // needs proper atomic add
							if (kernel_uchars[2])
								printf("# POS (%3d, %3d) comp %d facet %4d contributes %e to Dop (%3d)\n",
										x+pos->n, y+pos->n, c, f, fit_contribution, idop-1);
						} // end if (idop >= kernel_ints[2]....
				} // end for (idop = idop1; idop<=diop2.....
			} // end else
		}// end if (pos->cose[x][y] > 0.0 && pos->body[x][y] == body)
	}// end if (offset < (pos->xlim[1] - pos->xlim[0]) * (pos->ylim[1] - pos->ylim[0]))
}


int pos2doppler_cuda( struct par_t *par, struct photo_t *photo,
		double orbit_xoff, double orbit_yoff, double orbit_dopoff,
		struct doppler_t *doppler, int body, int set, int frm, int v)
{
	int j, j1, j2, badradar, any_overflow, *kernel_ints;
	double dopfact, w[3], lookfact, sdev_sq, variance, dop_extra, dopfactor,
			*kernel_doubles;
	float *fitf, *fit_overflow;
	struct dopfrm_t *frame;
	struct pos_t *pos;
	union radscat_t *radar;
	unsigned char *radtype, *kernel_uchars;

	/*	Allocate device memory for data passed to GPU	*/
	cudaCalloc1((void**)&frame, 			sizeof(struct dopfrm_t),	1);
	cudaCalloc1((void**)&pos, 			sizeof(struct pos_t), 		1);
	cudaCalloc1((void**)&fit_overflow,	sizeof(double),   MAXOVERFLOW);
	cudaCalloc1((void**)&radar, 			sizeof(union radscat_t), photo->nradlaws);
	cudaCalloc1((void**)&radtype,		sizeof(unsigned char),   photo->nradlaws);
	cudaCalloc1((void**)&kernel_ints, 	sizeof(int), 				9);
	cudaCalloc1((void**)&kernel_doubles, sizeof(double), 			9);
	cudaCalloc1((void**)&kernel_uchars, 	sizeof(unsigned char), 		4);
	cudaCalloc1((void**)&fitf, 			sizeof(float), doppler->frame[frm].ndop);

	/*	Variable initialization section	*/
	frame = &doppler->frame[frm];
	pos = &frame->pos;
	radtype = photo->radtype;
	radar = photo->radar;
	frame->idoplim[0] = frame->ndop + 999999;
	frame->idoplim[1] = -999999;
	frame->doplim[0] =  HUGENUMBER;
	frame->doplim[1] = -HUGENUMBER;
	badradar = 0;
	frame->badradar_logfactor = 0.0;

	kernel_ints[0] = par->sinc2width;
	kernel_ints[1] = par->nsinc2;
	kernel_ints[2] = par->map_doplim[0];
	kernel_ints[3] = par->map_doplim[1];
	kernel_ints[4] = body;
	kernel_ints[5] = doppler->iradlaw;
	kernel_ints[6] = 0;		// init for the old idop0[0] // Note that idop0[0] = idop0
	kernel_ints[7] = 0;		// init for the old idop0[1] // and idop[1] = any_overflow
	kernel_ints[8] = frame->ndop;

	kernel_uchars[0] = par->action;
	kernel_uchars[1] = par->map_mode;
	kernel_uchars[2] = par->map_verbose;
	kernel_uchars[3] = par->map_overflow;

	kernel_doubles[0] = 0.0;		// just for initialization (this is ax)
	kernel_doubles[1] = 0.0;		// just for initialization (this is ay)
	kernel_doubles[2] = orbit_xoff;
	kernel_doubles[3] = orbit_yoff;
	kernel_doubles[4] = 0.0;		// initialization of dopshift
	kernel_doubles[5] = 0.0;		// initialization of dopdiff_max
	kernel_doubles[6] = 0.0;		// initialization of dopdiff_bl
	kernel_doubles[7] = 0.0;		// initialization of xincr
	kernel_doubles[8] = 0.0;		// initialization of yincr
	/*	End variable initialization section	*/

	/*  Get w, the apparent spin vector in observer coordinates  */
	cotrans( w, frame->view[v].oe, frame->view[v].spin, 1);

	/*	Compute the Doppler bin increment per plane-of-sky pixel westward (ax) and
	 *  northward (ay); these values are scaled by the "dopscale" parameter for this
	 *  dataset.  Then compute km2Hz, the Doppler increment (Hz) per km perpendicular
	 *  to the projected spin axis in the plane of the sky.     */

	dopfact = doppler->dopscale.val * KM2HZFACT * pos->km_per_pixel
			* doppler->Ftx / doppler->dop_per_bin;
	kernel_doubles[0] = -w[1]*dopfact;		// this is ax
	kernel_doubles[1] =  w[0]*dopfact;		// this is ay
	frame->view[v].km2Hz = sqrt(kernel_doubles[0]*kernel_doubles[0]
	          + kernel_doubles[1]*kernel_doubles[1]) * doppler->dop_per_bin
			/ pos->km_per_pixel;

	/* Compute the absolute value of the difference between the maximum (or minimum)
	 * Doppler on any given POS pixel's edge and the Doppler at its center            */

	if (w[0] != 0.0 || w[1] != 0.0)
		dop_extra = frame->view[v].km2Hz * 0.5 * pos->km_per_pixel
		* sqrt(w[0]*w[0] + w[1]*w[1]) / MAX( fabs(w[0]), fabs(w[1]));
	else
		dop_extra = 0.0;

	/*  We may be evaluating the sinc^2 Doppler response function at more than one
	 * 	point per POS pixel.  xincr and yincr are the Doppler bin increments between
	 * 	adjacent evaluation points in the x and y directions.  dopdiff_bl is the Doppler
	 * 	bin difference between the bottom-leftmost (southeasternmost) evaluation point
	 * 	and the pixel center.  dopdiff_max is the maximum positive Doppler bin difference
	 * 	between any evaluation point and the pixel center.                        */

	kernel_doubles[7] = kernel_doubles[0] / kernel_ints[1];	// ax / par->nsinc2 and doubles[7] = xincr
	kernel_doubles[8] = kernel_doubles[1] / kernel_ints[1];	// ay / par->nsinc2 and doubles[8] = yincr
	kernel_doubles[6] = -(kernel_ints[1] - 1)*(kernel_doubles[7] + kernel_doubles[8])/2;	// kernel_doubles[6] = dopdiff_bl
	kernel_doubles[5] = (kernel_ints[1] - 1)*(fabs(kernel_doubles[7]) + fabs(kernel_doubles[8]))/2;  // kernel_doubles[5] = dopdiff_max
	if (2*kernel_doubles[5] + kernel_ints[0] + 1 > MAXBINS) {
		badradar = 1;
		frame->badradar_logfactor += log((2*kernel_doubles[5] + kernel_ints[0] + 1) / MAXBINS);
		if (par->warn_badradar) {
			printf("\nWARNING in pos2doppler.c for set %2d frame %2d:\n", set, frm);
			printf("        sinc^2 function evaluated at %d Doppler bins, must be <= %d\n",
					(int) ceil(2*kernel_doubles[5] + kernel_ints[0] + 1), MAXBINS);
			fflush(stdout);
		}
	}

	/*  Get the COM Doppler bin, corrected for ephemeris drift and adjusted for orbital motion */
	kernel_doubles[4] = frame->dopcom_vig + frame->view[v].dopoff + orbit_dopoff;	// [4] = dopshift

	/*	Now calculate launch parameters for the kernel and launch it. 		 */
	maxThreadsPerBlock = 128;
	int npixels = (pos->xlim[1] - pos->xlim[0]) * (pos->ylim[1] - pos->ylim[0]);
	int Bx = floor((maxThreadsPerBlock - 1 + npixels) / maxThreadsPerBlock);
	dim3 BLK, THD;
	BLK.x = Bx; BLK.y = 1; BLK.y = 1;					// Grid-block dimension for the 1-D case
	THD.x = maxThreadsPerBlock; THD.y = 1; THD.z = 1;	// Thread block dimensions
	char *location = "pos2doppler_cuda.cu line 492";

	/*	Call the kernel 	*/
	pos2doppler_pixel_kernel<<<BLK,THD>>>(pos, frame, fit_overflow,	radar, radtype,
			kernel_ints, kernel_uchars, kernel_doubles, fitf);

	checkErrorAfterKernelLaunch(location);
	deviceSyncAfterKernelLaunch(location);

	/*	Copy the float array of fit values back to doubles in frame->fit[]	*/
	int i;
	for (i=0; i<=frame->ndop; i++)
		frame->fit[i] = (double)fitf[i];

	/*	Convert the model's floating-point Doppler limits from floating-point bin numbers to Hz, and
	 * 	widen the limits to account for nonzero POS pixel width. kernel_doubles[4] = dopshift                 */
	frame->doplim[0] = (frame->doplim[0] - kernel_doubles[4])*doppler->dop_per_bin - dop_extra;
	frame->doplim[1] = (frame->doplim[1] - kernel_doubles[4])*doppler->dop_per_bin + dop_extra;

	/*	Calculate the overflow contributions to chi squared:
	    o2 = obs^2 contribution, m2 = model^2 contribution.
        Also compute the summed cross section and the mean Doppler bin
    	for the overflow region, for use with the "delcorinit" action   */
	frame->overflow_o2 = 0.0;
	frame->overflow_m2 = 0.0;
	frame->overflow_xsec = 0.0;
	frame->overflow_dopmean = 0.0;
	sdev_sq = frame->sdev*frame->sdev;
	variance = sdev_sq;
	lookfact = (frame->nlooks > 0.0) ? 1.0/frame->nlooks : 0.0;
	any_overflow = kernel_ints[7]; //idop0[1];	// re-assign the any_overflow flag from the pointer array

	if (any_overflow) {
		j1 = MAX( frame->idoplim[0] + kernel_ints[6], 0);
		j2 = MIN( frame->idoplim[1] + kernel_ints[6], MAXOVERFLOW - 1);
		for (j=j1; j<=j2; j++) {
			if (fit_overflow[j] != 0.0) {
				if (par->speckle)
					variance = sdev_sq + lookfact*(double)(fit_overflow[j]*fit_overflow[j]);
				frame->overflow_o2 += 1.0;
				frame->overflow_m2 += (double)(fit_overflow[j]*fit_overflow[j])/variance;
				frame->overflow_xsec += (double)fit_overflow[j];
				frame->overflow_dopmean += (j - kernel_ints[6])*double(fit_overflow[j]);
			}
		}
		if (frame->overflow_xsec != 0.0)
			frame->overflow_dopmean /= frame->overflow_xsec;

		/*  Print a warning if the model extends even beyond the overflow spectrum  */
		if ( ((frame->idoplim[0] + kernel_ints[6]) < 0) ||
				((frame->idoplim[1] + kernel_ints[6]) >= MAXOVERFLOW) ) {
			badradar = 1;
			dopfactor = (MAX(frame->idoplim[1] + kernel_ints[6],MAXOVERFLOW)
					- MIN(frame->idoplim[0] + kernel_ints[6],0)) / (1.0*MAXOVERFLOW);
			frame->badradar_logfactor += log(dopfactor);
			if (par->warn_badradar) {
				printf("\nWARNING in pos2doppler.c for set %2d frame %2d:\n", set, frm);
				printf("        model Doppler spectrum extends too far beyond the data spectrum\n");
				printf("             data:  bins %2d to %2d\n", 1, kernel_ints[8]);
				printf("            model:  bins %2d to %2d\n",
						frame->idoplim[0], frame->idoplim[1]);
				fflush(stdout);
			}
		}
	}

	// start debug code
	//int count =

	// Free up the cudaMallocManaged-allocated structs & variables
	//cudaFree(frame);
	//	cudaFree(pos);
	//	cudaFree(fit_overflow);
	//	cudaFree(parc);
	//	cudaFree(idop0);

	return badradar;
}
