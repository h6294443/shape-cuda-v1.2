extern "C" {
#include "../shape/head.h"
}
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

/*  Loop through all POS pixels within the rectangular plane-of-sky region
 *  spanned by the model; for each such pixel which isn't blank sky, compute
 *  the cross-section contributions to pixels in the model delay-Doppler
 *  frame. Note that functions posclr and posvis flag blank-sky pixels by
 *  assigning "cose" = cos(scattering angle) = 0. Only compute contributions
 *  from POS pixels that project onto the right body, in case this is the
 *  "orbit" action (for which this routine is called twice, once for each of
 *  the two orbiting bodies). */

__global__ void pos2deldop_pixel_kernel(

		struct pos_t *pos,
		struct deldopfrm_t *frame,
		union radscat_t *radar,
		unsigned char *radtype,

		float *deldoplim,
		float **fit_overflow,
		float **fitf,
		float **map_fitf,
		float **map_posf,

		int *kernel_ints,
		double *kernel_doubles,
		unsigned char *kernel_uchars)
{
	/*	Find the facet index with thread and block indices.  First are global
	 * 	thread indices in x and y directions.  Then total offset.  Finally i
	 * 	and j which index i1 <= i <= i2 and j1 <= j <= j2	*/
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int x = pos->xlim[0] + (offset % (pos->xlim[1]-pos->xlim[0]));
	int y = pos->ylim[0] + (offset / (pos->xlim[1]-pos->xlim[0]));
	double tmp, arg_sample, arg_bl, arg_left, sinc2_mean, sinc2arg, amp,
	sumweights, fit_contribution, dop_contribution[MAXBINS],
	del_contribution[MAXBINS], codefactor;
	int c, f, i, j, k, m, in_bounds, idop, idop1, idop2, idel1, idel2,
		m_min, m_max, idop_min, idop_max, idel_min, idel_max, idel;
	float delPOS, dopPOS;

	if (offset<(pos->xlim[0]-pos->xlim[1])*(pos->ylim[0]-pos->ylim[1])){

		if (pos->cose[x][y] > 0.0 && pos->body[x][y] == kernel_ints[3]) {

			/* Get the (floating-point) delay and Doppler bin of the POS pixel
			 * center: delPOS and dopPOS. Also get the minimum and maximum
			 * (integer) delay and Doppler bins to which pixel contributes
			 * power: idel_min and idel_max, idop_min and idop_max. Strictly
			 * speaking, each POS pixel contributes power to *all* Doppler
			 * columns, but here we're zeroing out the sinc^2 response function
			 * beyond the nearest sinc2width columns. Actually, if nsinc2 > 1,
			 * we'll distribute power to *at least* sinc2width Doppler bins:
			 * For pixels which span multiple bins we'll err on the side of
			 * computing more contributions rather than fewer.       */

			delPOS = pos->z[x][y]*kernel_doubles[9] + kernel_doubles[8];
			idel_min = (int) floor(delPOS - kernel_doubles[11]) + 1;
			idel_max = (int)  ceil(delPOS + kernel_doubles[11]) - 1;
			dopPOS = kernel_doubles[0]*(x - kernel_doubles[2]) +
			  		 kernel_doubles[1]*(y - kernel_doubles[3]) +
			  		 kernel_doubles[4];
			idop_min = (int) floor(dopPOS - kernel_doubles[5] + 1 -
					 kernel_ints[0]/2.0);
			idop_max = (int) floor(dopPOS + kernel_doubles[5] +
					kernel_ints[0]/2.0);

			/* For the short code, sensitivity drops as we move away from DC.
			 * (This variation is slow, so we can just evaluate the response
			 * at the center of the POS pixel.) Note that the SINC2 macro
			 * multiplies its argument by pi.        */
			codefactor = (kernel_ints[7] == SHORT) ? SINC2( (dopPOS -
					kernel_doubles[6])/kernel_ints[6] ) : 1.0;

			/* Update rectangular delay-Doppler region (row/column numbers)
			 * with nonzero power according to model. Old code commented out.
			 * New code uses atomicMin and atomicMax  */
			/*if (idel_min < frame->idellim[0])	frame->idellim[0] = idel_min;
			if (idel_max > frame->idellim[1])	frame->idellim[1] = idel_max;
			if (idop_min < frame->idoplim[0])	frame->idoplim[0] = idop_min;
			if (idop_max > frame->idoplim[1])	frame->idoplim[1] = idop_max;*/
			atomicMin(&frame->idellim[0], idel_min);
			atomicMax(&frame->idellim[1], idel_max);
			atomicMin(&frame->idoplim[0], idop_min);
			atomicMax(&frame->idoplim[1], idop_max);

			/* Update the model's floating-point delay-Doppler limits, as de-
			 * termined prior to convolution with delay and Doppler response
			 * functions. At this point in the code, dellim and doplim are
			 * pairs of floating-point row and column numbers which apply to
			 * POS pixel centers; when the loop over POS pixels is finished we
			 * will convert them to usec and Hz, and will widen the Doppler
			 * limits to account for nonzero POS pixel width.
			 * Old standard code commented out.  New CUDA code uses a float
			 * copy of frame->dellim[] and frame->doplim[] and atomics.   */
			/*frame->dellim[0] = MIN( frame->dellim[0], delPOS);
			frame->dellim[1] = MAX( frame->dellim[1], delPOS);
			frame->doplim[0] = MIN( frame->doplim[0], dopPOS);
			frame->doplim[1] = MAX( frame->doplim[1], dopPOS);*/
			atomicMinf(&deldoplim[0], (float)delPOS);
			atomicMaxf(&deldoplim[1], (float)delPOS);
			atomicMinf(&deldoplim[2], (float)dopPOS);
			atomicMaxf(&deldoplim[3], (float)dopPOS);

			/* Check whether or not all delay-Doppler pixels which will re-
			 * ceive power from this POS pixel fall within the data frame;
			 * if not, initialize the "overflow" image if necessary.       */
			if ( (idel_min >= 1) && (idel_max <= kernel_ints[12]) &&
					(idop_min >= 1) && (idop_max <= kernel_ints[13]))
				in_bounds = 1;
			else {
				in_bounds = 0;
				if (!kernel_ints[5]) {
					kernel_ints[5] = 1;
					for (i=0; i<MAXOVERFLOW; i++)
						for (j=0; j<MAXOVERFLOW; j++) // To-Do:  this may need
							fit_overflow[i][j] = 0.0; // atomic attention in the future
					/*  Center the COM in the overflow image: pixel[idel][idop]
					 *  in the fit frame corresponds to
					 *  pixel[idel+idel0][idop+idop0] in the fit_overflow frame.  */
					kernel_ints[14] = MAXOVERFLOW/2 -
							(int)floor(kernel_doubles[8]+0.5);
					kernel_ints[15] = MAXOVERFLOW/2 -
							(int)floor(kernel_doubles[4]+0.5);
				}
			}

			/* Loop through all delay bins to which this POS pixel contributes
			 * power (cross section), and compute the delay response function
			 * for each bin                 */
			for (idel=idel_min; idel<=idel_max; idel++) {
				if (kernel_ints[7] != LONG_ORIG) {

					/* Get the delay response function for image row idel: sum
					 * the triangle-function contributions from each sample
					 * per baud, then divide the sum by spb and square.
					 * The triangle function for sample m  (0 <= m <= spb-1)
					 * has unit height and a half-width of spb/stride image
					 * rows, and is centered [-const2 + m/stride] rows later
					 * than the row center (idel).
					 * In the code block below, the arguments to macros TRI
					 * and TRI2 have been divided by half-width spb/stride,
					 * since those two macros are defined to give nonzero va-
					 * lues for arguments between -1 and 1.  Each argument,
					 * then, is just
					 * 	(delPOS - [triangle-function center]) / half-width.
					 * Do the two most common cases (spb = 1 or 2) without
					 * loops in order to gain a bit of speed.  For the other
					 * cases, set m_min and m_max so as not to waste time
					 * computing contributions that are zero.             */

					/* NOTE: The atomics below are commented out because I don't think
					 * they are needed.  del_contribution[] is local to each thread.  */
					switch (kernel_ints[8]) {
					case 1:
						del_contribution[idel-idel_min] = TRI2( delPOS - idel );
						break;
					case 2:
						arg_sample = (delPOS - (idel - kernel_doubles[12])) /
							kernel_ints[10];
						del_contribution[idel-idel_min] = TRI(arg_sample)
										+ TRI( arg_sample - 0.5 );
						del_contribution[idel-idel_min] *=
								del_contribution[idel-idel_min]/4;
						break;
					default:
						del_contribution[idel-idel_min] = 0.0;
						m_min = MAX( (int) floor((delPOS - idel -
								kernel_doubles[12])*kernel_ints[11]) , 0 );
						m_max = MIN( (int) ceil((delPOS - idel +
								kernel_doubles[11])*kernel_ints[11]) ,
								kernel_ints[8] ) - 1;
						arg_sample = (delPOS - (idel - kernel_doubles[12])) /
								kernel_ints[10] - m_min*kernel_doubles[10];
						for (m=m_min; m<=m_max; m++) {
							del_contribution[idel-idel_min] += TRI( arg_sample );
							arg_sample -= kernel_doubles[10];
						}
						del_contribution[idel-idel_min] *=
								del_contribution[idel-idel_min]/kernel_ints[9];
						break;
					}
				} else {

					/* Long code with original (Harmon) reduction method: data
					 * for each sample per baud are reduced separately,as if
					 * datataking were just one sample per baud; then the im-
					 * age rows for spb/stride samples are interleaved.  */
					del_contribution[idel-idel_min] = TRI2( (delPOS - idel)/
							kernel_ints[10] );
				}
			}

			/* Now include sinc^2 factor for Doppler mismatching: Take the mean
			 * of nsinc2^2 points interior to POS pixel. Do two most common ca-
			 * ses (nsinc2=1 or 2) without loops. Note the SINC2 macro
			 * multiplies its argument by pi*/

			for (idop=idop_min; idop<=idop_max; idop++) {
				switch (kernel_ints[1]) {
				case 1:
					sinc2_mean = SINC2( dopPOS - idop );
					break;
				case 2:
					arg_bl = dopPOS + kernel_doubles[7] - idop;   /* bl = bottom left */
					sinc2_mean = ( SINC2( arg_bl ) +
							SINC2( arg_bl+kernel_doubles[13] ) +
							SINC2( arg_bl+kernel_doubles[14] ) +
							SINC2( arg_bl+kernel_doubles[13] +
									kernel_doubles[14] ) ) / 4;
					break;
				default:
					arg_left = dopPOS + kernel_doubles[7] - idop;
					sinc2_mean = 0.0;
					for (i=0; i<kernel_ints[1]; i++) {
						sinc2arg = arg_left;
						for (j=0; j<kernel_ints[1]; j++) {
							sinc2_mean += SINC2( sinc2arg );
							sinc2arg += kernel_doubles[13];
						}
						arg_left += kernel_doubles[14];
					}
					sinc2_mean /= kernel_ints[2];
					break;
				}
				k = MIN( idop - idop_min, MAXBINS);
				dop_contribution[k] = sinc2_mean;
			}

			/* Compute the sum of delay-Doppler weighting factors  */
			sumweights = 0.0;
			for (idel=idel_min; idel<=idel_max; idel++)
				for (idop=idop_min; idop<=idop_max; idop++) {
					k = MIN( idop - idop_min, MAXBINS);
					sumweights += del_contribution[idel-idel_min]*
							dop_contribution[k];
				}

			/* The radar cross section within this plane-of-sky pixel is
			 *
			 *    [differential radar scattering law]*[POS pixel area in km^2]
			 *
			 * The differential radar scattering law (function radlaw
			 *  = d[cross section]/d[area] ) includes a sec(theta) factor to
			 * account for the fact that the POS pixel area is projected area
			 * rather than physical area on the target surface.      */

			amp = radlaw_cuda(radar, radtype, kernel_ints[4], pos->cose[x][y],
					pos->comp[x][y], pos->f[x][y])
			    		  * pos->km_per_pixel * pos->km_per_pixel
			    		  * codefactor / sumweights;

			/* Add this POS pixel's power contributions to the model delay-
			 * Doppler frame if NONE of those contributions fall outside the
			 * frame limits.                                   */

			if (in_bounds) {

				/* Add the cross-section contributions to the model frame  */
				for (idel=idel_min; idel<=idel_max; idel++)
					for (idop=idop_min; idop<=idop_max; idop++) {
						k = MIN( idop - idop_min, MAXBINS);
						fit_contribution = amp * del_contribution[idel-idel_min]
						                                          * dop_contribution[k];
						atomicAdd(&fitf[idel][idop], fit_contribution);
						//frame->fit[idel][idop] += fit_contribution;

						if (kernel_uchars[0] == MAP) {
							if (kernel_uchars[1] == MAPMODE_DELDOP) {
								if (map_fitf[idel][idop] > 0.0) {
									//if (frame->map_fit[idel][idop] > 0.0) { //To-Do: trouble maybe
									atomicAdd(&map_posf[x][y], fit_contribution);
									//frame->map_pos[x][y] += fit_contribution;
									c = pos->comp[x][y];
									f = pos->f[x][y];
									/* The following is not allocated anywhere.  To-Do.	*/
									frame->map_facet_power[c][f] += fit_contribution;
									if (kernel_uchars[2])
										printf("# POS (%3d, %3d) comp %d facet %4d contributes %e to d-D (%3d, %3d)\n",
												x+pos->n, y+pos->n, c, f, fit_contribution, idel-1, idop-1);
								}
							} else if (kernel_uchars[1] == MAPMODE_POS) {
								if (map_posf[x][y] > 0.0) {
									//if (frame->map_pos[x][y] > 0.0) {
									atomicAdd(&map_fitf[idel][idop], fit_contribution);
									//frame->map_fit[idel][idop] += fit_contribution;
									c = pos->comp[x][y];
									f = pos->f[x][y];
									/* The following is not allocated anywhere.  To-Do.	*/
									frame->map_facet_power[c][f] += fit_contribution;
									if (kernel_uchars[2])
										printf("# POS (%3d, %3d) comp %d facet %4d contributes %e to d-D (%3d, %3d)\n",
												x+pos->n, y+pos->n, c, f, fit_contribution, idel-1, idop-1);
								}
							} else {
								if (map_posf[x][y] > 0.0) {
									//if (frame->map_pos[x][y] > 0.0) {
									atomicAdd(&map_fitf[idel][idop], fit_contribution);
									//frame->map_fit[idel][idop] += fit_contribution;
									if (kernel_uchars[2]) {
										c = pos->comp[x][y];
										f = pos->f[x][y];
										printf("# POS (%3d, %3d) comp %d facet %4d contributes %e to d-D (%3d, %3d)\n",
												x+pos->n, y+pos->n, c, f, fit_contribution, idel-1, idop-1);
									}
								}
							}
						}
					}

			} else {

				/*  Add the cross-section contributions to the "overflow" image  */
				if (kernel_uchars[0] == MAP && kernel_uchars[1] != MAPMODE_DELDOP)
					if (frame->map_pos[x][y] > 0.0)
						kernel_uchars[3] = 1; // To-Do: See if this is needed somewhere down the line

				idel1 = MAX( idel_min, -kernel_ints[14]);
				idel2 = MIN( idel_max, -kernel_ints[14] + MAXOVERFLOW - 1);
				idop1 = MAX( idop_min, -kernel_ints[15]);
				idop2 = MIN( idop_max, -kernel_ints[15] + MAXOVERFLOW - 1);
				for (idel=idel1; idel<=idel2; idel++)
					for (idop=idop1; idop<=idop2; idop++) {
						k = MIN( idop - idop_min, MAXBINS);
						fit_contribution = amp * del_contribution[idel-idel_min]
						                                          * dop_contribution[k];
						//fit_overflow[idel+idel0][idop+idop0] += fit_contribution;
						atomicAdd(&fit_overflow[idel+kernel_ints[14]][idop+kernel_ints[15]], (float)fit_contribution);
						if (kernel_uchars[0] == MAP && kernel_uchars[1] == MAPMODE_DELDOP)
							if (idel >= kernel_ints[16] && idel <= kernel_ints[17] &&
									idop >= kernel_ints[18] && idop <= kernel_ints[19]    ) {
								atomicAdd(&map_posf[x][y], fit_contribution);
								//frame->map_pos[x][y] += fit_contribution;
								c = pos->comp[x][y];
								f = pos->f[x][y];
								/* The following is not allocated anywhere.  To-Do.	*/
								frame->map_facet_power[c][f] += fit_contribution;
								if (kernel_uchars[2])
									printf("# POS (%3d, %3d) comp %d facet %4d contributes %e to d-D (%3d, %3d)\n",
											x+pos->n, y+pos->n, c, f, fit_contribution, idel-1, idop-1);
							}
					}
			}
		}  /* if cos(scattering angle) > 0 and POS pixel projects onto the right body */
	}  /* go to the next POS image pixel */
}



int pos2deldop_cuda(struct par_t *par, struct photo_t *photo, double orbit_xoff,
		double orbit_yoff, double orbit_dopoff, struct deldop_t *deldop, int body,
		int set, int frm, int v)
{
	int codemethod,	dopfftlen, spb, stride, spb_sq, nsinc2_sq, i,i1, i2, j, j1, j2,
	any_overflow, idel0, idop0, spb_over_stride, ndel, ndop,	badradar;
	int *kernel_ints;

	double dopDC_vig, delfact, dopfact, w[3], ax, ay, delshift, dopshift,
	const1, const2, xincr, yincr, dopdiff_bl, dopdiff_max, one_over_spb,
	lookfact, sdev_sq, variance, dop_extra, delfactor, dopfactor;
	double *kernel_doubles;
	float *deldoplim, **fitf, **map_fitf, **map_posf, **map_facet_powerf;
	static float **fit_overflow;// Using a float instead of double for atomics
	struct deldopfrm_t *frame;
	struct pos_t *pos;
	union radscat_t *radar;
	unsigned char *radtype, *kernel_uchars;

	/* Allocate pointers in CUDA Unified Memory	*/
	cudaCalloc1((void**)&frame, 			  sizeof(struct deldopfrm_t),		1);
	cudaCalloc1((void**)&pos, 			  sizeof(struct pos_t),				1);
	cudaCalloc1((void**)&kernel_ints, 	  sizeof(int), 			   		   10);
	cudaCalloc1((void**)&kernel_doubles,   sizeof(double), 		   		   16);
	cudaCalloc1((void**)&deldoplim,	 	  sizeof(float), 					4);
	cudaCalloc1((void**)&kernel_ints, 	  sizeof(int), 					   20);
	cudaCalloc1((void**)&kernel_doubles,   sizeof(double), 				   15);
	cudaCalloc1((void**)&kernel_uchars, 	  sizeof(unsigned char), 			4);
	cudaCalloc1((void**)&map_posf,		  sizeof(float*), 		   2*pos->n+1);
	cudaCalloc1((void**)&fit_overflow, 	  sizeof(float*), 	  		  MAXBINS);
	cudaCalloc1((void**)&map_fitf,		  sizeof(float*),deldop->frame[frm].ndel);
	cudaCalloc1((void**)&fitf, 			  sizeof(float*),deldop->frame[frm].ndel);
	cudaCalloc1((void**)&radar, 			  sizeof(union radscat_t),photo->nradlaws);
	cudaCalloc1((void**)&radtype,		  sizeof(unsigned char),  photo->nradlaws);
	/* Note that allocation for map_facet_powerf is skipped for now until I
	 * figure out where that is allocated in standard code.	 */

	/* Fix pointer indexing in map_posf[x][y] outer loop.  Then loop through
	 * outer loop and allocate inner loops in CUDA Unified Memory.  Lastly,
	 * fix the pointer indexing in the inner loop addressing.	*/
	map_posf -= -pos->n;
	for (i=0; i<MAXBINS; i++)
		cudaCalloc1((void**)&fit_overflow[i], sizeof(float), MAXBINS);
	for (i=0; i<2*pos->n+1; i++){
		cudaCalloc1((void**)&map_posf[i], sizeof(float), 2*pos->n+1);
		map_posf[i] -= -pos->n;
	}
	for (i=0; i<deldop->frame[frm].ndel; i++){
		cudaCalloc1((void**)&fitf[i],    sizeof(float),deldop->frame[frm].ndop);
		cudaCalloc1((void**)&map_fitf[i],sizeof(float),deldop->frame[frm].ndop);
	}

	/* Initialize variables to avoid compilation warnings  */
	idel0 = idop0 = 0;
	any_overflow = 0;
	frame = &deldop->frame[frm];
	pos = &frame->pos;
	radtype = photo->radtype;
	radar = photo->radar;
	ndel = frame->ndel;
	ndop = frame->ndop;
	badradar = 0;

	frame->badradar_logfactor = 0.0;
	frame->idellim[0] = ndel + 999999;
	frame->idellim[1] = -999999;
	frame->idoplim[0] = ndop + 999999;
	frame->idoplim[1] = -999999;
	frame->dellim[0] = deldoplim[0] =  HUGENUMBER;
	frame->dellim[1] = deldoplim[1] = -HUGENUMBER;
	frame->doplim[0] = deldoplim[2] =  HUGENUMBER;
	frame->doplim[1] = deldoplim[3] = -HUGENUMBER;

	/* Get parameters related to data sampling and data reduction; then
	 * compute two more (both in units of delay bins = image rows):
	 *  	const1: half of the base width of the delay response function
	 *  	const2: half the delay difference between the first and last
	 *  	image rows within each baud  */

	codemethod = deldop->codemethod;
	spb = deldop->spb;
	stride = deldop->stride;
	dopfftlen = deldop->dopfftlen;
	dopDC_vig = frame->dopDC_vig;
	spb_over_stride = spb/stride;
	one_over_spb = 1.0/spb;
	spb_sq = pow(spb,2);
	if (codemethod != LONG_ORIG) {
		const1 = (3*spb - 1)/(2.0*stride);
		const2 = (spb - 1)/(2.0*stride);
	} else {
		const1 = (double) spb_over_stride;
		const2 = 0.0;  /* not used with this code + reduction method */
	}

	/*  Converts from km towards radar to delay bins  */
	delfact = -KM2US/deldop->del_per_pixel;

	/*  Get w, the apparent spin vector in observer coordinates  */
	cotrans( w, frame->view[v].oe, frame->view[v].spin, 1);

	/*  Compute the Doppler bin increment per plane-of-sky pixel westward (ax)
	    and northward (ay); these values are scaled by the "dopscale" parameter
	    for this dataset.  Then compute km2Hz, the Doppler increment (Hz) per
	    km perpendicular to the projected spin axis in the plane of the sky.     */

	dopfact = deldop->dopscale.val * KM2HZFACT * pos->km_per_pixel
			* deldop->Ftx / deldop->dop_per_pixel;
	ax = -w[1]*dopfact;
	ay =  w[0]*dopfact;
	frame->view[v].km2Hz = sqrt(ax*ax + ay*ay) * deldop->dop_per_pixel
			/ pos->km_per_pixel;

	/*  Compute the absolute value of the difference between the maximum (or minimum)
	 *  Doppler on any given POS pixel's edge and the Doppler at its center             */

	if (w[0] != 0.0 || w[1] != 0.0)
		dop_extra = frame->view[v].km2Hz * 0.5 * pos->km_per_pixel
		* sqrt(w[0]*w[0] + w[1]*w[1]) / MAX( fabs(w[0]), fabs(w[1]));
	else
		dop_extra = 0.0;

	/*  We may be evaluating the sinc^2 Doppler response function at more than one point
	 *  per POS pixel.  xincr and yincr are the Doppler bin increments between adjacent
	 *  evaluation points in the x and y directions.  dopdiff_bl is the Doppler bin difference
	 *  between the bottom-leftmost (southeasternmost) evaluation point and the pixel center.
	 *  dopdiff_max is the maximum positive Doppler bin difference between any evaluation point
	 *  and the pixel center.                                                     */

	nsinc2_sq = pow(par->nsinc2,2);
	xincr = ax / par->nsinc2;
	yincr = ay / par->nsinc2;
	dopdiff_bl = -(par->nsinc2 - 1)*(xincr + yincr)/2;
	dopdiff_max = (par->nsinc2 - 1)*(fabs(xincr) + fabs(yincr))/2;
	if (2*dopdiff_max + par->sinc2width + 1 > MAXBINS) {
		badradar = 1;
		frame->badradar_logfactor += log((2*dopdiff_max + par->sinc2width + 1) / MAXBINS);
		if (par->warn_badradar) {
			printf("\nWARNING in pos2deldop.c for set %2d frame %2d:\n", set, frm);
			printf("        sinc^2 function evaluated at %d Doppler bins, must be <= %d\n",
					(int) ceil(2*dopdiff_max + par->sinc2width + 1), MAXBINS);
			fflush(stdout);
		}
	}

	/* Get the COM delay and Doppler bins, corrected for ephemeris drift and adjusted for
	 * orbital motion; the delay adjustment for orbital motion is done implicitly (the
	 * posvis routine has already adjusted the "z" values for all POS pixels), whereas the
	 * Doppler adjustment for orbital motion must be done here explicitly.                  */
	delshift = frame->delcom_vig + frame->view[v].deloff;
	dopshift = frame->dopcom_vig + frame->view[v].dopoff + orbit_dopoff;

	/* Load kernel parameter arrays	*/
	kernel_ints[0] = par->sinc2width;
	kernel_ints[1] = par->nsinc2;
	kernel_ints[2] = nsinc2_sq;
	kernel_ints[3] = body;
	kernel_ints[4] = deldop->iradlaw;
	kernel_ints[5] = any_overflow;
	kernel_ints[6] = dopfftlen;
	kernel_ints[7] = codemethod;
	kernel_ints[8] = spb;
	kernel_ints[9] = spb_sq;
	kernel_ints[10] = spb_over_stride;
	kernel_ints[11] = stride;
	kernel_ints[12] = ndel;
	kernel_ints[13] = ndop;
	kernel_ints[14] = idel0;
	kernel_ints[15] = idop0;
	kernel_ints[16] = par->map_dellim[0];
	kernel_ints[17] = par->map_dellim[1];
	kernel_ints[18] = par->map_doplim[0];
	kernel_ints[19] = par->map_doplim[1];

	kernel_uchars[0] = par->action;
	kernel_uchars[1] = par->map_mode;
	kernel_uchars[2] = par->map_verbose;
	kernel_uchars[3] = par->map_overflow;

	kernel_doubles[0] = ax;
	kernel_doubles[1] = ay;
	kernel_doubles[2] = orbit_xoff;
	kernel_doubles[3] = orbit_yoff;
	kernel_doubles[4] = dopshift;
	kernel_doubles[5] = dopdiff_max;
	kernel_doubles[6] = dopDC_vig;
	kernel_doubles[7] = dopdiff_bl;
	kernel_doubles[8] = delshift;
	kernel_doubles[9] = delfact;
	kernel_doubles[10] = one_over_spb;
	kernel_doubles[11] = const1;
	kernel_doubles[12] = const2;
	kernel_doubles[13] = xincr;
	kernel_doubles[14] = yincr;

	/* Set kernel launch parameters	*/
	maxThreadsPerBlock = 128;
	int npixels = (pos->xlim[1]-pos->xlim[0]) * (pos->ylim[1]-pos->ylim[0]);
	int Bx = floor((maxThreadsPerBlock - 1 + npixels) / maxThreadsPerBlock);
	dim3 BLK, THD;
	BLK.x = Bx; BLK.y = 1; BLK.y = 1;					// Grid-block dimension for the 1-D case
	THD.x = maxThreadsPerBlock; THD.y = 1; THD.z = 1;	// Thread block dimensions
	char *location = "pos2deldop_cuda.cu line 597";

	/* Launch kernel	 */
	pos2deldop_pixel_kernel<<<BLK,THD>>>(pos, frame, radar, radtype,
			deldoplim, fit_overflow, fitf, map_fitf, map_posf,
			kernel_ints, kernel_doubles, kernel_uchars);

	/* Post-launch checks	 */
	checkErrorAfterKernelLaunch(location);
	deviceSyncAfterKernelLaunch(location);

	/* Convert and copy the float deldop limits back to frame	*/
	frame->dellim[0] = (double)deldoplim[0];
	frame->dellim[1] = (double)deldoplim[1];
	frame->doplim[0] = (double)deldoplim[2];
	frame->doplim[1] = (double)deldoplim[3];

	/*	Debug code 	mixed with regular code */
	//double summapfit = 0, sumframefit = 0, summappos = 0;	//debug

	/* Copy the float arrays of float arrays back to doubles	*/
	for (i=frame->idellim[0]; i<frame->idellim[1]; i++)
		for (j=frame->idoplim[0]; j<frame->idoplim[1]; j++){
			frame->fit[i][j] = fitf[i][j];
			//sumframefit += frame->fit[i][j];			//debug
			if (par->action == MAP){
				frame->map_fit[i][j] = map_fitf[i][j];
				//summapfit += frame->map_fit[i][j];	//debug
			}
		}
	if (par->action == MAP)
		for (i=-pos->n; i<=pos->n; i++)
			for (j=-pos->n; j<=pos->n; j++){
				frame->map_pos[i][j] = map_posf[i][j];
				//summappos += map_posf[i][j];			//debug
			}

//	if(sumframefit>0.0)
//		printf("\nsum of frame->fit[][]: %f\n", sumframefit);
//	if(summapfit>0.0)
//		printf("\sum of frame->map_fit[][]: %f\n", summapfit);
//	if(summappos>0.0)
//		printf("\sum of frame->map_pos[][]: %f\n", summappos);

	/*	End debug	*/


	/* Convert the model's floating-point delay-Doppler limits from floating-
	 * point row and column numbers to usec and Hz, and widen the Doppler limits
	 * to account for nonzero POS pixel width  */
	frame->dellim[0] = (frame->dellim[0] - delshift)*deldop->del_per_pixel;
	frame->dellim[1] = (frame->dellim[1] - delshift)*deldop->del_per_pixel;
	frame->doplim[0] = (frame->doplim[0] - dopshift)*deldop->dop_per_pixel
			- dop_extra;
	frame->doplim[1] = (frame->doplim[1] - dopshift)*deldop->dop_per_pixel
			+ dop_extra;

	/* Calculate the overflow contributions to chi squared:
	 * 	 o2 = obs^2 contribution, m2 = model^2 contribution.
	 *
	 * Also compute the summed cross section and the mean delay and Doppler
	 * bins for the overflow region, for use with the "delcorinit" action    */

	frame->overflow_o2 = 0.0;
	frame->overflow_m2 = 0.0;
	frame->overflow_xsec = 0.0;
	frame->overflow_delmean = 0.0;
	frame->overflow_dopmean = 0.0;
	sdev_sq = frame->sdev*frame->sdev;
	variance = sdev_sq;
	lookfact = (frame->nlooks > 0.0) ? 1.0/frame->nlooks : 0.0;
	if (kernel_ints[5]) {
		i1 = MAX( frame->idellim[0] + idel0, 0);
		i2 = MIN( frame->idellim[1] + idel0, MAXOVERFLOW - 1);
		j1 = MAX( frame->idoplim[0] + idop0, 0);
		j2 = MIN( frame->idoplim[1] + idop0, MAXOVERFLOW - 1);
		for (i=i1; i<=i2; i++)
			for (j=j1; j<=j2; j++) {
				if (fit_overflow[i][j] != 0.0) {
					if (par->speckle)
						variance = sdev_sq + lookfact*fit_overflow[i][j]*fit_overflow[i][j];
					frame->overflow_o2 += 1.0;
					frame->overflow_m2 += fit_overflow[i][j]*fit_overflow[i][j]/variance;
					frame->overflow_xsec += fit_overflow[i][j];
					frame->overflow_delmean += (i - idel0)*fit_overflow[i][j];
					frame->overflow_dopmean += (j - idop0)*fit_overflow[i][j];
				}
			}
		if (frame->overflow_xsec != 0.0) {
			frame->overflow_delmean /= frame->overflow_xsec;
			frame->overflow_dopmean /= frame->overflow_xsec;
		}

		/* Print a warning if the model extends even beyond the overflow image  */

		if ( ((frame->idellim[0] + idel0) < 0)            ||
				((frame->idellim[1] + idel0) >= MAXOVERFLOW) ||
				((frame->idoplim[0] + idop0) < 0)            ||
				((frame->idoplim[1] + idop0) >= MAXOVERFLOW)    ) {
			badradar = 1;
			delfactor = (MAX( frame->idellim[1] + idel0, MAXOVERFLOW)
					- MIN( frame->idellim[0] + idel0, 0)         )
	                						  / (1.0*MAXOVERFLOW);
			dopfactor = (MAX( frame->idoplim[1] + idop0, MAXOVERFLOW)
					- MIN( frame->idoplim[0] + idop0, 0)         )
	                						  / (1.0*MAXOVERFLOW);
			frame->badradar_logfactor += log(delfactor*dopfactor);
			if (par->warn_badradar) {
				printf("\nWARNING in pos2deldop.c for set %2d frame %2d:\n", set, frm);
				printf("        model delay-Doppler image extends too far beyond the data image\n");
				printf("             data:  rows %2d to %2d , cols %2d to %2d\n", 1, ndel, 1, ndop);
				printf("            model:  rows %2d to %2d , cols %2d to %2d\n",
						frame->idellim[0], frame->idellim[1],
						frame->idoplim[0], frame->idoplim[1]);
				fflush(stdout);
			}
		}
	}



	return badradar;
}
