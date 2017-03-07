/*****************************************************************************************
 posvis.c

 Fill in the portion of a plane-of-sky image due to a particular model component: Assign
 each relevant POS pixel a z-value in observer coordinates (distance from the origin
 towards Earth) and a value of cos(scattering angle).

 Return 1 if any portion of this component lies outside the specified POS window,
 0 otherwise.

 If the "src" argument is true, the "observer" is the Sun rather than Earth, and
 "plane-of-sky" becomes "projection as viewed from the Sun."

 Modified 2014 February 20 by CM:
 Allow facets that partly project outside the POS frame to contribute to the POS frame
 (thus avoiding see-through "holes" in the model at the edge of a POS image)

 Modified 2010 May 18 by CM:
 Bug fix: When checking if a POS pixel hasn't already been assigned
 values during a previous call to posvis for a different component,
 check for fac[i][j] < 0 rather than cosa[i][j] == 0.0, since for
 bistatic situations the latter condition will also be true for
 pixels centered on Earth-facing facets that don't face the Sun

 Modified 2009 July 2 by CM:
 Eliminate the check that facets are "active": this term is now being
 interpreted to mean "not lying interior to the model," so the
 check is unnecessary and the determination of active vs. inactive
 status is inaccurate for half-exposed facets at the intersections
 between model components

 Modified 2009 April 3 by CM:
 Compute the "posbnd_logfactor" parameter: if the model extends beyond
 the POS frame, posbnd_logfactor is set to the logarithm of the
 ratio of the area that would have been required to "contain" the
 entire model divided by the area of the actual POS frame
 Work with floating-point pixel numbers (imin_dbl, etc.), at least
 initially, in case the sky rendering for a model with illegal
 parameters would involve huge pixel numbers that exceed the
 limits for valid integers

 Modified 2007 August 4 by CM:
 Add "orbit_offset" and "body" parameters and remove "facet" parameter
 Add body, bodyill, comp, and compill matrices for POS frames

 Modified 2006 June 21 by CM:
 For POS renderings, change res to km_per_pixel

 Modified 2005 September 19 by CM:
 Allow for roundoff error when determining which POS pixels project
 onto each model facet

 Modified 2005 June 27 by CM:
 Renamed "round" function to "iround" to avoid conflicts

 Modified 2005 June 22 by CM:
 Slightly modified some comments

 Modified 2005 January 25 by CM:
 Take care of unused and uninitialized variables

 Modified 2004 December 19 by CM:
 Added more comments
 Put update of rectangular POS area into "POSrect" routine and applied it
 even to facets which lie outside the POS frame

 Modified 2004 Feb 11 by CM:
 Added comments

 Modified 2003 May 5 by CM:
 Removed redundant coordinate transformation of the unit normal n
 for the no-pvs_smoothing case
 *****************************************************************************************/
extern "C" {
#include "head.h"
#include <limits.h>
}
__device__ int pvsoutbnd=0, pvs_nf, pvs_smooth, pvs_n;
__device__ struct vertices_t *_verts;

/* Note that the following two custom atomic functions are declared in each
 * file they are needed .  As static __device__ functions, this is the only
 * way to handle them. */
__device__ static float atomicMinf(float* address, float val) {
	int* address_as_i = (int*) address;
	int old = *address_as_i, assumed;
	do {
		assumed = old;
		old = ::atomicCAS(address_as_i, assumed,
				__float_as_int(::fminf(val, __int_as_float(assumed))));
	} while (assumed != old);
	return __int_as_float(old);
}
__device__ static float atomicMaxf(float* address, float val) {
	int* address_as_i = (int*) address;
	int old = *address_as_i, assumed;
	do {
		assumed = old;
		old = ::atomicCAS(address_as_i, assumed,
				__float_as_int(::fmaxf(val, __int_as_float(assumed))));
	} while (assumed != old);
	return __int_as_float(old);
}
__device__ void dev_POSrect3(struct pos_t *pos, int src, float imin_dbl, float
		imax_dbl, float jmin_dbl, float jmax_dbl, float4 *ijminmax_overall, int frm)	{
	int n, imin, imax, jmin, jmax;
	n = pvs_n;

	/* Update the POS region that contains the target without
	 * regard to whether or not it extends beyond the POS frame */
	atomicMinf(&ijminmax_overall[frm].w, imin_dbl);
	atomicMaxf(&ijminmax_overall[frm].x, imax_dbl);
	atomicMinf(&ijminmax_overall[frm].y, jmin_dbl);
	atomicMaxf(&ijminmax_overall[frm].z, jmax_dbl);

	/*  Update the subset of the POS frame that contains the target  */
	imin = (imin_dbl < INT_MIN) ? INT_MIN : (int) imin_dbl;
	imax = (imax_dbl > INT_MAX) ? INT_MAX : (int) imax_dbl;
	jmin = (jmin_dbl < INT_MIN) ? INT_MIN : (int) jmin_dbl;
	jmax = (jmax_dbl > INT_MAX) ? INT_MAX : (int) jmax_dbl;

	/* Make sure it's smaller than n */
	imin = MAX(imin,-n);
	imax = MIN(imax, n);
	jmin = MAX(jmin,-n);
	jmax = MIN(jmax, n);

	n = pos->n;
	if (src) {
		atomicMin(&pos->xlim2[0], imin);
		atomicMax(&pos->xlim2[1], imax);
		atomicMin(&pos->ylim2[0], jmin);
		atomicMax(&pos->ylim2[1], jmax);
	} else {
		atomicMin(&pos->xlim[0], imin);
		atomicMax(&pos->xlim[1], imax);
		atomicMin(&pos->ylim[0], jmin);
		atomicMax(&pos->ylim[1], jmax);
	}
}
__global__ void posvis_init_af_krnl(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, int set, int nframes, struct pos_t **pos,
		float4 *ijminmax_overall, int n) {
	/* nframes-threaded kernel to init variables for all frames in set */
	int frame = threadIdx.x;

	if (frame < nframes) {
		switch (ddat->set[set].type) {
		case DELAY:
			pos[frame] = &ddat->set[set].desc.deldop.frame[frame].pos;
			break;
		case DOPPLER:
			pos[frame] = &ddat->set[set].desc.doppler.frame[frame].pos;
			break;
		case POS:
			pos[frame] = &ddat->set[set].desc.poset.frame[frame].pos;
			break;
		case LGHTCRV:
			pos[frame] = &ddat->set[set].desc.lghtcrv.rend[frame].pos; //frame = i
			break;
		}
		/*  Initialize variables  */
		pos[frame]->posbnd_logfactor = 0.0;
		ijminmax_overall[frame].w = ijminmax_overall[frame].y = HUGENUMBER;
		ijminmax_overall[frame].x = ijminmax_overall[frame].z = -HUGENUMBER;

		_verts = &dmod->shape.comp[0].real;
		pvs_nf = _verts->nf;
		pvs_smooth = dpar->pos_smooth;
		pvs_n = pos[frame]->n;
	}
}

__global__ void posvis_facet_af_krnl(struct pos_t **pos, int src, int body,
		int comp, float3 orbit_offs, int total_size, int frame_size,
		int nframes, float4 *ijminmax_overall) {
	/* (nf * nframes)-threaded kernel */

	int total_offset = blockIdx.x * blockDim.x + threadIdx.x;
	int frm = total_offset / frame_size;	/* Frame number */
	int f   = total_offset % frame_size;	/* facet number */
	int pn, k, i, i1, i2, j, j1, j2, imin, imax, jmin, jmax;
	double n[3], v0[3], v1[3], v2[3], x[3], s, t, z, den;
	float imin_dbl, imax_dbl, jmin_dbl, jmax_dbl;
	int3 fidx;

	/* Variables that used to be in the init kernel or __device__ */
	double oa[3][3], usrc[3];
	//__syncthreads();
	if ((f < frame_size) && (frm < nframes)) {

		/* Take care of what was previously in init kernel */
		//if (f == 0) {
			dev_mtrnsps(oa, pos[frm]->ae);
			pn = pos[frm]->n;

			if (src) {
				/* We're viewing the model from the sun: at the center of each pixel
				 * in the projected view, we want cos(incidence angle), distance from
				 * the COM towards the sun, and the facet number.                */
				dev_mmmul(oa, pos[frm]->se, oa); /* oa takes ast into sun coords           */
			} else {
				/* We're viewing the model from Earth: at the center of each POS pixel
				 * we want cos(scattering angle), distance from the COM towards Earth,
				 * and the facet number.  For bistatic situations (lightcurves) we also
					 want cos(incidence angle) and the unit vector towards the source.     */
				dev_mmmul(oa, pos[frm]->oe, oa); /* oa takes ast into obs coords */
				if (pos[frm]->bistatic) {
					usrc[0] = usrc[1] = 0.0; /* unit vector towards source */
					usrc[2] = 1.0;
					dev_cotrans3(usrc, pos[frm]->se, usrc, -1);
					dev_cotrans3(usrc, pos[frm]->oe, usrc, 1); /* in observer coordinates */
				}
			}
//		}
//		__syncthreads();

		fidx.x = _verts->f[f].v[0];
		fidx.y = _verts->f[f].v[1];
		fidx.z = _verts->f[f].v[2];

		/* Get the normal to this facet in body-fixed (asteroid) coordinates
		 * and convert it to observer coordinates     */
		for (i = 0; i <= 2; i++)
			n[i] = _verts->f[f].n[i];

		dev_cotrans3(n, oa, n, 1);

		/* Consider this facet further only if its normal points somewhat
		 * towards the observer rather than away         */
		if (n[2] > 0.0) {
			/* Convert the three sets of vertex coordinates from body to ob-
			 * server coordinates; orbit_offset is the center-of-mass offset
			 * (in observer coordinates) for this model at this frame's epoch
			 * due to orbital motion, in case the model is half of a binary
			 * system.  */
			dev_cotrans3(v0, oa, _verts->v[fidx.x].x, 1);
			dev_cotrans3(v1, oa, _verts->v[fidx.y].x, 1);
			dev_cotrans3(v2, oa, _verts->v[fidx.z].x, 1);
			for (i = 0; i <= 2; i++) {
				v0[i] += orbit_offs.x;
				v1[i] += orbit_offs.y;
				v2[i] += orbit_offs.z;
			}

			/* Find rectangular region (in POS pixels) containing the projected
			 * facet - use floats in case model has illegal parameters and the
			 * pixel numbers exceed the limits for valid integers                         */
			imin_dbl = floor(MIN(v0[0],MIN(v1[0],v2[0])) / pos[frm]->km_per_pixel
							- SMALLVAL + 0.5);
			imax_dbl = floor(MAX(v0[0],MAX(v1[0],v2[0])) / pos[frm]->km_per_pixel
							+ SMALLVAL + 0.5);
			jmin_dbl = floor(MIN(v0[1],MIN(v1[1],v2[1])) / pos[frm]->km_per_pixel
							- SMALLVAL + 0.5);
			jmax_dbl = floor(MAX(v0[1],MAX(v1[1],v2[1])) / pos[frm]->km_per_pixel
							+ SMALLVAL + 0.5);
			imin = (imin_dbl < INT_MIN) ? INT_MIN : (int) imin_dbl;
			imax = (imax_dbl > INT_MAX) ? INT_MAX : (int) imax_dbl;
			jmin = (jmin_dbl < INT_MIN) ? INT_MIN : (int) jmin_dbl;
			jmax = (jmax_dbl > INT_MAX) ? INT_MAX : (int) jmax_dbl;

			/*  Set the outbnd flag if the facet extends beyond the POS window  */
			if ((imin < (-pn)) || (imax > pn) || (jmin < (-pn))	|| (jmax > pn))
				pvsoutbnd = 1;

			/* Figure out if facet projects at least partly within POS window;
			 * if it does, look at each "contained" POS pixel and get the
			 * z-coordinate and cos(scattering angle)           */
			i1 = MAX(imin, -pn);		j1 = MAX(jmin, -pn);
			i2 = MIN(imax,  pn);		j2 = MIN(jmax,  pn);

			if (i1 > pn || i2 < -pn || j1 > pn || j2 < -pn) {

				/* Facet is entirely outside the POS frame: just keep track of
				 * changed POS region     */
				dev_POSrect3(pos[frm], src, imin_dbl, imax_dbl, jmin_dbl, jmax_dbl,
						ijminmax_overall, frm);

			} else {

				dev_POSrect3(pos[frm], src, (float)i1, (float)i2, (float)j1,
						(float)j2, ijminmax_overall, frm);

				/* Facet is at least partly within POS frame: find all POS
				 * pixels whose centers project onto this facet  */
				int pxa;	// pixel address in 1D float array zzf etc
				for (i = i1; i <= i2; i++) {
					x[0] = i * pos[frm]->km_per_pixel;
					for (j = j1; j <= j2; j++) {
						x[1] = j * pos[frm]->km_per_pixel;

						/* Calculate the pixel address for 1D arrays */
						pxa = (j+pn)*(2*pn+1)+(i+pn);

						/* Compute parameters s(x,y) and t(x,y) which define a
						 * facet's surface as
						 *         z = z0 + s*(z1-z0) + t*(z2-z1)
						 * where z0, z1, and z2 are the z-coordinates at the
						 * vertices. The conditions 0 <= s <= 1 and
						 * 0 <= t <= s require the POS pixel center to be
						 * "within" the (projected) perimeter of facet f.    */
						den = 1	/ ((v1[0] - v0[0]) * (v2[1] - v1[1])
								 - (v2[0] - v1[0]) * (v1[1] - v0[1]));
						s = ((x[0] - v0[0]) * (v2[1] - v1[1])
						  - (v2[0] - v1[0]) * (x[1] - v0[1])) * den;

						if ((s >= -SMALLVAL) && (s <= 1.0 + SMALLVAL)) {

							t = ((v1[0] - v0[0]) * (x[1] - v0[1])
							    - (x[0] - v0[0]) * (v1[1]- v0[1])) * den;
							if ((t >= -SMALLVAL) && (t <= s + SMALLVAL)) {

								/* Compute z-coordinate of pixel center: its
								 * distance measured from the origin towards
								 * Earth.    */
								z = v0[2] + s*(v1[2]-v0[2]) + t*(v2[2]-v1[2]);

								/* If fac[i][j] is >= 0, pixel [i][j] was al-
								 * ready assigned values during a previous call
								 * to posvis for a different model component.
								 * If so, override only if the current component
								 * is blocking our view of (i.e., is closer to
								 * us than) the previous one.   */

								/* Following line replaces the previous if check
								 * for z > zz[i][j]
								 * atomicMaxf returns the value that was sitting
								 * at zzf[pxa] at time of call.  So if that value
								 * matches the z we compared to*/
								float old;
								if (src)
									old = atomicMaxf(&pos[frm]->zill_s[pxa], z);
								else
									old = atomicMaxf(&pos[frm]->z_s[pxa], z);

								if (old < z || pos[frm]->fill[i][j] < 0 ||
										pos[frm]->f[i][j] < 0) {

									/* Next line assigns distance of POS pixel
									 * center from COM towards Earth; that is,
									 * by changing zz,it changes pos->z or
									 * pos->zill                */
									/* following line is a first time z calc
									 * for this pixel  */
									if ( (pos[frm]->fill[i][j] < 0) || (pos[frm]->f[i][j] < 0)){
										if (src)	atomicExch(&pos[frm]->zill_s[pxa], z);
										else 		atomicExch(&pos[frm]->z_s[pxa], z);
									}

									if (pvs_smooth) {

										/* Get pvs_smoothed version of facet unit
										 * normal: Take the linear combination
										 * of the three vertex normals; trans-
										 * form from body to observer coordina-
										 * tes; and make sure that it points
										 * somewhat in our direction.         */
										for (k = 0; k <= 2; k++)
											n[k] =	_verts->v[fidx.x].n[k]
											 + s * (_verts->v[fidx.y].n[k] - _verts->v[fidx.x].n[k])
											 + t * (_verts->v[fidx.z].n[k]	- _verts->v[fidx.y].n[k]);
										dev_cotrans3(n, oa, n, 1);
										dev_normalize(n);
									}

									/* Determine scattering and/or incidence
									 * angles. Next lines change pos->cose/
									 * cosill. If bistatic (lightcurves), where
									 * we are viewing from Earth (src = 0),
									 * pos->cosi is also changed.                 */
									if (n[2] > 0.0) {
										if (src) atomicExch(&pos[frm]->cosill_s[pxa], n[2]);
										else	 atomicExch(&pos[frm]->cose_s[pxa], n[2]);
										if ((!src) && (pos[frm]->bistatic)) {
											float temp = (float)dev_dot(n,usrc);
											atomicExch(&pos[frm]->cosi_s[pxa], temp);
											if (pos[frm]->cosi_s[pxa] <= 0.0)
												pos[frm]->cose_s[pxa] = 0.0;
										}
									}

									/*  Keep track of the changed POS region  */
//									dev_POSrect2(pos, src, (double) i,
//											(double) i, (double) j, (double) j);

									/* Next lines change pos->body/bodyill,
									 * pos->comp/compill, pos->f/fill          */
									if (src) {
										pos[frm]->bodyill[i][j] = body;
										pos[frm]->compill[i][j] = comp;
										pos[frm]->fill[i][j] = f;
									} else {
										pos[frm]->body[i][j] = body;
										pos[frm]->comp[i][j] = comp;
										pos[frm]->f[i][j] = f;
									}

								} /* end if (no other facet yet blocks this facet from view) */
							} /* end if 0 <= t <= s (facet center is "in" this POS pixel) */
						} /* end if 0 <= s <= 1 */
					} /* end j-loop over POS rows */
				} /* end i-loop over POS columns */
			} /* end else of if (i1 > pos->n || i2 < -pos->n || j1 > pos->n || j2 < -pos->n) */
		} /* End if (n[2] > 0.0) */
	} /* end if (f < nf) */
}
__global__ void posvis_set_logfactor_af_krnl(struct pos_t **pos,
		int nframes, float4 *ijminmax_overall) {
	/* nframes-threaded kernel, but multiple streams */
	/*  If the target extends beyond the POS frame, set pos->posbnd_logfactor equal
	 to the logarithm of the ratio of the number of pixels in a frame extended to
	 include the entire target, divided by the number of pixels in the actual frame  */
	/* ijminmax_overall.w = imin_overall
	 * ijminmax_overall.x = imax_overall
	 * ijminmax_overall.y = jmin_overall
	 * ijminmax_overall.z = jmax_overall	 */
	int frm = threadIdx.x;
	if (frm < nframes) {
		if (pvsoutbnd) {
			double xfactor, yfactor;
			xfactor = (MAX(ijminmax_overall[frm].x,pos[frm]->n) - MIN(ijminmax_overall[frm].w,
					-pos[frm]->n)+ 1) / (2 * pos[frm]->n + 1);
			yfactor = (MAX(ijminmax_overall[frm].z,pos[frm]->n) - MIN(ijminmax_overall[frm].y,
					-pos[frm]->n)+ 1) / (2 * pos[frm]->n + 1);
			pos[frm]->posbnd_logfactor = log(xfactor * yfactor);
		}
	}
}

__host__ int posvis_af(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, float orbit_offset[3], int set, int nframes,
		int src, int body, int comp) {
	int nf, outbnd, n, nThreads;
	dim3 BLK,THD;

	struct pos_t **pos;
	float4 *ijminmax_overall;
	float3 orbit_offs;
	orbit_offs.x = orbit_offset[0];
	orbit_offs.y = orbit_offset[1];
	orbit_offs.z = orbit_offset[2];

	cudaCalloc((void**)&pos, sizeof(struct pos_t*), nframes);
	cudaCalloc((void**)&ijminmax_overall, sizeof(float4), nframes);

	/* Launch init kernel */
	THD.x = nframes;
	posvis_init_af_krnl<<<1,THD>>>(dpar, dmod, ddat, set,
		nframes, pos, ijminmax_overall, n);
	checkErrorAfterKernelLaunch("posvis_init_streams_krnl");
	gpuErrchk(cudaMemcpyFromSymbol(&nf, pvs_nf, sizeof(int), 0,
			cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&n, pvs_n, sizeof(int), 0,
			cudaMemcpyDeviceToHost));

	/* Configure and launch the facet kernel. Which kernel gets launched
	 * depends on flags DynProc (dynamic processing) and POSVIS_SEPARATE
	 * (xlim and ylim get calculated via a separate parallel reduction in an
	 * additional kernel after the facet kernel) */
	nThreads = nf * nframes;
	BLK.x = floor((maxThreadsPerBlock - 1 + nThreads) / maxThreadsPerBlock);
	THD.x = maxThreadsPerBlock;

	posvis_facet_af_krnl<<<BLK,THD>>>(pos, src,body,comp, orbit_offs,
			nThreads, nf, nframes, ijminmax_overall);
	checkErrorAfterKernelLaunch("posvis_facet_af_krnl");
	deviceSyncAfterKernelLaunch("posvis_facet_af_krnl");
	gpuErrchk(cudaMemcpyFromSymbol(&outbnd, pvsoutbnd, sizeof(outbnd), 0,
		cudaMemcpyDeviceToHost));

	/* Launch single-threaded kernel to set the pos logfactor */
	posvis_set_logfactor_af_krnl<<<1,1>>>(pos, nframes, ijminmax_overall);
	checkErrorAfterKernelLaunch("posvis_set_logfactor_af_krnl");

	/* Free temp arrays */
	cudaFree(ijminmax_overall);
	cudaFree(pos);

	return outbnd;
}
