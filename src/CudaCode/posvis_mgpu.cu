/*****************************************************************************************
 posvis_mgpu.cu

 Fill in the portion of a plane-of-sky image due to a particular model component: Assign
 each relevant POS pixel a z-value in observer coordinates (distance from the origin
 towards Earth) and a value of cos(scattering angle).

 Return 1 if any portion of this component lies outside the specified POS window,
 0 otherwise.

 If the "src" argument is true, the "observer" is the Sun rather than Earth, and
 "plane-of-sky" becomes "projection as viewed from the Sun."

 Modified 2017 May 27 by ME:
 This is a modified version of the posvis_gpu routine (which performs the
 calculations originally done by the posvis routine on an Nvidia GPU).
 This modified version uses two GPUs instead of just one.  This requires the
 presence of a dual-GPU card or multiple cards. Any Nvidia card can be paired
 up, but they must both meeth the minimum compute capability 3.5. Further it is
 extremely helpful to have both cards matched in performance, i.e. use two
 cards of the same make and model.

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
#include "../shape/head.h"
#include <limits.h>
}

__device__ int mposvis_streams_outbnd, mpvst_smooth;

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

__global__ void mposvis_init_krnl(
		struct par_t *dpar,
		struct pos_t **pos,
		float4 *ijminmax_overall,
		float3 *oa,
		float3 *usrc,
		int *outbndarr,
		int c,
		int f,
		int start,
		int end,
		int src) {

	/* Single-threaded, streamed kernel */
	if (threadIdx.x == 0) {
		if (f == start) {
			mposvis_streams_outbnd = 0;
			mpvst_smooth = dpar->pos_smooth;
		}
		ijminmax_overall[f].w = ijminmax_overall[f].y = HUGENUMBER;
		ijminmax_overall[f].x = ijminmax_overall[f].z = -HUGENUMBER;
		pos[f]->posbnd_logfactor = 0.0;

		dev_mtrnsps3(oa, pos[f]->ae, f);

		if (src) {
			/* We're viewing the model from the sun: at the center of each pixel
			 * in the projected view, we want cos(incidence angle), distance from
			 * the COM towards the sun, and the facet number.                */
			dev_mmmul3(oa, pos[f]->se, oa, f); /* oa takes ast into sun coords           */
		} else {
			/* We're viewing the model from Earth: at the center of each POS pixel
			 * we want cos(scattering angle), distance from the COM towards Earth,
			 * and the facet number.  For bistatic situations (lightcurves) we also
									 want cos(incidence angle) and the unit vector towards the source.     */
			dev_mmmul3(oa, pos[f]->oe, oa, f); /* oa takes ast into obs coords */
			if (pos[f]->bistatic) {
				usrc[f].x = usrc[f].y = 0.0; /* unit vector towards source */
				usrc[f].z = 1.0;
				dev_cotrans9(&usrc[f], pos[f]->se, usrc[f], -1);
				dev_cotrans9(&usrc[f], pos[f]->oe, usrc[f], 1); /* in observer coordinates */
			}
		}
		outbndarr[f] = 0;
	}
}
__global__ void mposvis_facet_krnl(
		struct pos_t **pos,
		struct vertices_t **verts,
		float4 *ijminmax_overall,
		float3 orbit_offs,
		float3 *oa,
		float3 *usrc,
		int src,
		int body,
		int comp,
		int nfacets,
		int frm,
		int smooth,
		int *outbndarr) {
	/* (nf * nframes)-threaded kernel.  This version eliminates as much double
	 * math as possible */

	int f = blockIdx.x * blockDim.x + threadIdx.x;
	int pxa, i, i1, i2, j, j1, j2, imin, imax, jmin, jmax;
	float imin_dbl, imax_dbl, jmin_dbl, jmax_dbl, old, s, t, z, den;
	int3 fidx;
	float3 n, v0, v1, v2, tv0, tv1, tv2, x;
	__shared__ int pn;
	__shared__ float kmpxl;

	if (threadIdx.x == 0) {
		pn = pos[frm]->n;
		kmpxl = __double2float_rn(pos[frm]->km_per_pixel);
	}

	if (f < nfacets) {
		/* The following section transfers vertex coordinates from double[3]
		 * storage to float3		 */
		fidx.x = verts[0]->f[f].v[0];
		fidx.y = verts[0]->f[f].v[1];
		fidx.z = verts[0]->f[f].v[2];
		tv0.x = __double2float_rn( verts[0]->v[fidx.x].x[0]);
		tv0.y = __double2float_rn(verts[0]->v[fidx.x].x[1]);
		tv0.z = __double2float_rn(verts[0]->v[fidx.x].x[2]);
		tv1.x = __double2float_rn(verts[0]->v[fidx.y].x[0]);
		tv1.y = __double2float_rn(verts[0]->v[fidx.y].x[1]);
		tv1.z = __double2float_rn(verts[0]->v[fidx.y].x[2]);
		tv2.x = __double2float_rn(verts[0]->v[fidx.z].x[0]);
		tv2.y = __double2float_rn(verts[0]->v[fidx.z].x[1]);
		tv2.z = __double2float_rn(verts[0]->v[fidx.z].x[2]);
		v0.x = v0.y = v0.z = v1.x = v1.y = v1.z = v2.x = v2.y = v2.z = 0.0;

		/* Get the normal to this facet in body-fixed (asteroid) coordinates
		 * and convert it to observer coordinates     */
		n.x = verts[0]->f[f].n[0];
		n.y = verts[0]->f[f].n[1];
		n.z = verts[0]->f[f].n[2];

		dev_cotrans8(&n, oa, n, 1, frm);

		/* Consider this facet further only if its normal points somewhat
		 * towards the observer rather than away         */
		if (n.z > 0.0) {
			/* Convert the three sets of vertex coordinates from body to ob-
			 * server coordinates; orbit_offset is the center-of-mass offset
			 * (in observer coordinates) for this model at this frame's epoch
			 * due to orbital motion, in case the model is half of a binary
			 * system.  */
			dev_cotrans8(&v0, oa, tv0, 1, frm);
			dev_cotrans8(&v1, oa, tv1, 1, frm);
			dev_cotrans8(&v2, oa, tv2, 1, frm);

			v0.x += orbit_offs.x;
			v0.y += orbit_offs.x;
			v0.z += orbit_offs.x;
			v1.x += orbit_offs.y;
			v1.y += orbit_offs.y;
			v1.z += orbit_offs.y;
			v2.x += orbit_offs.z;
			v2.y += orbit_offs.z;
			v2.z += orbit_offs.z;

			/* Find rectangular region (in POS pixels) containing the projected
			 * facet - use floats in case model has illegal parameters and the
			 * pixel numbers exceed the limits for valid integers                         */
			imin_dbl = floor(MIN(v0.x,MIN(v1.x,v2.x)) / kmpxl
							- SMALLVAL + 0.5);
			imax_dbl = floor(MAX(v0.x,MAX(v1.x,v2.x)) / kmpxl
							+ SMALLVAL + 0.5);
			jmin_dbl = floor(MIN(v0.y,MIN(v1.y,v2.y)) / kmpxl
							- SMALLVAL + 0.5);
			jmax_dbl = floor(MAX(v0.y,MAX(v1.y,v2.y)) / kmpxl
							+ SMALLVAL + 0.5);
			imin = (imin_dbl < INT_MIN) ? INT_MIN : (int) imin_dbl;
			imax = (imax_dbl > INT_MAX) ? INT_MAX : (int) imax_dbl;
			jmin = (jmin_dbl < INT_MIN) ? INT_MIN : (int) jmin_dbl;
			jmax = (jmax_dbl > INT_MAX) ? INT_MAX : (int) jmax_dbl;

			/*  Set the outbnd flag if the facet extends beyond the POS window  */
			if ((imin < (-pn)) || (imax > pn) || (jmin < (-pn))	|| (jmax > pn)) {
				mposvis_streams_outbnd = 1;
				outbndarr[f] = 1;
			}

			/* Figure out if facet projects at least partly within POS window;
			 * if it does, look at each "contained" POS pixel and get the
			 * z-coordinate and cos(scattering angle)           */
			i1 = MAX(imin, -pn);		j1 = MAX(jmin, -pn);
			i2 = MIN(imax,  pn);		j2 = MIN(jmax,  pn);

			if (i1 > pn || i2 < -pn || j1 > pn || j2 < -pn) {

				/* Facet is entirely outside the POS frame: just keep track of
				 * changed POS region     */
				dev_POSrect_gpu(pos, src, imin_dbl, imax_dbl, jmin_dbl, jmax_dbl,
						ijminmax_overall, frm);

			} else {

				dev_POSrect_gpu(pos, src, __double2float_rn(i1),
						__double2float_rn(i2), __double2float_rn(j1),
						__double2float_rn(j2), ijminmax_overall, frm);

				/* Facet is at least partly within POS frame: find all POS
				 * pixels whose centers project onto this facet  */
				for (i = i1; i <= i2; i++) {
					x.x = i * kmpxl;
					for (j = j1; j <= j2; j++) {
						x.y = j * kmpxl;

						/* Calculate the pixel address for 1D arrays */
						pxa = (j+pn) * (2*pn + 1) + (i+pn);

						/* Compute parameters s(x,y) and t(x,y) which define a
						 * facet's surface as
						 *         z = z0 + s*(z1-z0) + t*(z2-z1)
						 * where z0, z1, and z2 are the z-coordinates at the
						 * vertices. The conditions 0 <= s <= 1 and
						 * 0 <= t <= s require the POS pixel center to be
						 * "within" the (projected) perimeter of facet f.    */
						den = 1	/ ((v1.x - v0.x) * (v2.y - v1.y)
								 - (v2.x - v1.x) * (v1.y - v0.y));
						s = ((x.x - v0.x) * (v2.y - v1.y)
						  - (v2.x - v1.x) * (x.y - v0.y)) * den;

						if ((s >= -SMALLVAL) && (s <= 1.0 + SMALLVAL)) {

							t = ((v1.x - v0.x) * (x.y - v0.y)
							    - (x.x- v0.x) * (v1.y- v0.y)) * den;
							if ((t >= -SMALLVAL) && (t <= s + SMALLVAL)) {

								/* Compute z-coordinate of pixel center: its
								 * distance measured from the origin towards
								 * Earth.    */
								z = v0.z + s*(v1.z-v0.z) + t*(v2.z-v1.z);

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

									if (mpvst_smooth) {
										/* Assign temp. normal components as float3 */
										tv0.x = __double2float_rn(verts[0]->v[fidx.x].n[0]);
										tv0.y = __double2float_rn(verts[0]->v[fidx.x].n[1]);
										tv0.z = __double2float_rn(verts[0]->v[fidx.x].n[2]);
										tv1.x = __double2float_rn(verts[0]->v[fidx.y].n[0]);
										tv1.y = __double2float_rn(verts[0]->v[fidx.y].n[1]);
										tv1.z = __double2float_rn(verts[0]->v[fidx.y].n[2]);
										tv2.x = __double2float_rn(verts[0]->v[fidx.z].n[0]);
										tv2.y = __double2float_rn(verts[0]->v[fidx.z].n[1]);
										tv2.z = __double2float_rn(verts[0]->v[fidx.z].n[2]);

										/* Get pvs_smoothed version of facet unit
										 * normal: Take the linear combination
										 * of the three vertex normals; trans-
										 * form from body to observer coordina-
										 * tes; and make sure that it points
										 * somewhat in our direction.         */

										n.x = tv0.x + s * (tv1.x - tv0.x) + t * (tv2.x - tv1.x);
										n.y = tv0.y + s * (tv1.y - tv0.y) + t * (tv2.y - tv1.y);
										n.z = tv0.z + s * (tv1.z - tv0.z) + t * (tv2.z - tv1.z);

										dev_cotrans8(&n, oa, n, 1, frm);
										dev_normalize2(n);
									}

									/* Determine scattering and/or incidence
									 * angles. Next lines change pos->cose/
									 * cosill. If bistatic (lightcurves), where
									 * we are viewing from Earth (src = 0),
									 * pos->cosi is also changed.                 */
									if (n.z > 0.0) {
										if (src)
											atomicExch(&pos[frm]->cosill_s[pxa], n.z);
										else
											atomicExch(&pos[frm]->cose_s[pxa], n.z);
										if ((!src) && (pos[frm]->bistatic)) {
											float temp = dev_dot4(n,usrc[frm]);
											atomicExch(&pos[frm]->cosi_s[pxa], temp);
											if (pos[frm]->cosi_s[pxa] <= 0.0)
												pos[frm]->cose_s[pxa] = 0.0;
										}
									}

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
__global__ void mposvis_outbnd_krnl(struct pos_t **pos, int posn,
		int *outbndarr, float4 *ijminmax_overall, int f) {
	/* Single-threaded, streamed kernel */
	double xfactor, yfactor;
	if (threadIdx.x == 0) {
		if (outbndarr[f]) {
			/* ijminmax_overall.w = imin_overall
			 * ijminmax_overall.x = imax_overall
			 * ijminmax_overall.y = jmin_overall
			 * ijminmax_overall.z = jmax_overall	 */
			xfactor = (MAX( ijminmax_overall[f].x,  posn) -
					MIN( ijminmax_overall[f].w, -posn) + 1) / (2*posn+1);
			yfactor = (MAX( ijminmax_overall[f].z,  posn) -
					MIN( ijminmax_overall[f].y, -posn) + 1) / (2*posn+1);
			pos[f]->posbnd_logfactor = log(xfactor*yfactor);
		}
	}
}

/* The posvis_mgpu function performs the same tasks as posvis_gpu. The principal
 * difference is this:
 *
 * 	- posvis_mgpu alternates gpus as it goes through frames of one set.
 * 	- this is done with cudaSetDevice() to set the desired gpu before the
 * 	  intended instruction.
 * 	- additionally, streams are associated with the gpu they were created on.
 * 	  Therefore, the streams used must go with the gpu that is current.
 */
__host__ int posvis_mgpu(
		struct par_t *dpar,
		struct mod_t *dmod,
		struct dat_t *ddat,
		struct pos_t **pos0,
		struct pos_t **pos1,
		struct vertices_t **verts,
		float3 orbit_offset,
		int *posn0,
		int *posn1,
		int *outbndarr0,
		int *outbndarr1,
		int set,
		int nframes,
		int src,
		int nf,
		int body,
		int comp,
		unsigned char type,
		cudaStream_t *gpu0_stream,	/* streams owned by gpu0    	*/
		cudaStream_t *gpu1_stream)	/* streams owned by gpu1		*/
{
	int outbnd, smooth, start, end, frames_alloc;
	dim3 BLK,THD;
	cudaEvent_t start01, stop01;
	float milliseconds;
	float4 *ijminmax_overall;
	float3 *oa, *usrc;

	/* Launch parameters for the facet_streams kernel */
	THD.x = maxThreadsPerBlock;
	BLK.x = floor((THD.x - 1 + nf) / THD.x);

	/* Fix the frame offset if either of the two datasets is a lightcurve set */
	if (type == LGHTCRV) {
		start = 1;	/* fixes the lightcurve offsets */
		end = nframes + 1;
		frames_alloc = nframes + 1;
	} else {
		start = 0;
		end = nframes;
		frames_alloc = nframes;	}

	int oasize = frames_alloc*3;

	/* Allocate temporary arrays/structs */
	gpuErrchk(cudaMalloc((void**)&ijminmax_overall, sizeof(float4) * frames_alloc));
	gpuErrchk(cudaMalloc((void**)&oa, sizeof(float3) * oasize));
	gpuErrchk(cudaMalloc((void**)&usrc, sizeof(float3) * frames_alloc));

	if (TIMING) {
		/* Create the timer events */
		cudaEventCreate(&start01);
		cudaEventCreate(&stop01);
		cudaEventRecord(start01);
	}

	for (int f=start; f<end; f++) {

		gpuErrchk(cudaSetDevice(GPU0));
		/* Initialize via single-thread kernel first */
		mposvis_init_krnl<<<1,1,0,gpu0_stream[f-start]>>>(dpar,
				pos, ijminmax_overall, oa, usrc, outbndarr, comp, f, start,
				end, src);
//		checkErrorAfterKernelLaunch("mposvis_init_krnl on device 0");
		/* Now the main facet kernel */
		mposvis_facet_krnl<<<BLK,THD, 0, gpu0_stream[f-start]>>>(pos, verts,
				ijminmax_overall, orbit_offset, oa, usrc,	src, body, comp,
				nf, f, smooth, outbndarr);

		/* Take care of any posbnd flags */
		mposvis_outbnd_krnl<<<1,1,0,gpu0_stream[f-start]>>>(pos, posn[f],
				outbndarr, ijminmax_overall, f);

		/* Switch devices */
		gpuErrchk(cudaSetDevice(GPU1));

		/* Go to next frame and check if we're still in the set */
		f++;
		if (f==end)	break;


		/* Initialize via single-thread kernel first */
		mposvis_init_krnl<<<1,1,0,gpu1_stream[f-start]>>>(dpar,
				pos, ijminmax_overall, oa, usrc, outbndarr, comp, f, start,
				end, src);

		/* Now the main facet kernel */
		mposvis_facet_krnl<<<BLK,THD, 0, gpu1_stream[f-start]>>>(pos, verts,
				ijminmax_overall, orbit_offset, oa, usrc,	src, body, comp,
				nf, f, smooth, outbndarr);

		/* Take care of any posbnd flags */
		mposvis_outbnd_krnl<<<1,1,0,gpu1_stream[f-start]>>>(pos, posn[f],
				outbndarr, ijminmax_overall, f);
	}

	if (TIMING) {
		cudaEventRecord(stop01);
		cudaEventSynchronize(stop01);
		milliseconds = 0;
		cudaEventElapsedTime(&milliseconds, start01, stop01);
		printf("%i facets in posvis_cuda_2 in %3.3f ms with %i frames.\n", nf, milliseconds, nframes);
	}

	gpuErrchk(cudaMemcpyFromSymbol(&outbnd, mposvis_streams_outbnd, sizeof(outbnd), 0,
			cudaMemcpyDeviceToHost));

	/* Free temp arrays, destroy streams and timers, as applicable */
	cudaFree(ijminmax_overall);
	cudaFree(oa);
	cudaFree(usrc);

	if (TIMING) {
		cudaEventDestroy(start01);
		cudaEventDestroy(stop01);
	}

	/* Make sure the original GPU is set back to current device */
	gpuErrchk(cudaSetDevice(GPU0));
	return outbnd;
}

/* The posvis_mgpu function performs the same tasks as posvis_gpu. The principal
 * difference is this:
 *
 * 	- posvis_gpu operates on one dataset at a time and streams the frame
 * 	  calculation kernels,
 * 	- posvis_mgpu operates on two datasets at a time by assigning one dataset
 * 	  to gpu0 and the other dataset to gpu1,
 * 	- posvis_mgpu must calculate parameters needed for the facet kernel for
 * 	  two sets at a time,
 * 	- gpu is switched with cudaSetDevice(),
 * 	- Streams belong to the gpu that was active when they were created. So we
 * 	  have to use two sets of streams as input arguments.
 */
//__host__ int posvis_mgpu(
//		struct par_t *dpar,
//		struct mod_t *dmod,
//		struct dat_t *ddat,
//		struct pos_t **pos0,		/* all pos frames for gpu0 dataset */
//		struct pos_t **pos1,		/* all pos frames for gpu1 dataset */
//		struct vertices_t **verts,
//		float3 orbit_offset,
//		int *posn,
//		int src,
//		int nf,
//		int body,
//		int comp,
//		int *outbndarr0,			/* for dataset on gpu0 			*/
//		int *outbndarr1,			/* for dataset on gpu1 			*/
//		int set0,					/* dataset for gpu0    			*/
//		int set1,					/* dataset for gpu1	   			*/
//		int nframes0,				/* frames in dataset for gpu0 	*/
//		int nframes1,				/* frames in dataset for gpu1 	*/
//		unsigned char type0,		/* type of dataset for gpu0 	*/
//		unsigned char type1,		/* type of dataset for gpu1 	*/
//		cudaStream_t *gpu0_stream,	/* streams owned by gpu0    	*/
//		cudaStream_t *gpu1_stream)	/* streams owned by gpu1		*/
//{
//
//
//	int outbnd0, outbnd1, smooth, start0, start1, end0, end1, frames_alloc0,
//		frames_alloc1;
//	dim3 BLK,THD;
//	cudaEvent_t start01, stop01;
//	float milliseconds;
//	float4 *ijminmax_overall0, *ijminmax_overall1;
//	float3 *oa0, *oa1, *usrc0, *usrc1;
//
//	/* Launch parameters for the facet_streams kernel */
//	THD.x = maxThreadsPerBlock;
//	BLK.x = floor((THD.x - 1 + nf) / THD.x);
//
//	/* Fix the frame offset if either of the two datasets is a lightcurve set */
//	if (type0 == LGHTCRV) {
//		start0 = 1;	/* fixes the lightcurve offsets */
//		end0 = nframes0 + 1;
//		frames_alloc0 = nframes0 + 1;
//	} else {
//		start0 = 0;
//		end0 = nframes0;
//		frames_alloc0 = nframes0;	}
//	if (type1 == LGHTCRV) {
//		start1 = 1; /* fixes lightcurve offset in gpu1 dataset */
//		end1 = nframes1 + 1;
//		frames_alloc1 = nframes1 + 1;
//	} else {
//		start1= 0;
//		end1 = nframes1;
//		frames_alloc1 = nframes1;	}
//
//	int oasize0 = frames_alloc0*3;
//	int oasize1 = frames_alloc1*3;
//
//
//	/* Allocate temporary arrays/structs */
//	gpuErrchk(cudaMalloc((void**)&ijminmax_overall0, sizeof(float4) * frames_alloc0));
//	gpuErrchk(cudaMalloc((void**)&ijminmax_overall1, sizeof(float4) * frames_alloc1));
//	gpuErrchk(cudaMalloc((void**)&oa0, sizeof(float3) * oasize0));
//	gpuErrchk(cudaMalloc((void**)&oa1, sizeof(float3) * oasize1));
//	gpuErrchk(cudaMalloc((void**)&usrc0, sizeof(float3) * frames_alloc0));
//	gpuErrchk(cudaMalloc((void**)&usrc1, sizeof(float3) * frames_alloc1));
//
//	if (TIMING) {
//		/* Create the timer events */
//		cudaEventCreate(&start01);
//		cudaEventCreate(&stop01);
//		cudaEventRecord(start01);
//	}
//
//	for (int f=start; f<end; f++) {
//
//		/* Initialize via single-thread kernel first */
//		posvis_init_krnl<<<1,1,0,posvis_stream[f-start]>>>(dpar,
//				pos, ijminmax_overall, oa, usrc, outbndarr, comp, f, start,
//				end, src);
//
//		/* Now the main facet kernel */
//		posvis_facet_krnl<<<BLK,THD, 0, posvis_stream[f-start]>>>(pos, verts,
//				ijminmax_overall, orbit_offset, oa, usrc,	src, body, comp,
//				nf, f, smooth, outbndarr);
//
//		/* Take care of any posbnd flags */
//		posvis_outbnd_krnl<<<1,1,0,posvis_stream[f-start]>>>(pos, posn[f],
//				outbndarr, ijminmax_overall, f);
//	}
//
//	if (TIMING) {
//		cudaEventRecord(stop01);
//		cudaEventSynchronize(stop01);
//		milliseconds = 0;
//		cudaEventElapsedTime(&milliseconds, start01, stop01);
//		printf("%i facets in posvis_cuda_2 in %3.3f ms with %i frames.\n", nf, milliseconds, nframes);
//	}
//	checkErrorAfterKernelLaunch("The three posvis_cuda_streams2 kernels");
//	gpuErrchk(cudaMemcpyFromSymbol(&outbnd, mposvis_streams_outbnd, sizeof(outbnd), 0,
//			cudaMemcpyDeviceToHost));
//
//	/* Free temp arrays, destroy streams and timers, as applicable */
//
//	cudaFree(ijminmax_overall0);
//	cudaFree(ijminmax_overall1);
//	cudaFree(oa0);
//	cudaFree(oa1);
//	cudaFree(usrc0);
//	cudaFree(usrc1);
//
//
//	if (TIMING) {
//		cudaEventDestroy(start01);
//		cudaEventDestroy(stop01);
//	}
//	return outbnd;
//}
