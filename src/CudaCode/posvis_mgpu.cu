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
		int hf,
		int start,
		int end,
		int src) {

	/* Single-threaded, streamed, multi-GPU kernel. GPU0 handles even frames,
	 * GPU1 handles odd frames */
	if (threadIdx.x == 0) {

		ijminmax_overall[hf].w = ijminmax_overall[hf].y = HUGENUMBER;
		ijminmax_overall[hf].x = ijminmax_overall[hf].z = -HUGENUMBER;
		pos[hf]->posbnd_logfactor = 0.0;

		dev_mtrnsps3(oa, pos[hf]->ae, hf);

		if (src) {
			/* We're viewing the model from the sun: at the center of each pixel
			 * in the projected view, we want cos(incidence angle), distance from
			 * the COM towards the sun, and the facet number.                */
			dev_mmmul3(oa, pos[hf]->se, oa, hf); /* oa takes ast into sun coords           */
		} else {
			/* We're viewing the model from Earth: at the center of each POS pixel
			 * we want cos(scattering angle), distance from the COM towards Earth,
			 * and the facet number.  For bistatic situations (lightcurves) we also
									 want cos(incidence angle) and the unit vector towards the source.     */
			dev_mmmul3(oa, pos[hf]->oe, oa, hf); /* oa takes ast into obs coords */
			if (pos[hf]->bistatic) {
				usrc[hf].x = usrc[hf].y = 0.0; /* unit vector towards source */
				usrc[hf].z = 1.0;
				dev_cotrans9(&usrc[hf], pos[hf]->se, usrc[hf], -1);
				dev_cotrans9(&usrc[hf], pos[hf]->oe, usrc[hf], 1); /* in observer coordinates */
			}
		}
		outbndarr[hf] = 0;
	}
}

__global__ void posvis_mgpu_init_krnl(
		struct par_t *dpar,
		struct pos_t **pos,
		float4 *ijminmax_overall,
		float3 *oa,
		float3 *usrc,
		int *outbndarr,
		int c,
		int size,
		int oddflg,
		int src) {

	/* nfrm_half0/nfrm_half1-threaded kernel for multi-GPU operation. GPU0
	 * handles even frames, GPU1 handles odd frames */
	int hf = blockIdx.x * blockDim.x + threadIdx.x;
	int f = 2*hf + oddflg;

	if (hf < size) {

		ijminmax_overall[hf].w = ijminmax_overall[hf].y = HUGENUMBER;
		ijminmax_overall[hf].x = ijminmax_overall[hf].z = -HUGENUMBER;
		pos[hf]->posbnd_logfactor = 0.0;

		dev_mtrnsps3(oa, pos[hf]->ae, hf);

		if (src) {
			/* We're viewing the model from the sun: at the center of each pixel
			 * in the projected view, we want cos(incidence angle), distance from
			 * the COM towards the sun, and the facet number.                */
			dev_mmmul3(oa, pos[hf]->se, oa, hf); /* oa takes ast into sun coords           */
		} else {
			/* We're viewing the model from Earth: at the center of each POS pixel
			 * we want cos(scattering angle), distance from the COM towards Earth,
			 * and the facet number.  For bistatic situations (lightcurves) we also
									 want cos(incidence angle) and the unit vector towards the source.     */
			dev_mmmul3(oa, pos[hf]->oe, oa, hf); /* oa takes ast into obs coords */
			if (pos[hf]->bistatic) {
				usrc[hf].x = usrc[hf].y = 0.0; /* unit vector towards source */
				usrc[hf].z = 1.0;
				dev_cotrans9(&usrc[hf], pos[hf]->se, usrc[hf], -1);
				dev_cotrans9(&usrc[hf], pos[hf]->oe, usrc[hf], 1); /* in observer coordinates */
			}
		}
		outbndarr[hf] = 0;
	}
}

__global__ void mposvis_facet_krnl(
		struct par_t *dpar,
		struct mod_t *dmod,
		struct pos_t **pos,
		float4 *ijminmax_overall,
		float3 orbit_offs,
		float3 *oa,
		float3 *usrc,
		int src,
		int body,
		int comp,
		int nfacets,
		int frm,
		int hf,
		int *outbndarr) {

	/* nf-thread kernel.  It is streamed and GPU0 does even frames while GPU1
	 * does odd frames. hf is the half-frame array index for those arrays that
	 * needed to be split in half and assigned to the responsible gpu. */

	int f = blockIdx.x * blockDim.x + threadIdx.x;
	int pxa, i, i1, i2, j, j1, j2, imin, imax, jmin, jmax;
	float imin_dbl, imax_dbl, jmin_dbl, jmax_dbl, old, s, t, z, den;
	int3 fidx;
	float3 n, v0, v1, v2, tv0, tv1, tv2, x;
	__shared__ int pn;
	__shared__ float kmpxl;

	if (threadIdx.x == 0) {
		pn = pos[hf]->n;
		kmpxl = __double2float_rn(pos[hf]->km_per_pixel);
	}
	__syncthreads();
	if (f < nfacets) {
		/* The following section transfers vertex coordinates from double[3]
		 * storage to float3		 */
		fidx.x = dmod->shape.comp[0].real.f[f].v[0];
		fidx.y = dmod->shape.comp[0].real.f[f].v[1];
		fidx.z = dmod->shape.comp[0].real.f[f].v[2];
		tv0.x = __double2float_rn(dmod->shape.comp[0].real.v[fidx.x].x[0]);
		tv0.y = __double2float_rn(dmod->shape.comp[0].real.v[fidx.x].x[1]);
		tv0.z = __double2float_rn(dmod->shape.comp[0].real.v[fidx.x].x[2]);
		tv1.x = __double2float_rn(dmod->shape.comp[0].real.v[fidx.y].x[0]);
		tv1.y = __double2float_rn(dmod->shape.comp[0].real.v[fidx.y].x[1]);
		tv1.z = __double2float_rn(dmod->shape.comp[0].real.v[fidx.y].x[2]);
		tv2.x = __double2float_rn(dmod->shape.comp[0].real.v[fidx.z].x[0]);
		tv2.y = __double2float_rn(dmod->shape.comp[0].real.v[fidx.z].x[1]);
		tv2.z = __double2float_rn(dmod->shape.comp[0].real.v[fidx.z].x[2]);
		v0.x = v0.y = v0.z = v1.x = v1.y = v1.z = v2.x = v2.y = v2.z = 0.0;

		/* Get the normal to this facet in body-fixed (asteroid) coordinates
		 * and convert it to observer coordinates     */
		n.x = dmod->shape.comp[0].real.f[f].n[0];
		n.y = dmod->shape.comp[0].real.f[f].n[1];
		n.z = dmod->shape.comp[0].real.f[f].n[2];

		dev_cotrans8(&n, oa, n, 1, hf);

		/* Consider this facet further only if its normal points somewhat
		 * towards the observer rather than away         */
		if (n.z > 0.0) {
			/* Convert the three sets of vertex coordinates from body to ob-
			 * server coordinates; orbit_offset is the center-of-mass offset
			 * (in observer coordinates) for this model at this frame's epoch
			 * due to orbital motion, in case the model is half of a binary
			 * system.  */
			dev_cotrans8(&v0, oa, tv0, 1, hf);
			dev_cotrans8(&v1, oa, tv1, 1, hf);
			dev_cotrans8(&v2, oa, tv2, 1, hf);

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
				outbndarr[hf] = 1;
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
						ijminmax_overall, hf);

			} else {

				dev_POSrect_gpu(pos,src, __int2float_rn(i1), __int2float_rn(i2),
						__int2float_rn(j1),	__int2float_rn(j2),
						ijminmax_overall, hf);

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
									old = atomicMaxf(&pos[hf]->zill_s[pxa], z);
								else
									old = atomicMaxf(&pos[hf]->z_s[pxa], z);

								if (old < z || pos[hf]->fill[i][j] < 0 ||
										pos[hf]->f[i][j] < 0) {

									/* Next line assigns distance of POS pixel
									 * center from COM towards Earth; that is,
									 * by changing zz,it changes pos->z or
									 * pos->zill                */
									/* following line is a first time z calc
									 * for this pixel  */
									if ( (pos[hf]->fill[i][j] < 0) || (pos[hf]->f[i][j] < 0)){
										if (src)	atomicExch(&pos[hf]->zill_s[pxa], z);
										else 		atomicExch(&pos[hf]->z_s[pxa], z);
									}

									if (dpar->pos_smooth) {
										/* Assign temp. normal components as float3 */
										tv0.x = __double2float_rn(dmod->shape.comp[0].real.v[fidx.x].n[0]);
										tv0.y = __double2float_rn(dmod->shape.comp[0].real.v[fidx.x].n[1]);
										tv0.z = __double2float_rn(dmod->shape.comp[0].real.v[fidx.x].n[2]);
										tv1.x = __double2float_rn(dmod->shape.comp[0].real.v[fidx.y].n[0]);
										tv1.y = __double2float_rn(dmod->shape.comp[0].real.v[fidx.y].n[1]);
										tv1.z = __double2float_rn(dmod->shape.comp[0].real.v[fidx.y].n[2]);
										tv2.x = __double2float_rn(dmod->shape.comp[0].real.v[fidx.z].n[0]);
										tv2.y = __double2float_rn(dmod->shape.comp[0].real.v[fidx.z].n[1]);
										tv2.z = __double2float_rn(dmod->shape.comp[0].real.v[fidx.z].n[2]);

										/* Get pvs_smoothed version of facet unit
										 * normal: Take the linear combination
										 * of the three vertex normals; trans-
										 * form from body to observer coordina-
										 * tes; and make sure that it points
										 * somewhat in our direction.         */

										n.x = tv0.x + s * (tv1.x - tv0.x) + t * (tv2.x - tv1.x);
										n.y = tv0.y + s * (tv1.y - tv0.y) + t * (tv2.y - tv1.y);
										n.z = tv0.z + s * (tv1.z - tv0.z) + t * (tv2.z - tv1.z);

										dev_cotrans8(&n, oa, n, 1, hf);
										dev_normalize2(n);
									}

									/* Determine scattering and/or incidence
									 * angles. Next lines change pos->cose/
									 * cosill. If bistatic (lightcurves), where
									 * we are viewing from Earth (src = 0),
									 * pos->cosi is also changed.                 */
									if (n.z > 0.0) {
										if (src)
											atomicExch(&pos[hf]->cosill_s[pxa], n.z);
										else
											atomicExch(&pos[hf]->cose_s[pxa], n.z);
										if ((!src) && (pos[hf]->bistatic)) {
											float temp = dev_dot4(n,usrc[hf]);
											atomicExch(&pos[hf]->cosi_s[pxa], temp);
											if (pos[hf]->cosi_s[pxa] <= 0.0)
												pos[hf]->cose_s[pxa] = 0.0;
										}
									}

									/* Next lines change pos->body/bodyill,
									 * pos->comp/compill, pos->f/fill          */
									if (src) {
										pos[hf]->bodyill[i][j] = body;
										pos[hf]->compill[i][j] = comp;
										pos[hf]->fill[i][j] = f;
									} else {
										pos[hf]->body[i][j] = body;
										pos[hf]->comp[i][j] = comp;
										pos[hf]->f[i][j] = f;
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
		int *outbndarr, float4 *ijminmax_overall, int f, int hf) {
	/* Single-threaded, streamed, multi-GPU kernel. Odd frames go to GPU1,
	 * even frames go to GPU0 */
	double xfactor, yfactor;
	if (threadIdx.x == 0) {
		if (outbndarr[hf]) {
			/* ijminmax_overall.w = imin_overall
			 * ijminmax_overall.x = imax_overall
			 * ijminmax_overall.y = jmin_overall
			 * ijminmax_overall.z = jmax_overall	 */
			xfactor = (MAX( ijminmax_overall[hf].x,  posn) -
					MIN( ijminmax_overall[hf].w, -posn) + 1) / (2*posn+1);
			yfactor = (MAX( ijminmax_overall[hf].z,  posn) -
					MIN( ijminmax_overall[hf].y, -posn) + 1) / (2*posn+1);
			pos[hf]->posbnd_logfactor = log(xfactor*yfactor);
		}
	}
}

__global__ void posvis_mgpu_outbnd_krnl(struct pos_t **pos, int *outbndarr,
		float4 *ijminmax_overall, int size) {
	/* nfrm_half0/nfrm_half1-threaded multi-GPU kernel. Odd frames go to GPU1,
	 * even frames go to GPU0 */
	int hf = blockIdx.x * blockDim.x + threadIdx.x;
	int posn;
	double xfactor, yfactor;
	if (hf < size) {
		posn = pos[hf]->n;
		if (outbndarr[hf]) {
			/* ijminmax_overall.w = imin_overall
			 * ijminmax_overall.x = imax_overall
			 * ijminmax_overall.y = jmin_overall
			 * ijminmax_overall.z = jmax_overall	 */
			xfactor = (MAX( ijminmax_overall[hf].x,  posn) -
					MIN( ijminmax_overall[hf].w, -posn) + 1) / (2*posn+1);
			yfactor = (MAX( ijminmax_overall[hf].z,  posn) -
					MIN( ijminmax_overall[hf].y, -posn) + 1) / (2*posn+1);
			pos[hf]->posbnd_logfactor = log(xfactor*yfactor);
		}
	}
}

/* The posvis_mgpu function performs the same tasks as posvis_gpu. The principal
 * difference is this:
 *
 * 	- posvis_mgpu alternates gpus as it goes through frames of one set.
 * 	- this is done with cudaSetDevice() to set the desired gpu before the
 * 	  intended instruction.
 * 	- each gpu can perform write operations only on the constructs allocated
 * 	  on that device.
 * 	- Consequently, all arrays must be be split in half and allocated to the
 * 	  GPU device responsible for it.
 * 	- additionally, streams are associated with the gpu they were created on.
 * 	  Therefore, the streams used must go with the gpu that is current.
 */
__host__ int posvis_mgpu(
		struct par_t *dpar,
		struct mod_t *dmod,
		struct dat_t *ddat,
		struct pos_t **pos0,
		struct pos_t **pos1,
		float3 orbit_offset,
		int *hposn0,
		int *hposn1,
		int *outbndarr0,
		int *outbndarr1,
		int set,
		int nfrm_alloc,
		int src,
		int nf,
		int body,
		int comp,
		unsigned char type,
		cudaStream_t *gpu0_stream,	/* streams owned by gpu0    	*/
		cudaStream_t *gpu1_stream)	/* streams owned by gpu1		*/
{
	int outbnd=0, start, end;
	dim3 BLK,BLK_half0, BLK_half1, THD;
	cudaEvent_t start01, stop01;
	float milliseconds;
	float4 *ijminmax_overall0, *ijminmax_overall1;
	float3 *oa0, *oa1, *usrc0, *usrc1;

	/* Fix the frame offset if either of the two datasets is a lightcurve set */
	if (type == LGHTCRV)
		start = 1;	/* fixes the lightcurve offsets */
	else start = 0;

	int nfrm_half0 = nfrm_alloc/2 + nfrm_alloc%2;
	int nfrm_half1 = nfrm_alloc/2;
	int oasize0 = nfrm_half0 * 3;
	int oasize1 = nfrm_half1 * 3;
	int hf = 0;	/* Half-frame counter for the split arrays */
	int houtbnd0[nfrm_half0], houtbnd1[nfrm_half1];

	/* Launch parameters for the facet_streams kernel */
	THD.x = maxThreadsPerBlock;
	BLK.x = floor((THD.x - 1 + nf) / THD.x);
	BLK_half0.x = floor((THD.x - 1 + nfrm_half0)/THD.x);
	BLK_half1.x = floor((THD.x - 1 + nfrm_half1)/THD.x);

	/* GPU 0 allocation */
	gpuErrchk(cudaSetDevice(GPU0));	/* Make sure we're really on GPU0 */
	gpuErrchk(cudaMalloc((void**)&ijminmax_overall0, sizeof(float4) * nfrm_half0));
	gpuErrchk(cudaMalloc((void**)&oa0, sizeof(float3) * oasize0));
	gpuErrchk(cudaMalloc((void**)&usrc0, sizeof(float3) * nfrm_half0));

	/* GPU 1 allocation, then switch back to GPU0 */
	gpuErrchk(cudaSetDevice(GPU1));
	gpuErrchk(cudaMalloc((void**)&ijminmax_overall1, sizeof(float4) * nfrm_half1));
	gpuErrchk(cudaMalloc((void**)&oa1, sizeof(float3) * oasize1));
	gpuErrchk(cudaMalloc((void**)&usrc1, sizeof(float3) * nfrm_half1));

	if (TIMING) {
		/* Create the timer events */
		cudaEventCreate(&start01);
		cudaEventCreate(&stop01);
		cudaEventRecord(start01);
	}

	gpuErrchk(cudaSetDevice(GPU0));
	posvis_mgpu_init_krnl<<<BLK_half0,THD,0,gpu0_stream[0]>>>(dpar, pos0,
			ijminmax_overall0, oa0, usrc0, outbndarr0, 0, nfrm_half0, 0, src);
	gpuErrchk(cudaSetDevice(GPU1));
	posvis_mgpu_init_krnl<<<BLK_half1,THD,0,gpu1_stream[0]>>>(dpar, pos1,
			ijminmax_overall1, oa1, usrc1, outbndarr1, 0, nfrm_half1, 1, src);
	checkErrorAfterKernelLaunch("posvis_mgpu_init_krnl");

	hf = 0;
	gpuErrchk(cudaSetDevice(GPU0));
	for (int f=start; f<nfrm_alloc; f++) {
		/* Now the main facet kernel */
		mposvis_facet_krnl<<<BLK,THD, 0, gpu0_stream[hf]>>>(dpar, dmod, pos0,
				ijminmax_overall0, orbit_offset, oa0, usrc0, src, body, comp,
				nf, f, hf, outbndarr0);
		f++;	hf++;
		if (f>=nfrm_alloc) break;
	}

	hf = 0;
	gpuErrchk(cudaSetDevice(GPU1));
	for (int f=1+start; f<nfrm_alloc; f++) {
		/* Now the main facet kernel */
		mposvis_facet_krnl<<<BLK,THD, 0, gpu1_stream[hf]>>>(dpar, dmod, pos1,
				ijminmax_overall1, orbit_offset, oa1, usrc1, src, body, comp,
				nf, f, hf, outbndarr1);

		f++;	hf++;
		if (f>=nfrm_alloc)	break;
	}

	gpuErrchk(cudaSetDevice(GPU0));
	posvis_mgpu_outbnd_krnl<<<BLK_half0,THD,0,gpu0_stream[0]>>>(pos0,
			outbndarr0, ijminmax_overall0, nfrm_half0);
	gpuErrchk(cudaSetDevice(GPU1));
	posvis_mgpu_outbnd_krnl<<<BLK_half1,THD,0,gpu1_stream[0]>>>(pos1,
			outbndarr1, ijminmax_overall1, nfrm_half1);
	checkErrorAfterKernelLaunch("posvis_mgpu_outbnd_krnl");


	if (TIMING) {
		cudaEventRecord(stop01);
		cudaEventSynchronize(stop01);
		milliseconds = 0;
		cudaEventElapsedTime(&milliseconds, start01, stop01);
		printf("%i facets in posvis_cuda_2 in %3.3f ms with %i frames.\n",
				nf, milliseconds, nfrm_alloc);
	}
	/* Now asynchronously copy the outbnd arrays from each gpu back to host
	 * arrays.  If any of the frames were out of bounds, we'll return 1.   */
//	printf("nrm_half0: %i;     nfrmf_half1: %i\n", nfrm_half0, nfrm_half1);
//	fflush( stdout);
	gpuErrchk(cudaSetDevice(GPU0));
	gpuErrchk(cudaMemcpy(houtbnd0, outbndarr0, sizeof(int)*nfrm_half0,
			cudaMemcpyDeviceToHost));
	gpuErrchk(cudaSetDevice(GPU1));
	gpuErrchk(cudaMemcpy(houtbnd1, outbndarr1, sizeof(int)*nfrm_half1,
			cudaMemcpyDeviceToHost));

	for (int f=0; f<nfrm_half1; f++) {
		if (houtbnd0[f] || houtbnd1[f])
			outbnd = 1;
	}
	/* Check last frame in houtbnd0 if nframes is uneven */
	if (houtbnd0[nfrm_half0-1])
		outbnd = 1;


//	/* Debug start	 */
//	cudaSetDevice(GPU0);
//	/* Frame 0 */
//	hf = 0;
//	int npxls = (2*hposn0[hf]+1)*(2*hposn0[hf]+1);
//	dbg_print_pos_arrays_full(pos0, hf, (2*hf), npxls,  hposn0[hf]);
//	/* Frame 2 */
//	hf = 1;
//	npxls = (2*hposn0[hf]+1)*(2*hposn0[hf]+1);
//	dbg_print_pos_arrays_full(pos0, hf, (2*hf), npxls, hposn0[hf]);
////	/* Frame 4 */
////	hf = 2;
////	npxls = (2*hposn0[hf]+1)*(2*hposn0[hf]+1);
////	dbg_print_pos_arrays_full(pos0, hf, npxls, hposn0[hf]);
//
//	/* Frame 1 */
//	hf = 0;
//	npxls = (2*hposn1[hf]+1)*(2*hposn1[hf]+1);
//	cudaSetDevice(GPU1);
//	dbg_print_pos_arrays_full(pos1, hf, (2*hf+1), npxls, hposn1[hf]);
//	/* Frame 3 */
//	hf = 1;
//	npxls = (2*hposn1[hf]+1)*(2*hposn1[hf]+1);
//	cudaSetDevice(GPU1);
//	dbg_print_pos_arrays_full(pos1, hf, (2*hf+1), npxls, hposn1[hf]);
////	/* Frame 5 */
////	hf = 2;
////	npxls = (2*hposn1[hf]+1)*(2*hposn1[hf]+1);
////	cudaSetDevice(GPU1);
////	dbg_print_pos_arrays_full(pos1, hf, npxls, hposn1[hf]);
//
//	/* Debug end */



	/* Free GPU0 memory */
	gpuErrchk(cudaSetDevice(GPU0));
	cudaFree(ijminmax_overall0);
	cudaFree(oa0);
	cudaFree(usrc0);

	/* Free GPU1 memory */
	gpuErrchk(cudaSetDevice(GPU1));
	cudaFree(ijminmax_overall1);
	cudaFree(oa1);
	cudaFree(usrc1);

	/* Back to GPU0 and destroy timing events if they were used. */
	gpuErrchk(cudaSetDevice(GPU0));
	if (TIMING) {
		cudaEventDestroy(start01);
		cudaEventDestroy(stop01);	}

	return outbnd;
}
