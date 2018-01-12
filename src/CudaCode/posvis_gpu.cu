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
#include "../shape/head.h"
#include <limits.h>
}
__device__ int posvis_streams_outbnd, pvst_smooth, dbg_cntr=0;

/* Note that the following custom atomic functions must be declared in each
 * file it is needed (consequence of being a static device function) */

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
__device__ static float atomicMax64(double* address, double val)
{
	unsigned long long* address_as_i = (unsigned long long*) address;
	unsigned long long old = *address_as_i, assumed;
	do {
		assumed = old;
		old = ::atomicCAS(address_as_i, assumed,
				__double_as_longlong(::fmaxf(val, __longlong_as_double(assumed))));
	} while (assumed != old);
	return __longlong_as_double(old);
}
__device__ static float atomicMin64(double* address, double val)
{
	unsigned long long* address_as_i = (unsigned long long*) address;
	unsigned long long old = *address_as_i, assumed;
	do {
		assumed = old;
		old = ::atomicCAS(address_as_i, assumed,
				__double_as_longlong(::fminf(val, __longlong_as_double(assumed))));
	} while (assumed != old);
	return __longlong_as_double(old);
}

__global__ void posvis_init_krnl32(
		struct par_t *dpar,
		struct pos_t **pos,
		float4 *ijminmax_overall,
		float3 *oa,
		float3 *usrc,
		int *outbndarr,
		int c,
		int start,
		int src,
		int size,
		int set,
		int src_override) {

	/* nfrm_alloc-threaded */
	int f = blockIdx.x * blockDim.x + threadIdx.x + start;

	if (f < size) {
		if (f == start) {
			posvis_streams_outbnd = 0;
			pvst_smooth = dpar->pos_smooth;
			if (src_override)	pvst_smooth = 0;
			dbg_cntr=0;
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

__global__ void posvis_init_krnl64(
		struct par_t *dpar,
		struct pos_t **pos,
		double4 *ijminmax_overall,
		double3 *oa,
		double3 *usrc,
		int *outbndarr,
		int c,
		int start,
		int src,
		int size,
		int set,
		int src_override) {

	/* nfrm_alloc-threaded */
	int f = blockIdx.x * blockDim.x + threadIdx.x + start;

	if (f < size) {
		if (f == start) {
			posvis_streams_outbnd = 0;
			pvst_smooth = dpar->pos_smooth;
			if (src_override)	pvst_smooth = 0;
		}
		ijminmax_overall[f].w = ijminmax_overall[f].y = HUGENUMBER;
		ijminmax_overall[f].x = ijminmax_overall[f].z = -HUGENUMBER;
		pos[f]->posbnd_logfactor = 0.0;

		dev_mtrnsps2(oa, pos[f]->ae, f);
		if (src) {
			/* We're viewing the model from the sun: at the center of each pixel
			 * in the projected view, we want cos(incidence angle), distance from
			 * the COM towards the sun, and the facet number.                */
			dev_mmmul2(oa, pos[f]->se, oa, f); /* oa takes ast into sun coords           */
		} else {
			/* We're viewing the model from Earth: at the center of each POS pixel
			 * we want cos(scattering angle), distance from the COM towards Earth,
			 * and the facet number.  For bistatic situations (lightcurves) we also
									 want cos(incidence angle) and the unit vector towards the source.     */
			dev_mmmul2(oa, pos[f]->oe, oa, f); /* oa takes ast into obs coords */
			if (pos[f]->bistatic) {
				usrc[f].x = usrc[f].y = 0.0; /* unit vector towards source */
				usrc[f].z = 1.0;
				dev_cotrans1(&usrc[f], pos[f]->se, usrc[f], -1);
				dev_cotrans1(&usrc[f], pos[f]->oe, usrc[f], 1); /* in observer coordinates */
			}
		}
		outbndarr[f] = 0;
	}
}

__global__ void posvis_facet_krnl32(
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
		int *outbndarr,
		int set) {
	/* (nf * nframes)-threaded kernel.  This version eliminates as much double
	 * math as possible */

	int f = blockIdx.x * blockDim.x + threadIdx.x;
	int pxa, i, i1, i2, j, j1, j2, imin, imax, jmin, jmax, pn, fac;
	float imin_dbl, imax_dbl, jmin_dbl, jmax_dbl, old, s, t, z, den, kmpxl;
	int3 fidx;
	float3 n, v0, v1, v2, tv0, tv1, tv2, x;

	if (f < nfacets) {
		pn = pos[frm]->n;
		kmpxl = __double2float_rn(pos[frm]->km_per_pixel);
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
		n.x = __double2float_rn(verts[0]->f[f].n[0]);
		n.y = __double2float_rn(verts[0]->f[f].n[1]);
		n.z = __double2float_rn(verts[0]->f[f].n[2]);

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

			v0.x += orbit_offs.x;	v0.y += orbit_offs.x;	v0.z += orbit_offs.x;
			v1.x += orbit_offs.y;	v1.y += orbit_offs.y;	v1.z += orbit_offs.y;
			v2.x += orbit_offs.z;	v2.y += orbit_offs.z;	v2.z += orbit_offs.z;

			/* Find rectangular region (in POS pixels) containing the projected
			 * facet - use floats in case model has illegal parameters and the
			 * pixel numbers exceed the limits for valid integers                         */
			imin_dbl = floor(MIN(v0.x,MIN(v1.x,v2.x)) / kmpxl - SMALLVAL + 0.5);
			imax_dbl = floor(MAX(v0.x,MAX(v1.x,v2.x)) / kmpxl + SMALLVAL + 0.5);
			jmin_dbl = floor(MIN(v0.y,MIN(v1.y,v2.y)) / kmpxl - SMALLVAL + 0.5);
			jmax_dbl = floor(MAX(v0.y,MAX(v1.y,v2.y)) / kmpxl + SMALLVAL + 0.5);

			imin = (imin_dbl < INT_MIN) ? INT_MIN : (int) imin_dbl;
			imax = (imax_dbl > INT_MAX) ? INT_MAX : (int) imax_dbl;
			jmin = (jmin_dbl < INT_MIN) ? INT_MIN : (int) jmin_dbl;
			jmax = (jmax_dbl > INT_MAX) ? INT_MAX : (int) jmax_dbl;

			/*  Set the outbnd flag if the facet extends beyond the POS window  */
			if ((imin < (-pn)) || (imax > pn) || (jmin < (-pn))	|| (jmax > pn)) {
				posvis_streams_outbnd = 1;
				atomicExch(&outbndarr[frm], 1);
			}

			/* Figure out if facet projects at least partly within POS window;
			 * if it does, look at each "contained" POS pixel and get the
			 * z-coordinate and cos(scattering angle)           */
			i1 = MAX(imin, -pn);		j1 = MAX(jmin, -pn);
			i2 = MIN(imax,  pn);		j2 = MIN(jmax,  pn);

			if (i1 > pn || i2 < -pn || j1 > pn || j2 < -pn) {

				/* Facet is entirely outside the POS frame: just keep track of
				 * changed POS region     */
				dev_POSrect_gpu32(pos, src, imin_dbl, imax_dbl, jmin_dbl, jmax_dbl,
						ijminmax_overall, frm);

			} else {

				dev_POSrect_gpu32(pos, src, (float)i1, (float)i2, (float)j1,
						(float)j2, ijminmax_overall, frm);

				/* Facet is at least partly within POS frame: find all POS
				 * pixels whose centers project onto this facet  */
				for (i = i1; i <= i2; i++) {
					x.x = i * kmpxl;
					for (j = j1; j <= j2; j++) {

						x.y = j * kmpxl;
						if (src)	fac = pos[frm]->fill[i][j];
						if (!src)	fac = pos[frm]->f[i][j];

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
							    - (x.x- v0.x) * (v1.y- v0.y)) * den;atomicAdd(&dbg_cntr, 1);
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

								if (src)	old = atomicMaxf(&pos[frm]->zill_s[pxa], z);
								else		old = atomicMaxf(&pos[frm]->z_s[pxa], z);

								if (old < z ) {//|| fac < 0) {

									/* Next line assigns distance of POS pixel
									 * center from COM towards Earth; that is,
									 * by changing zz,it changes pos->z or
									 * pos->zill                */
									/* following line is a first time z calc
									 * for this pixel  */
//									if (fac < 0){
//										atomicAdd(&dbg_cntr, 1);
//										if (src)	atomicExch(&pos[frm]->zill_s[pxa], z);
//										else 		atomicExch(&pos[frm]->z_s[pxa], z);
//									}
//
//									if (pvst_smooth) {
//										/* Assign temp. normal components as float3 */
//										tv0.x = __double2float_rn(verts[0]->v[fidx.x].n[0]);
//										tv0.y = __double2float_rn(verts[0]->v[fidx.x].n[1]);
//										tv0.z = __double2float_rn(verts[0]->v[fidx.x].n[2]);
//										tv1.x = __double2float_rn(verts[0]->v[fidx.y].n[0]);
//										tv1.y = __double2float_rn(verts[0]->v[fidx.y].n[1]);
//										tv1.z = __double2float_rn(verts[0]->v[fidx.y].n[2]);
//										tv2.x = __double2float_rn(verts[0]->v[fidx.z].n[0]);
//										tv2.y = __double2float_rn(verts[0]->v[fidx.z].n[1]);
//										tv2.z = __double2float_rn(verts[0]->v[fidx.z].n[2]);
//
//										/* Get pvs_smoothed version of facet unit
//										 * normal: Take the linear combination
//										 * of the three vertex normals; trans-
//										 * form from body to observer coordina-
//										 * tes; and make sure that it points
//										 * somewhat in our direction.         */
//										n.x = tv0.x + s * (tv1.x - tv0.x) + t * (tv2.x - tv1.x);
//										n.y = tv0.y + s * (tv1.y - tv0.y) + t * (tv2.y - tv1.y);
//										n.z = tv0.z + s * (tv1.z - tv0.z) + t * (tv2.z - tv1.z);
//										dev_cotrans8(&n, oa, n, 1, frm);
//										dev_normalize2(&n);
//									}

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
											float temp = dev_dot_f3(n,usrc[frm]);
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

__global__ void posvis_facet_krnl32mod(
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
		int *outbndarr,
		int set) {
	/* (nf * nframes)-threaded kernel.  This version eliminates as much double
	 * math as possible */

	int f = blockIdx.x * blockDim.x + threadIdx.x;
	int pxa, i, i1, i2, j, j1, j2, imin, imax, jmin, jmax, pn, span;
	float imin_dbl, imax_dbl, jmin_dbl, jmax_dbl, old, s, t, z, den, kmpxl;
	int3 fidx;
	float3 n, n1n0, v0, v1, v2, tv0, tv1, tv2, x;

	if (f < nfacets) {
		pn = pos[frm]->n;
		span = 2*pn+1;
		kmpxl = __double2float_rn(pos[frm]->km_per_pixel);
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
		n.x = __double2float_rn(verts[0]->f[f].n[0]);
		n.y = __double2float_rn(verts[0]->f[f].n[1]);
		n.z = __double2float_rn(verts[0]->f[f].n[2]);

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

			v0.x += orbit_offs.x;	v0.y += orbit_offs.x;	v0.z += orbit_offs.x;
			v1.x += orbit_offs.y;	v1.y += orbit_offs.y;	v1.z += orbit_offs.y;
			v2.x += orbit_offs.z;	v2.y += orbit_offs.z;	v2.z += orbit_offs.z;

			/* Find rectangular region (in POS pixels) containing the projected
			 * facet - use floats in case model has illegal parameters and the
			 * pixel numbers exceed the limits for valid integers                         */
			imin_dbl = floor(MIN(v0.x,MIN(v1.x,v2.x)) / kmpxl - SMALLVAL + 0.5);
			imax_dbl = floor(MAX(v0.x,MAX(v1.x,v2.x)) / kmpxl + SMALLVAL + 0.5);
			jmin_dbl = floor(MIN(v0.y,MIN(v1.y,v2.y)) / kmpxl - SMALLVAL + 0.5);
			jmax_dbl = floor(MAX(v0.y,MAX(v1.y,v2.y)) / kmpxl + SMALLVAL + 0.5);

			imin = (imin_dbl < INT_MIN) ? INT_MIN : (int) imin_dbl;
			imax = (imax_dbl > INT_MAX) ? INT_MAX : (int) imax_dbl;
			jmin = (jmin_dbl < INT_MIN) ? INT_MIN : (int) jmin_dbl;
			jmax = (jmax_dbl > INT_MAX) ? INT_MAX : (int) jmax_dbl;

			/*  Set the outbnd flag if the facet extends beyond the POS window  */
			if ((imin < (-pn)) || (imax > pn) || (jmin < (-pn))	|| (jmax > pn)) {
				posvis_streams_outbnd = 1;
				atomicExch(&outbndarr[frm], 1);
			}

			/* Figure out if facet projects at least partly within POS window;
			 * if it does, look at each "contained" POS pixel and get the
			 * z-coordinate and cos(scattering angle)           */
			i1 = MAX(imin, -pn);		j1 = MAX(jmin, -pn);
			i2 = MIN(imax,  pn);		j2 = MIN(jmax,  pn);

			if (i1 > pn || i2 < -pn || j1 > pn || j2 < -pn) {

				/* Facet is entirely outside the POS frame: just keep track of
				 * changed POS region     */
				dev_POSrect_gpu32(pos, src, imin_dbl, imax_dbl, jmin_dbl, jmax_dbl,
						ijminmax_overall, frm);

			} else {

//				dev_POSrect_gpu32(pos, src, (float)i1, (float)i2, (float)j1,
//						(float)j2, ijminmax_overall, frm);

				/* Assign vertex normals if smoothing is enabled */
				if (pvst_smooth) {
					/* Assign temp. normal components as float3 */
					tv0.x = verts[0]->v[fidx.x].n[0];	tv0.y = verts[0]->v[fidx.x].n[1];
					tv0.z = verts[0]->v[fidx.x].n[2];	tv1.x = verts[0]->v[fidx.y].n[0];
					tv1.y = verts[0]->v[fidx.y].n[1];	tv1.z = verts[0]->v[fidx.y].n[2];
					tv2.x = verts[0]->v[fidx.z].n[0];	tv2.y = verts[0]->v[fidx.z].n[1];
					tv2.z = verts[0]->v[fidx.z].n[2];
					n1n0.x = tv1.x - tv0.x;	n1n0.y = tv1.y - tv0.y;	n1n0.z = tv1.z - tv0.z;
					tv2.x -= tv1.x;	tv2.y -= tv1.y; tv2.z -= tv1.z;
				}

				/* Precalculate s and t components */
				float a, b, c, d, e, h, ti, tj, si, sj, si0, sj0, ti0, tj0, sz, tz;
				a = i1*kmpxl - v0.x;
				b = v2.y - v1.y;
				c = v2.x - v1.x;
				d = j1*kmpxl - v0.y;
				e = v1.x - v0.x;
				h = v1.y - v0.y;
				den = e*b - c*h;
				ti = -h*kmpxl/den;
				tj = e*kmpxl/den;
				si = b*kmpxl/den;
				sj = -c*kmpxl/den;
				si0 = (a*b - c*d)/den;
				ti0 = (e*d -a*h)/den;
				sz = v1.z - v0.z;
				tz = v2.z - v1.z;

				/* Facet is at least partly within POS frame: find all POS
				 * pixels whose centers project onto this facet  */
				for (i=i1; i<=i2; i++) {

					sj0 = si0;	/* Initialize this loop's base sj0, tj0 */
					tj0 = ti0;

					for (j=j1; j<=j2; j++) {
						/* Calculate the pixel address for 1D arrays */
						pxa = (j+pn) * span + (i+pn);
						s = sj0;
						t = tj0;

						if ((s >= -SMALLVAL) && (s <= 1.0 + SMALLVAL)) {// &&
							if(	(t >= -SMALLVAL) && (t <= s + SMALLVAL))	{

								/* Compute z-coordinate of pixel center: its
								 * distance measured from the origin towards
								 * Earth.    */
								z = v0.z + s*sz + t*tz;
//								if (src)	fac = pos[frm]->fill[i][j];
//								if (!src)	fac = pos[frm]->f[i][j];

								/* Following line replaces the previous if check
								 * for z > zz[i][j]
								 * atomicMaxf returns the value that was sitting
								 * at zzf[pxa] at time of call.  So if that value
								 * matches the z we compared to*/
								if (src)	old = atomicMaxf(&pos[frm]->zill_s[pxa], z);
								else		old = atomicMaxf(&pos[frm]->z_s[pxa], z);

								if (old < z){
									/* Next line assigns distance of POS pixel
									 * center from COM towards Earth; that is,
									 * by changing zz,it changes pos->z or
									 * pos->zill                */
									/* following line is a first time z calc
									 * for this pixel  */
//									if (fac < 0){
									//										atomicAdd(&dbg_cntr, 1);
									//										if (src)
									//											atomicExch((unsigned long long int*)&pos[frm]->zill[i][j], __double_as_longlong(z));
									//
									//										else
									//											atomicExch((unsigned long long int*)&pos[frm]->z[i][j], __double_as_longlong(z));
									//									}
									//
									if (pvst_smooth) {

										/* Get pvs_smoothed version of facet unit
										 * normal: Take the linear combination
										 * of the three vertex normals; trans-
										 * form from body to observer coordina-
										 * tes; and make sure that it points
										 * somewhat in our direction.         */
										n.x = tv0.x + s * n1n0.x + t * tv2.x;
										n.y = tv0.y + s * n1n0.y + t * tv2.y;
										n.z = tv0.z + s * n1n0.z + t * tv2.z;
										dev_cotrans8(&n, oa, n, 1, frm);
										dev_normalize2(&n);
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
											float temp = dev_dot_f3(n,usrc[frm]);
											atomicExch(&pos[frm]->cosi_s[pxa], temp);
											if (pos[frm]->cosi_s[pxa] <= 0.0)
												pos[frm]->cose_s[pxa] = 0.0;
										}
									}
									dev_POSrect_gpu32(pos, src, (float)i, (float)i, (float)j,
											(float)j, ijminmax_overall, frm);
									/* Next lines change pos->body/bodyill,
									 * pos->comp/compill, pos->f/fill          */
									if (src) {
										pos[frm]->bodyill[i][j] = body;
										pos[frm]->compill[i][j] = comp;
										pos[frm]->fill[i][j] = f;
									} else {
										pos[frm]->body[i][j] = body;
										pos[frm]->comp[i][j] = comp;
										atomicExch(&pos[frm]->f[i][j], f);
									}

								} /* end if (no other facet yet blocks this facet from view) */
							} /* end if 0 <= t <= s (facet center is "in" this POS pixel) */
						} /* end if 0 <= s <= 1 */

						sj0 += sj;
						tj0 += tj;
					} /* end j-loop over POS rows */
					/* Modify s and t step-wise for the next i-iteration of the pixel loop */
					si0 += si;
					ti0 += ti;

				} /* end i-loop over POS columns */
			} /* end else of if (i1 > pos->n || i2 < -pos->n || j1 > pos->n || j2 < -pos->n) */
		} /* End if (n[2] > 0.0) */
	} /* end if (f < nf) */
}

__global__ void posvis_facet_krnl64(
		struct pos_t **pos,
		struct vertices_t **verts,
		double4 *ijminmax_overall,
		double3 orbit_offs,
		double3 *oa,
		double3 *usrc,
		int src,
		int body,
		int comp,
		int nfacets,
		int frm,
		int smooth,
		int *outbndarr,
		int set) {

	int f = blockIdx.x * blockDim.x + threadIdx.x;
	int i, i1, i2, j, j1, j2, imin, imax, jmin, jmax, pn, fac, span;
	double imin_dbl, imax_dbl, jmin_dbl, jmax_dbl, old, s, t, z, den, kmpxl;
	int3 fidx;
	double3 n, v0, v1, v2, tv0, tv1, tv2, x, n1n0;

	if (f < nfacets) {
		pn = pos[frm]->n;
		kmpxl = pos[frm]->km_per_pixel;
		span = 2*pn + 1;

/* The following section transfers vertex coordinates from double[3]
		 * storage to float3		 */
		fidx.x = verts[0]->f[f].v[0];
		fidx.y = verts[0]->f[f].v[1];
		fidx.z = verts[0]->f[f].v[2];
		tv0.x = verts[0]->v[fidx.x].x[0];
		tv0.y = verts[0]->v[fidx.x].x[1];
		tv0.z = verts[0]->v[fidx.x].x[2];
		tv1.x = verts[0]->v[fidx.y].x[0];
		tv1.y = verts[0]->v[fidx.y].x[1];
		tv1.z = verts[0]->v[fidx.y].x[2];
		tv2.x = verts[0]->v[fidx.z].x[0];
		tv2.y = verts[0]->v[fidx.z].x[1];
		tv2.z = verts[0]->v[fidx.z].x[2];
		v0.x = v0.y = v0.z = v1.x = v1.y = v1.z = v2.x = v2.y = v2.z = 0.0;

		/* Get the normal to this facet in body-fixed (asteroid) coordinates
		 * and convert it to observer coordinates     */
		n.x = verts[0]->f[f].n[0];
		n.y = verts[0]->f[f].n[1];
		n.z = verts[0]->f[f].n[2];
		dev_cotrans3(&n, oa, n, 1, frm);

		/* Consider this facet further only if its normal points somewhat
		 * towards the observer rather than away         */
		if (n.z > 0.0) {

			/* Convert the three sets of vertex coordinates from body to ob-
			 * server coordinates; orbit_offset is the center-of-mass offset
			 * (in observer coordinates) for this model at this frame's epoch
			 * due to orbital motion, in case the model is half of a binary
			 * system.  */
			dev_cotrans3(&v0, oa, tv0, 1, frm);
			dev_cotrans3(&v1, oa, tv1, 1, frm);
			dev_cotrans3(&v2, oa, tv2, 1, frm);

			v0.x += orbit_offs.x;	v0.y += orbit_offs.x;	v0.z += orbit_offs.x;
			v1.x += orbit_offs.y;	v1.y += orbit_offs.y;	v1.z += orbit_offs.y;
			v2.x += orbit_offs.z;	v2.y += orbit_offs.z;	v2.z += orbit_offs.z;

			/* Find rectangular region (in POS pixels) containing the projected
			 * facet - use floats in case model has illegal parameters and the
			 * pixel numbers exceed the limits for valid integers                         */
			imin_dbl = floor(MIN(v0.x,MIN(v1.x,v2.x)) / kmpxl - SMALLVAL + 0.5);
			imax_dbl = floor(MAX(v0.x,MAX(v1.x,v2.x)) / kmpxl + SMALLVAL + 0.5);
			jmin_dbl = floor(MIN(v0.y,MIN(v1.y,v2.y)) / kmpxl - SMALLVAL + 0.5);
			jmax_dbl = floor(MAX(v0.y,MAX(v1.y,v2.y)) / kmpxl + SMALLVAL + 0.5);

			imin = (imin_dbl < INT_MIN) ? INT_MIN : (int) imin_dbl;
			imax = (imax_dbl > INT_MAX) ? INT_MAX : (int) imax_dbl;
			jmin = (jmin_dbl < INT_MIN) ? INT_MIN : (int) jmin_dbl;
			jmax = (jmax_dbl > INT_MAX) ? INT_MAX : (int) jmax_dbl;

			/*  Set the outbnd flag if the facet extends beyond the POS window  */
			if ((imin < (-pn)) || (imax > pn) || (jmin < (-pn))	|| (jmax > pn)) {
				posvis_streams_outbnd = 1;
				atomicExch(&outbndarr[frm], 1);
			}

			/* Figure out if facet projects at least partly within POS window;
			 * if it does, look at each "contained" POS pixel and get the
			 * z-coordinate and cos(scattering angle)           */
			i1 = MAX(imin, -pn);		j1 = MAX(jmin, -pn);
			i2 = MIN(imax,  pn);		j2 = MIN(jmax,  pn);

			if (i1 > pn || i2 < -pn || j1 > pn || j2 < -pn) {

				/* Facet is entirely outside the POS frame: just keep track of
				 * changed POS region     */
				dev_POSrect_gpu64(pos, src, imin_dbl, imax_dbl, jmin_dbl, jmax_dbl,
						ijminmax_overall, frm);

			} else {

//				dev_POSrect_gpu64(pos, src, (double)i1, (double)i2, (double)j1,
//						(double)j2, ijminmax_overall, frm);

				/* Assign vertex normals if smoothing is enabled */
				if (pvst_smooth) {
					/* Assign temp. normal components as float3 */
					tv0.x = verts[0]->v[fidx.x].n[0];	tv0.y = verts[0]->v[fidx.x].n[1];
					tv0.z = verts[0]->v[fidx.x].n[2];	tv1.x = verts[0]->v[fidx.y].n[0];
					tv1.y = verts[0]->v[fidx.y].n[1];	tv1.z = verts[0]->v[fidx.y].n[2];
					tv2.x = verts[0]->v[fidx.z].n[0];	tv2.y = verts[0]->v[fidx.z].n[1];
					tv2.z = verts[0]->v[fidx.z].n[2];
					n1n0.x = tv1.x - tv0.x;	n1n0.y = tv1.y - tv0.y;	n1n0.z = tv1.z - tv0.z;
					tv2.x -= tv1.x;	tv2.y -= tv1.y; tv2.z -= tv1.z;
				}

				/* Facet is at least partly within POS frame: find all POS
				 * pixels whose centers project onto this facet  */
				for (i = i1; i <= i2; i++) {
					x.x = i * kmpxl;

					for (j = j1; j <= j2; j++) {
						x.y = j * kmpxl;

//						/* Compute parameters s(x,y) and t(x,y) which define a
//						 * facet's surface as
//						 *         z = z0 + s*(z1-z0) + t*(z2-z1)
//						 * where z0, z1, and z2 are the z-coordinates at the
//						 * vertices. The conditions 0 <= s <= 1 and
//						 * 0 <= t <= s require the POS pixel center to be
//						 * "within" the (projected) perimeter of facet f.    */
						den = 1	/ ((v1.x - v0.x) * (v2.y - v1.y)
								 - (v2.x - v1.x) * (v1.y - v0.y));
						s = ((x.x - v0.x) * (v2.y - v1.y)
						  - (v2.x - v1.x) * (x.y - v0.y)) * den;

						if ((s >= -SMALLVAL) && (s <= 1.0 + SMALLVAL)) {

							t = ((v1.x - v0.x) * (x.y - v0.y)
							    - (x.x- v0.x) * (v1.y- v0.y)) * den;
							if ((t >= -SMALLVAL) && (t <= s + SMALLVAL)) {
//								atomicAdd(&dbg_cntr, 1);

								/* Compute z-coordinate of pixel center: its
								 * distance measured from the origin towards
								 * Earth.    */
								z = v0.z + s*(v1.z-v0.z) + t*(v2.z-v1.z);
								if (src)	fac = pos[frm]->fill[i][j];
								if (!src)	fac = pos[frm]->f[i][j];

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

//								if (src)	old = atomicMax64(&pos[frm]->zill[i][j], z);
								//else
								old = atomicMax64(&pos[frm]->z[i][j], z);

								if (old < z){// || fac < 0) {
//									atomicAdd(&dbg_cntr, 1);
									/* Next line assigns distance of POS pixel
									 * center from COM towards Earth; that is,
									 * by changing zz,it changes pos->z or
									 * pos->zill                */
									/* following line is a first time z calc
									 * for this pixel  */
//									if (fac < 0){
//										atomicAdd(&dbg_cntr, 1);
//										if (src)
//											atomicExch((unsigned long long int*)&pos[frm]->zill[i][j], __double_as_longlong(z));
//
//										else
//											atomicExch((unsigned long long int*)&pos[frm]->z[i][j], __double_as_longlong(z));
//									}
//
//									if (pvst_smooth) {
//										/* Assign temp. normal components as float3 */
//										tv0.x = verts[0]->v[fidx.x].n[0];
//										tv0.y = verts[0]->v[fidx.x].n[1];
//										tv0.z = verts[0]->v[fidx.x].n[2];
//										tv1.x = verts[0]->v[fidx.y].n[0];
//										tv1.y = verts[0]->v[fidx.y].n[1];
//										tv1.z = verts[0]->v[fidx.y].n[2];
//										tv2.x = verts[0]->v[fidx.z].n[0];
//										tv2.y = verts[0]->v[fidx.z].n[1];
//										tv2.z = verts[0]->v[fidx.z].n[2];
//
//										/* Get pvs_smoothed version of facet unit
//										 * normal: Take the linear combination
//										 * of the three vertex normals; trans-
//										 * form from body to observer coordina-
//										 * tes; and make sure that it points
//										 * somewhat in our direction.         */
//										n.x = tv0.x + s * (tv1.x - tv0.x) + t * (tv2.x - tv1.x);
//										n.y = tv0.y + s * (tv1.y - tv0.y) + t * (tv2.y - tv1.y);
//										n.z = tv0.z + s * (tv1.z - tv0.z) + t * (tv2.z - tv1.z);
//
//										dev_cotrans3(&n, oa, n, 1, frm);
//										dev_normalize3(&n);
//									}

									/* Determine scattering and/or incidence
									 * angles. Next lines change pos->cose/
									 * cosill. If bistatic (lightcurves), where
									 * we are viewing from Earth (src = 0),
									 * pos->cosi is also changed.                 */
									if (n.z > 0.0) {

										if (src)
											atomicExch((unsigned long long int*)&pos[frm]->cosill[i][j],
													__double_as_longlong(n.z));
										else
											atomicExch((unsigned long long int*)&pos[frm]->cose[i][j],
													__double_as_longlong(n.z));

										if ((!src) && (pos[frm]->bistatic)) {

											double temp = dev_dot_d3(n,usrc[frm]);
											atomicExch((unsigned long long int*)&pos[frm]->cosi[i][j],
													__double_as_longlong(temp));
											if (pos[frm]->cosi[i][j] <= 0.0)
												pos[frm]->cose[i][j] = 0.0;
										}
									}
									dev_POSrect_gpu64(pos, src, (double)i1, (double)i2, (double)j1,
											(double)j2, ijminmax_overall, frm);
									/* Next lines change pos->body/bodyill,
									 * pos->comp/compill, pos->f/fill          */
									if (src) {
										pos[frm]->bodyill[i][j] = body;
										pos[frm]->compill[i][j] = comp;
										pos[frm]->fill[i][j] = f;
									} else {
										pos[frm]->body[i][j] = body;
										pos[frm]->comp[i][j] = comp;
										atomicExch(&pos[frm]->f[i][j], f);
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

__global__ void posvis_facet_krnl64mod(
		struct pos_t **pos,
		struct vertices_t **verts,
		double4 *ijminmax_overall,
		double3 orbit_offs,
		double3 *oa,
		double3 *usrc,
		int src,
		int body,
		int comp,
		int nfacets,
		int frm,
		int smooth,
		int *outbndarr,
		int set) {

	int f = blockIdx.x * blockDim.x + threadIdx.x;
	int i, i1, i2, j, j1, j2, imin, imax, jmin, jmax, pn, span;
	double imin_dbl, imax_dbl, jmin_dbl, jmax_dbl, old, s, t, z, den, kmpxl;
	int3 fidx;
	double3 n, v0, v1, v2, tv0, tv1, tv2, n1n0;

	if (f < nfacets) {
		pn = pos[frm]->n;
		kmpxl = pos[frm]->km_per_pixel;
		span = 2*pn + 1;

/* The following section transfers vertex coordinates from double[3]
		 * storage to float3		 */
		fidx.x = verts[0]->f[f].v[0];
		fidx.y = verts[0]->f[f].v[1];
		fidx.z = verts[0]->f[f].v[2];
		tv0.x = verts[0]->v[fidx.x].x[0];
		tv0.y = verts[0]->v[fidx.x].x[1];
		tv0.z = verts[0]->v[fidx.x].x[2];
		tv1.x = verts[0]->v[fidx.y].x[0];
		tv1.y = verts[0]->v[fidx.y].x[1];
		tv1.z = verts[0]->v[fidx.y].x[2];
		tv2.x = verts[0]->v[fidx.z].x[0];
		tv2.y = verts[0]->v[fidx.z].x[1];
		tv2.z = verts[0]->v[fidx.z].x[2];
		v0.x = v0.y = v0.z = v1.x = v1.y = v1.z = v2.x = v2.y = v2.z = 0.0;

		/* Get the normal to this facet in body-fixed (asteroid) coordinates
		 * and convert it to observer coordinates     */
		n.x = verts[0]->f[f].n[0];
		n.y = verts[0]->f[f].n[1];
		n.z = verts[0]->f[f].n[2];
		dev_cotrans3(&n, oa, n, 1, frm);

		/* Consider this facet further only if its normal points somewhat
		 * towards the observer rather than away         */
		if (n.z > 0.0) {

			/* Convert the three sets of vertex coordinates from body to ob-
			 * server coordinates; orbit_offset is the center-of-mass offset
			 * (in observer coordinates) for this model at this frame's epoch
			 * due to orbital motion, in case the model is half of a binary
			 * system.  */
			dev_cotrans3(&v0, oa, tv0, 1, frm);
			dev_cotrans3(&v1, oa, tv1, 1, frm);
			dev_cotrans3(&v2, oa, tv2, 1, frm);

			v0.x += orbit_offs.x;	v0.y += orbit_offs.x;	v0.z += orbit_offs.x;
			v1.x += orbit_offs.y;	v1.y += orbit_offs.y;	v1.z += orbit_offs.y;
			v2.x += orbit_offs.z;	v2.y += orbit_offs.z;	v2.z += orbit_offs.z;

			/* Find rectangular region (in POS pixels) containing the projected
			 * facet - use floats in case model has illegal parameters and the
			 * pixel numbers exceed the limits for valid integers                         */
			imin_dbl = floor(MIN(v0.x,MIN(v1.x,v2.x)) / kmpxl - SMALLVAL + 0.5);
			imax_dbl = floor(MAX(v0.x,MAX(v1.x,v2.x)) / kmpxl + SMALLVAL + 0.5);
			jmin_dbl = floor(MIN(v0.y,MIN(v1.y,v2.y)) / kmpxl - SMALLVAL + 0.5);
			jmax_dbl = floor(MAX(v0.y,MAX(v1.y,v2.y)) / kmpxl + SMALLVAL + 0.5);

			imin = (imin_dbl < INT_MIN) ? INT_MIN : (int) imin_dbl;
			imax = (imax_dbl > INT_MAX) ? INT_MAX : (int) imax_dbl;
			jmin = (jmin_dbl < INT_MIN) ? INT_MIN : (int) jmin_dbl;
			jmax = (jmax_dbl > INT_MAX) ? INT_MAX : (int) jmax_dbl;

			/*  Set the outbnd flag if the facet extends beyond the POS window  */
			if ((imin < (-pn)) || (imax > pn) || (jmin < (-pn))	|| (jmax > pn)) {
				posvis_streams_outbnd = 1;
				atomicExch(&outbndarr[frm], 1);
			}

			/* Figure out if facet projects at least partly within POS window;
			 * if it does, look at each "contained" POS pixel and get the
			 * z-coordinate and cos(scattering angle)           */
			i1 = MAX(imin, -pn);		j1 = MAX(jmin, -pn);
			i2 = MIN(imax,  pn);		j2 = MIN(jmax,  pn);

			if (i1 > pn || i2 < -pn || j1 > pn || j2 < -pn) {

				/* Facet is entirely outside the POS frame: just keep track of
				 * changed POS region     */
				dev_POSrect_gpu64(pos, src, imin_dbl, imax_dbl, jmin_dbl, jmax_dbl,
						ijminmax_overall, frm);

			} else {

//				dev_POSrect_gpu64(pos, src, (double)i1, (double)i2, (double)j1,
//						(double)j2, ijminmax_overall, frm);

				/* Assign vertex normals if smoothing is enabled */
				if (pvst_smooth) {
					/* Assign temp. normal components as float3 */
					tv0.x = verts[0]->v[fidx.x].n[0];	tv0.y = verts[0]->v[fidx.x].n[1];
					tv0.z = verts[0]->v[fidx.x].n[2];	tv1.x = verts[0]->v[fidx.y].n[0];
					tv1.y = verts[0]->v[fidx.y].n[1];	tv1.z = verts[0]->v[fidx.y].n[2];
					tv2.x = verts[0]->v[fidx.z].n[0];	tv2.y = verts[0]->v[fidx.z].n[1];
					tv2.z = verts[0]->v[fidx.z].n[2];
					n1n0.x = tv1.x - tv0.x;	n1n0.y = tv1.y - tv0.y;	n1n0.z = tv1.z - tv0.z;
					tv2.x -= tv1.x;	tv2.y -= tv1.y; tv2.z -= tv1.z;
				}

				/* Precalculate s and t components */
				double a, b, c, d, e, h, ti, tj, si, sj, si0, sj0, ti0, tj0, sz, tz;
				a = i1*kmpxl - v0.x;
				b = v2.y - v1.y;
				c = v2.x - v1.x;
				d = j1*kmpxl - v0.y;
				e = v1.x - v0.x;
				h = v1.y - v0.y;
				den = e*b - c*h;
				ti = -h*kmpxl/den;
				tj = e*kmpxl/den;
				si = b*kmpxl/den;
				sj = -c*kmpxl/den;
				si0 = (a*b - c*d)/den;
				ti0 = (e*d -a*h)/den;
				sz = v1.z - v0.z;
				tz = v2.z - v1.z;

				/* Facet is at least partly within POS frame: find all POS
				 * pixels whose centers project onto this facet  */
				for (i = i1; i <= i2; i++) {

					sj0 = si0;	/* Initialize this loop's base sj0, tj0 */
					tj0 = ti0;

					for (j = j1; j <= j2; j++) {

						s = sj0;
						t = tj0;

						if ((s >= -SMALLVAL) && (s <= 1.0 + SMALLVAL)) {// &&
							if(	(t >= -SMALLVAL) && (t <= s + SMALLVAL))	{

							/* Compute z-coordinate of pixel center: its
							 * distance measured from the origin towards
							 * Earth.    */
							z = v0.z + s*sz + t*tz;
//							if (src)	fac = pos[frm]->fill[i][j];
//							if (!src)	fac = pos[frm]->f[i][j];

							/* Following line replaces the previous if check
							 * for z > zz[i][j]
							 * atomicMaxf returns the value that was sitting
							 * at zzf[pxa] at time of call.  So if that value
							 * matches the z we compared to*/
							if (src)	old = atomicMax64(&pos[frm]->zill[i][j], z);
							if (!src)	old = atomicMax64(&pos[frm]->z[i][j], z);
							if (old < z){
								/* Next line assigns distance of POS pixel
								 * center from COM towards Earth; that is,
								 * by changing zz,it changes pos->z or
								 * pos->zill                */
								/* following line is a first time z calc
								 * for this pixel  */
//									if (fac < 0){
//										atomicAdd(&dbg_cntr, 1);
//										if (src)
//											atomicExch((unsigned long long int*)&pos[frm]->zill[i][j], __double_as_longlong(z));
//
//										else
//											atomicExch((unsigned long long int*)&pos[frm]->z[i][j], __double_as_longlong(z));
//									}
//
								if (pvst_smooth) {

									/* Get pvs_smoothed version of facet unit
									 * normal: Take the linear combination
									 * of the three vertex normals; trans-
									 * form from body to observer coordina-
									 * tes; and make sure that it points
									 * somewhat in our direction.         */
									n.x = tv0.x + s * n1n0.x + t * tv2.x;
									n.y = tv0.y + s * n1n0.y + t * tv2.y;
									n.z = tv0.z + s * n1n0.z + t * tv2.z;
									dev_cotrans3(&n, oa, n, 1, frm);
									dev_normalize3(&n);
								}

								/* Determine scattering and/or incidence
								 * angles. Next lines change pos->cose/
								 * cosill. If bistatic (lightcurves), where
								 * we are viewing from Earth (src = 0),
								 * pos->cosi is also changed.                 */
								if (n.z > 0.0) {
									if (src)
										atomicExch((unsigned long long int*)&pos[frm]->cosill[i][j],
												__double_as_longlong(n.z));
									else
										atomicExch((unsigned long long int*)&pos[frm]->cose[i][j],
												__double_as_longlong(n.z));

									if ((!src) && (pos[frm]->bistatic)) {

										double temp = dev_dot_d3(n,usrc[frm]);
										atomicExch((unsigned long long int*)&pos[frm]->cosi[i][j],
												__double_as_longlong(temp));
										if (pos[frm]->cosi[i][j] <= 0.0)
											pos[frm]->cose[i][j] = 0.0;
									}
								}
								dev_POSrect_gpu64(pos, src, (double)i, (double)i, (double)j,
										(double)j, ijminmax_overall, frm);
								/* Next lines change pos->body/bodyill,
								 * pos->comp/compill, pos->f/fill          */
								if (src) {
									pos[frm]->bodyill[i][j] = body;
									pos[frm]->compill[i][j] = comp;
									pos[frm]->fill[i][j] = f;
								} else {
									pos[frm]->body[i][j] = body;
									pos[frm]->comp[i][j] = comp;
									atomicExch(&pos[frm]->f[i][j], f);
								}

							} /* end if (no other facet yet blocks this facet from view) */
							} /* end if 0 <= t <= s (facet center is "in" this POS pixel) */
						} /* end if 0 <= s <= 1 */

						sj0 += sj;
						tj0 += tj;
					} /* end j-loop over POS rows */
					/* Modify s and t step-wise for the next i-iteration of the pixel loop */
					si0 += si;
					ti0 += ti;

				} /* end i-loop over POS columns */
			} /* end else of if (i1 > pos->n || i2 < -pos->n || j1 > pos->n || j2 < -pos->n) */
		} /* End if (n[2] > 0.0) */
	} /* end if (f < nf) */
}

__global__ void posvis_facet_krnl64modb(
		struct pos_t **pos,
		struct vertices_t **verts,
		double4 *ijminmax_overall_gm,
		double3 orbit_offs,
		double3 *oa_gm,
		double3 *usrc_gm,
		int src,
		int nfacets,
		int frm,
		int smooth,
		int *outbndarr,
		int set) {

	/* This kernel uses shared memory for oa and usrc and also ijminmax_overall.
	 * However, this doesn't work completely - ijminmax_overall will get written
	 * to and there are multiple thread blocks.  Abandonded for now (12/28/2017)
	 */

	int f = blockIdx.x * blockDim.x + threadIdx.x;
	int i, i1, i2, j, j1, j2, imin, imax, jmin, jmax;
	double imin_dbl, imax_dbl, jmin_dbl, jmax_dbl, old, s, t, z, den;
	__shared__ double kmpxl, oa_sh[3][3];
	__shared__ double3 usrc_sh;
	__shared__ double4 ijminmax_overall_sh;
	__shared__ int pn,bistatic;
	int3 fidx;
	double3 n, v0, v1, v2, tv0, tv1, tv2, n1n0;

	if (threadIdx.x==0) {
		pn = pos[frm]->n;
		kmpxl = pos[frm]->km_per_pixel;
		bistatic = pos[frm]->bistatic;

		ijminmax_overall_sh.w = ijminmax_overall_sh.x =
				ijminmax_overall_sh.y = ijminmax_overall_sh.z = 0.0f;

		/* Load oa for this frame into shared memory */
		oa_sh[0][0] = oa_gm[3*frm].x;	oa_sh[0][1] = oa_gm[3*frm].y;	oa_sh[0][2] = oa_gm[3*frm].z;
		oa_sh[1][0] = oa_gm[3*frm+1].x;	oa_sh[1][1] = oa_gm[3*frm+1].y;	oa_sh[1][2] = oa_gm[3*frm+1].z;
		oa_sh[2][0] = oa_gm[3*frm+2].x;	oa_sh[2][1] = oa_gm[3*frm+2].y;	oa_sh[2][2] = oa_gm[3*frm+2].z;

		/* Load usrc for this frame into shared memory if it's a bistatic situation */
		if (bistatic) {
			usrc_sh.x =	usrc_gm[frm].x; usrc_sh.y = usrc_gm[frm].y;	usrc_sh.z = usrc_gm[frm].z;
		}
	}
	__syncthreads();

	if (f < nfacets) {

		/* The following section transfers vertex coordinates from double[3]
		 * storage to float3		 */
		fidx.x = verts[0]->f[f].v[0];
		fidx.y = verts[0]->f[f].v[1];
		fidx.z = verts[0]->f[f].v[2];
		tv0.x = verts[0]->v[fidx.x].x[0];
		tv0.y = verts[0]->v[fidx.x].x[1];
		tv0.z = verts[0]->v[fidx.x].x[2];
		tv1.x = verts[0]->v[fidx.y].x[0];
		tv1.y = verts[0]->v[fidx.y].x[1];
		tv1.z = verts[0]->v[fidx.y].x[2];
		tv2.x = verts[0]->v[fidx.z].x[0];
		tv2.y = verts[0]->v[fidx.z].x[1];
		tv2.z = verts[0]->v[fidx.z].x[2];
		v0.x = v0.y = v0.z = v1.x = v1.y = v1.z = v2.x = v2.y = v2.z = 0.0;

		/* Get the normal to this facet in body-fixed (asteroid) coordinates
		 * and convert it to observer coordinates     */
		n.x = verts[0]->f[f].n[0];
		n.y = verts[0]->f[f].n[1];
		n.z = verts[0]->f[f].n[2];
		dev_cotrans1(&n, oa_sh, n, 1);

		/* Consider this facet further only if its normal points somewhat
		 * towards the observer rather than away         */
		if (n.z > 0.0) {

			/* Convert the three sets of vertex coordinates from body to ob-
			 * server coordinates; orbit_offset is the center-of-mass offset
			 * (in observer coordinates) for this model at this frame's epoch
			 * due to orbital motion, in case the model is half of a binary
			 * system.  */
			dev_cotrans1(&v0, oa_sh, tv0, 1);
			dev_cotrans1(&v1, oa_sh, tv1, 1);
			dev_cotrans1(&v2, oa_sh, tv2, 1);

			v0.x += orbit_offs.x;	v0.y += orbit_offs.x;	v0.z += orbit_offs.x;
			v1.x += orbit_offs.y;	v1.y += orbit_offs.y;	v1.z += orbit_offs.y;
			v2.x += orbit_offs.z;	v2.y += orbit_offs.z;	v2.z += orbit_offs.z;

			/* Find rectangular region (in POS pixels) containing the projected
			 * facet - use floats in case model has illegal parameters and the
			 * pixel numbers exceed the limits for valid integers                         */
			imin_dbl = floor(MIN(v0.x,MIN(v1.x,v2.x)) / kmpxl - SMALLVAL + 0.5);
			imax_dbl = floor(MAX(v0.x,MAX(v1.x,v2.x)) / kmpxl + SMALLVAL + 0.5);
			jmin_dbl = floor(MIN(v0.y,MIN(v1.y,v2.y)) / kmpxl - SMALLVAL + 0.5);
			jmax_dbl = floor(MAX(v0.y,MAX(v1.y,v2.y)) / kmpxl + SMALLVAL + 0.5);

			imin = (imin_dbl < INT_MIN) ? INT_MIN : (int) imin_dbl;
			imax = (imax_dbl > INT_MAX) ? INT_MAX : (int) imax_dbl;
			jmin = (jmin_dbl < INT_MIN) ? INT_MIN : (int) jmin_dbl;
			jmax = (jmax_dbl > INT_MAX) ? INT_MAX : (int) jmax_dbl;

			/*  Set the outbnd flag if the facet extends beyond the POS window  */
			if ((imin < (-pn)) || (imax > pn) || (jmin < (-pn))	|| (jmax > pn)) {
				posvis_streams_outbnd = 1;
				atomicExch(&outbndarr[frm], 1);
			}

			/* Figure out if facet projects at least partly within POS window;
			 * if it does, look at each "contained" POS pixel and get the
			 * z-coordinate and cos(scattering angle)           */
			i1 = MAX(imin, -pn);		j1 = MAX(jmin, -pn);
			i2 = MIN(imax,  pn);		j2 = MIN(jmax,  pn);

			if (i1 > pn || i2 < -pn || j1 > pn || j2 < -pn) {

				/* Facet is entirely outside the POS frame: just keep track of
				 * changed POS region     */
				dev_POSrect_gpu64_shared(imin_dbl, imax_dbl, jmin_dbl, jmax_dbl,
						&ijminmax_overall_sh, pn);
			} else {

				dev_POSrect_gpu64_shared((double)i1, (double)i2, (double)j1,
						(double)j2, &ijminmax_overall_sh, pn);

				/* Assign vertex normals if smoothing is enabled */
				if (pvst_smooth) {
					/* Assign temp. normal components as float3 */
					tv0.x = verts[0]->v[fidx.x].n[0];	tv0.y = verts[0]->v[fidx.x].n[1];
					tv0.z = verts[0]->v[fidx.x].n[2];	tv1.x = verts[0]->v[fidx.y].n[0];
					tv1.y = verts[0]->v[fidx.y].n[1];	tv1.z = verts[0]->v[fidx.y].n[2];
					tv2.x = verts[0]->v[fidx.z].n[0];	tv2.y = verts[0]->v[fidx.z].n[1];
					tv2.z = verts[0]->v[fidx.z].n[2];
					n1n0.x = tv1.x - tv0.x;	n1n0.y = tv1.y - tv0.y;	n1n0.z = tv1.z - tv0.z;
					tv2.x -= tv1.x;	tv2.y -= tv1.y; tv2.z -= tv1.z;
				}

				/* Precalculate s and t components */
				double a, b, c, d, e, h, ti, tj, si, sj, si0, sj0, ti0, tj0, sz, tz;
				a = i1*kmpxl - v0.x;
				b = v2.y - v1.y;
				c = v2.x - v1.x;
				d = j1*kmpxl - v0.y;
				e = v1.x - v0.x;
				h = v1.y - v0.y;
				den = e*b - c*h;
				ti = -h*kmpxl/den;
				tj = e*kmpxl/den;
				si = b*kmpxl/den;
				sj = -c*kmpxl/den;
				si0 = (a*b - c*d)/den;
				ti0 = (e*d -a*h)/den;
				sz = v1.z - v0.z;
				tz = v2.z - v1.z;

				/* Facet is at least partly within POS frame: find all POS
				 * pixels whose centers project onto this facet  */
				for (i = i1; i <= i2; i++) {

					sj0 = si0;	/* Initialize this loop's base sj0, tj0 */
					tj0 = ti0;

					for (j = j1; j <= j2; j++) {

						s = sj0;
						t = tj0;

						if ((s >= -SMALLVAL) && (s <= 1.0 + SMALLVAL)) {// &&
							if(	(t >= -SMALLVAL) && (t <= s + SMALLVAL))	{

								/* Compute z-coordinate of pixel center: its
								 * distance measured from the origin towards
								 * Earth.    */
								z = v0.z + s*sz + t*tz;

								/* Following line replaces the previous if check
								 * for z > zz[i][j]
								 * atomicMaxf returns the value that was sitting
								 * at zzf[pxa] at time of call.  So if that value
								 * matches the z we compared to*/
								if (src)	old = atomicMax64(&pos[frm]->zill[i][j], z);
								if (!src)	old = atomicMax64(&pos[frm]->z[i][j], z);

								if (old < z){

									if (pvst_smooth) {

										/* Get pvs_smoothed version of facet unit
										 * normal: Take the linear combination
										 * of the three vertex normals; trans-
										 * form from body to observer coordina-
										 * tes; and make sure that it points
										 * somewhat in our direction.         */
										n.x = tv0.x + s * n1n0.x + t * tv2.x;
										n.y = tv0.y + s * n1n0.y + t * tv2.y;
										n.z = tv0.z + s * n1n0.z + t * tv2.z;
										dev_cotrans1(&n, oa_sh, n, 1);
										dev_normalize3(&n);
									}

									/* Determine scattering and/or incidence
									 * angles. Next lines change pos->cose/
									 * cosill. If bistatic (lightcurves), where
									 * we are viewing from Earth (src = 0),
									 * pos->cosi is also changed.                 */
									if (n.z > 0.0) {

										if (src)
											atomicExch((unsigned long long int*)&pos[frm]->cosill[i][j],
													__double_as_longlong(n.z));
										else
											atomicExch((unsigned long long int*)&pos[frm]->cose[i][j],
													__double_as_longlong(n.z));

										if (!src && bistatic) {

											double temp = dev_dot_d3(n,usrc_sh);
											atomicExch((unsigned long long int*)&pos[frm]->cosi[i][j],
													__double_as_longlong(temp));
											if (pos[frm]->cosi[i][j] <= 0.0)
												pos[frm]->cose[i][j] = 0.0;
										}
									} else {
										if (src)
											atomicExch((unsigned long long int*)&pos[frm]->cosill[i][j],
													__double_as_longlong(0.0));
										else
											atomicExch((unsigned long long int*)&pos[frm]->cose[i][j],
													__double_as_longlong(0.0));

										if (!src && bistatic)
											atomicExch((unsigned long long int*)&pos[frm]->cosi[i][j],
													__double_as_longlong(0.0));
									}
									if (src) {
//										pos[frm]->bodyill[i][j] = body;
//										pos[frm]->compill[i][j] = comp;
										pos[frm]->fill[i][j] = f;
									} else {
//										pos[frm]->body[i][j] = body;
//										pos[frm]->comp[i][j] = comp;
										pos[frm]->f[i][j] = f;
									}
								} /* end if (no other facet yet blocks this facet from view) */
							} /* end if 0 <= t <= s (facet center is "in" this POS pixel) */
						} /* end if 0 <= s <= 1 */

						sj0 += sj;
						tj0 += tj;
					} /* end j-loop over POS rows */
					/* Modify s and t step-wise for the next i-iteration of the pixel loop */
					si0 += si;
					ti0 += ti;

				} /* end i-loop over POS columns */
			} /* end else of if (i1 > pos->n || i2 < -pos->n || j1 > pos->n || j2 < -pos->n) */
		} /* End if (n[2] > 0.0) */
	} /* end if (f < nf) */

	__syncthreads();

	/* Now write the POS frame window limits from shared mem back to global mem */
	if (threadIdx.x==0) {
//		ijminmax_overall_gm[frm].w = ijminmax_overall_sh.w;
//		ijminmax_overall_gm[frm].x = ijminmax_overall_sh.x;
//		ijminmax_overall_gm[frm].y = ijminmax_overall_sh.y;
//		ijminmax_overall_gm[frm].z = ijminmax_overall_sh.z;

		/* Do atomic min/max here because we have multiple blocks with shared mem */
		atomicMin64(&ijminmax_overall_gm[frm].w, ijminmax_overall_sh.w);
		atomicMax64(&ijminmax_overall_gm[frm].x, ijminmax_overall_sh.x);
		atomicMin64(&ijminmax_overall_gm[frm].y, ijminmax_overall_sh.y);
		atomicMax64(&ijminmax_overall_gm[frm].z, ijminmax_overall_sh.z);

//		/*  Update the subset of the POS frame that contains the target  */
//		/* imin_dbl - ijminmax_overall[frm].w
//		 * imax_dbl - ijminmax_overall[frm].x
//		 * jmin_dbl - ijminmax_overall[frm].y
//		 * jmax_dbl - ijminmax_overall[frm].z
//		 */
//		int imin = (ijminmax_overall_gm.w < INT_MIN) ? INT_MIN : (int) ijminmax_overall_gm[frm].w;
//		int imax = (ijminmax_overall_gm.x > INT_MAX) ? INT_MAX : (int) ijminmax_overall_sh.x;
//		int jmin = (ijminmax_overall_gm.y < INT_MIN) ? INT_MIN : (int) ijminmax_overall_sh.y;
//		int jmax = (ijminmax_overall_gm.z > INT_MAX) ? INT_MAX : (int) ijminmax_overall_sh.z;
//
//		/* Make sure it's smaller than n */
//		imin = MAX(imin,-pn);
//		imax = MIN(imax, pn);
//		jmin = MAX(jmin,-pn);
//		jmax = MIN(jmax, pn);
//
//		if (src) {
//			atomicMin(&pos[frm]->xlim2[0], imin);
//			atomicMax(&pos[frm]->xlim2[1], imax);
//			atomicMin(&pos[frm]->ylim2[0], jmin);
//			atomicMax(&pos[frm]->ylim2[1], jmax);
//		} else {
//			atomicMin(&pos[frm]->xlim[0], imin);
//			atomicMax(&pos[frm]->xlim[1], imax);
//			atomicMin(&pos[frm]->ylim[0], jmin);
//			atomicMax(&pos[frm]->ylim[1], jmax);
//		}
	}
}

__global__ void posvis_facet_krnl64af(
		struct pos_t **pos,
		struct vertices_t **verts,
		double4 *ijminmax_overall_gm,
		double3 orbit_offs,
		double3 *oa_gm,
		double3 *usrc_gm,
		int src,
		int nfacets,
		int start,
		int smooth,
		int *outbndarr,
		int set,
		int blocks) {

	/* This kernel uses shared memory for oa and usrc and also ijminmax_overall.
	 * However, this doesn't work completely - ijminmax_overall will get written
	 * to and there are multiple thread blocks.
	 * This particular version processes all frames in a set with one kernel
	 * call.  Each frame gets a configurable number of blocks (via input argument
	 * 'blocks') within which to perform a grid-stride loop over model facets.
	 * Larger models with higher facet count/vertex count should use a larger
	 * number of blocks
	 */

	/* First, calculate which frame this is - based on current blockIdx.x and
	 * how many blocks per frame we are using (via argument 'blocks')	 */
	int frm = blockIdx.x/blocks + start;
	int f, i, i1, i2, j, j1, j2, imin, imax, jmin, jmax;
	double imin_dbl, imax_dbl, jmin_dbl, jmax_dbl, old, s, t, z, den;
	__shared__ double kmpxl, oa_sh[3][3];
	__shared__ double3 usrc_sh;
	__shared__ double4 ijminmax_overall_sh;
	__shared__ int pn, bistatic, fr;
	int3 fidx;
	double3 n, v0, v1, v2, tv0, tv1, tv2, n1n0;


	/* Variables used by every thread get loaded into shared memory */
	if (threadIdx.x==0) {
		pn = pos[frm]->n;
		kmpxl = pos[frm]->km_per_pixel;
		bistatic = pos[frm]->bistatic;

		/* Calculate the thread offset depending which block this is in the
		 * blocks/frame unit		 */
		fr = (blockIdx.x % blocks) * blockDim.x;

		/* Reset the shared mem array for keeping track of floating point POS
		 * region		 */
		ijminmax_overall_sh.w = ijminmax_overall_sh.x =
				ijminmax_overall_sh.y = ijminmax_overall_sh.z = 0.0f;

		/* Load oa for this frame into shared memory */
		oa_sh[0][0] = oa_gm[3*frm].x;	oa_sh[0][1] = oa_gm[3*frm].y;	oa_sh[0][2] = oa_gm[3*frm].z;
		oa_sh[1][0] = oa_gm[3*frm+1].x;	oa_sh[1][1] = oa_gm[3*frm+1].y;	oa_sh[1][2] = oa_gm[3*frm+1].z;
		oa_sh[2][0] = oa_gm[3*frm+2].x;	oa_sh[2][1] = oa_gm[3*frm+2].y;	oa_sh[2][2] = oa_gm[3*frm+2].z;

		/* Load usrc for this frame into shared memory if it's a bistatic situation */
		if (bistatic) {
			usrc_sh.x =	usrc_gm[frm].x; usrc_sh.y = usrc_gm[frm].y;	usrc_sh.z = usrc_gm[frm].z;
		}
	}
	__syncthreads();

	for (f=threadIdx.x+fr; f<nfacets; f+=(blocks*blockDim.x)) {

		/* The following section transfers vertex coordinates from double[3]
		 * storage to float3		 */
		fidx.x = verts[0]->f[f].v[0];
		fidx.y = verts[0]->f[f].v[1];
		fidx.z = verts[0]->f[f].v[2];
		tv0.x = verts[0]->v[fidx.x].x[0];
		tv0.y = verts[0]->v[fidx.x].x[1];
		tv0.z = verts[0]->v[fidx.x].x[2];
		tv1.x = verts[0]->v[fidx.y].x[0];
		tv1.y = verts[0]->v[fidx.y].x[1];
		tv1.z = verts[0]->v[fidx.y].x[2];
		tv2.x = verts[0]->v[fidx.z].x[0];
		tv2.y = verts[0]->v[fidx.z].x[1];
		tv2.z = verts[0]->v[fidx.z].x[2];
		v0.x = v0.y = v0.z = v1.x = v1.y = v1.z = v2.x = v2.y = v2.z = 0.0;

		/* Get the normal to this facet in body-fixed (asteroid) coordinates
		 * and convert it to observer coordinates     */
		n.x = verts[0]->f[f].n[0];
		n.y = verts[0]->f[f].n[1];
		n.z = verts[0]->f[f].n[2];
		dev_cotrans1(&n, oa_sh, n, 1);

		/* Consider this facet further only if its normal points somewhat
		 * towards the observer rather than away         */
		if (n.z > 0.0) {

			/* Convert the three sets of vertex coordinates from body to ob-
			 * server coordinates; orbit_offset is the center-of-mass offset
			 * (in observer coordinates) for this model at this frame's epoch
			 * due to orbital motion, in case the model is half of a binary
			 * system.  */
			dev_cotrans1(&v0, oa_sh, tv0, 1);
			dev_cotrans1(&v1, oa_sh, tv1, 1);
			dev_cotrans1(&v2, oa_sh, tv2, 1);

			v0.x += orbit_offs.x;	v0.y += orbit_offs.x;	v0.z += orbit_offs.x;
			v1.x += orbit_offs.y;	v1.y += orbit_offs.y;	v1.z += orbit_offs.y;
			v2.x += orbit_offs.z;	v2.y += orbit_offs.z;	v2.z += orbit_offs.z;

			/* Find rectangular region (in POS pixels) containing the projected
			 * facet - use floats in case model has illegal parameters and the
			 * pixel numbers exceed the limits for valid integers                         */
			imin_dbl = floor(MIN(v0.x,MIN(v1.x,v2.x)) / kmpxl - SMALLVAL + 0.5);
			imax_dbl = floor(MAX(v0.x,MAX(v1.x,v2.x)) / kmpxl + SMALLVAL + 0.5);
			jmin_dbl = floor(MIN(v0.y,MIN(v1.y,v2.y)) / kmpxl - SMALLVAL + 0.5);
			jmax_dbl = floor(MAX(v0.y,MAX(v1.y,v2.y)) / kmpxl + SMALLVAL + 0.5);

			imin = (imin_dbl < INT_MIN) ? INT_MIN : (int) imin_dbl;
			imax = (imax_dbl > INT_MAX) ? INT_MAX : (int) imax_dbl;
			jmin = (jmin_dbl < INT_MIN) ? INT_MIN : (int) jmin_dbl;
			jmax = (jmax_dbl > INT_MAX) ? INT_MAX : (int) jmax_dbl;

			/*  Set the outbnd flag if the facet extends beyond the POS window  */
			if ((imin < (-pn)) || (imax > pn) || (jmin < (-pn))	|| (jmax > pn)) {
				posvis_streams_outbnd = 1;
				atomicExch(&outbndarr[frm], 1);
			}

			/* Figure out if facet projects at least partly within POS window;
			 * if it does, look at each "contained" POS pixel and get the
			 * z-coordinate and cos(scattering angle)           */
			i1 = MAX(imin, -pn);		j1 = MAX(jmin, -pn);
			i2 = MIN(imax,  pn);		j2 = MIN(jmax,  pn);

			if (i1 > pn || i2 < -pn || j1 > pn || j2 < -pn) {

				/* Facet is entirely outside the POS frame: just keep track of
				 * changed POS region     */
				dev_POSrect_gpu64_shared(imin_dbl, imax_dbl, jmin_dbl, jmax_dbl,
						&ijminmax_overall_sh, pn);
			} else {

				dev_POSrect_gpu64_shared((double)i1, (double)i2, (double)j1,
						(double)j2, &ijminmax_overall_sh, pn);

				/* Assign vertex normals if smoothing is enabled */
				if (pvst_smooth) {
					/* Assign temp. normal components as float3 */
					tv0.x = verts[0]->v[fidx.x].n[0];	tv0.y = verts[0]->v[fidx.x].n[1];
					tv0.z = verts[0]->v[fidx.x].n[2];	tv1.x = verts[0]->v[fidx.y].n[0];
					tv1.y = verts[0]->v[fidx.y].n[1];	tv1.z = verts[0]->v[fidx.y].n[2];
					tv2.x = verts[0]->v[fidx.z].n[0];	tv2.y = verts[0]->v[fidx.z].n[1];
					tv2.z = verts[0]->v[fidx.z].n[2];
					n1n0.x = tv1.x - tv0.x;	n1n0.y = tv1.y - tv0.y;	n1n0.z = tv1.z - tv0.z;
					tv2.x -= tv1.x;	tv2.y -= tv1.y; tv2.z -= tv1.z;
				}

				/* Precalculate s and t components */
				double a, b, c, d, e, h, ti, tj, si, sj, si0, sj0, ti0, tj0, sz, tz;
				a = i1*kmpxl - v0.x;
				b = v2.y - v1.y;
				c = v2.x - v1.x;
				d = j1*kmpxl - v0.y;
				e = v1.x - v0.x;
				h = v1.y - v0.y;
				den = e*b - c*h;
				ti = -h*kmpxl/den;
				tj = e*kmpxl/den;
				si = b*kmpxl/den;
				sj = -c*kmpxl/den;
				si0 = (a*b - c*d)/den;
				ti0 = (e*d -a*h)/den;
				sz = v1.z - v0.z;
				tz = v2.z - v1.z;

				/* Facet is at least partly within POS frame: find all POS
				 * pixels whose centers project onto this facet  */
				for (i = i1; i <= i2; i++) {

					sj0 = si0;	/* Initialize this loop's base sj0, tj0 */
					tj0 = ti0;

					for (j = j1; j <= j2; j++) {

						s = sj0;
						t = tj0;

						if ((s >= -SMALLVAL) && (s <= 1.0 + SMALLVAL)) {// &&
							if(	(t >= -SMALLVAL) && (t <= s + SMALLVAL))	{

								/* Compute z-coordinate of pixel center: its
								 * distance measured from the origin towards
								 * Earth.    */
								z = v0.z + s*sz + t*tz;

								/* Following line replaces the previous if check
								 * for z > zz[i][j]
								 * atomicMaxf returns the value that was sitting
								 * at zzf[pxa] at time of call.  So if that value
								 * matches the z we compared to*/
								if (src)	old = atomicMax64(&pos[frm]->zill[i][j], z);
								if (!src)	old = atomicMax64(&pos[frm]->z[i][j], z);
								if (old < z){

									if (pvst_smooth) {

										/* Get pvs_smoothed version of facet unit
										 * normal: Take the linear combination
										 * of the three vertex normals; trans-
										 * form from body to observer coordina-
										 * tes; and make sure that it points
										 * somewhat in our direction.         */
										n.x = tv0.x + s * n1n0.x + t * tv2.x;
										n.y = tv0.y + s * n1n0.y + t * tv2.y;
										n.z = tv0.z + s * n1n0.z + t * tv2.z;
										dev_cotrans1(&n, oa_sh, n, 1);
										dev_normalize3(&n);
									}

									/* Determine scattering and/or incidence
									 * angles. Next lines change pos->cose/
									 * cosill. If bistatic (lightcurves), where
									 * we are viewing from Earth (src = 0),
									 * pos->cosi is also changed.                 */
									if (n.z > 0.0) {
										if (src)
											atomicExch((unsigned long long int*)&pos[frm]->cosill[i][j],
													__double_as_longlong(n.z));
										else
											atomicExch((unsigned long long int*)&pos[frm]->cose[i][j],
													__double_as_longlong(n.z));

										if (!src && bistatic) {

											double temp = dev_dot_d3(n,usrc_sh);
											atomicExch((unsigned long long int*)&pos[frm]->cosi[i][j],
													__double_as_longlong(temp));
											if (pos[frm]->cosi[i][j] <= 0.0)
												pos[frm]->cose[i][j] = 0.0;
										}
									}

								} /* end if (no other facet yet blocks this facet from view) */
							} /* end if 0 <= t <= s (facet center is "in" this POS pixel) */
						} /* end if 0 <= s <= 1 */

						sj0 += sj;
						tj0 += tj;
					} /* end j-loop over POS rows */
					/* Modify s and t step-wise for the next i-iteration of the pixel loop */
					si0 += si;
					ti0 += ti;

				} /* end i-loop over POS columns */
			} /* end else of if (i1 > pos->n || i2 < -pos->n || j1 > pos->n || j2 < -pos->n) */
		} /* End if (n[2] > 0.0) */
	} /* end if (f < nf) */

	__syncthreads();

	/* Now write the POS frame window limits from shared mem back to global mem.
	 * Note that each threadblock has its own shared memory and it cannot be
	 * accessed by other blocks.  However, every block will perform the atomic
	 * operations below with its own shared memory against the global memory
	 * array ijminmax_overall_gm, accessible by every thread block.  This solves
	 * the access issue and should hopefully result in the right results. */
	if (threadIdx.x==0) {
		/* Do atomic min/max here because we have multiple blocks with shared mem */
		atomicMin64(&ijminmax_overall_gm[frm].w, ijminmax_overall_sh.w);
		atomicMax64(&ijminmax_overall_gm[frm].x, ijminmax_overall_sh.x);
		atomicMin64(&ijminmax_overall_gm[frm].y, ijminmax_overall_sh.y);
		atomicMax64(&ijminmax_overall_gm[frm].z, ijminmax_overall_sh.z);
	}
}

__global__ void posvis_outbnd_krnl32(struct pos_t **pos,
		int *outbndarr, float4 *ijminmax_overall, int size, int start) {
	/* nfrm_alloc-threaded kernel */
	int posn, f = blockIdx.x * blockDim.x + threadIdx.x + start;
	double xfactor, yfactor;
	if (f <size) {

		//		printf("dbg_cntr in posvis_gpu32 = %i\n", dbg_cntr);

		if (outbndarr[f]) {
			/* ijminmax_overall.w = imin_overall
			 * ijminmax_overall.x = imax_overall
			 * ijminmax_overall.y = jmin_overall
			 * ijminmax_overall.z = jmax_overall	 */
			posn = pos[f]->n;
			xfactor = (MAX( ijminmax_overall[f].x,  posn) -
					MIN( ijminmax_overall[f].w, -posn) + 1) / (2*posn+1);
			yfactor = (MAX( ijminmax_overall[f].z,  posn) -
					MIN( ijminmax_overall[f].y, -posn) + 1) / (2*posn+1);
			pos[f]->posbnd_logfactor = log(xfactor*yfactor);
		}
	}
}

__global__ void posvis_outbnd_krnl64(struct pos_t **pos,
		int *outbndarr, double4 *ijminmax_overall, int size, int start) {
	/* nfrm_alloc-threaded kernel */
	int posn, f = blockIdx.x * blockDim.x + threadIdx.x + start;
	double xfactor, yfactor;

	if (f <size) {

//		printf("dbg_cntr in posvis_gpu64 = %i\n", dbg_cntr);

		if (outbndarr[f]) {
			/* ijminmax_overall.w = imin_overall
			 * ijminmax_overall.x = imax_overall
			 * ijminmax_overall.y = jmin_overall
			 * ijminmax_overall.z = jmax_overall	 */
			posn = pos[f]->n;
			xfactor = (MAX( ijminmax_overall[f].x,  posn) -
					MIN( ijminmax_overall[f].w, -posn) + 1) / (2*posn+1);
			yfactor = (MAX( ijminmax_overall[f].z,  posn) -
					MIN( ijminmax_overall[f].y, -posn) + 1) / (2*posn+1);
			pos[f]->posbnd_logfactor = log(xfactor*yfactor);
		}
	}
}

__global__ void posvis_outbnd_krnl64modb(struct pos_t **pos, int *outbndarr,
		double4 *ijminmax_overall, int size, int start, int src) {
	/* nfrm_alloc-threaded kernel */
	int posn, f = blockIdx.x * blockDim.x + threadIdx.x + start;
	double xfactor, yfactor;
	int pn, imin, imax, jmin, jmax;

	if (f <size) {

		/* First calculate each frame's pos xlim and ylim values that define
		 * the pos window containing the model asteroid 		 */
		pn = pos[f]->n;
		imin = (ijminmax_overall[f].w < INT_MIN) ? INT_MIN : (int) ijminmax_overall[f].w;
		imax = (ijminmax_overall[f].x > INT_MAX) ? INT_MAX : (int) ijminmax_overall[f].x;
		jmin = (ijminmax_overall[f].y < INT_MIN) ? INT_MIN : (int) ijminmax_overall[f].y;
		jmax = (ijminmax_overall[f].z > INT_MAX) ? INT_MAX : (int) ijminmax_overall[f].z;

		/* Make sure it's smaller than n */
		imin = MAX(imin,-pn);
		imax = MIN(imax, pn);
		jmin = MAX(jmin,-pn);
		jmax = MIN(jmax, pn);

		if (src) {
			atomicMin(&pos[f]->xlim2[0], imin);
			atomicMax(&pos[f]->xlim2[1], imax);
			atomicMin(&pos[f]->ylim2[0], jmin);
			atomicMax(&pos[f]->ylim2[1], jmax);
		} else {
			atomicMin(&pos[f]->xlim[0], imin);
			atomicMax(&pos[f]->xlim[1], imax);
			atomicMin(&pos[f]->ylim[0], jmin);
			atomicMax(&pos[f]->ylim[1], jmax);
		}

		/* Now take care of out of bounds business */

		if (outbndarr[f]) {
			/* ijminmax_overall.w = imin_overall
			 * ijminmax_overall.x = imax_overall
			 * ijminmax_overall.y = jmin_overall
			 * ijminmax_overall.z = jmax_overall	 */
			posn = pos[f]->n;
			xfactor = (MAX( ijminmax_overall[f].x,  posn) -
					MIN( ijminmax_overall[f].w, -posn) + 1) / (2*posn+1);
			yfactor = (MAX( ijminmax_overall[f].z,  posn) -
					MIN( ijminmax_overall[f].y, -posn) + 1) / (2*posn+1);
			pos[f]->posbnd_logfactor = log(xfactor*yfactor);
		}
	}
}

__host__ int posvis_gpu32(
		struct par_t *dpar,
		struct mod_t *dmod,
		struct dat_t *ddat,
		struct pos_t **pos,
		struct vertices_t **verts,
		float3 orbit_offset,
		int *posn,
		int *outbndarr,
		int set,
		int nfrm_alloc,
		int src,
		int nf,
		int body, int comp, unsigned char type, cudaStream_t *pv_stream,
		int src_override) {

	int f, outbnd, smooth, start;
	dim3 BLK,THD, BLKfrm, THD64;
	float4 *ijminmax_overall;
	float3 *oa, *usrc;

	/* Launch parameters for the facet_streams kernel */
	THD.x = maxThreadsPerBlock;	THD64.x = 64;
	BLK.x = floor((THD.x - 1 + nf) / THD.x);
	BLKfrm.x = floor((THD64.x - 1 + nfrm_alloc)/THD64.x);

	/* Set up the offset addressing for lightcurves if this is a lightcurve */
	if (type == LGHTCRV)	start = 1;	/* fixes the lightcurve offsets */
	else 					start = 0;

	int oasize = nfrm_alloc*3;
	/* Allocate temporary arrays/structs */
	gpuErrchk(cudaMalloc((void**)&ijminmax_overall, sizeof(float4) * nfrm_alloc));
	gpuErrchk(cudaMalloc((void**)&oa, sizeof(float3) * oasize));
	gpuErrchk(cudaMalloc((void**)&usrc, sizeof(float3) * nfrm_alloc));

	posvis_init_krnl32<<<BLKfrm,THD64>>>(dpar, pos, ijminmax_overall, oa, usrc,
			outbndarr, comp, start, src, nfrm_alloc, set, src_override);
	checkErrorAfterKernelLaunch("posvis_init_krnl32");

	for (f=start; f<nfrm_alloc; f++) {
		/* Now the main facet kernel */
		posvis_facet_krnl32mod<<<BLK,THD, 0, pv_stream[f-start]>>>(pos, verts,
				ijminmax_overall, orbit_offset, oa, usrc,	src, body, comp,
				nf, f, smooth, outbndarr, set);
	}
	checkErrorAfterKernelLaunch("posvis_facet_krnl32");

	/* Take care of any posbnd flags */
	posvis_outbnd_krnl32<<<BLKfrm,THD64>>>(pos,
			outbndarr, ijminmax_overall, nfrm_alloc, start);
	checkErrorAfterKernelLaunch("posvis_outbnd_krnl32");
	gpuErrchk(cudaMemcpyFromSymbol(&outbnd, posvis_streams_outbnd, sizeof(int), 0,
			cudaMemcpyDeviceToHost));

//	int n = 200;
//	int npixels = 401*401;
//	f = 0;
//	dbg_print_pos_arrays_full(pos, f, npixels, n);

	/* Free temp arrays, destroy streams and timers, as applicable */

	cudaFree(ijminmax_overall);
	cudaFree(oa);
	cudaFree(usrc);

	return outbnd;
}

__host__ int posvis_gpu64(
		struct par_t *dpar,
		struct mod_t *dmod,
		struct dat_t *ddat,
		struct pos_t **pos,
		struct vertices_t **verts,
		double3 orbit_offset,
		int *posn,
		int *outbndarr,
		int set,
		int nfrm_alloc,
		int src,
		int nf,
		int body, int comp, unsigned char type, cudaStream_t *pv_stream,
		int src_override) {

	int f, outbnd, smooth, start, blocks;
	dim3 BLK,THD, BLKfrm, THD64, BLKaf, THDaf;
	double4 *ijminmax_overall;
	double3 *oa, *usrc;

	/* Launch parameters for the facet_streams kernel */
	THD.x = maxThreadsPerBlock;	THD64.x = 64;
	BLK.x = floor((THD.x - 1 + nf) / THD.x);
	BLKfrm.x = floor((THD64.x - 1 + nfrm_alloc)/THD64.x);
	THDaf.x = 384;
	blocks = 12;


	/* Set up the offset addressing for lightcurves if this is a lightcurve */
	if (type == LGHTCRV)	start = 1;	/* fixes the lightcurve offsets */
	else 					start = 0;
	BLKaf.x = blocks * (nfrm_alloc-start);
	int oasize = nfrm_alloc*3;
	/* Allocate temporary arrays/structs */
	gpuErrchk(cudaMalloc((void**)&ijminmax_overall, sizeof(double4) * nfrm_alloc));
	gpuErrchk(cudaMalloc((void**)&oa, sizeof(double3) * oasize));
//	gpuErrchk(cudaMalloc((void**)&usrc, sizeof(double3) * nfrm_alloc));
	cudaCalloc((void**)&usrc, sizeof(double3), nfrm_alloc);

	posvis_init_krnl64<<<BLKfrm,THD64>>>(dpar, pos, ijminmax_overall, oa, usrc,
			outbndarr, comp, start, src, nfrm_alloc, set, src_override);
	checkErrorAfterKernelLaunch("posvis_init_krnl64");


	for (f=start; f<nfrm_alloc; f++) {
		/* Now the main facet kernel */
//		posvis_facet_krnl64mod<<<BLK,THD, 0, pv_stream[f-start]>>>(pos, verts,
//				ijminmax_overall, orbit_offset, oa, usrc,	src, body, comp,
//				nf, f, smooth, outbndarr, set);
//		posvis_facet_krnl64moda<<<BLK,THD, 0, pv_stream[f-start]>>>(pos, verts,
//				ijminmax_overall, orbit_offset, oa, usrc,src, nf, f, smooth, outbndarr, set);
		posvis_facet_krnl64modb<<<BLK,THD, 0, pv_stream[f-start]>>>(pos, verts,
				ijminmax_overall, orbit_offset, oa, usrc,src, nf, f, smooth, outbndarr, set);
	}
//	posvis_facet_krnl64af<<<BLKaf,THDaf>>>(pos, verts, ijminmax_overall,
//			orbit_offset, oa, usrc, src, nf, start, smooth, outbndarr, set,
//			blocks);
	checkErrorAfterKernelLaunch("posvis_facet_krnl64");

	/* Take care of any posbnd flags */
	posvis_outbnd_krnl64modb<<<BLKfrm,THD64>>>(pos,
			outbndarr, ijminmax_overall, nfrm_alloc, start, src);
	checkErrorAfterKernelLaunch("posvis_outbnd_krnl64");
	gpuErrchk(cudaMemcpyFromSymbol(&outbnd, posvis_streams_outbnd, sizeof(int), 0,
			cudaMemcpyDeviceToHost));

//	int n = 200;
//	int npixels = 401*401;
//	f = 0;
//	dbg_print_pos_arrays_full64(pos, f, npixels, n);

	/* Free temp arrays, destroy streams and timers, as applicable */

	cudaFree(ijminmax_overall);
	cudaFree(oa);
	cudaFree(usrc);

	return outbnd;
}

