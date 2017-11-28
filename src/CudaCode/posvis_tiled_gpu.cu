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
__device__ int posvis_tiled_outbnd, posvis_tiled_smooth;

/* Note that the following custom atomic functions must be declared in each
 * file it is needed (consequence of being a static device function) */
__global__ void posvis_tiled_init_krnl64(
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
			posvis_tiled_outbnd = 0;
			posvis_tiled_smooth = dpar->pos_smooth;
			if (src_override)	posvis_tiled_smooth = 0;
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

__global__ void transform_facet_normals_krnl64(
		struct mod_t *dmod,
		struct pos_t **pos,
		struct vertices_t **verts,
		double4 *ijminmax_overall,
		double3 orbit_offs,
		double3 *oa,
		double3 *usrc,
		int *outbndarr,
		int nf,
		int *nvf,	/* The # of visible facets */
		int frm,
		int src)
{
	/* This kernel launches nf threads and transforms each facet normal with
	 * oa[frm] and stores the result back to dmod if n.z > 0.0.
	 * It also determines and stores the facet and global model bounding box
	 * via i1,i2,j1,j2 and xlim/ylim. These also get stored back to the dmod.
	 * The kernel counts the number of visible facets (where n.z > 0.0)
	 */
	int f = blockIdx.x * blockDim.x + threadIdx.x;	/* Facet index */

	/* Declare kernel variables */
	int pn = pos[frm]->n;
	int imin, jmin, imax, jmax, i1, i2, j1, j2;
	int3 fidx;
	double kmpxl = pos[frm]->km_per_pixel;
	double imin_dbl, jmin_dbl, imax_dbl, jmax_dbl;
	double3 n;
	double3 v0, v1, v2;

	/* Check f is within bounds */
	if (f < nf) {

		pn = pos[frm]->n;
		kmpxl = __double2float_rn(pos[frm]->km_per_pixel);

		/* Get vertex indices of the three vertices making up the facet */
		fidx.x = verts[0]->f[f].v[0]; fidx.y = verts[0]->f[f].v[1];
		fidx.z = verts[0]->f[f].v[2];

		/* Copy each vertex over to thread register memory */
		v0.x = verts[0]->v[fidx.x].x[0];	v0.y = verts[0]->v[fidx.x].x[1];
		v0.z = verts[0]->v[fidx.x].x[2];	v1.x = verts[0]->v[fidx.y].x[0];
		v1.y = verts[0]->v[fidx.y].x[1];	v1.z = verts[0]->v[fidx.y].x[2];
		v2.x = verts[0]->v[fidx.z].x[0];	v2.y = verts[0]->v[fidx.z].x[1];
		v2.z = verts[0]->v[fidx.z].x[2];

		/* Get the normal to this facet in body-fixed (asteroid) coordinates
		 * and convert it to observer coordinates     */
		n.x = verts[0]->f[f].n[0];
		n.y = verts[0]->f[f].n[1];
		n.z = verts[0]->f[f].n[2];
		dev_cotrans3(&n, oa, n, 1, frm);

		/* Check if this facet is visible - is the facet normal pointing
		 * roughly at the observer?				 */
		if (n.z > 0.0) {

			/* First, store the transformed normal back to the model and increase
			 * visible facet counter */
			verts[0]->f[f].nt = n;
			atomicAdd(&nvf[frm], 1);

			/* Convert the 3 vertex coordinates from body to observer
			 * coordinates; orbit_offset is the center-of-mass offset
			 * (in observer coordinates) for this model at this frame's
			 * epoch due to orbital motion, in case the model is half of
			 *  a binary system.  */
			dev_cotrans3(&v0, oa, v0, 1, frm);
			dev_cotrans3(&v1, oa, v1, 1, frm);
			dev_cotrans3(&v2, oa, v2, 1, frm);
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
				posvis_tiled_outbnd = 1;
				atomicExch(&outbndarr[frm], 1);
			}

			/* Figure out if facet projects at least partly within POS window;
			 * if it does, look at each "contained" POS pixel and get the
			 * z-coordinate and cos(scattering angle)           */
			i1 = MAX(imin, -pn);	i2 = MIN(imax, pn);
			j1 = MAX(jmin, -pn);	j2 = MIN(jmax, pn);

			verts[0]->f[f].ilim.x = i1;
			verts[0]->f[f].jlim.x = j1;
			verts[0]->f[f].ilim.y = i2;
			verts[0]->f[f].jlim.y = j2;

			/* Now keep track of the global region */
			if (i1 > pn || i2 < -pn || j1 > pn || j2 < -pn) {

				/* Facet is entirely outside the POS frame: just keep track of
				 * changed POS region     */
				dev_POSrect_gpu64(pos, src, imin_dbl, imax_dbl, jmin_dbl,
						jmax_dbl, ijminmax_overall, frm);

			} else {
				dev_POSrect_gpu64(pos, src, (double)i1, (double)i2,
						(double)j1, (double)j2, ijminmax_overall, frm);
			}
		}

	}
}


__host__ int posvis_tiled_gpu64(
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

	int f, outbnd, smooth, start, overlap, specific_span;
	dim3 BLK,THD, BLKfrm, THD64;
	double4 *ijminmax_overall;
	double3 *oa, *usrc;
	int *nvf, *xspan, *yspan, *n_tiles_x, *n_tiles_y, *n_tiles;
	int oasize = nfrm_alloc*3;

	/* To-Do:  Calculate these spans at program launch from max shared memory
	 * per thread block	 */
	int span_r64 = 55;		/* These four spans are the specific maximum tile */
	int span_r32 = 78;		/* sides depending on FP32/FP64 mode and data     */
	int span_lc64 = 45;		/* type - lightcurves need one more pos array     */
	int span_lc32 = 64;		/* than radar.									  */


	/* Launch parameters for the facet_streams kernel */
	THD.x = maxThreadsPerBlock;	THD64.x = 64;
	BLK.x = floor((THD.x - 1 + nf) / THD.x);
	BLKfrm.x = floor((THD64.x - 1 + nfrm_alloc)/THD64.x);

	/* Set up the offset addressing for lightcurves if this is a lightcurve */
	if (type == LGHTCRV) {
		start = 1;	/* fixes the lightcurve offsets */
		specific_span = span_lc64;
	}
	else {
		start = 0;
		specific_span = span_r64;
	}

	/* Allocate temporary arrays/structs */
	cudaCalloc1((void**)&ijminmax_overall, 	sizeof(double4), nfrm_alloc);
	cudaCalloc1((void**)&oa, 				sizeof(double3), oasize);
	cudaCalloc1((void**)&usrc, 				sizeof(double3), nfrm_alloc);
	cudaCalloc1((void**)&nvf, 				sizeof(int), 	 nfrm_alloc);
	cudaCalloc1((void**)&xspan, 			sizeof(int), 	 nfrm_alloc);
	cudaCalloc1((void**)&yspan, 			sizeof(int), 	 nfrm_alloc);
	cudaCalloc1((void**)&n_tiles, 			sizeof(int), 	 nfrm_alloc);
	cudaCalloc1((void**)&n_tiles_x, 		sizeof(int), 	 nfrm_alloc);
	cudaCalloc1((void**)&n_tiles_y, 		sizeof(int), 	 nfrm_alloc);

	/* Initialize/pre-calculate values for rasterization */
	posvis_tiled_init_krnl64<<<BLKfrm,THD64>>>(dpar, pos, ijminmax_overall, oa, usrc,
			outbndarr, comp, start, src, nfrm_alloc, set, src_override);
	checkErrorAfterKernelLaunch("posvis_tiled_init_krnl64");

	/* Transform facet normals and determine bounding for facets and pos */
	for (f=start; f<nfrm_alloc; f++) {
		/* Now the main facet kernel */
		transform_facet_normals_krnl64<<<BLK,THD,0,pv_stream[1]>>>(dmod, pos,
				verts, ijminmax_overall, orbit_offset, oa, usrc, outbndarr, nf,
				nvf, f, src);
	}
	checkErrorAfterKernelLaunch("posvis_facet_krnl64");
	cudaDeviceSynchronize();
	for (f=start; f<nfrm_alloc; f++) {
		printf("set[%i] frame[%i] xlim[0] = %i\n", set, f, pos[f]->xlim[0]);
		printf("set[%i] frame[%i] xlim[0] = %i\n", set, f, pos[f]->xlim[1]);
		printf("set[%i] frame[%i] ylim[1] = %i\n", set, f, pos[f]->ylim[0]);
		printf("set[%i] frame[%i] ylim[1] = %i\n", set, f, pos[f]->ylim[1]);
		printf("set[%i] frame[%i] # of visible facets = %i\n", set, f, nvf[f]);
	}

	/* Now calculate the tiling parameters to cover the POS view */
	for (f=start; f<nfrm_alloc; f++) {
		xspan[f] = pos[f]->xlim[1] - pos[f]->xlim[0] + 1;
		yspan[f] = pos[f]->ylim[1] - pos[f]->ylim[0] + 1;
		n_tiles_x[f] = (xspan[f]/specific_span) + 1;
		n_tiles_y[f] = (yspan[f]/specific_span) + 1;
		n_tiles[f] = n_tiles_x[f] * n_tiles_y[f];

		printf("xspan[%i] = %i\n", f, xspan[f]);
		printf("yspan[%i] = %i\n", f, yspan[f]);
		printf("n_tiles[%i] = %i\n", f, n_tiles[f]);
		printf("n_tiles_x[%i] = %i\n", f, n_tiles_x[f]);
		printf("n_tiles_y[%i] = %i\n", f, n_tiles_y[f]);

	}


//	posvis_facet_krnl64<<<BLK,THD, 0, pv_stream[1]>>>(pos, verts,
//					ijminmax_overall, orbit_offset, oa, usrc,	src, body, comp,
//					nf, 1, smooth, outbndarr, set);
//	dbg_krnl_psvs<<<1,1>>>();
//	gpuErrchk(cudaMemcpy(dbg_hn, dbg_n, sizeof(double3)*nf, cudaMemcpyDeviceToHost));
//	dbg_print_facet_normals_dbl3(dbg_hn, nf, "FP64_nrmls.csv");

//	/* Take care of any posbnd flags */
//	posvis_outbnd_krnl64<<<BLKfrm,THD64>>>(pos,
//			outbndarr, ijminmax_overall, nfrm_alloc, start);
//	checkErrorAfterKernelLaunch("posvis_outbnd_krnl64");
//	gpuErrchk(cudaMemcpyFromSymbol(&outbnd, posvis_streams_outbnd, sizeof(int), 0,
//			cudaMemcpyDeviceToHost));

//	int n = 200;
//	int npixels = 401*401;
//	f = 0;
//	dbg_print_pos_arrays_full64(pos, f, npixels, n);

	/* Free temp arrays, destroy streams and timers, as applicable */

//	cudaFree(ijminmax_overall);
//	cudaFree(oa);
//	cudaFree(usrc);

	return outbnd;
}



//__device__ static float atomicMaxf(float* address, float val) {
//	int* address_as_i = (int*) address;
//	int old = *address_as_i, assumed;
//	do {
//		assumed = old;
//		old = ::atomicCAS(address_as_i, assumed,
//				__float_as_int(::fmaxf(val, __int_as_float(assumed))));
//	} while (assumed != old);
//	return __int_as_float(old);
//}
//__device__ static float atomicMax64(double* address, double val)
//{
//	unsigned long long* address_as_i = (unsigned long long*) address;
//	unsigned long long old = *address_as_i, assumed;
//	do {
//		assumed = old;
//		old = ::atomicCAS(address_as_i, assumed,
//				__double_as_longlong(::fmaxf(val, __longlong_as_double(assumed))));
//	} while (assumed != old);
//	return __longlong_as_double(old);
//}
//
//__global__ void posvis_init_krnl32(
//		struct par_t *dpar,
//		struct pos_t **pos,
//		float4 *ijminmax_overall,
//		float3 *oa,
//		float3 *usrc,
//		int *outbndarr,
//		int c,
//		int start,
//		int src,
//		int size,
//		int set,
//		int src_override) {
//
//	/* nfrm_alloc-threaded */
//	int f = blockIdx.x * blockDim.x + threadIdx.x + start;
//
//	if (f < size) {
//		if (f == start) {
//			posvis_streams_outbnd = 0;
//			pvst_smooth = dpar->pos_smooth;
//			if (src_override)	pvst_smooth = 0;
//			dbg_cntr=0;
//		}
//		ijminmax_overall[f].w = ijminmax_overall[f].y = HUGENUMBER;
//		ijminmax_overall[f].x = ijminmax_overall[f].z = -HUGENUMBER;
//		pos[f]->posbnd_logfactor = 0.0;
//
//		dev_mtrnsps3(oa, pos[f]->ae, f);
//		if (src) {
//			/* We're viewing the model from the sun: at the center of each pixel
//			 * in the projected view, we want cos(incidence angle), distance from
//			 * the COM towards the sun, and the facet number.                */
//			dev_mmmul3(oa, pos[f]->se, oa, f); /* oa takes ast into sun coords           */
//		} else {
//			/* We're viewing the model from Earth: at the center of each POS pixel
//			 * we want cos(scattering angle), distance from the COM towards Earth,
//			 * and the facet number.  For bistatic situations (lightcurves) we also
//									 want cos(incidence angle) and the unit vector towards the source.     */
//			dev_mmmul3(oa, pos[f]->oe, oa, f); /* oa takes ast into obs coords */
//			if (pos[f]->bistatic) {
//				usrc[f].x = usrc[f].y = 0.0; /* unit vector towards source */
//				usrc[f].z = 1.0;
//				dev_cotrans9(&usrc[f], pos[f]->se, usrc[f], -1);
//				dev_cotrans9(&usrc[f], pos[f]->oe, usrc[f], 1); /* in observer coordinates */
//			}
//		}
//		outbndarr[f] = 0;
//	}
//}
//

//
//__global__ void posvis_facet_krnl32(
//		struct pos_t **pos,
//		struct vertices_t **verts,
//		float4 *ijminmax_overall,
//		float3 orbit_offs,
//		float3 *oa,
//		float3 *usrc,
//		int src,
//		int body,
//		int comp,
//		int nfacets,
//		int frm,
//		int smooth,
//		int *outbndarr,
//		int set) {
//	/* (nf * nframes)-threaded kernel.  This version eliminates as much double
//	 * math as possible */
//
//	int f = blockIdx.x * blockDim.x + threadIdx.x;
//	int pxa, i, i1, i2, j, j1, j2, imin, imax, jmin, jmax, pn, fac;
//	float imin_dbl, imax_dbl, jmin_dbl, jmax_dbl, old, s, t, z, den, kmpxl;
//	int3 fidx;
//	float3 n, v0, v1, v2, tv0, tv1, tv2, x;
//
//	if (f < nfacets) {
//		pn = pos[frm]->n;
//		kmpxl = __double2float_rn(pos[frm]->km_per_pixel);
//		/* The following section transfers vertex coordinates from double[3]
//		 * storage to float3		 */
//		fidx.x = verts[0]->f[f].v[0];
//		fidx.y = verts[0]->f[f].v[1];
//		fidx.z = verts[0]->f[f].v[2];
//		tv0.x = __double2float_rn( verts[0]->v[fidx.x].x[0]);
//		tv0.y = __double2float_rn(verts[0]->v[fidx.x].x[1]);
//		tv0.z = __double2float_rn(verts[0]->v[fidx.x].x[2]);
//		tv1.x = __double2float_rn(verts[0]->v[fidx.y].x[0]);
//		tv1.y = __double2float_rn(verts[0]->v[fidx.y].x[1]);
//		tv1.z = __double2float_rn(verts[0]->v[fidx.y].x[2]);
//		tv2.x = __double2float_rn(verts[0]->v[fidx.z].x[0]);
//		tv2.y = __double2float_rn(verts[0]->v[fidx.z].x[1]);
//		tv2.z = __double2float_rn(verts[0]->v[fidx.z].x[2]);
//		v0.x = v0.y = v0.z = v1.x = v1.y = v1.z = v2.x = v2.y = v2.z = 0.0;
//
//		/* Get the normal to this facet in body-fixed (asteroid) coordinates
//		 * and convert it to observer coordinates     */
//		n.x = __double2float_rn(verts[0]->f[f].n[0]);
//		n.y = __double2float_rn(verts[0]->f[f].n[1]);
//		n.z = __double2float_rn(verts[0]->f[f].n[2]);
//
//		dev_cotrans8(&n, oa, n, 1, frm);
//
//		/* Consider this facet further only if its normal points somewhat
//		 * towards the observer rather than away         */
//		if (n.z > 0.0) {
//			/* Convert the three sets of vertex coordinates from body to ob-
//			 * server coordinates; orbit_offset is the center-of-mass offset
//			 * (in observer coordinates) for this model at this frame's epoch
//			 * due to orbital motion, in case the model is half of a binary
//			 * system.  */
//			dev_cotrans8(&v0, oa, tv0, 1, frm);
//			dev_cotrans8(&v1, oa, tv1, 1, frm);
//			dev_cotrans8(&v2, oa, tv2, 1, frm);
//
//			v0.x += orbit_offs.x;	v0.y += orbit_offs.x;	v0.z += orbit_offs.x;
//			v1.x += orbit_offs.y;	v1.y += orbit_offs.y;	v1.z += orbit_offs.y;
//			v2.x += orbit_offs.z;	v2.y += orbit_offs.z;	v2.z += orbit_offs.z;
//
//			/* Find rectangular region (in POS pixels) containing the projected
//			 * facet - use floats in case model has illegal parameters and the
//			 * pixel numbers exceed the limits for valid integers                         */
//			imin_dbl = floor(MIN(v0.x,MIN(v1.x,v2.x)) / kmpxl - SMALLVAL + 0.5);
//			imax_dbl = floor(MAX(v0.x,MAX(v1.x,v2.x)) / kmpxl + SMALLVAL + 0.5);
//			jmin_dbl = floor(MIN(v0.y,MIN(v1.y,v2.y)) / kmpxl - SMALLVAL + 0.5);
//			jmax_dbl = floor(MAX(v0.y,MAX(v1.y,v2.y)) / kmpxl + SMALLVAL + 0.5);
//
//			imin = (imin_dbl < INT_MIN) ? INT_MIN : (int) imin_dbl;
//			imax = (imax_dbl > INT_MAX) ? INT_MAX : (int) imax_dbl;
//			jmin = (jmin_dbl < INT_MIN) ? INT_MIN : (int) jmin_dbl;
//			jmax = (jmax_dbl > INT_MAX) ? INT_MAX : (int) jmax_dbl;
//
//			/*  Set the outbnd flag if the facet extends beyond the POS window  */
//			if ((imin < (-pn)) || (imax > pn) || (jmin < (-pn))	|| (jmax > pn)) {
//				posvis_streams_outbnd = 1;
//				atomicExch(&outbndarr[frm], 1);
//			}
//
//			/* Figure out if facet projects at least partly within POS window;
//			 * if it does, look at each "contained" POS pixel and get the
//			 * z-coordinate and cos(scattering angle)           */
//			i1 = MAX(imin, -pn);		j1 = MAX(jmin, -pn);
//			i2 = MIN(imax,  pn);		j2 = MIN(jmax,  pn);
//
//			if (i1 > pn || i2 < -pn || j1 > pn || j2 < -pn) {
//
//				/* Facet is entirely outside the POS frame: just keep track of
//				 * changed POS region     */
//				dev_POSrect_gpu32(pos, src, imin_dbl, imax_dbl, jmin_dbl, jmax_dbl,
//						ijminmax_overall, frm);
//
//			} else {
//
//				dev_POSrect_gpu32(pos, src, (float)i1, (float)i2, (float)j1,
//						(float)j2, ijminmax_overall, frm);
//
//				/* Facet is at least partly within POS frame: find all POS
//				 * pixels whose centers project onto this facet  */
//				for (i = i1; i <= i2; i++) {
//					x.x = i * kmpxl;
//					for (j = j1; j <= j2; j++) {
//
//						x.y = j * kmpxl;
//						if (src)	fac = pos[frm]->fill[i][j];
//						if (!src)	fac = pos[frm]->f[i][j];
//
//						/* Calculate the pixel address for 1D arrays */
//						pxa = (j+pn) * (2*pn + 1) + (i+pn);
//
//						/* Compute parameters s(x,y) and t(x,y) which define a
//						 * facet's surface as
//						 *         z = z0 + s*(z1-z0) + t*(z2-z1)
//						 * where z0, z1, and z2 are the z-coordinates at the
//						 * vertices. The conditions 0 <= s <= 1 and
//						 * 0 <= t <= s require the POS pixel center to be
//						 * "within" the (projected) perimeter of facet f.    */
//						den = 1	/ ((v1.x - v0.x) * (v2.y - v1.y)
//								 - (v2.x - v1.x) * (v1.y - v0.y));
//						s = ((x.x - v0.x) * (v2.y - v1.y)
//						  - (v2.x - v1.x) * (x.y - v0.y)) * den;
//
//						if ((s >= -SMALLVAL) && (s <= 1.0 + SMALLVAL)) {
//
//							t = ((v1.x - v0.x) * (x.y - v0.y)
//							    - (x.x- v0.x) * (v1.y- v0.y)) * den;atomicAdd(&dbg_cntr, 1);
//							if ((t >= -SMALLVAL) && (t <= s + SMALLVAL)) {
//								/* Compute z-coordinate of pixel center: its
//								 * distance measured from the origin towards
//								 * Earth.    */
//								z = v0.z + s*(v1.z-v0.z) + t*(v2.z-v1.z);
//
//								/* If fac[i][j] is >= 0, pixel [i][j] was al-
//								 * ready assigned values during a previous call
//								 * to posvis for a different model component.
//								 * If so, override only if the current component
//								 * is blocking our view of (i.e., is closer to
//								 * us than) the previous one.   */
//
//								/* Following line replaces the previous if check
//								 * for z > zz[i][j]
//								 * atomicMaxf returns the value that was sitting
//								 * at zzf[pxa] at time of call.  So if that value
//								 * matches the z we compared to*/
//
//								if (src)	old = atomicMaxf(&pos[frm]->zill_s[pxa], z);
//								else		old = atomicMaxf(&pos[frm]->z_s[pxa], z);
//
//								if (old < z ) {//|| fac < 0) {
//
//									/* Next line assigns distance of POS pixel
//									 * center from COM towards Earth; that is,
//									 * by changing zz,it changes pos->z or
//									 * pos->zill                */
//									/* following line is a first time z calc
//									 * for this pixel  */
////									if (fac < 0){
////										atomicAdd(&dbg_cntr, 1);
////										if (src)	atomicExch(&pos[frm]->zill_s[pxa], z);
////										else 		atomicExch(&pos[frm]->z_s[pxa], z);
////									}
////
////									if (pvst_smooth) {
////										/* Assign temp. normal components as float3 */
////										tv0.x = __double2float_rn(verts[0]->v[fidx.x].n[0]);
////										tv0.y = __double2float_rn(verts[0]->v[fidx.x].n[1]);
////										tv0.z = __double2float_rn(verts[0]->v[fidx.x].n[2]);
////										tv1.x = __double2float_rn(verts[0]->v[fidx.y].n[0]);
////										tv1.y = __double2float_rn(verts[0]->v[fidx.y].n[1]);
////										tv1.z = __double2float_rn(verts[0]->v[fidx.y].n[2]);
////										tv2.x = __double2float_rn(verts[0]->v[fidx.z].n[0]);
////										tv2.y = __double2float_rn(verts[0]->v[fidx.z].n[1]);
////										tv2.z = __double2float_rn(verts[0]->v[fidx.z].n[2]);
////
////										/* Get pvs_smoothed version of facet unit
////										 * normal: Take the linear combination
////										 * of the three vertex normals; trans-
////										 * form from body to observer coordina-
////										 * tes; and make sure that it points
////										 * somewhat in our direction.         */
////										n.x = tv0.x + s * (tv1.x - tv0.x) + t * (tv2.x - tv1.x);
////										n.y = tv0.y + s * (tv1.y - tv0.y) + t * (tv2.y - tv1.y);
////										n.z = tv0.z + s * (tv1.z - tv0.z) + t * (tv2.z - tv1.z);
////										dev_cotrans8(&n, oa, n, 1, frm);
////										dev_normalize2(&n);
////									}
//
//									/* Determine scattering and/or incidence
//									 * angles. Next lines change pos->cose/
//									 * cosill. If bistatic (lightcurves), where
//									 * we are viewing from Earth (src = 0),
//									 * pos->cosi is also changed.                 */
//									if (n.z > 0.0) {
//
//										if (src)
//											atomicExch(&pos[frm]->cosill_s[pxa], n.z);
//										else
//											atomicExch(&pos[frm]->cose_s[pxa], n.z);
//										if ((!src) && (pos[frm]->bistatic)) {
//											float temp = dev_dot_f3(n,usrc[frm]);
//											atomicExch(&pos[frm]->cosi_s[pxa], temp);
//											if (pos[frm]->cosi_s[pxa] <= 0.0)
//												pos[frm]->cose_s[pxa] = 0.0;
//										}
//									}
//
//									/* Next lines change pos->body/bodyill,
//									 * pos->comp/compill, pos->f/fill          */
//									if (src) {
//										pos[frm]->bodyill[i][j] = body;
//										pos[frm]->compill[i][j] = comp;
//										pos[frm]->fill[i][j] = f;
//									} else {
//										pos[frm]->body[i][j] = body;
//										pos[frm]->comp[i][j] = comp;
//										pos[frm]->f[i][j] = f;
//									}
//
//								} /* end if (no other facet yet blocks this facet from view) */
//							} /* end if 0 <= t <= s (facet center is "in" this POS pixel) */
//						} /* end if 0 <= s <= 1 */
//					} /* end j-loop over POS rows */
//				} /* end i-loop over POS columns */
//			} /* end else of if (i1 > pos->n || i2 < -pos->n || j1 > pos->n || j2 < -pos->n) */
//		} /* End if (n[2] > 0.0) */
//	} /* end if (f < nf) */
//}
//
//__global__ void posvis_facet_krnl32mod(
//		struct pos_t **pos,
//		struct vertices_t **verts,
//		float4 *ijminmax_overall,
//		float3 orbit_offs,
//		float3 *oa,
//		float3 *usrc,
//		int src,
//		int body,
//		int comp,
//		int nfacets,
//		int frm,
//		int smooth,
//		int *outbndarr,
//		int set) {
//	/* (nf * nframes)-threaded kernel.  This version eliminates as much double
//	 * math as possible */
//
//	int f = blockIdx.x * blockDim.x + threadIdx.x;
//	int pxa, i, i1, i2, j, j1, j2, imin, imax, jmin, jmax, pn, span;
//	float imin_dbl, imax_dbl, jmin_dbl, jmax_dbl, old, s, t, z, den, kmpxl;
//	int3 fidx;
//	float3 n, n1n0, v0, v1, v2, tv0, tv1, tv2, x;
//
//	if (f < nfacets) {
//		pn = pos[frm]->n;
//		span = 2*pn+1;
//		kmpxl = __double2float_rn(pos[frm]->km_per_pixel);
//		/* The following section transfers vertex coordinates from double[3]
//		 * storage to float3		 */
//		fidx.x = verts[0]->f[f].v[0];
//		fidx.y = verts[0]->f[f].v[1];
//		fidx.z = verts[0]->f[f].v[2];
//		tv0.x = __double2float_rn( verts[0]->v[fidx.x].x[0]);
//		tv0.y = __double2float_rn(verts[0]->v[fidx.x].x[1]);
//		tv0.z = __double2float_rn(verts[0]->v[fidx.x].x[2]);
//		tv1.x = __double2float_rn(verts[0]->v[fidx.y].x[0]);
//		tv1.y = __double2float_rn(verts[0]->v[fidx.y].x[1]);
//		tv1.z = __double2float_rn(verts[0]->v[fidx.y].x[2]);
//		tv2.x = __double2float_rn(verts[0]->v[fidx.z].x[0]);
//		tv2.y = __double2float_rn(verts[0]->v[fidx.z].x[1]);
//		tv2.z = __double2float_rn(verts[0]->v[fidx.z].x[2]);
//		v0.x = v0.y = v0.z = v1.x = v1.y = v1.z = v2.x = v2.y = v2.z = 0.0;
//
//		/* Get the normal to this facet in body-fixed (asteroid) coordinates
//		 * and convert it to observer coordinates     */
//		n.x = __double2float_rn(verts[0]->f[f].n[0]);
//		n.y = __double2float_rn(verts[0]->f[f].n[1]);
//		n.z = __double2float_rn(verts[0]->f[f].n[2]);
//
//		dev_cotrans8(&n, oa, n, 1, frm);
//
//		/* Consider this facet further only if its normal points somewhat
//		 * towards the observer rather than away         */
//		if (n.z > 0.0) {
//			/* Convert the three sets of vertex coordinates from body to ob-
//			 * server coordinates; orbit_offset is the center-of-mass offset
//			 * (in observer coordinates) for this model at this frame's epoch
//			 * due to orbital motion, in case the model is half of a binary
//			 * system.  */
//			dev_cotrans8(&v0, oa, tv0, 1, frm);
//			dev_cotrans8(&v1, oa, tv1, 1, frm);
//			dev_cotrans8(&v2, oa, tv2, 1, frm);
//
//			v0.x += orbit_offs.x;	v0.y += orbit_offs.x;	v0.z += orbit_offs.x;
//			v1.x += orbit_offs.y;	v1.y += orbit_offs.y;	v1.z += orbit_offs.y;
//			v2.x += orbit_offs.z;	v2.y += orbit_offs.z;	v2.z += orbit_offs.z;
//
//			/* Find rectangular region (in POS pixels) containing the projected
//			 * facet - use floats in case model has illegal parameters and the
//			 * pixel numbers exceed the limits for valid integers                         */
//			imin_dbl = floor(MIN(v0.x,MIN(v1.x,v2.x)) / kmpxl - SMALLVAL + 0.5);
//			imax_dbl = floor(MAX(v0.x,MAX(v1.x,v2.x)) / kmpxl + SMALLVAL + 0.5);
//			jmin_dbl = floor(MIN(v0.y,MIN(v1.y,v2.y)) / kmpxl - SMALLVAL + 0.5);
//			jmax_dbl = floor(MAX(v0.y,MAX(v1.y,v2.y)) / kmpxl + SMALLVAL + 0.5);
//
//			imin = (imin_dbl < INT_MIN) ? INT_MIN : (int) imin_dbl;
//			imax = (imax_dbl > INT_MAX) ? INT_MAX : (int) imax_dbl;
//			jmin = (jmin_dbl < INT_MIN) ? INT_MIN : (int) jmin_dbl;
//			jmax = (jmax_dbl > INT_MAX) ? INT_MAX : (int) jmax_dbl;
//
//			/*  Set the outbnd flag if the facet extends beyond the POS window  */
//			if ((imin < (-pn)) || (imax > pn) || (jmin < (-pn))	|| (jmax > pn)) {
//				posvis_streams_outbnd = 1;
//				atomicExch(&outbndarr[frm], 1);
//			}
//
//			/* Figure out if facet projects at least partly within POS window;
//			 * if it does, look at each "contained" POS pixel and get the
//			 * z-coordinate and cos(scattering angle)           */
//			i1 = MAX(imin, -pn);		j1 = MAX(jmin, -pn);
//			i2 = MIN(imax,  pn);		j2 = MIN(jmax,  pn);
//
//			if (i1 > pn || i2 < -pn || j1 > pn || j2 < -pn) {
//
//				/* Facet is entirely outside the POS frame: just keep track of
//				 * changed POS region     */
//				dev_POSrect_gpu32(pos, src, imin_dbl, imax_dbl, jmin_dbl, jmax_dbl,
//						ijminmax_overall, frm);
//
//			} else {
//
////				dev_POSrect_gpu32(pos, src, (float)i1, (float)i2, (float)j1,
////						(float)j2, ijminmax_overall, frm);
//
//				/* Assign vertex normals if smoothing is enabled */
//				if (pvst_smooth) {
//					/* Assign temp. normal components as float3 */
//					tv0.x = verts[0]->v[fidx.x].n[0];	tv0.y = verts[0]->v[fidx.x].n[1];
//					tv0.z = verts[0]->v[fidx.x].n[2];	tv1.x = verts[0]->v[fidx.y].n[0];
//					tv1.y = verts[0]->v[fidx.y].n[1];	tv1.z = verts[0]->v[fidx.y].n[2];
//					tv2.x = verts[0]->v[fidx.z].n[0];	tv2.y = verts[0]->v[fidx.z].n[1];
//					tv2.z = verts[0]->v[fidx.z].n[2];
//					n1n0.x = tv1.x - tv0.x;	n1n0.y = tv1.y - tv0.y;	n1n0.z = tv1.z - tv0.z;
//					tv2.x -= tv1.x;	tv2.y -= tv1.y; tv2.z -= tv1.z;
//				}
//
//				/* Precalculate s and t components */
//				float a, b, c, d, e, h, ti, tj, si, sj, si0, sj0, ti0, tj0, sz, tz;
//				a = i1*kmpxl - v0.x;
//				b = v2.y - v1.y;
//				c = v2.x - v1.x;
//				d = j1*kmpxl - v0.y;
//				e = v1.x - v0.x;
//				h = v1.y - v0.y;
//				den = e*b - c*h;
//				ti = -h*kmpxl/den;
//				tj = e*kmpxl/den;
//				si = b*kmpxl/den;
//				sj = -c*kmpxl/den;
//				si0 = (a*b - c*d)/den;
//				ti0 = (e*d -a*h)/den;
//				sz = v1.z - v0.z;
//				tz = v2.z - v1.z;
//
//				/* Facet is at least partly within POS frame: find all POS
//				 * pixels whose centers project onto this facet  */
//				for (i=i1; i<=i2; i++) {
//
//					sj0 = si0;	/* Initialize this loop's base sj0, tj0 */
//					tj0 = ti0;
//
//					for (j=j1; j<=j2; j++) {
//						/* Calculate the pixel address for 1D arrays */
//						pxa = (j+pn) * span + (i+pn);
//						s = sj0;
//						t = tj0;
//
//						if ((s >= -SMALLVAL) && (s <= 1.0 + SMALLVAL)) {// &&
//							if(	(t >= -SMALLVAL) && (t <= s + SMALLVAL))	{
//
//								/* Compute z-coordinate of pixel center: its
//								 * distance measured from the origin towards
//								 * Earth.    */
//								z = v0.z + s*sz + t*tz;
////								if (src)	fac = pos[frm]->fill[i][j];
////								if (!src)	fac = pos[frm]->f[i][j];
//
//								/* Following line replaces the previous if check
//								 * for z > zz[i][j]
//								 * atomicMaxf returns the value that was sitting
//								 * at zzf[pxa] at time of call.  So if that value
//								 * matches the z we compared to*/
//								if (src)	old = atomicMaxf(&pos[frm]->zill_s[pxa], z);
//								else		old = atomicMaxf(&pos[frm]->z_s[pxa], z);
//
//								if (old < z){
//									/* Next line assigns distance of POS pixel
//									 * center from COM towards Earth; that is,
//									 * by changing zz,it changes pos->z or
//									 * pos->zill                */
//									/* following line is a first time z calc
//									 * for this pixel  */
////									if (fac < 0){
//									//										atomicAdd(&dbg_cntr, 1);
//									//										if (src)
//									//											atomicExch((unsigned long long int*)&pos[frm]->zill[i][j], __double_as_longlong(z));
//									//
//									//										else
//									//											atomicExch((unsigned long long int*)&pos[frm]->z[i][j], __double_as_longlong(z));
//									//									}
//									//
//									if (pvst_smooth) {
//
//										/* Get pvs_smoothed version of facet unit
//										 * normal: Take the linear combination
//										 * of the three vertex normals; trans-
//										 * form from body to observer coordina-
//										 * tes; and make sure that it points
//										 * somewhat in our direction.         */
//										n.x = tv0.x + s * n1n0.x + t * tv2.x;
//										n.y = tv0.y + s * n1n0.y + t * tv2.y;
//										n.z = tv0.z + s * n1n0.z + t * tv2.z;
//										dev_cotrans8(&n, oa, n, 1, frm);
//										dev_normalize2(&n);
//									}
//
//									/* Determine scattering and/or incidence
//									 * angles. Next lines change pos->cose/
//									 * cosill. If bistatic (lightcurves), where
//									 * we are viewing from Earth (src = 0),
//									 * pos->cosi is also changed.                 */
//									if (n.z > 0.0) {
//
//										if (src)
//											atomicExch(&pos[frm]->cosill_s[pxa], n.z);
//										else
//											atomicExch(&pos[frm]->cose_s[pxa], n.z);
//										if ((!src) && (pos[frm]->bistatic)) {
//											float temp = dev_dot_f3(n,usrc[frm]);
//											atomicExch(&pos[frm]->cosi_s[pxa], temp);
//											if (pos[frm]->cosi_s[pxa] <= 0.0)
//												pos[frm]->cose_s[pxa] = 0.0;
//										}
//									}
//									dev_POSrect_gpu32(pos, src, (float)i, (float)i, (float)j,
//											(float)j, ijminmax_overall, frm);
//									/* Next lines change pos->body/bodyill,
//									 * pos->comp/compill, pos->f/fill          */
//									if (src) {
//										pos[frm]->bodyill[i][j] = body;
//										pos[frm]->compill[i][j] = comp;
//										pos[frm]->fill[i][j] = f;
//									} else {
//										pos[frm]->body[i][j] = body;
//										pos[frm]->comp[i][j] = comp;
//										atomicExch(&pos[frm]->f[i][j], f);
//									}
//
//								} /* end if (no other facet yet blocks this facet from view) */
//							} /* end if 0 <= t <= s (facet center is "in" this POS pixel) */
//						} /* end if 0 <= s <= 1 */
//
//						sj0 += sj;
//						tj0 += tj;
//					} /* end j-loop over POS rows */
//					/* Modify s and t step-wise for the next i-iteration of the pixel loop */
//					si0 += si;
//					ti0 += ti;
//
//				} /* end i-loop over POS columns */
//			} /* end else of if (i1 > pos->n || i2 < -pos->n || j1 > pos->n || j2 < -pos->n) */
//		} /* End if (n[2] > 0.0) */
//	} /* end if (f < nf) */
//}
//
//__global__ void posvis_facet_krnl64(
//		struct pos_t **pos,
//		struct vertices_t **verts,
//		double4 *ijminmax_overall,
//		double3 orbit_offs,
//		double3 *oa,
//		double3 *usrc,
//		int src,
//		int body,
//		int comp,
//		int nfacets,
//		int frm,
//		int smooth,
//		int *outbndarr,
//		int set) {
//
//	int f = blockIdx.x * blockDim.x + threadIdx.x;
//	int i, i1, i2, j, j1, j2, imin, imax, jmin, jmax, pn, fac, span;
//	double imin_dbl, imax_dbl, jmin_dbl, jmax_dbl, old, s, t, z, den, kmpxl;
//	int3 fidx;
//	double3 n, v0, v1, v2, tv0, tv1, tv2, x, n1n0;
//
//	if (f < nfacets) {
//		pn = pos[frm]->n;
//		kmpxl = pos[frm]->km_per_pixel;
//		span = 2*pn + 1;
//
///* The following section transfers vertex coordinates from double[3]
//		 * storage to float3		 */
//		fidx.x = verts[0]->f[f].v[0];
//		fidx.y = verts[0]->f[f].v[1];
//		fidx.z = verts[0]->f[f].v[2];
//		tv0.x = verts[0]->v[fidx.x].x[0];
//		tv0.y = verts[0]->v[fidx.x].x[1];
//		tv0.z = verts[0]->v[fidx.x].x[2];
//		tv1.x = verts[0]->v[fidx.y].x[0];
//		tv1.y = verts[0]->v[fidx.y].x[1];
//		tv1.z = verts[0]->v[fidx.y].x[2];
//		tv2.x = verts[0]->v[fidx.z].x[0];
//		tv2.y = verts[0]->v[fidx.z].x[1];
//		tv2.z = verts[0]->v[fidx.z].x[2];
//		v0.x = v0.y = v0.z = v1.x = v1.y = v1.z = v2.x = v2.y = v2.z = 0.0;
//
//		/* Get the normal to this facet in body-fixed (asteroid) coordinates
//		 * and convert it to observer coordinates     */
//		n.x = verts[0]->f[f].n[0];
//		n.y = verts[0]->f[f].n[1];
//		n.z = verts[0]->f[f].n[2];
//		dev_cotrans3(&n, oa, n, 1, frm);
//
//		/* Consider this facet further only if its normal points somewhat
//		 * towards the observer rather than away         */
//		if (n.z > 0.0) {
//
//			/* Convert the three sets of vertex coordinates from body to ob-
//			 * server coordinates; orbit_offset is the center-of-mass offset
//			 * (in observer coordinates) for this model at this frame's epoch
//			 * due to orbital motion, in case the model is half of a binary
//			 * system.  */
//			dev_cotrans3(&v0, oa, tv0, 1, frm);
//			dev_cotrans3(&v1, oa, tv1, 1, frm);
//			dev_cotrans3(&v2, oa, tv2, 1, frm);
//
//			v0.x += orbit_offs.x;	v0.y += orbit_offs.x;	v0.z += orbit_offs.x;
//			v1.x += orbit_offs.y;	v1.y += orbit_offs.y;	v1.z += orbit_offs.y;
//			v2.x += orbit_offs.z;	v2.y += orbit_offs.z;	v2.z += orbit_offs.z;
//
//			/* Find rectangular region (in POS pixels) containing the projected
//			 * facet - use floats in case model has illegal parameters and the
//			 * pixel numbers exceed the limits for valid integers                         */
//			imin_dbl = floor(MIN(v0.x,MIN(v1.x,v2.x)) / kmpxl - SMALLVAL + 0.5);
//			imax_dbl = floor(MAX(v0.x,MAX(v1.x,v2.x)) / kmpxl + SMALLVAL + 0.5);
//			jmin_dbl = floor(MIN(v0.y,MIN(v1.y,v2.y)) / kmpxl - SMALLVAL + 0.5);
//			jmax_dbl = floor(MAX(v0.y,MAX(v1.y,v2.y)) / kmpxl + SMALLVAL + 0.5);
//
//			imin = (imin_dbl < INT_MIN) ? INT_MIN : (int) imin_dbl;
//			imax = (imax_dbl > INT_MAX) ? INT_MAX : (int) imax_dbl;
//			jmin = (jmin_dbl < INT_MIN) ? INT_MIN : (int) jmin_dbl;
//			jmax = (jmax_dbl > INT_MAX) ? INT_MAX : (int) jmax_dbl;
//
//			/*  Set the outbnd flag if the facet extends beyond the POS window  */
//			if ((imin < (-pn)) || (imax > pn) || (jmin < (-pn))	|| (jmax > pn)) {
//				posvis_streams_outbnd = 1;
//				atomicExch(&outbndarr[frm], 1);
//			}
//
//			/* Figure out if facet projects at least partly within POS window;
//			 * if it does, look at each "contained" POS pixel and get the
//			 * z-coordinate and cos(scattering angle)           */
//			i1 = MAX(imin, -pn);		j1 = MAX(jmin, -pn);
//			i2 = MIN(imax,  pn);		j2 = MIN(jmax,  pn);
//
//			if (i1 > pn || i2 < -pn || j1 > pn || j2 < -pn) {
//
//				/* Facet is entirely outside the POS frame: just keep track of
//				 * changed POS region     */
//				dev_POSrect_gpu64(pos, src, imin_dbl, imax_dbl, jmin_dbl, jmax_dbl,
//						ijminmax_overall, frm);
//
//			} else {
//
////				dev_POSrect_gpu64(pos, src, (double)i1, (double)i2, (double)j1,
////						(double)j2, ijminmax_overall, frm);
//
//				/* Assign vertex normals if smoothing is enabled */
//				if (pvst_smooth) {
//					/* Assign temp. normal components as float3 */
//					tv0.x = verts[0]->v[fidx.x].n[0];	tv0.y = verts[0]->v[fidx.x].n[1];
//					tv0.z = verts[0]->v[fidx.x].n[2];	tv1.x = verts[0]->v[fidx.y].n[0];
//					tv1.y = verts[0]->v[fidx.y].n[1];	tv1.z = verts[0]->v[fidx.y].n[2];
//					tv2.x = verts[0]->v[fidx.z].n[0];	tv2.y = verts[0]->v[fidx.z].n[1];
//					tv2.z = verts[0]->v[fidx.z].n[2];
//					n1n0.x = tv1.x - tv0.x;	n1n0.y = tv1.y - tv0.y;	n1n0.z = tv1.z - tv0.z;
//					tv2.x -= tv1.x;	tv2.y -= tv1.y; tv2.z -= tv1.z;
//				}
//
//				/* Facet is at least partly within POS frame: find all POS
//				 * pixels whose centers project onto this facet  */
//				for (i = i1; i <= i2; i++) {
//					x.x = i * kmpxl;
//
//					for (j = j1; j <= j2; j++) {
//						x.y = j * kmpxl;
//
////						/* Compute parameters s(x,y) and t(x,y) which define a
////						 * facet's surface as
////						 *         z = z0 + s*(z1-z0) + t*(z2-z1)
////						 * where z0, z1, and z2 are the z-coordinates at the
////						 * vertices. The conditions 0 <= s <= 1 and
////						 * 0 <= t <= s require the POS pixel center to be
////						 * "within" the (projected) perimeter of facet f.    */
//						den = 1	/ ((v1.x - v0.x) * (v2.y - v1.y)
//								 - (v2.x - v1.x) * (v1.y - v0.y));
//						s = ((x.x - v0.x) * (v2.y - v1.y)
//						  - (v2.x - v1.x) * (x.y - v0.y)) * den;
//
//						if ((s >= -SMALLVAL) && (s <= 1.0 + SMALLVAL)) {
//
//							t = ((v1.x - v0.x) * (x.y - v0.y)
//							    - (x.x- v0.x) * (v1.y- v0.y)) * den;
//							if ((t >= -SMALLVAL) && (t <= s + SMALLVAL)) {
////								atomicAdd(&dbg_cntr, 1);
//
//								/* Compute z-coordinate of pixel center: its
//								 * distance measured from the origin towards
//								 * Earth.    */
//								z = v0.z + s*(v1.z-v0.z) + t*(v2.z-v1.z);
//								if (src)	fac = pos[frm]->fill[i][j];
//								if (!src)	fac = pos[frm]->f[i][j];
//
//								/* If fac[i][j] is >= 0, pixel [i][j] was al-
//								 * ready assigned values during a previous call
//								 * to posvis for a different model component.
//								 * If so, override only if the current component
//								 * is blocking our view of (i.e., is closer to
//								 * us than) the previous one.   */
//
//								/* Following line replaces the previous if check
//								 * for z > zz[i][j]
//								 * atomicMaxf returns the value that was sitting
//								 * at zzf[pxa] at time of call.  So if that value
//								 * matches the z we compared to*/
//
////								if (src)	old = atomicMax64(&pos[frm]->zill[i][j], z);
//								//else
//								old = atomicMax64(&pos[frm]->z[i][j], z);
//
//								if (old < z){// || fac < 0) {
////									atomicAdd(&dbg_cntr, 1);
//									/* Next line assigns distance of POS pixel
//									 * center from COM towards Earth; that is,
//									 * by changing zz,it changes pos->z or
//									 * pos->zill                */
//									/* following line is a first time z calc
//									 * for this pixel  */
////									if (fac < 0){
////										atomicAdd(&dbg_cntr, 1);
////										if (src)
////											atomicExch((unsigned long long int*)&pos[frm]->zill[i][j], __double_as_longlong(z));
////
////										else
////											atomicExch((unsigned long long int*)&pos[frm]->z[i][j], __double_as_longlong(z));
////									}
////
////									if (pvst_smooth) {
////										/* Assign temp. normal components as float3 */
////										tv0.x = verts[0]->v[fidx.x].n[0];
////										tv0.y = verts[0]->v[fidx.x].n[1];
////										tv0.z = verts[0]->v[fidx.x].n[2];
////										tv1.x = verts[0]->v[fidx.y].n[0];
////										tv1.y = verts[0]->v[fidx.y].n[1];
////										tv1.z = verts[0]->v[fidx.y].n[2];
////										tv2.x = verts[0]->v[fidx.z].n[0];
////										tv2.y = verts[0]->v[fidx.z].n[1];
////										tv2.z = verts[0]->v[fidx.z].n[2];
////
////										/* Get pvs_smoothed version of facet unit
////										 * normal: Take the linear combination
////										 * of the three vertex normals; trans-
////										 * form from body to observer coordina-
////										 * tes; and make sure that it points
////										 * somewhat in our direction.         */
////										n.x = tv0.x + s * (tv1.x - tv0.x) + t * (tv2.x - tv1.x);
////										n.y = tv0.y + s * (tv1.y - tv0.y) + t * (tv2.y - tv1.y);
////										n.z = tv0.z + s * (tv1.z - tv0.z) + t * (tv2.z - tv1.z);
////
////										dev_cotrans3(&n, oa, n, 1, frm);
////										dev_normalize3(&n);
////									}
//
//									/* Determine scattering and/or incidence
//									 * angles. Next lines change pos->cose/
//									 * cosill. If bistatic (lightcurves), where
//									 * we are viewing from Earth (src = 0),
//									 * pos->cosi is also changed.                 */
//									if (n.z > 0.0) {
//
//										if (src)
//											atomicExch((unsigned long long int*)&pos[frm]->cosill[i][j],
//													__double_as_longlong(n.z));
//										else
//											atomicExch((unsigned long long int*)&pos[frm]->cose[i][j],
//													__double_as_longlong(n.z));
//
//										if ((!src) && (pos[frm]->bistatic)) {
//
//											double temp = dev_dot_d3(n,usrc[frm]);
//											atomicExch((unsigned long long int*)&pos[frm]->cosi[i][j],
//													__double_as_longlong(temp));
//											if (pos[frm]->cosi[i][j] <= 0.0)
//												pos[frm]->cose[i][j] = 0.0;
//										}
//									}
//									dev_POSrect_gpu64(pos, src, (double)i1, (double)i2, (double)j1,
//															(double)j2, ijminmax_overall, frm);
//									/* Next lines change pos->body/bodyill,
//									 * pos->comp/compill, pos->f/fill          */
//									if (src) {
//										pos[frm]->bodyill[i][j] = body;
//										pos[frm]->compill[i][j] = comp;
//										pos[frm]->fill[i][j] = f;
//									} else {
//										pos[frm]->body[i][j] = body;
//										pos[frm]->comp[i][j] = comp;
//										atomicExch(&pos[frm]->f[i][j], f);
//									}
//
//								} /* end if (no other facet yet blocks this facet from view) */
//							} /* end if 0 <= t <= s (facet center is "in" this POS pixel) */
//						} /* end if 0 <= s <= 1 */
//					} /* end j-loop over POS rows */
//				} /* end i-loop over POS columns */
//			} /* end else of if (i1 > pos->n || i2 < -pos->n || j1 > pos->n || j2 < -pos->n) */
//		} /* End if (n[2] > 0.0) */
//	} /* end if (f < nf) */
//}
//
//__global__ void posvis_facet_krnl64mod(
//		struct pos_t **pos,
//		struct vertices_t **verts,
//		double4 *ijminmax_overall,
//		double3 orbit_offs,
//		double3 *oa,
//		double3 *usrc,
//		int src,
//		int body,
//		int comp,
//		int nfacets,
//		int frm,
//		int smooth,
//		int *outbndarr,
//		int set) {
//
//	int f = blockIdx.x * blockDim.x + threadIdx.x;
//	int i, i1, i2, j, j1, j2, imin, imax, jmin, jmax, pn, span;
//	double imin_dbl, imax_dbl, jmin_dbl, jmax_dbl, old, s, t, z, den, kmpxl;
//	int3 fidx;
//	double3 n, v0, v1, v2, tv0, tv1, tv2, n1n0;
//
//	if (f < nfacets) {
//		pn = pos[frm]->n;
//		kmpxl = pos[frm]->km_per_pixel;
//		span = 2*pn + 1;
//
///* The following section transfers vertex coordinates from double[3]
//		 * storage to float3		 */
//		fidx.x = verts[0]->f[f].v[0];
//		fidx.y = verts[0]->f[f].v[1];
//		fidx.z = verts[0]->f[f].v[2];
//		tv0.x = verts[0]->v[fidx.x].x[0];
//		tv0.y = verts[0]->v[fidx.x].x[1];
//		tv0.z = verts[0]->v[fidx.x].x[2];
//		tv1.x = verts[0]->v[fidx.y].x[0];
//		tv1.y = verts[0]->v[fidx.y].x[1];
//		tv1.z = verts[0]->v[fidx.y].x[2];
//		tv2.x = verts[0]->v[fidx.z].x[0];
//		tv2.y = verts[0]->v[fidx.z].x[1];
//		tv2.z = verts[0]->v[fidx.z].x[2];
//		v0.x = v0.y = v0.z = v1.x = v1.y = v1.z = v2.x = v2.y = v2.z = 0.0;
//
//		/* Get the normal to this facet in body-fixed (asteroid) coordinates
//		 * and convert it to observer coordinates     */
//		n.x = verts[0]->f[f].n[0];
//		n.y = verts[0]->f[f].n[1];
//		n.z = verts[0]->f[f].n[2];
//		dev_cotrans3(&n, oa, n, 1, frm);
//
//		/* Consider this facet further only if its normal points somewhat
//		 * towards the observer rather than away         */
//		if (n.z > 0.0) {
//
//			/* Convert the three sets of vertex coordinates from body to ob-
//			 * server coordinates; orbit_offset is the center-of-mass offset
//			 * (in observer coordinates) for this model at this frame's epoch
//			 * due to orbital motion, in case the model is half of a binary
//			 * system.  */
//			dev_cotrans3(&v0, oa, tv0, 1, frm);
//			dev_cotrans3(&v1, oa, tv1, 1, frm);
//			dev_cotrans3(&v2, oa, tv2, 1, frm);
//
//			v0.x += orbit_offs.x;	v0.y += orbit_offs.x;	v0.z += orbit_offs.x;
//			v1.x += orbit_offs.y;	v1.y += orbit_offs.y;	v1.z += orbit_offs.y;
//			v2.x += orbit_offs.z;	v2.y += orbit_offs.z;	v2.z += orbit_offs.z;
//
//			/* Find rectangular region (in POS pixels) containing the projected
//			 * facet - use floats in case model has illegal parameters and the
//			 * pixel numbers exceed the limits for valid integers                         */
//			imin_dbl = floor(MIN(v0.x,MIN(v1.x,v2.x)) / kmpxl - SMALLVAL + 0.5);
//			imax_dbl = floor(MAX(v0.x,MAX(v1.x,v2.x)) / kmpxl + SMALLVAL + 0.5);
//			jmin_dbl = floor(MIN(v0.y,MIN(v1.y,v2.y)) / kmpxl - SMALLVAL + 0.5);
//			jmax_dbl = floor(MAX(v0.y,MAX(v1.y,v2.y)) / kmpxl + SMALLVAL + 0.5);
//
//			imin = (imin_dbl < INT_MIN) ? INT_MIN : (int) imin_dbl;
//			imax = (imax_dbl > INT_MAX) ? INT_MAX : (int) imax_dbl;
//			jmin = (jmin_dbl < INT_MIN) ? INT_MIN : (int) jmin_dbl;
//			jmax = (jmax_dbl > INT_MAX) ? INT_MAX : (int) jmax_dbl;
//
//			/*  Set the outbnd flag if the facet extends beyond the POS window  */
//			if ((imin < (-pn)) || (imax > pn) || (jmin < (-pn))	|| (jmax > pn)) {
//				posvis_streams_outbnd = 1;
//				atomicExch(&outbndarr[frm], 1);
//			}
//
//			/* Figure out if facet projects at least partly within POS window;
//			 * if it does, look at each "contained" POS pixel and get the
//			 * z-coordinate and cos(scattering angle)           */
//			i1 = MAX(imin, -pn);		j1 = MAX(jmin, -pn);
//			i2 = MIN(imax,  pn);		j2 = MIN(jmax,  pn);
//
//			if (i1 > pn || i2 < -pn || j1 > pn || j2 < -pn) {
//
//				/* Facet is entirely outside the POS frame: just keep track of
//				 * changed POS region     */
//				dev_POSrect_gpu64(pos, src, imin_dbl, imax_dbl, jmin_dbl, jmax_dbl,
//						ijminmax_overall, frm);
//
//			} else {
//
////				dev_POSrect_gpu64(pos, src, (double)i1, (double)i2, (double)j1,
////						(double)j2, ijminmax_overall, frm);
//
//				/* Assign vertex normals if smoothing is enabled */
//				if (pvst_smooth) {
//					/* Assign temp. normal components as float3 */
//					tv0.x = verts[0]->v[fidx.x].n[0];	tv0.y = verts[0]->v[fidx.x].n[1];
//					tv0.z = verts[0]->v[fidx.x].n[2];	tv1.x = verts[0]->v[fidx.y].n[0];
//					tv1.y = verts[0]->v[fidx.y].n[1];	tv1.z = verts[0]->v[fidx.y].n[2];
//					tv2.x = verts[0]->v[fidx.z].n[0];	tv2.y = verts[0]->v[fidx.z].n[1];
//					tv2.z = verts[0]->v[fidx.z].n[2];
//					n1n0.x = tv1.x - tv0.x;	n1n0.y = tv1.y - tv0.y;	n1n0.z = tv1.z - tv0.z;
//					tv2.x -= tv1.x;	tv2.y -= tv1.y; tv2.z -= tv1.z;
//				}
//
//				/* Precalculate s and t components */
//				double a, b, c, d, e, h, ti, tj, si, sj, si0, sj0, ti0, tj0, sz, tz;
//				a = i1*kmpxl - v0.x;
//				b = v2.y - v1.y;
//				c = v2.x - v1.x;
//				d = j1*kmpxl - v0.y;
//				e = v1.x - v0.x;
//				h = v1.y - v0.y;
//				den = e*b - c*h;
//				ti = -h*kmpxl/den;
//				tj = e*kmpxl/den;
//				si = b*kmpxl/den;
//				sj = -c*kmpxl/den;
//				si0 = (a*b - c*d)/den;
//				ti0 = (e*d -a*h)/den;
//				sz = v1.z - v0.z;
//				tz = v2.z - v1.z;
//
//				/* Facet is at least partly within POS frame: find all POS
//				 * pixels whose centers project onto this facet  */
//				for (i = i1; i <= i2; i++) {
//
//					sj0 = si0;	/* Initialize this loop's base sj0, tj0 */
//					tj0 = ti0;
//
//					for (j = j1; j <= j2; j++) {
//
//						s = sj0;
//						t = tj0;
//
//						if ((s >= -SMALLVAL) && (s <= 1.0 + SMALLVAL)) {// &&
//							if(	(t >= -SMALLVAL) && (t <= s + SMALLVAL))	{
//
//							/* Compute z-coordinate of pixel center: its
//							 * distance measured from the origin towards
//							 * Earth.    */
//							z = v0.z + s*sz + t*tz;
////							if (src)	fac = pos[frm]->fill[i][j];
////							if (!src)	fac = pos[frm]->f[i][j];
//
//							/* Following line replaces the previous if check
//							 * for z > zz[i][j]
//							 * atomicMaxf returns the value that was sitting
//							 * at zzf[pxa] at time of call.  So if that value
//							 * matches the z we compared to*/
//							if (src)	old = atomicMax64(&pos[frm]->zill[i][j], z);
//							if (!src)	old = atomicMax64(&pos[frm]->z[i][j], z);
//							if (old < z){
//								/* Next line assigns distance of POS pixel
//								 * center from COM towards Earth; that is,
//								 * by changing zz,it changes pos->z or
//								 * pos->zill                */
//								/* following line is a first time z calc
//								 * for this pixel  */
////									if (fac < 0){
////										atomicAdd(&dbg_cntr, 1);
////										if (src)
////											atomicExch((unsigned long long int*)&pos[frm]->zill[i][j], __double_as_longlong(z));
////
////										else
////											atomicExch((unsigned long long int*)&pos[frm]->z[i][j], __double_as_longlong(z));
////									}
////
//								if (pvst_smooth) {
//
//									/* Get pvs_smoothed version of facet unit
//									 * normal: Take the linear combination
//									 * of the three vertex normals; trans-
//									 * form from body to observer coordina-
//									 * tes; and make sure that it points
//									 * somewhat in our direction.         */
//									n.x = tv0.x + s * n1n0.x + t * tv2.x;
//									n.y = tv0.y + s * n1n0.y + t * tv2.y;
//									n.z = tv0.z + s * n1n0.z + t * tv2.z;
//									dev_cotrans3(&n, oa, n, 1, frm);
//									dev_normalize3(&n);
//								}
//
//								/* Determine scattering and/or incidence
//								 * angles. Next lines change pos->cose/
//								 * cosill. If bistatic (lightcurves), where
//								 * we are viewing from Earth (src = 0),
//								 * pos->cosi is also changed.                 */
//								if (n.z > 0.0) {
//									if (src)
//										atomicExch((unsigned long long int*)&pos[frm]->cosill[i][j],
//												__double_as_longlong(n.z));
//									else
//										atomicExch((unsigned long long int*)&pos[frm]->cose[i][j],
//												__double_as_longlong(n.z));
//
//									if ((!src) && (pos[frm]->bistatic)) {
//
//										double temp = dev_dot_d3(n,usrc[frm]);
//										atomicExch((unsigned long long int*)&pos[frm]->cosi[i][j],
//												__double_as_longlong(temp));
//										if (pos[frm]->cosi[i][j] <= 0.0)
//											pos[frm]->cose[i][j] = 0.0;
//									}
//								}
//								dev_POSrect_gpu64(pos, src, (double)i, (double)i, (double)j,
//										(double)j, ijminmax_overall, frm);
//								/* Next lines change pos->body/bodyill,
//								 * pos->comp/compill, pos->f/fill          */
//								if (src) {
//									pos[frm]->bodyill[i][j] = body;
//									pos[frm]->compill[i][j] = comp;
//									pos[frm]->fill[i][j] = f;
//								} else {
//									pos[frm]->body[i][j] = body;
//									pos[frm]->comp[i][j] = comp;
//									atomicExch(&pos[frm]->f[i][j], f);
//								}
//
//							} /* end if (no other facet yet blocks this facet from view) */
//							} /* end if 0 <= t <= s (facet center is "in" this POS pixel) */
//						} /* end if 0 <= s <= 1 */
//
//						sj0 += sj;
//						tj0 += tj;
//					} /* end j-loop over POS rows */
//					/* Modify s and t step-wise for the next i-iteration of the pixel loop */
//					si0 += si;
//					ti0 += ti;
//
//				} /* end i-loop over POS columns */
//			} /* end else of if (i1 > pos->n || i2 < -pos->n || j1 > pos->n || j2 < -pos->n) */
//		} /* End if (n[2] > 0.0) */
//	} /* end if (f < nf) */
//}
//
//__global__ void posvis_outbnd_krnl32(struct pos_t **pos,
//		int *outbndarr, float4 *ijminmax_overall, int size, int start) {
//	/* nfrm_alloc-threaded kernel */
//	int posn, f = blockIdx.x * blockDim.x + threadIdx.x + start;
//	double xfactor, yfactor;
//	if (f <size) {
//
////		printf("dbg_cntr in posvis_gpu32 = %i\n", dbg_cntr);
//
//		if (outbndarr[f]) {
//			/* ijminmax_overall.w = imin_overall
//			 * ijminmax_overall.x = imax_overall
//			 * ijminmax_overall.y = jmin_overall
//			 * ijminmax_overall.z = jmax_overall	 */
//			posn = pos[f]->n;
//			xfactor = (MAX( ijminmax_overall[f].x,  posn) -
//					MIN( ijminmax_overall[f].w, -posn) + 1) / (2*posn+1);
//			yfactor = (MAX( ijminmax_overall[f].z,  posn) -
//					MIN( ijminmax_overall[f].y, -posn) + 1) / (2*posn+1);
//			pos[f]->posbnd_logfactor = log(xfactor*yfactor);
//		}
//	}
//}
//
//__global__ void posvis_outbnd_krnl64(struct pos_t **pos,
//		int *outbndarr, double4 *ijminmax_overall, int size, int start) {
//	/* nfrm_alloc-threaded kernel */
//	int posn, f = blockIdx.x * blockDim.x + threadIdx.x + start;
//	double xfactor, yfactor;
//	if (f <size) {
//
////		printf("dbg_cntr in posvis_gpu64 = %i\n", dbg_cntr);
//
//		if (outbndarr[f]) {
//			/* ijminmax_overall.w = imin_overall
//			 * ijminmax_overall.x = imax_overall
//			 * ijminmax_overall.y = jmin_overall
//			 * ijminmax_overall.z = jmax_overall	 */
//			posn = pos[f]->n;
//			xfactor = (MAX( ijminmax_overall[f].x,  posn) -
//					MIN( ijminmax_overall[f].w, -posn) + 1) / (2*posn+1);
//			yfactor = (MAX( ijminmax_overall[f].z,  posn) -
//					MIN( ijminmax_overall[f].y, -posn) + 1) / (2*posn+1);
//			pos[f]->posbnd_logfactor = log(xfactor*yfactor);
//		}
//	}
//}
//
//__host__ int posvis_gpu32(
//		struct par_t *dpar,
//		struct mod_t *dmod,
//		struct dat_t *ddat,
//		struct pos_t **pos,
//		struct vertices_t **verts,
//		float3 orbit_offset,
//		int *posn,
//		int *outbndarr,
//		int set,
//		int nfrm_alloc,
//		int src,
//		int nf,
//		int body, int comp, unsigned char type, cudaStream_t *pv_stream,
//		int src_override) {
//
//	int f, outbnd, smooth, start;
//	dim3 BLK,THD, BLKfrm, THD64;
//	float4 *ijminmax_overall;
//	float3 *oa, *usrc;
//
//	/* Launch parameters for the facet_streams kernel */
//	THD.x = maxThreadsPerBlock;	THD64.x = 64;
//	BLK.x = floor((THD.x - 1 + nf) / THD.x);
//	BLKfrm.x = floor((THD64.x - 1 + nfrm_alloc)/THD64.x);
//
//	/* Set up the offset addressing for lightcurves if this is a lightcurve */
//	if (type == LGHTCRV)	start = 1;	/* fixes the lightcurve offsets */
//	else 					start = 0;
//
//	int oasize = nfrm_alloc*3;
//	/* Allocate temporary arrays/structs */
//	gpuErrchk(cudaMalloc((void**)&ijminmax_overall, sizeof(float4) * nfrm_alloc));
//	gpuErrchk(cudaMalloc((void**)&oa, sizeof(float3) * oasize));
//	gpuErrchk(cudaMalloc((void**)&usrc, sizeof(float3) * nfrm_alloc));
//
//	posvis_init_krnl32<<<BLKfrm,THD64>>>(dpar, pos, ijminmax_overall, oa, usrc,
//			outbndarr, comp, start, src, nfrm_alloc, set, src_override);
//	checkErrorAfterKernelLaunch("posvis_init_krnl32");
//
//	for (f=start; f<nfrm_alloc; f++) {
//		/* Now the main facet kernel */
//		posvis_facet_krnl32mod<<<BLK,THD, 0, pv_stream[f-start]>>>(pos, verts,
//				ijminmax_overall, orbit_offset, oa, usrc,	src, body, comp,
//				nf, f, smooth, outbndarr, set);
//	}
//	checkErrorAfterKernelLaunch("posvis_facet_krnl32");
//
//	/* Take care of any posbnd flags */
//	posvis_outbnd_krnl32<<<BLKfrm,THD64>>>(pos,
//			outbndarr, ijminmax_overall, nfrm_alloc, start);
//	checkErrorAfterKernelLaunch("posvis_outbnd_krnl32");
//	gpuErrchk(cudaMemcpyFromSymbol(&outbnd, posvis_streams_outbnd, sizeof(int), 0,
//			cudaMemcpyDeviceToHost));
//
////	int n = 200;
////	int npixels = 401*401;
////	f = 0;
////	dbg_print_pos_arrays_full(pos, f, npixels, n);
//
//	/* Free temp arrays, destroy streams and timers, as applicable */
//
//	cudaFree(ijminmax_overall);
//	cudaFree(oa);
//	cudaFree(usrc);
//
//	return outbnd;
//}
//
