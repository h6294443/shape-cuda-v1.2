///*****************************************************************************************
// posvis.c
//
// Fill in the portion of a plane-of-sky image due to a particular model component: Assign
// each relevant POS pixel a z-value in observer coordinates (distance from the origin
// towards Earth) and a value of cos(scattering angle).
//
// Return 1 if any portion of this component lies outside the specified POS window,
// 0 otherwise.
//
// If the "src" argument is true, the "observer" is the Sun rather than Earth, and
// "plane-of-sky" becomes "projection as viewed from the Sun."
//
// Modified 2014 February 20 by CM:
// Allow facets that partly project outside the POS frame to contribute to the POS frame
// (thus avoiding see-through "holes" in the model at the edge of a POS image)
//
// Modified 2010 May 18 by CM:
// Bug fix: When checking if a POS pixel hasn't already been assigned
// values during a previous call to posvis for a different component,
// check for fac[i][j] < 0 rather than cosa[i][j] == 0.0, since for
// bistatic situations the latter condition will also be true for
// pixels centered on Earth-facing facets that don't face the Sun
//
// Modified 2009 July 2 by CM:
// Eliminate the check that facets are "active": this term is now being
// interpreted to mean "not lying interior to the model," so the
// check is unnecessary and the determination of active vs. inactive
// status is inaccurate for half-exposed facets at the intersections
// between model components
//
// Modified 2009 April 3 by CM:
// Compute the "posbnd_logfactor" parameter: if the model extends beyond
// the POS frame, posbnd_logfactor is set to the logarithm of the
// ratio of the area that would have been required to "contain" the
// entire model divided by the area of the actual POS frame
// Work with floating-point pixel numbers (imin_dbl, etc.), at least
// initially, in case the sky rendering for a model with illegal
// parameters would involve huge pixel numbers that exceed the
// limits for valid integers
//
// Modified 2007 August 4 by CM:
// Add "orbit_offset" and "body" parameters and remove "facet" parameter
// Add body, bodyill, comp, and compill matrices for POS frames
//
// Modified 2006 June 21 by CM:
// For POS renderings, change res to km_per_pixel
//
// Modified 2005 September 19 by CM:
// Allow for roundoff error when determining which POS pixels project
// onto each model facet
//
// Modified 2005 June 27 by CM:
// Renamed "round" function to "iround" to avoid conflicts
//
// Modified 2005 June 22 by CM:
// Slightly modified some comments
//
// Modified 2005 January 25 by CM:
// Take care of unused and uninitialized variables
//
// Modified 2004 December 19 by CM:
// Added more comments
// Put update of rectangular POS area into "POSrect" routine and applied it
// even to facets which lie outside the POS frame
//
// Modified 2004 Feb 11 by CM:
// Added comments
//
// Modified 2003 May 5 by CM:
// Removed redundant coordinate transformation of the unit normal n
// for the no-pvs_smoothing case
// *****************************************************************************************/
//extern "C" {
//#include "../shape/head.h"
//#include <limits.h>
//}
//
//#define maxbins 100
//__device__ int posvis_tiled_outbnd, posvis_tiled_smooth;
//
///* Note that the following custom atomic functions must be declared in each
// * file it is needed (consequence of being a static device function) */
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
//__global__ void posvis_tiled_init_krnl64(
//		struct par_t *dpar,
//		struct pos_t **pos,
//		double4 *ijminmax_overall,
//		double3 *oa,
//		double3 *usrc,
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
//			posvis_tiled_outbnd = 0;
//			posvis_tiled_smooth = dpar->pos_smooth;
//			if (src_override)	posvis_tiled_smooth = 0;
//		}
//		ijminmax_overall[f].w = ijminmax_overall[f].y = HUGENUMBER;
//		ijminmax_overall[f].x = ijminmax_overall[f].z = -HUGENUMBER;
//		pos[f]->posbnd_logfactor = 0.0;
//
//		dev_mtrnsps2(oa, pos[f]->ae, f);
//		if (src) {
//			/* We're viewing the model from the sun: at the center of each pixel
//			 * in the projected view, we want cos(incidence angle), distance from
//			 * the COM towards the sun, and the facet number.                */
//			dev_mmmul2(oa, pos[f]->se, oa, f); /* oa takes ast into sun coords           */
//		} else {
//			/* We're viewing the model from Earth: at the center of each POS pixel
//			 * we want cos(scattering angle), distance from the COM towards Earth,
//			 * and the facet number.  For bistatic situations (lightcurves) we also
//									 want cos(incidence angle) and the unit vector towards the source.     */
//			dev_mmmul2(oa, pos[f]->oe, oa, f); /* oa takes ast into obs coords */
//			if (pos[f]->bistatic) {
//				usrc[f].x = usrc[f].y = 0.0; /* unit vector towards source */
//				usrc[f].z = 1.0;
//				dev_cotrans1(&usrc[f], pos[f]->se, usrc[f], -1);
//				dev_cotrans1(&usrc[f], pos[f]->oe, usrc[f], 1); /* in observer coordinates */
//			}
//		}
//		outbndarr[f] = 0;
//	}
//}
//
//__global__ void transform_facet_normals_krnl64a(
//		struct mod_t *dmod,
//		struct pos_t **pos,
//		struct vertices_t **verts,
//		double4 *ijminmax_overall,
//		double3 orbit_offs,
//		double3 *oa,
//		double3 *usrc,
//		int *outbndarr,
//		int nf,
//		int frm,
//		int src,
//		int blockSize)
//{
//	/* This kernel launches 256 threads, performs a grid-stride loop through
//	 * all model facets and transforms each facet normal with oa[frm] and stores
//	 * the result back to dmod if n.z > 0.0. It also determines and stores the
//	 * facet and global model bounding box via i1,i2,j1,j2 and xlim/ylim.
//	 * These quantities are stored in pos_facet_t structures inside each frame's
//	 * pos. 	 */
//
//	/* Declare kernel variables */
//	__shared__ int pn;
//	__shared__ double kmpxl;
//	int imin, jmin, imax, jmax, i1, i2, j1, j2;
//	int3 fidx;
//	double imin_dbl, jmin_dbl, imax_dbl, jmax_dbl;
//	double3 n, v0, v1, v2;
//
//	/* Initialize the shared variables (accessed by every thread) */
//	if (threadIdx.x==0) {
//		pn = pos[frm]->n;
//		kmpxl = pos[frm]->km_per_pixel;
//	}
//	__syncthreads();
//
//	/* Do a grid-stride loop on all facets */
//	for (int f=threadIdx.x; f<nf; f+=blockSize) {
//		/* Get vertex indices of the three vertices making up the facet */
//		fidx.x = verts[0]->f[f].v[0]; fidx.y = verts[0]->f[f].v[1];
//		fidx.z = verts[0]->f[f].v[2];
//
//		/* Copy each vertex over to thread register memory */
//		v0.x = verts[0]->v[fidx.x].x[0];	v0.y = verts[0]->v[fidx.x].x[1];
//		v0.z = verts[0]->v[fidx.x].x[2];	v1.x = verts[0]->v[fidx.y].x[0];
//		v1.y = verts[0]->v[fidx.y].x[1];	v1.z = verts[0]->v[fidx.y].x[2];
//		v2.x = verts[0]->v[fidx.z].x[0];	v2.y = verts[0]->v[fidx.z].x[1];
//		v2.z = verts[0]->v[fidx.z].x[2];
//
//		/* Get the normal to this facet in body-fixed (asteroid) coordinates
//		 * and convert it to observer coordinates     */
//		n.x = verts[0]->f[f].n[0];
//		n.y = verts[0]->f[f].n[1];
//		n.z = verts[0]->f[f].n[2];
//		dev_cotrans3(&n, oa, n, 1, frm);
//
//		/* Check if this facet is visible - is the facet normal pointing
//		 * roughly at the observer?				 */
//		if (n.z > 0.0) {
//			/* First, store the transformed normal back to the model and increase
//			 * visible facet counter */
//			pos[frm]->facet[f].nt.x = n.x;
//			pos[frm]->facet[f].nt.y = n.y;
//			pos[frm]->facet[f].nt.z = n.z;
//			/* Convert the 3 vertex coordinates from body to observer
//			 * coordinates; orbit_offset is the center-of-mass offset
//			 * (in observer coordinates) for this model at this frame's
//			 * epoch due to orbital motion, in case the model is half of
//			 *  a binary system.  */
//			dev_cotrans3(&v0, oa, v0, 1, frm);
//			dev_cotrans3(&v1, oa, v1, 1, frm);
//			dev_cotrans3(&v2, oa, v2, 1, frm);
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
//				posvis_tiled_outbnd = 1;
//				atomicExch(&outbndarr[frm], 1);
//			}
//
//			/* Figure out if facet projects at least partly within POS window;
//			 * if it does, look at each "contained" POS pixel and get the
//			 * z-coordinate and cos(scattering angle)           */
//			i1 = MAX(imin, -pn);	i2 = MIN(imax, pn);
//			j1 = MAX(jmin, -pn);	j2 = MIN(jmax, pn);
//
//			pos[frm]->facet[f].ilim.x = i1;
//			pos[frm]->facet[f].ilim.y = i2;
//			pos[frm]->facet[f].jlim.x = j1;
//			pos[frm]->facet[f].jlim.y = j2;
//			pos[frm]->facet[f].v0t = v0;
//			pos[frm]->facet[f].v1t = v1;
//			pos[frm]->facet[f].v2t = v2;
//
//			/* Now keep track of the global region */
//			if (i1 > pn || i2 < -pn || j1 > pn || j2 < -pn) {
//				/* Facet is entirely outside POS frame: just track POS region */
//				dev_POSrect_gpu64(pos, src, imin_dbl, imax_dbl, jmin_dbl,
//						jmax_dbl, ijminmax_overall, frm);
//
//			} else {
//				dev_POSrect_gpu64(pos, src, (double)i1, (double)i2,
//						(double)j1, (double)j2, ijminmax_overall, frm);
//			}
//		}
//		else {
//			/* The following makes a check in the bin_facets_krnl64 kernel easier */
//			pos[frm]->facet[f].nt.x = -1.0;
//			pos[frm]->facet[f].nt.y = -1.0;
//			pos[frm]->facet[f].nt.z = -1.0;
//		}
//	}
//}
//
//__global__ void transform_facet_normals_krnl64b(
//		struct mod_t *dmod,
//		struct pos_t **pos,
//		struct vertices_t **verts,
//		double4 *ijminmax_overall,
//		double3 orbit_offs,
//		double3 *oa_gm,
//		int *outbndarr,
//		int nf,
//		int src)
//{
//	/* This kernel launches nframes blocks of threads, performs a grid-stride loop through
//	 * all model facets and transforms each facet normal with oa[frm] and stores
//	 * the result back to dmod if n.z > 0.0. It also determines and stores the
//	 * facet and global model bounding box via i1,i2,j1,j2 and xlim/ylim.
//	 * These quantities are stored in pos_facet_t structures inside each frame's
//	 * pos. This kernel also uses shared memory for ijminmax_overall_sh, used
//	 * as temporary (faster) storage for pos window calculation.  Additionally,
//	 * the pos->xlim/ylim atomic operations have been moved to the very end of
//	 * this kernel to be processed just once instead of for every facet.	 */
//
//	/* Declare kernel variables */
//	__shared__ int pn;
//	__shared__ double kmpxl, oa_sh[3][3];
//	__shared__ double4 ijminmax_overall_sh;
//	int frm=blockIdx.x, imin, jmin, imax, jmax, i1, i2, j1, j2;
//	int3 fidx;
//	double imin_dbl, jmin_dbl, imax_dbl, jmax_dbl;
//	double3 n, v0, v1, v2;
//
//	/* Initialize the shared variables (accessed by every thread) */
//	if (threadIdx.x==0) {
//		pn = pos[frm]->n;
//		kmpxl = pos[frm]->km_per_pixel;
//		ijminmax_overall_sh.w = ijminmax_overall_sh.x =
//				ijminmax_overall_sh.y = ijminmax_overall_sh.z = 0.0f;
//
//		/* Load oa for this frame into shared memory */
//		oa_sh[0][0] = oa_gm[3*frm].x;	oa_sh[0][1] = oa_gm[3*frm].y;	oa_sh[0][2] = oa_gm[3*frm].z;
//		oa_sh[1][0] = oa_gm[3*frm+1].x;	oa_sh[1][1] = oa_gm[3*frm+1].y;	oa_sh[1][2] = oa_gm[3*frm+1].z;
//		oa_sh[2][0] = oa_gm[3*frm+2].x;	oa_sh[2][1] = oa_gm[3*frm+2].y;	oa_sh[2][2] = oa_gm[3*frm+2].z;
//	}
//	__syncthreads();
//
//	/* Do a grid-stride loop on all facets */
//	for (int f=threadIdx.x; f<nf; f+=blockDim.x) {
//		/* Get vertex indices of the three vertices making up the facet */
//		fidx.x = verts[0]->f[f].v[0]; fidx.y = verts[0]->f[f].v[1];
//		fidx.z = verts[0]->f[f].v[2];
//
//		/* Copy each vertex over to thread register memory */
//		v0.x = verts[0]->v[fidx.x].x[0];	v0.y = verts[0]->v[fidx.x].x[1];
//		v0.z = verts[0]->v[fidx.x].x[2];	v1.x = verts[0]->v[fidx.y].x[0];
//		v1.y = verts[0]->v[fidx.y].x[1];	v1.z = verts[0]->v[fidx.y].x[2];
//		v2.x = verts[0]->v[fidx.z].x[0];	v2.y = verts[0]->v[fidx.z].x[1];
//		v2.z = verts[0]->v[fidx.z].x[2];
//
//		/* Get the normal to this facet in body-fixed (asteroid) coordinates
//		 * and convert it to observer coordinates     */
//		n.x = verts[0]->f[f].n[0];
//		n.y = verts[0]->f[f].n[1];
//		n.z = verts[0]->f[f].n[2];
//		dev_cotrans1(&n, oa_sh, n, 1);
//
//		/* Check if this facet is visible - is the facet normal pointing
//		 * roughly at the observer?				 */
//		if (n.z > 0.0) {
//			/* First, store the transformed normal back to the model and increase
//			 * visible facet counter */
//			pos[frm]->facet[f].nt.x = n.x;
//			pos[frm]->facet[f].nt.y = n.y;
//			pos[frm]->facet[f].nt.z = n.z;
//			/* Convert the 3 vertex coordinates from body to observer
//			 * coordinates; orbit_offset is the center-of-mass offset
//			 * (in observer coordinates) for this model at this frame's
//			 * epoch due to orbital motion, in case the model is half of
//			 *  a binary system.  */
//			dev_cotrans1(&v0, oa_sh, v0, 1);
//			dev_cotrans1(&v1, oa_sh, v1, 1);
//			dev_cotrans1(&v2, oa_sh, v2, 1);
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
//				posvis_tiled_outbnd = 1;
//				atomicExch(&outbndarr[frm], 1);
//			}
//
//			/* Figure out if facet projects at least partly within POS window;
//			 * if it does, look at each "contained" POS pixel and get the
//			 * z-coordinate and cos(scattering angle)           */
//			i1 = MAX(imin, -pn);	i2 = MIN(imax, pn);
//			j1 = MAX(jmin, -pn);	j2 = MIN(jmax, pn);
//
//			pos[frm]->facet[f].ilim.x = i1;
//			pos[frm]->facet[f].ilim.y = i2;
//			pos[frm]->facet[f].jlim.x = j1;
//			pos[frm]->facet[f].jlim.y = j2;
//			pos[frm]->facet[f].v0t = v0;
//			pos[frm]->facet[f].v1t = v1;
//			pos[frm]->facet[f].v2t = v2;
//
//			/* Now keep track of the global region */
//			if (i1 > pn || i2 < -pn || j1 > pn || j2 < -pn) {
//				/* Facet is entirely outside POS frame: just track POS region */
//				dev_POSrect_gpu64_shared(imin_dbl,imax_dbl,	jmin_dbl,jmax_dbl,
//						&ijminmax_overall_sh, pn);
//
//			} else {
//				dev_POSrect_gpu64_shared((double)i1, (double)i2,
//						(double)j1, (double)j2, &ijminmax_overall_sh, pn);
//			}
//		}
//		else {
//			/* The following makes a check in the bin_facets_krnl64 kernel easier */
////			pos[frm]->facet[f].nt.x = -1.0;
////			pos[frm]->facet[f].nt.y = -1.0;
////			pos[frm]->facet[f].nt.z = -1.0;
//		}
//	}
//	__syncthreads();
//
//	/* Now write the POS frame window limits from shared mem back to global mem */
//	if (threadIdx.x==0) {
//		ijminmax_overall[frm].w = ijminmax_overall_sh.w;
//		ijminmax_overall[frm].x = ijminmax_overall_sh.x;
//		ijminmax_overall[frm].y = ijminmax_overall_sh.y;
//		ijminmax_overall[frm].z = ijminmax_overall_sh.z;
//
//		/*  Update the subset of the POS frame that contains the target  */
//		/* imin_dbl - ijminmax_overall[frm].w
//		 * imax_dbl - ijminmax_overall[frm].x
//		 * jmin_dbl - ijminmax_overall[frm].y
//		 * jmax_dbl - ijminmax_overall[frm].z
//		 */
//		int imin = (ijminmax_overall_sh.w < INT_MIN) ? INT_MIN : (int) ijminmax_overall_sh.w;
//		int imax = (ijminmax_overall_sh.x > INT_MAX) ? INT_MAX : (int) ijminmax_overall_sh.x;
//		int jmin = (ijminmax_overall_sh.y < INT_MIN) ? INT_MIN : (int) ijminmax_overall_sh.y;
//		int jmax = (ijminmax_overall_sh.z > INT_MAX) ? INT_MAX : (int) ijminmax_overall_sh.z;
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
//	}
//}
//
//__global__ void transform_facet_normals_krnl64c(
//		struct mod_t *dmod,
//		struct pos_t **pos,
//		struct vertices_t **verts,
//		double4 *ijminmax_overall,
//		double3 orbit_offs,
//		double3 *oa_gm,
//		int *outbndarr,
//		int nf,
//		int frm,
//		int src)
//{
//	/* This kernel launches 256 threads, performs a grid-stride loop through
//	 * all model facets and transforms each facet normal with oa[frm] and stores
//	 * the result back to dmod if n.z > 0.0. It also determines and stores the
//	 * facet and global model bounding box via i1,i2,j1,j2 and xlim/ylim.
//	 * These quantities are stored in pos_facet_t structures inside each frame's
//	 * pos. This kernel also uses shared memory for ijminmax_overall_sh, used
//	 * as temporary (faster) storage for pos window calculation.  Additionally,
//	 * the pos->xlim/ylim atomic operations have been moved to the very end of
//	 * this kernel to be processed just once instead of for every facet.	 */
//
//	/* Declare kernel variables */
//	__shared__ int pn;
//	__shared__ double kmpxl, oa_sh[3][3];
//	__shared__ double4 ijminmax_overall_sh;
//	int imin, jmin, imax, jmax, i1, i2, j1, j2;
//	int3 fidx;
//	double imin_dbl, jmin_dbl, imax_dbl, jmax_dbl;
//	double3 n, v0, v1, v2;
//
//	/* Initialize the shared variables (accessed by every thread) */
//	if (threadIdx.x==0) {
//		pn = pos[frm]->n;
//		kmpxl = pos[frm]->km_per_pixel;
//		ijminmax_overall_sh.w = ijminmax_overall_sh.x =
//				ijminmax_overall_sh.y = ijminmax_overall_sh.z = 0.0f;
//
//		/* Load oa for this frame into shared memory */
//		oa_sh[0][0] = oa_gm[3*frm].x;	oa_sh[0][1] = oa_gm[3*frm].y;	oa_sh[0][2] = oa_gm[3*frm].z;
//		oa_sh[1][0] = oa_gm[3*frm+1].x;	oa_sh[1][1] = oa_gm[3*frm+1].y;	oa_sh[1][2] = oa_gm[3*frm+1].z;
//		oa_sh[2][0] = oa_gm[3*frm+2].x;	oa_sh[2][1] = oa_gm[3*frm+2].y;	oa_sh[2][2] = oa_gm[3*frm+2].z;
//	}
//	__syncthreads();
//
//	/* Do a grid-stride loop on all facets */
//	for (int f=threadIdx.x; f<nf; f+=blockDim.x) {
//		/* Get vertex indices of the three vertices making up the facet */
//		fidx.x = verts[0]->f[f].v[0]; fidx.y = verts[0]->f[f].v[1];
//		fidx.z = verts[0]->f[f].v[2];
//
//		/* Copy each vertex over to thread register memory */
//		v0.x = verts[0]->v[fidx.x].x[0];	v0.y = verts[0]->v[fidx.x].x[1];
//		v0.z = verts[0]->v[fidx.x].x[2];	v1.x = verts[0]->v[fidx.y].x[0];
//		v1.y = verts[0]->v[fidx.y].x[1];	v1.z = verts[0]->v[fidx.y].x[2];
//		v2.x = verts[0]->v[fidx.z].x[0];	v2.y = verts[0]->v[fidx.z].x[1];
//		v2.z = verts[0]->v[fidx.z].x[2];
//
//		/* Get the normal to this facet in body-fixed (asteroid) coordinates
//		 * and convert it to observer coordinates     */
//		n.x = verts[0]->f[f].n[0];
//		n.y = verts[0]->f[f].n[1];
//		n.z = verts[0]->f[f].n[2];
//		dev_cotrans1(&n, oa_sh, n, 1);
//
//		/* Check if this facet is visible - is the facet normal pointing
//		 * roughly at the observer?				 */
//		if (n.z > 0.0) {
//			/* First, store the transformed normal back to the model and increase
//			 * visible facet counter */
//			pos[frm]->facet[f].nt.x = n.x;
//			pos[frm]->facet[f].nt.y = n.y;
//			pos[frm]->facet[f].nt.z = n.z;
//			/* Convert the 3 vertex coordinates from body to observer
//			 * coordinates; orbit_offset is the center-of-mass offset
//			 * (in observer coordinates) for this model at this frame's
//			 * epoch due to orbital motion, in case the model is half of
//			 *  a binary system.  */
//			dev_cotrans1(&v0, oa_sh, v0, 1);
//			dev_cotrans1(&v1, oa_sh, v1, 1);
//			dev_cotrans1(&v2, oa_sh, v2, 1);
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
//				posvis_tiled_outbnd = 1;
//				atomicExch(&outbndarr[frm], 1);
//			}
//
//			/* Figure out if facet projects at least partly within POS window;
//			 * if it does, look at each "contained" POS pixel and get the
//			 * z-coordinate and cos(scattering angle)           */
//			i1 = MAX(imin, -pn);	i2 = MIN(imax, pn);
//			j1 = MAX(jmin, -pn);	j2 = MIN(jmax, pn);
//
//			pos[frm]->facet[f].ilim.x = i1;
//			pos[frm]->facet[f].ilim.y = i2;
//			pos[frm]->facet[f].jlim.x = j1;
//			pos[frm]->facet[f].jlim.y = j2;
//			pos[frm]->facet[f].v0t = v0;
//			pos[frm]->facet[f].v1t = v1;
//			pos[frm]->facet[f].v2t = v2;
//
//			/* Now keep track of the global region */
//			if (i1 > pn || i2 < -pn || j1 > pn || j2 < -pn) {
//				/* Facet is entirely outside POS frame: just track POS region */
//				dev_POSrect_gpu64_shared(imin_dbl,imax_dbl,	jmin_dbl,jmax_dbl,
//						&ijminmax_overall_sh, pn);
//
//			} else {
//				dev_POSrect_gpu64_shared((double)i1, (double)i2,
//						(double)j1, (double)j2, &ijminmax_overall_sh, pn);
//			}
//		}
//		else {
//			/* The following makes a check in the bin_facets_krnl64 kernel easier */
//			pos[frm]->facet[f].nt.x = -1.0;
//			pos[frm]->facet[f].nt.y = -1.0;
//			pos[frm]->facet[f].nt.z = -1.0;
//		}
//	}
//	__syncthreads();
//
//	/* Now write the POS frame window limits from shared mem back to global mem */
//	if (threadIdx.x==0) {
//		ijminmax_overall[frm].w = ijminmax_overall_sh.w;
//		ijminmax_overall[frm].x = ijminmax_overall_sh.x;
//		ijminmax_overall[frm].y = ijminmax_overall_sh.y;
//		ijminmax_overall[frm].z = ijminmax_overall_sh.z;
//
//		/*  Update the subset of the POS frame that contains the target  */
//		/* imin_dbl - ijminmax_overall[frm].w
//		 * imax_dbl - ijminmax_overall[frm].x
//		 * jmin_dbl - ijminmax_overall[frm].y
//		 * jmax_dbl - ijminmax_overall[frm].z
//		 */
//		int imin = (ijminmax_overall_sh.w < INT_MIN) ? INT_MIN : (int) ijminmax_overall_sh.w;
//		int imax = (ijminmax_overall_sh.x > INT_MAX) ? INT_MAX : (int) ijminmax_overall_sh.x;
//		int jmin = (ijminmax_overall_sh.y < INT_MIN) ? INT_MIN : (int) ijminmax_overall_sh.y;
//		int jmax = (ijminmax_overall_sh.z > INT_MAX) ? INT_MAX : (int) ijminmax_overall_sh.z;
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
//	}
//}
//
//__global__ void bin_facets_krnl64a(struct pos_t **pos,
//		struct vertices_t **verts,
//		int ***facet_index,
//		int **entries,
//		int nf,
//		int frm,
//		int *n_tiles,
//		int *n_tiles_x,
//		int	*n_tiles_y,
//		int tile_size)
//{
//	/* This kernel is responsible for binning visible model facets according to
//	 * which screen tile they appear on. Each facet can belong to 1, 2, or 4
//	 * different facets.  (If the size of individual triangles should exceed
//	 * the tile size, this is no longer true.)
//	 * The kernel version has just one thread block with 1024 threads. It uses a grid-
//	 * stride loop to cover all facets
//	 */
//	int f, current_i, next_i, current_j, next_j, i1, i2, j1, j2, bi, bin, old_indx;
//	__shared__ int2 xlim, ylim; /* These are the global pos limits */
//	extern __shared__ int addr_index[];	/* Used for the facet_index entries */
//
//	/* Initialize shared variables that will be accessed by every thread */
//	if (threadIdx.x==0) {
//		xlim.x = pos[frm]->xlim[0];
//		xlim.y = pos[frm]->xlim[1];
//		ylim.x = pos[frm]->ylim[0];
//		ylim.y = pos[frm]->ylim[1];
//		for (bin=0; bin<n_tiles[frm]; bin++)
//			addr_index[bin] = 0;
//	}
//	__syncthreads();
//
//	/* Do grid-stride loop through all facets in model  */
//	for (f=threadIdx.x; f<nf; f+=blockDim.x) {
//		/* Weed out any facets not visible to observer */
//		if (pos[frm]->facet[f].nt.z > 0.0) {
//			bi = 0;	/* Bin index for the four facet bin entries*/
//			/* Copy facet limits into register memory for faster access */
//			i1 = pos[frm]->facet[f].ilim.x;
//			i2 = pos[frm]->facet[f].ilim.y;
//			j1 = pos[frm]->facet[f].jlim.x;
//			j2 = pos[frm]->facet[f].jlim.y;
//
//			/* Now check where the current facet lies, stepping through each
//			 * tile */
//			for (int k=0; k<n_tiles_y[frm]; k++) {
//				current_j = ylim.x + k * tile_size;
//				next_j = current_j + tile_size;
//				for (int n=0; n<n_tiles_x[frm]; n++) {
//					bin = k*n_tiles_x[frm] + n;
//					current_i = xlim.x + n * tile_size;
//					next_i = current_i + tile_size;
//
//					/* If i1 or i2 AND j1 or j2 fall into this tile, register it */
//					if ((i1>=current_i && i1<next_i)  || (i2>=current_i && i2<next_i)) {
//						if ((j1>=current_j && j1<next_j)  || (j2>=current_j && j2<next_j)) {
//							pos[frm]->facet[f].bin[bi] = bin;
//							old_indx = atomicAdd(&addr_index[bin], 1);
//							facet_index[frm][bin][old_indx] = f;
//							atomicAdd(&entries[frm][bin], 1);
//							bi++;
//						}
//					}
//				}
//			}
//		}
//	}
//}
//
//__global__ void bin_facets_krnl64b(struct pos_t **pos,
//		struct vertices_t **verts,
//		int ***facet_index,
//		int **entries,
//		int **addr_index,
//		int nf,
//		int frm,
//		int *n_tiles,
//		int *n_tiles_x,
//		int	*n_tiles_y,
//		int tile_size)
//{
//	/* This kernel is responsible for binning visible model facets according to
//	 * which screen tile they appear on. Each facet can belong to 1, 2, or 4
//	 * different facets.  (If the size of individual triangles should exceed
//	 * the tile size, this is no longer true.)
//	 * The kernel version uses nf-threads with as many thread blocks as it
//	 * takes, considering the previously defined maxThreadsPerBlock. Because of
//	 * this, the addr_index array is in global memory (instead of shared). */
//
//	int f, current_i, next_i, current_j, next_j, i1, i2, j1, j2, bi, bin, old_indx;
//	__shared__ int2 xlim, ylim; /* These are the global pos limits */
//
//	f = blockDim.x * blockIdx.x + threadIdx.x;
//
//	/* Initialize shared variables that will be accessed by every thread */
//	if (threadIdx.x==0) {
//		xlim.x = pos[frm]->xlim[0];
//		xlim.y = pos[frm]->xlim[1];
//		ylim.x = pos[frm]->ylim[0];
//		ylim.y = pos[frm]->ylim[1];
////		for (bin=0; bin<n_tiles[frm]; bin++)
////			addr_index[frm][bin] = 0;
//	}
//	__syncthreads();
//
//	/* Do grid-stride loop through all facets in model  */
//	if (f < nf) {
//		/* Weed out any facets not visible to observer */
//		if (pos[frm]->facet[f].nt.z > 0.0) {
//			bi = 0;	/* Bin index for the four facet bin entries*/
//			/* Copy facet limits into register memory for faster access */
//			i1 = pos[frm]->facet[f].ilim.x;
//			i2 = pos[frm]->facet[f].ilim.y;
//			j1 = pos[frm]->facet[f].jlim.x;
//			j2 = pos[frm]->facet[f].jlim.y;
//
//			/* Now check where the current facet lies, stepping through each
//			 * tile */
//			for (int k=0; k<n_tiles_y[frm]; k++) {
//				current_j = ylim.x + k * tile_size;
//				next_j = current_j + tile_size;
//				for (int n=0; n<n_tiles_x[frm]; n++) {
//					bin = k*n_tiles_x[frm] + n;
//					current_i = xlim.x + n * tile_size;
//					next_i = current_i + tile_size;
//
//					/* If i1 or i2 AND j1 or j2 fall into this tile, register it */
//					if ((i1>=current_i && i1<next_i)  || (i2>=current_i && i2<next_i)) {
//						if ((j1>=current_j && j1<next_j)  || (j2>=current_j && j2<next_j)) {
//							pos[frm]->facet[f].bin[bi] = bin;
//							old_indx = atomicAdd(&addr_index[frm][bin], 1);
//							facet_index[frm][bin][old_indx] = f;
//							atomicAdd(&entries[frm][bin], 1);
//							bi++;
//						}
//					}
//				}
//			}
//		}
//	}
//}
//
//__global__ void bin_facets_krnl64c(struct pos_t **pos,
//		struct vertices_t **verts,
//		int ***facet_index,
//		int **entries,
//		int nf,
//		int *n_tiles,
//		int *n_tiles_x,
//		int	*n_tiles_y,
//		int tile_size)
//{
//	/* This kernel is responsible for binning visible model facets according to
//	 * which screen tile they appear on. Each facet can belong to 1, 2, or 4
//	 * different facets.  (If the size of individual triangles should exceed
//	 * the tile size, this is no longer true.)
//	 * The kernel version uses nframes-thread blocks to cover all frames in one
//	 * go. Each block has its own __shared___ addr_index array and each of them
//	 * uses a grid-stride loop to cover all facets
//	 */
//	int f, frm, current_i, next_i, current_j, next_j, i1, i2, j1, j2, bi, bin, old_indx;
//	__shared__ int2 xlim, ylim; /* These are the global pos limits */
//	__shared__ int addr_index[160];	/* This allows for 40x40 tiles (at 32x32
//	tile size, allowing for a maximum POS resolution of 1280x1280 pixels  */
//	frm = blockIdx.x;
//
//	/* Initialize shared variables that will be accessed by every thread */
//	if (threadIdx.x==0) {
//		xlim.x = pos[frm]->xlim[0];
//		xlim.y = pos[frm]->xlim[1];
//		ylim.x = pos[frm]->ylim[0];
//		ylim.y = pos[frm]->ylim[1];
//		for (bin=0; bin<n_tiles[frm]; bin++)
//			addr_index[bin] = 0;
//	}
//	__syncthreads();
//
//	/* Do grid-stride loop through all facets in model  */
//	for (f=threadIdx.x; f<nf; f+=blockDim.x) {
//		/* Weed out any facets not visible to observer */
//		if (pos[frm]->facet[f].nt.z > 0.0) {
//			bi = 0;	/* Bin index for the four facet bin entries*/
//			/* Copy facet limits into register memory for faster access */
//			i1 = pos[frm]->facet[f].ilim.x;
//			i2 = pos[frm]->facet[f].ilim.y;
//			j1 = pos[frm]->facet[f].jlim.x;
//			j2 = pos[frm]->facet[f].jlim.y;
//
//			/* Now check where the current facet lies, stepping through each
//			 * tile */
//			for (int k=0; k<n_tiles_y[frm]; k++) {
//				current_j = ylim.x + k * tile_size;
//				next_j = current_j + tile_size;
//				for (int n=0; n<n_tiles_x[frm]; n++) {
//					bin = k*n_tiles_x[frm] + n;
//					current_i = xlim.x + n * tile_size;
//					next_i = current_i + tile_size;
//
//					/* If i1 or i2 AND j1 or j2 fall into this tile, register it */
//					if ((i1>=current_i && i1<next_i)  || (i2>=current_i && i2<next_i)) {
//						if ((j1>=current_j && j1<next_j)  || (j2>=current_j && j2<next_j)) {
//							pos[frm]->facet[f].bin[bi] = bin;
//							old_indx = atomicAdd(&addr_index[bin], 1);
//							facet_index[frm][bin][old_indx] = f;
//							atomicAdd(&entries[frm][bin], 1);
//							bi++;
//						}
//					}
//				}
//			}
//		}
//	}
//}
//
//__global__ void radar_raster_krnl64(struct pos_t **pos,
//		struct vertices_t **verts,
//		double3 *oa_gm,
//		int ***facet_index,
//		int **entries,
//		int nf,
//		int frm,
//		int *n_tiles,
//		int *n_tiles_x,
//		int	*n_tiles_y,
//		int tile_size,
//		int tile_x) {
//
//	/* This kernel performs the rasterization tile by tile.  Each thread block
//	 * is responsible for one tile. */
//	/* Determine which tile this thread block is responsible for and
//	 * which element of the thread block this thread is. 	 */
//	int bin = blockIdx.x;
//	int index = threadIdx.x;
//
//	/* Declare the shared memory variables and others */
//	__shared__ double pos_z[32][32];//[55][55];	/* One per thread block */
//	__shared__ double pos_cose[32][32];//[55][55];	/* One per thread block */
//	__shared__ int2 bn;
//	__shared__ int xlim, ylim, offsetx, offsety;
//	__shared__ double kmpxl;
//	__shared__ double oa_sh[3][3];
//	int i, j, ig, jg, i1, i2, j1, j2;	/* ig,jg are global indices */
//	int tile_i1, tile_i2, tile_j1, tile_j2, fct_indx;
//	double3 v0, v1, v2, n, tv0, tv1, tv2, n1n0;
//	int3 fidx;
//
//	/* Initialize the shared memory arrays with grid-stride loop */
//	for (int index=threadIdx.x; index<tile_size; index+=blockDim.x) {
//		i = index % tile_x;
//		j = index / tile_x;
//		pos_z[i][j] = -1e20;
//		pos_cose[i][j] = 0.0;
//	}
//	__syncthreads();
//
//	/* Load variables used by every thread (per block) to shared memory for
//	 * faster access */
//	if (threadIdx.x==0) {
//		xlim = pos[frm]->xlim[0];
//		ylim = pos[frm]->ylim[0];
//		kmpxl = pos[frm]->km_per_pixel;
//		bn.x = bin % n_tiles_x[frm];
//		bn.y = bin / n_tiles_x[frm];
//
//		/* Calculate the pixel offsets needed to go back and forth between
//		 * tiled POS space for this block's tile and global POS space 	 */
//		offsetx = xlim + tile_x * bn.x;
//		offsety = ylim + tile_x * bn.y;
//
//		/* Load oa for this frame into shared memory */
//		if (posvis_tiled_smooth) {
//			oa_sh[0][0] = oa_gm[3*frm].x;	oa_sh[0][1] = oa_gm[3*frm].y;	oa_sh[0][2] = oa_gm[3*frm].z;
//			oa_sh[1][0] = oa_gm[3*frm+1].x;	oa_sh[1][1] = oa_gm[3*frm+1].y;	oa_sh[1][2] = oa_gm[3*frm+1].z;
//			oa_sh[2][0] = oa_gm[3*frm+2].x;	oa_sh[2][1] = oa_gm[3*frm+2].y;	oa_sh[2][2] = oa_gm[3*frm+2].z;
//		}
//
//	}
//	__syncthreads();
//
//	/* Using grid-stride loop, step through all facet entries for each bin where
//	 * each thread block is responsible for one bin/tile	 */
//	for (index=threadIdx.x; index<entries[frm][bin]; index+=blockDim.x) {
//
//		/* Load facet index into registers */
//		fct_indx = facet_index[frm][bin][index];
//
//		/* Load transformed facet vertices into registers */
//		v0 = pos[frm]->facet[fct_indx].v0t;
//		v1 = pos[frm]->facet[fct_indx].v1t;
//		v2 = pos[frm]->facet[fct_indx].v2t;
//		n = pos[frm]->facet[fct_indx].nt;
//
//		/* Calculate and store the boundaries of this tile */
//		tile_i1 = offsetx;
//		tile_i2 = tile_i1 + tile_x;
//		tile_j1 = offsety;
//		tile_j2 = tile_j1 + tile_x;
//
//		/* Load this facet's boundaries and clamp them if needed, then
//		 * convert to local shared memory array addressing  */
//		i1 = max(pos[frm]->facet[fct_indx].ilim.x, tile_i1);
//		i2 = min(pos[frm]->facet[fct_indx].ilim.y, (tile_i2-1));
//		j1 = max(pos[frm]->facet[fct_indx].jlim.x, tile_j1);
//		j2 = min(pos[frm]->facet[fct_indx].jlim.y, (tile_j2-1));
//
//		/* Precalculate s and t components for the pixel loop */
//		double a, b, c, d, e, h, ti, tj, si, sj, si0, sj0, ti0, tj0, sz, tz, den, s, t, z, old;
//		a = i1*kmpxl - v0.x;
//		b = v2.y - v1.y;
//		c = v2.x - v1.x;
//		d = j1*kmpxl - v0.y;
//		e = v1.x - v0.x;
//		h = v1.y - v0.y;
//		den = e*b - c*h;
//		ti = -h*kmpxl/den;
//		tj = e*kmpxl/den;
//		si = b*kmpxl/den;
//		sj = -c*kmpxl/den;
//		si0 = (a*b - c*d)/den;
//		ti0 = (e*d -a*h)/den;
//		sz = v1.z - v0.z;
//		tz = v2.z - v1.z;
//
//		/* Now convert i1, i2, j1, j2 to shared-memory tile coordinates */
//		i1 -= (offsetx);
//		i2 -= (offsetx);
//		j1 -= (offsety);
//		j2 -= (offsety);
//
//		/* Pre-calculate some quantities for cosine smoothing if enabled */
//		if (posvis_tiled_smooth) {
//			/* Assign temp. normal components as float3 */
//			fidx.x = verts[0]->f[fct_indx].v[0];
//			fidx.y = verts[0]->f[fct_indx].v[1];
//			fidx.z = verts[0]->f[fct_indx].v[2];
//
//			tv0.x = verts[0]->v[fidx.x].n[0];	tv0.y = verts[0]->v[fidx.x].n[1];
//			tv0.z = verts[0]->v[fidx.x].n[2];	tv1.x = verts[0]->v[fidx.y].n[0];
//			tv1.y = verts[0]->v[fidx.y].n[1];	tv1.z = verts[0]->v[fidx.y].n[2];
//			tv2.x = verts[0]->v[fidx.z].n[0];	tv2.y = verts[0]->v[fidx.z].n[1];
//			tv2.z = verts[0]->v[fidx.z].n[2];
//			n1n0.x = tv1.x - tv0.x;	n1n0.y = tv1.y - tv0.y;	n1n0.z = tv1.z - tv0.z;
//			tv2.x -= tv1.x;	tv2.y -= tv1.y; tv2.z -= tv1.z;
//		}
//
//		/* Facet is at least partly within POS frame: find all POS
//		 * pixels whose centers project onto this facet  */
//		for (i=i1; i<=i2; i++) {
//
//			sj0 = si0;	/* Initialize this loop's base sj0, tj0 */
//			tj0 = ti0;
//
//			for (j=j1; j<=j2; j++) {
//
//				/* Calculate s and t parameters */
//				s = sj0;
//				t = tj0;
//
//				if ((s >= -SMALLVAL) && (s <= 1.0 + SMALLVAL)) {// &&
//
//					if(	(t >= -SMALLVAL) && (t <= s + SMALLVAL))	{
//
//						/* Compute z-coordinate of pixel center: its
//						 * distance measured from the origin towards
//						 * Earth.    */
//						z = v0.z + s*sz + t*tz;
//
//						/* Compare calculated z to stored shared memory z
//						 * array at this address and store the bigger value */
//						old = atomicMax64(&pos_z[i][j], z);
//
//						if (old < z){
//
//								if (posvis_tiled_smooth) {
//									/* Get pvs_smoothed version of facet unit
//									 * normal: Take the linear combination
//									 * of the three vertex normals; trans-
//									 * form from body to observer coordina-
//									 * tes; and make sure that it points
//									 * somewhat in our direction.         */
//									n.x = tv0.x + s * n1n0.x + t * tv2.x;
//									n.y = tv0.y + s * n1n0.y + t * tv2.y;
//									n.z = tv0.z + s * n1n0.z + t * tv2.z;
//									dev_cotrans1(&n, oa_sh, n, 1);
////									dev_cotrans3(&n, oa_gm, n, 1, frm);
//									dev_normalize3(&n);
//								}
//
//								/* Determine scattering angles.   */
//								if (n.z > 0.0) {
//									atomicExch((unsigned long long int*)&pos_cose[i][j],
//											__double_as_longlong(n.z));
//								}
//							/* Keeping track of facets may not be required.  */
////								atomicExch(&pos[frm]->f[i][j], f);
//
//						} /* end if (no other facet yet blocks this facet from view) */
//					} /* end if 0 <= t <= s (facet center is "in" this POS pixel) */
//				} /* end if 0 <= s <= 1 */
//
//				sj0 += sj;
//				tj0 += tj;
//			} /* end j-loop over POS rows */
//			/* Modify s and t step-wise for the next i-iteration of the pixel loop */
//			si0 += si;
//			ti0 += ti;
//
//		} /* end i-loop over POS columns */
//	}
//	__syncthreads();
//
//	/* Now write the shared memory array tiles into the global memory z buffer
//	 * and cosine array, again with a block-stride loop */
//	if (facet_index[frm][bin][0]!=0)
//		for (int index=threadIdx.x; index<tile_size; index+=blockDim.x) {
//			i = index % tile_x;
//			j = index / tile_x;
//			ig = i + offsetx;
//			jg = j + offsety;
//			if (pos_z[i][j]!=-1e20)
//				pos[frm]->z[ig][jg] = pos_z[i][j];
//			if (pos_cose[i][j]!=0.0)
//				pos[frm]->cose[ig][jg] = pos_cose[i][j];
//		}
//	__syncthreads();
//}
//
//__global__ void lightcurve_raster_krnl64(struct pos_t **pos,
//		struct vertices_t **verts,
//		double3 *oa,
//		double3 *usrc,
//		int ***facet_index,
//		int **entries,
//		int nf,
//		int frm,
//		int *n_tiles,
//		int *n_tiles_x,
//		int	*n_tiles_y,
//		int tile_size,
//		int tile_x,
//		int src) {
//
//	/* This kernel performs the rasterization tile by tile.  Each thread block
//	 * is responsible for one tile. */
//	/* Determine which tile this thread block is responsible for and
//	 * which element of the thread block this thread is. 	 */
//	int bin = blockIdx.x;
//	int index = threadIdx.x;
//
//	/* Declare the shared memory variables and others */
//	__shared__ double pos_z[32][32];//[55][55];	/* One per thread block */
//	__shared__ double pos_cose[32][32];//[55][55];	/* One per thread block */
//	__shared__ double pos_cosi[32][32];
//	__shared__ int2 bn;
//	__shared__ int xlim, ylim, offsetx, offsety, bistatic;
//	__shared__ double kmpxl;
//	int i, j, ig, jg, i1, i2, j1, j2;	/* ig,jg are global indices */
//	int tile_i1, tile_i2, tile_j1, tile_j2, fct_indx;
//	double3 v0, v1, v2, n, tv0, tv1, tv2, n1n0;
//	int3 fidx;
//
//	/* Initialize the shared memory arrays with grid-stride loop */
//	for (int index=threadIdx.x; index<tile_size; index+=blockDim.x) {
//		i = index % tile_x;
//		j = index / tile_x;
//		pos_z[i][j] = -1e20;
//		pos_cose[i][j] = 0.0;
//	}
//	__syncthreads();
//
//	/* Load variables used by every thread (per block) to shared memory for
//	 * faster access */
//	if (threadIdx.x==0) {
//		bistatic = pos[frm]->bistatic;
//		xlim = pos[frm]->xlim[0];
//		ylim = pos[frm]->ylim[0];
//		kmpxl = pos[frm]->km_per_pixel;
//		bn.x = bin % n_tiles_x[frm];
//		bn.y = bin / n_tiles_x[frm];
//
//		/* Calculate the pixel offsets needed to go back and forth between
//		 * tiled POS space for this block's tile and global POS space 	 */
//		offsetx = xlim + tile_x * bn.x;
//		offsety = ylim + tile_x * bn.y;
//	}
//	__syncthreads();
//
//	/* Using grid-stride loop, step through all facet entries for each bin where
//	 * each thread block is responsible for one bin/tile	 */
//	for (index=threadIdx.x; index<entries[frm][bin]; index+=blockDim.x) {
//
//		/* Load facet index into registers */
//		fct_indx = facet_index[frm][bin][index];
//
//		/* Load transformed facet vertices into registers */
//		v0 = pos[frm]->facet[fct_indx].v0t;
//		v1 = pos[frm]->facet[fct_indx].v1t;
//		v2 = pos[frm]->facet[fct_indx].v2t;
//		n = pos[frm]->facet[fct_indx].nt;
//
//		/* Calculate and store the boundaries of this tile */
//		tile_i1 = offsetx;
//		tile_i2 = tile_i1 + tile_x;
//		tile_j1 = offsety;
//		tile_j2 = tile_j1 + tile_x;
//
//		/* Load this facet's boundaries and clamp them if needed, then
//		 * convert to local shared memory array addressing  */
//		i1 = max(pos[frm]->facet[fct_indx].ilim.x, tile_i1);
//		i2 = min(pos[frm]->facet[fct_indx].ilim.y, (tile_i2-1));
//		j1 = max(pos[frm]->facet[fct_indx].jlim.x, tile_j1);
//		j2 = min(pos[frm]->facet[fct_indx].jlim.y, (tile_j2-1));
//
//		/* Pre-calculate s and t components for the pixel loop */
//		double a, b, c, d, e, h, ti, tj, si, sj, si0, sj0, ti0, tj0, sz, tz, den, s, t, z, old;
//		a = i1*kmpxl - v0.x;
//		b = v2.y - v1.y;
//		c = v2.x - v1.x;
//		d = j1*kmpxl - v0.y;
//		e = v1.x - v0.x;
//		h = v1.y - v0.y;
//		den = e*b - c*h;
//		ti = -h*kmpxl/den;
//		tj = e*kmpxl/den;
//		si = b*kmpxl/den;
//		sj = -c*kmpxl/den;
//		si0 = (a*b - c*d)/den;
//		ti0 = (e*d -a*h)/den;
//		sz = v1.z - v0.z;
//		tz = v2.z - v1.z;
//
//		/* Now convert i1, i2, j1, j2 to shared-memory tile coordinates */
//		i1 -= (offsetx);
//		i2 -= (offsetx);
//		j1 -= (offsety);
//		j2 -= (offsety);
//
//		/* Pre-calculate some quantities for cosine smoothing if enabled */
//		if (posvis_tiled_smooth) {
//			/* Assign temp. normal components as float3 */
//			fidx.x = verts[0]->f[fct_indx].v[0];
//			fidx.y = verts[0]->f[fct_indx].v[1];
//			fidx.z = verts[0]->f[fct_indx].v[2];
//
//			tv0.x = verts[0]->v[fidx.x].n[0];	tv0.y = verts[0]->v[fidx.x].n[1];
//			tv0.z = verts[0]->v[fidx.x].n[2];	tv1.x = verts[0]->v[fidx.y].n[0];
//			tv1.y = verts[0]->v[fidx.y].n[1];	tv1.z = verts[0]->v[fidx.y].n[2];
//			tv2.x = verts[0]->v[fidx.z].n[0];	tv2.y = verts[0]->v[fidx.z].n[1];
//			tv2.z = verts[0]->v[fidx.z].n[2];
//			n1n0.x = tv1.x - tv0.x;	n1n0.y = tv1.y - tv0.y;	n1n0.z = tv1.z - tv0.z;
//			tv2.x -= tv1.x;	tv2.y -= tv1.y; tv2.z -= tv1.z;
//		}
//
//		/* Facet is at least partly within POS frame: find all POS
//		 * pixels whose centers project onto this facet  */
//		for (i=i1; i<=i2; i++) {
//
//			sj0 = si0;	/* Initialize this loop's base sj0, tj0 */
//			tj0 = ti0;
//
//			for (j=j1; j<=j2; j++) {
//
//				/* Calculate s and t parameters */
//				s = sj0;
//				t = tj0;
//
//				if ((s >= -SMALLVAL) && (s <= 1.0 + SMALLVAL)) {// &&
//
//					if(	(t >= -SMALLVAL) && (t <= s + SMALLVAL))	{
//
//						/* Compute z-coordinate of pixel center: its
//						 * distance measured from the origin towards
//						 * Earth.    */
//						z = v0.z + s*sz + t*tz;
//
//						/* Compare calculated z to stored shared memory z
//						 * array at this address and store the bigger value */
//						old = atomicMax64(&pos_z[i][j], z);
//
//						if (old < z){
//								if (posvis_tiled_smooth) {
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
//								/* Determine scattering angles.   */
//								if (n.z > 0.0) {
//									if (src)
//										atomicExch((unsigned long long int*)&pos[frm]->cosill[i][j],
//												__double_as_longlong(n.z));
//									else
//										atomicExch((unsigned long long int*)&pos[frm]->cose[i][j],
//												__double_as_longlong(n.z));
//
//									if ((!src) && (bistatic)) {
//
//										double temp = dev_dot_d3(n,usrc[frm]);
//										atomicExch((unsigned long long int*)&pos[frm]->cosi[i][j],
//												__double_as_longlong(temp));
//										if (pos[frm]->cosi[i][j] <= 0.0)
//											pos[frm]->cose[i][j] = 0.0;
//									}
//
//								}
//							/* Keeping track of facets may not be required.  */
////								atomicExch(&pos[frm]->f[i][j], f);
//
//						} /* end if (no other facet yet blocks this facet from view) */
//					} /* end if 0 <= t <= s (facet center is "in" this POS pixel) */
//				} /* end if 0 <= s <= 1 */
//
//				sj0 += sj;
//				tj0 += tj;
//			} /* end j-loop over POS rows */
//			/* Modify s and t step-wise for the next i-iteration of the pixel loop */
//			si0 += si;
//			ti0 += ti;
//
//		} /* end i-loop over POS columns */
//	}
//	__syncthreads();
//
//	/* Now write the shared memory array tiles into the global memory z buffer
//	 * and cosine array, again with a block-stride loop */
//	if (facet_index[frm][bin][0]!=0)
//		for (int index=threadIdx.x; index<tile_size; index+=blockDim.x) {
//			i = index % tile_x;
//			j = index / tile_x;
//			ig = i + offsetx;
//			jg = j + offsety;
//			if (pos_z[i][j]!=-1e20)
//				pos[frm]->z[ig][jg] = pos_z[i][j];
//			if (pos_cose[i][j]!=0.0)
//				pos[frm]->cose[ig][jg] = pos_cose[i][j];
//		}
//	__syncthreads();
//}
//
//__global__ void posvis_outbnd_tiled_krnl64(struct pos_t **pos,
//		int *outbndarr, double4 *ijminmax_overall, int size, int start) {
//	/* nfrm_alloc-threaded kernel */
//	int posn, f = blockIdx.x * blockDim.x + threadIdx.x + start;
//	double xfactor, yfactor;
//	if (f <size) {
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
//__host__ int posvis_tiled_gpu64(
//		struct par_t *dpar,
//		struct mod_t *dmod,
//		struct dat_t *ddat,
//		struct pos_t **pos,
//		struct vertices_t **verts,
//		double3 orbit_offset,
//		int *posn,
//		int *outbndarr,
//		int set,
//		int nfrm_alloc,
//		int src,
//		int nf,
//		int body, int comp, unsigned char type, cudaStream_t *pv_stream,
//		int src_override) {
//
//	dim3 BLK,THD, BLKfrm, THD64, *BLKtile, BLKaf, THDaf;
//	double4 *ijminmax_overall;
//	double3 *oa, *usrc;
//	int ***facet_index, *xspan, *yspan, *n_tiles_x, *n_tiles_y, *n_tiles, **entries;
//	int f, outbnd, start, sharedMem, oasize, span, tile_size, **addr_index;
//
//	oasize=nfrm_alloc*3;
//	span=32;
//	tile_size=span*span;
//
//	/* Launch parameters for the facet_streams kernel */
//	THD.x = 256;	THD64.x = 64;
//	BLK.x = floor((THD.x - 1 + nf) / THD.x);
//	BLKfrm.x = floor((THD64.x - 1 + nfrm_alloc)/THD64.x);
//	THDaf.x = 1024;
//	BLKaf.x = nfrm_alloc;
//
//	/* Set up the offset addressing for lightcurves if this is a lightcurve */
//	if (type == LGHTCRV) 	start = 1;	/* fixes the lightcurve offsets */
//	else					start = 0;
//
//	/* Allocate temporary arrays/structs */
//	cudaCalloc1((void**)&ijminmax_overall, 	sizeof(double4), nfrm_alloc);
//	cudaCalloc1((void**)&oa, 				sizeof(double3), oasize);
//	cudaCalloc1((void**)&usrc, 				sizeof(double3), nfrm_alloc);
//	cudaCalloc1((void**)&xspan, 			sizeof(int), 	 nfrm_alloc);
//	cudaCalloc1((void**)&yspan, 			sizeof(int), 	 nfrm_alloc);
//	cudaCalloc1((void**)&n_tiles, 			sizeof(int), 	 nfrm_alloc);
//	cudaCalloc1((void**)&n_tiles_x, 		sizeof(int), 	 nfrm_alloc);
//	cudaCalloc1((void**)&n_tiles_y, 		sizeof(int), 	 nfrm_alloc);
//	/* Allocate the frame portion of the facet index triple pointer and
//	 * the bin entries counter */
//	cudaCalloc((void**)&facet_index,		sizeof(int**),	 nfrm_alloc);
//	cudaCalloc1((void**)&entries, 			sizeof(int*), 	 nfrm_alloc);
////	cudaCalloc1((void**)&addr_index,		sizeof(int*), 	 nfrm_alloc);
//	cudaCalloc((void**)&BLKtile, 			sizeof(dim3), 	 nfrm_alloc);
//
//	/* Initialize/pre-calculate values for rasterization */
//	posvis_tiled_init_krnl64<<<BLKfrm,THD64>>>(dpar, pos, ijminmax_overall, oa, usrc,
//			outbndarr, comp, start, src, nfrm_alloc, set, src_override);
//	checkErrorAfterKernelLaunch("posvis_tiled_init_krnl64");
//
//	/* Transform facet normals and determine bounding for facets and pos.  */
//	for (f=start; f<nfrm_alloc; f++)
//		transform_facet_normals_krnl64c<<<1,THD,0,pv_stream[f]>>>(dmod, pos,
//				verts, ijminmax_overall, orbit_offset, oa, outbndarr, nf,
//				f, src);
//	checkErrorAfterKernelLaunch("transform_facet_normals_krnl64a");
////	/* Transform facet normals and determine bounding for facets and pos.  */
////	for (f=start; f<nfrm_alloc; f++)
////		transform_facet_normals_krnl64b<<<BLK,THD,0,pv_stream[f]>>>(dmod, pos,
////				verts, ijminmax_overall, orbit_offset, oa, usrc, outbndarr, nf,
////				f, src);
////	checkErrorAfterKernelLaunch("transform_facet_normals_krnl64b");
////	transform_facet_normals_krnl64b<<<BLKaf,THD>>>(dmod, pos, verts, ijminmax_overall,
////			orbit_offset, oa, outbndarr, nf, src);
////	checkErrorAfterKernelLaunch("transform_facet_normals_krnl64a");
//
//	for (f=start; f<nfrm_alloc; f++)
//		cudaStreamSynchronize(pv_stream[f]);
//
//	/* Now calculate the tiling parameters to cover the POS view */
//	for (f=start; f<nfrm_alloc; f++) {
//		xspan[f] = pos[f]->xlim[1] - pos[f]->xlim[0] + 1;
//		yspan[f] = pos[f]->ylim[1] - pos[f]->ylim[0] + 1;
//		n_tiles_x[f] = (xspan[f]/span) + 1;
//		n_tiles_y[f] = (yspan[f]/span) + 1;
//		n_tiles[f] = n_tiles_x[f] * n_tiles_y[f];
//		BLKtile[f].x = n_tiles[f];
//		BLKtile[f].y = BLKtile[f].z = 1;
//
//		/* Now allocate the tiles section of the facet index and then step
//		 * through each tile section to allocate enough space for 1024
//		 * facet indices.  This is the maximum number of entries allowable
//		 * per thread block		 */
//		/* Allocate the entries array to keep track of how many facets each bin holds */
////		cudaCalloc((void**)&addr_index[f], 	sizeof(int), n_tiles[f]);
//		cudaCalloc((void**)&entries[f], 	sizeof(int), n_tiles[f]);
//		cudaCalloc((void**)&facet_index[f], sizeof(int*),n_tiles[f]);
//		for (int ti=0; ti<n_tiles[f]; ti++)
//			cudaCalloc((void**)&facet_index[f][ti], sizeof(int), 4*1024);
//	}
//
//	bin_facets_krnl64c<<<BLKaf,THDaf>>>(pos, verts, facet_index,
//			entries, nf, n_tiles, n_tiles_x, n_tiles_y, span);
//	checkErrorAfterKernelLaunch("bin_facets_krnl64");
//
//	/* Now we bin the triangles into the tiles */
//	for (f=start; f<nfrm_alloc; f++) {
////		sharedMem = sizeof(int)*n_tiles[f];
////		bin_facets_krnl64a<<<1,THD,sharedMem,pv_stream[f]>>>(pos, verts, facet_index,
////				entries, nf, f, n_tiles, n_tiles_x, n_tiles_y, span);
////		bin_facets_krnl64b<<<BLK,THD,sharedMem,pv_stream[f]>>>(pos, verts, facet_index,
////				entries, addr_index, nf, f, n_tiles, n_tiles_x, n_tiles_y, span);
//		radar_raster_krnl64<<<BLKtile[f],THD,0,pv_stream[f]>>>(pos,	verts, oa,
//				facet_index, entries, nf, f, n_tiles, n_tiles_x,
//				n_tiles_y, tile_size, span);
//	}
//	checkErrorAfterKernelLaunch("bin_facets_krnl64");
//
//	for (f=start; f<nfrm_alloc; f++)
//		cudaStreamSynchronize(pv_stream[f]);
//
//	/* Take care of any posbnd flags */
//	posvis_outbnd_tiled_krnl64<<<BLKfrm,THD64>>>(pos,
//			outbndarr, ijminmax_overall, nfrm_alloc, start);
//	checkErrorAfterKernelLaunch("posvis_outbnd_krnl64");
//	gpuErrchk(cudaMemcpyFromSymbol(&outbnd, posvis_tiled_outbnd, sizeof(int), 0,
//			cudaMemcpyDeviceToHost));
//
////	int n = 75;
////	int npixels = 151*151;
////	f = 0;
////	dbg_print_pos_arrays_full64(pos, 0, npixels, n);
////	dbg_print_pos_arrays_full64(pos, 1, npixels, n);
////	dbg_print_pos_arrays_full64(pos, 2, npixels, n);
////	dbg_print_pos_arrays_full64(pos, 3, npixels, n);
//
//	/* Free temp arrays, destroy streams and timers, as applicable */
//	cudaFree(oa);
//	cudaFree(usrc);
//	cudaFree(xspan);
//	cudaFree(yspan);
//	cudaFree(n_tiles);
//	cudaFree(entries);
//	cudaFree(BLKtile);
//	cudaFree(n_tiles_x);
//	cudaFree(n_tiles_y);
////	cudaFree(addr_index);
//	cudaFree(facet_index);
//	cudaFree(ijminmax_overall);
//
//	return outbnd;
//}
