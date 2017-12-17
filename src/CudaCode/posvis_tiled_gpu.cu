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

#define maxbins 100
__device__ int posvis_tiled_outbnd, posvis_tiled_smooth;
__device__ int dbg_cntr1=0;

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
		dbg_cntr1=0;
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
		int src,
		int blockSize)
{
	/* This kernel launches nf threads and transforms each facet normal with
	 * oa[frm] and stores the result back to dmod if n.z > 0.0.
	 * It also determines and stores the facet and global model bounding box
	 * via i1,i2,j1,j2 and xlim/ylim. These also get stored back to the dmod.
	 * The kernel counts the number of visible facets (where n.z > 0.0)
	 */

	/* Declare kernel variables */
	__shared__ int pn;
	__shared__ double kmpxl;
	int imin, jmin, imax, jmax, i1, i2, j1, j2;
	int3 fidx;
	double imin_dbl, jmin_dbl, imax_dbl, jmax_dbl;
	double3 n;
	double3 v0, v1, v2;

	if (threadIdx.x==0) {
		pn = pos[frm]->n;
		kmpxl = pos[frm]->km_per_pixel;
	}

	__syncthreads();

	/* Check f is within bounds */
	for (int f=threadIdx.x; f<nf; f+=blockSize) {
		if (frm==0) atomicAdd(&dbg_cntr1, 1);
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
//			atomicAdd(&dbg_cntr1, 1);
			/* First, store the transformed normal back to the model and increase
			 * visible facet counter */
			verts[0]->f[f].nt.x = n.x;	verts[0]->f[f].nt.y = n.y;	verts[0]->f[f].nt.z = n.z;
			//atomicAdd(&nvf[frm], 1);

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

//			if (j2!=0) printf("j2=%i\n", j2);

			verts[0]->f[f].ilim.x = i1;
			verts[0]->f[f].jlim.x = j1;
			verts[0]->f[f].ilim.y = i2;
			verts[0]->f[f].jlim.y = j2;
			verts[0]->f[f].v0t.x = v0.x;	verts[0]->f[f].v0t.y = v0.y;	verts[0]->f[f].v0t.z = v0.z;
			verts[0]->f[f].v1t.x = v1.x;	verts[0]->f[f].v1t.y = v1.y;	verts[0]->f[f].v1t.z = v1.z;
			verts[0]->f[f].v2t.x = v2.x;	verts[0]->f[f].v2t.y = v2.y;	verts[0]->f[f].v2t.z = v2.z;

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
//			printf("facet [%i] i1=%i, i2=%i, j1=%i, j2=%i\n",
//					verts[0]->f[f].ilim.x, verts[0]->f[f].ilim.y, verts[0]->f[f].jlim.x, verts[0]->f[f].jlim.y);

		}
		else {
			/* The following makes a check in the bin_facets_krnl64 kernel easier */
			verts[0]->f[f].nt.x = -1.0;
			verts[0]->f[f].nt.y = -1.0;
			verts[0]->f[f].nt.z = -1.0;
		}
	}
}

__global__ void bin_facets_krnl64(struct pos_t **pos,
		struct vertices_t **verts,
		int ***facet_index,
		int **entries,
		int nf,
		int frm,
		int *n_tiles,
		int *n_tiles_x,
		int	*n_tiles_y,
		int tile_size)
{
	/* This kernel is responsible for binning visible model facets according to
	 * which screen tile they appear on. Each facet can belong to 1, 2, or 4
	 * different facets.  (If the size of individual triangles should exceed
	 * the tile size, this is no longer true.)
	 * The kernel has just one thread block with 1024 threads. It uses a grid-
	 * stride loop to cover all facets
	 */
	int f, current_i, next_i, current_j, next_j, i1, i2, j1, j2, bi, bin, old_indx;
	__shared__ int2 xlim, ylim; /* These are the global pos limits */
	extern __shared__ int addr_index[];	/* Used for the facet_index entries */
//	__shared__ int addr_index[maxbins];
	if (threadIdx.x==0) {
		xlim.x = pos[frm]->xlim[0];
		xlim.y = pos[frm]->xlim[1];
		ylim.x = pos[frm]->ylim[0];
		ylim.y = pos[frm]->ylim[1];
		for (bin=0; bin<n_tiles[frm]; bin++)
			addr_index[bin] = 0;
	}
	__syncthreads();

	/* Check that the thread/facet number is smaller than the # of facets
	 * and that it's a visible facet  */
	for (f=threadIdx.x; f<nf; f+=blockDim.x) {
		/* Weed out any facets not visible to observer */
		if (verts[0]->f[f].nt.z > 0.0) {

			bi = 0;	/* Bin index for the four facet bin entries*/
			/* Copy facet limits into register memory for faster access */
			i1 = verts[0]->f[f].ilim.x;
			i2 = verts[0]->f[f].ilim.y;
			j1 = verts[0]->f[f].jlim.x;
			j2 = verts[0]->f[f].jlim.y;

			/* Now check where the current facet lies, stepping through each
			 * tile */
			for (int k=0; k<n_tiles_y[frm]; k++) {
				current_j = ylim.x + k * tile_size;
				next_j = current_j + tile_size;
				for (int n=0; n<n_tiles_x[frm]; n++) {
					bin = k*n_tiles_x[frm] + n;
					current_i = xlim.x + n * tile_size;
					next_i = current_i + tile_size;

					/* If i1 or i2 AND j1 or j2 fall into this tile, register it */
					if ((i1>=current_i && i1<next_i)  || (i2>=current_i && i2<next_i)) {
						if ((j1>=current_j && j1<next_j)  || (j2>=current_j && j2<next_j)) {
							verts[0]->f[f].bin[bi] = bin;
							old_indx = atomicAdd(&addr_index[bin], 1);
							facet_index[frm][bin][old_indx] = f;
							atomicAdd(&entries[frm][bin], 1);
							bi++;

						}
					}
				}
			}
		}
	}
}

__global__ void radar_raster_krnl64(struct pos_t **pos,
		struct vertices_t **verts,
		int ***facet_index,
		int **entries,
		int nf,
		int frm,
		int *n_tiles,
		int *n_tiles_x,
		int	*n_tiles_y,
		int tile_size,
		int specific_span) {

	/* This kernel performs the rasterization tile by tile.  Each thread block
	 * is responsible for one tile. */

	/* Declare the shared memory arrays and initialize with block-stride loop,
	 * then synchronize all threads in the current thread block */
	__shared__ double pos_z[32][32];//[55][55];	/* One per thread block */
	__shared__ double pos_cose[32][32];//[55][55];	/* One per thread block */
	int i, j, ig, jg, i1, i2, j1, j2;	/* ig,jg are global indices */
	int tile_i1, tile_i2, tile_j1, tile_j2, fct_indx;
	double3 v0, v1, v2, n;

	for (int index=threadIdx.x; index<tile_size; index+=blockDim.x) {
		i = index % specific_span;
		j = index / specific_span;
		pos_z[i][j] = -1e20;
		pos_cose[i][j] = 0.0;
	}
	__syncthreads();


	/* Determine which tile this thread block is responsible for and
	 * which element of the thread block this thread is. 	 */
	int bin = blockIdx.x;
	int index = threadIdx.x;
	int2 bn; bn.x = bin % n_tiles_x[frm]; bn.y = bin / n_tiles_x[frm];
	__shared__ int xlim, ylim;
	__shared__ double kmpxl;

	/* Load the pos limits to shared memory for faster access */
	if (threadIdx.x==0) {
		xlim = pos[frm]->xlim[0];
		ylim = pos[frm]->ylim[0];
		kmpxl = pos[frm]->km_per_pixel;
	}
	__syncthreads();
	/* Check that we are within bounds on both bin counter and # of entries in
	 * that bin	*/
	if (bin < n_tiles[frm]) {
		for (index=threadIdx.x; index<entries[frm][bin]; index+=blockDim.x) {

			/* Load facet index into registers */
			fct_indx = facet_index[frm][bin][index];
//			printf("fct_index in thread %i, bin %i, index %i = %i\n", threadIdx.x, bin, index, fct_indx);

			/* Load transformed facet vertices into registers */
			v0 = verts[0]->f[fct_indx].v0t;
			v1 = verts[0]->f[fct_indx].v1t;
			v2 = verts[0]->f[fct_indx].v2t;
			n  = verts[0]->f[fct_indx].nt;

			/* Calculate and store the boundaries of this tile */
			tile_i1 = xlim + bn.x * specific_span;
			tile_i2 = tile_i1 + specific_span;
			tile_j1 = ylim + bn.y * specific_span;
			tile_j2 = tile_j1 + specific_span;

			/* Load this facet's boundaries and clamp them if needed, then
			 * convert to local shared memory array addressing  */
			i1 = max(verts[0]->f[fct_indx].ilim.x, tile_i1);
			i2 = min(verts[0]->f[fct_indx].ilim.y, tile_i2);
			j1 = max(verts[0]->f[fct_indx].jlim.x, tile_j1);
			j2 = min(verts[0]->f[fct_indx].jlim.y, tile_j2);

//			if (frm==0) {
//				printf("facet_i1=%i, facet_i2=%i, facet_j1=%i, facet_j2=%i, i1=%i, i2=%i, j1=%i, j2=%i, tile_i1=%i, tile_i2=%i, tile_j1=%i, tile_j2=%i\n",
//						verts[0]->f[fct_indx].ilim.x, verts[0]->f[fct_indx].ilim.y, verts[0]->f[fct_indx].jlim.x, verts[0]->f[fct_indx].jlim.y,
//						i1, i2, j1, j2, tile_i1, tile_i2, tile_j1, tile_j2);
//			}

			/* Precalculate s and t components for the pixel loop */
			double a, b, c, d, e, h, ti, tj, si, sj, si0, sj0, ti0, tj0, sz, tz, den, s, t, z, old;
//			int pxa;
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

			/* Now convert i1, i2, j1, j2 to shared-memory tile coordinates */
			i1 -= (xlim + specific_span * bn.x);
			i2 -= (xlim + specific_span * bn.x);
			j1 -= (ylim + specific_span * bn.y);
			j2 -= (ylim + specific_span * bn.y);

			/* Facet is at least partly within POS frame: find all POS
			 * pixels whose centers project onto this facet  */
			for (i=i1; i<=i2; i++) {

				sj0 = si0;	/* Initialize this loop's base sj0, tj0 */
				tj0 = ti0;

				for (j=j1; j<=j2; j++) {

					/* Calculate local pixel address for shared memory arrays */
					//pxa = (j+pn) * span + (i+pn);
					s = sj0;
					t = tj0;

					if ((s >= -SMALLVAL) && (s <= 1.0 + SMALLVAL)) {// &&
//						if (frm==0) atomicAdd(&dbg_cntr1,1);
						if(	(t >= -SMALLVAL) && (t <= s + SMALLVAL))	{

							/* Compute z-coordinate of pixel center: its
							 * distance measured from the origin towards
							 * Earth.    */
							z = v0.z + s*sz + t*tz;
//printf("z in facet %i = %3.8g\n", fct_indx, z);
							/* Compare calculated z to stored shared memory z
							 * array at this address and store the bigger value */
							old = atomicMax64(&pos_z[i][j], z);

							if (old < z){

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
//									dev_cotrans8(&n, oa, n, 1, frm);
//									dev_normalize2(&n);
//								}
//
								/* Determine scattering angles.   */
//								if (n.z > 0.0) {
//									if(bin==0 )
//										printf("n.z in facet %i = %3.8g\n", fct_indx, n.z);
//									pos_cose[i][j]=0.0;
//									atomicExch((unsigned long long int*)&pos_cose[i][j],
//											__double_as_longlong(n.z));
//								}

								/* Keeping track of facets may not be required.  */
//								atomicExch(&pos[frm]->f[i][j], f);

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
		}
	}
	__syncthreads();

	/* Now write the shared memory array tiles into the global memory z buffer
	 * and cosine array, again with a block-stride loop */
	for (int index=threadIdx.x; index<tile_size; index+=blockDim.x) {

		i = index % specific_span;
		j = index / specific_span;
		ig = i + xlim + specific_span * bn.x;
		jg = j + ylim + specific_span * bn.y;
		pos[frm]->z[ig][jg] = pos_z[i][j];
		pos[frm]->cose[ig][jg] = pos_cose[i][j];

//		if (frm==0) {
//			atomicAdd(&dbg_cntr1,1);
//			printf("i=%i, ig=%i, j=%i, jg=%i, pos_z[i][j]=%3.8g\n", i, ig, j, jg, pos_z[i][j]);
//		}
	}
	__syncthreads();
}

__global__ void posvis_outbnd_tiled_krnl64(struct pos_t **pos,
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

__global__ void print_dbg_cntr_krnl(int flag) {

	if (threadIdx.x==0) {
		if (flag==0) {
			printf("# of facets with n.z > 0.0 is %i\n", dbg_cntr1);
			dbg_cntr1=0;
		}
		else
			printf("dbg_cntr=%i\n", dbg_cntr1);
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

	int f, outbnd, smooth, start, overlap, specific_span, bin;
	dim3 BLK,THD, BLKfrm, THD64, *BLKtile;
	double4 *ijminmax_overall;
	double3 *oa, *usrc;
	int *nvf, *xspan, *yspan, *n_tiles_x, *n_tiles_y, *n_tiles, **entries;
	int oasize = nfrm_alloc*3;
	int sharedMem;

	/* The following triple pointer is used to store model facet indices. They
	 * are organized by tiles/bins and frames. This is essential for the shared
	 * memory bucket rasterization	 */
	int ***facet_index;	/* will be addressed facet_index[frame][bin][index]

	/* To-Do:  Calculate these spans at program launch from max shared memory
	 * per thread block	 */
	int span_r64 = 32;//55;		/* These four spans are the specific maximum tile */
	int span_r32 = 32;//78;		/* sides depending on FP32/FP64 mode and data     */
	int span_lc64 = 45;		/* type - lightcurves need one more pos array     */
	int span_lc32 = 63;		/* than radar.									  */
	int tile_size = span_r64*span_r64;

	/* Launch parameters for the facet_streams kernel */
	THD.x = 256;	THD64.x = 64;
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
	//cudaCalloc1((void**)&nvf, 				sizeof(int), 	 nfrm_alloc);
	cudaCalloc1((void**)&xspan, 			sizeof(int), 	 nfrm_alloc);
	cudaCalloc1((void**)&yspan, 			sizeof(int), 	 nfrm_alloc);
	cudaCalloc1((void**)&n_tiles, 			sizeof(int), 	 nfrm_alloc);
	cudaCalloc1((void**)&n_tiles_x, 		sizeof(int), 	 nfrm_alloc);
	cudaCalloc1((void**)&n_tiles_y, 		sizeof(int), 	 nfrm_alloc);
	/* Allocate the frame portion of the facet index triple pointer and
	 * the bin entries counter */
	cudaCalloc((void**)&facet_index,		sizeof(int**),	 nfrm_alloc);
	cudaCalloc1((void**)&entries, 			sizeof(int*), 	 nfrm_alloc);
	cudaCalloc((void**)&BLKtile, 			sizeof(dim3), 	 nfrm_alloc);

	/* Initialize/pre-calculate values for rasterization */
	posvis_tiled_init_krnl64<<<BLKfrm,THD64>>>(dpar, pos, ijminmax_overall, oa, usrc,
			outbndarr, comp, start, src, nfrm_alloc, set, src_override);
	checkErrorAfterKernelLaunch("posvis_tiled_init_krnl64");

	/* Transform facet normals and determine bounding for facets and pos.  */
	for (f=start; f<nfrm_alloc; f++)
		/* Now the main facet kernel */
		transform_facet_normals_krnl64<<<1,THD,0,pv_stream[f]>>>(dmod, pos,
				verts, ijminmax_overall, orbit_offset, oa, usrc, outbndarr, nf,
				nvf, f, src, THD.x);

	checkErrorAfterKernelLaunch("transform_facet_normals_krnl64");
	cudaDeviceSynchronize();

//	print_dbg_cntr_krnl<<<1,1>>>(0);

	/* Now calculate the tiling parameters to cover the POS view */
	for (f=start; f<nfrm_alloc; f++) {
//		maxentries[f] = 0;
		xspan[f] = pos[f]->xlim[1] - pos[f]->xlim[0] + 1;
		yspan[f] = pos[f]->ylim[1] - pos[f]->ylim[0] + 1;
		n_tiles_x[f] = (xspan[f]/specific_span) + 1;
		n_tiles_y[f] = (yspan[f]/specific_span) + 1;
		n_tiles[f] = n_tiles_x[f] * n_tiles_y[f];
		BLKtile[f].x = n_tiles[f];
		BLKtile[f].y = BLKtile[f].z = 1;

		/* Now allocate the tiles section of the facet index and then step
		 * through each tile section to allocate enough space for 1024
		 * facet indices.  This is the maximum number of entries allowable
		 * per thread block		 */
		/* Allocate the entries array to keep track of how many facets each bin holds */
		cudaCalloc((void**)&entries[f], sizeof(int), n_tiles[f]);
		cudaCalloc((void**)&facet_index[f], sizeof(int*), n_tiles[f]);
		for (int ti=0; ti<n_tiles[f]; ti++)
			cudaCalloc((void**)&facet_index[f][ti], sizeof(int), 4*1024);
	}

	/* Now we bin the triangles into the tiles */
	for (f=start; f<nfrm_alloc; f++) {
		sharedMem = sizeof(int)*n_tiles[f];

		bin_facets_krnl64<<<1,THD,sharedMem,pv_stream[f]>>>(pos, verts, facet_index,
				entries, nf, f, n_tiles, n_tiles_x, n_tiles_y, specific_span);
		radar_raster_krnl64<<<BLKtile[f],THD,0,pv_stream[f]>>>(pos,
				verts, facet_index, entries, nf, f, n_tiles, n_tiles_x,
				n_tiles_y, tile_size, specific_span);
	}
	checkErrorAfterKernelLaunch("bin_facets_krnl64");

	cudaDeviceSynchronize();
	print_dbg_cntr_krnl<<<1,1>>>(1);
	checkErrorAfterKernelLaunch("print_dbg_cntr_krnl");

	/* Now check that there are not more than 1024 facet entries per bin.
	 * Warning only for now.
	 * Also determine the maximum number of facets in any bin in a frame */
//	for (f=start; f<nfrm_alloc; f++) {
//		for (bin=0; bin<n_tiles[f]; bin++) {
//			maxentries[f] = max(maxentries[f], entries[f][bin]);
//			printf("%i facets in bin %i in set[%i] frame[%i]\n", entries[f][bin], bin, set, f);
//		}
//		THDtile[f].x = maxentries[f];
//	}


//	/* Now the main rasterization kernel.  */
//	for (f=start; f<nfrm_alloc; f++) {
//		radar_raster_krnl64<<<BLKtile[f],THDtile[f],0,pv_stream[f]>>>(pos,
//				verts, facet_index, entries, nf, f, n_tiles, n_tiles_x,
//				n_tiles_y, tile_size, specific_span);
//	}
//	checkErrorAfterKernelLaunch("radar_raster_krnl_64");
//
//	cudaDeviceSynchronize();

//	posvis_facet_krnl64<<<BLK,THD, 0, pv_stream[1]>>>(pos, verts,
//					ijminmax_overall, orbit_offset, oa, usrc,	src, body, comp,
//					nf, 1, smooth, outbndarr, set);
//	dbg_krnl_psvs<<<1,1>>>();
//	gpuErrchk(cudaMemcpy(dbg_hn, dbg_n, sizeof(double3)*nf, cudaMemcpyDeviceToHost));
//	dbg_print_facet_normals_dbl3(dbg_hn, nf, "FP64_nrmls.csv");

	/* Take care of any posbnd flags */
	posvis_outbnd_tiled_krnl64<<<BLKfrm,THD64>>>(pos,
			outbndarr, ijminmax_overall, nfrm_alloc, start);
	checkErrorAfterKernelLaunch("posvis_outbnd_krnl64");
	gpuErrchk(cudaMemcpyFromSymbol(&outbnd, posvis_tiled_outbnd, sizeof(int), 0,
			cudaMemcpyDeviceToHost));

	int n = 75;
	int npixels = 151*151;
	f = 0;
	dbg_print_pos_arrays_full64(pos, 0, npixels, n);
	dbg_print_pos_arrays_full64(pos, 1, npixels, n);
	dbg_print_pos_arrays_full64(pos, 2, npixels, n);
	dbg_print_pos_arrays_full64(pos, 3, npixels, n);

	/* Free temp arrays, destroy streams and timers, as applicable */

	cudaFree(ijminmax_overall);
	cudaFree(oa);
	cudaFree(usrc);

	cudaFree(facet_index);
	cudaFree(entries);
	cudaFree(BLKtile);

	return outbnd;
}
