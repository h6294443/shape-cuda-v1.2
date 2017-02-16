/*****************************************************************************************
                                                                                 posvis.c

Fill in the portion of a plane-of-sky image due to a particular model component: Assign
each relevant POS pixel a z-value in observer coordinates (distance from the origin
towards Earth) and a value of cos(scattering angle).

Return 1 if any portion of this component lies outside the specified POS window,
0 otherwise.

If the "src" argument is true, the "observer" is the Sun rather than Earth, and
"plane-of-sky" becomes "projection as viewed from the Sun."

This is CUDA-enabled (and CUDA-only!) version of posvis
 *****************************************************************************************/
extern "C" {
#include "../shape/head.h"
#include "../shape/shape-cuda.h"
//#include <limits.h>
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

__device__ void dev_POSrect( struct pos_t *pos, int xylim[4], int xylim2[4], int src,
		double imin_dbl, double imax_dbl, double jmin_dbl, double jmax_dbl,
		double *max_overall, int n)
{
	/* Keep track of the changed POS region. What follows is the code for the
	 * POSrect() function, which cannot be called from this kernel.
	 * Update the POS region that contains the target without regard whether or
	 * not it extends beyond the POS frame*/

	int imin, imax, jmin, jmax;

	/*  Update POS region that contains target without regard to whether or not
	 * 	it extends beyond the POS frame  */
	max_overall[0] = MIN(max_overall[0], imin_dbl);
	max_overall[1] = MAX(max_overall[1], imax_dbl);
	max_overall[2] = MIN(max_overall[2], jmin_dbl);
	max_overall[3] = MAX(max_overall[3], jmax_dbl);

	/*  Update the subset of the POS frame that contains the target  */
	imin = (imin_dbl < INT_MIN) ? INT_MIN : (int) imin_dbl;
	imax = (imax_dbl > INT_MAX) ? INT_MAX : (int) imax_dbl;
	jmin = (jmin_dbl < INT_MIN) ? INT_MIN : (int) jmin_dbl;
	jmax = (jmax_dbl > INT_MAX) ? INT_MAX : (int) jmax_dbl;

	imin = MAX(imin, -n);
	imax = MIN(imax,  n);
	jmin = MAX(jmin, -n);
	jmax = MIN(jmax,  n);
//
//	if (src) {
//		if (imin < pos->xlim2[0])
//			pos->xlim2[0] = MAX( imin, -n);
//		if (imax > pos->xlim2[1])
//			pos->xlim2[1] = MIN( imax,  n);
//		if (jmin < pos->ylim2[0])
//			pos->ylim2[0] = MAX( jmin, -n);
//		if (jmax > pos->ylim2[1])
//			pos->ylim2[1] = MIN( jmax,  n);
//	} else {
//		if (imin < pos->xlim[0])
//			pos->xlim[0] = MAX( imin, -n);
//		if (imax > pos->xlim[1])
//			pos->xlim[1] = MIN( imax,  n);
//		if (jmin < pos->ylim[0])
//			pos->ylim[0] = MAX( jmin, -n);
//		if (jmax > pos->ylim[1])
//			pos->ylim[1] = MIN( jmax,  n);
//	}

	if (src) {
//		if (imin < xylim2[0])	xylim2[0] = MAX( imin, -n);
//		if (imax > xylim2[1])	xylim2[1] = MIN( imax,  n);
//		if (jmin < xylim2[2])	xylim2[2] = MAX( jmin, -n);
//		if (jmax > xylim2[3])	xylim2[3] = MIN( jmax,  n);
		atomicMin(&pos->xlim2[0], imin);
		atomicMax(&pos->xlim2[1], imax);
		atomicMin(&pos->ylim2[0], jmin);
		atomicMax(&pos->ylim2[1], jmax);
	} else {
		atomicMin(&pos->xlim[0], imin);
		atomicMax(&pos->xlim[1], imax);
		atomicMin(&pos->ylim[0], jmin);
		atomicMax(&pos->ylim[1], jmax);
//		if (imin < xylim[0])	xylim[0] = imin;
//		if (imax > xylim[1])	xylim[1] = imax;
//		if (jmin < xylim[2])	xylim[2] = jmin;
//		if (jmax > xylim[3])	xylim[3] = jmax;
	}
}


/*	The following __device__ kernel is associated with the posvis_cuda() routine and the
 * 	__global__ posvis_facets kernel.  The kernel accomplishes roughly the same thing as
 * 	the two nested for-loops (i1 < i < i2, j1 < j < j2) in the original posvis single CPU
 * 	routine.
 * 	This kernel is called from another kernel - posvis_facets.  As such, it implements
 * 	dynamic parallelism and requires CUDA Compute Capability 3.5 to work.
 *
 * 	Input arguments:
 * 		*verts - 		Pointer to address of vertices substructure of *mod.
 * 		*pos - 			Pointer to address of pos substructure of *dat.
 * 		npixels - 		Total number of pixels
 * 		body - 			current body #
 * 		comp - 			current component #
 * 		v0x,v0y,v0z - 	vertex 0
 * 		v1x,v1y,v1z - 	vertex 1
 * 		v2x,v2y,v2z - 	vertex 2
 * 		nx,ny,nz - 		normal vextor
 * 		imin_overall, imax_overall, jmin_overall, jmax_overall - used for POS sub-frame
 * 			boundaries.
 * 		smooth - 	the smoothing flag
 * 		src - 			Source flag.  If src=1, obs = sun.
 * 		**bod - 		Double pointer to current body of model, bod[i][j].
 * 		**cmp - 		Double pointer to current component of model, mod[i][j].
 * 		**fac - 		Double pointer to current facet of component, f[i][j].
 * 		**zz - 			Double pointer to obs z coord backwards along LOS, pos->z[i][j].
 * 		facet - 		Facet number for current facet.
 * 		i1,j1 - 		pos window starting coordinates
 * 		fno - 			facet number for this pos pixel
 * 		*oa - 			the oa[3][3] matrix created in the calling function
 * 		**cosa - 		Double pointer to incident angle cose[i][j].
 * 		**cosb - 		Double pointer to departing angle cosi[i][j].
 * 		usrc[3] - 		unit vector towards source		*/

__global__ void  posvis_pixels(struct vertices_t *verts, struct pos_t *pos,
		int npixels, int body, int comp,
		double v0x, double v0y, double v0z,
		double v1x, double v1y, double v1z,
		double v2x, double v2y, double v2z,
		double nx,  double ny,  double nz,
		double *max_overall, double *oa, double **cosa, double **cosb, double usrc[3],
		int smooth, int src, int **bod,	int **cmp, int **fac, double **zz,
		int i1, int j1, int fno)
{
	/*	Find the facet index with thread and block indices.  First are global
	 * 	thread indices in x and y directions.  Then total offset.  Finally i
	 * 	and j which index i1 <= i <= i2 and j1 <= j <= j2	*/
	int x_index = blockIdx.x * blockDim.x + threadIdx.x;
	int y_index = blockIdx.y * blockDim.y + threadIdx.y;
	int offset = y_index * blockDim.x + x_index;
	int i = x_index + i1;		// gets row index of POS space
	int j = y_index + j1;		// gets column index of POS space
	double s, t, z, den;		// the surface parameters plus a helper var (den)
	int k;						// for loop iterator
	double x[3];
	double v0[3] = {v0x, v0y, v0z};	// reconstruct the three vertex coordinates
	double v1[3] = {v1x, v1y, v1z};	// from the constituent parts.
	double v2[3] = {v2x, v2y, v2z};
	double n[3] = {nx, ny, nz};		// same for the normal vector for this facet

	if (offset < npixels){
		/*	Check if j = j1 which indicates the 0th column which means an extra computation
		 * 	step for just that thread where j = j1	*/
		if (j == j1)	x[0] = i*pos->km_per_pixel;		// determine physical x loc

		x[1] = j*pos->km_per_pixel;						// determine physical y loc

		/*	Now compute parameters s(x,y) and t(x,y) which define a facet surface:
		 *
		 * 			z = z0 + s*(z1-z0) + t*(z2-z1)
		 *
		 * 	where z0, z1, and z2 are the z-coordinates at the vertices.  The
		 * 	conditions 0 <= s <= 1 and 0 <= t <= s require the POS pixel center
		 * 	to be 'within' the (projected) perimeter of facet f.				*/

		den = 1/( (v1[0]-v0[0])*(v2[1]-v1[1]) - (v2[0]-v1[0])*(v1[1]-v0[1]) );
		s = ( (x[0]-v0[0]) * (v2[1]-v1[1]) - (v2[0]-v1[0]) * (x[1]-v0[1]) ) * den;

		/*	Now check if 0 <= s <= 1 */
		if ( (s >= -SMALLVAL) && (s <= 1.0+SMALLVAL) ) {

			t = ( (v1[0]-v0[0])*(x[1]-v0[1]) - (x[0]-v0[0])*(v1[1]-v0[1]) ) * den;

			/*Now check if 0 <= t <= s*/
			if ( (t >= -SMALLVAL) && (t <= s+SMALLVAL) ) {

				/*	Compute the z coordinate of the center of this pixel: its distance measured
				 * 	from the origin towards Earth.*/
				z = v0[2] + s*(v1[2] - v0[2]) + t*(v2[2] - v1[2]);

				/*  If fac[i][j] is >= 0, pixel [i][j] was already assigned
				 *  values during a previous call to posvis for a different
				 *  model component.  If so, override only if the current
				 *  component is blocking our view of (i.e., is closer to us
				 *  than) the previous one.*/
				if ( (z > zz[i][j]) || (fac[i][j] < 0) ) {

					/*	Next line assigns distance of POS pixel center from
					 * 	COM towards Earth; that is, by changing zz, it changes
					 * 	pos->z or pos->zill */
					zz[i][j] = z;

					if (smooth) {

						/*  Get smoothed version of this facet's unit normal:
				        Take the appropriate linear combination of the
				        three vertex normals; transform from body to
				        observer coordinates; and make sure that it
				        points somewhat in our direction.*/

						for (k=0; k<=2; k++)
							n[k] = verts->v[verts->f[fno].v[0]].n[k]
							                                      + s*(verts->v[verts->f[fno].v[1]].n[k]- verts->v[verts->f[fno].v[0]].n[k])
							                                      + t*(verts->v[verts->f[fno].v[2]].n[k]- verts->v[verts->f[fno].v[1]].n[k]);

						/*	Perform coordinate transformation n = oa * n
						 * 	This code performs the same as cotrans(), which
						 * 	cannot be called from a CUDA kernel.		*/
						dev_cotrans1(n, oa, n, 1);
						dev_normalize(n);
					}

					/*	Determine scattering and/or incidence angles - The following
					 * 	lines change cosa, thus changing pos->cose or pos->cosill.
					 * 	For bistatic situations (lightcurves) where we are viewing from
					 * 	Earth (src = 0), they also change cosb,thus changing pos->cosi.*/
					if (n[2] > 0.0) {
						cosa[i][j] = n[2];
						if ((!src) && (pos->bistatic)) {
							cosb[i][j] = dev_dot( n, usrc);
							if (cosb[i][j] <= 0.0)
								cosa[i][j] = 0.0;
						}
					}

					//dev_POSrect(pos, src, (double) i, (double) i, (double) j, (double) j, max_overall);

					/*  Next lines change pos->body or pos->bodyill, pos->comp or pos->compill,
					 *  and pos->f or pos->fill: the body, component, and facet numbers at the
					 *  center of this POS pixel.*/

					bod[i][j] = body;
					cmp[i][j] = comp;
					fac[i][j] = offset;		// f in original code
				}   //end if (no other facet yet blocks this facet from view)*/
			}   //end if 0 <= t <= s (facet center is "in" this POS pixel)
		}   //end if 0 <= s <= 1*/
	}
	__syncthreads();
}


__global__ void posvis_facets(struct vertices_t *verts, struct pos_t *pos,
		int nf, double *oa, double *orbit_offset, double *max_overall,
		int smooth, int body, int comp, int **bod, int **cmp, int **fac, int *outbnd, int src,
		int maxThreadsPerBlock, double **cosa, double **cosb, double **zz, double usrc[3])
{
	/*	Find the facet index with thread and block indices.	Offset is calculated for a one-
	 * 	dimensional thread string. In other words, this kernel should be launched with (Bx,0,0)
	 * 	and (Tx,0,0)*/
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int imin, imax, jmin, jmax, i, i1, i2, j1, j2;
	double imin_dbl, imax_dbl, jmin_dbl, jmax_dbl, v0[3], v1[3], v2[3], n[3];
	int fac_ind0 = verts->f[offset].v[0];
	int fac_ind1 = verts->f[offset].v[1];
	int fac_ind2 = verts->f[offset].v[2];

	if (offset < nf) {
		/*	Get normal for this facet (in body-fixed ast coord.) and
		 * 	transform it to obs coord.*/

		for (i=0; i<=2; i++)
			n[i] = verts->f[offset].n[i];

		dev_cotrans1( n, oa, n, 1);

		/*	Consider facet further only if the z-component of the normal is >0.
		 * 	This ensures it is sort-of pointed at the observer.*/
		if (n[2] > 0.0)
		{
			/*	First convert the three vertex coordinates that make up this facet
			 * 	from body to observer coordinates.  Then add the orbit offset
			 * 	(center-of-mass offset) in obs coord for this model at this frame's
			 * 	epoch due to orbital motion, in case the model is half of a binary
			 * 	system.*/

			dev_cotrans1( v0, oa, verts->v[fac_ind0].x, 1);
			dev_cotrans1( v1, oa, verts->v[fac_ind1].x, 1);
			dev_cotrans1( v2, oa, verts->v[fac_ind2].x, 1);

			for (i=0; i<=2; i++) {
				v0[i] += orbit_offset[i];
				v1[i] += orbit_offset[i];
				v2[i] += orbit_offset[i];
			}

			/*	Find the rectangular region (in POS pixels) that contains the
			 * 	projected facet.  Uses doubles in case the model has illegal
			 * 	parameters and the pixel numbers exceed the limits for valid
			 * 	integers.*/

			imin_dbl = floor(MIN(v0[0],MIN(v1[0],v2[0]))/pos->km_per_pixel - SMALLVAL+0.5);
			imax_dbl = floor(MAX(v0[0],MAX(v1[0],v2[0]))/pos->km_per_pixel + SMALLVAL+0.5);
			jmin_dbl = floor(MIN(v0[1],MIN(v1[1],v2[1]))/pos->km_per_pixel - SMALLVAL+0.5);
			jmax_dbl = floor(MAX(v0[1],MAX(v1[1],v2[1]))/pos->km_per_pixel	+ SMALLVAL+0.5);
			imin = (imin_dbl < INT_MIN) ? INT_MIN : (int) imin_dbl;
			imax = (imax_dbl > INT_MAX) ? INT_MAX : (int) imax_dbl;
			jmin = (jmin_dbl < INT_MIN) ? INT_MIN : (int) jmin_dbl;
			jmax = (jmax_dbl > INT_MAX) ? INT_MAX : (int) jmax_dbl;

			/*  Set the outbnd flag if the facet extends beyond the POS window*/

			if ((imin < -pos->n) || (imax > pos->n) ||	(jmin < -pos->n) || (jmax > pos->n))
				outbnd[0] = 1;

			/* Figure out if the facet projects at least partly within the POS
			 * 	window; if it does, look at each "contained" POS pixel and get
			 * 	the z-coordinate and cos(scattering angle)*/

			i1 = MAX( imin, -pos->n);
			i2 = MIN( imax,  pos->n);
			j1 = MAX( jmin, -pos->n);
			j2 = MIN( jmax,  pos->n);

			/*	Check if this facet is entirely outside the POS frame.	If that's
			 * 	the case, just keep track of the changed POS region. Original code
			 * 	called POSrect.  Here, we execute the code directly.*/

			if (i1 > pos->n || i2 < -pos->n || j1 > pos->n || j2 < -pos->n) {

				/*	Update the POS region that contains the target without regard
				 * 	to whether or not it extends beyond the POS frame*/

				//dev_POSrect(pos, src, imin_dbl, imax_dbl, jmin_dbl, jmax_dbl, max_overall);

			} else {
				/*	This path means that the facet is at least partially within the
				 * 	POS frame.  Must now find all POS pixels whose centers project
				 * 	onto this facet under consideration right now by this kernel.
				 *
				 * 	This means calling a second kernel that will open up child-grids
				 * 	of (i2-i1)*(j2-j1) threads for every facet.  This is dynamic
				 * 	parallelism in action.

			 	First figure out the child-kernel launch parameters
				 */
				int npixels = (i2-i1) * (j2-j1);
				int Bx = (maxThreadsPerBlock - 1 + npixels)/maxThreadsPerBlock;
				dim3 BLK, THD;
				BLK.x = Bx; BLK.y = 1; BLK.z = 1;
				THD.x = maxThreadsPerBlock; THD.y = 1; THD.z = 1;

				/*	Synchronize all threads in this block before continuing*/
				__syncthreads();

				//kernel launch here
				posvis_pixels <<<BLK, THD >>>(verts, pos,
						npixels, body, comp,
						v0[0], v0[1], v0[2],		/* components of v0 */
						v1[0], v1[1], v1[2],		/* components of v1 */
						v2[0], v2[1], v2[2],		/* components of v2 */
						n[0],  n[1],  n[2],			/* components of n  */
						max_overall,				/* used for pixel windowing */
						oa, cosa, cosb, usrc,
						smooth, src, bod, cmp, fac, zz, i1, j1, offset);

				cudaError_t code = cudaGetLastError();
				if (code != cudaSuccess)
					printf("GPUassert: %s %s %d\n", cudaGetErrorString(code), __FILE__, __LINE__);
				cudaDeviceSynchronize();

			}
		}
	}
	else {
		__syncthreads();
		cudaDeviceSynchronize();
	}
}

__global__ void posvis_facets_ndp(struct vertices_t *verts, struct pos_t *pos,
		int nf, double *oa, double *orbit_offset, double *max_overall,
		int smooth, int body, int comp, int **bod, int **cmp, int **fac, int *outbnd, int src,
		int maxThreadsPerBlock, float **fcosa, float **fcosb, float **fzz, float fusrc[3])
{
	/*	Find the facet index with thread and block indices.	Offset is calculated for a one-
	 * 	dimensional thread string. In other words, this kernel should be launched with (Bx,0,0)
	 * 	and (Tx,0,0)*/
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int imin, imax, jmin, jmax, i, i1, i2, j, j1, j2, k, flag = 0;
	double imin_dbl, imax_dbl, jmin_dbl, jmax_dbl, v0[3], v1[3], v2[3], n[3], x[3];
	double s, t, den;
	float z;
	int fac_ind0 = verts->f[offset].v[0];
	int fac_ind1 = verts->f[offset].v[1];
	int fac_ind2 = verts->f[offset].v[2];
	int xylim[4], xylim2[4];		// this is the temporary, per-facet xlim and ylim.

	if (offset < nf) { // make sure we stay within bounds of # of facets

//		/* 	transfer xlim and ylim from pos	*/
//		xylim[0] = pos->xlim[0];
//		xylim[1] = pos->xlim[1];
//		xylim[2] = pos->ylim[0];
//		xylim[3] = pos->ylim[1];
//		//if (src) {
//		xylim2[0] = pos->xlim2[0];
//		xylim2[1] = pos->xlim2[1];
//		xylim2[2] = pos->ylim2[0];
//		xylim2[3] = pos->ylim2[1];

		/* Get normal for facet (in body-fixed asteroid coord.), transform to
		 * observer coord.*/
		for (i=0; i<=2; i++)
			n[i] = verts->f[offset].n[i];
		dev_cotrans1( n, oa, n, 1);

		/* Consider facet further only if z-component of normal is >0. This
		 * ensures it is sort-of pointed at the observer.*/
		if (n[2] > 0.0)
		{
			/*	1st convert the 3 vertex coords for facet from body to obs co-
			 * 	ordinates. Then add orbit offset (center-of-mass offset) in obs
			 * 	coord for this model at this frame's epoch due to orbital mo-
			 * 	tion, in case the model is half of a binary system.*/
			dev_cotrans1(v0, oa, verts->v[fac_ind0].x, 1);
			dev_cotrans1(v1, oa, verts->v[fac_ind1].x, 1);
			dev_cotrans1(v2, oa, verts->v[fac_ind2].x, 1);

			for (i=0; i<=2; i++) {
				v0[i] += orbit_offset[i];
				v1[i] += orbit_offset[i];
				v2[i] += orbit_offset[i];
			}

			/* Find rectangular region (in POS pixels) containing the projected
			 * facet. Uses doubles in case model has illegal parameters and pi-
			 * xel numbers exceed the limits for valid integers.*/
			imin_dbl = floor ( MIN( v0[0], MIN( v1[0], v2[0] )) /
						pos->km_per_pixel - SMALLVAL+0.5);
			imax_dbl = floor ( MAX( v0[0], MAX( v1[0], v2[0] )) /
						pos->km_per_pixel + SMALLVAL+0.5);
			jmin_dbl = floor ( MIN( v0[1], MIN( v1[1], v2[1] )) /
						pos->km_per_pixel - SMALLVAL+0.5);
			jmax_dbl = floor ( MAX( v0[1], MAX( v1[1], v2[1] )) /
						pos->km_per_pixel + SMALLVAL+0.5);

			imin = (imin_dbl < INT_MIN) ? INT_MIN : (int) imin_dbl;
			imax = (imax_dbl > INT_MAX) ? INT_MAX : (int) imax_dbl;
			jmin = (jmin_dbl < INT_MIN) ? INT_MIN : (int) jmin_dbl;
			jmax = (jmax_dbl > INT_MAX) ? INT_MAX : (int) jmax_dbl;

			/*  Set outbnd flag if facet extends beyond POS window	*/
			if ((imin < -pos->n) || (imax > pos->n) ||
					(jmin < -pos->n) || (jmax > pos->n))
				atomicExch(&outbnd[0], 1);	// To-Do: May be able to remove atomics here

			/* Figure out if facet projects at least partly within POS window;
			 * if it does, look at each "contained" POS pixel and get z-coor-
			 * dinate and cos(scattering angle)*/
			i1 = MAX( imin, -pos->n);
			i2 = MIN( imax,  pos->n);
			j1 = MAX( jmin, -pos->n);
			j2 = MIN( jmax,  pos->n);

			/*	Check if facet is entirely outside POS frame. If so, just keep track of changed
			 * 	POS region. Original code called POSrect.  Changed to dev_POSrect.	*/
			if (i1 > pos->n || i2 < -pos->n || j1 > pos->n || j2 < -pos->n) {

				/*	Update POS region that contains the target without regard to whether or not
				 * 	it extends beyond the POS frame	*/
				dev_POSrect(pos, xylim, xylim2, src, imin_dbl, imax_dbl, jmin_dbl, jmax_dbl, max_overall, pos->n);

			} else {
				/*	This path means that the facet is at least partially within the POS frame.
				 * 	Must now find all POS pixels whose centers project onto this facet under
				 * 	consideration right now by this kernel. This is the non-dynamic processing
				 * 	kernel (NDP), so we do 2 nested for-loops here to get through the pixels. */

				for (i=i1; i<=i2; i++) {
					x[0] = i*pos->km_per_pixel;
					for (j=j1; j<=j2; j++) {
						x[1] = j*pos->km_per_pixel;

						/*  Compute parameters s(x,y) and t(x,y) which define a facet's surface as
						 *
						 *                     z = z0 + s*(z1-z0) + t*(z2-z1)
						 *
						 *	where z0, z1, and z2 are the z-coordinates at the vertices. The
						 *	conditions 0 <= s <= 1 and 0 <= t <= s require the POS pixel center to
						 *	be "within" the (projected) perimeter of facet f.                   */

						s = ( (x[0]-v0[0])*(v2[1]-v1[1]) - (v2[0]-v1[0])*(x[1]-v0[1]) )
				            		  * (den = 1/( (v1[0]-v0[0])*(v2[1]-v1[1]) - (v2[0]-v1[0])*(v1[1]-v0[1]) ));
						if ( (s >= -SMALLVAL) && (s <= 1.0+SMALLVAL) ) {
							t = ( (v1[0]-v0[0])*(x[1]-v0[1]) - (x[0]-v0[0])*(v1[1]-v0[1]) ) * den;
							if ( (t >= -SMALLVAL) && (t <= s+SMALLVAL) ) {

								/*	Compute the z coordinate of the center of this pixel:
 			                      	its distance measured from the origin towards Earth.    */
								z = v0[2] + s*(v1[2] - v0[2]) + t*(v2[2] - v1[2]);

								/*  If fac[i][j] is >= 0, pixel [i][j] was already assigned values
								 * 	during a previous call to posvis for a different model component.
								 * 	If so, override only if the current component is blocking our view
								 * 	of (i.e., is closer to us than) the previous one.   */
								float old = atomicMaxf(&fzz[i][j], z);
								if (old == z) flag = 1;

								//if ( (z > fzz[i][j]) || (fac[i][j] < 0) ) {
								if (flag || fac[i][j] < 0) {
									/*  Next line assigns distance of POS pixel center from COM towards
									 * 	Earth; that is, by changing zz, it changes pos->z or pos->zill */
									//zz[i][j] = z;
									if (!flag) atomicExch(&fzz[i][j], z);

									if (smooth) {

										/*  Get smoothed version of this facet's unit normal: Take the
										 * 	appropriate linear combination of the three vertex normals;
										 * 	transform from body to observer coordinates; and make sure
										 * 	that it points somewhat in our direction.                  */
										for (k=0; k<=2; k++)
											n[k] = verts->v[verts->f[offset].v[0]].n[k]
									        + s*(  verts->v[verts->f[offset].v[1]].n[k]
                                            - verts->v[verts->f[offset].v[0]].n[k])
                                            + t*(  verts->v[verts->f[offset].v[2]].n[k]
                                            - verts->v[verts->f[offset].v[1]].n[k]);

										dev_cotrans1( n, oa, n, 1);
										dev_normalize( n);
									}

									/*  Determine scattering and/or incidence angles. The following
									 * 	lines change cosa, thus changing pos->cose or pos->cosill.
									 * 	For bistatic situations (lightcurves) where we are viewing from
									 * 	Earth (src = 0), they also change cosb,thus changing pos->cosi. */
									if (n[2] > 0.0) {
										fcosa[i][j] = n[2];
										//atomicExch(&fcosa[i][j], n[2]);
										if ((!src) && (pos->bistatic)) {
											//fcosb[i][j] = dev_dot( n, usrc);
											float temp = (float)n[0]*fusrc[0] + (float)n[1]*fusrc[1] + (float)n[2]*fusrc[2];
											fcosb[i][j] = temp;
											//atomicExch(&fcosb[i][j], temp);
											if (fcosb[i][j] <= 0.0)
												fcosa[i][j] = 0.0;
										}
									}

									/*  Keep track of the changed POS region  */
									dev_POSrect( pos, xylim, xylim2, src, (double) i, (double) i, (double) j, (double) j,
											max_overall, pos->n);

									/*  Next lines change pos->body or pos->bodyill,pos->comp or
									 * 	pos->compill, and pos->f or pos->fill: the body, component,
									 * 	and facet numbers at the center of this POS pixel           */
									/*atomicExch(&bod[i][j], body);	//*/bod[i][j] = body;
									/*atomicExch(&cmp[i][j], comp);	//*/cmp[i][j] = comp;
									/*atomicExch(&fac[i][j], offset);	//*/fac[i][j] = offset;

								}  /* end if (no other facet yet blocks this facet from view) */
							}  /* end if 0 <= t <= s (facet center is "in" this POS pixel) */
						}  /* end if 0 <= s <= 1 */
					}  /* end j-loop over POS rows */
				}  /* end i-loop over POS columns */
			}
		}
		/* 	transfer xlim and ylim from pos	*/
//		pos->xlim[0] = xylim[0];
//		pos->xlim[1] = xylim[1];
//		pos->ylim[0] = xylim[2];
//		pos->ylim[1] = xylim[3];
//		pos->xlim2[0] = xylim2[0];
//		pos->xlim2[1] = xylim2[1];
//		pos->ylim2[0] = xylim2[2];
//		pos->ylim2[1] = xylim2[3];

		__syncthreads();
		//cudaDeviceSynchronize();
	}
	else {
		__syncthreads();
		//cudaDeviceSynchronize();
	}
}

extern "C"{
int posvis_cuda( struct vertices_t *verts, double orbit_offset[3], struct pos_t *pos,
		int smooth, int src, int body, int comp)
{
	int i, j, *outbnd, **bod, **cmp, **fac;
	double xfactor, yfactor, usrc[3], **cosa, **cosb, **zz;
	float fusrc[3];
	double *dev_orbit_offset, *oa, *max_overall;
	float **fzz, **fcosa, **fcosb;

	/*	Allocate unified memory	*/
	//npixels = (pos->xlim[1]-pos->xlim[0])*(pos->ylim[1]-pos->ylim[0]);
	cudaCalloc((void**)&outbnd, 			sizeof(int), 		      1);
	cudaCalloc((void**)&max_overall,		sizeof(double), 	      4);
	cudaCalloc((void**)&oa, 				sizeof(double), 	      9);
	cudaCalloc((void**)&dev_orbit_offset,	sizeof(double),    	      3);
	cudaCalloc((void**)&bod, 				sizeof(int*), 	 2*pos->n+1);
	cudaCalloc((void**)&cmp, 				sizeof(int*),  	 2*pos->n+1);
	cudaCalloc((void**)&fac, 				sizeof(int*), 	 2*pos->n+1);
	cudaCalloc((void**)&cosa, 				sizeof(double*), 2*pos->n+1);
	cudaCalloc((void**)&cosb, 				sizeof(double*), 2*pos->n+1);
	cudaCalloc((void**)&zz, 				sizeof(double*), 2*pos->n+1);
	cudaCalloc((void**)&fzz,				sizeof(float*),	 2*pos->n+1);
	cudaCalloc((void**)&fcosa,				sizeof(float*),	 2*pos->n+1);
	cudaCalloc((void**)&fcosb,				sizeof(float*),  2*pos->n+1);
	cudaCalloc((void**)&usrc,				sizeof(double),		      3);
	cudaCalloc((void**)&fusrc,				sizeof(float),		      3);

	/*	Fix pointer offset for float arrays	*/
	fzz -= -pos->n;
	fcosa -= -pos->n;
	fcosb -= -pos->n;

	/* Inner loop of allocating the single pointers along the outer loop */
	for (i=-pos->n; i<=pos->n; i++) {
		/* First allocate with cuda managed memory */
		cudaCalloc((void**)&fzz[i],		sizeof(float), (2*pos->n+1));
		cudaCalloc((void**)&fcosa[i],	sizeof(float), (2*pos->n+1));
		cudaCalloc((void**)&fcosb[i],	sizeof(float), (2*pos->n+1));
		/* Now offset indexing for the inner loop */
		fzz[i] 		-= -pos->n;
		fcosa[i] 	-= -pos->n;
		fcosb[i] 	-= -pos->n;
	}

	//  Initialize variables
	outbnd[0] = 0;
	pos->posbnd_logfactor = 0.0;
	max_overall[0] =  HUGENUMBER;	// imin_overall
	max_overall[1] = -HUGENUMBER;	// imax_overall
	max_overall[2] =  HUGENUMBER;	// jmin_overall
	max_overall[3] = -HUGENUMBER;	// jmax_overall
	mtrnsps_cuda(oa, pos->ae);		// transpose pos->ae into oa

	if (src) {
		/* We're viewing the model from the sun: at the center of each pixel in
		 * the projected view, we want cos(incidence angle), distance from the
		 * COM towards the sun, and the facet number.*/
		mmmul_cuda( oa, pos->se, oa);	//	oa takes ast into sun coords
		cosa = pos->cosill;            	//	cos(incidence angle)
		cosb = NULL;                    //	<avoid compilation warnings>
		zz = pos->zill;                	//	distance towards sun
		bod = pos->bodyill;            	//	body number at projected pixel center
		cmp = pos->compill;            	//	component number at projected pixel center
		fac = pos->fill;               	//	facet number at projected pixel center
		fcosb = NULL;
		/* Copy zill and cosill into their float equivalents	*/
		for (j=-pos->n; j<=pos->n; j++) {
			for (i=-pos->n; i<=pos->n; i++){
				fzz[i][j] = (float)pos->zill[i][j];
				fcosa[i][j] = (float)pos->cosill[i][j];
			}
		}
	} else {
		/* We're viewing model from Earth: at the center of each POS pixel we
		 * want cos(scattering angle), distance from COM towards Earth, and fa-
		 * cet number. If bistatic (lightcurves) we also want cos(incidence an-
		 * gle) and unit vector towards source.*/
		mmmul_cuda( oa, pos->oe, oa);   //	oa takes ast into obs coords
		cosa = pos->cose;              	//	scattering angle
		cosb = pos->cosi;              	//	incident angle
		zz = pos->z;                   	// 	observer z-coordinate (backwards along LOS)
		bod = pos->body;               	//	body number at projected pixel center
		cmp = pos->comp;               	//	component number at POS pixel center
		fac = pos->f;                  	//	facet number at POS pixel center
		if (pos->bistatic) {
			usrc[0] = usrc[1] = 0.0;   	//	unit vector towards source
			usrc[2] = 1.0;
			cotrans_cuda( usrc, pos->se, usrc, -1);
			cotrans_cuda( usrc, pos->oe, usrc,  1);    //	in observer coordinates
			for (i=0; i<3; i++)
				fusrc[i] = (float)usrc[i];
		}

		/* Copy zz, cosa, and cosb into their float equivalents	*/
		for (j=-pos->n; j<=pos->n; j++) {
			for (i=-pos->n; i<=pos->n; i++){
				fzz[i][j] = 0;//(float)pos->z[i][j];
				fcosa[i][j] = (float)pos->cose[i][j];
				fcosb[i][j] = (float)pos->cosi[i][j];;
			}
		}
	}

	/* For all comments below this point: If the "src" argument is true, "observer" is
	 * actually "source" (sun), "Earth" is actually "source" and "POS" is actually
	 * "projection as viewed from the source"		*/
	/* This kernel performs the same tasks as the facets loop in the original posvis */
	/*	Now calculate launch parameters for the facets kernel and launch it. 		 */
	maxThreadsPerBlock = 512;
	int Bx = floor((maxThreadsPerBlock - 1 + verts->nf) / maxThreadsPerBlock);
	dim3 BLK, THD;
	BLK.x = Bx; BLK.y = 1; BLK.y = 1;					// Grid-block dimension for the 1-D case
	THD.x = maxThreadsPerBlock; THD.y = 1; THD.z = 1;	// Thread block dimensions

	//kernel launch here
	if (!DYNPROC)
		posvis_facets_ndp <<<BLK, THD >>>(verts, pos, verts->nf, oa, dev_orbit_offset,
				max_overall, smooth, body, comp, bod, cmp, fac, outbnd, src,
				maxThreadsPerBlock, fcosa, fcosb, fzz, fusrc);
	else
		posvis_facets <<<BLK, THD >>>(verts, pos, verts->nf, oa, dev_orbit_offset,
				max_overall, smooth, body, comp, bod, cmp, fac, outbnd, src,
				maxThreadsPerBlock, cosa, cosb, zz, usrc);
	char *location = "posvis_cuda.cu line 692";
	checkErrorAfterKernelLaunch(location);
	deviceSyncAfterKernelLaunch(location);

	/*	Copy fzz, fcosa, fcosb back to their double equivalents	*/
	for (j=-pos->n; j<=pos->n; j++) {
		for (i=-pos->n; i<=pos->n; i++){
			pos->z[i][j] = (double)fzz[i][j];
			pos->cose[i][j] = (double)fcosa[i][j];
			pos->cosi[i][j] = (double)fcosb[i][j];
		}
	}
	/* If the target extends beyond the POS frame, set pos->posbnd_logfactor equal
     to the logarithm of the ratio of the number of pixels in a frame extended to
     include the entire target, divided by the number of pixels in the actual frame	*/

	if (outbnd[0]) {
		xfactor = (MAX(max_overall[1],pos->n) - MIN(max_overall[0],-pos->n)+1) / (2*pos->n+1);
		yfactor = (MAX(max_overall[3],pos->n) - MIN(max_overall[2],-pos->n)+1) / (2*pos->n+1);
		pos->posbnd_logfactor = log(xfactor*yfactor);
		printf("OUTBOUND!\n\n");
	}

	/*------------------------------------------------------------------------------------------*/
	/* START DEBUGGING CODE SECTION - WRITES RESULTS TO CSV	*/
	/*------------------------------------------------------------------------------------------*/

//	npixels = (pos->xlim[1]-pos->xlim[0])*(pos->ylim[1]-pos->ylim[0]);
//	float sumz = 0.0;
//	float sumcsa = 0.0;
//	float sumcsb = 0.0;
//	int sumbdy = 0;
//	int sumcmp = 0;
//	int sumfct = 0;
//
//	for (i = -pos->n; i <= pos->n; i++)
//		for (j = -pos->n; j <= pos->n; j++){
//			if (fzz[i][j] >= 0) sumz += fzz[i][j];//pos->z[i][j];}
//			sumcsa += fcosa[i][j];//pos->cose[i][j];
//			sumcsb += fcosb[i][j];//pos->cosi[i][j];
//			sumbdy += pos->body[i][j];
//			sumcmp += pos->comp[i][j];
//			sumfct += pos->f[i][j];
//		}
//	printf("pos->z: %f\n", sumz);
//	printf("pos->cosa: %f\n", sumcsa);
//	//printf("pos->cosb: %f\n", sumcsb);
//	//printf("pos->body: %i\n", sumbdy);
//	printf("pos->comp: %i\n", sumcmp);
//	printf("pos->f: %i\n", sumfct);
//	int cntr = 0;
//	for (i=0; i<2; i++){
//		printf("xlim[%i] = %i\n", i, pos->xlim[i]);
//		printf("ylim[%i] = %i\n", i, pos->ylim[i]);
//	}
//	for (i = pos->xlim[0]; i <= pos->xlim[1]; i++)
//		for (j = pos->ylim[0]; j <= pos->ylim[1]; j++)
//			if (pos->z[i][j] != HUGENUMBER && pos->z[i][j] != -HUGENUMBER){
//				//printf("pos->z[%i][%i] = %g\n", i, j, pos->z[i][j]);
//				cntr++;
//			}
//	printf("Total pixels > 0 pwr: %i\n", cntr);
//	printf("npixels in posvis_cuda: %i\n", npixels);

	//	FILE *fp_z;
	//	FILE *fp_f;
	//	FILE *fp_cose;
	//	FILE *fp_cosi;
	//	char *filename_z, *filename_f, *filename_cose, *filename_cosi;
	//	filename_z = 	"dbg_z_cuda.csv";
	//	filename_f = 	"dbg_f_cuda.csv";
	//	filename_cose = "dbg_cose_cuda.csv";
	//	filename_cosi = "dbg_cosi_cuda.csv";
	//
	//	//printf("\n %sfile created",filename_z);
	//	//printf("\n\nFilename: %s",filename_z);
	//	fp_z = 		fopen(filename_z,		"w+");
	//	fp_f = 		fopen(filename_f, 		"w+");
	//	fp_cose = 	fopen(filename_cose, 	"w+");
	//	fp_cosi = 	fopen(filename_cosi, 	"w+");
	//
	//	fprintf(fp_z, 		"i/j , ");
	//	fprintf(fp_f, 		"i/j , ");
	//	fprintf(fp_cose, 	"i/j , ");
	//	fprintf(fp_cosi, 	"i/j , ");
	//
	//	// The following loop writes the top row of -i to + i
	//	for (i = -pos->n; i <= pos->n; i++){
	//		fprintf(fp_z, 		"%i , ", i);
	//		fprintf(fp_f, 		"%i , ", i);
	//		fprintf(fp_cose, 	"%i , ", i);
	//		fprintf(fp_cosi, 	"%i , ", i);
	//	}
	//
	//	// The following loops write the pos->z data content into the file
	//	for (j = -pos->n; j <= pos->n; j++){
	//		fprintf(fp_z, 		"\n%i , ", j);
	//		fprintf(fp_f, 		"\n%i , ", j);
	//		fprintf(fp_cose, 	"\n%i , ", j);
	//		fprintf(fp_cosi, 	"\n%i , ", j);
	//
	//		for (i = -pos->n; i <= pos->n; i++){
	//			fprintf(fp_z, 		"%f , ", pos->z[i][j]);
	//			fprintf(fp_f, 		"%i , ", pos->f[i][j]);
	//			fprintf(fp_cose, 	"%f , ", pos->cose[i][j]);
	//			fprintf(fp_cosi, 	"%f , ", pos->cosi[i][j]);
	//		}
	//	}
	//
	//	fclose(fp_z);
	//	fclose(fp_f);
	//	fclose(fp_cose);
	//	fclose(fp_cosi);
	/*------------------------------------------------------------------------------------------*/
	/*	END DEBUGGING SECTION	*/
	/*------------------------------------------------------------------------------------------*/
	// Release any cudaMallocManaged-allocated objects
	cudaFree(max_overall);
	cudaFree(oa);
	cudaFree(dev_orbit_offset);
//		cudaFree(bod);
//		cudaFree(cmp);
//		cudaFree(fac);
//		cudaFree(cosa);
//		cudaFree(cosb);
//		cudaFree(zz);

	return outbnd[0];
}
}
