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
 for the no-smoothing case
 *****************************************************************************************/
extern "C" {
#include "head.h"
#include <limits.h>
}
__device__ double oa[3][3], usrc[3], orbit_offset[3];
__device__ float *cosa, *cosb, *zz;
__device__ float imin_overall, imax_overall, jmin_overall, jmax_overall;
__device__ int **bod, **cmp, **fac, posvis_nf, doutbnd, smooth, psvs_n;
__device__ struct pos_t *pos;
__device__ struct vertices_t *verts;

__device__ int dbg_pxl_occ=0;

//__device__ void dev_cotrans3(double y[3], double a[3][3], double x[3],
//		int dir) {
//	double t[3];
//	int i, j;
//
//	if (dir == 1)
//		for (i = 0; i <= 2; i++) {
//			t[i] = 0.0;
//			for (j = 0; j <= 2; j++)
//				t[i] += a[i][j] * x[j];
//		}
//	if (dir == (-1))
//		for (i = 0; i <= 2; i++) {
//			t[i] = 0.0;
//			for (j = 0; j <= 2; j++)
//				t[i] += a[j][i] * x[j];
//		}
//	for (i = 0; i <= 2; i++)
//		y[i] = t[i];
//}
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
/* Update the rectangular plane-of-sky region that "contains" the target */
/* the i/jmin/max_overall variables have been removed here as arguments
 * because they are declared as __device__ variables at file scope instead  */
__device__ void dev_POSrect2(struct pos_t *pos, int src, float imin_dbl,
		float imax_dbl, float jmin_dbl, float jmax_dbl)	{
	int n, imin, imax, jmin, jmax;
	n = psvs_n;

	/*  Update the POS region that contains the target without
	 regard to whether or not it extends beyond the POS frame  */

	atomicMinf(&imin_overall, imin_dbl);
	atomicMaxf(&imax_overall, imax_dbl);
	atomicMinf(&jmin_overall, jmin_dbl);
	atomicMaxf(&jmax_overall, jmax_dbl);

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

	/* The below code is for debugging use only. It drops the use of atomics.*/
//	if (src) {
//		pos->xlim2[0] = min(pos->xlim2[0], imin);
//		pos->xlim2[1] = max(pos->xlim2[1], imax);
//		pos->ylim2[0] = min(pos->ylim2[0], jmin);
//		pos->ylim2[1] = max(pos->ylim2[1], jmax);
//	} else {
//		pos->xlim[0] = min(pos->xlim[0], imin);
//		pos->xlim[1] = max(pos->xlim[1], imax);
//		pos->ylim[0] = min(pos->ylim[0], jmin);
//		pos->ylim[1] = max(pos->ylim[1], jmax);
//	}

}
__global__ void posvis_init_krnl(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, int set, int frame, int src,
		double orbit_offset0, double orbit_offset1, double orbit_offset2) {
	/* Single-threaded kernel to init variables */
	if (threadIdx.x == 0) {
		switch (ddat->set[set].type) {
		case DELAY:
			pos = &ddat->set[set].desc.deldop.frame[frame].pos;
			break;
		case DOPPLER:
			pos = &ddat->set[set].desc.doppler.frame[frame].pos;
			break;
		case POS:
			pos = &ddat->set[set].desc.poset.frame[frame].pos;
			break;
		case LGHTCRV:
			pos = &ddat->set[set].desc.lghtcrv.rend[frame].pos; //frame = i
			break;
		}
		/*  Initialize variables  */
		doutbnd = 0;
		pos->posbnd_logfactor = 0.0;
		imin_overall = jmin_overall = HUGENUMBER;
		imax_overall = jmax_overall = -HUGENUMBER;
		dev_mtrnsps(oa, pos->ae);

		if (src) {

			/* We're viewing the model from the sun: at the center of each pixel
			 * in the projected view, we want cos(incidence angle), distance from
			 * the COM towards the sun, and the facet number.                */
			dev_mmmul(oa, pos->se, oa); /* oa takes ast into sun coords           */
			cosa = pos->cosill_s; /* cos(incidence angle)                   */
			cosb = NULL; /* <avoid compilation warnings>           */
			zz = pos->zill_s; /* distance towards sun                   */
			bod = pos->bodyill; /* body number at projected pixel center  */
			cmp = pos->compill; /* component number at projected pixel center */
			fac = pos->fill; /* facet number at projected pixel center */

		} else {

			/* We're viewing the model from Earth: at the center of each POS pixel
			 * we want cos(scattering angle), distance from the COM towards Earth,
			 * and the facet number.  For bistatic situations (lightcurves) we also
			 want cos(incidence angle) and the unit vector towards the source.     */
			dev_mmmul(oa, pos->oe, oa); /* oa takes ast into obs coords */
			cosa = pos->cose_s; /* scattering angle */
			cosb = pos->cosi_s; /* incident angle */
			zz = pos->z_s; /* observer z-coordinate (backwards along LOS) */
			bod = pos->body; /* body number at projected pixel center  */
			cmp = pos->comp; /* component number at POS pixel center */
			fac = pos->f; /* facet number at POS pixel center */
			if (pos->bistatic) {
				usrc[0] = usrc[1] = 0.0; /* unit vector towards source */
				usrc[2] = 1.0;
				dev_cotrans3(usrc, pos->se, usrc, -1);
				dev_cotrans3(usrc, pos->oe, usrc, 1); /* in observer coordinates */
			}
		}
		/* transfer the host orbit_offset[3] values to the device copy */
		orbit_offset[0] = orbit_offset0;
		orbit_offset[1] = orbit_offset1;
		orbit_offset[2] = orbit_offset2;

		verts = &dmod->shape.comp[0].real;
		posvis_nf = verts->nf;
		smooth = dpar->pos_smooth;
		psvs_n = pos->n;

		dbg_pxl_occ = 0;
	}
}
__global__ void posvis_fct_pxl_krnl(int i1, int i2, int j1, int j2,
		double v0x, double v0y, double v0z, double v1x, double v1y, double v1z,
		double v2x, double v2y, double v2z, double nx, double ny, double nz,
		int src, int body, int comp, int f, int fvx, int fvy, int fvz) {

	/* Multi-threaded kernel to be launched from facet kernel */
	int ispan = i2-i1+1;		int jspan = j2-j1+1;
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int npixels = ispan*jspan;

	int i = offset % ispan + i1;
	int j = offset / ispan + j1;

	unsigned int pxa; // Pixel address in 1D float arrays (zzf, cosaf, etc)
	unsigned int k, f_v[3];
	double den, s, t, z, v0[3], v1[3], v2[3], n[3], x[3];

	v0[0]=v0x; 	v0[1]=v0y; 	v0[2]=v0z;
	v1[0]=v1x; 	v1[1]=v1y; 	v1[2]=v1z;
	v2[0]=v2x; 	v2[1]=v2y; 	v2[2]=v2z;
	n[0]=nx;   	n[1]=ny;   	n[2]=nz;
	f_v[0]=fvx; f_v[1]=fvy; f_v[2]=fvz;

	if (offset < npixels){
		/* Calculate (x,y) position of pixel and linear pixel address */
		x[0] = i * pos->km_per_pixel;
		x[1] = j * pos->km_per_pixel;
		pxa = (j+pos->n)*(2*pos->n+1)+(i+pos->n);

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
				z = v0[2] + s * (v1[2] - v0[2])  + t * (v2[2] - v1[2]);

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

				float old = atomicMaxf(&zz[pxa], z);
				//	if ((z > zz[i][j]) || (fac[i][j] < 0)) {
				if (old < z || fac[i][j] < 0) {

					/* Next line assigns distance of POS pixel
					 * center from COM towards Earth; that is,
					 * by changing zz,it changes pos->z or
					 * pos->zill                */
					/* following line is a first time z calc
					 * for this pixel  */
					if (fac[i][j] < 0) atomicExch(&zz[pxa], z);
					//									zz[i][j] = z;
					if (smooth) {

						/* Get smoothed version of facet unit
						 * normal: Take the linear combination
						 * of the three vertex normals; trans-
						 * form from body to observer coordina-
						 * tes; and make sure that it points
						 * somewhat in our direction.         */
						for (k = 0; k <= 2; k++)
							n[k] =	verts->v[f_v[0]].n[k]
							 + s * (verts->v[f_v[1]].n[k] - verts->v[f_v[0]].n[k])
							 + t * (verts->v[f_v[2]].n[k] - verts->v[f_v[1]].n[k]);
						dev_cotrans3(n, oa, n, 1);
						dev_normalize(n);
					}

					/* Determine scattering and/or incidence
					 * angles. Next lines change pos->cose/
					 * cosill. If bistatic (lightcurves), where
					 * we are viewing from Earth (src = 0),
					 * pos->cosi is also changed.                 */
					if (n[2] > 0.0) {
						atomicExch(&cosa[pxa], n[2]);
						if ((!src) && (pos->bistatic)) {
							float temp = (float)dev_dot(n,usrc);
							atomicExch(&cosb[pxa], temp);
							if (cosb[pxa] <= 0.0)
								cosa[pxa] = 0.0;
						}
					}

					/*  Keep track of the changed POS region  */
					dev_POSrect2(pos, src, (double) i,
							(double) i, (double) j, (double) j);

					/* Next lines change pos->body/bodyill,
					 * pos->comp/compill, pos->f/fill          */
					bod[i][j] = body;
					cmp[i][j] = comp;
					fac[i][j] = f;

				} /* end if (no other facet yet blocks this facet from view) */
			} /* end if 0 <= t <= s (facet center is "in" this POS pixel) */
		} /* end if 0 <= s <= 1 */
	}
}
__global__ void posvis_fct_dynp_krnl(int src, int body, int comp) {
	/* nf-threaded kernel */
	/* Variables used:
	 * 	f = offset and facet index
	 *	i,j,k = iterators
	 *	i1,i2,j1,j2 = used for pos bracketing
	 *	imin,imax,jmin,jmax = used for pos bracketing
	 *	imin_dbl,imax_dbl,jmin_dbl,jmax_dbl = same as the last two
	 *	xfactor,yfactor = for calculating pos->posbnd_logfactor
	 *	n[3] = normal for this facet/thread
	 *	v0[3],v1[3],v2[3] = the three vertices making up this facet
	 *	s,t = s(x,y) and t(x,y) used in surface calculation
	 *	z = distance to pixel center
	 *	den = denominator used in the surface calculation
	 */
	int f = blockIdx.x * blockDim.x + threadIdx.x;
	int i, i1, i2, j1, j2, imin, imax, jmin, jmax;
	double n[3], v0[3], v1[3], v2[3];
	float imin_dbl, imax_dbl, jmin_dbl, jmax_dbl;
	int f_v[3];//0, f_v1, f_v2;

	if (f < posvis_nf) {
		/* 1st, put the indices for this facet's 3 vertices into f_v[]. Then
		 * get this facet's normal in body-fixed (asteroid) coordinates and
		 * convert it to observer coordinates     */
		for (i=0; i<3; i++) {
			f_v[i] = verts->f[f].v[i];
			n[i] = verts->f[f].n[i];
		}
		dev_cotrans3(n, oa, n, 1);

		/* Consider facet further only if normal points towards observer  */
		if (n[2] > 0.0) {
			/* Convert 3 sets of vertex coordinates from body to observer
			 * coordinates; orbit_offset is the center-of-mass offset (in obs.
			 * co.) for model at frame epoch due to orbital motion, in case
			 * model is half of a binary system.  */
			dev_cotrans3(v0, oa, verts->v[f_v[0]].x, 1);
			dev_cotrans3(v1, oa, verts->v[f_v[1]].x, 1);
			dev_cotrans3(v2, oa, verts->v[f_v[2]].x, 1);
			for (i = 0; i <= 2; i++) {
				v0[i] += orbit_offset[i];
				v1[i] += orbit_offset[i];
				v2[i] += orbit_offset[i];
			}

			/* Find rectangular region (in POS pixels) containing the projected
			 * facet - use floats in case model has illegal parameters and the
			 * pixel numbers exceed the limits for valid integers                         */
			imin_dbl = floor(MIN(v0[0],MIN(v1[0],v2[0])) / pos->km_per_pixel
							- SMALLVAL + 0.5);
			imax_dbl = floor(MAX(v0[0],MAX(v1[0],v2[0])) / pos->km_per_pixel
							+ SMALLVAL + 0.5);
			jmin_dbl = floor(MIN(v0[1],MIN(v1[1],v2[1])) / pos->km_per_pixel
							- SMALLVAL + 0.5);
			jmax_dbl = floor(MAX(v0[1],MAX(v1[1],v2[1])) / pos->km_per_pixel
							+ SMALLVAL + 0.5);
			imin = (imin_dbl < INT_MIN) ? INT_MIN : (int) imin_dbl;
			imax = (imax_dbl > INT_MAX) ? INT_MAX : (int) imax_dbl;
			jmin = (jmin_dbl < INT_MIN) ? INT_MIN : (int) jmin_dbl;
			jmax = (jmax_dbl > INT_MAX) ? INT_MAX : (int) jmax_dbl;

			/*  Set the outbnd flag if the facet extends beyond the POS window  */
			if ((imin < (-pos->n)) || (imax > pos->n) || (jmin < (-pos->n))
					|| (jmax > pos->n))
				doutbnd = 1;

			/* Figure out if facet projects at least partly within POS window;
			 * if it does, look at each "contained" POS pixel and get the
			 * z-coordinate and cos(scattering angle)           */
			i1 = MAX(imin, -pos->n);
			i2 = MIN(imax, pos->n);
			j1 = MAX(jmin, -pos->n);
			j2 = MIN(jmax, pos->n);

			if (i1 > pos->n || i2 < -pos->n || j1 > pos->n || j2 < -pos->n) {

				/* Facet is entirely outside the POS frame: just keep track of
				 * changed POS region     */
				dev_POSrect2(pos, src, imin_dbl, imax_dbl, jmin_dbl, jmax_dbl);
			} else {

				/* Facet is at least partly within POS frame: find all POS
				 * pixels whose centers project onto this facet  */
				/* Here we use dynamic parallelism to launch a child grid of
				 * threads for the pixel loops */
				dim3 BLK, THD;
				unsigned int dmaxThreadsPerBlock = 32;
				unsigned int npixels, ispan, jspan;
				ispan = i2-i1 + 1;			jspan = j2-j1 + 1;
				npixels = ispan * jspan;
				int Bx = (dmaxThreadsPerBlock - 1 + npixels)/dmaxThreadsPerBlock;
				BLK.x = Bx;
				THD.x = dmaxThreadsPerBlock;

				if (npixels > 0) {
				posvis_fct_pxl_krnl<<<BLK,THD>>>(i1,i2,j1,j2, v0[0],v0[1],v0[2],
						v1[0],v1[1],v1[2], v2[0],v2[1],v2[2], n[0],n[1],n[2],
						src,body,comp,f,f_v[0], f_v[1],f_v[2]);
				}
			} /* end else of if (i1 > pos->n || i2 < -pos->n || j1 > pos->n || j2 < -pos->n) */
		} /* End if (n[2] > 0.0) */
	} /* end if (f < nf) */
}
__global__ void posvis_facet_krnl(int src, int body, int comp) {
	/* nf-threaded kernel */
	/* Variables used:
	 * 	f = offset and facet index
	 *	i,j,k = iterators
	 *	i1,i2,j1,j2 = used for pos bracketing
	 *	imin,imax,jmin,jmax = used for pos bracketing
	 *	imin_dbl,imax_dbl,jmin_dbl,jmax_dbl = same as the last two
	 *	xfactor,yfactor = for calculating pos->posbnd_logfactor
	 *	n[3] = normal for this facet/thread
	 *	v0[3],v1[3],v2[3] = the three vertices making up this facet
	 *	s,t = s(x,y) and t(x,y) used in surface calculation
	 *	z = distance to pixel center
	 *	den = denominator used in the surface calculation
	 */
	int f = blockIdx.x * blockDim.x + threadIdx.x;
	int k, i, i1, i2, j, j1, j2, imin, imax, jmin, jmax;
	double n[3], v0[3], v1[3], v2[3], x[3], s, t, z, den;
	float imin_dbl, imax_dbl, jmin_dbl, jmax_dbl;
	int f_v0, f_v1, f_v2;
	__shared__ int posn, span;
	__shared__ float km_per_pixel;
	posn = pos->n;
	span = 2*posn+1;
	km_per_pixel = (float)pos->km_per_pixel;

	if (f < posvis_nf) {

		f_v0 = verts->f[f].v[0];
		f_v1 = verts->f[f].v[1];
		f_v2 = verts->f[f].v[2];

		/* Get the normal to this facet in body-fixed (asteroid) coordinates
		 * and convert it to observer coordinates     */
		for (i = 0; i <= 2; i++)
			n[i] = verts->f[f].n[i];

		dev_cotrans3(n, oa, n, 1);

		/* Consider this facet further only if its normal points somewhat
		 * towards the observer rather than away         */
		if (n[2] > 0.0) {
			/* Convert the three sets of vertex coordinates from body to ob-
			 * server coordinates; orbit_offset is the center-of-mass offset
			 * (in observer coordinates) for this model at this frame's epoch
			 * due to orbital motion, in case the model is half of a binary
			 * system.  */
			dev_cotrans3(v0, oa, verts->v[f_v0].x, 1);
			dev_cotrans3(v1, oa, verts->v[f_v1].x, 1);
			dev_cotrans3(v2, oa, verts->v[f_v2].x, 1);
			for (i = 0; i <= 2; i++) {
				v0[i] += orbit_offset[i];
				v1[i] += orbit_offset[i];
				v2[i] += orbit_offset[i];
			}

			/* Find rectangular region (in POS pixels) containing the projected
			 * facet - use floats in case model has illegal parameters and the
			 * pixel numbers exceed the limits for valid integers                         */
			imin_dbl = floor(MIN(v0[0],MIN(v1[0],v2[0])) / km_per_pixel
							- SMALLVAL + 0.5);
			imax_dbl = floor(MAX(v0[0],MAX(v1[0],v2[0])) / km_per_pixel
							+ SMALLVAL + 0.5);
			jmin_dbl = floor(MIN(v0[1],MIN(v1[1],v2[1])) / km_per_pixel
							- SMALLVAL + 0.5);
			jmax_dbl = floor(MAX(v0[1],MAX(v1[1],v2[1])) / km_per_pixel
							+ SMALLVAL + 0.5);
			imin = (imin_dbl < INT_MIN) ? INT_MIN : (int) imin_dbl;
			imax = (imax_dbl > INT_MAX) ? INT_MAX : (int) imax_dbl;
			jmin = (jmin_dbl < INT_MIN) ? INT_MIN : (int) jmin_dbl;
			jmax = (jmax_dbl > INT_MAX) ? INT_MAX : (int) jmax_dbl;

			/*  Set the outbnd flag if the facet extends beyond the POS window  */
			if ((imin < (-posn)) || (imax > posn) || (jmin < (-posn))
					|| (jmax > posn))
				doutbnd = 1;

			/* Figure out if facet projects at least partly within POS window;
			 * if it does, look at each "contained" POS pixel and get the
			 * z-coordinate and cos(scattering angle)           */
			i1 = MAX(imin, -posn);
			i2 = MIN(imax, posn);
			j1 = MAX(jmin, -posn);
			j2 = MIN(jmax, posn);

			if (i1 > posn || i2 < -posn || j1 > posn || j2 < -posn) {

				/* Facet is entirely outside the POS frame: just keep track of
				 * changed POS region     */
				dev_POSrect2(pos, src, imin_dbl, imax_dbl, jmin_dbl, jmax_dbl);
			} else {
				dev_POSrect2(pos, src, imin_dbl, imax_dbl, jmin_dbl, jmax_dbl);
				/* Facet is at least partly within POS frame: find all POS
				 * pixels whose centers project onto this facet  */
				int pxa;	// pixel address in 1D float array zzf etc
				for (i = i1; i <= i2; i++) {

					x[0] = i * km_per_pixel;
					for (j = j1; j <= j2; j++) {

						x[1] = j * km_per_pixel;

						/* Calculate the pixel address for 1D arrays */
						pxa = (j+posn)*span + (i+posn);

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
								z = v0[2] + s * (v1[2] - v0[2])
										+ t * (v2[2] - v1[2]);

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

								float old = atomicMaxf(&zz[pxa], z);

//								if ((z > zz[i][j]) || (fac[i][j] < 0)) {
								if (old < z || fac[i][j] < 0) {
atomicAdd(&dbg_pxl_occ, 1);
									/* Next line assigns distance of POS pixel
									 * center from COM towards Earth; that is,
									 * by changing zz,it changes pos->z or
									 * pos->zill                */
									/* following line is a first time z calc
									 * for this pixel  */
									if (fac[i][j] < 0) atomicExch(&zz[pxa], z);
//									zz[i][j] = z;
									if (smooth) {

										/* Get smoothed version of facet unit
										 * normal: Take the linear combination
										 * of the three vertex normals; trans-
										 * form from body to observer coordina-
										 * tes; and make sure that it points
										 * somewhat in our direction.         */
										for (k = 0; k <= 2; k++)
											n[k] =	verts->v[f_v0].n[k]
											 + s * (verts->v[f_v1].n[k] - verts->v[f_v0].n[k])
											 + t * (verts->v[f_v2].n[k]	- verts->v[f_v1].n[k]);
										dev_cotrans3(n, oa, n, 1);
										dev_normalize(n);
									}

									/* Determine scattering and/or incidence
									 * angles. Next lines change pos->cose/
									 * cosill. If bistatic (lightcurves), where
									 * we are viewing from Earth (src = 0),
									 * pos->cosi is also changed.                 */
									if (n[2] > 0.0) {
										atomicExch(&cosa[pxa], n[2]);
										if ((!src) && (pos->bistatic)) {
											float temp = (float)dev_dot(n,usrc);
											atomicExch(&cosb[pxa], temp);
											if (cosb[pxa] <= 0.0)
												cosa[pxa] = 0.0;
										}
									}

									/* Next lines change pos->body/bodyill,
									 * pos->comp/compill, pos->f/fill          */
									bod[i][j] = body;
									cmp[i][j] = comp;
									fac[i][j] = f;

								} /* end if (no other facet yet blocks this facet from view) */
							} /* end if 0 <= t <= s (facet center is "in" this POS pixel) */
						} /* end if 0 <= s <= 1 */
					} /* end j-loop over POS rows */
				} /* end i-loop over POS columns */
			} /* end else of if (i1 > pos->n || i2 < -pos->n || j1 > pos->n || j2 < -pos->n) */
		} /* End if (n[2] > 0.0) */
	} /* end if (f < nf) */
}
__global__ void posvis_facet_sep_krnl(int src, int body, int comp,
		float *iminflt, float *imaxflt, float *jminflt, float *jmaxflt) {
	/* nf-threaded kernel */
	/* Variables used:
	 * 	f = offset and facet index
	 *	i,j,k = iterators
	 *	i1,i2,j1,j2 = used for pos bracketing
	 *	imin,imax,jmin,jmax = used for pos bracketing
	 *	imin_dbl,imax_dbl,jmin_dbl,jmax_dbl = same as the last two
	 *	xfactor,yfactor = for calculating pos->posbnd_logfactor
	 *	n[3] = normal for this facet/thread
	 *	v0[3],v1[3],v2[3] = the three vertices making up this facet
	 *	s,t = s(x,y) and t(x,y) used in surface calculation
	 *	z = distance to pixel center
	 *	den = denominator used in the surface calculation
	 */
	int f = blockIdx.x * blockDim.x + threadIdx.x;
	int k, i, i1, i2, j, j1, j2, imin, imax, jmin, jmax;
	double n[3], v0[3], v1[3], v2[3], x[3], s, t, z, den;
	int f_v0, f_v1, f_v2;

	if (f < posvis_nf) {

		f_v0 = verts->f[f].v[0];
		f_v1 = verts->f[f].v[1];
		f_v2 = verts->f[f].v[2];

		/* Get the normal to this facet in body-fixed (asteroid) coordinates
		 * and convert it to observer coordinates     */
		for (i = 0; i <= 2; i++)
			n[i] = verts->f[f].n[i];

		dev_cotrans3(n, oa, n, 1);

		/* Consider this facet further only if its normal points somewhat
		 * towards the observer rather than away         */
		if (n[2] > 0.0) {
			/* Convert the three sets of vertex coordinates from body to ob-
			 * server coordinates; orbit_offset is the center-of-mass offset
			 * (in observer coordinates) for this model at this frame's epoch
			 * due to orbital motion, in case the model is half of a binary
			 * system.  */
			dev_cotrans3(v0, oa, verts->v[f_v0].x, 1);
			dev_cotrans3(v1, oa, verts->v[f_v1].x, 1);
			dev_cotrans3(v2, oa, verts->v[f_v2].x, 1);
			for (i = 0; i <= 2; i++) {
				v0[i] += orbit_offset[i];
				v1[i] += orbit_offset[i];
				v2[i] += orbit_offset[i];
			}

			/* Find rectangular region (in POS pixels) containing the projected
			 * facet - use floats in case model has illegal parameters and the
			 * pixel numbers exceed the limits for valid integers                         */
			iminflt[f] = floor(MIN(v0[0],MIN(v1[0],v2[0])) / pos->km_per_pixel
							- SMALLVAL + 0.5);
			imaxflt[f] = floor(MAX(v0[0],MAX(v1[0],v2[0])) / pos->km_per_pixel
							+ SMALLVAL + 0.5);
			jminflt[f] = floor(MIN(v0[1],MIN(v1[1],v2[1])) / pos->km_per_pixel
							- SMALLVAL + 0.5);
			jmaxflt[f] = floor(MAX(v0[1],MAX(v1[1],v2[1])) / pos->km_per_pixel
							+ SMALLVAL + 0.5);
			imin = (iminflt[f] < INT_MIN) ? INT_MIN : (int) iminflt[f];
			imax = (imaxflt[f] > INT_MAX) ? INT_MAX : (int) imaxflt[f];
			jmin = (jminflt[f] < INT_MIN) ? INT_MIN : (int) jminflt[f];
			jmax = (jmaxflt[f] > INT_MAX) ? INT_MAX : (int) jmaxflt[f];

			/* Set the outbnd flag if the facet extends beyond the POS window */
			if ((imin < (-pos->n)) || (imax > pos->n) || (jmin < (-pos->n))
					|| (jmax > pos->n))
				doutbnd = 1;

			/* Figure out if facet projects at least partly within POS window;
			 * if it does, look at each "contained" POS pixel and get the
			 * z-coordinate and cos(scattering angle)           */
			i1 = MAX(imin, -pos->n);		i2 = MIN(imax, pos->n);
			j1 = MAX(jmin, -pos->n);		j2 = MIN(jmax, pos->n);

			if (i1 > pos->n || i2 < -pos->n || j1 > pos->n || j2 < -pos->n) {

				/* Might need to set a flag here and exit the kernel */
				/* Facet is entirely outside the POS frame: just keep track of
				 * changed POS region     */
				//dev_POSrect2(pos, src, imin_dbl, imax_dbl, jmin_dbl, jmax_dbl);
			} else {

				/* Facet is at least partly within POS frame: find all POS
				 * pixels whose centers project onto this facet  */
				int pxa;	// pixel address in 1D float array zzf etc
				for (i = i1; i <= i2; i++) {
					x[0] = i * pos->km_per_pixel;
					for (j = j1; j <= j2; j++) {
						x[1] = j * pos->km_per_pixel;

						/* Calculate the pixel address for 1D arrays */
						pxa = (j+pos->n)*(2*pos->n+1)+(i+pos->n);

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
								z = v0[2] + s * (v1[2] - v0[2])
										+ t * (v2[2] - v1[2]);

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

								float old = atomicMaxf(&zz[pxa], z);

//								if ((z > zz[i][j]) || (fac[i][j] < 0)) {
								if (old < z || fac[i][j] < 0) {

									/* Next line assigns distance of POS pixel
									 * center from COM towards Earth; that is,
									 * by changing zz,it changes pos->z or
									 * pos->zill                */
									/* following line is a first time z calc
									 * for this pixel  */
									if (fac[i][j] < 0) atomicExch(&zz[pxa], z);
//									zz[i][j] = z;
									if (smooth) {

										/* Get smoothed version of facet unit
										 * normal: Take the linear combination
										 * of the three vertex normals; trans-
										 * form from body to observer coordina-
										 * tes; and make sure that it points
										 * somewhat in our direction.         */
										for (k = 0; k <= 2; k++)
											n[k] =	verts->v[f_v0].n[k]
											 + s * (verts->v[f_v1].n[k] - verts->v[f_v0].n[k])
											 + t * (verts->v[f_v2].n[k]	- verts->v[f_v1].n[k]);
										dev_cotrans3(n, oa, n, 1);
										dev_normalize(n);
									}

									/* Determine scattering and/or incidence
									 * angles. Next lines change pos->cose/
									 * cosill. If bistatic (lightcurves), where
									 * we are viewing from Earth (src = 0),
									 * pos->cosi is also changed.                 */
									if (n[2] > 0.0) {
										atomicExch(&cosa[pxa], n[2]);
										if ((!src) && (pos->bistatic)) {
											float temp = (float)dev_dot(n,usrc);
											atomicExch(&cosb[pxa], temp);
											if (cosb[pxa] <= 0.0)
												cosa[pxa] = 0.0;
										}
									}

									/* This now happens in a separate kernel */
									/*  Keep track of the changed POS region  */
//									dev_POSrect2(pos, src, (double) i,
//											(double) i, (double) j, (double) j);

									/* Next lines change pos->body/bodyill,
									 * pos->comp/compill, pos->f/fill          */
									bod[i][j] = body;
									cmp[i][j] = comp;
									fac[i][j] = f;

								} /* end if (no other facet yet blocks this facet from view) */
							} /* end if 0 <= t <= s (facet center is "in" this POS pixel) */
						} /* end if 0 <= s <= 1 */
					} /* end j-loop over POS rows */
				} /* end i-loop over POS columns */
			} /* end else of if (i1 > pos->n || i2 < -pos->n || j1 > pos->n || j2 < -pos->n) */
		} /* End if (n[2] > 0.0) */
	} /* end if (f < nf) */
}
__global__ void posvis_set_logfactor() {
	/* Single-threaded kernel */
	/*  If the target extends beyond the POS frame, set pos->posbnd_logfactor equal
	 to the logarithm of the ratio of the number of pixels in a frame extended to
	 include the entire target, divided by the number of pixels in the actual frame  */
	if (threadIdx.x == 0) {
		if (doutbnd) {
			double xfactor, yfactor;
			xfactor = (MAX(imax_overall,pos->n) - MIN(imin_overall, -pos->n)
							+ 1) / (2 * pos->n + 1);
			yfactor = (MAX(jmax_overall,pos->n) - MIN(jmin_overall, -pos->n)
							+ 1) / (2 * pos->n + 1);
			pos->posbnd_logfactor = log(xfactor * yfactor);
		}
	}
}
__global__ void posvis_set_logfactor_sep_krnl(float *minmax_overall) {
	/* Single-threaded kernel */
	/*  If the target extends beyond the POS frame, set pos->posbnd_logfactor equal
	 to the logarithm of the ratio of the number of pixels in a frame extended to
	 include the entire target, divided by the number of pixels in the actual frame  */
	/* minmax_overall[0] = imin_overall
	 * minmax_overall[1] = imax_overall
	 * minmax_overall[2] = jmin_overall
	 * minmax_overall[3] = jmax_overall
	 */
	if (threadIdx.x == 0) {
		if (doutbnd) {
			double xfactor, yfactor;
			xfactor = (MAX(minmax_overall[1],pos->n) - MIN(minmax_overall[0], -pos->n)
							+ 1) / (2 * pos->n + 1);
			yfactor = (MAX(minmax_overall[3],pos->n) - MIN(minmax_overall[2], -pos->n)
							+ 1) / (2 * pos->n + 1);
			pos->posbnd_logfactor = log(xfactor * yfactor);
		}
	}
}

__host__ int posvis_cuda_2(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, double horbit_offset[3], int set, int frame,
		int src, int body, int comp) {
	int nf, outbnd, n;
	dim3 BLK,THD;
	cudaEvent_t start, stop;
	float milliseconds;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	int dbg_occ = 0;

	/* The xxxxflt arrays below are for the parallel reduction to find xlim
	 * and ylim for each frame's POS. The minmax_overall array keeps the
	 * absolute min/maxes for i and j regardless of POS frame limits - that's
	 * the old imin_overall/imax_overall/jmin_overall/jmax_overall values  */
	float *iminflt, *imaxflt, *jminflt, *jmaxflt, *minmax_overall;

	/* Launch init kernel */
	posvis_init_krnl<<<1,1>>>(dpar, dmod, ddat, set, frame, src,
			horbit_offset[0], horbit_offset[1],horbit_offset[2]);
	checkErrorAfterKernelLaunch("posvis_init_krnl, line ");
	gpuErrchk(cudaMemcpyFromSymbol(&nf, posvis_nf, sizeof(int), 0,
				cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&n, psvs_n, sizeof(int), 0,
				cudaMemcpyDeviceToHost));

	if (POSVIS_SEPARATE) {
		cudaCalloc1((void**)&minmax_overall, sizeof(float), 4);
		cudaCalloc1((void**)&iminflt, 		sizeof(float), nf);
		cudaCalloc1((void**)&imaxflt, 		sizeof(float), nf);
		cudaCalloc1((void**)&jminflt, 		sizeof(float), nf);
		cudaCalloc1((void**)&jmaxflt, 		sizeof(float), nf);
	}

	/* Configure and launch the facet kernel. Which kernel gets launched
	 * depends on flags DynProc (dynamic processing) and POSVIS_SEPARATE
	 * (xlim and ylim get calculated via a separate parallel reduction in an
	 * additional kernel after the facet kernel) */
	//npixels = (2*n+1)*(2*n+1);
	BLK.x = floor((maxThreadsPerBlock - 1 + nf) / maxThreadsPerBlock);
	THD.x = maxThreadsPerBlock;
	if (!DYNPROC) {
		if (POSVIS_SEPARATE) {
			/* Launch facet kernel that doesn't calculate/set pos->xlim/ylim */
			posvis_facet_sep_krnl<<<BLK,THD>>>(src, body, comp, iminflt,
					imaxflt, jminflt, jmaxflt);
			checkErrorAfterKernelLaunch("posvis_facet_sep_krnl");
			gpuErrchk(cudaMemcpyFromSymbol(&outbnd, doutbnd, sizeof(outbnd), 0,
							cudaMemcpyDeviceToHost));

			/* Now that the imin/imax/jmin/jmax arrays are loaded, proceed
			 * to the parallel reduction to find the max/mins in them and
			 * update pos->xlim and pos->ylim accordingly			 */
			compute_xlim_ylim(ddat, nf,	set, frame, src, iminflt, imaxflt,
					jminflt, jmaxflt, minmax_overall);

			/* Call modified posvis_set_logfactor kernel */
			posvis_set_logfactor_sep_krnl<<<1,1>>>(minmax_overall);
			checkErrorAfterKernelLaunch("posvis_set_logfactor_sep_krnl");

			cudaFree(iminflt);
			cudaFree(imaxflt);
			cudaFree(jminflt);
			cudaFree(jmaxflt);
			cudaFree(minmax_overall);
		}
		else {

			if (TIMING)	cudaEventRecord(start);

			/* Launch the non-dynamic processing facet kernel for
			 * Nvidia GPUs w/Compute Capability of less than 3.5 (slower) */
			posvis_facet_krnl<<<BLK,THD>>>(src, body, comp);

			if (TIMING) {
				cudaEventRecord(stop);
				cudaEventSynchronize(stop);
				milliseconds = 0;
				cudaEventElapsedTime(&milliseconds, start, stop);
				gpuErrchk(cudaMemcpyFromSymbol(&dbg_occ, dbg_pxl_occ, sizeof(int), 0,
								cudaMemcpyDeviceToHost));
				printf("%i facets in posvis_facet_krnl completed in %3.3f ms with %i occurrences.\n", nf, milliseconds, dbg_occ);
			}

			checkErrorAfterKernelLaunch("posvis_facet_krnl, line ");
			gpuErrchk(cudaMemcpyFromSymbol(&outbnd, doutbnd, sizeof(outbnd), 0,
					cudaMemcpyDeviceToHost));
		}
	}
	if (DYNPROC) {
		/* Configure and launch the dynamic processing facet kernel for Nvidia
		 * GPUs with Compute Capability of 3.5 or higher (faster).   */
		posvis_fct_dynp_krnl<<<BLK,THD>>>(src, body, comp);
		checkErrorAfterKernelLaunch("posvis_fct_dynp_krnl, line ");
		gpuErrchk(cudaMemcpyFromSymbol(&outbnd, doutbnd, sizeof(outbnd), 0,
				cudaMemcpyDeviceToHost));

		/* Launch single-threaded kernel to set the pos logfactor */
		posvis_set_logfactor<<<1,1>>>();
		checkErrorAfterKernelLaunch("posvis_set_logfactor_krnl, line ");
	}

	return outbnd;
}

