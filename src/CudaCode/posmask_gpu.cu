extern "C" {
#include "../shape/head.h"
}

__global__ void posmask_init_krnl(struct pos_t **pos, double3 *so,
		float *pixels_per_km, int f) {
	/* This single-threaded kernel performs the first few tasks (outside the
	 * pixel loop) of routine posmask.	 */
	if (threadIdx.x == 0) {
		dev_mtrnsps2(so, pos[f]->oe, f);
		dev_mmmul2(so, pos[f]->se, so, f);
		pixels_per_km[f] = 1/pos[f]->km_per_pixel;
	}
}
__global__ void posmask_init_krnl2(struct pos_t **pos, double3 *so,
		float *pixels_per_km, int size) {
	/* nfrm_alloc-threaded kernel */
	int f = blockIdx.x * blockDim.x + threadIdx.x;
	if (f < size) {
		dev_mtrnsps2(so, pos[f]->oe, f);
		dev_mmmul2(so, pos[f]->se, so, f);
		pixels_per_km[f] = 1/pos[f]->km_per_pixel;
	}
}
__global__ void posmask_init_mgpu_krnl(struct pos_t **pos, double3 *so,
		float *pixels_per_km, int size, int oddflg) {

	/* nfrm_half0/nfrm_half1-threaded kernel for dual-GPU operation. It
	 * performs the first few tasks (outside the pixel loop) of routine
	 * posmask.	 */

	int hf = blockIdx.x * blockDim.x + threadIdx.x;
	int f = 2*hf + oddflg + 1;

	if (hf < size) {
		dev_mtrnsps2(so, pos[hf]->oe, hf);
		dev_mmmul2(so, pos[hf]->se, so, hf);
		pixels_per_km[hf] = 1/pos[hf]->km_per_pixel;
	}
}

__global__ void posmask_krnl(
		struct par_t *dpar,
		struct pos_t **pos,
		double3 *so,
		float *pixels_per_km,
		int *posn,
		int nThreads,
		int xspan,
		int f)
{
	/* multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int n = posn[f];
	int i = offset % xspan - n;
	int j = offset / xspan - n;
	float tol = dpar->mask_tol;
	float kmpxl = (float)pos[f]->km_per_pixel;
	int im, jm, i1, j1, i2, j2, i_sign, j_sign, pxa, pxa1, pxa2;
	float3 xk;
	float i0_f, j0_f, zill, t, u, bignum;
	bignum = 0.99*HUGENUMBER;  /* z = -HUGENUMBER for blank-sky pixels */

	/*  Loop through all POS pixels  */
	if (offset < nThreads) {

		if (pos[f]->cose_s[offset] != 0.0) {     /* if there's something there */
			xk.x = i*kmpxl;     /* calculate 3D position */
			xk.y = j*kmpxl;
			xk.z = pos[f]->z_s[offset];

			/* Given the observer coordinates x of of POS pixel (i,j), find
			 * which pixel (im,jm) this corresponds to in the projected view as
			 * seen from the source (sun)             */
			dev_cotrans7(&xk, so, xk, 1, f);	/* go into source coordinates */
			i0_f = xk.x*pixels_per_km[f];     /* unrounded (double precision) */
			j0_f = xk.y*pixels_per_km[f];
			im = dev_vp_iroundf(i0_f);            /* center of nearest pixel in mask */
			jm = dev_vp_iroundf(j0_f);
			pxa = (jm+posn[f])*xspan + (im+posn[f]); /* The float single pointer address */

			/* If center of projected pixel "seen" from source (as determined
			 * by routine posvis) lies within the boundaries of the mask,
			 * projects onto model rather than onto blank space, and represents
			 * a body, component, and facet different from those seen in the
			 * POS, calculate distance from mask pixel to source and compare to
			 * distance from POS pixel to source.                             */
			if (fabs(i0_f) < n && fabs(j0_f) < n
					&& pos[f]->zill_s[pxa] > -bignum
					&&(pos[f]->f[i][j]    != pos[f]->fill[im][jm]    ||
							pos[f]->comp[i][j] != pos[f]->compill[im][jm] ||
							pos[f]->body[i][j] != pos[f]->bodyill[im][jm]    )) {

				/* Rather than using distance towards source of CENTER of mask
				 * pixel, use bilinear interpolation to get distance towards
				 * source where the line between source and POS pixel's center
				 * intersects the mask pixel.                                */
				i1 = (int) floor( i0_f);
				j1 = (int) floor( j0_f);
				pxa = (j1+posn[f])*xspan + (i1+posn[f]);
				pxa1 = (j1+posn[f]+1)*xspan + (i1+posn[f]);

				if (pos[f]->zill_s[pxa]     > -bignum &&
						pos[f]->zill_s[pxa+1]   > -bignum &&
						pos[f]->zill_s[pxa1]   > -bignum &&
						pos[f]->zill_s[pxa1+1] > -bignum    ) {

					/* Do standard bilinear interpolation: None of the four
					 * surrounding "grid square" pixels in the mask is
					 * blank sky                           */
					t = i0_f - i1;
					u = j0_f - j1;
					zill = (1 - t)*(1 - u)*pos[f]->zill_s[pxa]
					      + t*(1 - u)*pos[f]->zill_s[pxa+1]
					      + t*u*pos[f]->zill_s[pxa1+1]
					      + (1 - t)*u*pos[f]->zill_s[pxa1];
				} else {
					/* The following code block is a kludge: One or more of the
					 * four surrounding "grid square" pixels in mask is blank
					 * sky, so standard bilinear interpolation won't work.  */
					pxa = (jm+posn[f])*xspan + (im+posn[f]);
					zill = pos[f]->zill_s[pxa];

					i_sign = (i0_f >= im) ? 1 : -1;
					i2 = im + i_sign;
					pxa1 = (jm+posn[f])*xspan + (i2+posn[f]);

					if (abs(i2) <= n && pos[f]->zill_s[pxa1] > -bignum) {
						zill += fabs(i0_f - im)
           				  * (pos[f]->zill_s[pxa1] - pos[f]->zill_s[pxa]);
					} else {
						i2 = im - i_sign;
						pxa1 = (jm+posn[f])*xspan + (i2+posn[f]);
						if (abs(i2) <= n && pos[f]->zill_s[pxa1] > -bignum)
							zill -= fabs(i0_f - im)
							* (pos[f]->zill_s[pxa1] - pos[f]->zill_s[pxa]);
					}

					j_sign = (j0_f >= jm) ? 1 : -1;
					j2 = jm + j_sign;
					pxa1 = (j2+posn[f])*xspan + (im+posn[f]);

					if (abs(j2) <= n && pos[f]->zill_s[pxa1] > -bignum) {
						zill += fabs(j0_f - jm)
                        	  * (pos[f]->zill_s[pxa1] - pos[f]->zill_s[pxa]);
					} else {
						j2 = jm - j_sign;
						pxa1 = (j2+posn[f])*xspan + (im+posn[f]);

						if (abs(j2) <= n && pos[f]->zill_s[pxa1] > -bignum)
							zill -= fabs(j0_f - jm)
							* (pos[f]->zill_s[pxa1] - pos[f]->zill_s[pxa]);
					}
				}
				/* If interpolated point within mask pixel is at least tol km
				 * closer to source than is the center of POS pixel, the facet
				 * represented by the mask pixel is shadowing the POS pixel:
				 * represent this by setting
				 * 		cos(scattering angle) = 0.0 for the POS pixel.      */
				if (zill - xk.z > tol)
					pos[f]->cose_s[offset] = 0.0;
			}
		}
	}
}
