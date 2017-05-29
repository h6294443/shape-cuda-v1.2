/*****************************************************************************************
                                                                              mkparlist.c

Takes mod_t and dat_t files and makes a list of pointers to the free parameters.

Modified 2016 Dec 20 by ME:
	Adapted for use with CUDA.  This version of mkparlist is used exclusively
	for structures residing in device (GPU) memory.

Modified 2014 Aug 22 by SN:
    Add spin.lib_amp, spin.lib_freq, spin.lib_phase for taking into account librations.

Modified 2014 February 12 by CM:
    Implement multiple radar and optical scattering laws

Modified 2013 May 20 by CM:
    Implement ovoid shape components

Modified 2012 March 23 by CM:
    Add the "dopscale" parameter for delay-Doppler and Doppler datasets

Modified 2012 March 14 by CM:
    Have separate code blocks for delay-Doppler vs. Doppler datasets when handling the
        delay correction polynomial coefficients

Modified 2011 September 2 by CM:
    Add the "harmlambert" and "inholambert" optical scattering laws

Modified 2011 August 7 by CM:
    Add spin impulses

Modified 2010 June 1 by CM:
    Change "scalefactor" parameter from a scalar to a 3-component vector
        for harmonic and vertex shape structures

Modified 2010 April 27 by CM:
    Add "tabular" radar scattering law

Modified 2008 August 10 by CM:
    Change "delcor_step" and "delcor_abstol" parameters to be vectors
        rather than scalars
    Change parameter type for each component's linear offsets in the
        mod file from SHAPEPAR to SIZEPAR

Modified 2006 October 1 by CM:
    Replace ellipsoid diameters D with two_a, a_over_b, and b_over_c
    Add "scalefactor" to harmonic and vertex shape structures
    Add "ratio_step" "ratio_tol" and "ratio_abstol" fit parameters
    Add SIZEPAR parameters

Modified 2006 March 6 by PT:
    Add "spin.omegadot" parameter for changing spin rate

Modified 2005 September 7 by CM:
    Add "harmlommel" "harmhapke" and "harmkaas" optical scattering laws
    Add "harmcosine" radar scattering law

Modified 2005 August 1 by CM:
    Add "inhokaas" optical scattering law

Modified 2005 July 20 by CM:
    Add "hagfors" and "cosine_qs" and "gauss+cosine" and "hagfors+cosine"
        and "cosine+cosine" and inhomogeneous "inhocosine" radar
        scattering laws
    Eliminate "flat" radar scattering law

Modified 2005 July 4 by CM:
    Add inhomogeneous "inholommel" and "inhohapke" optical scattering laws

Modified 2005 February 24 by CM:
    Add POS case for data parameters so that the horizontal and vertical
        offsets for plane-of-sky datasets can be fitted

Modified 2004 August 13 by CM:
    Assign "fparabstol" pointers to handle absolute fitting tolerances

Modified 2004 April 29 by CM:
    For Kaasalainen scattering law, switch from weighting factor "c"
        (ranging from 0 to infinity) to "wt" (ranging from 0 to 1)

Modified 2004 March 13 by CM:
    Add code to handle harmonic components

Modified 2004 February 25 by CM:
    Add Kaasalainen "Lambert + Lommel-Seeliger" scattering law parameters

Modified 2003 April 3 by CM:
    1) Fit ellipsoid diameters using length_step and length_tol
       rather than angle_step and angle_tol
    2) Fit photometric parameters using photo_step and photo_tol
       rather than length_step and length_tol
    3) Fit moments of inertia using inertia_step and inertia_tol
       rather than spin_step and spin_tol
    4) Fit delay correction polynomial coefficients using
       delcor_step and delcor_tol rather than spin_step and spin_tol
 *****************************************************************************************/
extern "C" {
#include "../shape/head.h"
}

__device__ int p;

__global__ void mpl_comp_krnl(struct par_t *dpar, struct mod_t *dmod,
		double **fpntr, double *fparstep, double *fpartol, double *fparabstol,
		int *fpartype) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		p = -1;
		int i = 0, j, k; /* component index - always zero */
		for (j=0; j<=2; j++) {     /* check linear offsets */
			if (dmod->shape.comp[i].off[j].state == 'f') {
				fpntr[++p] = &dmod->shape.comp[i].off[j].val;
				fparstep[p] = dpar->length_step;
				fpartol[p] = dpar->length_tol;
				fparabstol[p] = dpar->length_abstol;
				fpartype[p] = SIZEPAR;
			}
		}
		for (j=0; j<=2; j++) {    /* check angular offsets */
			if (dmod->shape.comp[i].rot[j].state == 'f') {
				fpntr[++p] = &dmod->shape.comp[i].rot[j].val;
				fparstep[p] = dpar->angle_step;
				fpartol[p] = dpar->angle_tol;
				fparabstol[p] = dpar->angle_abstol;
				fpartype[p] = SHAPEPAR;
			}
		}

		switch (dmod->shape.comp[i].type) {
		case ELLIPSE:
			if (dmod->shape.comp[i].desc.ell.two_a.state == 'f') {
				fpntr[++p] = &dmod->shape.comp[i].desc.ell.two_a.val;
				fparstep[p] = dpar->length_step;
				fpartol[p] = dpar->length_tol;
				fparabstol[p] = dpar->length_abstol;
				fpartype[p] = SIZEPAR;
			}
			if (dmod->shape.comp[i].desc.ell.a_over_b.state == 'f') {
				fpntr[++p] = &dmod->shape.comp[i].desc.ell.a_over_b.val;
				fparstep[p] = dpar->ratio_step;
				fpartol[p] = dpar->ratio_tol;
				fparabstol[p] = dpar->ratio_abstol;
				fpartype[p] = SHAPEPAR;
			}
			if (dmod->shape.comp[i].desc.ell.b_over_c.state == 'f') {
				fpntr[++p] = &dmod->shape.comp[i].desc.ell.b_over_c.val;
				fparstep[p] = dpar->ratio_step;
				fpartol[p] = dpar->ratio_tol;
				fparabstol[p] = dpar->ratio_abstol;
				fpartype[p] = SHAPEPAR;
			}
			break;
		case OVOID:
			if (dmod->shape.comp[i].desc.ovoid.two_a.state == 'f') {
				fpntr[++p] = &dmod->shape.comp[i].desc.ovoid.two_a.val;
				fparstep[p] = dpar->length_step;
				fpartol[p] = dpar->length_tol;
				fparabstol[p] = dpar->length_abstol;
				fpartype[p] = SIZEPAR;
			}
			if (dmod->shape.comp[i].desc.ovoid.a_over_b.state == 'f') {
				fpntr[++p] = &dmod->shape.comp[i].desc.ovoid.a_over_b.val;
				fparstep[p] = dpar->ratio_step;
				fpartol[p] = dpar->ratio_tol;
				fparabstol[p] = dpar->ratio_abstol;
				fpartype[p] = SHAPEPAR;
			}
			if (dmod->shape.comp[i].desc.ovoid.b_over_c.state == 'f') {
				fpntr[++p] = &dmod->shape.comp[i].desc.ovoid.b_over_c.val;
				fparstep[p] = dpar->ratio_step;
				fpartol[p] = dpar->ratio_tol;
				fparabstol[p] = dpar->ratio_abstol;
				fpartype[p] = SHAPEPAR;
			}
			if (dmod->shape.comp[i].desc.ovoid.k.state == 'f') {
				fpntr[++p] = &dmod->shape.comp[i].desc.ovoid.k.val;
				fparstep[p] = dpar->ratio_step;
				fpartol[p] = dpar->ratio_tol;
				fparabstol[p] = dpar->ratio_abstol;
				fpartype[p] = SHAPEPAR;
			}
			break;
		case HARMONIC:
			for (j=0; j<=2; j++)
				if (dmod->shape.comp[i].desc.har.scalefactor[j].state == 'f') {
					fpntr[++p] = &dmod->shape.comp[i].desc.har.scalefactor[j].val;
					fparstep[p] = dpar->ratio_step;
					fpartol[p] = dpar->ratio_tol;
					fparabstol[p] = dpar->ratio_abstol;
					fpartype[p] = SIZEPAR;
				}
			for (j=0; j<=dmod->shape.comp[i].desc.har.nhar; j++) {
				if (dmod->shape.comp[i].desc.har.a[j][0].state == 'f') {
					fpntr[++p] = &dmod->shape.comp[i].desc.har.a[j][0].val;
					fparstep[p] = dpar->length_step;
					fpartol[p] = dpar->length_tol;
					fparabstol[p] = dpar->length_abstol;
					fpartype[p] = SHAPEPAR;
				}
				for (k=1; k<=j; k++) {
					if (dmod->shape.comp[i].desc.har.a[j][k].state == 'f') {
						fpntr[++p] = &dmod->shape.comp[i].desc.har.a[j][k].val;
						fparstep[p] = dpar->length_step;
						fpartol[p] = dpar->length_tol;
						fparabstol[p] = dpar->length_abstol;
						fpartype[p] = SHAPEPAR;
					}
					if (dmod->shape.comp[i].desc.har.b[j][k].state == 'f') {
						fpntr[++p] = &dmod->shape.comp[i].desc.har.b[j][k].val;
						fparstep[p] = dpar->length_step;
						fpartol[p] = dpar->length_tol;
						fparabstol[p] = dpar->length_abstol;
						fpartype[p] = SHAPEPAR;
					}
				}
			}
			break;
		case VERTEX:
			for (j=0; j<=2; j++)
				if (dmod->shape.comp[i].desc.ver.scalefactor[j].state == 'f') {
					fpntr[++p] = &dmod->shape.comp[i].desc.ver.scalefactor[j].val;
					fparstep[p] = dpar->ratio_step;
					fpartol[p] = dpar->ratio_tol;
					fparabstol[p] = dpar->ratio_abstol;
					fpartype[p] = SIZEPAR;
				}
			for (j=0; j<dmod->shape.comp[i].desc.ver.nv; j++) {
				if (dmod->shape.comp[i].desc.ver.v[j].r.state == 'f') {
					fpntr[++p] = &dmod->shape.comp[i].desc.ver.v[j].r.val;
					fparstep[p] = dpar->length_step;
					fpartol[p] = dpar->length_tol;
					fparabstol[p] = dpar->length_abstol;
					fpartype[p] = SHAPEPAR;
				}
			}
			break;
		default:
			printf("mkparlist_cuda.cu: can't do that type of model yet\n");
		}
	}
}
__global__ void mpl_rad_krnl(struct par_t *dpar, struct mod_t *dmod,
		double **fpntr, double *fparstep, double *fpartol, double *fparabstol,
		int *fpartype) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		int ilaw = 0;	/* component index = always zero for now */
		int i, c, f, m, l;
		switch (dmod->photo.radtype[ilaw]) {
		case COSINELAW_DIFF:
			if (dmod->photo.radar[ilaw].RC.R.state == 'f') {
				fpntr[++p] = &dmod->photo.radar[ilaw].RC.R.val;
				fparstep[p] = dpar->photo_step;
				fpartol[p] = dpar->photo_tol;
				fparabstol[p] = dpar->photo_abstol;
				fpartype[p] = PHOTOPAR;
			}
			if (dmod->photo.radar[ilaw].RC.C.state == 'f') {
				fpntr[++p] = &dmod->photo.radar[ilaw].RC.C.val;
				fparstep[p] = dpar->photo_step;
				fpartol[p] = dpar->photo_tol;
				fparabstol[p] = dpar->photo_abstol;
				fpartype[p] = PHOTOPAR;
			}
			break;
		case TABULARLAW:
			for (i=0; i<dmod->photo.radar[ilaw].tabular.n; i++) {
				if (dmod->photo.radar[ilaw].tabular.rho[i].state == 'f') {
					fpntr[++p] = &dmod->photo.radar[ilaw].tabular.rho[i].val;
					fparstep[p] = dpar->photo_step;
					fpartol[p] = dpar->photo_tol;
					fparabstol[p] = dpar->photo_abstol;
					fpartype[p] = PHOTOPAR;
				}
			}
			break;
		case HARMCOSINE_DIFF:
			for (l=0; l<=dmod->photo.radar[ilaw].harmcosine.R.nhar; l++)
				for (m=0; m<=l; m++) {
					if (dmod->photo.radar[ilaw].harmcosine.R.a[l][m].state == 'f') {
						fpntr[++p] = &dmod->photo.radar[ilaw].harmcosine.R.a[l][m].val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
					if (m > 0 && dmod->photo.radar[ilaw].harmcosine.R.b[l][m].state == 'f') {
						fpntr[++p] = &dmod->photo.radar[ilaw].harmcosine.R.b[l][m].val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
				}
			for (l=0; l<=dmod->photo.radar[ilaw].harmcosine.C.nhar; l++)
				for (m=0; m<=l; m++) {
					if (dmod->photo.radar[ilaw].harmcosine.C.a[l][m].state == 'f') {
						fpntr[++p] = &dmod->photo.radar[ilaw].harmcosine.C.a[l][m].val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
					if (m > 0 && dmod->photo.radar[ilaw].harmcosine.C.b[l][m].state == 'f') {
						fpntr[++p] = &dmod->photo.radar[ilaw].harmcosine.C.b[l][m].val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
				}
			break;
		case GAUSSIANLAW :
		case HAGFORSLAW  :
		case COSINELAW_QS:
			if (dmod->photo.radar[ilaw].quasispec.R.state == 'f') {
				fpntr[++p] = &dmod->photo.radar[ilaw].quasispec.R.val;
				fparstep[p] = dpar->photo_step;
				fpartol[p] = dpar->photo_tol;
				fparabstol[p] = dpar->photo_abstol;
				fpartype[p] = PHOTOPAR;
			}
			if (dmod->photo.radar[ilaw].quasispec.C.state == 'f') {
				fpntr[++p] = &dmod->photo.radar[ilaw].quasispec.C.val;
				fparstep[p] = dpar->photo_step;
				fpartol[p] = dpar->photo_tol;
				fparabstol[p] = dpar->photo_abstol;
				fpartype[p] = PHOTOPAR;
			}
			break;
		case GAUSSIAN_COSINE:
		case HAGFORS_COSINE :
		case COSINE_COSINE  :
			if (dmod->photo.radar[ilaw].hybrid.qs.R.state == 'f') {
				fpntr[++p] = &dmod->photo.radar[ilaw].hybrid.qs.R.val;
				fparstep[p] = dpar->photo_step;
				fpartol[p] = dpar->photo_tol;
				fparabstol[p] = dpar->photo_abstol;
				fpartype[p] = PHOTOPAR;
			}
			if (dmod->photo.radar[ilaw].hybrid.qs.C.state == 'f') {
				fpntr[++p] = &dmod->photo.radar[ilaw].hybrid.qs.C.val;
				fparstep[p] = dpar->photo_step;
				fpartol[p] = dpar->photo_tol;
				fparabstol[p] = dpar->photo_abstol;
				fpartype[p] = PHOTOPAR;
			}
			if (dmod->photo.radar[ilaw].hybrid.diff.R.state == 'f') {
				fpntr[++p] = &dmod->photo.radar[ilaw].hybrid.diff.R.val;
				fparstep[p] = dpar->photo_step;
				fpartol[p] = dpar->photo_tol;
				fparabstol[p] = dpar->photo_abstol;
				fpartype[p] = PHOTOPAR;
			}
			if (dmod->photo.radar[ilaw].hybrid.diff.C.state == 'f') {
				fpntr[++p] = &dmod->photo.radar[ilaw].hybrid.diff.C.val;
				fparstep[p] = dpar->photo_step;
				fpartol[p] = dpar->photo_tol;
				fparabstol[p] = dpar->photo_abstol;
				fpartype[p] = PHOTOPAR;
			}
			break;
		case INHOCOSINE_DIFF:
			if (dmod->photo.radar[ilaw].inhocosine.global.R.state == 'f') {
				fpntr[++p] = &dmod->photo.radar[ilaw].inhocosine.global.R.val;
				fparstep[p] = dpar->photo_step;
				fpartol[p] = dpar->photo_tol;
				fparabstol[p] = dpar->photo_abstol;
				fpartype[p] = PHOTOPAR;
			}
			if (dmod->photo.radar[ilaw].inhocosine.global.C.state == 'f') {
				fpntr[++p] = &dmod->photo.radar[ilaw].inhocosine.global.C.val;
				fparstep[p] = dpar->photo_step;
				fpartol[p] = dpar->photo_tol;
				fparabstol[p] = dpar->photo_abstol;
				fpartype[p] = PHOTOPAR;
			}
			for (c=0; c<dmod->shape.ncomp; c++)
				for (f=0; f<dmod->shape.comp[c].real.nf; f++) {
					if (dmod->photo.radar[ilaw].inhocosine.local[c][f].R.state == 'f') {
						fpntr[++p] = &dmod->photo.radar[ilaw].inhocosine.local[c][f].R.val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
					if (dmod->photo.radar[ilaw].inhocosine.local[c][f].C.state == 'f') {
						fpntr[++p] = &dmod->photo.radar[ilaw].inhocosine.local[c][f].C.val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
				}
			break;
		case NOLAW:
			break;
		default:
			printf("mkparlist_cuda.cu: can't do that radar law yet\n");
		}
	}
}
__global__ void mpl_photo_krnl(struct par_t *dpar, struct mod_t *dmod,
		double **fpntr, double *fparstep, double *fpartol, double *fparabstol,
		int *fpartype) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		int c, f, l, m;
		int ilaw = 0;
		switch (dmod->photo.opttype[ilaw]) {
		case NOLAW:
			break;
		case GEOMETRICAL:
		case LAMBERTLAW:
		case LOMMEL:
			if (dmod->photo.optical[ilaw].R.R.state == 'f') {
				fpntr[++p] = &dmod->photo.optical[ilaw].R.R.val;
				fparstep[p] = dpar->photo_step;
				fpartol[p] = dpar->photo_tol;
				fparabstol[p] = dpar->photo_abstol;
				fpartype[p] = PHOTOPAR;
			}
			break;
		case HARMLAMBERT:
		case HARMLOMMEL:
			for (l=0; l<=dmod->photo.optical[ilaw].harmR.R.nhar; l++)
				for (m=0; m<=l; m++) {
					if (dmod->photo.optical[ilaw].harmR.R.a[l][m].state == 'f') {
						fpntr[++p] = &dmod->photo.optical[ilaw].harmR.R.a[l][m].val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
					if (m > 0 && dmod->photo.optical[ilaw].harmR.R.b[l][m].state == 'f') {
						fpntr[++p] = &dmod->photo.optical[ilaw].harmR.R.b[l][m].val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
				}
			break;
		case INHOLAMBERT:
		case INHOLOMMEL:
			if (dmod->photo.optical[ilaw].inhoR.global.R.state == 'f') {
				fpntr[++p] = &dmod->photo.optical[ilaw].inhoR.global.R.val;
				fparstep[p] = dpar->photo_step;
				fpartol[p] = dpar->photo_tol;
				fparabstol[p] = dpar->photo_abstol;
				fpartype[p] = PHOTOPAR;
			}
			for (c=0; c<dmod->shape.ncomp; c++)
				for (f=0; f<dmod->shape.comp[c].real.nf; f++)
					if (dmod->photo.optical[ilaw].inhoR.local[c][f].R.state == 'f') {
						fpntr[++p] = &dmod->photo.optical[ilaw].inhoR.local[c][f].R.val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
			break;
		case HAPKE:
			if (dmod->photo.optical[ilaw].hapke.w.state == 'f') {
				fpntr[++p] = &dmod->photo.optical[ilaw].hapke.w.val;
				fparstep[p] = dpar->photo_step;
				fpartol[p] = dpar->photo_tol;
				fparabstol[p] = dpar->photo_abstol;
				fpartype[p] = PHOTOPAR;
			}
			if (dmod->photo.optical[ilaw].hapke.h.state == 'f') {
				fpntr[++p] = &dmod->photo.optical[ilaw].hapke.h.val;
				fparstep[p] = dpar->photo_step;
				fpartol[p] = dpar->photo_tol;
				fparabstol[p] = dpar->photo_abstol;
				fpartype[p] = PHOTOPAR;
			}
			if (dmod->photo.optical[ilaw].hapke.B0.state == 'f') {
				fpntr[++p] = &dmod->photo.optical[ilaw].hapke.B0.val;
				fparstep[p] = dpar->photo_step;
				fpartol[p] = dpar->photo_tol;
				fparabstol[p] = dpar->photo_abstol;
				fpartype[p] = PHOTOPAR;
			}
			if (dmod->photo.optical[ilaw].hapke.g.state == 'f') {
				fpntr[++p] = &dmod->photo.optical[ilaw].hapke.g.val;
				fparstep[p] = dpar->photo_step;
				fpartol[p] = dpar->photo_tol;
				fparabstol[p] = dpar->photo_abstol;
				fpartype[p] = PHOTOPAR;
			}
			if (dmod->photo.optical[ilaw].hapke.theta.state == 'f') {
				fpntr[++p] = &dmod->photo.optical[ilaw].hapke.theta.val;
				fparstep[p] = dpar->photo_step;
				fpartol[p] = dpar->photo_tol;
				fparabstol[p] = dpar->photo_abstol;
				fpartype[p] = PHOTOPAR;
			}
			break;
		case HARMHAPKE:
			for (l=0; l<=dmod->photo.optical[ilaw].harmhapke.w.nhar; l++)
				for (m=0; m<=l; m++) {
					if (dmod->photo.optical[ilaw].harmhapke.w.a[l][m].state == 'f') {
						fpntr[++p] = &dmod->photo.optical[ilaw].harmhapke.w.a[l][m].val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
					if (m > 0 && dmod->photo.optical[ilaw].harmhapke.w.b[l][m].state == 'f') {
						fpntr[++p] = &dmod->photo.optical[ilaw].harmhapke.w.b[l][m].val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
				}
			for (l=0; l<=dmod->photo.optical[ilaw].harmhapke.h.nhar; l++)
				for (m=0; m<=l; m++) {
					if (dmod->photo.optical[ilaw].harmhapke.h.a[l][m].state == 'f') {
						fpntr[++p] = &dmod->photo.optical[ilaw].harmhapke.h.a[l][m].val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
					if (m > 0 && dmod->photo.optical[ilaw].harmhapke.h.b[l][m].state == 'f') {
						fpntr[++p] = &dmod->photo.optical[ilaw].harmhapke.h.b[l][m].val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
				}
			for (l=0; l<=dmod->photo.optical[ilaw].harmhapke.B0.nhar; l++)
				for (m=0; m<=l; m++) {
					if (dmod->photo.optical[ilaw].harmhapke.B0.a[l][m].state == 'f') {
						fpntr[++p] = &dmod->photo.optical[ilaw].harmhapke.B0.a[l][m].val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
					if (m > 0 && dmod->photo.optical[ilaw].harmhapke.B0.b[l][m].state == 'f') {
						fpntr[++p] = &dmod->photo.optical[ilaw].harmhapke.B0.b[l][m].val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
				}
			for (l=0; l<=dmod->photo.optical[ilaw].harmhapke.g.nhar; l++)
				for (m=0; m<=l; m++) {
					if (dmod->photo.optical[ilaw].harmhapke.g.a[l][m].state == 'f') {
						fpntr[++p] = &dmod->photo.optical[ilaw].harmhapke.g.a[l][m].val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
					if (m > 0 && dmod->photo.optical[ilaw].harmhapke.g.b[l][m].state == 'f') {
						fpntr[++p] = &dmod->photo.optical[ilaw].harmhapke.g.b[l][m].val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
				}
			for (l=0; l<=dmod->photo.optical[ilaw].harmhapke.theta.nhar; l++)
				for (m=0; m<=l; m++) {
					if (dmod->photo.optical[ilaw].harmhapke.theta.a[l][m].state == 'f') {
						fpntr[++p] = &dmod->photo.optical[ilaw].harmhapke.theta.a[l][m].val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
					if (m > 0 && dmod->photo.optical[ilaw].harmhapke.theta.b[l][m].state == 'f') {
						fpntr[++p] = &dmod->photo.optical[ilaw].harmhapke.theta.b[l][m].val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
				}
			break;
		case INHOHAPKE:
			if (dmod->photo.optical[ilaw].inhohapke.global.w.state == 'f') {
				fpntr[++p] = &dmod->photo.optical[ilaw].inhohapke.global.w.val;
				fparstep[p] = dpar->photo_step;
				fpartol[p] = dpar->photo_tol;
				fparabstol[p] = dpar->photo_abstol;
				fpartype[p] = PHOTOPAR;
			}
			if (dmod->photo.optical[ilaw].inhohapke.global.h.state == 'f') {
				fpntr[++p] = &dmod->photo.optical[ilaw].inhohapke.global.h.val;
				fparstep[p] = dpar->photo_step;
				fpartol[p] = dpar->photo_tol;
				fparabstol[p] = dpar->photo_abstol;
				fpartype[p] = PHOTOPAR;
			}
			if (dmod->photo.optical[ilaw].inhohapke.global.B0.state == 'f') {
				fpntr[++p] = &dmod->photo.optical[ilaw].inhohapke.global.B0.val;
				fparstep[p] = dpar->photo_step;
				fpartol[p] = dpar->photo_tol;
				fparabstol[p] = dpar->photo_abstol;
				fpartype[p] = PHOTOPAR;
			}
			if (dmod->photo.optical[ilaw].inhohapke.global.g.state == 'f') {
				fpntr[++p] = &dmod->photo.optical[ilaw].inhohapke.global.g.val;
				fparstep[p] = dpar->photo_step;
				fpartol[p] = dpar->photo_tol;
				fparabstol[p] = dpar->photo_abstol;
				fpartype[p] = PHOTOPAR;
			}
			if (dmod->photo.optical[ilaw].inhohapke.global.theta.state == 'f') {
				fpntr[++p] = &dmod->photo.optical[ilaw].inhohapke.global.theta.val;
				fparstep[p] = dpar->photo_step;
				fpartol[p] = dpar->photo_tol;
				fparabstol[p] = dpar->photo_abstol;
				fpartype[p] = PHOTOPAR;
			}
			for (c=0; c<dmod->shape.ncomp; c++)
				for (f=0; f<dmod->shape.comp[c].real.nf; f++) {
					if (dmod->photo.optical[ilaw].inhohapke.local[c][f].w.state == 'f') {
						fpntr[++p] = &dmod->photo.optical[ilaw].inhohapke.local[c][f].w.val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
					if (dmod->photo.optical[ilaw].inhohapke.local[c][f].h.state == 'f') {
						fpntr[++p] = &dmod->photo.optical[ilaw].inhohapke.local[c][f].h.val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
					if (dmod->photo.optical[ilaw].inhohapke.local[c][f].B0.state == 'f') {
						fpntr[++p] = &dmod->photo.optical[ilaw].inhohapke.local[c][f].B0.val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
					if (dmod->photo.optical[ilaw].inhohapke.local[c][f].g.state == 'f') {
						fpntr[++p] = &dmod->photo.optical[ilaw].inhohapke.local[c][f].g.val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
					if (dmod->photo.optical[ilaw].inhohapke.local[c][f].theta.state == 'f') {
						fpntr[++p] = &dmod->photo.optical[ilaw].inhohapke.local[c][f].theta.val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
				}
			break;
		case KAASALAINEN:
			if (dmod->photo.optical[ilaw].kaas.R.state == 'f') {
				fpntr[++p] = &dmod->photo.optical[ilaw].kaas.R.val;
				fparstep[p] = dpar->photo_step;
				fpartol[p] = dpar->photo_tol;
				fparabstol[p] = dpar->photo_abstol;
				fpartype[p] = PHOTOPAR;
			}
			if (dmod->photo.optical[ilaw].kaas.wt.state == 'f') {
				fpntr[++p] = &dmod->photo.optical[ilaw].kaas.wt.val;
				fparstep[p] = dpar->photo_step;
				fpartol[p] = dpar->photo_tol;
				fparabstol[p] = dpar->photo_abstol;
				fpartype[p] = PHOTOPAR;
			}
			if (dmod->photo.optical[ilaw].kaas.A0.state == 'f') {
				fpntr[++p] = &dmod->photo.optical[ilaw].kaas.A0.val;
				fparstep[p] = dpar->photo_step;
				fpartol[p] = dpar->photo_tol;
				fparabstol[p] = dpar->photo_abstol;
				fpartype[p] = PHOTOPAR;
			}
			if (dmod->photo.optical[ilaw].kaas.D.state == 'f') {
				fpntr[++p] = &dmod->photo.optical[ilaw].kaas.D.val;
				fparstep[p] = dpar->photo_step;
				fpartol[p] = dpar->photo_tol;
				fparabstol[p] = dpar->photo_abstol;
				fpartype[p] = PHOTOPAR;
			}
			if (dmod->photo.optical[ilaw].kaas.k.state == 'f') {
				fpntr[++p] = &dmod->photo.optical[ilaw].kaas.k.val;
				fparstep[p] = dpar->photo_step;
				fpartol[p] = dpar->photo_tol;
				fparabstol[p] = dpar->photo_abstol;
				fpartype[p] = PHOTOPAR;
			}
			break;
		case HARMKAAS:
			for (l=0; l<=dmod->photo.optical[ilaw].harmkaas.R.nhar; l++)
				for (m=0; m<=l; m++) {
					if (dmod->photo.optical[ilaw].harmkaas.R.a[l][m].state == 'f') {
						fpntr[++p] = &dmod->photo.optical[ilaw].harmkaas.R.a[l][m].val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
					if (m > 0 && dmod->photo.optical[ilaw].harmkaas.R.b[l][m].state == 'f') {
						fpntr[++p] = &dmod->photo.optical[ilaw].harmkaas.R.b[l][m].val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
				}
			for (l=0; l<=dmod->photo.optical[ilaw].harmkaas.wt.nhar; l++)
				for (m=0; m<=l; m++) {
					if (dmod->photo.optical[ilaw].harmkaas.wt.a[l][m].state == 'f') {
						fpntr[++p] = &dmod->photo.optical[ilaw].harmkaas.wt.a[l][m].val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
					if (m > 0 && dmod->photo.optical[ilaw].harmkaas.wt.b[l][m].state == 'f') {
						fpntr[++p] = &dmod->photo.optical[ilaw].harmkaas.wt.b[l][m].val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
				}
			for (l=0; l<=dmod->photo.optical[ilaw].harmkaas.A0.nhar; l++)
				for (m=0; m<=l; m++) {
					if (dmod->photo.optical[ilaw].harmkaas.A0.a[l][m].state == 'f') {
						fpntr[++p] = &dmod->photo.optical[ilaw].harmkaas.A0.a[l][m].val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
					if (m > 0 && dmod->photo.optical[ilaw].harmkaas.A0.b[l][m].state == 'f') {
						fpntr[++p] = &dmod->photo.optical[ilaw].harmkaas.A0.b[l][m].val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
				}
			for (l=0; l<=dmod->photo.optical[ilaw].harmkaas.D.nhar; l++)
				for (m=0; m<=l; m++) {
					if (dmod->photo.optical[ilaw].harmkaas.D.a[l][m].state == 'f') {
						fpntr[++p] = &dmod->photo.optical[ilaw].harmkaas.D.a[l][m].val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
					if (m > 0 && dmod->photo.optical[ilaw].harmkaas.D.b[l][m].state == 'f') {
						fpntr[++p] = &dmod->photo.optical[ilaw].harmkaas.D.b[l][m].val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
				}
			for (l=0; l<=dmod->photo.optical[ilaw].harmkaas.k.nhar; l++)
				for (m=0; m<=l; m++) {
					if (dmod->photo.optical[ilaw].harmkaas.k.a[l][m].state == 'f') {
						fpntr[++p] = &dmod->photo.optical[ilaw].harmkaas.k.a[l][m].val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
					if (m > 0 && dmod->photo.optical[ilaw].harmkaas.k.b[l][m].state == 'f') {
						fpntr[++p] = &dmod->photo.optical[ilaw].harmkaas.k.b[l][m].val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
				}
			break;
		case INHOKAAS:
			if (dmod->photo.optical[ilaw].inhokaas.global.R.state == 'f') {
				fpntr[++p] = &dmod->photo.optical[ilaw].inhokaas.global.R.val;
				fparstep[p] = dpar->photo_step;
				fpartol[p] = dpar->photo_tol;
				fparabstol[p] = dpar->photo_abstol;
				fpartype[p] = PHOTOPAR;
			}
			if (dmod->photo.optical[ilaw].inhokaas.global.wt.state == 'f') {
				fpntr[++p] = &dmod->photo.optical[ilaw].inhokaas.global.wt.val;
				fparstep[p] = dpar->photo_step;
				fpartol[p] = dpar->photo_tol;
				fparabstol[p] = dpar->photo_abstol;
				fpartype[p] = PHOTOPAR;
			}
			if (dmod->photo.optical[ilaw].inhokaas.global.A0.state == 'f') {
				fpntr[++p] = &dmod->photo.optical[ilaw].inhokaas.global.A0.val;
				fparstep[p] = dpar->photo_step;
				fpartol[p] = dpar->photo_tol;
				fparabstol[p] = dpar->photo_abstol;
				fpartype[p] = PHOTOPAR;
			}
			if (dmod->photo.optical[ilaw].inhokaas.global.D.state == 'f') {
				fpntr[++p] = &dmod->photo.optical[ilaw].inhokaas.global.D.val;
				fparstep[p] = dpar->photo_step;
				fpartol[p] = dpar->photo_tol;
				fparabstol[p] = dpar->photo_abstol;
				fpartype[p] = PHOTOPAR;
			}
			if (dmod->photo.optical[ilaw].inhokaas.global.k.state == 'f') {
				fpntr[++p] = &dmod->photo.optical[ilaw].inhokaas.global.k.val;
				fparstep[p] = dpar->photo_step;
				fpartol[p] = dpar->photo_tol;
				fparabstol[p] = dpar->photo_abstol;
				fpartype[p] = PHOTOPAR;
			}
			for (c=0; c<dmod->shape.ncomp; c++)
				for (f=0; f<dmod->shape.comp[c].real.nf; f++) {
					if (dmod->photo.optical[ilaw].inhokaas.local[c][f].R.state == 'f') {
						fpntr[++p] = &dmod->photo.optical[ilaw].inhokaas.local[c][f].R.val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
					if (dmod->photo.optical[ilaw].inhokaas.local[c][f].wt.state == 'f') {
						fpntr[++p] = &dmod->photo.optical[ilaw].inhokaas.local[c][f].wt.val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
					if (dmod->photo.optical[ilaw].inhokaas.local[c][f].A0.state == 'f') {
						fpntr[++p] = &dmod->photo.optical[ilaw].inhokaas.local[c][f].A0.val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
					if (dmod->photo.optical[ilaw].inhokaas.local[c][f].D.state == 'f') {
						fpntr[++p] = &dmod->photo.optical[ilaw].inhokaas.local[c][f].D.val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
					if (dmod->photo.optical[ilaw].inhokaas.local[c][f].k.state == 'f') {
						fpntr[++p] = &dmod->photo.optical[ilaw].inhokaas.local[c][f].k.val;
						fparstep[p] = dpar->photo_step;
						fpartol[p] = dpar->photo_tol;
						fparabstol[p] = dpar->photo_abstol;
						fpartype[p] = PHOTOPAR;
					}
				}
			break;
		default:
			printf("mkparlist_cuda.c: can't do that optical law yet\n");
		}
	}
}
__global__ void mpl_spin_krnl(struct par_t *dpar, struct mod_t *dmod,
		double **fpntr, double *fparstep, double *fpartol, double *fparabstol,
		int *fpartype) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		int i, k;
		for (i=0; i<=2; i++) {
			if (dmod->spin.angle[i].state == 'f') {
				fpntr[++p] = &dmod->spin.angle[i].val;
				fparstep[p] = dpar->angle_step;
				fpartol[p] = dpar->angle_tol;
				fparabstol[p] = dpar->angle_abstol;
				fpartype[p] = SPINPAR;
			}
		}
		for (i=0; i<=2; i++) {
			if (dmod->spin.omega[i].state == 'f') {
				fpntr[++p] = &dmod->spin.omega[i].val;
				fparstep[p] = dpar->spin_step;
				fpartol[p] = dpar->spin_tol;
				fparabstol[p] = dpar->spin_abstol;
				fpartype[p] = SPINPAR;
			}
		}
		for (i=0; i<=2; i++) {
			if (dmod->spin.inertia[i].state == 'f') {
				fpntr[++p] = &dmod->spin.inertia[i].val;
				fparstep[p] = dpar->inertia_step;
				fpartol[p] = dpar->inertia_tol;
				fparabstol[p] = dpar->inertia_abstol;
				fpartype[p] = SPINPAR;
			}
		}
		for (i=0; i<=2; i++) {
			if (dmod->spin.omegadot[i].state == 'f') {
				fpntr[++p] = &dmod->spin.omegadot[i].val;
				fparstep[p] = dpar->spindot_step;
				fpartol[p] = dpar->spindot_tol;
				fparabstol[p] = dpar->spindot_abstol;
				fpartype[p] = SPINPAR;
			}
		}
		for (k=0; k<dmod->spin.n_impulse; k++)
			for (i=0; i<=2; i++) {
				if (dmod->spin.impulse[k][i].state == 'f') {
					fpntr[++p] = &dmod->spin.impulse[k][i].val;
					fparstep[p] = dpar->spin_step;
					fpartol[p] = dpar->spin_tol;
					fparabstol[p] = dpar->spin_abstol;
					fpartype[p] = SPINPAR;
				}
			}

		if (dmod->spin.lib_amp.state == 'f') {
			fpntr[++p] = &dmod->spin.lib_amp.val;
			fparstep[p] = dpar->lib_amp_step;
			fpartol[p] = dpar->lib_amp_tol;
			fparabstol[p] = dpar->lib_amp_abstol;
			fpartype[p] = SPINPAR;
		}

		if (dmod->spin.lib_freq.state == 'f') {
			fpntr[++p] = &dmod->spin.lib_freq.val;
			fparstep[p] = dpar->lib_freq_step;
			fpartol[p] = dpar->lib_freq_tol;
			fparabstol[p] = dpar->lib_freq_abstol;
			fpartype[p] = SPINPAR;
		}

		if (dmod->spin.lib_phase.state == 'f') {
			fpntr[++p] = &dmod->spin.lib_phase.val;
			fparstep[p] = dpar->lib_phase_step;
			fpartol[p] = dpar->lib_phase_tol;
			fparabstol[p] = dpar->lib_phase_abstol;
			fpartype[p] = SPINPAR;
		}
	}

}
__global__ void mpl_dat_krnl(struct par_t *dpar, struct dat_t *ddat,
		double **fpntr, double *fparstep, double *fpartol, double *fparabstol,
		int *fpartype, int nsets) {
	/* nsets-threaded kernel */
	int s = blockIdx.x * blockDim.x + threadIdx.x;
	int i, j;

	if (s < nsets) {
		for (i=0; i<=2; i++) {
			if (ddat->set[s].angleoff[i].state == 'f') {
				fpntr[++p] = &ddat->set[s].angleoff[i].val;
				fparstep[p] = dpar->angle_step;
				fpartol[p] = dpar->angle_tol;
				fparabstol[p] = dpar->angle_abstol;
				fpartype[p] = SPINPAR;
			}
		}
		for (i=0; i<=2; i++) {
			if (ddat->set[s].omegaoff[i].state == 'f') {
				fpntr[++p] = &ddat->set[s].omegaoff[i].val;
				fparstep[p] = dpar->spin_step;
				fpartol[p] = dpar->spin_tol;
				fparabstol[p] = dpar->spin_abstol;
				fpartype[p] = SPINPAR;
			}
		}

		switch (ddat->set[s].type) {
		case DELAY:
			for (i=0; i<=ddat->set[s].desc.deldop.delcor.n; i++) {
				if (ddat->set[s].desc.deldop.delcor.a[i].state == 'f') {
					j = (i < MAXDELCORPAR) ? i : 0;
					fpntr[++p] = &ddat->set[s].desc.deldop.delcor.a[i].val;
					fparstep[p] = dpar->delcor_step[j];
					fpartol[p] = dpar->delcor_tol;
					fparabstol[p] = dpar->delcor_abstol[j];
					fpartype[p] = DELCORPAR;
				}
			}
			if (ddat->set[s].desc.deldop.dopscale.state == 'f') {
				fpntr[++p] = &ddat->set[s].desc.deldop.dopscale.val;
				fparstep[p] = dpar->ratio_step;
				fpartol[p] = dpar->ratio_tol;
				fparabstol[p] = dpar->ratio_abstol;
				fpartype[p] = DOPSCALEPAR;
			}
			break;
		case DOPPLER:
			for (i=0; i<=ddat->set[s].desc.doppler.delcor.n; i++) {
				if (ddat->set[s].desc.doppler.delcor.a[i].state == 'f') {
					j = (i < MAXDELCORPAR) ? i : 0;
					fpntr[++p] = &ddat->set[s].desc.doppler.delcor.a[i].val;
					fparstep[p] = dpar->delcor_step[j];
					fpartol[p] = dpar->delcor_tol;
					fparabstol[p] = dpar->delcor_abstol[j];
					fpartype[p] = DELCORPAR;
				}
			}
			if (ddat->set[s].desc.doppler.dopscale.state == 'f') {
				fpntr[++p] = &ddat->set[s].desc.doppler.dopscale.val;
				fparstep[p] = dpar->ratio_step;
				fpartol[p] = dpar->ratio_tol;
				fparabstol[p] = dpar->ratio_abstol;
				fpartype[p] = DOPSCALEPAR;
			}
			break;
		case POS:
			for (i=0; i<ddat->set[s].desc.poset.nframes; i++) {
				for (j=0; j<=1; j++) {
					if (ddat->set[s].desc.poset.frame[i].off[j].state == 'f') {
						fpntr[++p] = &ddat->set[s].desc.poset.frame[i].off[j].val;
						fparstep[p] = dpar->xyoff_step;
						fpartol[p] = dpar->xyoff_tol;
						fparabstol[p] = dpar->xyoff_abstol;
						fpartype[p] = XYOFFPAR;
					}
				}
			}
			break;
		case LGHTCRV:
			break;
		default:
			printf("mkparlist_cuda.cu: can't handle that type of data yet\n");
		}

	}
}

__host__ void mkparlist_gpu(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, double *fparstep, double *fpartol,
		double *fparabstol, int *fpartype, double **fpntr,
		int nfpar, int nsets)
{
	dim3 BLK,THD;
	THD.x = maxThreadsPerBlock;

	/* Shape parameters - single component only */
	//for (i=0; i<dmod->shape.ncomp; i++) { /* read each component */
	/* Launch first parameter kernel */
	mpl_comp_krnl<<<1,1>>>(dpar, dmod, fpntr, fparstep, fpartol,
			fparabstol,	fpartype);
	checkErrorAfterKernelLaunch("mpl_comp_krnl (mkparlist_cuda.cu)");

	/* Photometric parameters - only one radlaw at a time for now */
	//for (ilaw=0; ilaw<dmod->photo.nradlaws; ilaw++) {
	/* Launch photometric kernel */
	mpl_rad_krnl<<<1,1>>>(dpar, dmod, fpntr, fparstep, fpartol,
			fparabstol,	fpartype);
	checkErrorAfterKernelLaunch("mpl_rad_krnl (mkparlist_cuda.cu)");

	/* Photometric parameters - only one optlaw at a time */
	//for (ilaw=0; ilaw<dmod->photo.noptlaws; ilaw++) {
	mpl_photo_krnl<<<1,1>>>(dpar, dmod, fpntr, fparstep, fpartol,
			fparabstol,	fpartype);
	checkErrorAfterKernelLaunch("mpl_photo_krnl (mkparlist_cuda.cu)");

	/* Spin parameters  */
	mpl_spin_krnl<<<1,1>>>(dpar, dmod, fpntr, fparstep, fpartol,
			fparabstol,	fpartype);
	checkErrorAfterKernelLaunch("mpl_spin_krnl (mkparlist_cuda.cu)");

	/* Data parameters (i.e., those in the obs file, other than the "calfact"
	 * parameters which are computed analytically)
	 * Launching nsets threads here */
	BLK.x = floor((THD.x - 1 + nsets)/THD.x);
	mpl_dat_krnl<<<1,1>>>(dpar, ddat, fpntr, fparstep, fpartol,
			fparabstol,	fpartype, nsets);
	checkErrorAfterKernelLaunch("mpl_dat_krnl (mkparlist_cuda.cu)");

}
