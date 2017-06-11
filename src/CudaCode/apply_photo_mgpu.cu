/*****************************************************************************************
                                                                            apply_photo.c

For each plane-of-sky pixel, compute the model's scattered optical power per unit
projected (POS) area per unit solid angle per unit incident flux, and then sum these
values over the entire POS.  (The POS pixel area is multiplied in elsewhere.)

The expressions given here differ from the bidirectional reflectance functions defined by,
say, Hapke 1993: bidirectional reflectance includes an extra factor of
cos(scattering angle), since it is defined per unit surface area rather than per unit
projected area.

Modified 2014 February 12 by CM:
    Implement multiple optical scatering laws

Modified 2011 September 2 by CM:
    Add the "harmlambert" and "inholambert" optical scattering laws

Modified 2007 August 4 by CM:
    Add body parameter for use with the "orbit" action: it denotes which
        orbiting body's optical power contributions are being computed
        on this call to the routine
    Don't zero out blank-sky and shadowed POS pixels in the sky rendering
        (the pos->b matrix): do it instead in the calling routine by
        having it call the posclr routine.  This way apply_photo can be
        called twice for the "orbit" action, once for each orbiting body.
    Add comp matrix for POS frames

Modified 2006 October 1 by CM:
    Add "intensityfactor" parameter: account for POS pixel area,
        1 AU Sun-target distance, and solar apparent magnitude here
        rather than after calling the routine

Modified 2006 September 1 by CM and MCN:
    For inhomogeneous laws, add check that facet number pos->f[i][j]
        is nonnegative

Modified 2005 September 7 by CM:
    Implement the "harmlommel" "harmhapke" and "harmkaas" optical
        scattering laws

Modified 2005 August 8 by CM:
    Implement the "inhokaas" optical scattering law
    Add some (cosi > 0) checks
    Move "sum == 0" check to the end

Modified 2005 July 4 by CM:
    Changed structure name for the INHOLOMMEL optical scattering law

Modified 2005 March 1 by CM:
    Add NOLAW case

Modified 2005 January 25 by CM:
    Eliminate unused variables

Modified 2004 April 29 by CM:
    Modify Kaasalainen scattering law to use "wt" as the relative
        weighting factor (0 = pure Lommel-Seeliger, 1 = pure Lambert)
        rather than "c" (which ranged from 0 to infinity)

Modified 2004 March 25 by CM:
    hapke routine now takes phase rather than cos(phase) as argument

Modified 2004 February 29 by CM:
    Added comments
    Added Kaasalainen "Lommel-Seeliger + Lambert" scattering law
    Eliminated "type" argument, since this routine was only being
       used to handle optical scattering.  (Radar scattering is
       instead handled by the "radlaw" routine.)
    Added "phase" argument (solar phase angle) so that we can compute
       the phase just once per calculated lightcurve point (in read_dat)
       rather than computing it every time we call apply_photo
*****************************************************************************************/
extern "C" {
#include "../shape/head.h"
}

#define TINY 1.0e-40

__device__ int mgpu_ilaw;
//__device__ double phasefunc, scale_lommsee, scale_lambert, intensityfactor, phase;
//__device__ float sum, phasefuncf, scale_lommseef, scale_lambertf;

__global__ void ap_init_mgpu_krnl(
		struct dat_t *ddat,
		struct mod_t *dmod,
		struct pos_t **pos,
		int set, int size, unsigned char *type, float *dsum,
		double *intensity_factor,
		double *phase_d,
		int oddflg) {

	/* nfrm_half0/nfrm_half1-threaded kernel for multi-GPU operation */
	int hf = blockIdx.x * blockDim.x + threadIdx.x + 1;
	int f = 2 * hf + oddflg + 1;
	if (hf < size) {
		if (hf==0 && oddflg==0) {
			mgpu_ilaw = ddat->set[set].desc.lghtcrv.ioptlaw;
			type[0] = dmod->photo.opttype[mgpu_ilaw];
		}
		dsum[hf] = 0.0;
		intensity_factor[hf] = (pos[hf]->km_per_pixel/AU) * (pos[hf]->km_per_pixel/AU);
		phase_d[hf] = ddat->set[set].desc.lghtcrv.solar_phase[f];
	}
}
__global__ void ap_init_mgpu_f_krnl(
		struct dat_t *ddat,
		struct mod_t *dmod,
		struct pos_t **pos,
		int set, int size, unsigned char *type, float *dsum,
		float *intensity_factor,
		float *phase_f,
		int oddflg) {

	/* nfrm_half0/nfrm_half1-threaded kernel for multi-GPU operation */
	int hf = blockIdx.x * blockDim.x + threadIdx.x;
	int f = 2 * hf + oddflg + 1;
	float temp;
	if (hf < size) {
		if (hf==0 && oddflg==0) {
			mgpu_ilaw = ddat->set[set].desc.lghtcrv.ioptlaw;
			type[0] = dmod->photo.opttype[mgpu_ilaw];
		}

		dsum[hf] = 0.0;
		temp = __double2float_rn(pos[hf]->km_per_pixel/AU);
		intensity_factor[hf] = temp*temp;
		phase_f[hf] = __double2float_rn(ddat->set[set].desc.lghtcrv.solar_phase[f]);
	}
}

__global__ void ap_lambertlaw_mgpu_f_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		float *intensity_factor,
		int4 *xylim,
		int nThreads,
		int body,
		int2 span,
		int hf) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[hf].w;
	int j = offset / span.x + xylim[hf].y;
	float scale;

	if (offset < nThreads) {
		scale = __double2float_rn(dmod->photo.optical[mgpu_ilaw].R.R.val/PIE);

		if (pos[hf]->cose_s[offset] > 0.0 && pos[hf]->cosi_s[offset] > 0.0
				&& pos[hf]->body[i][j] == body) {
			pos[hf]->b_s[offset] = intensity_factor[hf] * scale * pos[hf]->cosi_s[offset];
		}
	}
}

__global__ void ap_harmlambert_mgpu_f_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		float *intensity_factor,
		int hf) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[hf].w;
	int j = offset / span.x + xylim[hf].y;
	int c, f;
	float scale;

	if (offset < nThreads) {

		if (pos[hf]->cose_s[offset] > 0.0 && pos[hf]->cosi_s[offset] > 0.0
		 && pos[hf]->body[i][j] == body && pos[hf]->f[i][j] >= 0) {
			c = pos[hf]->comp[i][j];
			f = pos[hf]->f[i][j];
			scale = __double2float_rn(dmod->photo.optical[mgpu_ilaw].harmR.local[c][f].R.val)/PIE;
			pos[hf]->b_s[offset] = intensity_factor[hf] * scale *
					pos[hf]->cosi_s[offset];
		}
	}
}

__global__ void ap_inholambert_mgpu_f_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		float *intensity_factor,
		int hf) {
	/* Multi-threaded kernel for dual-GPU operation */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[hf].w;
	int j = offset / span.x + xylim[hf].y;
	int c, f;
	float scale;

	if (offset < nThreads) {

		if (pos[hf]->cose_s[offset] > 0.0 && pos[hf]->cosi_s[offset] > 0.0
		 && pos[hf]->body[i][j] == body && pos[hf]->f[i][j] >= 0) {
			c = pos[hf]->comp[i][j];
			f = pos[hf]->f[i][j];
			scale = __double2float_rn(dmod->photo.optical[mgpu_ilaw].inhoR.local[c][f].R.val)/PIE;
			pos[hf]->b_s[offset] = intensity_factor[hf] * scale * pos[hf]->cosi_s[offset];
		}
	}
}

__global__ void ap_lommel_mgpu_f_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		float *intensity_factor,
		int hf) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[hf].w;
	int j = offset / span.x + xylim[hf].y;
	float scale;

	if (offset < nThreads) {
		scale = __double2float_rn(dmod->photo.optical[mgpu_ilaw].R.R.val)/(4*PIE);
		if (pos[hf]->cose_s[offset] > 0.0 && pos[hf]->cosi_s[offset] > 0.0
				&& pos[hf]->body[i][j] == body) {
			pos[hf]->b_s[offset] = intensity_factor[hf] * scale * pos[hf]->cosi_s[offset]
			  / (pos[hf]->cosi_s[offset] + pos[hf]->cose_s[offset]);
		}
	}
}

__global__ void ap_harmlommel_mgpu_f_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		float *intensity_factor,
		int hf) {

	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[hf].w;
	int j = offset / span.x + xylim[hf].y;
	int c, f;
	float scale;

	if (offset < nThreads) {
		if (pos[hf]->cose_s[offset] > 0.0 && pos[hf]->cosi_s[offset] > 0.0
		 && pos[hf]->body[i][j] == body && pos[hf]->f[i][j] >= 0) {
			c = pos[hf]->comp[i][j];
			f = pos[hf]->f[i][j];
			scale = __double2float_rn(dmod->photo.optical[mgpu_ilaw].harmR.local[c][f].R.val)/(4*PIE);
			pos[hf]->b_s[offset] = intensity_factor[hf] * scale * pos[hf]->cosi_s[offset]
			   / (pos[hf]->cosi_s[offset] + pos[hf]->cose_s[offset]);
		}
	}
}

__global__ void ap_inholommel_mgpu_f_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		float *intensity_factor,
		int hf) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[hf].w;
	int j = offset / span.x + xylim[hf].y;
	int c, f;
	float scale;

	if (offset < nThreads) {
		if (pos[hf]->cose_s[offset] > 0.0 && pos[hf]->cosi_s[offset] > 0.0
		 && pos[hf]->body[i][j] == body && pos[hf]->f[i][j] >= 0) {
			c = pos[hf]->comp[i][j];
			f = pos[hf]->f[i][j];
			scale = __double2float_rn(dmod->photo.optical[mgpu_ilaw].inhoR.local[c][f].R.val)/(4*PIE);
			pos[hf]->b_s[offset] = intensity_factor[hf] * scale * pos[hf]->cosi_s[offset]
			   / (pos[hf]->cosi_s[offset] + pos[hf]->cose_s[offset]);
		}
	}
}

__global__ void ap_geometrical_mgpu_f_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		float *intensity_factor,
		int hf) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[hf].w;
	int j = offset / span.x + xylim[hf].y;

	if (offset < nThreads) {
		if (pos[hf]->cose_s[offset] > 0.0 && pos[hf]->cosi_s[offset] > 0.0
		 && pos[hf]->body[i][j] == body) {
			pos[hf]->b_s[offset] = intensity_factor[hf] * dmod->photo.optical[mgpu_ilaw].R.R.val;
		}
	}
}

__global__ void ap_hapke_mgpu_f_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		float *intensity_factor,
		float *phase_f,
		int hf) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[hf].w;
	int j = offset / span.x + xylim[hf].y;

	if (offset < nThreads) {
		if (pos[hf]->cose_s[offset] > 0.0 && pos[hf]->cosi_s[offset] > 0.0
		 && pos[hf]->body[i][j] == body) {
			pos[hf]->b_s[offset] = intensity_factor[hf]
					* dev_hapke_f(pos[hf]->cosi_s[offset],
							pos[hf]->cose_s[offset],
							phase_f[hf],
							__double2float_rn(dmod->photo.optical[mgpu_ilaw].hapke.w.val),
							__double2float_rn(dmod->photo.optical[mgpu_ilaw].hapke.h.val),
							__double2float_rn(dmod->photo.optical[mgpu_ilaw].hapke.B0.val),
							__double2float_rn(dmod->photo.optical[mgpu_ilaw].hapke.g.val),
							__double2float_rn(dmod->photo.optical[mgpu_ilaw].hapke.theta.val));
		}
	}
}

__global__ void ap_harmhapke_mgpu_f_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		float *intensity_factor,
		float *phase_f,
		int hf) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[hf].w;
	int j = offset / span.x + xylim[hf].y;
	int c, f;

	if (offset < nThreads) {
		if (pos[hf]->cose_s[offset] > 0.0 && pos[hf]->cosi_s[offset] > 0.0
	     && pos[hf]->body[i][j] == body && pos[hf]->f[i][j] >= 0) {
			c = pos[hf]->comp[i][j];
			f = pos[hf]->f[i][j];
			pos[hf]->b_s[offset] = intensity_factor[hf]
		    *dev_hapke_f(pos[hf]->cosi_s[offset], pos[hf]->cose_s[offset],
		                   		phase_f[hf],
		__double2float_rn(dmod->photo.optical[mgpu_ilaw].harmhapke.local[c][f].w.val),
		__double2float_rn(dmod->photo.optical[mgpu_ilaw].harmhapke.local[c][f].h.val),
		__double2float_rn(dmod->photo.optical[mgpu_ilaw].harmhapke.local[c][f].B0.val),
		__double2float_rn(dmod->photo.optical[mgpu_ilaw].harmhapke.local[c][f].g.val),
		__double2float_rn(dmod->photo.optical[mgpu_ilaw].harmhapke.local[c][f].theta.val));
		}
	}
}

__global__ void ap_inhohapke_mgpu_f_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		float *intensity_factor,
		float *phase_f,
		int hf) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[hf].w;
	int j = offset / span.x + xylim[hf].y;
	int c, f;

	if (offset < nThreads) {
		if (pos[hf]->cose_s[offset] > 0.0 && pos[hf]->cosi_s[offset] > 0.0
		 && pos[hf]->body[i][j] == body && pos[hf]->f[i][j] >= 0) {
			c = pos[hf]->comp[i][j];
			f = pos[hf]->f[i][j];
			pos[hf]->b_s[offset] = intensity_factor[hf]
					* dev_hapke_f(pos[hf]->cosi_s[offset], pos[hf]->cose_s[offset],
							phase_f[hf],
							__double2float_rn(dmod->photo.optical[mgpu_ilaw].inhohapke.local[c][f].w.val),
							__double2float_rn(dmod->photo.optical[mgpu_ilaw].inhohapke.local[c][f].h.val),
							__double2float_rn(dmod->photo.optical[mgpu_ilaw].inhohapke.local[c][f].B0.val),
							__double2float_rn(dmod->photo.optical[mgpu_ilaw].inhohapke.local[c][f].g.val),
							__double2float_rn(dmod->photo.optical[mgpu_ilaw].inhohapke.local[c][f].theta.val));
		}
	}
}

__global__ void ap_kaas_init_mgpu_f_krnl(
		struct mod_t *dmod,
		float *phasefuncd,
		float *phase_d,
		float *scale_lommsee,
		float *scale_lambert,
		int nfrm_alloc) {
	/* nframes-threaded kernel */
	int hf = blockIdx.x * blockDim.x + threadIdx.x;
	int mgpu_ilaw = 0;

	if (hf < nfrm_alloc) {
		phasefuncd[hf] = dmod->photo.optical[mgpu_ilaw].kaas.A0.val
				* exp( -phase_d[hf] / dmod->photo.optical[mgpu_ilaw].kaas.D.val)
		+ dmod->photo.optical[mgpu_ilaw].kaas.k.val * phase_d[hf] + 1;

		scale_lommsee[hf] = (1 - dmod->photo.optical[mgpu_ilaw].kaas.wt.val)
				* phasefuncd[hf] * dmod->photo.optical[mgpu_ilaw].kaas.R.val/(4*PIE);
		scale_lambert[hf] = dmod->photo.optical[mgpu_ilaw].kaas.wt.val
				* phasefuncd[hf] * dmod->photo.optical[mgpu_ilaw].kaas.R.val/PIE;
	}
}

__global__ void ap_kaas_mgpu_f_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		float *intensity_factor,
		float *phase_d,
		float *phasefuncd,
		float *scale_lommsee,
		float *scale_lambert,
		int hf) {

	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[hf].w;
	int j = offset / span.x + xylim[hf].y;
	int n = pos[hf]->n;
	int pos_spn = 2*n+1;
	int pxa = (j+n)*pos_spn + (i+n);

	if (offset < nThreads) {
		if (pos[hf]->cose_s[pxa] > 0.0 && pos[hf]->cosi_s[pxa] > 0.0
				&& pos[hf]->body[i][j] == body) {
			pos[hf]->b_s[pxa] = intensity_factor[hf] * pos[hf]->cosi_s[pxa]
			     *(scale_lommsee[hf] / (pos[hf]->cosi_s[pxa] + pos[hf]->cose_s[pxa])
			    + scale_lambert[hf]);
		}
	}
}

__global__ void ap_harmkaas_mgpu_f_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		float *intensity_factor,
		float *phase_f,
		int hf) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[hf].w;
	int j = offset / span.x + xylim[hf].y;
	int c, f;
	float phasefuncf, scale_lommseef, scale_lambertf;

	if (offset < nThreads) {
		if (pos[hf]->cose_s[offset] > 0.0 && pos[hf]->cosi_s[offset] > 0.0
		 && pos[hf]->body[i][j] == body && pos[hf]->f[i][j] >= 0) {
			c = pos[hf]->comp[i][j];
			f = pos[hf]->f[i][j];
			phasefuncf = __double2float_rn(dmod->photo.optical[mgpu_ilaw].harmkaas.local[c][f].A0.val)
			* exp( -phase_f[hf] / __double2float_rn(dmod->photo.optical[mgpu_ilaw].harmkaas.local[c][f].D.val))
			+ __double2float_rn(dmod->photo.optical[mgpu_ilaw].harmkaas.local[c][f].k.val) * phase_f[hf] + 1;

			scale_lommseef = (1 - __double2float_rn(dmod->photo.optical[mgpu_ilaw].harmkaas.local[c][f].wt.val))
		    * phasefuncf * __double2float_rn(dmod->photo.optical[mgpu_ilaw].harmkaas.local[c][f].R.val) / (4*PIE);
			scale_lambertf = __double2float_rn(dmod->photo.optical[mgpu_ilaw].harmkaas.local[c][f].wt.val)
			* phasefuncf * __double2float_rn(dmod->photo.optical[mgpu_ilaw].harmkaas.local[c][f].R.val) / PIE;
			 pos[hf]->b_s[offset] = intensity_factor[hf] * pos[hf]->cosi_s[offset] * (scale_lommseef /
					( pos[hf]->cosi_s[offset] +  pos[hf]->cose_s[offset]) + scale_lambertf);
		}
	}
}

__global__ void ap_inhokaas_mgpu_f_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		float *intensity_factor,
		float *phase_f,
		int hf) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[hf].w;
	int j = offset / span.x + xylim[hf].y;
	int c, f;
	float phasefuncf, scale_lommseef, scale_lambertf;

	if (offset < nThreads) {
		if (pos[hf]->cose_s[offset] > 0.0 && pos[hf]->cosi_s[offset] > 0.0
		 && pos[hf]->body[i][j] == body && pos[hf]->f[i][j] >= 0) {
			c = pos[hf]->comp[i][j];
			f = pos[hf]->f[i][j];
			phasefuncf = __double2float_rn(dmod->photo.optical[mgpu_ilaw].inhokaas.local[c][f].A0.val)
			* exp( -phase_f[hf] / __double2float_rn(dmod->photo.optical[mgpu_ilaw].inhokaas.local[c][f].D.val))
			+ __double2float_rn(dmod->photo.optical[mgpu_ilaw].inhokaas.local[c][f].k.val) * phase_f[hf] + 1;
			scale_lommseef = (1 - __double2float_rn(dmod->photo.optical[mgpu_ilaw].inhokaas.local[c][f].wt.val))
		    * phasefuncf * __double2float_rn(dmod->photo.optical[mgpu_ilaw].inhokaas.local[c][f].R.val) / (4*PIE);
			scale_lambertf = __double2float_rn(dmod->photo.optical[mgpu_ilaw].inhokaas.local[c][f].wt.val)
			* phasefuncf * __double2float_rn(dmod->photo.optical[mgpu_ilaw].inhokaas.local[c][f].R.val) / PIE;
			pos[hf]->b_s[offset] = intensity_factor[hf] * pos[hf]->cosi_s[offset] * (scale_lommseef /
					(pos[hf]->cosi_s[offset] + pos[hf]->cose_s[offset]) + scale_lambertf);
		}
	}
}


__host__ void apply_photo_mgpu_f(
		struct mod_t *dmod,
		struct dat_t *ddat,
		struct pos_t **pos0,
		struct pos_t **pos1,
		int4 *xylim0,
		int4 *xylim1,
		int2 *span,
		dim3 *BLKpx,
		int *nThreads,
		int body,
		int set,
		int nfrm_alloc,	/* This already includes the additional +1 for lghtcrv */
		int nfrm_half0,
		int nfrm_half1,
		int *nThreadspx,
		cudaStream_t *gpu0_stream,
		cudaStream_t *gpu1_stream)
{
	unsigned char *type, *htype;
	int f=0, hf=0;
	float *dsum0, *dsum1;
	double *hsum, *sum;
	float *intensity_factor0, *intensity_factor1, *phase_f0, *phase_f1,
		*phasefuncf0, *phasefuncf1, *scale_lommsee0, *scale_lommsee1,
		*scale_lambert0, *scale_lambert1;
	dim3 BLK, THD, BLK_half0, BLK_half1, THD64;

	/* Host allocations */
	htype = (unsigned char *) malloc(2*sizeof(unsigned char));
	hsum = (double *) malloc(nfrm_alloc*sizeof(double));

	/* Allocate GPU0 memory */
	gpuErrchk(cudaSetDevice(GPU0));
	gpuErrchk(cudaMalloc((void**)&type, sizeof(unsigned char) * 2));
	gpuErrchk(cudaMalloc((void**)&sum, sizeof(double)*nfrm_alloc));
	gpuErrchk(cudaMalloc((void**)&dsum0, sizeof(float)*nfrm_half0));
	gpuErrchk(cudaMalloc((void**)&intensity_factor0, sizeof(float)*nfrm_half0));
	gpuErrchk(cudaMalloc((void**)&phase_f0, sizeof(float)*nfrm_half0));
	gpuErrchk(cudaMalloc((void**)&phasefuncf0, sizeof(float)*nfrm_half0));
	gpuErrchk(cudaMalloc((void**)&scale_lommsee0, sizeof(float)*nfrm_half0));
	gpuErrchk(cudaMalloc((void**)&scale_lambert0, sizeof(float)*nfrm_half0));

	/* Allocate GPU1 memory */
	gpuErrchk(cudaSetDevice(GPU1));
	gpuErrchk(cudaMalloc((void**)&dsum1, sizeof(float)*nfrm_half1));
	gpuErrchk(cudaMalloc((void**)&intensity_factor1, sizeof(float)*nfrm_half1));
	gpuErrchk(cudaMalloc((void**)&phase_f1, sizeof(float)*nfrm_half1));
	gpuErrchk(cudaMalloc((void**)&phasefuncf1, sizeof(float)*nfrm_half1));
	gpuErrchk(cudaMalloc((void**)&scale_lommsee1, sizeof(float)*nfrm_half1));
	gpuErrchk(cudaMalloc((void**)&scale_lambert1, sizeof(float)*nfrm_half1));

	gpuErrchk(cudaSetDevice(GPU0));

	/* Assign pos addresses and get type, first for GPU0, then GPU1 */
	THD.x = maxThreadsPerBlock;
	THD64.x = 64;
	BLK.x = floor((THD.x - 1 + nfrm_alloc) / THD.x);
	BLK_half0.x = floor((THD64.x - 1 + nfrm_half0)/THD64.x);
	BLK_half1.x = floor((THD64.x - 1 + nfrm_half1)/THD64.x);
	ap_init_mgpu_f_krnl<<<BLK_half0,THD64,0,gpu0_stream[0]>>>(ddat, dmod, pos0,
			set, nfrm_half0, type, dsum0, intensity_factor0, phase_f0, 0);
	gpuErrchk(cudaSetDevice(GPU1));
	ap_init_mgpu_f_krnl<<<BLK_half0,THD64,0,gpu1_stream[0]>>>(ddat, dmod, pos1,
				set, nfrm_half1, type, dsum1, intensity_factor1, phase_f1, 0);
	checkErrorAfterKernelLaunch("ap_init_mgpu_f_krnl");
	gpuErrchk(cudaMemcpy(htype, type, sizeof(unsigned char) *2,
			cudaMemcpyDeviceToHost));
	gpuErrchk(cudaSetDevice(GPU0));

	switch (htype[0]) {
	case LAMBERTLAW:
		/* Launch Lambert Law kernel */
		hf = 0;
		for (f=1; f<nfrm_alloc; f++) {
			gpuErrchk(cudaSetDevice(GPU0));
			ap_lambertlaw_mgpu_f_krnl<<<BLKpx[f],THD,0,gpu0_stream[hf]>>>(dmod,pos0,
					intensity_factor0, xylim0, nThreads[f], body, span[f], hf);
			gpuErrchk(cudaSetDevice(GPU1));
			f++;	if (f>=nfrm_alloc)	break;
			ap_lambertlaw_mgpu_f_krnl<<<BLKpx[f],THD,0,gpu1_stream[hf]>>>(dmod,pos1,
					intensity_factor1, xylim1, nThreads[f], body, span[f], hf);
			hf++;
		}
		checkErrorAfterKernelLaunch("ap_lambertlaw_mgpu_f_krnl");
		break;
	case HARMLAMBERT:
		/* Launch the HarmLambert kernel */
		hf = 0;
		for (f=1; f<nfrm_alloc; f++) {
			gpuErrchk(cudaSetDevice(GPU0));
			ap_harmlambert_mgpu_f_krnl<<<BLKpx[f],THD,0,gpu0_stream[hf]>>>(dmod,
					pos0,nThreads[f],body,xylim0,span[f],intensity_factor0,hf);
			gpuErrchk(cudaSetDevice(GPU1));
			f++;	if (f>=nfrm_alloc)	break;
			ap_harmlambert_mgpu_f_krnl<<<BLKpx[f],THD,0,gpu1_stream[hf]>>>(dmod,
					pos1,nThreads[f],body,xylim1,span[f],intensity_factor1,hf);
			hf++;
		}
		checkErrorAfterKernelLaunch("ap_harmlambert_mgpu_f_krnl");
		break;
	case INHOLAMBERT:
		/* Launch the Inhomogeneous Lambert kernel */
		hf = 0;
		for (f=1; f<nfrm_alloc; f++) {
			gpuErrchk(cudaSetDevice(GPU0));
			ap_inholambert_mgpu_f_krnl<<<BLKpx[f],THD,0,gpu0_stream[hf]>>>(dmod,
					pos0,nThreads[f],body,xylim0,span[f],intensity_factor0,hf);
			gpuErrchk(cudaSetDevice(GPU1));
			f++;	if (f>=nfrm_alloc)	break;
			ap_inholambert_mgpu_f_krnl<<<BLKpx[f],THD,0,gpu1_stream[hf]>>>(dmod,
					pos1,nThreads[f],body,xylim1,span[f],intensity_factor0,hf);
			hf++;
		}
		checkErrorAfterKernelLaunch("ap_inholambert_mgpu_f_krnl");
		break;
	case LOMMEL:
		/* Launch the Lommel kernel */
		hf = 0;
		for (f=1; f<nfrm_alloc; f++) {
			gpuErrchk(cudaSetDevice(GPU0));
			ap_lommel_mgpu_f_krnl<<<BLKpx[f],THD,0,gpu0_stream[hf]>>>(dmod,pos0,
					nThreads[f], body, xylim0, span[f], intensity_factor0, hf);
			gpuErrchk(cudaSetDevice(GPU1));
			f++;	if (f>=nfrm_alloc)	break;
			ap_lommel_mgpu_f_krnl<<<BLKpx[f],THD,0,gpu1_stream[hf]>>>(dmod,pos1,
					nThreads[f], body, xylim1, span[f], intensity_factor1, hf);
			hf++;
		}
		checkErrorAfterKernelLaunch("ap_lommel_mgpu_f_krnl");
		break;
	case HARMLOMMEL:
		/* Launch the HarmLommel kernel */
		hf = 0;
		for (f=1; f<nfrm_alloc; f++) {
			gpuErrchk(cudaSetDevice(GPU0));
			ap_harmlommel_mgpu_f_krnl<<<BLKpx[f],THD,0,gpu0_stream[hf]>>>(dmod,
					pos0,nThreads[f],body,xylim0,span[f],intensity_factor0,hf);
			gpuErrchk(cudaSetDevice(GPU1));
			f++;	if (f>=nfrm_alloc)	break;
			ap_harmlommel_mgpu_f_krnl<<<BLKpx[f],THD,0,gpu1_stream[hf]>>>(dmod,
					pos1,nThreads[f],body,xylim1,span[f],intensity_factor1,hf);
			hf++;
		}
		checkErrorAfterKernelLaunch("ap_harmlommel_mgpu_f_krnl");
		break;
	case INHOLOMMEL:
		/* Launch the Inhomogeneous Lommel kernel */
		hf = 0;
		for (f=1; f<nfrm_alloc; f++) {
			gpuErrchk(cudaSetDevice(GPU0));
			ap_inholommel_mgpu_f_krnl<<<BLKpx[f],THD,0,gpu0_stream[hf]>>>(dmod,
					pos0,nThreads[f],body,xylim0,span[f],intensity_factor0,hf);
			gpuErrchk(cudaSetDevice(GPU1));
			f++;	if (f>=nfrm_alloc)	break;
			ap_inholommel_mgpu_f_krnl<<<BLKpx[f],THD,0,gpu1_stream[hf]>>>(dmod,
					pos1,nThreads[f],body,xylim1,span[f],intensity_factor1,hf);
			hf++;
		}
		checkErrorAfterKernelLaunch("ap_inholommel_mgpu_f_krnl");
		break;
	case GEOMETRICAL:
		/* Launch the Geometrical law kernel */
		hf = 0;
		for (f=1; f<nfrm_alloc; f++) {
			gpuErrchk(cudaSetDevice(GPU0));
			ap_geometrical_mgpu_f_krnl<<<BLKpx[f],THD,0,gpu0_stream[hf]>>>(dmod,
					pos0,nThreads[f],body,xylim0,span[f],intensity_factor0,hf);
			gpuErrchk(cudaSetDevice(GPU1));
			f++;	if (f>=nfrm_alloc)	break;
			ap_geometrical_mgpu_f_krnl<<<BLKpx[f],THD,0,gpu1_stream[hf]>>>(dmod,
					pos1,nThreads[f],body,xylim1,span[f],intensity_factor1,hf);
			hf++;
		}
		checkErrorAfterKernelLaunch("ap_geometrical_mgpu_f_krnl");
		break;
	case HAPKE:
		/* Launch the Hapke kernel */
		hf = 0;
		for (f=1; f<nfrm_alloc; f++) {
			gpuErrchk(cudaSetDevice(GPU0));
			ap_hapke_mgpu_f_krnl<<<BLKpx[f],THD,0,gpu0_stream[hf]>>>(dmod,pos0,
					nThreads[f], body, xylim0, span[f],
					intensity_factor0, phase_f0, hf);
			gpuErrchk(cudaSetDevice(GPU1));
			f++;	if (f>=nfrm_alloc)	break;
			ap_hapke_mgpu_f_krnl<<<BLKpx[f],THD,0,gpu1_stream[hf]>>>(dmod,pos1,
					nThreads[f], body, xylim1, span[f],
					intensity_factor1, phase_f1, hf);
			hf++;
		}
		checkErrorAfterKernelLaunch("ap_hapke_mgpu_f_krnl");
		break;
	case HARMHAPKE:
		/* Launch the HarmHapke kernel */
		hf = 0;
		for (f=1; f<nfrm_alloc; f++) {
			gpuErrchk(cudaSetDevice(GPU0));
			ap_harmhapke_mgpu_f_krnl<<<BLKpx[f],THD,0,gpu0_stream[hf]>>>(dmod,
					pos0, nThreads[f], body, xylim0, span[f], intensity_factor0,
					phase_f0, hf);
			gpuErrchk(cudaSetDevice(GPU1));
			f++;	if (f>=nfrm_alloc)	break;
			ap_harmhapke_mgpu_f_krnl<<<BLKpx[f],THD,0,gpu1_stream[hf]>>>(dmod,
					pos1, nThreads[f], body, xylim1, span[f], intensity_factor1,
					phase_f1, hf);
			hf++;
		}
		checkErrorAfterKernelLaunch("ap_harmhapke_mgpu_krnl");
		break;
	case INHOHAPKE:
		/* Launch the Inhomogeneous Hapke kernel */
		hf = 0;
		for (f=1; f<nfrm_alloc; f++) {
			gpuErrchk(cudaSetDevice(GPU0));
			ap_inhohapke_mgpu_f_krnl<<<BLKpx[f],THD,0,gpu0_stream[hf]>>>(dmod,
					pos0, nThreads[f], body, xylim0, span[f],
					intensity_factor0,	phase_f0, hf);
			gpuErrchk(cudaSetDevice(GPU1));
			f++;	if (f>=nfrm_alloc)	break;
			ap_inhohapke_mgpu_f_krnl<<<BLKpx[f],THD,0,gpu1_stream[hf]>>>(dmod,
					pos1, nThreads[f], body, xylim1, span[f],
					intensity_factor1,	phase_f1, hf);
			hf++;
		}
		checkErrorAfterKernelLaunch("ap_inhohapke_mgpu_krnl");
		break;
	case KAASALAINEN:
		/* Set up needed arrays first, on both GPUs */
		gpuErrchk(cudaSetDevice(GPU0));
		ap_kaas_init_mgpu_f_krnl<<<BLK_half0,THD64,0,gpu0_stream[0]>>>(dmod,
				phasefuncf0,phase_f0,scale_lommsee0,scale_lambert0,nfrm_half0);
		gpuErrchk(cudaSetDevice(GPU1));
		ap_kaas_init_mgpu_f_krnl<<<BLK_half1,THD64,0,gpu1_stream[0]>>>(dmod,
				phasefuncf1,phase_f1,scale_lommsee1,scale_lambert1,nfrm_half1);
		checkErrorAfterKernelLaunch("ap_kaas_init_mgpu_f_krnl");

		/* Launch the main Kaasalainen kernel */
		hf = 0;
		for (f=1; f<nfrm_alloc; f++){
			gpuErrchk(cudaSetDevice(GPU0));
			ap_kaas_mgpu_f_krnl<<<BLKpx[f],THD,0,gpu0_stream[hf]>>>(dmod, pos0,
					nThreads[f],body,xylim0,span[f],intensity_factor0,phase_f0,
					phasefuncf0, scale_lommsee0, scale_lambert0, hf);
			gpuErrchk(cudaSetDevice(GPU1));
			f++;	if (f>=nfrm_alloc)	break;
			ap_kaas_mgpu_f_krnl<<<BLKpx[f],THD,0,gpu1_stream[hf]>>>(dmod, pos1,
					nThreads[f],body,xylim1,span[f],intensity_factor1,phase_f1,
					phasefuncf1, scale_lommsee1, scale_lambert1, hf);
			hf++;
		}
		checkErrorAfterKernelLaunch("ap_kaas_mgpu_f_krnl");

		break;
	case HARMKAAS:
		/* Launch the HarmKaas kernel */
		hf = 0;
		for (f=1; f<nfrm_alloc; f++) {
			gpuErrchk(cudaSetDevice(GPU0));
			ap_harmkaas_mgpu_f_krnl<<<BLKpx[f],THD,0,gpu0_stream[hf]>>>(dmod,
					pos0, nThreads[f], body, xylim0, span[f], intensity_factor0,
					phase_f0, hf);
			gpuErrchk(cudaSetDevice(GPU1));
			f++;	if (f>=nfrm_alloc)	break;
			ap_harmkaas_mgpu_f_krnl<<<BLKpx[f],THD,0,gpu1_stream[hf]>>>(dmod,
					pos1, nThreads[f], body, xylim1, span[f], intensity_factor1,
					phase_f1, hf);
			hf++;
		}
		checkErrorAfterKernelLaunch("ap_harmkaas_mgpu_krnl");
		break;
	case INHOKAAS:
		/* Launch the InhoKaas kernel */
		hf = 0;
		for (f=1; f<nfrm_alloc; f++) {
			gpuErrchk(cudaSetDevice(GPU0));
			ap_inhokaas_mgpu_f_krnl<<<BLKpx[f],THD,0,gpu0_stream[hf]>>>(dmod,
					pos0, nThreads[f], body, xylim0, span[f], intensity_factor0,
					phase_f0, hf);
			gpuErrchk(cudaSetDevice(GPU0));
			f++;	if (f>=nfrm_alloc)	break;
			ap_inhokaas_mgpu_f_krnl<<<BLKpx[f],THD,0,gpu1_stream[hf]>>>(dmod,
					pos1, nThreads[f], body, xylim1, span[f], intensity_factor1,
					phase_f1, hf);
			hf++;
		}
		checkErrorAfterKernelLaunch("ap_inhokaas_mgpu_krnl");
		break;
	case NOLAW:
		bailout("apply_photo.c: can't set optical scattering law = \"none\" when optical data are used\n");
		break;
	default:
		bailout("apply_photo.c: can't handle that optical scattering law yet\n");
	}

	/* Call a streamed parallel reduction which calculates the sums of pos->b
	 * for all frames in a dataset (up to 4 simultaneously)	 */
	sum_brightness_mgpu(ddat, pos0, pos1, nfrm_alloc, nfrm_half0, nfrm_half1,
			nThreadspx[1], set, gpu0_stream, gpu1_stream);
	//sum_brightness_streams(ddat, pos, nframes, nThreadspx[1], 1, set, ap_stream);


	/* Host deallocations */
	free(htype);
	free(hsum);

	/* GPU0 deallocations */
	gpuErrchk(cudaSetDevice(GPU0));
	gpuErrchk(cudaFree(type));
	gpuErrchk(cudaFree(sum));
	gpuErrchk(cudaFree(dsum0));
	gpuErrchk(cudaFree(intensity_factor0));
	gpuErrchk(cudaFree(phase_f0));
	gpuErrchk(cudaFree(phasefuncf0));
	gpuErrchk(cudaFree(scale_lommsee0));
	gpuErrchk(cudaFree(scale_lambert0));

	/* Allocate GPU1 memory */
	gpuErrchk(cudaSetDevice(GPU1));
	gpuErrchk(cudaFree(dsum1));
	gpuErrchk(cudaFree(intensity_factor1));
	gpuErrchk(cudaFree(phase_f1));
	gpuErrchk(cudaFree(phasefuncf1));
	gpuErrchk(cudaFree(scale_lommsee1));
	gpuErrchk(cudaFree(scale_lambert1));

	gpuErrchk(cudaSetDevice(GPU0));

}


#undef TINY
