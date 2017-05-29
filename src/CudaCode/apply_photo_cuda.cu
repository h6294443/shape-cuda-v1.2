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

__device__ int ap_ilaw;
__device__ double phasefunc, scale_lommsee, scale_lambert, intensityfactor, phase;
//__device__ float sum, phasefuncf, scale_lommseef, scale_lambertf;

__global__ void ap_init_streams_krnl(
		struct dat_t *ddat,
		struct mod_t *dmod,
		struct pos_t **pos,
		int set, int nframes, unsigned char *type, float *dsum,
		double *intensity_factor,
		double *phase_d) {
	/* nframes-threaded kernel */
	int f = blockIdx.x * blockDim.x + threadIdx.x + 1;
	if (f <= nframes) {
		ap_ilaw = ddat->set[set].desc.lghtcrv.ioptlaw;
		type[0] = dmod->photo.opttype[ap_ilaw];
		dsum[f] = 0.0;
		intensity_factor[f] = (pos[f]->km_per_pixel/AU) * (pos[f]->km_per_pixel/AU);
		phase_d[f] = ddat->set[set].desc.lghtcrv.solar_phase[f];
	}
}
__global__ void ap_init_streams_f_krnl(
		struct dat_t *ddat,
		struct mod_t *dmod,
		struct pos_t **pos,
		int set, int nframes, unsigned char *type, float *dsum,
		float *intensity_factor,
		float *phase_f) {
	/* nframes-threaded kernel */
	int f = blockIdx.x * blockDim.x + threadIdx.x + 1;
	float temp;
	if (f <= nframes) {
		ap_ilaw = ddat->set[set].desc.lghtcrv.ioptlaw;
		type[0] = dmod->photo.opttype[ap_ilaw];
		dsum[f] = 0.0;
		temp = __double2float_rn(pos[f]->km_per_pixel/AU);
		intensity_factor[f] = temp*temp;
		phase_f[f] = __double2float_rn(ddat->set[set].desc.lghtcrv.solar_phase[f]);
	}
}
__global__ void ap_lambertlaw_streams_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		double *intensity_factor,
		int4 *xylim,
		int nThreads,
		int body,
		int2 span,
		int f) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[f].w;
	int j = offset / span.x + xylim[f].y;
	double scale;

	if (offset < nThreads) {
		scale = dmod->photo.optical[ap_ilaw].R.R.val/PIE;

		if (pos[f]->cose_s[offset] > 0.0 && pos[f]->cosi_s[offset] > 0.0
				&& pos[f]->body[i][j] == body) {
			pos[f]->b_d[offset] = intensity_factor[f] * scale * pos[f]->cosi_s[offset];
		}
	}
}
__global__ void ap_lambertlaw_streams_f_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		float *intensity_factor,
		int4 *xylim,
		int nThreads,
		int body,
		int2 span,
		int f) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[f].w;
	int j = offset / span.x + xylim[f].y;
	float scale;

	if (offset < nThreads) {
		scale = __double2float_rn(dmod->photo.optical[ap_ilaw].R.R.val/PIE);

		if (pos[f]->cose_s[offset] > 0.0 && pos[f]->cosi_s[offset] > 0.0
				&& pos[f]->body[i][j] == body) {
			pos[f]->b_s[offset] = intensity_factor[f] * scale * pos[f]->cosi_s[offset];
		}
	}
}
__global__ void ap_harmlambert_streams_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		double *intensity_factor,
		int frm) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[frm].w;
	int j = offset / span.x + xylim[frm].y;
	int c, f;
	double scale;

	if (offset < nThreads) {

		if (pos[frm]->cose_s[offset] > 0.0 && pos[frm]->cosi_s[offset] > 0.0
		 && pos[frm]->body[i][j] == body && pos[frm]->f[i][j] >= 0) {
			c = pos[frm]->comp[i][j];
			f = pos[frm]->f[i][j];
			scale = dmod->photo.optical[ap_ilaw].harmR.local[c][f].R.val/PIE;
			pos[frm]->b_d[offset] = intensity_factor[frm] * scale *
					pos[frm]->cosi_s[offset];
		}
	}
}
__global__ void ap_harmlambert_streams_f_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		float *intensity_factor,
		int frm) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[frm].w;
	int j = offset / span.x + xylim[frm].y;
	int c, f;
	float scale;

	if (offset < nThreads) {

		if (pos[frm]->cose_s[offset] > 0.0 && pos[frm]->cosi_s[offset] > 0.0
		 && pos[frm]->body[i][j] == body && pos[frm]->f[i][j] >= 0) {
			c = pos[frm]->comp[i][j];
			f = pos[frm]->f[i][j];
			scale = __double2float_rn(dmod->photo.optical[ap_ilaw].harmR.local[c][f].R.val)/PIE;
			pos[frm]->b_s[offset] = intensity_factor[frm] * scale *
					pos[frm]->cosi_s[offset];
		}
	}
}
__global__ void ap_inholambert_streams_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		double *intensity_factor,
		int frm) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[frm].w;
	int j = offset / span.x + xylim[frm].y;
	int c, f;
	double scale;

	if (offset < nThreads) {

		if (pos[frm]->cose_s[offset] > 0.0 && pos[frm]->cosi_s[offset] > 0.0
		 && pos[frm]->body[i][j] == body && pos[frm]->f[i][j] >= 0) {

			c = pos[frm]->comp[i][j];
			f = pos[frm]->f[i][j];

			scale = dmod->photo.optical[ap_ilaw].inhoR.local[c][f].R.val/PIE;
			pos[frm]->b_d[offset] = intensity_factor[frm] * scale * pos[frm]->cosi_s[offset];
		}
	}
}
__global__ void ap_inholambert_streams_f_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		float *intensity_factor,
		int frm) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[frm].w;
	int j = offset / span.x + xylim[frm].y;
	int c, f;
	float scale;

	if (offset < nThreads) {

		if (pos[frm]->cose_s[offset] > 0.0 && pos[frm]->cosi_s[offset] > 0.0
		 && pos[frm]->body[i][j] == body && pos[frm]->f[i][j] >= 0) {
			c = pos[frm]->comp[i][j];
			f = pos[frm]->f[i][j];
			scale = __double2float_rn(dmod->photo.optical[ap_ilaw].inhoR.local[c][f].R.val)/PIE;
			pos[frm]->b_s[offset] = intensity_factor[frm] * scale * pos[frm]->cosi_s[offset];
		}
	}
}
__global__ void ap_lommel_streams_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		double *intensity_factor,
		int frm) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[frm].w;
	int j = offset / span.x + xylim[frm].y;
	double scale;

	if (offset < nThreads) {

		scale = dmod->photo.optical[ap_ilaw].R.R.val/(4*PIE);

		if (pos[frm]->cose_s[offset] > 0.0 && pos[frm]->cosi_s[offset] > 0.0
				&& pos[frm]->body[i][j] == body) {
			pos[frm]->b_d[offset] = intensity_factor[frm] * scale * pos[frm]->cosi_s[offset]
			                    / (pos[frm]->cosi_s[offset] + pos[frm]->cose_s[offset]);
		}
	}
}
__global__ void ap_lommel_streams_f_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		float *intensity_factor,
		int frm) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[frm].w;
	int j = offset / span.x + xylim[frm].y;
	float scale;

	if (offset < nThreads) {
		scale = __double2float_rn(dmod->photo.optical[ap_ilaw].R.R.val)/(4*PIE);
		if (pos[frm]->cose_s[offset] > 0.0 && pos[frm]->cosi_s[offset] > 0.0
				&& pos[frm]->body[i][j] == body) {
			pos[frm]->b_s[offset] = intensityfactor * scale * pos[frm]->cosi_s[offset]
			  / (pos[frm]->cosi_s[offset] + pos[frm]->cose_s[offset]);
		}
	}
}
__global__ void ap_harmlommel_streams_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		double *intensity_factor,
		int frm) {

	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[frm].w;
	int j = offset / span.x + xylim[frm].y;
	int c, f;
	double scale;

	if (offset < nThreads) {
		if (pos[frm]->cose_s[offset] > 0.0 && pos[frm]->cosi_s[offset] > 0.0
		 && pos[frm]->body[i][j] == body && pos[frm]->f[i][j] >= 0) {

			c = pos[frm]->comp[i][j];
			f = pos[frm]->f[i][j];

			scale = dmod->photo.optical[ap_ilaw].harmR.local[c][f].R.val/(4*PIE);

			pos[frm]->b_d[offset] = intensity_factor[frm] * scale * pos[frm]->cosi_s[offset]
			   / (pos[frm]->cosi_s[offset] + pos[frm]->cose_s[offset]);
		}
	}
}
__global__ void ap_harmlommel_streams_f_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		float *intensity_factor,
		int frm) {

	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[frm].w;
	int j = offset / span.x + xylim[frm].y;
	int c, f;
	float scale;

	if (offset < nThreads) {
		if (pos[frm]->cose_s[offset] > 0.0 && pos[frm]->cosi_s[offset] > 0.0
		 && pos[frm]->body[i][j] == body && pos[frm]->f[i][j] >= 0) {
			c = pos[frm]->comp[i][j];
			f = pos[frm]->f[i][j];
			scale = __double2float_rn(dmod->photo.optical[ap_ilaw].harmR.local[c][f].R.val)/(4*PIE);
			pos[frm]->b_s[offset] = intensityfactor * scale * pos[frm]->cosi_s[offset]
			   / (pos[frm]->cosi_s[offset] + pos[frm]->cose_s[offset]);
		}
	}
}
__global__ void ap_inholommel_streams_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		double *intensity_factor,
		int frm) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[frm].w;
	int j = offset / span.x + xylim[frm].y;
	int c, f;
	double scale;

	if (offset < nThreads) {

		if (pos[frm]->cose_s[offset] > 0.0 && pos[frm]->cosi_s[offset] > 0.0
		 && pos[frm]->body[i][j] == body && pos[frm]->f[i][j] >= 0) {

			c = pos[frm]->comp[i][j];
			f = pos[frm]->f[i][j];

			scale = dmod->photo.optical[ap_ilaw].inhoR.local[c][f].R.val/(4*PIE);

			pos[frm]->b_d[offset] = intensity_factor[frm] * scale * pos[frm]->cosi_s[offset]
			   / (pos[frm]->cosi_s[offset] + pos[frm]->cose_s[offset]);
		}
	}
}
__global__ void ap_inholommel_streams_f_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		float *intensity_factor,
		int frm) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[frm].w;
	int j = offset / span.x + xylim[frm].y;
	int c, f;
	float scale;

	if (offset < nThreads) {
		if (pos[frm]->cose_s[offset] > 0.0 && pos[frm]->cosi_s[offset] > 0.0
		 && pos[frm]->body[i][j] == body && pos[frm]->f[i][j] >= 0) {
			c = pos[frm]->comp[i][j];
			f = pos[frm]->f[i][j];
			scale = __double2float_rn(dmod->photo.optical[ap_ilaw].inhoR.local[c][f].R.val)/(4*PIE);
			pos[frm]->b_s[offset] = intensityfactor * scale * pos[frm]->cosi_s[offset]
			   / (pos[frm]->cosi_s[offset] + pos[frm]->cose_s[offset]);
		}
	}
}
__global__ void ap_geometrical_streams_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		double *intensity_factor,
		int frm) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[frm].w;
	int j = offset / span.x + xylim[frm].y;

	if (offset < nThreads) {
		if (pos[frm]->cose_s[offset] > 0.0 && pos[frm]->cosi_s[offset] > 0.0
		 && pos[frm]->body[i][j] == body) {

			pos[frm]->b_d[offset] = intensity_factor[frm] * dmod->photo.optical[ap_ilaw].R.R.val;
		}
	}
}
__global__ void ap_geometrical_streams_f_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		float *intensity_factor,
		int frm) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[frm].w;
	int j = offset / span.x + xylim[frm].y;

	if (offset < nThreads) {
		if (pos[frm]->cose_s[offset] > 0.0 && pos[frm]->cosi_s[offset] > 0.0
		 && pos[frm]->body[i][j] == body) {
			pos[frm]->b_s[offset] = intensityfactor * dmod->photo.optical[ap_ilaw].R.R.val;
		}
	}
}
__global__ void ap_hapke_streams_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		double *intensity_factor,
		double *phase_d,
		int frm) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[frm].w;
	int j = offset / span.x + xylim[frm].y;

	if (offset < nThreads) {
		if (pos[frm]->cose_s[offset] > 0.0 && pos[frm]->cosi_s[offset] > 0.0
		 && pos[frm]->body[i][j] == body) {

			pos[frm]->b_d[offset] = intensity_factor[frm]
					* dev_hapke(pos[frm]->cosi_s[offset], pos[frm]->cose_s[offset],
							phase_d[frm],
							dmod->photo.optical[ap_ilaw].hapke.w.val,
							dmod->photo.optical[ap_ilaw].hapke.h.val,
							dmod->photo.optical[ap_ilaw].hapke.B0.val,
							dmod->photo.optical[ap_ilaw].hapke.g.val,
							dmod->photo.optical[ap_ilaw].hapke.theta.val);
		}
	}
}
__global__ void ap_hapke_streams_f_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		float *intensity_factor,
		float *phase_f,
		int frm) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[frm].w;
	int j = offset / span.x + xylim[frm].y;

	if (offset < nThreads) {
		if (pos[frm]->cose_s[offset] > 0.0 && pos[frm]->cosi_s[offset] > 0.0
		 && pos[frm]->body[i][j] == body) {
			pos[frm]->b_s[offset] = intensity_factor[frm]
					* dev_hapke_f(pos[frm]->cosi_s[offset],
							pos[frm]->cose_s[offset],
							phase_f[frm],
							__double2float_rn(dmod->photo.optical[ap_ilaw].hapke.w.val),
							__double2float_rn(dmod->photo.optical[ap_ilaw].hapke.h.val),
							__double2float_rn(dmod->photo.optical[ap_ilaw].hapke.B0.val),
							__double2float_rn(dmod->photo.optical[ap_ilaw].hapke.g.val),
							__double2float_rn(dmod->photo.optical[ap_ilaw].hapke.theta.val));

		}
	}
}
__global__ void ap_harmhapke_streams_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		double *intensity_factor,
		double *phase_d,
		int frm) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[frm].w;
	int j = offset / span.x + xylim[frm].y;
	int c, f;

	if (offset < nThreads) {
		if (pos[frm]->cose_s[offset] > 0.0 && pos[frm]->cosi_s[offset] > 0.0
	     && pos[frm]->body[i][j] == body && pos[frm]->f[i][j] >= 0) {

			c = pos[frm]->comp[i][j];
			f = pos[frm]->f[i][j];
			pos[frm]->b_d[offset] = intensity_factor[frm]
					* dev_hapke(pos[frm]->cosi_s[offset], pos[frm]->cose_s[offset],
							phase_d[frm],
							dmod->photo.optical[ap_ilaw].harmhapke.local[c][f].w.val,
							dmod->photo.optical[ap_ilaw].harmhapke.local[c][f].h.val,
							dmod->photo.optical[ap_ilaw].harmhapke.local[c][f].B0.val,
							dmod->photo.optical[ap_ilaw].harmhapke.local[c][f].g.val,
							dmod->photo.optical[ap_ilaw].harmhapke.local[c][f].theta.val);
		}
	}
}
__global__ void ap_harmhapke_streams_f_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		float *intensity_factor,
		float *phase_f,
		int frm) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[frm].w;
	int j = offset / span.x + xylim[frm].y;
	int c, f;

	if (offset < nThreads) {
		if (pos[frm]->cose_s[offset] > 0.0 && pos[frm]->cosi_s[offset] > 0.0
	     && pos[frm]->body[i][j] == body && pos[frm]->f[i][j] >= 0) {
			c = pos[frm]->comp[i][j];
			f = pos[frm]->f[i][j];
			pos[frm]->b_s[offset] = intensity_factor[frm]
			            *dev_hapke_f(pos[frm]->cosi_s[offset], pos[frm]->cose_s[offset],
			            phase_f[frm],
			            __double2float_rn(dmod->photo.optical[ap_ilaw].harmhapke.local[c][f].w.val),
			            __double2float_rn(dmod->photo.optical[ap_ilaw].harmhapke.local[c][f].h.val),
			            __double2float_rn(dmod->photo.optical[ap_ilaw].harmhapke.local[c][f].B0.val),
			            __double2float_rn(dmod->photo.optical[ap_ilaw].harmhapke.local[c][f].g.val),
			            __double2float_rn(dmod->photo.optical[ap_ilaw].harmhapke.local[c][f].theta.val));
		}
	}
}
__global__ void ap_inhohapke_streams_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		double *intensity_factor,
		double *phase_d,
		int frm) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[frm].w;
	int j = offset / span.x + xylim[frm].y;
	int c, f;

	if (offset < nThreads) {
		if (pos[frm]->cose_s[offset] > 0.0 && pos[frm]->cosi_s[offset] > 0.0
		 && pos[frm]->body[i][j] == body && pos[frm]->f[i][j] >= 0) {

			c = pos[frm]->comp[i][j];
			f = pos[frm]->f[i][j];

			pos[frm]->b_d[offset] = intensity_factor[frm]
					* dev_hapke(pos[frm]->cosi_s[offset], pos[frm]->cose_s[offset],
							phase_d[frm],
							dmod->photo.optical[ap_ilaw].inhohapke.local[c][f].w.val,
							dmod->photo.optical[ap_ilaw].inhohapke.local[c][f].h.val,
							dmod->photo.optical[ap_ilaw].inhohapke.local[c][f].B0.val,
							dmod->photo.optical[ap_ilaw].inhohapke.local[c][f].g.val,
							dmod->photo.optical[ap_ilaw].inhohapke.local[c][f].theta.val);
		}
	}
}
__global__ void ap_inhohapke_streams_f_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		float *intensity_factor,
		float *phase_f,
		int frm) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[frm].w;
	int j = offset / span.x + xylim[frm].y;
	int c, f;

	if (offset < nThreads) {
		if (pos[frm]->cose_s[offset] > 0.0 && pos[frm]->cosi_s[offset] > 0.0
		 && pos[frm]->body[i][j] == body && pos[frm]->f[i][j] >= 0) {
			c = pos[frm]->comp[i][j];
			f = pos[frm]->f[i][j];
			pos[frm]->b_s[offset] = intensity_factor[frm]
					* dev_hapke_f(pos[frm]->cosi_s[offset], pos[frm]->cose_s[offset],
							phase_f[frm],
							__double2float_rn(dmod->photo.optical[ap_ilaw].inhohapke.local[c][f].w.val),
							__double2float_rn(dmod->photo.optical[ap_ilaw].inhohapke.local[c][f].h.val),
							__double2float_rn(dmod->photo.optical[ap_ilaw].inhohapke.local[c][f].B0.val),
							__double2float_rn(dmod->photo.optical[ap_ilaw].inhohapke.local[c][f].g.val),
							__double2float_rn(dmod->photo.optical[ap_ilaw].inhohapke.local[c][f].theta.val));
		}
	}
}
__global__ void ap_kaas_init2_krnl(
		struct mod_t *dmod,
		double *phasefuncd,
		double *phase_d,
		double *scale_lommsee,
		double *scale_lambert,
		int nframes) {
	/* nframes-threaded kernel */
	int frm = blockIdx.x * blockDim.x + threadIdx.x + 1;
	int ap_ilaw = 0;

	if (frm <= nframes) {
		phasefuncd[frm] = dmod->photo.optical[ap_ilaw].kaas.A0.val
				* exp( -phase_d[frm] / dmod->photo.optical[ap_ilaw].kaas.D.val)
		+ dmod->photo.optical[ap_ilaw].kaas.k.val * phase_d[frm] + 1;

		scale_lommsee[frm] = (1 - dmod->photo.optical[ap_ilaw].kaas.wt.val)
				* phasefuncd[frm] * dmod->photo.optical[ap_ilaw].kaas.R.val/(4*PIE);
		scale_lambert[frm] = dmod->photo.optical[ap_ilaw].kaas.wt.val
				* phasefuncd[frm] * dmod->photo.optical[ap_ilaw].kaas.R.val/PIE;
	}
}
__global__ void ap_kaas_init2f_krnl(
		struct mod_t *dmod,
		float *phasefuncd,
		float *phase_d,
		float *scale_lommsee,
		float *scale_lambert,
		int nframes) {
	/* nframes-threaded kernel */
	int frm = blockIdx.x * blockDim.x + threadIdx.x + 1;
	int ap_ilaw = 0;

	if (frm <= nframes) {
		phasefuncd[frm] = dmod->photo.optical[ap_ilaw].kaas.A0.val
				* exp( -phase_d[frm] / dmod->photo.optical[ap_ilaw].kaas.D.val)
		+ dmod->photo.optical[ap_ilaw].kaas.k.val * phase_d[frm] + 1;

		scale_lommsee[frm] = (1 - dmod->photo.optical[ap_ilaw].kaas.wt.val)
				* phasefuncd[frm] * dmod->photo.optical[ap_ilaw].kaas.R.val/(4*PIE);
		scale_lambert[frm] = dmod->photo.optical[ap_ilaw].kaas.wt.val
				* phasefuncd[frm] * dmod->photo.optical[ap_ilaw].kaas.R.val/PIE;
	}
}
__global__ void ap_kaas_streams2_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		double *intensity_factor,
		double *phase_d,
		double *phasefuncd,
		double *scale_lommsee,
		double *scale_lambert,
		int frm,
		int nframes) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[frm].w;
	int j = offset / span.x + xylim[frm].y;
	int n = pos[frm]->n;
	int pos_spn = 2*n+1;
	int pxa = (j+n)*pos_spn + (i+n);

	if (offset < nThreads) {
		if (pos[frm]->cose_s[pxa] > 0.0 && pos[frm]->cosi_s[pxa] > 0.0
				&& pos[frm]->body[i][j] == body) {

			pos[frm]->b_d[pxa] = intensity_factor[frm] * pos[frm]->cosi_s[pxa]
			   *(scale_lommsee[frm] / (pos[frm]->cosi_s[pxa] + pos[frm]->cose_s[pxa])
					   + scale_lambert[frm]);
			}
	}
}
__global__ void ap_kaas_streams2f_krnl(
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
				int frm,
				int nframes) {
			/* Multi-threaded kernel */
			int offset = blockIdx.x * blockDim.x + threadIdx.x;
			int i = offset % span.x + xylim[frm].w;
			int j = offset / span.x + xylim[frm].y;
			int n = pos[frm]->n;
			int pos_spn = 2*n+1;
			int pxa = (j+n)*pos_spn + (i+n);

			if (offset < nThreads) {
				if (pos[frm]->cose_s[pxa] > 0.0 && pos[frm]->cosi_s[pxa] > 0.0
						&& pos[frm]->body[i][j] == body) {

					pos[frm]->b_s[pxa] = intensity_factor[frm] * pos[frm]->cosi_s[pxa]
			       *(scale_lommsee[frm] / (pos[frm]->cosi_s[pxa] + pos[frm]->cose_s[pxa])
			  		  + scale_lambert[frm]);
				}
			}
}
__global__ void ap_harmkaas_streams_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		double *intensity_factor,
		double *phase_d,
		int frm) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[frm].w;
	int j = offset / span.x + xylim[frm].y;
	int c, f;
	double phasefuncd, scale_lommsee, scale_lambert;

	if (offset < nThreads) {
		if (pos[frm]->cose_s[offset] > 0.0 && pos[frm]->cosi_s[offset] > 0.0
		 && pos[frm]->body[i][j] == body && pos[frm]->f[i][j] >= 0) {

			c = pos[frm]->comp[i][j];
			f = pos[frm]->f[i][j];

			phasefuncd = dmod->photo.optical[ap_ilaw].harmkaas.local[c][f].A0.val
			* exp( -phase_d[frm] / dmod->photo.optical[ap_ilaw].harmkaas.local[c][f].D.val)
			+ dmod->photo.optical[ap_ilaw].harmkaas.local[c][f].k.val * phase_d[frm] + 1;

			scale_lommsee = (1 - dmod->photo.optical[ap_ilaw].harmkaas.local[c][f].wt.val)
		    * phasefuncd * dmod->photo.optical[ap_ilaw].harmkaas.local[c][f].R.val / (4*PIE);

			scale_lambert = dmod->photo.optical[ap_ilaw].harmkaas.local[c][f].wt.val
			* phasefuncd * dmod->photo.optical[ap_ilaw].harmkaas.local[c][f].R.val / PIE;

			 pos[frm]->b_d[offset] = intensity_factor[frm] * pos[frm]->cosi_s[offset] * (scale_lommsee /
					( pos[frm]->cosi_s[offset] +  pos[frm]->cose_s[offset]) + scale_lambert);
		}
	}
}
__global__ void ap_harmkaas_streams_f_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		float *intensity_factor,
		float *phase_f,
		int frm) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[frm].w;
	int j = offset / span.x + xylim[frm].y;
	int c, f;
	float phasefuncf, scale_lommseef, scale_lambertf;

	if (offset < nThreads) {
		if (pos[frm]->cose_s[offset] > 0.0 && pos[frm]->cosi_s[offset] > 0.0
		 && pos[frm]->body[i][j] == body && pos[frm]->f[i][j] >= 0) {
			c = pos[frm]->comp[i][j];
			f = pos[frm]->f[i][j];
			phasefuncf = __double2float_rn(dmod->photo.optical[ap_ilaw].harmkaas.local[c][f].A0.val)
			* exp( -phase_f[frm] / __double2float_rn(dmod->photo.optical[ap_ilaw].harmkaas.local[c][f].D.val))
			+ __double2float_rn(dmod->photo.optical[ap_ilaw].harmkaas.local[c][f].k.val) * phase_f[frm] + 1;

			scale_lommseef = (1 - __double2float_rn(dmod->photo.optical[ap_ilaw].harmkaas.local[c][f].wt.val))
		    * phasefuncf * __double2float_rn(dmod->photo.optical[ap_ilaw].harmkaas.local[c][f].R.val) / (4*PIE);
			scale_lambertf = __double2float_rn(dmod->photo.optical[ap_ilaw].harmkaas.local[c][f].wt.val)
			* phasefuncf * __double2float_rn(dmod->photo.optical[ap_ilaw].harmkaas.local[c][f].R.val) / PIE;
			 pos[frm]->b_s[offset] = intensity_factor[frm] * pos[frm]->cosi_s[offset] * (scale_lommseef /
					( pos[frm]->cosi_s[offset] +  pos[frm]->cose_s[offset]) + scale_lambertf);
		}
	}
}
__global__ void ap_inhokaas_streams_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		double *intensity_factor,
		double *phase_d,
		int frm) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[frm].w;
	int j = offset / span.x + xylim[frm].y;
	int c, f;
	double phasefuncd, scale_lommsee, scale_lambert;

	if (offset < nThreads) {
		if (pos[frm]->cose_s[offset] > 0.0 && pos[frm]->cosi_s[offset] > 0.0
		 && pos[frm]->body[i][j] == body && pos[frm]->f[i][j] >= 0) {

			c = pos[frm]->comp[i][j];
			f = pos[frm]->f[i][j];

			phasefuncd = dmod->photo.optical[ap_ilaw].inhokaas.local[c][f].A0.val
			* exp( -phase_d[frm] / dmod->photo.optical[ap_ilaw].inhokaas.local[c][f].D.val)
			+ dmod->photo.optical[ap_ilaw].inhokaas.local[c][f].k.val * phase_d[frm] + 1;

			scale_lommsee = (1 - dmod->photo.optical[ap_ilaw].inhokaas.local[c][f].wt.val)
		    * phasefuncd * dmod->photo.optical[ap_ilaw].inhokaas.local[c][f].R.val / (4*PIE);

			scale_lambert = dmod->photo.optical[ap_ilaw].inhokaas.local[c][f].wt.val
			* phasefuncd * dmod->photo.optical[ap_ilaw].inhokaas.local[c][f].R.val / PIE;

			pos[frm]->b_d[offset] = intensity_factor[frm] * pos[frm]->cosi_s[offset] * (scale_lommsee /
					(pos[frm]->cosi_s[offset] + pos[frm]->cose_s[offset]) + scale_lambert);
		}
	}
}
__global__ void ap_inhokaas_streams_f_krnl(
		struct mod_t *dmod,
		struct pos_t **pos,
		int nThreads,
		int body,
		int4 *xylim,
		int2 span,
		float *intensity_factor,
		float *phase_f,
		int frm) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int i = offset % span.x + xylim[frm].w;
	int j = offset / span.x + xylim[frm].y;
	int c, f;
	float phasefuncf, scale_lommseef, scale_lambertf;

	if (offset < nThreads) {
		if (pos[frm]->cose_s[offset] > 0.0 && pos[frm]->cosi_s[offset] > 0.0
		 && pos[frm]->body[i][j] == body && pos[frm]->f[i][j] >= 0) {
			c = pos[frm]->comp[i][j];
			f = pos[frm]->f[i][j];
			phasefuncf = __double2float_rn(dmod->photo.optical[ap_ilaw].inhokaas.local[c][f].A0.val)
			* exp( -phase_f[frm] / __double2float_rn(dmod->photo.optical[ap_ilaw].inhokaas.local[c][f].D.val))
			+ __double2float_rn(dmod->photo.optical[ap_ilaw].inhokaas.local[c][f].k.val) * phase_f[frm] + 1;
			scale_lommseef = (1 - __double2float_rn(dmod->photo.optical[ap_ilaw].inhokaas.local[c][f].wt.val))
		    * phasefuncf * __double2float_rn(dmod->photo.optical[ap_ilaw].inhokaas.local[c][f].R.val) / (4*PIE);
			scale_lambertf = __double2float_rn(dmod->photo.optical[ap_ilaw].inhokaas.local[c][f].wt.val)
			* phasefuncf * __double2float_rn(dmod->photo.optical[ap_ilaw].inhokaas.local[c][f].R.val) / PIE;
			pos[frm]->b_s[offset] = intensity_factor[frm] * pos[frm]->cosi_s[offset] * (scale_lommseef /
					(pos[frm]->cosi_s[offset] + pos[frm]->cose_s[offset]) + scale_lambertf);
		}
	}
}

__host__ void apply_photo_cuda_streams(
		struct mod_t *dmod,
		struct dat_t *ddat,
		struct pos_t **pos,
		int4 *xylim,
		int2 *span,
		dim3 *BLKpx,
		int *nThreads,
		int body,
		int set,
		int nframes,
		int *nThreadspx,
		cudaStream_t *ap_stream)
{
	unsigned char *type, *htype;
	int f;
	float *dsum;
	double *hsum, *sum, *intensity_factor, *phase_d, *phasefuncd,
		   *scale_lommsee, *scale_lambert;
	dim3 BLK, THD;

	gpuErrchk(cudaMalloc((void**)&type, sizeof(unsigned char) * 2));
	gpuErrchk(cudaMalloc((void**)&sum, sizeof(double) * (nframes+1)));
	gpuErrchk(cudaMalloc((void**)&dsum, sizeof(float)*(nframes+1)));
	gpuErrchk(cudaMalloc((void**)&intensity_factor, sizeof(double)*(nframes+1)));
	gpuErrchk(cudaMalloc((void**)&phase_d, sizeof(double)*(nframes+1)));
	gpuErrchk(cudaMalloc((void**)&scale_lommsee, sizeof(double)*(nframes+1)));
	gpuErrchk(cudaMalloc((void**)&scale_lambert, sizeof(double)*(nframes+1)));
	htype = (unsigned char *) malloc(2*sizeof(unsigned char));
	hsum = (double *) malloc((nframes+1)*sizeof(double));

	/* Launch single-thread kernel to assign pos address and get type */
	THD.x = maxThreadsPerBlock;
	BLK.x = floor((THD.x - 1 + nframes) / THD.x);
	ap_init_streams_krnl<<<BLK,THD>>>(ddat, dmod, pos, set, nframes, type, dsum,
			intensity_factor, phase_d);
	checkErrorAfterKernelLaunch("ap_init_streams_krnl");
	gpuErrchk(cudaMemcpy(htype, type, sizeof(unsigned char) *2,
			cudaMemcpyDeviceToHost));

	switch (htype[0]) {
	case LAMBERTLAW:
		/* Launch Lambert Law kernel */
		for (f=1; f<=nframes; f++)
			ap_lambertlaw_streams_krnl<<<BLKpx[f],THD,0,ap_stream[f-1]>>>(dmod,pos,
					intensity_factor, xylim, nThreads[f], body, span[f], f);
		checkErrorAfterKernelLaunch("ap_lambertlaw_streams_krnl");
		break;
	case HARMLAMBERT:
		/* Launch the HarmLambert kernel */
		for (f=1; f<=nframes; f++)
			ap_harmlambert_streams_krnl<<<BLKpx[f],THD,0,ap_stream[f-1]>>>(
					dmod, pos, nThreads[f], body, xylim, span[f],
					intensity_factor, f);
		checkErrorAfterKernelLaunch("ap_harmlambert_streams_krnl");
		break;
	case INHOLAMBERT:
		/* Launch the Inhomogeneous Lambert kernel */
		for (f=1; f<=nframes; f++)
			ap_inholambert_streams_krnl<<<BLKpx[f],THD,0,ap_stream[f-1]>>>(dmod,
					pos, nThreads[f], body, xylim, span[f], intensity_factor,f);
		checkErrorAfterKernelLaunch("ap_inholambert_streams_krnl");
		break;
	case LOMMEL:
		/* Launch the Lommel kernel */
		for (f=1; f<=nframes; f++)
			ap_lommel_streams_krnl<<<BLKpx[f],THD,0,ap_stream[f-1]>>>(dmod, pos,
					nThreads[f], body, xylim, span[f], intensity_factor, f);
		checkErrorAfterKernelLaunch("ap_lommel_streams_krnl");
		break;
	case HARMLOMMEL:
		/* Launch the HarmLommel kernel */
		for (f=1; f<=nframes; f++)
			ap_harmlommel_streams_krnl<<<BLKpx[f],THD,0,ap_stream[f-1]>>>(dmod,
					pos, nThreads[f], body, xylim, span[f], intensity_factor, f);
		checkErrorAfterKernelLaunch("ap_harmlommel_streams_krnl");
		break;
	case INHOLOMMEL:
		/* Launch the Inhomogeneous Lommel kernel */
		for (f=1; f<=nframes; f++)
			ap_inholommel_streams_krnl<<<BLKpx[f],THD,0,ap_stream[f-1]>>>(dmod,
					pos, nThreads[f], body, xylim, span[f], intensity_factor, f);
		checkErrorAfterKernelLaunch("ap_inholommel_streams_krnl");
		break;
	case GEOMETRICAL:
		/* Launch the Geometrical law kernel */
		for (f=1; f<=nframes; f++)
			ap_geometrical_streams_krnl<<<BLKpx[f],THD,0,ap_stream[f-1]>>>(dmod,
					pos, nThreads[f], body, xylim, span[f], intensity_factor, f);
		checkErrorAfterKernelLaunch("ap_geometrical_streams_krnl");
		break;
	case HAPKE:
		/* Launch the Hapke kernel */
		for (f=1; f<=nframes; f++)
			ap_hapke_streams_krnl<<<BLKpx[f],THD,0,ap_stream[f-1]>>>(dmod, pos,
					nThreads[f], body, xylim, span[f], intensity_factor,
					phase_d, f);
		checkErrorAfterKernelLaunch("ap_hapke_streams_krnl");
		break;
	case HARMHAPKE:
		/* Launch the HarmHapke kernel */
		for (f=1; f<=nframes; f++)
			ap_harmhapke_streams_krnl<<<BLKpx[f],THD,0,ap_stream[f-1]>>>(dmod,
					pos, nThreads[f], body, xylim, span[f],	intensity_factor,
					phase_d, f);
		checkErrorAfterKernelLaunch("ap_harmhapke_streams_krnl");
		break;
	case INHOHAPKE:
		/* Launch the Inhomogeneous Hapke kernel */
		for (f=1; f<=nframes; f++)
			ap_inhohapke_streams_krnl<<<BLKpx[f],THD,0,ap_stream[f-1]>>>(dmod,
					pos, nThreads[f], body, xylim, span[f], intensity_factor,
					phase_d, f);
		checkErrorAfterKernelLaunch("ap_inhohapke_streams_krnl");
		break;
	case KAASALAINEN:
		/* Launch single-thread kernel to init Kaas */
		gpuErrchk(cudaMalloc((void**)&phasefuncd, sizeof(double)*(nframes+1)));
		ap_kaas_init2_krnl<<<BLK,THD>>>(dmod, phasefuncd, phase_d,
				scale_lommsee, scale_lambert, nframes);
		checkErrorAfterKernelLaunch("ap_kaas_init2_krnl");

		/* Launch the main Kaasalainen kernel */
		for (f=1; f<=nframes; f++){
			ap_kaas_streams2_krnl<<<BLKpx[f],THD,0,ap_stream[f-1]>>>(dmod, pos,
					nThreads[f], body, xylim, span[f], intensity_factor,
					phase_d, phasefuncd, scale_lommsee, scale_lambert, f, nframes);
		}
		checkErrorAfterKernelLaunch("ap_kaas_streams_krnl");
		cudaFree(phasefuncd);
		break;
	case HARMKAAS:
		/* Launch the HarmKaas kernel */
		for (f=1; f<=nframes; f++)
			ap_harmkaas_streams_krnl<<<BLKpx[f],THD,0,ap_stream[f-1]>>>(dmod, pos,
					nThreads[f], body, xylim, span[f], intensity_factor, phase_d, f);
		checkErrorAfterKernelLaunch("ap_harmkaas_streams_krnl");
		break;
	case INHOKAAS:
		/* Launch the HarmKaas kernel */
		for (f=1; f<=nframes; f++)
			ap_inhokaas_streams_krnl<<<BLKpx[f],THD,0,ap_stream[f-1]>>>(dmod,
					pos, nThreads[f], body, xylim, span[f], intensity_factor,
					phase_d, f);
		checkErrorAfterKernelLaunch("ap_inhokaas_streams_krnl");
		break;
	case NOLAW:
		bailout("apply_photo.c: can't set optical scattering law = \"none\" when optical data are used\n");
		break;
	default:
		bailout("apply_photo.c: can't handle that optical scattering law yet\n");
	}

	/* Call a streamed parallel reduction which calculates the sums of pos->b
	 * for all frames in a dataset (up to 4 simultaneously)	 */
	sum_brightness_streams(ddat, pos, nframes, nThreadspx[1], 1, set, ap_stream);

	cudaFree(dsum);
	cudaFree(sum);
	cudaFree(type);
	free(htype);
	free(hsum);
}

__host__ void apply_photo_cuda_streams_f(
		struct mod_t *dmod,
		struct dat_t *ddat,
		struct pos_t **pos,
		int4 *xylim,
		int2 *span,
		dim3 *BLKpx,
		int *nThreads,
		int body,
		int set,
		int nframes,
		int *nThreadspx,
		cudaStream_t *ap_stream)
{
	unsigned char *type, *htype;
	int f;
	float *dsum;
	double *hsum, *sum;
	float *intensity_factor, *phase_f, *phasefuncf, *scale_lommsee, *scale_lambert;
	dim3 BLK, THD;

	gpuErrchk(cudaMalloc((void**)&type, sizeof(unsigned char) * 2));
	gpuErrchk(cudaMalloc((void**)&sum, sizeof(double) * (nframes+1)));
	gpuErrchk(cudaMalloc((void**)&dsum, sizeof(float)*(nframes+1)));
	gpuErrchk(cudaMalloc((void**)&intensity_factor, sizeof(float)*(nframes+1)));
	gpuErrchk(cudaMalloc((void**)&phase_f, sizeof(float)*(nframes+1)));
	gpuErrchk(cudaMalloc((void**)&scale_lommsee, sizeof(float)*(nframes+1)));
	gpuErrchk(cudaMalloc((void**)&scale_lambert, sizeof(float)*(nframes+1)));
	htype = (unsigned char *) malloc(2*sizeof(unsigned char));
	hsum = (double *) malloc((nframes+1)*sizeof(double));

	/* Launch single-thread kernel to assign pos address and get type */
	THD.x = maxThreadsPerBlock;
	BLK.x = floor((THD.x - 1 + nframes) / THD.x);
	ap_init_streams_f_krnl<<<BLK,THD>>>(ddat, dmod, pos, set, nframes, type, dsum,
			intensity_factor, phase_f);
	checkErrorAfterKernelLaunch("ap_init_streams_f_krnl");
	gpuErrchk(cudaMemcpy(htype, type, sizeof(unsigned char) *2,
			cudaMemcpyDeviceToHost));

	switch (htype[0]) {
	case LAMBERTLAW:
		/* Launch Lambert Law kernel */
		for (f=1; f<=nframes; f++)
			ap_lambertlaw_streams_f_krnl<<<BLKpx[f],THD,0,ap_stream[f-1]>>>(dmod,pos,
					intensity_factor, xylim, nThreads[f], body, span[f], f);
		checkErrorAfterKernelLaunch("ap_lambertlaw_streams_f_krnl");
		break;
	case HARMLAMBERT:
		/* Launch the HarmLambert kernel */
		for (f=1; f<=nframes; f++)
			ap_harmlambert_streams_f_krnl<<<BLKpx[f],THD,0,ap_stream[f-1]>>>(
					dmod, pos, nThreads[f], body, xylim, span[f],
					intensity_factor, f);
		checkErrorAfterKernelLaunch("ap_harmlambert_streams_f_krnl");
		break;
	case INHOLAMBERT:
		/* Launch the Inhomogeneous Lambert kernel */
		for (f=1; f<=nframes; f++)
			ap_inholambert_streams_f_krnl<<<BLKpx[f],THD,0,ap_stream[f-1]>>>(dmod,
					pos, nThreads[f], body, xylim, span[f], intensity_factor,f);
		checkErrorAfterKernelLaunch("ap_inholambert_streams_f_krnl");
		break;
	case LOMMEL:
		/* Launch the Lommel kernel */
		for (f=1; f<=nframes; f++)
			ap_lommel_streams_f_krnl<<<BLKpx[f],THD,0,ap_stream[f-1]>>>(dmod, pos,
					nThreads[f], body, xylim, span[f], intensity_factor, f);
		checkErrorAfterKernelLaunch("ap_lommel_streams_f_krnl");
		break;
	case HARMLOMMEL:
		/* Launch the HarmLommel kernel */
		for (f=1; f<=nframes; f++)
			ap_harmlommel_streams_f_krnl<<<BLKpx[f],THD,0,ap_stream[f-1]>>>(dmod,
					pos, nThreads[f], body, xylim, span[f], intensity_factor, f);
		checkErrorAfterKernelLaunch("ap_harmlommel_streams_f_krnl");
		break;
	case INHOLOMMEL:
		/* Launch the Inhomogeneous Lommel kernel */
		for (f=1; f<=nframes; f++)
			ap_inholommel_streams_f_krnl<<<BLKpx[f],THD,0,ap_stream[f-1]>>>(dmod,
					pos, nThreads[f], body, xylim, span[f], intensity_factor, f);
		checkErrorAfterKernelLaunch("ap_inholommel_streams_f_krnl");
		break;
	case GEOMETRICAL:
		/* Launch the Geometrical law kernel */
		for (f=1; f<=nframes; f++)
			ap_geometrical_streams_f_krnl<<<BLKpx[f],THD,0,ap_stream[f-1]>>>(dmod,
					pos, nThreads[f], body, xylim, span[f], intensity_factor, f);
		checkErrorAfterKernelLaunch("ap_geometrical_streams_f_krnl");
		break;
	case HAPKE:
		/* Launch the Hapke kernel */
		for (f=1; f<=nframes; f++)
			ap_hapke_streams_f_krnl<<<BLKpx[f],THD,0,ap_stream[f-1]>>>(dmod, pos,
					nThreads[f], body, xylim, span[f], intensity_factor,
					phase_f, f);
		checkErrorAfterKernelLaunch("ap_hapke_streams_f_krnl");
		break;
	case HARMHAPKE:
		/* Launch the HarmHapke kernel */
		for (f=1; f<=nframes; f++)
			ap_harmhapke_streams_f_krnl<<<BLKpx[f],THD,0,ap_stream[f-1]>>>(dmod,
					pos, nThreads[f], body, xylim, span[f],	intensity_factor,
					phase_f, f);
		checkErrorAfterKernelLaunch("ap_harmhapke_streams_krnl");
		break;
	case INHOHAPKE:
		/* Launch the Inhomogeneous Hapke kernel */
		for (f=1; f<=nframes; f++)
			ap_inhohapke_streams_f_krnl<<<BLKpx[f],THD,0,ap_stream[f-1]>>>(dmod,
					pos, nThreads[f], body, xylim, span[f], intensity_factor,
					phase_f, f);
		checkErrorAfterKernelLaunch("ap_inhohapke_streams_krnl");
		break;
	case KAASALAINEN:
		/* Launch single-thread kernel to init Kaas */
		gpuErrchk(cudaMalloc((void**)&phasefuncf, sizeof(float)*(nframes+1)));
		ap_kaas_init2f_krnl<<<BLK,THD>>>(dmod, phasefuncf, phase_f,
				scale_lommsee, scale_lambert, nframes);
		checkErrorAfterKernelLaunch("ap_kaas_init2f_krnl");

		/* Launch the main Kaasalainen kernel */
		for (f=1; f<=nframes; f++){
			ap_kaas_streams2f_krnl<<<BLKpx[f],THD,0,ap_stream[f-1]>>>(dmod, pos,
					nThreads[f], body, xylim, span[f], intensity_factor,
					phase_f, phasefuncf, scale_lommsee, scale_lambert, f, nframes);
		}
		checkErrorAfterKernelLaunch("ap_kaas_streams2f_krnl");
		cudaFree(phasefuncf);
		break;
	case HARMKAAS:
		/* Launch the HarmKaas kernel */
		for (f=1; f<=nframes; f++)
			ap_harmkaas_streams_f_krnl<<<BLKpx[f],THD,0,ap_stream[f-1]>>>(dmod, pos,
					nThreads[f], body, xylim, span[f], intensity_factor, phase_f, f);
		checkErrorAfterKernelLaunch("ap_harmkaas_streams_krnl");
		break;
	case INHOKAAS:
		/* Launch the HarmKaas kernel */
		for (f=1; f<=nframes; f++)
			ap_inhokaas_streams_f_krnl<<<BLKpx[f],THD,0,ap_stream[f-1]>>>(dmod,
					pos, nThreads[f], body, xylim, span[f], intensity_factor,
					phase_f, f);
		checkErrorAfterKernelLaunch("ap_inhokaas_streams_krnl");
		break;
	case NOLAW:
		bailout("apply_photo.c: can't set optical scattering law = \"none\" when optical data are used\n");
		break;
	default:
		bailout("apply_photo.c: can't handle that optical scattering law yet\n");
	}

	/* Call a streamed parallel reduction which calculates the sums of pos->b
	 * for all frames in a dataset (up to 4 simultaneously)	 */
	sum_brightness_streams(ddat, pos, nframes, nThreadspx[1], 1, set, ap_stream);

	cudaFree(dsum);
	cudaFree(sum);
	cudaFree(type);
	cudaFree(intensity_factor);
	cudaFree(phase_f);
	cudaFree(scale_lommsee);
	cudaFree(scale_lambert);
	free(htype);
	free(hsum);
}


#undef TINY
