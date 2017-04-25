/*______________________________________________________________________________*
 * This is the header file for shape-cuda.										*
 * 																				*
 * Written: Thursday, July 7, 2016 by Matt Engels								*
 *______________________________________________________________________________*/

// includes, cuda
#include <cuda.h>
#include <cuda_runtime_api.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <math.h>


/*	Macro to call the error check inline function for all things Cuda */
#define gpuErrchk(ans) do{ gpuAssert((ans), __FILE__, __LINE__); }while(0)

#define cudaCalloc(A, B, C) \
    do { \
    	gpuErrchk(cudaMallocManaged(A, B*C, cudaMemAttachGlobal)); \
		gpuErrchk(cudaMemset(*A, 0, B*C)); \
    } while (0)

#define cudaCalloc1(A, B, C) \
    do { \
    	gpuErrchk(cudaMallocManaged(A, B*C, cudaMemAttachGlobal)); \
    } while (0)

#define MAXIMP 10

/* Flags */
extern int CUDA;			/* Use CUDA 									*/
extern int STREAMS;			/* Use CUDA streams 							*/
extern int STREAMS2;		/* Slightly different iteration of STREAMS		*/
extern int GPU;				/* Which GPU to use 							*/
extern int DYNPROC; 		/* Use dynamic processing (launch child kernels)*/
extern int POSVIS_SEPARATE;	/* Use separate xlim/ylim calculation 			*/
extern int AF;				/* Process all frames in a set at once 			*/
extern int TIMING;			/* Time execution of certain kernels 			*/
extern int FLOAT;			/* Uses singles (floats) instead of doubles in much of the calculation kernels */


extern int maxThreadsPerBlock;
extern struct par_t *dev_par;
extern struct mod_t *dev_mod;
extern struct dat_t *dev_dat;
extern float *dfit, *ddfit;		/* pointers for doppler fit and delay doppler fit */

extern double *fparstep;		/* par->fparstep 	*/
extern double *fpartol;			/* par->partol 		*/
extern double *fparabstol;		/* par->fparabstol 	*/
extern double **fpntr;			/* par->pntr 		*/
extern int *fpartype;			/* par->fpartype	*/

void CUDACount();
void allocate_CUDA_structs(struct par_t par, struct mod_t mod, struct dat_t dat);
void checkErrorAfterKernelLaunch(char *location);
void cotrans_cuda( double y[3], double a[3][3], double x[3], int dir);
void deviceSyncAfterKernelLaunch(char *location);
void mmmul_cuda( double *x, double y[3][3], double *z);
void mtrnsps_cuda( double *a, double b[3][3]);
void pickGPU(int gpuid);

int posvis_cuda( struct vertices_t *verts, double orbit_offset[3], struct pos_t *pos,
		int smooth, int src, int body, int comp);
int pos2deldop_cuda(struct par_t *par, struct photo_t *photo, double orbit_xoff,
		double orbit_yoff, double orbit_dopoff, struct deldop_t *deldop, int body,
		int set, int frm, int v);
int pos2doppler_cuda( struct par_t *par, struct photo_t *photo,
        double orbit_xoff, double orbit_yoff, double orbit_dopoff,
        struct doppler_t *doppler, int body, int set, int frm, int v);

__host__ double apply_photo_cuda(struct mod_t *dmod, struct dat_t *ddat,
		int body, int set, int frm);
__host__ void apply_photo_cuda_streams(struct mod_t *dmod, struct dat_t *ddat,
		struct pos_t **pos, int4 *xylim, int2 *span, dim3 *BLKpx, int *nThreads,
		int body, int set, int nframes, int *nThreadsPx, cudaStream_t *ap_stream);
__host__ void apply_photo_cuda_streams_f(struct mod_t *dmod, struct dat_t *ddat,
		struct pos_t **pos,	int4 *xylim, int2 *span, dim3 *BLKpx, int *nThreads,
		int body, int set, int nframes, cudaStream_t *ap_stream);
__host__ double bestfit_CUDA(struct par_t *dpar, struct mod_t *dmod, struct dat_t
		*ddat, struct par_t *par, struct mod_t *mod, struct dat_t *dat);
__host__ double bestfit_CUDA2(struct par_t *dpar, struct mod_t *dmod, struct
		dat_t *ddat, struct par_t *par, struct mod_t *mod, struct dat_t *dat);
__host__ double brent_abs_streams(double ax,double bx,double cx,double (*f)(double,
		struct vertices_t**, unsigned char*, unsigned char*, int*, int*, int*, int,
		int, cudaStream_t*), double tol,double abstol,double *xmin, struct
		vertices_t **verts, unsigned char *htype, unsigned char *dtype,
		int *nframes, int *nviews, int *lc_n, int nsets, int nf,
		cudaStream_t *bf_stream);
__host__ void c2af_deldop_add_o2_m2(float **temp_o2, float **temp_m2,
		float **temp_om, int size, int nframes);
__host__ void calc_fits_cuda(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat);
__host__ void calc_fits_cuda_af(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat);
__host__ void calc_fits_cuda_streams(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat);
__host__ void calc_fits_cuda_streams2(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, struct vertices_t **verts, int *nviews, int *nframes,
		int *lc_n, unsigned char *type, int nsets, int nf, cudaStream_t *cf_stream);
__host__ double chi2_cuda(struct par_t *dpar, struct dat_t *ddat, int list_breakdown);
__host__ double chi2_cuda_af(struct par_t *dpar,struct dat_t *ddat,
		int list_breakdown, int nsets);
__host__ double chi2_cuda_streams(struct par_t *dpar, struct dat_t *ddat,
		int list_breakdown, int nsets);
__host__ double chi2_cuda_streams2(struct par_t *dpar, struct dat_t *ddat,
		unsigned char *htype, unsigned char *dtype, int *nframes, int *lc_n,
		int list_breakdown,	int nsets, cudaStream_t *c2s_stream);
__host__ double chi2_lghtcrv_cuda_lghtcrv(struct par_t *dpar, struct dat_t *ddat,
		int s, int list_breakdown, double *chi2_all_lghtcrv, int nframes,int lc_n);
__host__ void compute_dv_dcom_dI_reduction(float *dv, float *dcom0, float
		*dcom1, float *dcom2, float *dI00, float *dI01, float *dI02, float
		*dI10, float *dI11, float *dI12, float *dI20, float *dI21, float *dI22,
		 int c, int size, struct mod_t *dmod);
__host__ float compute_doppler_xsec(struct dat_t *ddat, int ndop,
		int set, int frm);
__host__ float compute_deldop_xsec_all_frames(struct dat_t *ddat, int ndel, int ndop,
		int set, int nframes);
__host__ float compute_deldop_xsec_pr6(struct dat_t *ddat, int ndel, int ndop, int set, int frm);
__host__ float compute_deldop_xsec_snglkrnl(struct dat_t *ddat, int ndel, int ndop,
		int set, int frm);
__host__ float compute_pos_zmax(struct dat_t *ddat, int size,
		int set, int frm);
__host__ float compute_pos_zmax_all_frames(struct dat_t *ddat, int frame_size,
		int set, int nframes);
__host__ float compute_pos_zmax_all_frames_2(struct dat_t *ddat, int size,
		int set, int nframes);
__host__ float compute_model_area(struct mod_t *dmod, int c, int size);
__host__ float compute_model_area1(struct mod_t *dmod, int c, int size);
__host__ void compute_xlim_ylim(struct dat_t *ddat, int size, int set, int frm,
		int src, float *iminflt, float *imaxflt, float *jminflt,float *jmaxflt,
		float *minmax_overall);
__host__ void deldopoffs_cuda(struct dat_t *ddat, int s);
__host__ void deldopoffs_cuda_f(struct dat_t *ddat, int s, int nframes);
__host__ void dopoffs_cuda(struct dat_t *ddat, int s);
__host__ void dopoffs_cuda_f(struct dat_t *ddat, int s, int nframes);
__host__ void dvdI_reduce_single(struct mod_t *dmod, float *dv, float *dcom0,
		float *dcom1, float *dcom2, float *dI00, float *dI01, float *dI02,
		float *dI10, float *dI11, float *dI12, float *dI20, float *dI21,
		float *dI22, int size, int c);
__host__ void gpuAssert(cudaError_t code, const char *file, int line);
__host__ void mnbrak_streams(double *ax,double *bx,double *cx,double *fa,double *fb,
	    double *fc,double (*func)(double, struct vertices_t**, unsigned char*,
	    		unsigned char*, int*, int*, int*, int, int, cudaStream_t*), struct
	    		vertices_t **verts, unsigned char *htype, unsigned char *dtype,
	    		int *nframes, int *nviews, int *lc_n, int nsets, int nf,
	    		cudaStream_t *bf_stream );
__host__ void mkparlist_cuda(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, double *fparstep, double *fpartol,
		double *fparabstol, int *fpartype, double **fpntr);
__host__ void mkparlist_cuda2(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, double *fparstep, double *fpartol,
		double *fparabstol, int *fpartype, double **fpntr,
		int nfpar, int nsets);
__host__ double penalties_cuda(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat);
__host__ int pos2deldop_cuda_2(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, double orbit_xoff, double orbit_yoff, double
		orbit_dopoff, int body, int set, int frm, int v);
__host__ int pos2deldop_cuda_af(struct par_t *dpar, struct mod_t *dmod, struct
		dat_t *ddat, double orbit_xoff, double orbit_yoff, double
		orbit_dopoff, int body, int set, int nframes, int v);
__host__ int pos2deldop_cuda_streams(struct par_t *dpar, struct mod_t *dmod, struct
		dat_t *ddat, struct pos_t **pos, int *ndel, int *ndop, double orbit_xoff,
		double orbit_yoff, double orbit_dopoff, int body, int set, int nframes,
		int v, int *badradararr, cudaStream_t *p2d_stream);
__host__ int pos2deldop_cuda_streams_f(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, struct pos_t **pos, int *ndel, int *ndop,
		float3 orbit_xydopoff, int body,
		int set, int nframes, int v, int *badradararr, cudaStream_t *p2d_stream);
__host__ int pos2doppler_cuda_2( struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, double orbit_xoff, double orbit_yoff, double
		orbit_dopoff, int body, int set, int frm, int v);
__host__ int pos2doppler_cuda_af( struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, double orbit_xoff, double orbit_yoff, double
		orbit_dopoff, int body, int set, int nframes, int v);
__host__ int pos2doppler_cuda_streams(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, struct pos_t **pos, double orbit_xoff, double
		orbit_yoff, double orbit_dopoff, int *ndop, int body, int set,
		int nframes, int v,	int *badradararr, cudaStream_t *pds_stream);
__host__ int pos2doppler_cuda_streams_f(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, struct pos_t **pos, float3 orbit_xydopoff, int *ndop,
		int body, int set, int nframes, int v, int *badradararr, cudaStream_t *pds_stream);
__host__ int posvis_cuda_2(struct par_t *dpar, struct mod_t *dmod, struct
		dat_t *ddat, double orbit_offset[3], int set, int frame, int src,
		int body, int comp);
__host__ int posvis_af(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, float orbit_offset[3], int set, int nframes,
		int src, int body, int comp);
__host__ int posvis_cuda_streams(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, double orbit_offset[3], int set, int nframes,
		int src, int body, int comp, int *outndarr, cudaStream_t *posvis_stream);
__host__ int posvis_cuda_streams2(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, struct pos_t **pos, struct vertices_t **verts,
		float3 orbit_offset, int *posn, int *outbndarr, int set, int nframes,
		int src, int nf, int body, int comp, unsigned char type,
		cudaStream_t *posvis_stream);
__host__ void realize_delcor_cuda(struct dat_t *ddat, double delta_delcor0,
		int delcor0_mode, int nsets);
__host__ void realize_delcor_cuda_f(struct dat_t *ddat, double delta_delcor0,
		int delcor0_mode, int nsets, unsigned char *type, int *nframes);
__host__ void realize_dopscale_cuda(struct par_t *dpar, struct dat_t *ddat,
		double dopscale_factor, int dopscale_mode);
__host__ void realize_dopscale_cuda_streams(struct par_t *dpar, struct dat_t
		*ddat, double dopscale_factor, int dopscale_mode, int nsets,
		unsigned char *dtype);
__host__ void realize_mod_cuda( struct par_t *dpar, struct mod_t *dmod,
		unsigned char type);
__host__ void realize_photo_cuda( struct par_t *dpar, struct mod_t *dmod,
        double radalb_factor, double optalb_factor, int albedo_mode);
__host__ void realize_spin_cuda( struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, int nsets);
__host__ void realize_spin_cuda_af( struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, int nsets);
__host__ void realize_spin_cuda_streams( struct par_t *dpar,
		struct mod_t *dmod, struct dat_t *ddat, int nsets);
__host__ void realize_spin_cuda_streams2(struct par_t *dpar,struct mod_t *dmod,
		struct dat_t *ddat, unsigned char *htype, int *nframes, int *nviews,
		int nsets, cudaStream_t *rs_stream);
__host__ void realize_spin_cuda_streams2f(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, unsigned char *htype, int *nframes, int *nviews,
		int nsets, cudaStream_t *rs_stream);
__host__ void realize_xyoff_cuda( struct dat_t *ddat);
__host__ void realize_xyoff_cuda_streams( struct dat_t *ddat, int nsets,
		unsigned char *dtype);
__host__ void show_deldoplim_cuda(struct dat_t *dat, struct dat_t *ddat);
__host__ void show_deldoplim_cuda_streams(struct dat_t *ddat,
		unsigned char *type, int nsets, int *nframes, int maxframes);
__host__ void sum_2_double_arrays(double *a, double *b, double *absum, int size);
__host__ double sum_brightness(struct pos_t **pos, int i, int size);
__host__ float sum_brightness_f(struct pos_t **pos, int i, int size);
__host__ void vary_params_cuda(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, int action, double *deldop_zmax,
		double *rad_xsec, double *opt_brightness, double *cos_subradarlat,
		int nsets);
__host__ void vary_params_af( struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, int action, double *deldop_zmax, double
		*rad_xsec, double *opt_brightness, double *cos_subradarlat, int nsets);
__host__ void vary_params_cuda_streams( struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, int action, double *deldop_zmax, double
		*rad_xsec, double *opt_brightness, double *cos_subradarlat, int nsets);
__host__ void vary_params_cuda_streams2( struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, int action, double *deldop_zmax, double
		*rad_xsec, double *opt_brightness, double *cos_subradarlat, int nsets);
__host__ void vary_params_cuda_streams3(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, int action, double *deldop_zmax, double *rad_xsec,
		double *opt_brightness, double *cos_subradarlat, int *hnframes, int *hlc_n,
		int *nviews, struct vertices_t **verts, unsigned char *htype,
		unsigned char *dtype, int nf, int nsets, cudaStream_t *vp_stream);
__host__ void vary_params_cuda_streams3f(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, int action, double *deldop_zmax, double *rad_xsec,
		double *opt_brightness, double *cos_subradarlat, int *hnframes, int *hlc_n,
		int *nviews, struct vertices_t **verts, unsigned char *htype,
		unsigned char *dtype, int nf, int nsets, cudaStream_t *vp_stream);

__device__ int cubic_realroots_cuda( double *coeff, double *realroot);
__device__ int cubic_realroots_f_cuda( float4 coeff, float3 *realroot);
__device__ void dev_bsstep(double *y, double *dydx, int nv, double *xx, double htry, double eps,
		double *yscal, double *hdid, double *hnext, void (*derivs)(double,double *,double *));
__device__ double dev_cel(double qqc, double pp, double aa, double bb);
__device__ void dev_cotrans1( double y[3], double *a, double x[3], int dir);
__device__ void dev_cotrans2( double y[3], double a[3][3], double x[3], int dir);
__device__ void dev_cotrans3(float y[3], float3 *a, float x[3], int dir, int frm);
__device__ void dev_cotrans4(float3 *y, double a[3][3], double x[3], int dir, int f);
__device__ void dev_cotrans5(double3 *y, double a[3][3], double3 x, int dir);
__device__ void dev_cotrans6(double y[3], double3 *a, double x[3], int dir, int f);
__device__ void dev_cotrans7(float3 *y, double3 *a, float3 x, int dir, int frm);
__device__ void dev_cotrans8(float3 *y, float3 *a, float3 x, int dir, int frm);
__device__ void dev_cotrans9(float3 *y, double a[3][3], float3 x, int dir);
__device__ double dev_cross( double z[3], double x[3], double y[3]);
__device__ float dev_cross_f( double z[3], float3 x, float3 y);
__device__ double dev_dot( double x[3], double y[3]);
__device__ double dev_dot2( double x[3], double3 y);
__device__ double dev_dot3( float3 x, double3 y);
__device__ float dev_dot4(float3 x, float3 y);
__device__ void dev_euler2mat( double m[3][3], double phi, double theta, double psi);
__device__ void dev_facmom( double fv0[3], double fv1[3], double fv2[3], double fn[3],
        double *dv, double dvr[3], double dI[3][3]);
__device__ void dev_facmom_f( float3 fv0, float3 fv1, float3 fv2, float3 fn,
             float *dv, float3 *dvr, float3 dI[3]);
__device__ double dev_facnrm( struct vertices_t verts, int fi);
__device__ float dev_facnrm_f( struct vertices_t verts, int fi);
__device__ int dev_gamma_trans(float *datum, double gamma);
__device__ int dev_gamma_trans_f(float *datum, float gamma);
__device__ double dev_gammln(double xx);
__device__ double dev_hapke( double cosi, double cose, double phase,
        double w, double h, double B0, double g, double theta);
__device__ double dev_hapke_f(float cosi, float cose, float phase,
		float w, float h, float B0, float g, float theta);
__device__ void dev_inteuler( struct spin_t spin, double t[], double impulse[][3], int n,
		double w[3], double m[3][3], unsigned char pa, unsigned char method, double int_abstol);
__device__ int dev_vp_iround(double x);
__device__ int dev_vp_iroundf(float x);
__device__ void dev_mat2euler( double m[3][3], double *phi, double *theta, double *psi);
__device__ void dev_mmid( double *y, double *dydx, int nvar1, double xs, double htot,
		int nstep, double *yout, void (*dev_derivs)( double, double *, double *));
__device__ void dev_mmmul(double x[3][3], double y[3][3], double z[3][3]);
__device__ void dev_mmmul2(double3 *x, double y[3][3], double3 *z, int f);
__device__ void dev_mmmul3(float3 *x, double y[3][3], float3 *z, int frm);
__device__ void dev_mmmul4( float x[3][3], double y[3][3], float z[3][3]);
__device__ void dev_mtrnsps( double a[3][3], double b[3][3]);
__device__ void dev_mtrnsps2(double3 *a, double b[3][3], int f);
__device__ void dev_mtrnsps3(float3 *a, double b[3][3], int frm);
__device__ void dev_mtrnsps4( float a[3][3], double b[3][3]);
__device__ double dev_normalize(double *u);
__device__ float dev_normalize2(float3 u);
__device__ void dev_odeint( double *ystart, int nvar, double x1, double x2, double eps,
	double h1, double hmin, int *nok, int *nbad, void (*derivs)(double,double *,double *),
	void (*drkqc)(double *,double *,int,double *,double,double,double
			*,double *,double *,void (*)(double,double *,double *)));
__device__ double dev_plgndr(int l,int m,double x);
__device__ void dev_POSrect2(struct pos_t *pos, int src, float imin_dbl,
		float imax_dbl, float jmin_dbl, float jmax_dbl);
__device__ double dev_radlaw( struct photo_t *photo, int ilaw, double cosinc, int c, int f);
__device__ double dev_radlaw_f( struct photo_t *photo, int ilaw, float cosinc, int c, int f);
__device__ void dev_realize_impulse(struct spin_t spin, double t,double t_integrate[], double impulse[][3], int *n_integrate, int s, int f, int k);
__device__ void dev_rzextr( int iest, double xest, double *yest, double *yz, double *dy);
__device__ double radlaw_cuda(union radscat_t *radar, unsigned char *radtype,
		int ilaw, double cosinc, int c, int f);

__global__ void c2_add_chi2_krnl(struct dat_t *ddat, int s);
__global__ void c2_set_chi2_krnl(struct dat_t *ddat, double chi2, int s);
__global__ void cf_init_seen_flags_krnl(struct mod_t *dmod, int nf);
__global__ void cf_set_final_pars_krnl(struct par_t *dpar, struct dat_t *ddat);
__global__ void clrvect_krnl(struct dat_t *ddat, int s, int f, int nThreads);
__global__ void clrvect_af_krnl(struct dat_t *ddat, int s, int nframes,	int nThreads, int frame_size);
__global__ void euler2mat_krnl( double m[3][3], double phi, double theta, double psi);
__global__ void euler2mat_realize_mod_krnl(struct mod_t *dmod);
__global__ void get_types_krnl(struct dat_t *ddat, unsigned char *dtype);
__global__ void lghtcrv_spline_streams_krnl(struct dat_t *ddat, int set, double
		yp1, double ypn, double *u, int ncalc);
__global__ void lghtcrv_spline_streams_test_krnl(struct dat_t *ddat, int set, double
		yp1, double ypn, double *u, int ncalc);
__global__ void lghtcrv_spline_streams_f_krnl(struct dat_t *ddat, int set,
		float yp1, float ypn, float *u, int ncalc);
__global__ void lghtcrv_splint_streams3_krnl(struct dat_t *ddat, int set, int n, int ncalc);
__global__ void lghtcrv_splint_streams3f_krnl(struct dat_t *ddat, int set, int ncalc);
__global__ void lghtcrv_splint_streams3_test_krnl(struct dat_t *ddat, int set, int n, int ncalc);
__global__ void posclr_streams_krnl(struct pos_t **pos, int *posn, int f);
__global__ void posmask_universal_krnl(struct par_t *dpar, struct pos_t *pos,
		int nThreads, int xspan);
__global__ void posmask_init_streams_krnl(struct pos_t **pos, double3 *so,
		double *pixels_per_km, int f);
__global__ void posmask_init_streams_f_krnl(struct pos_t **pos, float3 *so,
		float *pixels_per_km, int f);
__global__ void posmask_init_streams2_krnl(struct pos_t **pos, double3 *so,
		float *pixels_per_km, int f);
__global__ void posmask_streams_krnl(struct par_t *dpar, struct pos_t **pos,
		double3 *so, double *pixels_per_km,	int *posn, int nThreads, int xspan,	int f);
__global__ void posmask_streams_f_krnl(struct par_t *dpar, struct pos_t **pos,
		float3 *so,	float *pixels_per_km, int nThreads,	int xspan,	int f);
__global__ void posmask_streams2_krnl(struct par_t *dpar,struct pos_t **pos,double3 *so,
		float *pixels_per_km, int *posn, int nThreads, int xspan, int f);
__global__ void realize_angleoff_krnl(struct dat_t *ddat);
__global__ void realize_omegaoff_krnl(struct dat_t *ddat);
__global__ void update_spin_angle_krnl(struct mod_t *dmod);
__global__ void dbg_vertex_nrmls_krnl(struct mod_t *dmod, int *nafnas);

__host__ void dbg_print_fit(struct dat_t *ddat, int s, int f, char *filename);
void dbg_print_fit_host(struct dat_t *ddat, int s, int f);
__host__ void dbg_print_deldop_fit(struct dat_t *ddat, int s, int , char *filename);
void dbg_print_deldop_fit_host(struct dat_t *ddat, int s, int f);
__host__ void dbg_print_RandC(struct mod_t *dmod);
void dbg_print_RandC_host(struct mod_t *dmod);
__host__ void dbg_print_array(float *data, int x, int y);
__host__ void dbg_print_array1D(float *data, int size);
__host__ void dbg_sum_up_pos(struct dat_t *ddat, int s, int f);
void dbg_sum_up_pos_host(struct dat_t *ddat, int s, int f);
void reduce(int size, int threads, int blocks, int whichKernel, float *d_idata, float *d_odata);
void maxz(int size, int threads, int blocks, int whichKernel, float *d_idata, float *d_odata);
int isPow2(unsigned int x);
unsigned int nextPow2(unsigned int x);
float2 getNumBlocksAndThreads(int n, int maxBlocks, int maxThreads);
__host__ void dbg_check_array_for_content(float *in, int size);
__host__ void dbg_print_array1(float *in, int size);
__host__ void dbg_print_pos_z(struct dat_t *ddat, int set, int frm, int n, char *filename);
__host__ void dbg_print_pos_cose_s(struct dat_t *ddat, int set, int frm, int n, char *filename);
__host__ void dbg_print_pos_z_af(struct dat_t *ddat, int set, int n);
__host__ void dbg_print_cose_af(struct dat_t *ddat, int set, int n);
__host__ void dbg_print_array1D_dbl(double *data, int size, int offset, char *filename);
__host__ void dbg_print_lghtcrv_arrays(struct dat_t *ddat, int set, int n, char *filename);
__host__ void dbg_print_lghtcrv_arrays_host(struct lghtcrv_t *lghtcrv, int set, int n, char *filename);
__host__ void dbg_print_lghtcrv_xyy2(struct dat_t *ddat, int set, int ncalc, char *filename);
__host__ void dbg_print_lghtcrv_xyy2_host(struct lghtcrv_t *lghtcrv, int set, int ncalc, char *filename);
__host__ void dbg_print_lghtcrv_pos_arrays(struct dat_t *ddat, int set, int f, int npixels, int n);
__host__ void dbg_print_lghtcrv_pos_arrays_host(struct lghtcrv_t *lghtcrv, int f, int set);
__host__ void dbg_print_pos_arrays2(struct pos_t **pos, int f, int npixels, int n);
__host__ void dbg_print_pos_arrays2_host(struct pos_t *pos);
__host__ void dbg_print_pos_z_host(struct pos_t *pos, char *fn);
__host__ void dbg_print_pos_arrays_full(struct pos_t **pos, int f, int npixels, int n);
__host__ void dbg_print_pos_arrays_full_host(struct pos_t *pos);
__host__ void dbg_print_facet_normals_host(struct mod_t *mod, char *fn);
__host__ void dbg_print_facet_normals(struct mod_t *dmod, int nf, char *fn);
__host__ void dbg_print_posfacets(struct pos_t **pos, int f, int posn, char *filename);
__host__ void dbg_print_posfacets_host(struct pos_t *pos, char *fn);
__host__ void dbg_print_pos_bd(struct pos_t **pos, int f, int npixels, int n);
__host__ void dbg_print_lc_fit(struct dat_t *ddat, int s, char *filename_fit, int n);
__host__ void dbg_print_lc_fit_host(struct lghtcrv_t *lghtcrv, char *filename_fit, int n);
