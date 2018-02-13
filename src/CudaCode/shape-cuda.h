/*______________________________________________________________________________*
 * This is the header file for shape-cuda.										*
 * 																				*
 * Written: Thursday, July 7, 2016 by Matt Engels								*
 *______________________________________________________________________________*/

// includes, cuda
#include <cuda.h>
#include <cuda_runtime_api.h>

/*  Define macros */
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
extern int GPU0, GPU1;		/* Which GPU to use 							*/
extern int TIMING;			/* Time execution of certain kernels 			*/
extern int MGPU;			/* Switch for dual-gpu mode */
extern int PIN;				/* Use pinned host memory instead of GPU memory */
extern int EXP;				/* Experimental flag - currently used for the a
							 * tiled, shared memory posvis implementation */
extern int MFS;				/* Multi-frame set flag. Use only for combining multiple
							 * single-frame sets into a single multi-frame set */

/* Structures */
extern struct par_t *dev_par;
extern struct mod_t *dev_mod;
extern struct dat_t *dev_dat;

extern int maxThreadsPerBlock;
extern double *fparstep;		/* par->fparstep 	*/
extern double *fpartol;			/* par->partol 		*/
extern double *fparabstol;		/* par->fparabstol 	*/
extern double **fpntr;			/* par->pntr 		*/
extern int *fpartype;			/* par->fpartype	*/

/* Device Functions */
void CUDACount(int showCUDAInfo);
void allocate_CUDA_structs(struct par_t par, struct mod_t mod, struct dat_t dat);
void checkErrorAfterKernelLaunch(const char *location);
void deviceSyncAfterKernelLaunch(const char *location);
void pickGPU(int gpuid);

__host__ void apply_photo_gpu(struct mod_t *dmod, struct dat_t *ddat,
		struct pos_t **pos, int4 *xylim, int2 *span, dim3 *BLKpx_bbox,
		int *nThreads, int body, int set, int nframes, int maxthds,
		int4 maxxylim,cudaStream_t *ap_stream);

__host__ double bestfit_gpu(struct par_t *dpar, struct mod_t *dmod, struct
		dat_t *ddat, struct par_t *par, struct mod_t *mod, struct dat_t *dat);

__host__ double bestfit_MFS_gpu(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, struct par_t *par, struct mod_t *mod,
		struct dat_t *dat);

__host__ double bestfit_gpu_pthreads(struct par_t *dpar, struct par_t *dpar1,
		struct mod_t *dmod, struct mod_t *dmod1, struct dat_t *ddat, struct
		dat_t *ddat1, struct par_t *par, struct par_t *par1, struct mod_t *mod,
		struct mod_t *mod1, struct dat_t *dat, struct dat_t *dat1, pthread_t
		thread1, pthread_t thread2);

__host__ double brent_abs_gpu(double ax,double bx,double cx,double (*f)(double,
		struct vertices_t**, unsigned char*, unsigned char*, int*, int*, int*, int,
		int, cudaStream_t*), double tol,double abstol,double *xmin, struct
		vertices_t **verts, unsigned char *htype, unsigned char *dtype,
		int *nframes, int *nviews, int *lc_n, int nsets, int nf,
		cudaStream_t *bf_stream);

__host__ double brent_abs_MFS_gpu(double ax,double bx,double cx,double (*f)(double,
		struct vertices_t**, int*, int, int, cudaStream_t*), double tol,double
		abstol,double *xmin, struct vertices_t **verts, int *nviews, int nsets,
		int nf, cudaStream_t *bf_stream);

__host__ double brent_abs_pthreads(double ax, double bx, double cx, double (*f)
		(double, struct vertices_t**, struct vertices_t**, unsigned char*,
		unsigned char*, unsigned char*, int*, int*, int*, int*, int, int, int,
		pthread_t, pthread_t, cudaStream_t*, cudaStream_t*), double tol, double
		abstol, double *xmin, struct vertices_t **verts0, struct vertices_t
		**verts1, unsigned char *htype, unsigned char *dtype0, unsigned char
		*dtype1, int *nframes, int *nviews, int *lc_n, int *GPUID, int nsets,
		int nf, int max_frames, pthread_t thread1, pthread_t thread2,
		cudaStream_t *gpu0_stream, cudaStream_t *gpu1_stream);

__host__ void calc_fits_gpu(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, struct vertices_t **verts, int *nviews, int
		*nframes, int *lc_n, unsigned char *type, int nsets, int nf,
		cudaStream_t *cf_stream, int max_frames);

__host__ void calc_fits_MFS_gpu(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, struct vertices_t **verts, int *nviews,
		int nsets, int nf, cudaStream_t *cf_stream);

__host__ void calc_fits_pthreads(struct par_t *dpar0, struct par_t *dpar1,
		struct mod_t *dmod0, struct mod_t *dmod1, struct dat_t *ddat0, struct
		dat_t *ddat1, struct vertices_t **verts0, struct vertices_t **verts1,
		int *nviews, int *nframes, int *lc_n, int *GPUID, unsigned char *type,
		int nsets, int nf, int max_frames, pthread_t thread1, pthread_t thread2,
		cudaStream_t *gpu0_stream, cudaStream_t *gpu1_stream);

__host__ double chi2_gpu(struct par_t *dpar, struct dat_t *ddat, unsigned char
		*htype, unsigned char *dtype, int *nframes, int *lc_n, int
		list_breakdown,	int nsets, cudaStream_t *c2s_stream, int max_frames);

__host__ double chi2_MFS_gpu(struct par_t *dpar, struct dat_t *ddat, int
		list_breakdown,	int nsets, cudaStream_t *c2s_stream);

__host__ double chi2_pthreads(struct par_t *dpar0, struct par_t *dpar1, struct
		dat_t *ddat0, struct dat_t *ddat1, unsigned char *htype, unsigned char
		*dtype0, unsigned char *dtype1, int *hnframes, int *hlc_n, int *GPUID,
		int list_breakdown, int nsets, int max_frames, pthread_t thread1,
		pthread_t thread2,cudaStream_t *gpu0_stream,cudaStream_t *gpu1_stream);

__host__ void deldopoffs_gpu(struct dat_t *ddat, int s, int nframes);

__host__ void dopoffs_gpu(struct dat_t *ddat, int s, int nframes);

__host__ void gpuAssert(cudaError_t code, const char *file, int line);

__host__ void mnbrak_gpu(double *ax,double *bx,double *cx,double *fa,double *fb,
	    double *fc,double (*func)(double, struct vertices_t**, unsigned char*,
	    unsigned char*, int*, int*, int*, int, int, cudaStream_t*), struct
	    vertices_t **verts, unsigned char *htype, unsigned char *dtype,
	    int *nframes, int *nviews, int *lc_n, int nsets, int nf,
	    cudaStream_t *bf_stream );

__host__ void mnbrak_MFS_gpu(double *ax, double *bx, double *cx, double *fa,
		double *fb, double *fc,	double (*func)(double, struct vertices_t**,
				int*, int, int, cudaStream_t*), struct vertices_t **verts,
				int *nviews, int nsets, int nf,	cudaStream_t *bf_stream );

__host__ void mnbrak_pthreads(double *ax, double *bx, double *cx, double *fa,
		double *fb, double *fc, double (*func)(double, struct vertices_t**,
		struct vertices_t**, unsigned char*, unsigned char*, unsigned char*,
		int*, int*, int*, int*,	int, int, int, pthread_t, pthread_t,
		cudaStream_t*, cudaStream_t*), struct vertices_t **verts0, struct
		vertices_t **verts1, unsigned char *htype, unsigned char *dtype0,
		unsigned char *dtype1, int *nframes, int *nviews, int *lc_n, int *GPUID,
		int nsets, int nf, int max_frames, pthread_t thread1, pthread_t thread2,
		cudaStream_t *gpu0_stream, cudaStream_t *gpu1_stream);

__host__ void mkparlist_gpu(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, double *fparstep, double *fpartol,
		double *fparabstol, int *fpartype, double **fpntr,
		int nfpar, int nsets);

__host__ double penalties_gpu(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat);

__host__ void pos2deldop_gpu64(struct par_t *dpar, struct mod_t *dmod, struct
		dat_t *ddat, struct pos_t **pos, struct deldopfrm_t **frame, int4
		*xylim, int *ndel, int *ndop, double orbit_xoff, double orbit_yoff,
		double orbit_dopoff, int body, int set,	int nfrm_alloc,	int v, int
		*badradararr, cudaStream_t *p2d_stream);

__host__ void pos2deldop_MFS_gpu64(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, struct pos_t **pos, struct deldopfrm_t **frame,
		int4 *xylim, int *ndel, int *ndop, double orbit_xoff, double
		orbit_yoff, double orbit_dopoff, int body, int nsets, int v,
		int *badradararr, cudaStream_t *p2d_stream);

__host__ int pos2doppler_gpu64(struct par_t *dpar, struct mod_t *dmod, struct
		dat_t *ddat, struct pos_t **pos, struct dopfrm_t **frame, int4 *xylim,
		double orbit_xoff, double orbit_yoff, double orbit_dopoff, int *ndop,
		int body, int set, int nfrm_alloc, int v, int *badradararr,
		cudaStream_t *pds_stream);

__host__ int posvis_gpu64(struct par_t *dpar, struct mod_t *dmod, struct pos_t
		**pos, struct vertices_t **verts, double3
		orbit_offset, int *posn, int *outbndarr, int set, int nfrm_alloc, int
		src, int nf, int body, int comp, unsigned char type,
		cudaStream_t *pv_stream, int src_override);

__host__ int posvis_MFS_gpu64(struct par_t *dpar, struct mod_t *dmod, struct pos_t
		**pos, struct vertices_t **verts, double3 orbit_offset, int *posn, int
		*outbndarr, int nsets, int nf, int body, int comp, cudaStream_t *pv_stream);

__host__ int read_dat_gpu( struct par_t *par, struct mod_t *mod, struct dat_t *dat);

__host__ int read_dat_gpu_MFS( struct par_t *par, struct mod_t *mod, struct dat_t *dat);

__host__ int read_dat_mgpu( struct par_t *par, struct mod_t *mod,
		struct dat_t *dat, int gpuid);

__host__ void realize_delcor_gpu(struct dat_t *ddat, double delta_delcor0,
		int delcor0_mode, int nsets, int *nframes);

__host__ void realize_delcor_MFS_gpu(struct dat_t *ddat, double delta_delcor0, int delcor0_mode,
		int nsets);

__host__ void realize_delcor_pthreads(struct dat_t *ddat0, struct dat_t *ddat1,
		double delta_delcor0, int delcor0_mode, int nsets, int *nframes, int
		*GPUID, unsigned char *type, pthread_t thread1, pthread_t thread2);

__host__ void realize_dopscale_gpu(struct par_t *dpar, struct dat_t *ddat0,
		double dopscale_factor, int dopscale_mode,	int
		nsets, unsigned char *dtype);

__host__ void realize_dopscale_MFS_gpu(struct par_t *dpar, struct dat_t
		*ddat, double dopscale_factor, int dopscale_mode, int nsets);

__host__ void realize_dopscale_pthreads(struct par_t *dpar0, struct par_t
		*dpar1, struct dat_t *ddat0, struct dat_t *ddat1, double
		dopscale_factor, int dopscale_mode, int nsets, unsigned char *dtype0,
		unsigned char *dtype1, int *GPUID);

__host__ void realize_mod_gpu( struct par_t *dpar, struct mod_t *dmod,
		unsigned char type, int nf, cudaStream_t *rm_streams);

__host__ void realize_mod_pthread(struct par_t *dpar0, struct par_t *dpar1,
		struct mod_t *dmod0, struct mod_t *dmod1, unsigned char type, int nf,
		pthread_t thread1, pthread_t thread2, cudaStream_t *gpu0_stream,
		cudaStream_t *gpu1_stream);

__host__ void realize_photo_gpu( struct par_t *dpar, struct mod_t *dmod,
		double radalb_factor, double optalb_factor, int albedo_mode, int nf);

__host__ void realize_photo_pthread(struct par_t *dpar0, struct par_t *dpar1,
		struct mod_t *dmod0, struct mod_t *dmod1, double radalb_factor,
		double optalb_factor, int albedo_mode, int nf, pthread_t thread1,
		pthread_t thread2);

__host__ void realize_spin_gpu(struct par_t *dpar, struct mod_t *dmod, struct
		dat_t *ddat, unsigned char *htype, int *nframes, int *nviews, int nsets,
		cudaStream_t *rs_stream);

__host__ void realize_spin_MFS_gpu(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, int nsets);

__host__ void realize_spin_pthread(struct par_t *dpar0, struct par_t *dpar1,
		struct mod_t *dmod0, struct mod_t *dmod1, struct dat_t *ddat0, struct
		dat_t *ddat1,unsigned char *htype, int *nframes, int *GPUID, int
		*nviews, int nsets, pthread_t thread1, pthread_t thread2, cudaStream_t
		*gpu0_stream, cudaStream_t *gpu1_stream);

__host__ void realize_xyoff_gpu( struct dat_t *ddat, int nsets,
		unsigned char *dtype);

__host__ void realize_xyoff_pthreads(struct dat_t *ddat0, struct dat_t *ddat1,
		int nsets, unsigned char *dtype0, unsigned char *dtype1, int *GPUID);

__host__ int read_mod_pthread( struct par_t *par, struct mod_t *mod0, struct mod_t *mod1,
		pthread_t thread1, pthread_t thread2);

__host__ void show_deldoplim_gpu(struct dat_t *ddat,
		unsigned char *type, int nsets, int *nframes, int maxframes);

__host__ void show_deldoplim_MFS_gpu(struct dat_t *ddat, int nsets);

__host__ void show_deldoplim_pthread(struct dat_t *ddat0, struct dat_t *ddat1,
		unsigned char *type, int nsets, int *nframes, int maxframes, int *GPUID);

__host__ void vary_params_gpu64(struct par_t *dpar, struct mod_t *dmod, struct
		dat_t *ddat, int action, double *deldop_zmax, double *rad_xsec, double
		*opt_brightness, double *cos_subradarlat, int *hnframes, int *hlc_n,
		int *nviews, struct vertices_t **verts, unsigned char *htype, unsigned
		char *dtype, int nf, int nsets, cudaStream_t *vp_stream, int max_frames);

__host__ void vary_params_MFS_gpu64(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, int action, double *deldop_zmax, double *rad_xsec,
		double *opt_brightness, double *cos_subradarlat, int *nviews, struct
		vertices_t **verts, int nf, int nsets, cudaStream_t *vp_stream);

__host__ int posvis_tiled_gpu64(struct par_t *dpar, struct mod_t *dmod, struct
		dat_t *ddat, struct pos_t **pos, struct vertices_t **verts,	double3
		orbit_offset, int *posn, int *outbndarr, int set, int nfrm_alloc,
		int src, int nf, int body, int comp, unsigned char type, cudaStream_t
		*pv_stream,	int src_override);

__host__ void vary_params_pthreads(struct par_t *dpar0, struct par_t *dpar1,
		struct mod_t *dmod0, struct mod_t *dmod1, struct dat_t *ddat0, struct
		dat_t *ddat1, int action, double *deldop_zmax, double *rad_xsec, double
		*opt_brightness, double *cos_subradarlat, int *hnframes, int *hlc_n,
		int *nviews, int *GPUID, struct vertices_t **verts0, struct vertices_t
		**verts1, unsigned char *htype, unsigned char *dtype0, unsigned char
		*dtype1, int nf, int nsets, int max_frames, pthread_t thread1,pthread_t
		thread2, cudaStream_t *gpu0_stream, cudaStream_t *gpu1_stream);

/* CUDA kernels (not all) */
__global__ void add_offsets_to_euler_krnl(struct mod_t *dmod, struct dat_t
		*ddat, double3 *angle_omega_save, int s);

__global__ void cfs_init_devpar_krnl(struct par_t *dpar);

__global__ void cf_add_fit_store_krnl1(struct dat_t *ddat, float **fit_store,
		int nThreads, int s, int f, unsigned char type);

__global__ void cf_add_fit_store_krnl2(struct dat_t *ddat, int s, int f,
		float *overflow, unsigned char type);

__global__ void cf_init_seen_flags_krnl(struct mod_t *dmod, int nf);

__global__ void cf_finish_fit_store_krnl(struct dat_t *ddat, float **fit_store,
		int s, int f, int nThreads, unsigned char type);

__global__ void cf_set_final_pars_krnl(struct par_t *dpar, struct
		dat_t *ddat);

__global__ void posclr_krnl(struct pos_t **pos, int *posn, int f, int dblflg, int lcflg);

//__global__ void posclr_fast_krnl(struct pos_t **pos, int *posn, int dblflg, int lcflg);

__global__ void posclr_radar_krnl64(struct pos_t **pos, int *posn, int f);

//__global__ void posclr_radar_krnl64mod(struct pos_t **pos, int *posn, int f,
//		int xspan, int4 *xylim, int npixels);

__global__ void posclr_radar_krnl64af(struct pos_t **pos, int *posn,
		int *xspan, int4 *xylim, int *npixels, int factor);

//__global__ void posclr_lc_krnl64(struct pos_t **pos, int *posn, int f);

__global__ void posclr_lc_krnl64af(struct pos_t **pos, int *posn, int *xspan, int4 *xylim, int *npixels, int blocks);

//__global__ void calcfits_posclr_lc_krnl64(struct pos_t **pos, int *posn, int f);

//__global__ void calcfits_posclr_lc_krnl64mod(struct pos_t **pos, int *posn, int f,
//		int xspan, int4 *xylim, int npixels);

__global__ void clrvect_krnl64af(struct dat_t *ddat, int *size, int s, int blocks);

__global__ void clrvect_MFS_krnl(struct dat_t *ddat, int size, int s);

__global__ void posclr_mgpu_krnl(struct pos_t **pos, int *posn, int f, int hf,
		int bdflag);

__global__ void cf_mark_pixels_seen_krnl(struct par_t *dpar, struct mod_t *dmod,
		struct pos_t **pos, int4 *xylim, int npixels, int xspan, int f);

__global__ void posmask_krnl64(struct par_t *dpar,struct pos_t **pos,double3 *so,
		double *pixels_per_km, int *posn, int nThreads, int xspan, int f);

__global__ void realize_angleoff_krnl(struct dat_t *ddat, int gpuid);

__global__ void realize_omegaoff_krnl(struct dat_t *ddat, int gpuid);

__global__ void update_spin_angle_krnl(struct mod_t *dmod, double3 *angle_omega_save);

__global__ void posmask_init_krnl64(struct pos_t **pos, double3 *so,
		double *pixels_per_km, int size);

__global__ void set_ae_oe_MFS_krnl(struct dat_t *ddat, struct pos_t **pos, int nsets);

__device__ void dev_POSrect_gpu64(struct pos_t **pos, int src, double imin_dbl,
		double imax_dbl, double jmin_dbl, double jmax_dbl, double4
		*ijminmax_overall, int frm);

__device__ void dev_POSrect_gpu64_short(double imin_dbl, double imax_dbl, double
		jmin_dbl, double jmax_dbl, double4 *ijminmax_overall, int frm);

__device__ void dev_POSrect_gpu64_shared(double	imin_dbl, double imax_dbl,
		double jmin_dbl, double jmax_dbl, double4 *ijminmax_overall_sh, int n);

__device__ void dev_realize_impulse(struct spin_t spin, double t,
		double t_integrate[], double impulse[][3], int *n_integrate, int s, int f, int k);

__device__ void dev_splint_cfs(double *xa,double *ya,double *y2a,int n,double x,double *y);

__host__ void set_up_pos_pinned(struct par_t *par, struct dat_t *dat);
__host__ void set_up_pos_mgpu(struct par_t *par, struct dat_t *dat, int gpuid);
__host__ void set_up_pos_gpu( struct par_t *par, struct dat_t *dat, int nf);
