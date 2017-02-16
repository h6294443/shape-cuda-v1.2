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
#define MAXIMP 10

extern int CUDA;	/* Flag whether to use CUDA or not*/
extern int STREAMS;	/* Flag whether to use streams or not */
extern int GPU;		/* Which GPU will be used */
extern int DYNPROC;  /* Flag wether to use dynamic processing or not */
extern int POSVIS_SEPARATE;	/* Flag wether to use separate xlim/ylim calculation */
extern int AF;		/* Flag wether to process all frames in a set at once */
extern int maxThreadsPerBlock;
extern struct par_t *dev_par;
extern struct mod_t *dev_mod;
extern struct dat_t *dev_dat;
extern float *dfit, *ddfit;		// pointers for doppler fit and delay doppler fit

extern double *fparstep;		// par->fparstep
extern double *fpartol;			// par->partol
extern double *fparabstol;		// par->fparabstol
extern double **fpntr;			// par->pntr
extern int *fpartype;			// par->fpartype

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
__host__ double bestfit_CUDA(struct par_t *dpar, struct mod_t *dmod, struct dat_t
		*ddat, struct par_t *par, struct mod_t *mod, struct dat_t *dat);
__host__ void calc_fits_cuda(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat);
__host__ double chi2_cuda(struct par_t *dpar, struct dat_t *ddat, int list_breakdown);
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
__host__ void dopoffs_cuda(struct dat_t *ddat, int s);
__host__ void dvdI_reduce_single(struct mod_t *dmod, float *dv, float *dcom0,
		float *dcom1, float *dcom2, float *dI00, float *dI01, float *dI02,
		float *dI10, float *dI11, float *dI12, float *dI20, float *dI21,
		float *dI22, int size, int c);
__host__ void gpuAssert(cudaError_t code, const char *file, int line);
__host__ void mkparlist_cuda(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, double *fparstep, double *fpartol,
		double *fparabstol, int *fpartype, double **fpntr);
__host__ double penalties_cuda(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat);
__host__ int pos2deldop_cuda_2(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, double orbit_xoff, double orbit_yoff, double
		orbit_dopoff, int body, int set, int frm, int v);
__host__ int pos2deldop_cuda_af(struct par_t *dpar, struct mod_t *dmod, struct
		dat_t *ddat, double orbit_xoff, double orbit_yoff, double
		orbit_dopoff, int body, int set, int nframes, int v);
__host__ int pos2doppler_cuda_2( struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, double orbit_xoff, double orbit_yoff, double
		orbit_dopoff, int body, int set, int frm, int v);
__host__ int pos2doppler_cuda_af( struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, double orbit_xoff, double orbit_yoff, double
		orbit_dopoff, int body, int set, int nframes, int v);
__host__ int posvis_cuda_2(struct par_t *dpar, struct mod_t *dmod, struct
		dat_t *ddat, double orbit_offset[3], int set, int frame, int src,
		int body, int comp);
__host__ int posvis_af(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, double orbit_offset[3], int set, int nframes,
		int src, int body, int comp);
__host__ void realize_delcor_cuda(struct dat_t *ddat, double delta_delcor0,
		int delcor0_mode, int nsets);
__host__ void realize_dopscale_cuda(struct par_t *dpar, struct dat_t *ddat,
		double dopscale_factor, int dopscale_mode);
__host__ void realize_mod_cuda( struct par_t *dpar, struct mod_t *dmod,
		unsigned char type);
__host__ void realize_photo_cuda( struct par_t *dpar, struct mod_t *dmod,
        double radalb_factor, double optalb_factor, int albedo_mode);
__host__ void realize_spin_cuda( struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, int nsets);
__host__ void realize_xyoff_cuda( struct dat_t *ddat);
__host__ void show_deldoplim_cuda(struct dat_t *dat, struct dat_t *ddat);
__host__ void vary_params_cuda(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, int action, double *deldop_zmax,
		double *rad_xsec, double *opt_brightness, double *cos_subradarlat,
		int nsets);
__host__ void vary_params_af( struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, int action, double *deldop_zmax, double
		*rad_xsec, double *opt_brightness, double *cos_subradarlat, int nsets);

__device__ int cubic_realroots_cuda( double *coeff, double *realroot);
__device__ void dev_bsstep(double *y, double *dydx, int nv, double *xx, double htry, double eps,
		double *yscal, double *hdid, double *hnext, void (*derivs)(double,double *,double *));
__device__ double dev_cel(double qqc, double pp, double aa, double bb);
__device__ void dev_cotrans1( double y[3], double *a, double x[3], int dir);
__device__ void dev_cotrans2( double y[3], double a[3][3], double x[3], int dir);
__device__ void dev_cotrans3(double y[3], double a[3][3], double x[3],
		int dir);
__device__ void dev_cotrans4(float3 *y, double a[3][3], double x[3], int dir, int f);
__device__ double dev_cross( double z[3], double x[3], double y[3]);
__device__ double dev_dot( double x[3], double y[3]);
__device__ void dev_euler2mat( double m[3][3], double phi, double theta, double psi);
__device__ void dev_facmom( double fv0[3], double fv1[3], double fv2[3], double fn[3],
        double *dv, double dvr[3], double dI[3][3]);
__device__ double dev_facnrm( struct vertices_t verts, int fi);
__device__ int dev_gamma_trans(float *datum, double gamma);
__device__ double dev_gammln(double xx);
__device__ double dev_hapke( double cosi, double cose, double phase,
        double w, double h, double B0, double g, double theta);
__device__ void dev_inteuler( struct spin_t spin, double t[], double impulse[][3], int n,
		double w[3], double m[3][3], unsigned char pa, unsigned char method, double int_abstol);
__device__ void dev_mat2euler( double m[3][3], double *phi, double *theta, double *psi);
__device__ void dev_mmid( double *y, double *dydx, int nvar1, double xs, double htot,
		int nstep, double *yout, void (*dev_derivs)( double, double *, double *));
__device__ void dev_mmmul(double x[3][3], double y[3][3], double z[3][3]);
__device__ void dev_mtrnsps( double a[3][3], double b[3][3]);
__device__ double dev_normalize(double *u);
__device__ void dev_odeint( double *ystart, int nvar, double x1, double x2, double eps,
	double h1, double hmin, int *nok, int *nbad, void (*derivs)(double,double *,double *),
	void (*drkqc)(double *,double *,int,double *,double,double,double
			*,double *,double *,void (*)(double,double *,double *)));
__device__ double dev_plgndr(int l,int m,double x);
__device__ void dev_POSrect2(struct pos_t *pos, int src, float imin_dbl,
		float imax_dbl, float jmin_dbl, float jmax_dbl);
__device__ double dev_radlaw( struct photo_t *photo, int ilaw, double cosinc, int c, int f);
__device__ void dev_rzextr( int iest, double xest, double *yest, double *yz, double *dy);

__device__ double radlaw_cuda(union radscat_t *radar, unsigned char *radtype,
		int ilaw, double cosinc, int c, int f);

__global__ void clrvect_krnl(struct dat_t *ddat, int s, int f, int nThreads);
__global__ void euler2mat_krnl( double m[3][3], double phi, double theta, double psi);
__global__ void euler2mat_realize_mod_krnl(struct mod_t *dmod);



__host__ void dbg_print_fit(struct dat_t *ddat, int s, int f);
void dbg_print_fit_host(struct dat_t *ddat, int s, int f);
__host__ void dbg_print_deldop_fit(struct dat_t *ddat, int s, int f);
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
__host__ void dbg_print_pos_z(struct dat_t *ddat, int set, int frm, int n);
__host__ void dbg_print_pos_cose_s(struct dat_t *ddat, int set, int frm, int n);
__host__ void dbg_print_pos_z_af(struct dat_t *ddat, int set, int n);
__host__ void dbg_print_cose_af(struct dat_t *ddat, int set, int n);
