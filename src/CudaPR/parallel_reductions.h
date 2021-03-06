#include <cuda.h>
#include <cuda_runtime_api.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <math.h>
/* Common reduction functions */

int isPow2(unsigned int x);
unsigned int nextPow2(unsigned int x);

__global__ void deviceReduceWarpAtomicKernel(float *in, float* out, int N);
__global__ void device_reduce_block_atomic_kernel_f(float *in, float* out, int N);

float2 getNumBlocksAndThreads(int n, int maxBlocks, int maxThreads);

/* Specific reduction kernels */
__global__ void device_reduce_block_atomic_kernel_brt(struct pos_t **pos, double *out,
		int N, int f, int flt);
__global__ void device_reduce_block_atomic_kernel_ddf(struct dat_t *ddat, float *out,
		int N, int f, int s);
__global__ void device_sum_block_atomic_kernel(float *in, float* out, int N);

__global__ void set_idata_pntr_krnl(struct dat_t *ddat, float *d_idata,
		int set, int frm, int size);

/* Specific reduction functions */
__host__ double compute_deldop_xsec_gpu(struct dat_t *ddat, int nframes,
		int size, int set, cudaStream_t *sb_stream);

__host__ double compute_deldop_xsec_hyb(struct dat_t *ddat, int nframes,
		int size, int set, cudaStream_t *sb_stream);		
		
__host__ double compute_doppler_xsec(struct dat_t *ddat, int ndop,
		int set, int frm);

__host__ double compute_model_area(struct mod_t *dmod, int c, int size);

__host__ double compute_zmax_gpu(struct dat_t *ddat, struct pos_t **pos,
		int nframes, int size, int set, cudaStream_t *sb_stream);

__host__ double compute_zmax_MFS_gpu(struct dat_t *ddat, struct pos_t **pos,
		int nsets, int pos_size, cudaStream_t *sb_stream);
		
__host__ double compute_zmax_MFS_hyb(struct dat_t *ddat, struct pos_t **pos,
		int nsets, int pos_size, cudaStream_t *sb_stream);

__host__ void dvdI_reduce_gpu(struct mod_t *dmod, double *dv, double *dcom0,
		double *dcom1, double *dcom2, double *dI00, double *dI01, double *dI02,
		double *dI10, double *dI11, double *dI12, double *dI20, double *dI21,
		double *dI22, int size, int c, cudaStream_t *dv_streams);
		
__host__ void dvdI_reduce_hyb(struct mod_t *dmod, double *dv, double *dcom0,
		double *dcom1, double *dcom2, double *dI00, double *dI01, double *dI02,
		double *dI10, double *dI11, double *dI12, double *dI20, double *dI21,
		double *dI22, int size, int c, cudaStream_t *dv_streams);

__host__ double find_max_in_double_array(double *in, int size);

__host__ double find_min_in_double_array(double *in, int size);

__host__ void sum_brightness_gpu(struct dat_t *ddat, struct pos_t **pos,
		int nframes, int size, int flt, int set, int maxthds,
		int4 maxxylim, cudaStream_t *sb_stream);

__host__ double sum_double_array(double *a, int size);

__host__ void sum_2_double_arrays(double *a, double *b, double *absum, int size);

__host__ void sum_o2m2om_gpu(struct dat_t *ddat, double *o2, double *m2, double *om,
		int nframes, int size, int set, cudaStream_t *sb_stream);
		
__host__ void sum_o2m2om_hyb(struct dat_t *ddat, double *o2, double *m2, double *om,
		int nframes, int size, int set, cudaStream_t *sb_stream);
