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
__global__ void device_reduce_block_atomic_kernel_brt(struct pos_t **pos, double** out,
		int N, int f, int flt);
__global__ void device_reduce_block_atomic_kernel_ddf(struct dat_t *ddat, float** out,
		int N, int f, int s);
__global__ void device_sum_block_atomic_kernel(float *in, float* out, int N);

__global__ void set_idata_pntr_krnl(struct dat_t *ddat, float *d_idata,
		int set, int frm, int size);

/* Specific reduction functions */
__host__ float compute_deldop_xsec_gpu(struct dat_t *ddat, int nframes,
		int size, int set, cudaStream_t *sb_stream);
__host__ float compute_doppler_xsec(struct dat_t *ddat, int ndop,
		int set, int frm);;
__host__ float compute_model_area(struct mod_t *dmod, int c, int size);
__host__ float compute_zmax_gpu(struct dat_t *ddat, struct pos_t **pos,
		int nframes, int size, int set, cudaStream_t *sb_stream);
__host__ void dvdI_reduce_streams(struct mod_t *dmod, float *dv, float *dcom0,
		float *dcom1, float *dcom2, float *dI00, float *dI01, float *dI02,
		float *dI10, float *dI11, float *dI12, float *dI20, float *dI21,
		float *dI22, int size, int c, cudaStream_t *dv_streams);
__host__ double find_max_in_double_array(double *in, int size);
__host__ double find_min_in_double_array(double *in, int size);
__host__ void sum_brightness_streams(struct dat_t *ddat, struct pos_t **pos,
		int nframes, int size, int flt, int set, cudaStream_t *sb_stream);
__host__ double sum_double_array(double *a, int size);
__host__ void sum_2_double_arrays(double *a, double *b, double *absum, int size);



