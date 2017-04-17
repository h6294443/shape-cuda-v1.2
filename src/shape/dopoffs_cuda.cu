/*****************************************************************************************
                                                                                dopoffs.c

Takes the delay-correction polynomial for a Doppler dataset and figures out the COM
Doppler corrections (in units of Doppler bins) for each frame.

Modified 2015 June 3 by CM:
    Implement smearing for the "fit" and "write" actions

Modified 2006 June 21 by CM:
    Changed dopres to dop_per_bin

Modified 2003 April 26 by CM:
    Removed delay computation
*****************************************************************************************/
extern "C" {
#include "head.h"
}

__device__ int dopoffs_nframes;

__global__ void dopoffs_get_frames_krnl(struct dat_t *ddat, int s) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0)
		dopoffs_nframes = ddat->set[s].desc.doppler.nframes;
}
__global__ void dopoffs_cuda_krnl(struct dat_t *ddat, int s) {
	/* nframes-threaded kernel */
	int f = blockIdx.x * blockDim.x + threadIdx.x;
	int k, n;
	double dop, arg, x;

	if (f < dopoffs_nframes) {
		for (k=0; k<ddat->set[s].desc.doppler.nviews; k++) {
			x = 1.0;
			dop = 0.0;
			arg = ddat->set[s].desc.doppler.frame[f].view[k].t -
					ddat->set[s].desc.doppler.delcor.t0;
			for (n=1; n<=ddat->set[s].desc.doppler.delcor.n; n++) {
				dop += n*ddat->set[s].desc.doppler.delcor.a[n].val*x;
				x *= arg;
			}

			/*  dop has units of usec/day and there are 86400 sec/day  */
			ddat->set[s].desc.doppler.frame[f].view[k].dopoff =
					-dop*ddat->set[s].desc.doppler.Ftx
					/ (ddat->set[s].desc.doppler.dop_per_bin*86400.0);
		}
	}
}
__global__ void dopoffs_cuda_f_krnl(struct dat_t *ddat, int s, int nframes) {
	/* nframes-threaded kernel */
	int f = blockIdx.x * blockDim.x + threadIdx.x;
	int k, n;
	float dop, arg, x;

	if (f < nframes) {
		for (k=0; k<ddat->set[s].desc.doppler.nviews; k++) {
			x = 1.0;
			dop = 0.0;
			arg = __double2float_rn(ddat->set[s].desc.doppler.frame[f].view[k].t) -
					__double2float_rn(ddat->set[s].desc.doppler.delcor.t0);
			for (n=1; n<=ddat->set[s].desc.doppler.delcor.n; n++) {
				dop += n*(__double2float_rn(ddat->set[s].desc.doppler.delcor.a[n].val))*x;
				x *= arg;
			}

			/*  dop has units of usec/day and there are 86400 sec/day  */
			ddat->set[s].desc.doppler.frame[f].view[k].dopoff =
					-dop*ddat->set[s].desc.doppler.Ftx
					/ (ddat->set[s].desc.doppler.dop_per_bin*86400.0);
		}
	}

}
__host__ void dopoffs_cuda(struct dat_t *ddat, int s)
{
  int nframes = 0;
  dim3 BLK,THD;
  //ddat->set[s].desc.doppler

  /* Get frames and views */
  dopoffs_get_frames_krnl<<<1,1>>>(ddat, s);
  checkErrorAfterKernelLaunch("dopoffs_get_frames_krnl (dopoffs_cuda)");
  gpuErrchk(cudaMemcpyFromSymbol(&nframes, dopoffs_nframes, sizeof(int),
		  0, cudaMemcpyDeviceToHost));

  /* Launch nframes-threaded kernel */
  THD.x = nframes;
  dopoffs_cuda_krnl<<<BLK,THD>>>(ddat, s);
  checkErrorAfterKernelLaunch("dopoffs_cuda_krnl (dopoffs_cuda)");
}
__host__ void dopoffs_cuda_f(struct dat_t *ddat, int s, int nframes)
{
  dim3 BLK,THD;
  THD.x = maxThreadsPerBlock;
  BLK.x = floor((THD.x - 1 + nframes)/THD.x);

  /* Launch nframes-threaded kernel */
  dopoffs_cuda_f_krnl<<<BLK,THD>>>(ddat, s, nframes);
  checkErrorAfterKernelLaunch("dopoffs_cuda_f_krnl (dopoffs_cuda_f)");
}
