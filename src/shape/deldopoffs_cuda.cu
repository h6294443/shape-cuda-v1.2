/*****************************************************************************************
                                                                             deldopoffs.c

Takes the delay-correction polynomial for a delay-doppler set and figures out the COM
delay and doppler corrections (in units of image rows and columns) for each frame.

Modified 2015 June 3 by CM:
    Implement smearing for the "fit" and "write" actions

Modified 2006 June 21 by CM:
    Changed delres to del_per_pixel and dopres to dop_per_pixel
 *****************************************************************************************/
extern "C" {
#include "head.h"
}

__device__ int deldopoffs_nframes;

__global__ void deldopoffs_get_nframes_krnl(struct dat_t *ddat, int s) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0)
		deldopoffs_nframes = ddat->set[s].desc.deldop.nframes;
}
__global__ void deldopoffs_cuda_krnl(struct dat_t *ddat, int s) {
	/* nframes-threaded kernel */
	int f = blockIdx.x * blockDim.x + threadIdx.x;
	int k, n;
	double del, dop, arg, x;

	if (f < deldopoffs_nframes) {
		for (k=0; k<ddat->set[s].desc.deldop.nviews; k++) {
			x = 1.0;
			dop = 0.0;
			del = ddat->set[s].desc.deldop.delcor.a[0].val;
			arg = ddat->set[s].desc.deldop.frame[f].view[k].t -
					ddat->set[s].desc.deldop.delcor.t0;

			for (n=1; n<=ddat->set[s].desc.deldop.delcor.n; n++) {
				dop += n*ddat->set[s].desc.deldop.delcor.a[n].val*x;
				del +=   ddat->set[s].desc.deldop.delcor.a[n].val*(x*=arg);
			}

			/* del has units of usec */
			ddat->set[s].desc.deldop.frame[f].view[k].deloff =
					del/ddat->set[s].desc.deldop.del_per_pixel;

			/* dop has units of usec/day and there are 86400 sec/day */
			ddat->set[s].desc.deldop.frame[f].view[k].dopoff =
					-dop*ddat->set[s].desc.deldop.Ftx
					/ (ddat->set[s].desc.deldop.dop_per_pixel*86400.0);
		}
	}
}
__global__ void deldopoffs_cuda_f_krnl(struct dat_t *ddat, int s, int nframes) {
	/* nframes-threaded kernel */
	int f = blockIdx.x * blockDim.x + threadIdx.x;
	int k, n;
	float del, dop, arg, x;

	if (f < nframes) {
		for (k=0; k<ddat->set[s].desc.deldop.nviews; k++) {
			x = 1.0;
			dop = 0.0;
			del = __double2float_rn(ddat->set[s].desc.deldop.delcor.a[0].val);
			arg = __double2float_rn(ddat->set[s].desc.deldop.frame[f].view[k].t) -
					__double2float_rn(ddat->set[s].desc.deldop.delcor.t0);

			for (n=1; n<=ddat->set[s].desc.deldop.delcor.n; n++) {
				dop += n*(__double2float_rn(ddat->set[s].desc.deldop.delcor.a[n].val))*x;
				del +=   __double2float_rn(ddat->set[s].desc.deldop.delcor.a[n].val)*(x*=arg);
			}

			/* del has units of usec */
			ddat->set[s].desc.deldop.frame[f].view[k].deloff =
					del/ddat->set[s].desc.deldop.del_per_pixel;

			/* dop has units of usec/day and there are 86400 sec/day */
			ddat->set[s].desc.deldop.frame[f].view[k].dopoff =
					-dop*ddat->set[s].desc.deldop.Ftx
					/ (ddat->set[s].desc.deldop.dop_per_pixel*86400.0);
		}
	}
}
__host__ void deldopoffs_cuda(struct dat_t *ddat, int s)
{
	int nframes = 0;
	dim3 BLK,THD;

	/* Get # of frames first   */
	deldopoffs_get_nframes_krnl<<<1,1>>>(ddat, s);
	checkErrorAfterKernelLaunch("deldopoffs_get_frames_krnl (deldopoffs_cuda)");
	gpuErrchk(cudaMemcpyFromSymbol(&nframes, deldopoffs_nframes, sizeof(int),
			0, cudaMemcpyDeviceToHost));

	/* Launch nframes-threaded kernel */
	THD.x = nframes;
	deldopoffs_cuda_krnl<<<BLK,THD>>>(ddat, s);
	checkErrorAfterKernelLaunch("deldopoffs_cuda_krnl (deldopoffs_cuda)");
}
__host__ void deldopoffs_cuda_f(struct dat_t *ddat, int s, int nframes)
{
	dim3 BLK,THD;
	THD.x = maxThreadsPerBlock;
	BLK.x = floor((THD.x - 1 + nframes)/THD.x);

	/* Launch nframes-threaded kernel */
	deldopoffs_cuda_f_krnl<<<BLK,THD>>>(ddat, s, nframes);
	checkErrorAfterKernelLaunch("deldopoffs_cuda_f_krnl (deldopoffs_cuda_f)");
}
