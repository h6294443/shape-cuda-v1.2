extern "C" {
#include "../shape/head.h"
}

void checkErrorAfterKernelLaunch(const char *location) {
	cudaError_t cudaStatus;
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Kernel launch failed in %s: %s\n", location, cudaGetErrorString(cudaStatus));
	}
}
void deviceSyncAfterKernelLaunch(const char *location) {
	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaError_t cudaStatus;
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess)
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching the kernel in %s.\n", cudaStatus, location);
}
__host__ void gpuAssert(cudaError_t code, const char *file, int line)
{
	if (code != cudaSuccess)
	{
		fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		exit(code);
	}
}
void pickGPU(int gpuid) {
	cudaSetDeviceFlags(cudaDeviceMapHost);
	cudaError_t cudaStatus;
	cudaStatus = cudaSetDevice(gpuid);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
	}
}

__global__ void zero_fit_overflow_krnl32(struct deldop_t *deldop, int f, int size) {
	/* MAXOVERFLOW^2 - threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int x = offset % MAXOVERFLOW;
	int y = offset / MAXOVERFLOW;

	if (offset < size) {
		deldop->frame[f].fit_overflow32[x][y] = 0.0;
	}
}

__global__ void zero_fit_overflow_krnl64(struct deldop_t *deldop, int f, int size) {
	/* MAXOVERFLOW^2 - threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int x = offset % MAXOVERFLOW;
	int y = offset / MAXOVERFLOW;

	if (offset < size) {
		deldop->frame[f].fit_overflow64[x][y] = 0.0;
	}
}

__host__ void zero_fit_overflow(struct deldop_t *deldop, int f) {
	/* Wrapper function */
	dim3 THD, BLK;
	int threads = 0;
	threads = MAXOVERFLOW*MAXOVERFLOW;
	THD.x = maxThreadsPerBlock;
	BLK.x = floor((THD.x - 1 + threads) / THD.x);

	if (FP64)
		zero_fit_overflow_krnl64<<<BLK,THD>>>(deldop, f, threads);
	else
		zero_fit_overflow_krnl32<<<BLK,THD>>>(deldop, f, threads);
	checkErrorAfterKernelLaunch("zero_fit_overflow_krnl");
}
