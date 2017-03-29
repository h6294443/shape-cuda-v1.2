extern "C" {
#include "head.h"
}
#include <stdio.h>

#define MAX_BLOCK_DIM_SIZE 65535
#ifndef MIN
#define MIN(x,y) ((x < y) ? x : y)
#endif
//#ifndef _REDUCE_KERNEL_H_
//#define _REDUCE_KERNEL_H_

__device__ float reduction_sum_rad_xsec = 0.0;	// Used for the deldop xsec all frames fx
int isPow2(unsigned int x)
{
    return ((x&(x-1))==0);
}

unsigned int nextPow2(unsigned int x)
{
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return ++x;
}

__inline__ __device__
float warpReduceSum(float val) {
  for (int offset = warpSize/2; offset > 0; offset /= 2)
    val += __shfl_down(val, offset);
  return val;
}
__inline__ __device__
float warpReduceMax(float val) {
  for (int offset = warpSize/2; offset > 0; offset /= 2)
    val = fmaxf((__shfl_down(val, offset)),val);
  return val;
}
__inline__ __device__
float blockReduceSum(float val) {
  static __shared__ float shared[32];
  int lane=threadIdx.x%warpSize;
  int wid=threadIdx.x/warpSize;
  val=warpReduceSum(val);

  //write reduced value to shared memory
  if(lane==0) shared[wid]=val;
  __syncthreads();

  //ensure we only grab a value from shared memory if that warp existed
  val = (threadIdx.x<blockDim.x/warpSize) ? shared[lane] : float(0.0);
  if(wid==0) val=warpReduceSum(val);

  return val;
}
__inline__ __device__
float blockReduceMax(float val) {
  static __shared__ float shared[32];
  int lane=threadIdx.x%warpSize;
  int wid=threadIdx.x/warpSize;
  val=warpReduceMax(val);

  //write reduced value to shared memory
  if(lane==0) shared[wid]=val;
  __syncthreads();

  //ensure we only grab a value from shared memory if that warp existed
  val = (threadIdx.x<blockDim.x/warpSize) ? shared[lane] : float(0.0);
  if(wid==0) val=warpReduceMax(val);

  return val;
}
__global__ void deviceReduceWarpAtomicKernel(float *in, float* out, int N) {
	float sum = float(0.0);
  for(int i = blockIdx.x * blockDim.x + threadIdx.x;
      i < N;
      i += blockDim.x * gridDim.x) {
    sum += in[i];
  }
  sum = warpReduceSum(sum);
  if (threadIdx.x & (warpSize - 1) == 0)
    atomicAdd(out, sum);
}
__global__ void device_reduce_block_atomic_kernel(float *in, float* out, int N) {
  float sum=float(0.0);
  for(int i=blockIdx.x*blockDim.x+threadIdx.x;i<N;i+=blockDim.x*gridDim.x) {
    sum+=in[i];
  }
  sum=blockReduceSum(sum);
  if(threadIdx.x==0)
    atomicAdd(out,sum);
}
__global__ void device_sum_block_atomic_kernel(float *in, float* out, int N) {
	float maxz=float(0.0);
	for(int i=blockIdx.x*blockDim.x+threadIdx.x;i<N;i+=blockDim.x*gridDim.x) {
		maxz=fmaxf(in[i],maxz);
	}
	maxz=blockReduceMax(maxz);
	if(threadIdx.x==0)
		atomicExch(out,maxz);
}

/* Compute the number of threads and blocks to use for reduction kernel 6
 * For kernel 6, we observe the maximum specified number of blocks, because
 * each thread in that kernel can process a variable number of elements. */
float2 getNumBlocksAndThreads(int n, int maxBlocks, int maxThreads)
{
	//get device capability, to avoid block/grid size exceed the upper bound
	cudaDeviceProp prop;
	int device, threads, blocks;
	float2 xb_yt;
	gpuErrchk(cudaGetDevice(&device));
	gpuErrchk(cudaGetDeviceProperties(&prop, device));

	threads = (n < maxThreads*2) ? nextPow2((n + 1)/ 2) : maxThreads;
	blocks = (n + (threads * 2 - 1)) / (threads * 2);

	if ((float)threads*blocks > (float)prop.maxGridSize[0] * prop.maxThreadsPerBlock)
		printf("Array size for parallel reduction is too large!\n");

	if (blocks > prop.maxGridSize[0])
	{
		printf("Grid size <%d> exceeds the device capability <%d>, set block size as %d (original %d)\n",
				blocks, prop.maxGridSize[0], threads*2, threads);

		blocks /= 2;
		threads *= 2;
	}
	blocks = MIN(maxBlocks, blocks);
	xb_yt.x = blocks;
	xb_yt.y = threads;
	return xb_yt;
}

// Utility class used to avoid linker errors with extern
// unsized shared memory arrays with templated type
template<class T>
struct SharedMemory
{
    __device__ inline operator       T *()
    {
        extern __shared__ int __smem[];
        return (T *)__smem;
    }

    __device__ inline operator const T *() const
    {
        extern __shared__ int __smem[];
        return (T *)__smem;
    }
};

// specialize for double to avoid unaligned memory
// access compile errors
template<>
struct SharedMemory<double>
{
    __device__ inline operator       double *()
    {
        extern __shared__ double __smem_d[];
        return (double *)__smem_d;
    }

    __device__ inline operator const double *() const
    {
        extern __shared__ double __smem_d[];
        return (double *)__smem_d;
    }
};

/*
    Parallel sum reduction using shared memory
    - takes log(n) steps for n input elements
    - uses n threads
    - only works for power-of-2 arrays
*/

/*
    This version adds multiple elements per thread sequentially.  This reduces the overall
    cost of the algorithm while keeping the work complexity O(n) and the step complexity O(log n).
    (Brent's Theorem optimization)

    Note, this kernel needs a minimum of 64*sizeof(T) bytes of shared memory.
    In other words if blockSize <= 32, allocate 64*sizeof(T) bytes.
    If blockSize > 32, allocate blockSize*sizeof(T) bytes.
*/
template <class T, unsigned int blockSize, int nIsPow2>
__global__ void
reduce6(T *g_idata, T *g_odata, unsigned int n)
{
    T *sdata = SharedMemory<T>();

    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*blockSize*2 + threadIdx.x;
    unsigned int gridSize = blockSize*2*gridDim.x;

    T mySum = 0;

    // we reduce multiple elements per thread.  The number is determined by the
    // number of active thread blocks (via gridDim).  More blocks will result
    // in a larger gridSize and therefore fewer elements per thread
    while (i < n)
    {
        mySum += g_idata[i];

        // ensure we don't read out of bounds -- this is optimized away for powerOf2 sized arrays
        if (nIsPow2 || i + blockSize < n)
            mySum += g_idata[i+blockSize];

        i += gridSize;
    }

    // each thread puts its local sum into shared memory
    sdata[tid] = mySum;
    __syncthreads();


    // do reduction in shared mem
    if ((blockSize >= 512) && (tid < 256))
    {
        sdata[tid] = mySum = mySum + sdata[tid + 256];
    }

    __syncthreads();

    if ((blockSize >= 256) &&(tid < 128))
    {
            sdata[tid] = mySum = mySum + sdata[tid + 128];
    }

     __syncthreads();

    if ((blockSize >= 128) && (tid <  64))
    {
       sdata[tid] = mySum = mySum + sdata[tid +  64];
    }

    __syncthreads();

#if (__CUDA_ARCH__ >= 300 )
    if ( tid < 32 )
    {
        // Fetch final intermediate sum from 2nd warp
        if (blockSize >=  64) mySum += sdata[tid + 32];
        // Reduce final warp using shuffle
        for (int offset = warpSize/2; offset > 0; offset /= 2)
        {
            mySum += __shfl_down(mySum, offset);
        }
    }
#else
    // fully unroll reduction within a single warp
    if ((blockSize >=  64) && (tid < 32))
    {
        sdata[tid] = mySum = mySum + sdata[tid + 32];
    }

    __syncthreads();

    if ((blockSize >=  32) && (tid < 16))
    {
        sdata[tid] = mySum = mySum + sdata[tid + 16];
    }

    __syncthreads();

    if ((blockSize >=  16) && (tid <  8))
    {
        sdata[tid] = mySum = mySum + sdata[tid +  8];
    }

    __syncthreads();

    if ((blockSize >=   8) && (tid <  4))
    {
        sdata[tid] = mySum = mySum + sdata[tid +  4];
    }

    __syncthreads();

    if ((blockSize >=   4) && (tid <  2))
    {
        sdata[tid] = mySum = mySum + sdata[tid +  2];
    }

    __syncthreads();

    if ((blockSize >=   2) && ( tid <  1))
    {
        sdata[tid] = mySum = mySum + sdata[tid +  1];
    }

    __syncthreads();
#endif

    // write result for this block to global mem
    if (tid == 0) g_odata[blockIdx.x] = mySum;
}

/*
    This version adds multiple elements per thread sequentially.  This reduces the overall
    cost of the algorithm while keeping the work complexity O(n) and the step complexity O(log n).
    (Brent's Theorem optimization)

    Note, this kernel needs a minimum of 64*sizeof(T) bytes of shared memory.
    In other words if blockSize <= 32, allocate 64*sizeof(T) bytes.
    If blockSize > 32, allocate blockSize*sizeof(T) bytes.
*/
template <class T, unsigned int blockSize, int nIsPow2>
__global__ void
maxz6(T *g_idata, T *g_odata, unsigned int n)
{
    T *sdata = SharedMemory<T>();

    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*blockSize*2 + threadIdx.x;
    unsigned int gridSize = blockSize*2*gridDim.x;

    T myMax = 0;

    // we reduce multiple elements per thread.  The number is determined by the
    // number of active thread blocks (via gridDim).  More blocks will result
    // in a larger gridSize and therefore fewer elements per thread
    while (i < n)
    {
        myMax = fmaxf(g_idata[i], myMax);

        // ensure we don't read out of bounds -- this is optimized away for powerOf2 sized arrays
        if (nIsPow2 || i + blockSize < n)
            myMax = fmaxf(g_idata[i+blockSize], myMax);

        i += gridSize;
    }

    // each thread puts its local sum into shared memory
    sdata[tid] = myMax;
    __syncthreads();


    // do reduction in shared mem
    if ((blockSize >= 512) && (tid < 256))
    {
        sdata[tid] = myMax = fmaxf(myMax, sdata[tid + 256]);
    }

    __syncthreads();

    if ((blockSize >= 256) &&(tid < 128))
    {
        sdata[tid] = myMax = fmaxf(myMax, sdata[tid + 128]);
    }

     __syncthreads();

    if ((blockSize >= 128) && (tid <  64))
    {
       sdata[tid] = myMax = fmaxf(myMax, sdata[tid +  64]);
    }

    __syncthreads();

#if (__CUDA_ARCH__ >= 300 )
    if ( tid < 32 )
    {
        // Fetch final intermediate sum from 2nd warp
        if (blockSize >=  64) myMax = fmaxf(sdata[tid + 32], myMax);
        // Reduce final warp using shuffle
        for (int offset = warpSize/2; offset > 0; offset /= 2)
        {
            myMax = fmaxf(myMax, __shfl_down(myMax, offset));
        }
    }
#else
    // fully unroll reduction within a single warp
    if ((blockSize >=  64) && (tid < 32))
    {
        sdata[tid] = myMax = fmaxf(myMax, sdata[tid + 32]);
    }

    __syncthreads();

    if ((blockSize >=  32) && (tid < 16))
    {
        sdata[tid] = myMax = fmaxf(myMax, sdata[tid + 16]);
    }

    __syncthreads();

    if ((blockSize >=  16) && (tid <  8))
    {
        sdata[tid] = myMax = fmaxf(myMax, sdata[tid +  8]);
    }

    __syncthreads();

    if ((blockSize >=   8) && (tid <  4))
    {
        sdata[tid] = myMax = fmaxf(myMax, sdata[tid +  4]);
    }

    __syncthreads();

    if ((blockSize >=   4) && (tid <  2))
    {
        sdata[tid] = myMax = fmaxf(myMax, sdata[tid +  2]);
    }

    __syncthreads();

    if ((blockSize >=   2) && ( tid <  1))
    {
        sdata[tid] = myMax = fmaxf(myMax, sdata[tid +  1]);
    }

    __syncthreads();
#endif

    // write result for this block to global mem
    if (tid == 0)
    	g_odata[blockIdx.x] = myMax;
}
template <class T, unsigned int blockSize, int nIsPow2>
__global__ void
minz6(T *g_idata, T *g_odata, unsigned int n)
{
    T *sdata = SharedMemory<T>();

    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*blockSize*2 + threadIdx.x;
    unsigned int gridSize = blockSize*2*gridDim.x;

    T myMin = 0;

    // we reduce multiple elements per thread.  The number is determined by the
    // number of active thread blocks (via gridDim).  More blocks will result
    // in a larger gridSize and therefore fewer elements per thread
    while (i < n)
    {
        myMin = fminf(g_idata[i], myMin);

        // ensure we don't read out of bounds -- this is optimized away for powerOf2 sized arrays
        if (nIsPow2 || i + blockSize < n)
            myMin = fminf(g_idata[i+blockSize], myMin);

        i += gridSize;
    }

    // each thread puts its local sum into shared memory
    sdata[tid] = myMin;
    __syncthreads();


    // do reduction in shared mem
    if ((blockSize >= 512) && (tid < 256))
    {
        sdata[tid] = myMin = fminf(myMin, sdata[tid + 256]);
    }

    __syncthreads();

    if ((blockSize >= 256) &&(tid < 128))
    {
            sdata[tid] = myMin = fminf(myMin, sdata[tid + 128]);
    }

     __syncthreads();

    if ((blockSize >= 128) && (tid <  64))
    {
       sdata[tid] = myMin = fminf(myMin, sdata[tid +  64]);
    }

    __syncthreads();

#if (__CUDA_ARCH__ >= 300 )
    if ( tid < 32 )
    {
        // Fetch final intermediate sum from 2nd warp
        if (blockSize >=  64) myMin = fminf(sdata[tid + 32], myMin);
        // Reduce final warp using shuffle
        for (int offset = warpSize/2; offset > 0; offset /= 2)
        {
            myMin = fminf(myMin, __shfl_down(myMin, offset));
        }
    }
#else
    // fully unroll reduction within a single warp
    if ((blockSize >=  64) && (tid < 32))
    {
        sdata[tid] = myMin = fminf(myMin, sdata[tid + 32]);
    }

    __syncthreads();

    if ((blockSize >=  32) && (tid < 16))
    {
        sdata[tid] = myMin = fminf(myMin, sdata[tid + 16]);
    }

    __syncthreads();

    if ((blockSize >=  16) && (tid <  8))
    {
        sdata[tid] = myMin = fminf(myMin, sdata[tid +  8]);
    }

    __syncthreads();

    if ((blockSize >=   8) && (tid <  4))
    {
        sdata[tid] = myMin = fminf(myMin, sdata[tid +  4]);
    }

    __syncthreads();

    if ((blockSize >=   4) && (tid <  2))
    {
        sdata[tid] = myMin = fminf(myMin, sdata[tid +  2]);
    }

    __syncthreads();

    if ((blockSize >=   2) && ( tid <  1))
    {
        sdata[tid] = myMin = fminf(myMin, sdata[tid +  1]);
    }

    __syncthreads();
#endif

    // write result for this block to global mem
    if (tid == 0)
    	g_odata[blockIdx.x] = myMin;
}

////////////////////////////////////////////////////////////////////////////////
// Wrapper function for kernel launch
////////////////////////////////////////////////////////////////////////////////
template <class T>
void
reduce(int size, int threads, int blocks,
       int whichKernel, T *d_idata, T *d_odata)
{
    dim3 dimBlock(threads, 1, 1);
    dim3 dimGrid(blocks, 1, 1);

    // when there is only one warp per block, we need to allocate two warps
    // worth of shared memory so that we don't index shared memory out of bounds
    int smemSize = (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);

    // choose which of the optimized versions of reduction to launch
    switch (whichKernel)
    {
        case 6:
        default:
            if (isPow2(size))
            {
                switch (threads)
                {
                    case 512:
                        reduce6<T, 512, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
                        break;

                    case 256:
                        reduce6<T, 256, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
                        break;

                    case 128:
                        reduce6<T, 128, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
                        break;

                    case 64:
                        reduce6<T,  64, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
                        break;

                    case 32:
                        reduce6<T,  32, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
                        break;

                    case 16:
                        reduce6<T,  16, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
                        break;

                    case  8:
                        reduce6<T,   8, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
                        break;

                    case  4:
                        reduce6<T,   4, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
                        break;

                    case  2:
                        reduce6<T,   2, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
                        break;

                    case  1:
                        reduce6<T,   1, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
                        break;
                }
            }
            else
            {
                switch (threads)
                {
                    case 512:
                        reduce6<T, 512, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
                        break;

                    case 256:
                        reduce6<T, 256, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
                        break;

                    case 128:
                        reduce6<T, 128, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
                        break;

                    case 64:
                        reduce6<T,  64, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
                        break;

                    case 32:
                        reduce6<T,  32, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
                        break;

                    case 16:
                        reduce6<T,  16, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
                        break;

                    case  8:
                        reduce6<T,   8, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
                        break;

                    case  4:
                        reduce6<T,   4, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
                        break;

                    case  2:
                        reduce6<T,   2, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
                        break;

                    case  1:
                        reduce6<T,   1, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
                        break;
                }
            }
            break;
    }
}

// Instantiate the reduction function for 3 types
template void
reduce<int>(int size, int threads, int blocks,
            int whichKernel, int *d_idata, int *d_odata);

template void
reduce<float>(int size, int threads, int blocks,
              int whichKernel, float *d_idata, float *d_odata);

template void
reduce<double>(int size, int threads, int blocks,
               int whichKernel, double *d_idata, double *d_odata);

template <class T>
void
reducedI(int size, int threads, int blocks, int whichKernel, T *dv,
		T *dcom0, T *dcom1, T *dcom2, T *dI00, T *dI01, T *dI02, T *dI10,
		T *dI11, T *dI12, T *dI20, T *dI21, T *dI22, T *d_odata_dv,
		T *d_odata_dcom0, T *d_odata_dcom1, T *d_odata_dcom2, T *d_odata_dI00,
		T *d_odata_dI01, T *d_odata_dI02, T *d_odata_dI10, T *d_odata_dI11,
		T *d_odata_dI12, T *d_odata_dI20, T *d_odata_dI21, T *d_odata_dI22)
{
    dim3 dimBlock(threads, 1, 1);
    dim3 dimGrid(blocks, 1, 1);

    // when there is only one warp per block, we need to allocate two warps
    // worth of shared memory so that we don't index shared memory out of bounds
    int smemSize = (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);

    /* Create 13 streams */
    cudaStream_t stream01, stream02, stream03, stream04, stream05, stream06,
    	stream07, stream08, stream09, stream10, stream11, stream12, stream13;
    cudaStreamCreate(&stream01);	cudaStreamCreate(&stream02);
    cudaStreamCreate(&stream03);	cudaStreamCreate(&stream04);
    cudaStreamCreate(&stream05);	cudaStreamCreate(&stream06);
    cudaStreamCreate(&stream07);	cudaStreamCreate(&stream08);
    cudaStreamCreate(&stream09);	cudaStreamCreate(&stream10);
    cudaStreamCreate(&stream11);	cudaStreamCreate(&stream12);
    cudaStreamCreate(&stream13);
    // choose which of the optimized versions of reduction to launch
    switch (whichKernel)
    {
    case 6:
    default:
    	if (isPow2(size))
    	{
    		switch (threads)
    		{
    		case 512:
    			reduce6<T, 512, true><<< dimGrid, dimBlock, smemSize, stream01 >>>(dv, d_odata_dv, size);
    			reduce6<T, 512, true><<< dimGrid, dimBlock, smemSize, stream02 >>>(dcom0, d_odata_dcom0, size);
    			reduce6<T, 512, true><<< dimGrid, dimBlock, smemSize, stream03 >>>(dcom1, d_odata_dcom1, size);
    			reduce6<T, 512, true><<< dimGrid, dimBlock, smemSize, stream04 >>>(dcom2, d_odata_dcom2, size);
    			reduce6<T, 512, true><<< dimGrid, dimBlock, smemSize, stream05 >>>(dI00, d_odata_dI00, size);
    			reduce6<T, 512, true><<< dimGrid, dimBlock, smemSize, stream06 >>>(dI01, d_odata_dI01, size);
    			reduce6<T, 512, true><<< dimGrid, dimBlock, smemSize, stream07 >>>(dI02, d_odata_dI02, size);
    			reduce6<T, 512, true><<< dimGrid, dimBlock, smemSize, stream08 >>>(dI10, d_odata_dI10, size);
    			reduce6<T, 512, true><<< dimGrid, dimBlock, smemSize, stream09 >>>(dI11, d_odata_dI11, size);
    			reduce6<T, 512, true><<< dimGrid, dimBlock, smemSize, stream10 >>>(dI12, d_odata_dI12, size);
    			reduce6<T, 512, true><<< dimGrid, dimBlock, smemSize, stream11 >>>(dI20, d_odata_dI20, size);
    			reduce6<T, 512, true><<< dimGrid, dimBlock, smemSize, stream12 >>>(dI21, d_odata_dI21, size);
    			reduce6<T, 512, true><<< dimGrid, dimBlock, smemSize, stream13 >>>(dI22, d_odata_dI22, size);
    			break;

    		case 256:
    			reduce6<T, 256, true><<< dimGrid, dimBlock, smemSize, stream01 >>>(dv, d_odata_dv, size);
    			checkErrorAfterKernelLaunch("reduce6 in reducedI, dv (pow2) and 256 threads");
    			reduce6<T, 256, true><<< dimGrid, dimBlock, smemSize, stream02 >>>(dcom0, d_odata_dcom0, size);
    			checkErrorAfterKernelLaunch("reduce6 in reducedI, dcom0 (pow2) and 256 threads");
    			reduce6<T, 256, true><<< dimGrid, dimBlock, smemSize, stream03 >>>(dcom1, d_odata_dcom1, size);
    			checkErrorAfterKernelLaunch("reduce6 in reducedI, dcom1 (pow2) and 256 threads");
    			reduce6<T, 256, true><<< dimGrid, dimBlock, smemSize, stream04 >>>(dcom2, d_odata_dcom2, size);
    			checkErrorAfterKernelLaunch("reduce6 in reducedI, dcom2 (pow2) and 256 threads");
    			reduce6<T, 256, true><<< dimGrid, dimBlock, smemSize, stream05 >>>(dI00, d_odata_dI00, size);
    			checkErrorAfterKernelLaunch("reduce6 in reducedI, dI00 (pow2) and 256 threads");
    			reduce6<T, 256, true><<< dimGrid, dimBlock, smemSize, stream06 >>>(dI01, d_odata_dI01, size);
    			checkErrorAfterKernelLaunch("reduce6 in reducedI, dI01 (pow2) and 256 threads");
    			reduce6<T, 256, true><<< dimGrid, dimBlock, smemSize, stream07 >>>(dI02, d_odata_dI02, size);
    			checkErrorAfterKernelLaunch("reduce6 in reducedI, dI02 (pow2) and 256 threads");
    			reduce6<T, 256, true><<< dimGrid, dimBlock, smemSize, stream08 >>>(dI10, d_odata_dI10, size);
    			checkErrorAfterKernelLaunch("reduce6 in reducedI, dI10 (pow2) and 256 threads");
    			reduce6<T, 256, true><<< dimGrid, dimBlock, smemSize, stream09 >>>(dI11, d_odata_dI11, size);
    			checkErrorAfterKernelLaunch("reduce6 in reducedI, dI11 (pow2) and 256 threads");
    			reduce6<T, 256, true><<< dimGrid, dimBlock, smemSize, stream10 >>>(dI12, d_odata_dI12, size);
    			checkErrorAfterKernelLaunch("reduce6 in reducedI, dI12 (pow2) and 256 threads");
    			reduce6<T, 256, true><<< dimGrid, dimBlock, smemSize, stream11 >>>(dI20, d_odata_dI20, size);
    			checkErrorAfterKernelLaunch("reduce6 in reducedI, dI20 (pow2) and 256 threads");
    			reduce6<T, 256, true><<< dimGrid, dimBlock, smemSize, stream12 >>>(dI21, d_odata_dI21, size);
    			checkErrorAfterKernelLaunch("reduce6 in reducedI, dI21 (pow2) and 256 threads");
    			reduce6<T, 256, true><<< dimGrid, dimBlock, smemSize, stream13 >>>(dI22, d_odata_dI22, size);
    			checkErrorAfterKernelLaunch("reduce6 in reducedI, dI22 (pow2) and 256 threads");
    			break;

    		case 128:
    			reduce6<T, 128, true><<< dimGrid, dimBlock, smemSize, stream01 >>>(dv, d_odata_dv, size);
    			reduce6<T, 128, true><<< dimGrid, dimBlock, smemSize, stream02 >>>(dcom0, d_odata_dcom0, size);
    			reduce6<T, 128, true><<< dimGrid, dimBlock, smemSize, stream03 >>>(dcom1, d_odata_dcom1, size);
    			reduce6<T, 128, true><<< dimGrid, dimBlock, smemSize, stream04 >>>(dcom2, d_odata_dcom2, size);
    			reduce6<T, 128, true><<< dimGrid, dimBlock, smemSize, stream05 >>>(dI00, d_odata_dI00, size);
    			reduce6<T, 128, true><<< dimGrid, dimBlock, smemSize, stream06 >>>(dI01, d_odata_dI01, size);
    			reduce6<T, 128, true><<< dimGrid, dimBlock, smemSize, stream07 >>>(dI02, d_odata_dI02, size);
    			reduce6<T, 128, true><<< dimGrid, dimBlock, smemSize, stream08 >>>(dI10, d_odata_dI10, size);
    			reduce6<T, 128, true><<< dimGrid, dimBlock, smemSize, stream09 >>>(dI11, d_odata_dI11, size);
    			reduce6<T, 128, true><<< dimGrid, dimBlock, smemSize, stream10 >>>(dI12, d_odata_dI12, size);
    			reduce6<T, 128, true><<< dimGrid, dimBlock, smemSize, stream11 >>>(dI20, d_odata_dI20, size);
    			reduce6<T, 128, true><<< dimGrid, dimBlock, smemSize, stream12 >>>(dI21, d_odata_dI21, size);
    			reduce6<T, 128, true><<< dimGrid, dimBlock, smemSize, stream13 >>>(dI22, d_odata_dI22, size);
    			break;

    		case 64:
    			reduce6<T, 64, true><<< dimGrid, dimBlock, smemSize, stream01 >>>(dv, d_odata_dv, size);
    			reduce6<T, 64, true><<< dimGrid, dimBlock, smemSize, stream02 >>>(dcom0, d_odata_dcom0, size);
    			reduce6<T, 64, true><<< dimGrid, dimBlock, smemSize, stream03 >>>(dcom1, d_odata_dcom1, size);
    			reduce6<T, 64, true><<< dimGrid, dimBlock, smemSize, stream04 >>>(dcom2, d_odata_dcom2, size);
    			reduce6<T, 64, true><<< dimGrid, dimBlock, smemSize, stream05 >>>(dI00, d_odata_dI00, size);
    			reduce6<T, 64, true><<< dimGrid, dimBlock, smemSize, stream06 >>>(dI01, d_odata_dI01, size);
    			reduce6<T, 64, true><<< dimGrid, dimBlock, smemSize, stream07 >>>(dI02, d_odata_dI02, size);
    			reduce6<T, 64, true><<< dimGrid, dimBlock, smemSize, stream08 >>>(dI10, d_odata_dI10, size);
    			reduce6<T, 64, true><<< dimGrid, dimBlock, smemSize, stream09 >>>(dI11, d_odata_dI11, size);
    			reduce6<T, 64, true><<< dimGrid, dimBlock, smemSize, stream10 >>>(dI12, d_odata_dI12, size);
    			reduce6<T, 64, true><<< dimGrid, dimBlock, smemSize, stream11 >>>(dI20, d_odata_dI20, size);
    			reduce6<T, 64, true><<< dimGrid, dimBlock, smemSize, stream12 >>>(dI21, d_odata_dI21, size);
    			reduce6<T, 64, true><<< dimGrid, dimBlock, smemSize, stream13 >>>(dI22, d_odata_dI22, size);
    			break;

    		case 32:
    			reduce6<T, 32, true><<< dimGrid, dimBlock, smemSize, stream01 >>>(dv, d_odata_dv, size);
    			reduce6<T, 32, true><<< dimGrid, dimBlock, smemSize, stream02 >>>(dcom0, d_odata_dcom0, size);
    			reduce6<T, 32, true><<< dimGrid, dimBlock, smemSize, stream03 >>>(dcom1, d_odata_dcom1, size);
    			reduce6<T, 32, true><<< dimGrid, dimBlock, smemSize, stream04 >>>(dcom2, d_odata_dcom2, size);
    			reduce6<T, 32, true><<< dimGrid, dimBlock, smemSize, stream05 >>>(dI00, d_odata_dI00, size);
    			reduce6<T, 32, true><<< dimGrid, dimBlock, smemSize, stream06 >>>(dI01, d_odata_dI01, size);
    			reduce6<T, 32, true><<< dimGrid, dimBlock, smemSize, stream07 >>>(dI02, d_odata_dI02, size);
    			reduce6<T, 32, true><<< dimGrid, dimBlock, smemSize, stream08 >>>(dI10, d_odata_dI10, size);
    			reduce6<T, 32, true><<< dimGrid, dimBlock, smemSize, stream09 >>>(dI11, d_odata_dI11, size);
    			reduce6<T, 32, true><<< dimGrid, dimBlock, smemSize, stream10 >>>(dI12, d_odata_dI12, size);
    			reduce6<T, 32, true><<< dimGrid, dimBlock, smemSize, stream11 >>>(dI20, d_odata_dI20, size);
    			reduce6<T, 32, true><<< dimGrid, dimBlock, smemSize, stream12 >>>(dI21, d_odata_dI21, size);
    			reduce6<T, 32, true><<< dimGrid, dimBlock, smemSize, stream13 >>>(dI22, d_odata_dI22, size);
    			break;

    		case 16:
    			reduce6<T, 16, true><<< dimGrid, dimBlock, smemSize, stream01 >>>(dv, d_odata_dv, size);
    			reduce6<T, 16, true><<< dimGrid, dimBlock, smemSize, stream02 >>>(dcom0, d_odata_dcom0, size);
    			reduce6<T, 16, true><<< dimGrid, dimBlock, smemSize, stream03 >>>(dcom1, d_odata_dcom1, size);
    			reduce6<T, 16, true><<< dimGrid, dimBlock, smemSize, stream04 >>>(dcom2, d_odata_dcom2, size);
    			reduce6<T, 16, true><<< dimGrid, dimBlock, smemSize, stream05 >>>(dI00, d_odata_dI00, size);
    			reduce6<T, 16, true><<< dimGrid, dimBlock, smemSize, stream06 >>>(dI01, d_odata_dI01, size);
    			reduce6<T, 16, true><<< dimGrid, dimBlock, smemSize, stream07 >>>(dI02, d_odata_dI02, size);
    			reduce6<T, 16, true><<< dimGrid, dimBlock, smemSize, stream08 >>>(dI10, d_odata_dI10, size);
    			reduce6<T, 16, true><<< dimGrid, dimBlock, smemSize, stream09 >>>(dI11, d_odata_dI11, size);
    			reduce6<T, 16, true><<< dimGrid, dimBlock, smemSize, stream10 >>>(dI12, d_odata_dI12, size);
    			reduce6<T, 16, true><<< dimGrid, dimBlock, smemSize, stream11 >>>(dI20, d_odata_dI20, size);
    			reduce6<T, 16, true><<< dimGrid, dimBlock, smemSize, stream12 >>>(dI21, d_odata_dI21, size);
    			reduce6<T, 16, true><<< dimGrid, dimBlock, smemSize, stream13 >>>(dI22, d_odata_dI22, size);
    			break;

    		case  8:
    			reduce6<T, 8, true><<< dimGrid, dimBlock, smemSize, stream01 >>>(dv, d_odata_dv, size);
    			reduce6<T, 8, true><<< dimGrid, dimBlock, smemSize, stream02 >>>(dcom0, d_odata_dcom0, size);
    			reduce6<T, 8, true><<< dimGrid, dimBlock, smemSize, stream03 >>>(dcom1, d_odata_dcom1, size);
    			reduce6<T, 8, true><<< dimGrid, dimBlock, smemSize, stream04 >>>(dcom2, d_odata_dcom2, size);
    			reduce6<T, 8, true><<< dimGrid, dimBlock, smemSize, stream05 >>>(dI00, d_odata_dI00, size);
    			reduce6<T, 8, true><<< dimGrid, dimBlock, smemSize, stream06 >>>(dI01, d_odata_dI01, size);
    			reduce6<T, 8, true><<< dimGrid, dimBlock, smemSize, stream07 >>>(dI02, d_odata_dI02, size);
    			reduce6<T, 8, true><<< dimGrid, dimBlock, smemSize, stream08 >>>(dI10, d_odata_dI10, size);
    			reduce6<T, 8, true><<< dimGrid, dimBlock, smemSize, stream09 >>>(dI11, d_odata_dI11, size);
    			reduce6<T, 8, true><<< dimGrid, dimBlock, smemSize, stream10 >>>(dI12, d_odata_dI12, size);
    			reduce6<T, 8, true><<< dimGrid, dimBlock, smemSize, stream11 >>>(dI20, d_odata_dI20, size);
    			reduce6<T, 8, true><<< dimGrid, dimBlock, smemSize, stream12 >>>(dI21, d_odata_dI21, size);
    			reduce6<T, 8, true><<< dimGrid, dimBlock, smemSize, stream13 >>>(dI22, d_odata_dI22, size);
    			break;

    		case  4:
    			reduce6<T, 4, true><<< dimGrid, dimBlock, smemSize, stream01 >>>(dv, d_odata_dv, size);
    			reduce6<T, 4, true><<< dimGrid, dimBlock, smemSize, stream02 >>>(dcom0, d_odata_dcom0, size);
    			reduce6<T, 4, true><<< dimGrid, dimBlock, smemSize, stream03 >>>(dcom1, d_odata_dcom1, size);
    			reduce6<T, 4, true><<< dimGrid, dimBlock, smemSize, stream04 >>>(dcom2, d_odata_dcom2, size);
    			reduce6<T, 4, true><<< dimGrid, dimBlock, smemSize, stream05 >>>(dI00, d_odata_dI00, size);
    			reduce6<T, 4, true><<< dimGrid, dimBlock, smemSize, stream06 >>>(dI01, d_odata_dI01, size);
    			reduce6<T, 4, true><<< dimGrid, dimBlock, smemSize, stream07 >>>(dI02, d_odata_dI02, size);
    			reduce6<T, 4, true><<< dimGrid, dimBlock, smemSize, stream08 >>>(dI10, d_odata_dI10, size);
    			reduce6<T, 4, true><<< dimGrid, dimBlock, smemSize, stream09 >>>(dI11, d_odata_dI11, size);
    			reduce6<T, 4, true><<< dimGrid, dimBlock, smemSize, stream10 >>>(dI12, d_odata_dI12, size);
    			reduce6<T, 4, true><<< dimGrid, dimBlock, smemSize, stream11 >>>(dI20, d_odata_dI20, size);
    			reduce6<T, 4, true><<< dimGrid, dimBlock, smemSize, stream12 >>>(dI21, d_odata_dI21, size);
    			reduce6<T, 4, true><<< dimGrid, dimBlock, smemSize, stream13 >>>(dI22, d_odata_dI22, size);
    			break;

    		case  2:
    			reduce6<T, 2, true><<< dimGrid, dimBlock, smemSize, stream01 >>>(dv, d_odata_dv, size);
    			reduce6<T, 2, true><<< dimGrid, dimBlock, smemSize, stream02 >>>(dcom0, d_odata_dcom0, size);
    			reduce6<T, 2, true><<< dimGrid, dimBlock, smemSize, stream03 >>>(dcom1, d_odata_dcom1, size);
    			reduce6<T, 2, true><<< dimGrid, dimBlock, smemSize, stream04 >>>(dcom2, d_odata_dcom2, size);
    			reduce6<T, 2, true><<< dimGrid, dimBlock, smemSize, stream05 >>>(dI00, d_odata_dI00, size);
    			reduce6<T, 2, true><<< dimGrid, dimBlock, smemSize, stream06 >>>(dI01, d_odata_dI01, size);
    			reduce6<T, 2, true><<< dimGrid, dimBlock, smemSize, stream07 >>>(dI02, d_odata_dI02, size);
    			reduce6<T, 2, true><<< dimGrid, dimBlock, smemSize, stream08 >>>(dI10, d_odata_dI10, size);
    			reduce6<T, 2, true><<< dimGrid, dimBlock, smemSize, stream09 >>>(dI11, d_odata_dI11, size);
    			reduce6<T, 2, true><<< dimGrid, dimBlock, smemSize, stream10 >>>(dI12, d_odata_dI12, size);
    			reduce6<T, 2, true><<< dimGrid, dimBlock, smemSize, stream11 >>>(dI20, d_odata_dI20, size);
    			reduce6<T, 2, true><<< dimGrid, dimBlock, smemSize, stream12 >>>(dI21, d_odata_dI21, size);
    			reduce6<T, 2, true><<< dimGrid, dimBlock, smemSize, stream13 >>>(dI22, d_odata_dI22, size);
    			break;

    		case  1:
    			reduce6<T, 1, true><<< dimGrid, dimBlock, smemSize, stream01 >>>(dv, d_odata_dv, size);
    			reduce6<T, 1, true><<< dimGrid, dimBlock, smemSize, stream02 >>>(dcom0, d_odata_dcom0, size);
    			reduce6<T, 1, true><<< dimGrid, dimBlock, smemSize, stream03 >>>(dcom1, d_odata_dcom1, size);
    			reduce6<T, 1, true><<< dimGrid, dimBlock, smemSize, stream04 >>>(dcom2, d_odata_dcom2, size);
    			reduce6<T, 1, true><<< dimGrid, dimBlock, smemSize, stream05 >>>(dI00, d_odata_dI00, size);
    			reduce6<T, 1, true><<< dimGrid, dimBlock, smemSize, stream06 >>>(dI01, d_odata_dI01, size);
    			reduce6<T, 1, true><<< dimGrid, dimBlock, smemSize, stream07 >>>(dI02, d_odata_dI02, size);
    			reduce6<T, 1, true><<< dimGrid, dimBlock, smemSize, stream08 >>>(dI10, d_odata_dI10, size);
    			reduce6<T, 1, true><<< dimGrid, dimBlock, smemSize, stream09 >>>(dI11, d_odata_dI11, size);
    			reduce6<T, 1, true><<< dimGrid, dimBlock, smemSize, stream10 >>>(dI12, d_odata_dI12, size);
    			reduce6<T, 1, true><<< dimGrid, dimBlock, smemSize, stream11 >>>(dI20, d_odata_dI20, size);
    			reduce6<T, 1, true><<< dimGrid, dimBlock, smemSize, stream12 >>>(dI21, d_odata_dI21, size);
    			reduce6<T, 1, true><<< dimGrid, dimBlock, smemSize, stream13 >>>(dI22, d_odata_dI22, size);
    			break;
    		}
    	}
    	else
    	{
    		switch (threads)
    		{
    		case 512:
    			reduce6<T, 512, false><<< dimGrid, dimBlock, smemSize, stream01 >>>(dv, d_odata_dv, size);
    			reduce6<T, 512, false><<< dimGrid, dimBlock, smemSize, stream02 >>>(dcom0, d_odata_dcom0, size);
    			reduce6<T, 512, false><<< dimGrid, dimBlock, smemSize, stream03 >>>(dcom1, d_odata_dcom1, size);
    			reduce6<T, 512, false><<< dimGrid, dimBlock, smemSize, stream04 >>>(dcom2, d_odata_dcom2, size);
    			reduce6<T, 512, false><<< dimGrid, dimBlock, smemSize, stream05 >>>(dI00, d_odata_dI00, size);
    			reduce6<T, 512, false><<< dimGrid, dimBlock, smemSize, stream06 >>>(dI01, d_odata_dI01, size);
    			reduce6<T, 512, false><<< dimGrid, dimBlock, smemSize, stream07 >>>(dI02, d_odata_dI02, size);
    			reduce6<T, 512, false><<< dimGrid, dimBlock, smemSize, stream08 >>>(dI10, d_odata_dI10, size);
    			reduce6<T, 512, false><<< dimGrid, dimBlock, smemSize, stream09 >>>(dI11, d_odata_dI11, size);
    			reduce6<T, 512, false><<< dimGrid, dimBlock, smemSize, stream10 >>>(dI12, d_odata_dI12, size);
    			reduce6<T, 512, false><<< dimGrid, dimBlock, smemSize, stream11 >>>(dI20, d_odata_dI20, size);
    			reduce6<T, 512, false><<< dimGrid, dimBlock, smemSize, stream12 >>>(dI21, d_odata_dI21, size);
    			reduce6<T, 512, false><<< dimGrid, dimBlock, smemSize, stream13 >>>(dI22, d_odata_dI22, size);
    			break;

    		case 256:
    			reduce6<T, 256, false><<< dimGrid, dimBlock, smemSize, stream01 >>>(dv, d_odata_dv, size);
    			checkErrorAfterKernelLaunch("reduce6 in reducedI, dv and 256 threads");
    			reduce6<T, 256, false><<< dimGrid, dimBlock, smemSize, stream02 >>>(dcom0, d_odata_dcom0, size);
    			checkErrorAfterKernelLaunch("reduce6 in reducedI, dcom0 and 256 threads");
    			reduce6<T, 256, false><<< dimGrid, dimBlock, smemSize, stream03 >>>(dcom1, d_odata_dcom1, size);
    			checkErrorAfterKernelLaunch("reduce6 in reducedI, dcom1 and 256 threads");
    			reduce6<T, 256, false><<< dimGrid, dimBlock, smemSize, stream04 >>>(dcom2, d_odata_dcom2, size);
    			checkErrorAfterKernelLaunch("reduce6 in reducedI, dcom2 and 256 threads");
    			reduce6<T, 256, false><<< dimGrid, dimBlock, smemSize, stream05 >>>(dI00, d_odata_dI00, size);
    			checkErrorAfterKernelLaunch("reduce6 in reducedI, dI00 and 256 threads");
    			reduce6<T, 256, false><<< dimGrid, dimBlock, smemSize, stream06 >>>(dI01, d_odata_dI01, size);
    			checkErrorAfterKernelLaunch("reduce6 in reducedI, dI01 and 256 threads");
    			reduce6<T, 256, false><<< dimGrid, dimBlock, smemSize, stream07 >>>(dI02, d_odata_dI02, size);
    			checkErrorAfterKernelLaunch("reduce6 in reducedI, dI02 and 256 threads");
    			reduce6<T, 256, false><<< dimGrid, dimBlock, smemSize, stream08 >>>(dI10, d_odata_dI10, size);
    			checkErrorAfterKernelLaunch("reduce6 in reducedI, dI10 and 256 threads");
    			reduce6<T, 256, false><<< dimGrid, dimBlock, smemSize, stream09 >>>(dI11, d_odata_dI11, size);
    			checkErrorAfterKernelLaunch("reduce6 in reducedI, dI11 and 256 threads");
    			reduce6<T, 256, false><<< dimGrid, dimBlock, smemSize, stream10 >>>(dI12, d_odata_dI12, size);
    			checkErrorAfterKernelLaunch("reduce6 in reducedI, dI12 and 256 threads");
    			reduce6<T, 256, false><<< dimGrid, dimBlock, smemSize, stream11 >>>(dI20, d_odata_dI20, size);
    			checkErrorAfterKernelLaunch("reduce6 in reducedI, dI20 and 256 threads");
    			reduce6<T, 256, false><<< dimGrid, dimBlock, smemSize, stream12 >>>(dI21, d_odata_dI21, size);
    			checkErrorAfterKernelLaunch("reduce6 in reducedI, dI21 and 256 threads");
    			reduce6<T, 256, false><<< dimGrid, dimBlock, smemSize, stream13 >>>(dI22, d_odata_dI22, size);
    			checkErrorAfterKernelLaunch("reduce6 in reducedI, dI22 and 256 threads");
    			break;

    		case 128:
    			reduce6<T, 128, false><<< dimGrid, dimBlock, smemSize, stream01 >>>(dv, d_odata_dv, size);
    			reduce6<T, 128, false><<< dimGrid, dimBlock, smemSize, stream02 >>>(dcom0, d_odata_dcom0, size);
    			reduce6<T, 128, false><<< dimGrid, dimBlock, smemSize, stream03 >>>(dcom1, d_odata_dcom1, size);
    			reduce6<T, 128, false><<< dimGrid, dimBlock, smemSize, stream04 >>>(dcom2, d_odata_dcom2, size);
    			reduce6<T, 128, false><<< dimGrid, dimBlock, smemSize, stream05 >>>(dI00, d_odata_dI00, size);
    			reduce6<T, 128, false><<< dimGrid, dimBlock, smemSize, stream06 >>>(dI01, d_odata_dI01, size);
    			reduce6<T, 128, false><<< dimGrid, dimBlock, smemSize, stream07 >>>(dI02, d_odata_dI02, size);
    			reduce6<T, 128, false><<< dimGrid, dimBlock, smemSize, stream08 >>>(dI10, d_odata_dI10, size);
    			reduce6<T, 128, false><<< dimGrid, dimBlock, smemSize, stream09 >>>(dI11, d_odata_dI11, size);
    			reduce6<T, 128, false><<< dimGrid, dimBlock, smemSize, stream10 >>>(dI12, d_odata_dI12, size);
    			reduce6<T, 128, false><<< dimGrid, dimBlock, smemSize, stream11 >>>(dI20, d_odata_dI20, size);
    			reduce6<T, 128, false><<< dimGrid, dimBlock, smemSize, stream12 >>>(dI21, d_odata_dI21, size);
    			reduce6<T, 128, false><<< dimGrid, dimBlock, smemSize, stream13 >>>(dI22, d_odata_dI22, size);
    			break;

    		case 64:
    			reduce6<T, 64, false><<< dimGrid, dimBlock, smemSize, stream01 >>>(dv, d_odata_dv, size);
    			reduce6<T, 64, false><<< dimGrid, dimBlock, smemSize, stream02 >>>(dcom0, d_odata_dcom0, size);
    			reduce6<T, 64, false><<< dimGrid, dimBlock, smemSize, stream03 >>>(dcom1, d_odata_dcom1, size);
    			reduce6<T, 64, false><<< dimGrid, dimBlock, smemSize, stream04 >>>(dcom2, d_odata_dcom2, size);
    			reduce6<T, 64, false><<< dimGrid, dimBlock, smemSize, stream05 >>>(dI00, d_odata_dI00, size);
    			reduce6<T, 64, false><<< dimGrid, dimBlock, smemSize, stream06 >>>(dI01, d_odata_dI01, size);
    			reduce6<T, 64, false><<< dimGrid, dimBlock, smemSize, stream07 >>>(dI02, d_odata_dI02, size);
    			reduce6<T, 64, false><<< dimGrid, dimBlock, smemSize, stream08 >>>(dI10, d_odata_dI10, size);
    			reduce6<T, 64, false><<< dimGrid, dimBlock, smemSize, stream09 >>>(dI11, d_odata_dI11, size);
    			reduce6<T, 64, false><<< dimGrid, dimBlock, smemSize, stream10 >>>(dI12, d_odata_dI12, size);
    			reduce6<T, 64, false><<< dimGrid, dimBlock, smemSize, stream11 >>>(dI20, d_odata_dI20, size);
    			reduce6<T, 64, false><<< dimGrid, dimBlock, smemSize, stream12 >>>(dI21, d_odata_dI21, size);
    			reduce6<T, 64, false><<< dimGrid, dimBlock, smemSize, stream13 >>>(dI22, d_odata_dI22, size);
    			break;

    		case 32:
    			reduce6<T, 32, false><<< dimGrid, dimBlock, smemSize, stream01 >>>(dv, d_odata_dv, size);
    			reduce6<T, 32, false><<< dimGrid, dimBlock, smemSize, stream02 >>>(dcom0, d_odata_dcom0, size);
    			reduce6<T, 32, false><<< dimGrid, dimBlock, smemSize, stream03 >>>(dcom1, d_odata_dcom1, size);
    			reduce6<T, 32, false><<< dimGrid, dimBlock, smemSize, stream04 >>>(dcom2, d_odata_dcom2, size);
    			reduce6<T, 32, false><<< dimGrid, dimBlock, smemSize, stream05 >>>(dI00, d_odata_dI00, size);
    			reduce6<T, 32, false><<< dimGrid, dimBlock, smemSize, stream06 >>>(dI01, d_odata_dI01, size);
    			reduce6<T, 32, false><<< dimGrid, dimBlock, smemSize, stream07 >>>(dI02, d_odata_dI02, size);
    			reduce6<T, 32, false><<< dimGrid, dimBlock, smemSize, stream08 >>>(dI10, d_odata_dI10, size);
    			reduce6<T, 32, false><<< dimGrid, dimBlock, smemSize, stream09 >>>(dI11, d_odata_dI11, size);
    			reduce6<T, 32, false><<< dimGrid, dimBlock, smemSize, stream10 >>>(dI12, d_odata_dI12, size);
    			reduce6<T, 32, false><<< dimGrid, dimBlock, smemSize, stream11 >>>(dI20, d_odata_dI20, size);
    			reduce6<T, 32, false><<< dimGrid, dimBlock, smemSize, stream12 >>>(dI21, d_odata_dI21, size);
    			reduce6<T, 32, false><<< dimGrid, dimBlock, smemSize, stream13 >>>(dI22, d_odata_dI22, size);
    			break;

    		case 16:
    			reduce6<T, 16, false><<< dimGrid, dimBlock, smemSize, stream01 >>>(dv, d_odata_dv, size);
    			reduce6<T, 16, false><<< dimGrid, dimBlock, smemSize, stream02 >>>(dcom0, d_odata_dcom0, size);
    			reduce6<T, 16, false><<< dimGrid, dimBlock, smemSize, stream03 >>>(dcom1, d_odata_dcom1, size);
    			reduce6<T, 16, false><<< dimGrid, dimBlock, smemSize, stream04 >>>(dcom2, d_odata_dcom2, size);
    			reduce6<T, 16, false><<< dimGrid, dimBlock, smemSize, stream05 >>>(dI00, d_odata_dI00, size);
    			reduce6<T, 16, false><<< dimGrid, dimBlock, smemSize, stream06 >>>(dI01, d_odata_dI01, size);
    			reduce6<T, 16, false><<< dimGrid, dimBlock, smemSize, stream07 >>>(dI02, d_odata_dI02, size);
    			reduce6<T, 16, false><<< dimGrid, dimBlock, smemSize, stream08 >>>(dI10, d_odata_dI10, size);
    			reduce6<T, 16, false><<< dimGrid, dimBlock, smemSize, stream09 >>>(dI11, d_odata_dI11, size);
    			reduce6<T, 16, false><<< dimGrid, dimBlock, smemSize, stream10 >>>(dI12, d_odata_dI12, size);
    			reduce6<T, 16, false><<< dimGrid, dimBlock, smemSize, stream11 >>>(dI20, d_odata_dI20, size);
    			reduce6<T, 16, false><<< dimGrid, dimBlock, smemSize, stream12 >>>(dI21, d_odata_dI21, size);
    			reduce6<T, 16, false><<< dimGrid, dimBlock, smemSize, stream13 >>>(dI22, d_odata_dI22, size);
    			break;

    		case  8:
    			reduce6<T, 8, false><<< dimGrid, dimBlock, smemSize, stream01 >>>(dv, d_odata_dv, size);
    			reduce6<T, 8, false><<< dimGrid, dimBlock, smemSize, stream02 >>>(dcom0, d_odata_dcom0, size);
    			reduce6<T, 8, false><<< dimGrid, dimBlock, smemSize, stream03 >>>(dcom1, d_odata_dcom1, size);
    			reduce6<T, 8, false><<< dimGrid, dimBlock, smemSize, stream04 >>>(dcom2, d_odata_dcom2, size);
    			reduce6<T, 8, false><<< dimGrid, dimBlock, smemSize, stream05 >>>(dI00, d_odata_dI00, size);
    			reduce6<T, 8, false><<< dimGrid, dimBlock, smemSize, stream06 >>>(dI01, d_odata_dI01, size);
    			reduce6<T, 8, false><<< dimGrid, dimBlock, smemSize, stream07 >>>(dI02, d_odata_dI02, size);
    			reduce6<T, 8, false><<< dimGrid, dimBlock, smemSize, stream08 >>>(dI10, d_odata_dI10, size);
    			reduce6<T, 8, false><<< dimGrid, dimBlock, smemSize, stream09 >>>(dI11, d_odata_dI11, size);
    			reduce6<T, 8, false><<< dimGrid, dimBlock, smemSize, stream10 >>>(dI12, d_odata_dI12, size);
    			reduce6<T, 8, false><<< dimGrid, dimBlock, smemSize, stream11 >>>(dI20, d_odata_dI20, size);
    			reduce6<T, 8, false><<< dimGrid, dimBlock, smemSize, stream12 >>>(dI21, d_odata_dI21, size);
    			reduce6<T, 8, false><<< dimGrid, dimBlock, smemSize, stream13 >>>(dI22, d_odata_dI22, size);
    			break;

    		case  4:
    			reduce6<T, 4, false><<< dimGrid, dimBlock, smemSize, stream01 >>>(dv, d_odata_dv, size);
    			reduce6<T, 4, false><<< dimGrid, dimBlock, smemSize, stream02 >>>(dcom0, d_odata_dcom0, size);
    			reduce6<T, 4, false><<< dimGrid, dimBlock, smemSize, stream03 >>>(dcom1, d_odata_dcom1, size);
    			reduce6<T, 4, false><<< dimGrid, dimBlock, smemSize, stream04 >>>(dcom2, d_odata_dcom2, size);
    			reduce6<T, 4, false><<< dimGrid, dimBlock, smemSize, stream05 >>>(dI00, d_odata_dI00, size);
    			reduce6<T, 4, false><<< dimGrid, dimBlock, smemSize, stream06 >>>(dI01, d_odata_dI01, size);
    			reduce6<T, 4, false><<< dimGrid, dimBlock, smemSize, stream07 >>>(dI02, d_odata_dI02, size);
    			reduce6<T, 4, false><<< dimGrid, dimBlock, smemSize, stream08 >>>(dI10, d_odata_dI10, size);
    			reduce6<T, 4, false><<< dimGrid, dimBlock, smemSize, stream09 >>>(dI11, d_odata_dI11, size);
    			reduce6<T, 4, false><<< dimGrid, dimBlock, smemSize, stream10 >>>(dI12, d_odata_dI12, size);
    			reduce6<T, 4, false><<< dimGrid, dimBlock, smemSize, stream11 >>>(dI20, d_odata_dI20, size);
    			reduce6<T, 4, false><<< dimGrid, dimBlock, smemSize, stream12 >>>(dI21, d_odata_dI21, size);
    			reduce6<T, 4, false><<< dimGrid, dimBlock, smemSize, stream13 >>>(dI22, d_odata_dI22, size);
    			break;

    		case  2:
    			reduce6<T, 2, false><<< dimGrid, dimBlock, smemSize, stream01 >>>(dv, d_odata_dv, size);
    			reduce6<T, 2, false><<< dimGrid, dimBlock, smemSize, stream02 >>>(dcom0, d_odata_dcom0, size);
    			reduce6<T, 2, false><<< dimGrid, dimBlock, smemSize, stream03 >>>(dcom1, d_odata_dcom1, size);
    			reduce6<T, 2, false><<< dimGrid, dimBlock, smemSize, stream04 >>>(dcom2, d_odata_dcom2, size);
    			reduce6<T, 2, false><<< dimGrid, dimBlock, smemSize, stream05 >>>(dI00, d_odata_dI00, size);
    			reduce6<T, 2, false><<< dimGrid, dimBlock, smemSize, stream06 >>>(dI01, d_odata_dI01, size);
    			reduce6<T, 2, false><<< dimGrid, dimBlock, smemSize, stream07 >>>(dI02, d_odata_dI02, size);
    			reduce6<T, 2, false><<< dimGrid, dimBlock, smemSize, stream08 >>>(dI10, d_odata_dI10, size);
    			reduce6<T, 2, false><<< dimGrid, dimBlock, smemSize, stream09 >>>(dI11, d_odata_dI11, size);
    			reduce6<T, 2, false><<< dimGrid, dimBlock, smemSize, stream10 >>>(dI12, d_odata_dI12, size);
    			reduce6<T, 2, false><<< dimGrid, dimBlock, smemSize, stream11 >>>(dI20, d_odata_dI20, size);
    			reduce6<T, 2, false><<< dimGrid, dimBlock, smemSize, stream12 >>>(dI21, d_odata_dI21, size);
    			reduce6<T, 2, false><<< dimGrid, dimBlock, smemSize, stream13 >>>(dI22, d_odata_dI22, size);
    			break;

    		case  1:
    			reduce6<T, 1, false><<< dimGrid, dimBlock, smemSize, stream01 >>>(dv, d_odata_dv, size);
    			reduce6<T, 1, false><<< dimGrid, dimBlock, smemSize, stream02 >>>(dcom0, d_odata_dcom0, size);
    			reduce6<T, 1, false><<< dimGrid, dimBlock, smemSize, stream03 >>>(dcom1, d_odata_dcom1, size);
    			reduce6<T, 1, false><<< dimGrid, dimBlock, smemSize, stream04 >>>(dcom2, d_odata_dcom2, size);
    			reduce6<T, 1, false><<< dimGrid, dimBlock, smemSize, stream05 >>>(dI00, d_odata_dI00, size);
    			reduce6<T, 1, false><<< dimGrid, dimBlock, smemSize, stream06 >>>(dI01, d_odata_dI01, size);
    			reduce6<T, 1, false><<< dimGrid, dimBlock, smemSize, stream07 >>>(dI02, d_odata_dI02, size);
    			reduce6<T, 1, false><<< dimGrid, dimBlock, smemSize, stream08 >>>(dI10, d_odata_dI10, size);
    			reduce6<T, 1, false><<< dimGrid, dimBlock, smemSize, stream09 >>>(dI11, d_odata_dI11, size);
    			reduce6<T, 1, false><<< dimGrid, dimBlock, smemSize, stream10 >>>(dI12, d_odata_dI12, size);
    			reduce6<T, 1, false><<< dimGrid, dimBlock, smemSize, stream11 >>>(dI20, d_odata_dI20, size);
    			reduce6<T, 1, false><<< dimGrid, dimBlock, smemSize, stream12 >>>(dI21, d_odata_dI21, size);
    			reduce6<T, 1, false><<< dimGrid, dimBlock, smemSize, stream13 >>>(dI22, d_odata_dI22, size);
    			break;
    		}
    	}
    	break;
    }

    /* Destroy the streams */
    cudaStreamDestroy(stream01);		cudaStreamDestroy(stream02);
    cudaStreamDestroy(stream03);		cudaStreamDestroy(stream04);
    cudaStreamDestroy(stream05);		cudaStreamDestroy(stream06);
    cudaStreamDestroy(stream07);		cudaStreamDestroy(stream08);
    cudaStreamDestroy(stream09);		cudaStreamDestroy(stream10);
    cudaStreamDestroy(stream11);		cudaStreamDestroy(stream12);
    cudaStreamDestroy(stream12);
}
/* Instantiate, but just for floats for now */
template void
reducedI<float>(int size, int numThreads, int numBlocks, int whichKernel,
		float *dv, float *dcom0, float *dcom1, float *dcom2, float *dI00,
		float *dI01, float *dI02, float *dI10, float *dI11, float *dI12,
		float *dI20, float *dI21, float *dI22, float *d_odata_dv, float
		*d_odata_dcom0, float *d_odata_dcom1, float *d_odata_dcom2, float
		*d_odata_dI00, float *d_odata_dI01, float *d_odata_dI02, float
		*d_odata_dI10, float *d_odata_dI11,	float *d_odata_dI12, float
		*d_odata_dI20, float *d_odata_dI21, float *d_odata_dI22);

////////////////////////////////////////////////////////////////////////////////
// Wrapper function for kernel launch
////////////////////////////////////////////////////////////////////////////////
template <class T>
void
maxz(int size, int threads, int blocks,
		int whichKernel, T *d_idata, T *d_odata)
{
	dim3 dimBlock(threads, 1, 1);
	dim3 dimGrid(blocks, 1, 1);

	// when there is only one warp per block, we need to allocate two warps
	// worth of shared memory so that we don't index shared memory out of bounds
	int smemSize = (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);

	// choose which of the optimized versions of reduction to launch
	switch (whichKernel)
	{
	case 6:
	default:
		if (isPow2(size))
		{
			switch (threads)
			{
			case 512:
				maxz6<T, 512, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case 256:
				maxz6<T, 256, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case 128:
				maxz6<T, 128, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case 64:
				maxz6<T,  64, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case 32:
				maxz6<T,  32, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case 16:
				maxz6<T,  16, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case  8:
				maxz6<T,   8, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case  4:
				maxz6<T,   4, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case  2:
				maxz6<T,   2, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case  1:
				maxz6<T,   1, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;
			}
		}
		else
		{
			switch (threads)
			{
			case 512:
				maxz6<T, 512, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case 256:
				maxz6<T, 256, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case 128:
				maxz6<T, 128, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case 64:
				maxz6<T,  64, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case 32:
				maxz6<T,  32, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case 16:
				maxz6<T,  16, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case  8:
				maxz6<T,   8, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case  4:
				maxz6<T,   4, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case  2:
				maxz6<T,   2, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case  1:
				maxz6<T,   1, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;
			}
		}
		break;
	}
}
// Instantiate the maxz function for 3 types
template void
maxz<int>(int size, int threads, int blocks,
		int whichKernel, int *d_idata, int *d_odata);

template void
maxz<float>(int size, int threads, int blocks,
              int whichKernel, float *d_idata, float *d_odata);

template void
maxz<double>(int size, int threads, int blocks,
               int whichKernel, double *d_idata, double *d_odata);


////////////////////////////////////////////////////////////////////////////////
// Wrapper function for kernel launch
////////////////////////////////////////////////////////////////////////////////
template <class T>
void
minz(int size, int threads, int blocks,
		int whichKernel, T *d_idata, T *d_odata)
{
	dim3 dimBlock(threads, 1, 1);
	dim3 dimGrid(blocks, 1, 1);

	// when there is only one warp per block, we need to allocate two warps
	// worth of shared memory so that we don't index shared memory out of bounds
	int smemSize = (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);

	// choose which of the optimized versions of reduction to launch
	switch (whichKernel)
	{
	case 6:
	default:
		if (isPow2(size))
		{
			switch (threads)
			{
			case 512:
				minz6<T, 512, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case 256:
				minz6<T, 256, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case 128:
				minz6<T, 128, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case 64:
				minz6<T,  64, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case 32:
				minz6<T,  32, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case 16:
				minz6<T,  16, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case  8:
				minz6<T,   8, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case  4:
				minz6<T,   4, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case  2:
				minz6<T,   2, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case  1:
				minz6<T,   1, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;
			}
		}
		else
		{
			switch (threads)
			{
			case 512:
				minz6<T, 512, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case 256:
				minz6<T, 256, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case 128:
				minz6<T, 128, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case 64:
				minz6<T,  64, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case 32:
				minz6<T,  32, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case 16:
				minz6<T,  16, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case  8:
				minz6<T,   8, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case  4:
				minz6<T,   4, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case  2:
				minz6<T,   2, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;

			case  1:
				minz6<T,   1, false><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);
				break;
			}
		}
		break;
	}
}
// Instantiate the minz function for 3 types
template void
minz<int>(int size, int threads, int blocks,
		int whichKernel, int *d_idata, int *d_odata);

template void
minz<float>(int size, int threads, int blocks,
              int whichKernel, float *d_idata, float *d_odata);

template void
minz<double>(int size, int threads, int blocks,
               int whichKernel, double *d_idata, double *d_odata);


////////////////////////////////////////////////////////////////////////////////
// Wrapper function for kernel launch
////////////////////////////////////////////////////////////////////////////////
template <class T>
void
maxzexp(int size, int threads, int blocks,
		int whichKernel, T *d_idata1, T *d_idata2, T *d_idata3, T *d_idata4,
		T *d_odata1, T *d_odata2, T *d_odata3, T *d_odata4,
		cudaStream_t *stream1, cudaStream_t *stream2, cudaStream_t *stream3,
		cudaStream_t *stream4)
{
	dim3 dimBlock(threads, 1, 1);
	dim3 dimGrid(blocks, 1, 1);

	// when there is only one warp per block, we need to allocate two warps
	// worth of shared memory so that we don't index shared memory out of bounds
	int smemSize = (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);

	// choose which of the optimized versions of reduction to launch
	switch (whichKernel)
	{
	case 6:
	default:
		if (isPow2(size))
		{
			switch (threads)
			{
			case 512:
				maxz6<T, 512, true><<< dimGrid, dimBlock, smemSize, *stream1 >>>(d_idata1, d_odata1, size);
				maxz6<T, 512, true><<< dimGrid, dimBlock, smemSize, *stream2 >>>(d_idata2, d_odata2, size);
				minz6<T, 512, true><<< dimGrid, dimBlock, smemSize, *stream3 >>>(d_idata3, d_odata3, size);
				minz6<T, 512, true><<< dimGrid, dimBlock, smemSize, *stream4 >>>(d_idata4, d_odata4, size);
				break;

			case 256:
				maxz6<T, 256, true><<< dimGrid, dimBlock, smemSize, *stream1 >>>(d_idata1, d_odata1, size);
				maxz6<T, 256, true><<< dimGrid, dimBlock, smemSize, *stream2 >>>(d_idata2, d_odata2, size);
				minz6<T, 256, true><<< dimGrid, dimBlock, smemSize, *stream3 >>>(d_idata3, d_odata3, size);
				minz6<T, 256, true><<< dimGrid, dimBlock, smemSize, *stream4 >>>(d_idata4, d_odata4, size);
				break;

			case 128:
				maxz6<T, 128, true><<< dimGrid, dimBlock, smemSize, *stream1 >>>(d_idata1, d_odata1, size);
				maxz6<T, 128, true><<< dimGrid, dimBlock, smemSize, *stream2 >>>(d_idata2, d_odata2, size);
				minz6<T, 128, true><<< dimGrid, dimBlock, smemSize, *stream3 >>>(d_idata3, d_odata3, size);
				minz6<T, 128, true><<< dimGrid, dimBlock, smemSize, *stream4 >>>(d_idata4, d_odata4, size);
				break;

			case 64:
				maxz6<T,  64, true><<< dimGrid, dimBlock, smemSize, *stream1 >>>(d_idata1, d_odata1, size);
				maxz6<T,  64, true><<< dimGrid, dimBlock, smemSize, *stream2 >>>(d_idata2, d_odata2, size);
				minz6<T,  64, true><<< dimGrid, dimBlock, smemSize, *stream3 >>>(d_idata3, d_odata3, size);
				minz6<T,  64, true><<< dimGrid, dimBlock, smemSize, *stream4 >>>(d_idata4, d_odata4, size);
				break;

			case 32:
				maxz6<T,  32, true><<< dimGrid, dimBlock, smemSize, *stream1 >>>(d_idata1, d_odata1, size);
				maxz6<T,  32, true><<< dimGrid, dimBlock, smemSize, *stream2 >>>(d_idata2, d_odata2, size);
				minz6<T,  32, true><<< dimGrid, dimBlock, smemSize, *stream3 >>>(d_idata3, d_odata3, size);
				minz6<T,  32, true><<< dimGrid, dimBlock, smemSize, *stream4 >>>(d_idata4, d_odata4, size);
				break;

			case 16:
				maxz6<T,  16, true><<< dimGrid, dimBlock, smemSize, *stream1 >>>(d_idata1, d_odata1, size);
				maxz6<T,  16, true><<< dimGrid, dimBlock, smemSize, *stream2 >>>(d_idata2, d_odata2, size);
				minz6<T,  16, true><<< dimGrid, dimBlock, smemSize, *stream3 >>>(d_idata3, d_odata3, size);
				minz6<T,  16, true><<< dimGrid, dimBlock, smemSize, *stream4 >>>(d_idata4, d_odata4, size);
				break;

			case  8:
				maxz6<T,   8, true><<< dimGrid, dimBlock, smemSize, *stream1 >>>(d_idata1, d_odata1, size);
				maxz6<T,   8, true><<< dimGrid, dimBlock, smemSize, *stream2 >>>(d_idata2, d_odata2, size);
				minz6<T,   8, true><<< dimGrid, dimBlock, smemSize, *stream3 >>>(d_idata3, d_odata3, size);
				minz6<T,   8, true><<< dimGrid, dimBlock, smemSize, *stream4 >>>(d_idata4, d_odata4, size);
				break;

			case  4:
				maxz6<T,   4, true><<< dimGrid, dimBlock, smemSize, *stream1 >>>(d_idata1, d_odata1, size);
				maxz6<T,   4, true><<< dimGrid, dimBlock, smemSize, *stream2 >>>(d_idata2, d_odata2, size);
				minz6<T,   4, true><<< dimGrid, dimBlock, smemSize, *stream3 >>>(d_idata3, d_odata3, size);
				minz6<T,   4, true><<< dimGrid, dimBlock, smemSize, *stream4 >>>(d_idata4, d_odata4, size);
				break;

			case  2:
				maxz6<T,   2, true><<< dimGrid, dimBlock, smemSize, *stream1 >>>(d_idata1, d_odata1, size);
				maxz6<T,   2, true><<< dimGrid, dimBlock, smemSize, *stream2 >>>(d_idata2, d_odata2, size);
				minz6<T,   2, true><<< dimGrid, dimBlock, smemSize, *stream3 >>>(d_idata3, d_odata3, size);
				minz6<T,   2, true><<< dimGrid, dimBlock, smemSize, *stream4 >>>(d_idata4, d_odata4, size);
				break;

			case  1:
				maxz6<T,   1, true><<< dimGrid, dimBlock, smemSize, *stream1 >>>(d_idata1, d_odata1, size);
				maxz6<T,   1, true><<< dimGrid, dimBlock, smemSize, *stream2 >>>(d_idata2, d_odata2, size);
				minz6<T,   1, true><<< dimGrid, dimBlock, smemSize, *stream3 >>>(d_idata3, d_odata3, size);
				minz6<T,   1, true><<< dimGrid, dimBlock, smemSize, *stream4 >>>(d_idata4, d_odata4, size);
				break;
			}
		}
		else
		{
			switch (threads)
			{
			case 512:
				maxz6<T, 512, false><<< dimGrid, dimBlock, smemSize, *stream1 >>>(d_idata1, d_odata1, size);
				maxz6<T, 512, false><<< dimGrid, dimBlock, smemSize, *stream2 >>>(d_idata2, d_odata2, size);
				minz6<T, 512, false><<< dimGrid, dimBlock, smemSize, *stream3 >>>(d_idata3, d_odata3, size);
				minz6<T, 512, false><<< dimGrid, dimBlock, smemSize, *stream4 >>>(d_idata4, d_odata4, size);
				break;

			case 256:
				maxz6<T, 256, false><<< dimGrid, dimBlock, smemSize, *stream1 >>>(d_idata1, d_odata1, size);
				maxz6<T, 256, false><<< dimGrid, dimBlock, smemSize, *stream2 >>>(d_idata2, d_odata2, size);
				minz6<T, 256, false><<< dimGrid, dimBlock, smemSize, *stream3 >>>(d_idata3, d_odata3, size);
				minz6<T, 256, false><<< dimGrid, dimBlock, smemSize, *stream4 >>>(d_idata4, d_odata4, size);
				break;

			case 128:
				maxz6<T, 128, false><<< dimGrid, dimBlock, smemSize, *stream1 >>>(d_idata1, d_odata1, size);
				maxz6<T, 128, false><<< dimGrid, dimBlock, smemSize, *stream2 >>>(d_idata2, d_odata2, size);
				minz6<T, 128, false><<< dimGrid, dimBlock, smemSize, *stream3 >>>(d_idata3, d_odata3, size);
				minz6<T, 128, false><<< dimGrid, dimBlock, smemSize, *stream4 >>>(d_idata4, d_odata4, size);
				break;

			case 64:
				maxz6<T,  64, false><<< dimGrid, dimBlock, smemSize, *stream1 >>>(d_idata1, d_odata1, size);
				maxz6<T,  64, false><<< dimGrid, dimBlock, smemSize, *stream2 >>>(d_idata2, d_odata2, size);
				minz6<T,  64, false><<< dimGrid, dimBlock, smemSize, *stream3 >>>(d_idata3, d_odata3, size);
				minz6<T,  64, false><<< dimGrid, dimBlock, smemSize, *stream4 >>>(d_idata4, d_odata4, size);
				break;

			case 32:
				maxz6<T,  32, false><<< dimGrid, dimBlock, smemSize, *stream1 >>>(d_idata1, d_odata1, size);
				maxz6<T,  32, false><<< dimGrid, dimBlock, smemSize, *stream2 >>>(d_idata2, d_odata2, size);
				minz6<T,  32, false><<< dimGrid, dimBlock, smemSize, *stream3 >>>(d_idata3, d_odata3, size);
				minz6<T,  32, false><<< dimGrid, dimBlock, smemSize, *stream4 >>>(d_idata4, d_odata4, size);
				break;

			case 16:
				maxz6<T,  16, false><<< dimGrid, dimBlock, smemSize, *stream1 >>>(d_idata1, d_odata1, size);
				maxz6<T,  16, false><<< dimGrid, dimBlock, smemSize, *stream2 >>>(d_idata2, d_odata2, size);
				minz6<T,  16, false><<< dimGrid, dimBlock, smemSize, *stream3 >>>(d_idata3, d_odata3, size);
				minz6<T,  16, false><<< dimGrid, dimBlock, smemSize, *stream4 >>>(d_idata4, d_odata4, size);
				break;

			case  8:
				maxz6<T,   8, false><<< dimGrid, dimBlock, smemSize, *stream1 >>>(d_idata1, d_odata1, size);
				maxz6<T,   8, false><<< dimGrid, dimBlock, smemSize, *stream2 >>>(d_idata2, d_odata2, size);
				minz6<T,   8, false><<< dimGrid, dimBlock, smemSize, *stream3 >>>(d_idata3, d_odata3, size);
				minz6<T,   8, false><<< dimGrid, dimBlock, smemSize, *stream4 >>>(d_idata4, d_odata4, size);
				break;

			case  4:
				maxz6<T,   4, false><<< dimGrid, dimBlock, smemSize, *stream1 >>>(d_idata1, d_odata1, size);
				maxz6<T,   4, false><<< dimGrid, dimBlock, smemSize, *stream2 >>>(d_idata2, d_odata2, size);
				minz6<T,   4, false><<< dimGrid, dimBlock, smemSize, *stream3 >>>(d_idata3, d_odata3, size);
				minz6<T,   4, false><<< dimGrid, dimBlock, smemSize, *stream4 >>>(d_idata4, d_odata4, size);
				break;

			case  2:
				maxz6<T,   2, false><<< dimGrid, dimBlock, smemSize, *stream1 >>>(d_idata1, d_odata1, size);
				maxz6<T,   2, false><<< dimGrid, dimBlock, smemSize, *stream2 >>>(d_idata2, d_odata2, size);
				minz6<T,   2, false><<< dimGrid, dimBlock, smemSize, *stream3 >>>(d_idata3, d_odata3, size);
				minz6<T,   2, false><<< dimGrid, dimBlock, smemSize, *stream4 >>>(d_idata4, d_odata4, size);
				break;

			case  1:
				maxz6<T,   1, false><<< dimGrid, dimBlock, smemSize, *stream1 >>>(d_idata1, d_odata1, size);
				maxz6<T,   1, false><<< dimGrid, dimBlock, smemSize, *stream2 >>>(d_idata2, d_odata2, size);
				minz6<T,   1, false><<< dimGrid, dimBlock, smemSize, *stream3 >>>(d_idata3, d_odata3, size);
				minz6<T,   1, false><<< dimGrid, dimBlock, smemSize, *stream4 >>>(d_idata4, d_odata4, size);
				break;
			}
		}
		break;
	}
}
// Instantiate the maxz function for 3 types
template void
maxzexp<int>(int size, int threads, int blocks,	int whichKernel, int *d_idata1,
		int *d_idata2, int *d_idata3, int *d_idata4, int *d_odata1, int *d_odata2,
		int *d_odata3, int *d_odata4, cudaStream_t *stream1, cudaStream_t *stream2,
		cudaStream_t *stream3, cudaStream_t *stream4);

template void
maxzexp<float>(int size, int threads, int blocks, int whichKernel, float *d_idata1,
		float *d_idata2, float *d_idata3, float *d_idata4, float *d_odata1,
		float *d_odata2, float *d_odata3, float *d_odata4, cudaStream_t *stream1,
		cudaStream_t *stream2, cudaStream_t *stream3, cudaStream_t *stream4);

template void
maxzexp<double>(int size, int threads, int blocks, int whichKernel, double *d_idata1,
		double *d_idata2, double *d_idata3, double *d_idata4, double *d_odata1,
		double *d_odata2, double *d_odata3, double *d_odata4, cudaStream_t *stream1,
		cudaStream_t *stream2, cudaStream_t *stream3, cudaStream_t *stream4);


__global__ void set_idata_zmax_krnl(struct dat_t *ddat, float *d_idata,
		int set, int frm, int size) {
	/* MULTI-threaded kernel */
	int i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i < size){
		d_idata[i] = ddat->set[set].desc.deldop.frame[frm].pos.z_s[i];
	}
}
__global__ void set_idata_zmax_all_frames_krnl(struct dat_t *ddat,
		float *d_idata,	int set, int frm, int frame_size) {
	/* MULTI-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;

	if (offset < frame_size){
		d_idata[offset] = ddat->set[set].desc.deldop.frame[frm].pos.z_s[offset];
	}
}
__global__ void set_idata_pntr_krnl(struct dat_t *ddat, float *d_idata,
		int set, int frm, int size) {
	/* MULTI-threaded kernel */
	int i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i < size){

		switch (ddat->set[set].type) {
		case DELAY:
			d_idata[i] = ddat->set[set].desc.deldop.frame[frm].fit_s[i];
			break;
		case DOPPLER:
			d_idata[i] = ddat->set[set].desc.doppler.frame[frm].fit_s[i];
			break;
		}
	}
}
__global__ void set_idata_modarea_krnl(struct mod_t *dmod, float *d_idata,
		int c, int size) {
	/* MULTI-threaded kernel */
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < size)
		d_idata[i] = dmod->shape.comp[c].real.f[i].area;
}
__global__ void set_dv_dcom_di_krnl(float *d_odata_dv, float *d_odata_dcom0,
		float *d_odata_dcom1, float *d_odata_dcom2, float *d_odata_dI00,
		float *d_odata_dI01, float *d_odata_dI02, float *d_odata_dI10,
		float *d_odata_dI11, float *d_odata_dI12, float *d_odata_dI20,
		float *d_odata_dI21, float *d_odata_dI22, struct mod_t *dmod, int c) {
	/* Single threaded kernel to update the model with volume, COM, and inertia */
	/* Note that this kernel ignores multi-threaded models for now */
	if (threadIdx.x == 0) {
		dmod->shape.comp[c].volume = dmod->shape.volume = d_odata_dv[0];
		dmod->shape.comp[c].com[0] = dmod->shape.com[0] = d_odata_dcom0[0];
		dmod->shape.comp[c].com[1] = dmod->shape.com[1] = d_odata_dcom1[0];
		dmod->shape.comp[c].com[2] = dmod->shape.com[2] = d_odata_dcom2[0];
		dmod->shape.comp[c].inertia[0][0] = dmod->shape.inertia[0][0] = d_odata_dI00[0];
		dmod->shape.comp[c].inertia[0][1] = dmod->shape.inertia[0][1] = d_odata_dI01[0];
		dmod->shape.comp[c].inertia[0][2] = dmod->shape.inertia[0][2] = d_odata_dI02[0];
		dmod->shape.comp[c].inertia[1][0] = dmod->shape.inertia[1][0] = d_odata_dI10[0];
		dmod->shape.comp[c].inertia[1][1] = dmod->shape.inertia[1][1] = d_odata_dI11[0];
		dmod->shape.comp[c].inertia[1][2] = dmod->shape.inertia[1][2] = d_odata_dI12[0];
		dmod->shape.comp[c].inertia[2][0] = dmod->shape.inertia[2][0] = d_odata_dI20[0];
		dmod->shape.comp[c].inertia[2][1] = dmod->shape.inertia[2][1] = d_odata_dI21[0];
		dmod->shape.comp[c].inertia[2][2] = dmod->shape.inertia[2][2] = d_odata_dI22[0];
	}
}
__global__ void set_xlim_ylim_krnl(struct dat_t *ddat, int set, int frm, int src,
		float *d_odata_imax, float *d_odata_imin, float *d_odata_jmax,
		float *d_odata_jmin, float *minmax_overall) {
	/* Single-threaded kernel to update pos->xlim and ylim and also the overall
	 * model limits regardless of the POS limits */
	int n;
	if (threadIdx.x == 0) {
		/* First set the overall model limits, regardless of POS frame limits */
		minmax_overall[0] = d_odata_imin[0];
		minmax_overall[1] = d_odata_imax[0];
		minmax_overall[2] = d_odata_jmin[0];
		minmax_overall[3] = d_odata_jmax[0];

		/* Now set pos->xlim and pos->ylim depending on data type */
		switch(ddat->set[set].type) {
		case DELAY:
			n = ddat->set[set].desc.deldop.frame[frm].pos.n;
			if (src) {
				ddat->set[set].desc.deldop.frame[frm].pos.xlim2[0] = max((int)d_odata_imin[0], -n);
				ddat->set[set].desc.deldop.frame[frm].pos.xlim2[1] = min((int)d_odata_imax[0],  n);
				ddat->set[set].desc.deldop.frame[frm].pos.ylim2[0] = max((int)d_odata_jmin[0], -n);
				ddat->set[set].desc.deldop.frame[frm].pos.ylim2[1] = min((int)d_odata_jmax[0],  n);
			}
			else {
				ddat->set[set].desc.deldop.frame[frm].pos.xlim[0] = max((int)d_odata_imin[0], -n);
				ddat->set[set].desc.deldop.frame[frm].pos.xlim[1] = min((int)d_odata_imax[0],  n);
				ddat->set[set].desc.deldop.frame[frm].pos.ylim[0] = max((int)d_odata_jmin[0], -n);
				ddat->set[set].desc.deldop.frame[frm].pos.ylim[1] = min((int)d_odata_jmax[0],  n);
			}
			break;
		case DOPPLER:
			n = ddat->set[set].desc.doppler.frame[frm].pos.n;
			if (src) {
				ddat->set[set].desc.doppler.frame[frm].pos.xlim2[0] = max((int)d_odata_imin[0], -n);
				ddat->set[set].desc.doppler.frame[frm].pos.xlim2[1] = min((int)d_odata_imax[0],  n);
				ddat->set[set].desc.doppler.frame[frm].pos.ylim2[0] = max((int)d_odata_jmin[0], -n);
				ddat->set[set].desc.doppler.frame[frm].pos.ylim2[1] = min((int)d_odata_jmax[0],  n);
			}
			else {
				ddat->set[set].desc.doppler.frame[frm].pos.xlim[0] = max((int)d_odata_imin[0], -n);
				ddat->set[set].desc.doppler.frame[frm].pos.xlim[1] = min((int)d_odata_imax[0],  n);
				ddat->set[set].desc.doppler.frame[frm].pos.ylim[0] = max((int)d_odata_jmin[0], -n);
				ddat->set[set].desc.doppler.frame[frm].pos.ylim[1] = min((int)d_odata_jmax[0],  n);
			}
			break;
		case POS:
			//			if (src) {
			//				ddat->set[set].desc.poset.frame[frm].pos.xlim2[0] = max((int)d_odata_imin[0], -n);
			//				ddat->set[set].desc.poset.frame[frm].pos.xlim2[1] = min((int)d_odata_imax[0],  n);
			//				ddat->set[set].desc.poset.frame[frm].pos.ylim2[0] = max((int)d_odata_jmin[0], -n);
			//				ddat->set[set].desc.poset.frame[frm].pos.ylim2[1] = min((int)d_odata_jmax[0],  n);
			//			}
			//			else {
			//				ddat->set[set].desc.poset.frame[frm].pos.xlim[0] = max((int)d_odata_imin[0], -n);
			//				ddat->set[set].desc.poset.frame[frm].pos.xlim[1] = min((int)d_odata_imax[0],  n);
			//				ddat->set[set].desc.poset.frame[frm].pos.ylim[0] = max((int)d_odata_jmin[0], -n);
			//				ddat->set[set].desc.poset.frame[frm].pos.ylim[1] = min((int)d_odata_jmax[0],  n);
			//			}
			//			break;
		case LGHTCRV:
			n = ddat->set[set].desc.lghtcrv.rend[frm].pos.n;
			if (src) {
				ddat->set[set].desc.lghtcrv.rend[frm].pos.xlim2[0] = max((int)d_odata_imin[0], -n);
				ddat->set[set].desc.lghtcrv.rend[frm].pos.xlim2[1] = min((int)d_odata_imax[0],  n);
				ddat->set[set].desc.lghtcrv.rend[frm].pos.ylim2[0] = max((int)d_odata_jmin[0], -n);
				ddat->set[set].desc.lghtcrv.rend[frm].pos.ylim2[1] = min((int)d_odata_jmax[0],  n);
			}
			else {
				ddat->set[set].desc.lghtcrv.rend[frm].pos.xlim[0] = max((int)d_odata_imin[0], -n);
				ddat->set[set].desc.lghtcrv.rend[frm].pos.xlim[1] = min((int)d_odata_imax[0],  n);
				ddat->set[set].desc.lghtcrv.rend[frm].pos.ylim[0] = max((int)d_odata_jmin[0], -n);
				ddat->set[set].desc.lghtcrv.rend[frm].pos.ylim[1] = min((int)d_odata_jmax[0],  n);
			}
			break;
		}
	}
}
__global__ void zmax_all_frames_finalize_krnl(struct dat_t *ddat, float *zmax,
	int nframes, int s) {
		/* single-threaded kernel */
		float sum = 0.0;
		if (threadIdx.x == 0) {
			for (int f=0; f<nframes; f++)
				sum += (zmax[f] * ddat->set[s].desc.deldop.frame[f].weight);

			zmax[0] = sum;
		}

	}
__global__ void deldop_xsec_frame_finalize_krnl(struct dat_t *ddat,
	int s, int f, float value, float *xsec) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		xsec[f] = ddat->set[s].desc.deldop.frame[f].overflow_xsec;
		xsec[f] += value;
		xsec[f] *= ddat->set[s].desc.deldop.frame[f].cal.val;
	}
}
__global__ void deldop_xsec_set_finalize_krnl(struct dat_t *ddat,
	int s, float *xsec, int nframes) {
		/* Single-threaded kernel */
		if (threadIdx.x == 0) {
			for (int f=0; f<nframes; f++)
				reduction_sum_rad_xsec += xsec[f]*ddat->set[s].desc.deldop.frame[f].weight;
		}
	}
__global__ void c2af_set_data_krnl(float **input, float *d_idata, int frm, int frmsz) {
	/* frmsz-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	if (offset < frmsz) {
		d_idata[offset] = input[frm][offset];
	}
}
__host__ float compute_deldop_xsec_pr6(struct dat_t *ddat, int ndel, int ndop,
		int set, int frm) {
	/* Function calculates a delay-Doppler frame's radar cross section with
	 * Nvidia's reduction sample code (simplified and adapted for use with
	 * shape).  The function returns the cross section as a float 	 */

	int s, size = ndel*ndop;		// array size
	int maxThreads = maxThreadsPerBlock;		// max # of threads per block
	int maxBlocks = 2048;		// max # of blocks per grid
	int whichKernel = 6;		// id of reduction kernel
	int numBlocks = 0;			// initialize numBlocks
	int numThreads = 0;			// initialize numThreads
	float xsec = 0.0;			// radar cross section; return value
	float *d_odata;				// temp. float array for reduction output
	float *d_idata; 			// temp. float arrays for reduction input
	float *h_odata;				// the host output array
	float2 xblock_ythread;		// used for return value of getNumBlocksAndThreads

	dim3 BLK,THD;
	BLK.x = floor((maxThreadsPerBlock - 1 + size)/maxThreadsPerBlock);
	THD.x = maxThreadsPerBlock; // Thread block dimensions

	/* Allocate memory for d_idata, then set that pointer equal to the right
	 * data set and frame to the right deldop fit array	 */
	cudaCalloc((void**)&d_idata, sizeof(float), size);

	set_idata_pntr_krnl<<<BLK,THD>>>(ddat, d_idata, set, frm, size);
	checkErrorAfterKernelLaunch("set_idata_pntr_krnl in compute_deldop_xsec");

	/* Find number of blocks & threads needed for reduction call */
	xblock_ythread = getNumBlocksAndThreads(size, maxBlocks, maxThreads);
	numBlocks = xblock_ythread.x;
	numThreads = xblock_ythread.y;

	/* Allocate memory for d_odata and d_odata2 with enough elements to hold
	 * the reduction of each block during the first call */
	cudaCalloc((void**)&d_odata,  sizeof(float), numBlocks);
	h_odata = (float *) malloc(numBlocks*sizeof(float));

	/* Call reduction for first time */
	reduce<float>(size, numThreads, numBlocks, whichKernel, d_idata, d_odata);
	checkErrorAfterKernelLaunch("reduce<float> in compute_deldop_xsec");

	/* Reset d_idata for later use as buffer */
	cudaMemset(d_idata, 0, size*sizeof(float));

	/* Now sum partial block sums on GPU, using the reduce6<> kernel */
	s = numBlocks;

	while (s > 1)
	{
		int threads = 0, blocks = 0;
		xblock_ythread = getNumBlocksAndThreads(s, maxBlocks, maxThreads);
		blocks = xblock_ythread.x;
		threads = xblock_ythread.y;

		/* Copy the first d_odata back into d_idata2  */
		cudaMemcpy(d_idata, d_odata, s*sizeof(float), cudaMemcpyDeviceToDevice);

		reduce<float>(s, threads, blocks, whichKernel, d_idata, d_odata);

		if (whichKernel < 3)
			s = (s + threads - 1) / threads;
		else
			s = (s + (threads*2-1)) / (threads*2);
		if (s > 1)
			printf("s is bigger than one");
	}

	gpuErrchk(cudaMemcpy(h_odata, d_odata, 2*sizeof(float), cudaMemcpyDeviceToHost));
	xsec = h_odata[0];
	free(h_odata);
	cudaFree(d_odata);
	cudaFree(d_idata);
	return xsec;
}
__host__ float compute_deldop_xsec_snglkrnl(struct dat_t *ddat, int ndel, int ndop,
		int set, int frm) {
	/* Function calculates a delay-Doppler frame's radar cross section with
	 * Nvidia's reduction sample code (simplified and adapted for use with
	 * shape).  The function returns the cross section as a float 	 */

	int size = ndel*ndop;		// array size
	int maxThreads = maxThreadsPerBlock;		// max # of threads per block
	int maxBlocks = 2048;		// max # of blocks per grid
	int numBlocks = 0;			// initialize numBlocks
	int numThreads = 0;			// initialize numThreads
	float xsec = 0.0;			// radar cross section; return value
	float *d_odata;				// temp. float array for reduction output
	float *d_idata; 			// temp. float arrays for reduction input
	float2 xblock_ythread;		// used for return value of getNumBlocksAndThreads

	dim3 BLK,THD;
	BLK.x = floor((maxThreadsPerBlock - 1 + size)/maxThreadsPerBlock);
	THD.x = maxThreadsPerBlock; // Thread block dimensions

	/* Find number of blocks & threads needed for reduction call */
	xblock_ythread = getNumBlocksAndThreads(size, maxBlocks, maxThreads);
	numBlocks = xblock_ythread.x;
	numThreads = xblock_ythread.y;
	dim3 dimBlock(numThreads, 1, 1);
	dim3 dimGrid(numBlocks, 1, 1);

	/* Allocate memory for d_idata and d_odata */
	cudaCalloc((void**)&d_idata, sizeof(float), size);
	cudaCalloc((void**)&d_odata,  sizeof(float), numBlocks);

	set_idata_pntr_krnl<<<BLK,THD>>>(ddat, d_idata, set, frm, size);
	checkErrorAfterKernelLaunch("set_idata_pntr_krnl in compute_deldop_xsec");

	/* Call reduction  */
	device_reduce_block_atomic_kernel<<< dimGrid, dimBlock>>>(d_idata, d_odata, size);
	checkErrorAfterKernelLaunch("device_reduce_block_atomic_kernel (compute_deldop_xsec_snglkrnl)");
	deviceSyncAfterKernelLaunch("device_reduce_block_atomic_kernel (compute_deldop_xsec_snglkrnl)");

	xsec = d_odata[0];
	cudaFree(d_odata);
	cudaFree(d_idata);
	return xsec;
}
__host__ float compute_deldop_xsec_all_frames(struct dat_t *ddat, int ndel, int ndop,
		int set, int nframes) {
	/* Function calculates a delay-Doppler frame's radar cross section with
	 * Nvidia's reduction sample code (simplified and adapted for use with
	 * shape).  The function returns the cross section as a float 	 */

	int size = ndel*ndop;		// array size
	int maxThreads = maxThreadsPerBlock;		// max # of threads per block
	int maxBlocks = 2048;		// max # of blocks per grid
	int numBlocks = 0;			// initialize numBlocks
	int numThreads = 0;			// initialize numThreads
	float *xsec, xsec_set;   	// radar cross section; return value
	float *d_odata;				// temp. float array for reduction output
	float *d_idata; 			// temp. float arrays for reduction input
	float2 xblock_ythread;		// used for return value of getNumBlocksAndThreads

	dim3 BLK,THD;
	BLK.x = floor((maxThreadsPerBlock - 1 + size)/maxThreadsPerBlock);
	THD.x = maxThreadsPerBlock; // Thread block dimensions

	/* Find number of blocks & threads needed for reduction call */
	xblock_ythread = getNumBlocksAndThreads(size, maxBlocks, maxThreads);
	numBlocks = xblock_ythread.x;
	numThreads = xblock_ythread.y;
	dim3 dimBlock(numThreads, 1, 1);
	dim3 dimGrid(numBlocks, 1, 1);

	/* Allocate memory for d_idata and d_odata */
	cudaCalloc((void**)&d_idata, sizeof(float), size);
	cudaMallocManaged((void**)&d_odata, sizeof(float)*numBlocks, cudaMemAttachHost);
	//cudaMallocManaged((void**)&xsec, sizeof(float)*nframes, cudaMemAttachHost);
	cudaCalloc((void**)&xsec, sizeof(float), nframes);

	for (int i=0; i<nframes; i++)
		xsec[i]=0.0;

	/* Start loop through frames */
	for (int f=0; f<nframes; f++) {
		//set_idata_zmax_all_frames_krnl<<<BLK,THD>>>(ddat, d_idata,	set, f, size);
		set_idata_pntr_krnl<<<BLK,THD>>>(ddat, d_idata, set, f, size);
		checkErrorAfterKernelLaunch("set_idata_pntr_krnl in compute_deldop_xsec_all_frames");

		/* Call reduction  */
		device_reduce_block_atomic_kernel<<< dimGrid, dimBlock>>>(d_idata, d_odata, size);
		checkErrorAfterKernelLaunch("device_reduce_block_atomic_kernel (compute_deldop_xsec_all_frames)");
		//deviceSyncAfterKernelLaunch("device_reduce_block_atomic_kernel (compute_deldop_xsec_all_frames)");

		/* Calculate frames weighting */
		deldop_xsec_frame_finalize_krnl<<<1,1>>>(ddat, set, f, d_odata[0], xsec);
		checkErrorAfterKernelLaunch("deldop_xsec_frame_finalize_krnl");
		deviceSyncAfterKernelLaunch("");
	}
	/* Now finalize the set and return value */
	deldop_xsec_set_finalize_krnl<<<1,1>>>(ddat, set, xsec, nframes);
	checkErrorAfterKernelLaunch("deldop_xsec_set_finalize_krnl");
	deviceSyncAfterKernelLaunch("deldop_xsec_set_finalize_krnl");
	gpuErrchk(cudaMemcpyFromSymbol(&xsec_set, reduction_sum_rad_xsec,
			sizeof(float), 0, cudaMemcpyDeviceToHost));

	cudaFree(d_odata);
	cudaFree(d_idata);
	cudaFree(xsec);
	return xsec_set;
}
__host__ float compute_pos_zmax(struct dat_t *ddat, int size,
		int set, int frm) {
	/* Function calculates the zmax in a pos->z arrary with
	 * Nvidia's reduction sample code (simplified and adapted for use with
	 * shape).  The function returns the cross section as a float 	 */

	int s;						// array size
	int maxThreads = maxThreadsPerBlock;		// max # of threads per block
	int maxBlocks = 2048;		// max # of blocks per grid
	int whichKernel = 6;		// id of reduction kernel
	int numBlocks = 0;			// initialize numBlocks
	int numThreads = 0;			// initialize numThreads
	float zmax = 0.0;			// radar cross section; return value
	float *d_odata;				// temp. float array for reduction output
	float *d_idata; 			// temp. float arrays for reduction input
	float *h_odata;				// the host output array
	float2 xblock_ythread;		// used for return value of getNumBlocksAndThreads

	dim3 BLK,THD;
	BLK.x = floor((maxThreadsPerBlock - 1 + size)/maxThreadsPerBlock);
	THD.x = maxThreadsPerBlock; // Thread block dimensions

	/* Allocate memory for d_idata, then set that pointer equal to the right
	 * data set and frame to the right deldop fit array	 */
	cudaCalloc((void**)&d_idata, sizeof(float), size);

	set_idata_zmax_krnl<<<BLK,THD>>>(ddat, d_idata, set, frm, size);
	checkErrorAfterKernelLaunch("set_idata_zmax_krnl in compute_zmax");

	/* Find number of blocks & threads needed for reduction call */
	xblock_ythread = getNumBlocksAndThreads(size, maxBlocks, maxThreads);
	numBlocks = xblock_ythread.x;
	numThreads = xblock_ythread.y;

	/* Allocate memory for d_odata and d_odata2 with enough elements to hold
	 * the reduction of each block during the first call */
	cudaCalloc((void**)&d_odata,  sizeof(float), numBlocks);
	h_odata = (float *) malloc(numBlocks*sizeof(float));

	/* Call maxz for first time */
	maxz<float>(size, numThreads, numBlocks, whichKernel, d_idata, d_odata);
	checkErrorAfterKernelLaunch("maxz<float> in compute_zmax");

	/* Reset d_idata for later use as buffer */
	cudaMemset(d_idata, 0, size*sizeof(float));

	/* Now sum partial block sums on GPU, using the maxz6<> kernel */
	s = numBlocks;

	while (s > 1)
	{
		int threads = 0, blocks = 0;
		xblock_ythread = getNumBlocksAndThreads(s, maxBlocks, maxThreads);
		blocks = xblock_ythread.x;
		threads = xblock_ythread.y;

		/* Copy the first d_odata back into d_idata2  */
		cudaMemcpy(d_idata, d_odata, s*sizeof(float), cudaMemcpyDeviceToDevice);

		maxz<float>(s, threads, blocks, whichKernel, d_idata, d_odata);

		if (whichKernel < 3)
			s = (s + threads - 1) / threads;
		else
			s = (s + (threads*2-1)) / (threads*2);
		if (s > 1)
			printf("s is bigger than one");
	}

	gpuErrchk(cudaMemcpy(h_odata, d_odata, 2*sizeof(float), cudaMemcpyDeviceToHost));
	zmax = h_odata[0];
	free(h_odata);
	cudaFree(d_odata);
	cudaFree(d_idata);
	return zmax;
}
//__host__ float compute_pos_zmax_streams(struct dat_t *ddat, int size,
//		int set, int nframes) {
//	/* Function calculates the zmax in a pos->z array with Nvidia's reduction
//	 * sample code (simplified and adapted for use with shape). This version
//	 * uses cudaStreams to process all frames concurrently.  	 */
//
//	int s;						// array size
//	int maxThreads = maxThreadsPerBlock;		// max # of threads per block
//	int maxBlocks = 2048;		// max # of blocks per grid
//	int whichKernel = 6;		// id of reduction kernel
//	int numBlocks = 0;			// initialize numBlocks
//	int numThreads = 0;			// initialize numThreads
//	float zmax = 0.0;			// radar cross section; return value
//	float *d_odata;				// temp. float array for reduction output
//	float *d_idata; 			// temp. float arrays for reduction input
//	float *h_odata;				// the host output array
//	float2 xblock_ythread;		// used for return value of getNumBlocksAndThreads
//
//	dim3 BLK,THD;
//	BLK.x = floor((maxThreadsPerBlock - 1 + size)/maxThreadsPerBlock);
//	THD.x = maxThreadsPerBlock; // Thread block dimensions
//
//	cudaStream_t zmax_stream[nframes];
//
//
//	/* Allocate memory for d_idata, then set that pointer equal to the right
//	 * data set and frame to the right deldop fit array	 */
//	cudaCalloc((void**)&d_idata, sizeof(float), size);
//
//	set_idata_zmax_krnl<<<BLK,THD>>>(ddat, d_idata, set, frm, size);
//	checkErrorAfterKernelLaunch("set_idata_zmax_krnl in compute_zmax");
//
//	/* Find number of blocks & threads needed for reduction call */
//	xblock_ythread = getNumBlocksAndThreads(size, maxBlocks, maxThreads);
//	numBlocks = xblock_ythread.x;
//	numThreads = xblock_ythread.y;
//
//	/* Allocate memory for d_odata and d_odata2 with enough elements to hold
//	 * the reduction of each block during the first call */
//	cudaCalloc((void**)&d_odata,  sizeof(float), numBlocks);
//	h_odata = (float *) malloc(numBlocks*sizeof(float));
//
//	/* Call maxz for first time */
//	maxz<float>(size, numThreads, numBlocks, whichKernel, d_idata, d_odata);
//	checkErrorAfterKernelLaunch("maxz<float> in compute_zmax");
//
//	/* Reset d_idata for later use as buffer */
//	cudaMemset(d_idata, 0, size*sizeof(float));
//
//	/* Now sum partial block sums on GPU, using the maxz6<> kernel */
//	s = numBlocks;
//
//	while (s > 1)
//	{
//		int threads = 0, blocks = 0;
//		xblock_ythread = getNumBlocksAndThreads(s, maxBlocks, maxThreads);
//		blocks = xblock_ythread.x;
//		threads = xblock_ythread.y;
//
//		/* Copy the first d_odata back into d_idata2  */
//		cudaMemcpy(d_idata, d_odata, s*sizeof(float), cudaMemcpyDeviceToDevice);
//
//		maxz<float>(s, threads, blocks, whichKernel, d_idata, d_odata);
//
//		if (whichKernel < 3)
//			s = (s + threads - 1) / threads;
//		else
//			s = (s + (threads*2-1)) / (threads*2);
//		if (s > 1)
//			printf("s is bigger than one");
//	}
//
//	gpuErrchk(cudaMemcpy(h_odata, d_odata, 2*sizeof(float), cudaMemcpyDeviceToHost));
//	zmax = h_odata[0];
//	free(h_odata);
//	cudaFree(d_odata);
//	cudaFree(d_idata);
//	return zmax;
//}
__host__ float compute_pos_zmax_all_frames(struct dat_t *ddat, int frame_size, int set, int nframes) {
	/* Function calculates the zmax per frame and then the final zmax for the set.
	 * Code assumes that frame_size is the same for all frames in set */

	int s;						// array size
	int maxBlocks = 2048;		// max # of blocks per grid
	int numBlocks = 0;			// initialize numBlocks
	int numThreads = 0;			// initialize numThreads
	float *zmax, final;			// radar cross section (per frame)
	float *d_odata;				// temp. float array for reduction output
	float *d_idata; 			// temp. float arrays for reduction input
	float2 xblock_ythread;		// used for return value of getNumBlocksAndThreads

	/* Find number of blocks & threads needed for reduction call */
	xblock_ythread = getNumBlocksAndThreads(frame_size, maxBlocks, maxThreadsPerBlock);
	numBlocks = xblock_ythread.x;
	numThreads = xblock_ythread.y;
	dim3 dimBlock(numThreads, 1, 1);
	dim3 dimGrid(numBlocks, 1, 1);

	/* Need to calculate zmax per frame, then multiply by frame weight and
	 * add to sum_deldop_zmax, which is again weighted and then returned.
	 * Copy each frame's pos->z_s into a double pointer or maybe split into subsets
	 * of the main total_size array?	 */

	/* Allocate memory for d_idata and d_odata */
	cudaCalloc((void**)&d_idata, sizeof(float), frame_size);
	cudaMallocManaged((void**)&d_odata, sizeof(float)*numBlocks, cudaMemAttachHost);
	cudaMallocManaged((void**)&zmax, sizeof(float)*nframes, cudaMemAttachHost);

	/* Configure data copy kernel launch */
	dim3 BLK,THD;
	BLK.x = floor((maxThreadsPerBlock - 1 + frame_size)/maxThreadsPerBlock);
	THD.x = maxThreadsPerBlock; // Thread block dimensions

	/* Start loop through frames */
	for (int f=0; f<nframes; f++) {

		/* Copy input data into input array */
		set_idata_zmax_all_frames_krnl<<<BLK,THD>>>(ddat, d_idata, set,
				f, frame_size);
		checkErrorAfterKernelLaunch("set_idata_zmax_all_frames_krnl in compute_zmax_all_frames");
		deviceSyncAfterKernelLaunch("");
		/* Start debug */
		//dbg_print_array(d_idata, 151, 151);
		/* End debug */

		/* Call reduction  */
		device_sum_block_atomic_kernel<<< dimGrid, dimBlock>>>(d_idata, d_odata, frame_size);
		checkErrorAfterKernelLaunch("device_sum_block_atomic_kernel (compute_pos_zmax_all_frames)");
		deviceSyncAfterKernelLaunch("");

		/* Copy zmax for this frame from the output array into the zmax array for the right frame */
		zmax[f] = d_odata[0];
	}
	/* Now apply frame weighting factors and sum up all frames to get sum_deldop_zmax to return */
	zmax_all_frames_finalize_krnl<<<1,1>>>(ddat,zmax,nframes,set);
	checkErrorAfterKernelLaunch("zmax_all_frames_finalize_krnl");
	final = zmax[0];
	cudaFree(d_odata);
	cudaFree(d_idata);
	return final;
}
__host__ float compute_pos_zmax_all_frames_2(struct dat_t *ddat, int size,
		int set, int nframes) {
	/* Function calculates the zmax in a pos->z arrary with
	 * Nvidia's reduction sample code (simplified and adapted for use with
	 * shape).  The function returns the cross section as a float 	 */

	int s;								// array size
	int maxThreads = maxThreadsPerBlock;// max # of threads per block
	int maxBlocks = 2048;				// max # of blocks per grid
	int whichKernel = 6;				// id of reduction kernel
	int numBlocks = 0;					// initialize numBlocks
	int numThreads = 0;					// initialize numThreads
	float *zmax, zmaxfinal = 0.0;		// max z values; return value
	float *d_odata;						// temp. float array for reduction output
	float *d_idata; 					// temp. float arrays for reduction input
	float *h_odata;						// the host output array
	float2 xblock_ythread;				// used for return value of getNumBlocksAndThreads

	dim3 BLK,THD;
	BLK.x = floor((maxThreadsPerBlock - 1 + size)/maxThreadsPerBlock);
	THD.x = maxThreadsPerBlock; // Thread block dimensions

	/* Allocate memory for d_idata, then set that pointer equal to the right
	 * data set and frame to the right deldop fit array	 */
	cudaCalloc((void**)&d_idata, sizeof(float), size);
	//cudaMallocManaged((void**)&zmax, sizeof(float)*nframes, cudaMemAttachHost);
	cudaCalloc((void**)&zmax, sizeof(float), nframes);
	/* Find number of blocks & threads needed for reduction call */
	xblock_ythread = getNumBlocksAndThreads(size, maxBlocks, maxThreads);
	numBlocks = xblock_ythread.x;
	numThreads = xblock_ythread.y;

	/* Allocate memory for d_odata and d_odata2 with enough elements to hold
	 * the reduction of each block during the first call */
	cudaCalloc((void**)&d_odata,  sizeof(float), numBlocks);
	h_odata = (float *) malloc(numBlocks*sizeof(float));

	for (int frm=0; frm<nframes; frm++) {
	set_idata_zmax_krnl<<<BLK,THD>>>(ddat, d_idata, set, frm, size);
	checkErrorAfterKernelLaunch("set_idata_zmax_krnl in compute_zmax");

	/* Call maxz for first time */
	maxz<float>(size, numThreads, numBlocks, whichKernel, d_idata, d_odata);
	checkErrorAfterKernelLaunch("maxz<float> in compute_zmax");

	/* Reset d_idata for later use as buffer */
	cudaMemset(d_idata, 0, size*sizeof(float));

	/* Now sum partial block sums on GPU, using the maxz6<> kernel */
	s = numBlocks;

	while (s > 1)
	{
		int threads = 0, blocks = 0;
		xblock_ythread = getNumBlocksAndThreads(s, maxBlocks, maxThreads);
		blocks = xblock_ythread.x;
		threads = xblock_ythread.y;

		/* Copy the first d_odata back into d_idata2  */
		cudaMemcpy(d_idata, d_odata, s*sizeof(float), cudaMemcpyDeviceToDevice);

		maxz<float>(s, threads, blocks, whichKernel, d_idata, d_odata);

		if (whichKernel < 3)
			s = (s + threads - 1) / threads;
		else
			s = (s + (threads*2-1)) / (threads*2);
		if (s > 1)
			printf("s is bigger than one");
	}

	gpuErrchk(cudaMemcpy(h_odata, d_odata, 2*sizeof(float), cudaMemcpyDeviceToHost));
	deviceSyncAfterKernelLaunch("");
	zmax[frm] = h_odata[0];

	} /*End frame loop*/
	/* Now apply frame weighting factors and sum up all frames to get sum_deldop_zmax to return */
	zmax_all_frames_finalize_krnl<<<1,1>>>(ddat,zmax,nframes,set);
	checkErrorAfterKernelLaunch("zmax_all_frames_finalize_krnl");
	deviceSyncAfterKernelLaunch("zmax_all_frames_finalize_krnl");
	zmaxfinal = zmax[0];
	free(h_odata);
	cudaFree(d_odata);
	cudaFree(d_idata);
	cudaFree(zmax);
	return zmaxfinal;
}
__host__ float compute_doppler_xsec(struct dat_t *ddat, int ndop,
		int set, int frm) {
	/* Function calculates a delay-Doppler frame's radar cross section with
	 * Nvidia's reduction sample code (simplified and adapted for use with
	 * shape).  The function returns the cross section as a float 	 */

	int s, size=ndop;			// array size
	int maxThreads = maxThreadsPerBlock;		// max # of threads per block
	int maxBlocks = 2048;		// max # of blocks per grid
	int whichKernel = 6;		// id of reduction kernel
	int numBlocks = 0;			// initialize numBlocks
	int numThreads = 0;			// initialize numThreads
	float xsec = 0.0;			// radar cross section; return value
	float *d_odata;				// temp. float array for reduction output
	float *d_idata; 			// temp. float arrays for reduction input
	float *h_odata;				// the host output array
	float2 xblock_ythread;		// used for return value of getNumBlocksAndThreads

	dim3 BLK,THD;
	BLK.x = floor((maxThreadsPerBlock - 1 + size)/maxThreadsPerBlock);
	THD.x = maxThreadsPerBlock; // Thread block dimensions

	/* Allocate memory for d_idata, then set that pointer equal to the right
	 * data set and frame to the right deldop fit array	 */
	cudaCalloc((void**)&d_idata, sizeof(float), size);

	set_idata_pntr_krnl<<<BLK,THD>>>(ddat, d_idata, set, frm, size);
	checkErrorAfterKernelLaunch("set_idata_pntr_krnl in compute_doppler_xsec");

	/* Find number of blocks & threads needed for reduction call */
	xblock_ythread = getNumBlocksAndThreads(size, maxBlocks, maxThreads);
	numBlocks = xblock_ythread.x;
	numThreads = xblock_ythread.y;

	/* Allocate memory for d_odata and d_odata2 with enough elements to hold
	 * the reduction of each block during the first call */
	cudaCalloc((void**)&d_odata,  sizeof(float), numBlocks);
	h_odata = (float *) malloc(numBlocks*sizeof(float));

	/* Call reduction for first time */
	reduce<float>(size, numThreads, numBlocks, whichKernel, d_idata, d_odata);
	checkErrorAfterKernelLaunch("reduce<float> in compute_deldop_xsec");

	/* Reset d_idata for later use as buffer */
	gpuErrchk(cudaMemset(d_idata, 0, size*sizeof(float)));

	/* Now sum partial block sums on GPU, using the reduce6<> kernel */
	s = numBlocks;

	while (s > 1)
	{
		int threads = 0, blocks = 0;
		xblock_ythread = getNumBlocksAndThreads(s, maxBlocks, maxThreads);
		blocks = xblock_ythread.x;
		threads = xblock_ythread.y;

		/* Copy the first d_odata back into d_idata2  */
		cudaMemcpy(d_idata, d_odata, s*sizeof(float), cudaMemcpyDeviceToDevice);

		reduce<float>(s, threads, blocks, whichKernel, d_idata, d_odata);

		if (whichKernel < 3)
			s = (s + threads - 1) / threads;
		else
			s = (s + (threads*2-1)) / (threads*2);
		if (s > 1)
			printf("s is bigger than one");
	}

	gpuErrchk(cudaMemcpy(h_odata, d_odata, 1*sizeof(float), cudaMemcpyDeviceToHost));
	xsec = h_odata[0];
	free(h_odata);
	cudaFree(d_odata);
	cudaFree(d_idata);
	return xsec;
}
__host__ float compute_model_area1(struct mod_t *dmod, int c, int size) {
	/* Function calculates a delay-Doppler frame's radar cross section with
	 * Nvidia's reduction sample code (simplified and adapted for use with
	 * shape).  The function returns the cross section as a float 	 */

	int s;						// array size
	int maxThreads = maxThreadsPerBlock;		// max # of threads per block
	int maxBlocks = 2048;		// max # of blocks per grid
	int whichKernel = 6;		// id of reduction kernel
	int numBlocks = 0;			// initialize numBlocks
	int numThreads = 0;			// initialize numThreads
	float area = 0.0;			// radar cross section; return value
	float *d_odata;				// temp. float array for reduction output
	float *d_idata; 			// temp. float arrays for reduction input
	float *h_odata;				// the host output array
	float2 xblock_ythread;		// used for return value of getNumBlocksAndThreads

	dim3 BLK,THD;
	BLK.x = floor((maxThreadsPerBlock - 1 + size)/maxThreadsPerBlock);
	THD.x = maxThreadsPerBlock; // Thread block dimensions

	/* Allocate memory for d_idata, then set that pointer equal to the right
	 * data set and frame to the right deldop fit array	 */
	cudaCalloc((void**)&d_idata, sizeof(float), size);

	set_idata_modarea_krnl<<<BLK,THD>>>(dmod, d_idata, c, size);
	checkErrorAfterKernelLaunch("set_idata_modarea_krnl in compute_doppler_xsec");

	/* Find number of blocks & threads needed for reduction call */
	xblock_ythread = getNumBlocksAndThreads(size, maxBlocks, maxThreads);
	numBlocks = xblock_ythread.x;
	numThreads = xblock_ythread.y;

	/* Allocate memory for d_odata and d_odata2 with enough elements to hold
	 * the reduction of each block during the first call */
	cudaCalloc((void**)&d_odata,  sizeof(float), numBlocks);
	h_odata = (float *) malloc(numBlocks*sizeof(float));

	/* Call reduction for first time */
	reduce<float>(size, numThreads, numBlocks, whichKernel, d_idata, d_odata);
	checkErrorAfterKernelLaunch("reduce<float> in compute_deldop_xsec");

	/* Reset d_idata for later use as buffer */
	gpuErrchk(cudaMemset(d_idata, 0, size*sizeof(float)));

	/* Now sum partial block sums on GPU, using the reduce6<> kernel */
	s = numBlocks;

	while (s > 1)
	{
		int threads = 0, blocks = 0;
		xblock_ythread = getNumBlocksAndThreads(s, maxBlocks, maxThreads);
		blocks = xblock_ythread.x;
		threads = xblock_ythread.y;

		/* Copy the first d_odata back into d_idata2  */
		gpuErrchk(cudaMemcpy(d_idata, d_odata, s*sizeof(float),
				cudaMemcpyDeviceToDevice));

		reduce<float>(s, threads, blocks, whichKernel, d_idata, d_odata);

		if (whichKernel < 3)
			s = (s + threads - 1) / threads;
		else
			s = (s + (threads*2-1)) / (threads*2);
		if (s > 1)
			printf("s is bigger than one");
	}

	gpuErrchk(cudaMemcpy(h_odata, d_odata, 1*sizeof(float), cudaMemcpyDeviceToHost));
	area = h_odata[0];
	free(h_odata);
	cudaFree(d_odata);
	cudaFree(d_idata);
	return area;
}
__host__ float compute_model_area(struct mod_t *dmod, int c, int size) {
	/* Function calculates a delay-Doppler frame's radar cross section with
	 * Nvidia's reduction sample code (simplified and adapted for use with
	 * shape).  The function returns the cross section as a float 	 */

	int maxThreads = maxThreadsPerBlock;		// max # of threads per block
	int maxBlocks = 2048;		// max # of blocks per grid
	int numBlocks = 0;			// initialize numBlocks
	int numThreads = 0;			// initialize numThreads
	float area = 0.0;			// radar cross section; return value
	float *d_odata;				// temp. float array for reduction output
	float *d_idata; 			// temp. float arrays for reduction input
	float2 xblock_ythread;		// used for return value of getNumBlocksAndThreads

	dim3 BLK,THD;
	BLK.x = floor((maxThreadsPerBlock - 1 + size)/maxThreadsPerBlock);
	THD.x = maxThreadsPerBlock; // Thread block dimensions

	/* Find number of blocks & threads needed for reduction call */
	xblock_ythread = getNumBlocksAndThreads(size, maxBlocks, maxThreads);
	numBlocks = xblock_ythread.x;
	numThreads = xblock_ythread.y;
	dim3 dimBlock(numThreads, 1, 1);
	dim3 dimGrid(numBlocks, 1, 1);

	/* Allocate memory for d_idata and d_odata */
	cudaCalloc((void**)&d_idata, sizeof(float), size);
	cudaCalloc((void**)&d_odata,  sizeof(float), numBlocks);

	/* Load the d_idata array */
	set_idata_modarea_krnl<<<BLK,THD>>>(dmod, d_idata, c, size);
	checkErrorAfterKernelLaunch("set_idata_modarea_krnl in compute_doppler_xsec");

	/* Call reduction  */
	device_reduce_block_atomic_kernel<<< dimGrid, dimBlock>>>(d_idata, d_odata, size);
	checkErrorAfterKernelLaunch("device_reduce_block_atomic_kernel (compute_model_area)");
	deviceSyncAfterKernelLaunch("device_reduce_block_atomic_kernel (compute_model_area)");

	area = d_odata[0];
	cudaFree(d_odata);
	cudaFree(d_idata);
	return area;
}
__host__ void dvdI_reduce_single(struct mod_t *dmod, float *dv, float *dcom0,
		float *dcom1, float *dcom2, float *dI00, float *dI01, float *dI02,
		float *dI10, float *dI11, float *dI12, float *dI20, float *dI21,
		float *dI22, int size, int c)
{
	int maxThreads = maxThreadsPerBlock;		// max # of threads per block
	int maxBlocks = 2048;		// max # of blocks per grid
	int numBlocks = 0;			// initialize numBlocks
	int numThreads = 0;			// initialize numThreads
	float2 xblock_ythread;		// used for return value of getNumBlocksAndThreads
	/* Device output arrays */
	float *d_odata_dcom0, *d_odata_dcom1, *d_odata_dcom2, *d_odata_dI00,
	*d_odata_dI01, *d_odata_dI02, *d_odata_dI10, *d_odata_dI11,
	*d_odata_dI12, *d_odata_dI20, *d_odata_dI21, *d_odata_dI22,
	*d_odata_dv;

	/* Find number of blocks & threads needed for reduction call */
	xblock_ythread = getNumBlocksAndThreads(size, maxBlocks, maxThreads);
	numBlocks = xblock_ythread.x;
	numThreads = xblock_ythread.y;
	dim3 dimBlock(numThreads, 1, 1);
	dim3 dimGrid(numBlocks, 1, 1);

	/* Allocate memory for the device output arrays for first reduction */
	cudaCalloc((void**)&d_odata_dv,    	 sizeof(float), numBlocks);
	cudaCalloc((void**)&d_odata_dcom0, 	 sizeof(float), numBlocks);
	cudaCalloc((void**)&d_odata_dcom1, 	 sizeof(float), numBlocks);
	cudaCalloc((void**)&d_odata_dcom2, 	 sizeof(float), numBlocks);
	cudaCalloc((void**)&d_odata_dI00, 	 sizeof(float), numBlocks);
	cudaCalloc((void**)&d_odata_dI01, 	 sizeof(float), numBlocks);
	cudaCalloc((void**)&d_odata_dI02, 	 sizeof(float), numBlocks);
	cudaCalloc((void**)&d_odata_dI10, 	 sizeof(float), numBlocks);
	cudaCalloc((void**)&d_odata_dI11, 	 sizeof(float), numBlocks);
	cudaCalloc((void**)&d_odata_dI12, 	 sizeof(float), numBlocks);
	cudaCalloc((void**)&d_odata_dI20, 	 sizeof(float), numBlocks);
	cudaCalloc((void**)&d_odata_dI21, 	 sizeof(float), numBlocks);
	cudaCalloc((void**)&d_odata_dI22, 	 sizeof(float), numBlocks);

	cudaStream_t stream[13];
	for (int i=0; i<13; i++)
		gpuErrchk(cudaStreamCreate(&stream[i]));


	device_reduce_block_atomic_kernel<<< dimGrid, dimBlock, 0, stream[0] >>>(dv, d_odata_dv, size);
	checkErrorAfterKernelLaunch("device_reduce_block_atomic_kernel (dv)");
	device_reduce_block_atomic_kernel<<< dimGrid, dimBlock, 0, stream[1] >>>(dcom0, d_odata_dcom0, size);
	checkErrorAfterKernelLaunch("device_reduce_block_atomic_kernel (dcom0)");
	device_reduce_block_atomic_kernel<<< dimGrid, dimBlock, 0, stream[2] >>>(dcom1, d_odata_dcom1, size);
	checkErrorAfterKernelLaunch("device_reduce_block_atomic_kernel (dcom1)");
	device_reduce_block_atomic_kernel<<< dimGrid, dimBlock, 0, stream[3] >>>(dcom2, d_odata_dcom2, size);
	checkErrorAfterKernelLaunch("device_reduce_block_atomic_kernel (dcom2)");
	device_reduce_block_atomic_kernel<<< dimGrid, dimBlock, 0, stream[4] >>>(dI00, d_odata_dI00, size);
	checkErrorAfterKernelLaunch("device_reduce_block_atomic_kernel (dI00)");
	device_reduce_block_atomic_kernel<<< dimGrid, dimBlock, 0, stream[5] >>>(dI01, d_odata_dI01, size);
	checkErrorAfterKernelLaunch("device_reduce_block_atomic_kernel (dI01)");
	device_reduce_block_atomic_kernel<<< dimGrid, dimBlock, 0, stream[6] >>>(dI02, d_odata_dI02, size);
	checkErrorAfterKernelLaunch("device_reduce_block_atomic_kernel (dI02)");
	device_reduce_block_atomic_kernel<<< dimGrid, dimBlock, 0, stream[7] >>>(dI10, d_odata_dI10, size);
	checkErrorAfterKernelLaunch("device_reduce_block_atomic_kernel (dI10)");
	device_reduce_block_atomic_kernel<<< dimGrid, dimBlock, 0, stream[8] >>>(dI11, d_odata_dI11, size);
	checkErrorAfterKernelLaunch("device_reduce_block_atomic_kernel (dI11)");
	device_reduce_block_atomic_kernel<<< dimGrid, dimBlock, 0, stream[9] >>>(dI12, d_odata_dI12, size);
	checkErrorAfterKernelLaunch("device_reduce_block_atomic_kernel (dI12)");
	device_reduce_block_atomic_kernel<<< dimGrid, dimBlock, 0, stream[10] >>>(dI20, d_odata_dI20, size);
	checkErrorAfterKernelLaunch("device_reduce_block_atomic_kernel (dI20)");
	device_reduce_block_atomic_kernel<<< dimGrid, dimBlock, 0, stream[11] >>>(dI21, d_odata_dI21, size);
	checkErrorAfterKernelLaunch("device_reduce_block_atomic_kernel (dI21)");
	device_reduce_block_atomic_kernel<<< dimGrid, dimBlock, 0, stream[12] >>>(dI22, d_odata_dI22, size);
	checkErrorAfterKernelLaunch("device_reduce_block_atomic_kernel (dI22)");

	for (int i=0; i<13; i++)
		gpuErrchk(cudaStreamSynchronize(stream[i]));

	/* Copy and assign */
	set_dv_dcom_di_krnl<<<1,1>>>(d_odata_dv, d_odata_dcom0, d_odata_dcom1,
			d_odata_dcom2, d_odata_dI00, d_odata_dI01, d_odata_dI02, d_odata_dI10,
			d_odata_dI11, d_odata_dI12, d_odata_dI20, d_odata_dI21, d_odata_dI22,
			dmod, c);

	/* Free up the temporary arrays */
	cudaFree(d_odata_dv);
	cudaFree(d_odata_dcom0);	cudaFree(d_odata_dcom1);	cudaFree(d_odata_dcom2);
	cudaFree(d_odata_dI00);		cudaFree(d_odata_dI01);		cudaFree(d_odata_dI02);
	cudaFree(d_odata_dI10);		cudaFree(d_odata_dI11);		cudaFree(d_odata_dI12);
	cudaFree(d_odata_dI20);		cudaFree(d_odata_dI21);		cudaFree(d_odata_dI22);

	for (int i=0; i<13; i++)
		gpuErrchk(cudaStreamDestroy(stream[i]));

}
__host__ void compute_dv_dcom_dI_reduction(float *dv, float *dcom0, float
		*dcom1, float *dcom2, float *dI00, float *dI01, float *dI02, float
		*dI10, float *dI11, float *dI12, float *dI20, float *dI21, float *dI22,
		int c, int size, struct mod_t *dmod) {
	/* Function calculates the model's COM and Inertia tensors	 */

	int s;						// array size
	int maxThreads = maxThreadsPerBlock;		// max # of threads per block
	int maxBlocks = 2048;		// max # of blocks per grid
	int whichKernel = 6;		// id of reduction kernel
	int numBlocks = 0;			// initialize numBlocks
	int numThreads = 0;			// initialize numThreads

	/* Device output arrays */
	float *d_odata_dcom0, *d_odata_dcom1, *d_odata_dcom2, *d_odata_dI00,
	*d_odata_dI01, *d_odata_dI02, *d_odata_dI10, *d_odata_dI11,
	*d_odata_dI12, *d_odata_dI20, *d_odata_dI21, *d_odata_dI22,
	*d_odata_dv;

	float2 xblock_ythread;		// used for return value of getNumBlocksAndThreads
	//
	//	dim3 BLK,THD;
	//	BLK.x = floor((maxThreadsPerBlock - 1 + size)/maxThreadsPerBlock);
	//	THD.x = maxThreadsPerBlock; // Thread block dimensions

	/* Find number of blocks & threads needed for reduction call */
	xblock_ythread = getNumBlocksAndThreads(size, maxBlocks, maxThreads);
	numBlocks = xblock_ythread.x;
	numThreads = xblock_ythread.y;

	/* Allocate memory for the device output arrays for first reduction */
	cudaCalloc((void**)&d_odata_dv,    	 sizeof(float), numBlocks);
	cudaCalloc((void**)&d_odata_dcom0, 	 sizeof(float), numBlocks);
	cudaCalloc((void**)&d_odata_dcom1, 	 sizeof(float), numBlocks);
	cudaCalloc((void**)&d_odata_dcom2, 	 sizeof(float), numBlocks);
	cudaCalloc((void**)&d_odata_dI00, 	 sizeof(float), numBlocks);
	cudaCalloc((void**)&d_odata_dI01, 	 sizeof(float), numBlocks);
	cudaCalloc((void**)&d_odata_dI02, 	 sizeof(float), numBlocks);
	cudaCalloc((void**)&d_odata_dI10, 	 sizeof(float), numBlocks);
	cudaCalloc((void**)&d_odata_dI11, 	 sizeof(float), numBlocks);
	cudaCalloc((void**)&d_odata_dI12, 	 sizeof(float), numBlocks);
	cudaCalloc((void**)&d_odata_dI20, 	 sizeof(float), numBlocks);
	cudaCalloc((void**)&d_odata_dI21, 	 sizeof(float), numBlocks);
	cudaCalloc((void**)&d_odata_dI22, 	 sizeof(float), numBlocks);

	/* Call reductions for first time */
	if (STREAMS) {
	reducedI<float>(size, numThreads, numBlocks, whichKernel, dv, dcom0, dcom1, dcom2,
			dI00, dI01, dI02, dI10, dI11, dI12, dI20, dI21, dI22, d_odata_dv,
			d_odata_dcom0, d_odata_dcom1, d_odata_dcom2, d_odata_dI00, d_odata_dI01,
			d_odata_dI02, d_odata_dI10, d_odata_dI11, d_odata_dI12, d_odata_dI20,
			d_odata_dI21, d_odata_dI22);
	checkErrorAfterKernelLaunch("reducedI<float> in compute_dv_dcom_dI_reduction"); }

	else {
	reduce<float>(size, numThreads, numBlocks, whichKernel, dv, d_odata_dv);
	checkErrorAfterKernelLaunch("reduce<float> in compute_dv_dcom_dI_reduction");
	reduce<float>(size, numThreads, numBlocks, whichKernel, dcom0, d_odata_dcom0);
	checkErrorAfterKernelLaunch("reduce<float> in compute_dv_dcom_dI_reduction");
	reduce<float>(size, numThreads, numBlocks, whichKernel, dcom1, d_odata_dcom1);
	checkErrorAfterKernelLaunch("reduce<float> in compute_dv_dcom_dI_reduction");
	reduce<float>(size, numThreads, numBlocks, whichKernel, dcom2, d_odata_dcom2);
	checkErrorAfterKernelLaunch("reduce<float> in compute_dv_dcom_dI_reduction");
	reduce<float>(size, numThreads, numBlocks, whichKernel, dI00, d_odata_dI00);
	checkErrorAfterKernelLaunch("reduce<float> in compute_dv_dcom_dI_reduction");
	reduce<float>(size, numThreads, numBlocks, whichKernel, dI01, d_odata_dI01);
	checkErrorAfterKernelLaunch("reduce<float> in compute_dv_dcom_dI_reduction");
	reduce<float>(size, numThreads, numBlocks, whichKernel, dI02, d_odata_dI02);
	checkErrorAfterKernelLaunch("reduce<float> in compute_dv_dcom_dI_reduction");
	reduce<float>(size, numThreads, numBlocks, whichKernel, dI10, d_odata_dI10);
	checkErrorAfterKernelLaunch("reduce<float> in compute_dv_dcom_dI_reduction");
	reduce<float>(size, numThreads, numBlocks, whichKernel, dI11, d_odata_dI11);
	checkErrorAfterKernelLaunch("reduce<float> in compute_dv_dcom_dI_reduction");
	reduce<float>(size, numThreads, numBlocks, whichKernel, dI12, d_odata_dI12);
	checkErrorAfterKernelLaunch("reduce<float> in compute_dv_dcom_dI_reduction");
	reduce<float>(size, numThreads, numBlocks, whichKernel, dI20, d_odata_dI20);
	checkErrorAfterKernelLaunch("reduce<float> in compute_dv_dcom_dI_reduction");
	reduce<float>(size, numThreads, numBlocks, whichKernel, dI21, d_odata_dI21);
	checkErrorAfterKernelLaunch("reduce<float> in compute_dv_dcom_dI_reduction");
	reduce<float>(size, numThreads, numBlocks, whichKernel, dI22, d_odata_dI22);
	checkErrorAfterKernelLaunch("reduce<float> in compute_dv_dcom_dI_reduction");
	}

	/* Reset the orig. input arrays for later use as buffer */
	gpuErrchk(cudaMemset(dv, 	  0, size*sizeof(float)));
	gpuErrchk(cudaMemset(dcom0, 0, size*sizeof(float)));
	gpuErrchk(cudaMemset(dcom1, 0, size*sizeof(float)));
	gpuErrchk(cudaMemset(dcom2, 0, size*sizeof(float)));
	gpuErrchk(cudaMemset(dI00,  0, size*sizeof(float)));
	gpuErrchk(cudaMemset(dI01,  0, size*sizeof(float)));
	gpuErrchk(cudaMemset(dI02,  0, size*sizeof(float)));
	gpuErrchk(cudaMemset(dI10,  0, size*sizeof(float)));
	gpuErrchk(cudaMemset(dI11,  0, size*sizeof(float)));
	gpuErrchk(cudaMemset(dI12,  0, size*sizeof(float)));
	gpuErrchk(cudaMemset(dI20,  0, size*sizeof(float)));
	gpuErrchk(cudaMemset(dI21,  0, size*sizeof(float)));
	gpuErrchk(cudaMemset(dI22,  0, size*sizeof(float)));

	/* Now sum partial block sums on GPU, using the reduce6<> kernel */
	s = numBlocks;

	while (s > 1)
	{
		int threads = 0, blocks = 0;
		xblock_ythread = getNumBlocksAndThreads(s, maxBlocks, maxThreads);
		blocks = xblock_ythread.x;
		threads = xblock_ythread.y;

		/* Copy the d_odata_xx arrays back into the zeroed-out input arrays */
		gpuErrchk(cudaMemcpy(dv, 	d_odata_dv,    s*sizeof(float), cudaMemcpyDeviceToDevice));
		gpuErrchk(cudaMemcpy(dcom0, d_odata_dcom0, s*sizeof(float), cudaMemcpyDeviceToDevice));
		gpuErrchk(cudaMemcpy(dcom1, d_odata_dcom1, s*sizeof(float), cudaMemcpyDeviceToDevice));
		gpuErrchk(cudaMemcpy(dcom2, d_odata_dcom2, s*sizeof(float), cudaMemcpyDeviceToDevice));
		gpuErrchk(cudaMemcpy(dI00,  d_odata_dI00,  s*sizeof(float), cudaMemcpyDeviceToDevice));
		gpuErrchk(cudaMemcpy(dI01,  d_odata_dI01,  s*sizeof(float), cudaMemcpyDeviceToDevice));
		gpuErrchk(cudaMemcpy(dI02,  d_odata_dI02,  s*sizeof(float), cudaMemcpyDeviceToDevice));
		gpuErrchk(cudaMemcpy(dI10,  d_odata_dI10,  s*sizeof(float), cudaMemcpyDeviceToDevice));
		gpuErrchk(cudaMemcpy(dI11,  d_odata_dI11,  s*sizeof(float), cudaMemcpyDeviceToDevice));
		gpuErrchk(cudaMemcpy(dI12,  d_odata_dI12,  s*sizeof(float), cudaMemcpyDeviceToDevice));
		gpuErrchk(cudaMemcpy(dI20,  d_odata_dI20,  s*sizeof(float), cudaMemcpyDeviceToDevice));
		gpuErrchk(cudaMemcpy(dI21,  d_odata_dI21,  s*sizeof(float), cudaMemcpyDeviceToDevice));
		gpuErrchk(cudaMemcpy(dI22,  d_odata_dI22,  s*sizeof(float), cudaMemcpyDeviceToDevice));

		/* Call all reductions again until s = 1 */
		if (STREAMS) {
			reducedI<float>(s, threads, blocks, whichKernel, dv, dcom0, dcom1, dcom2,
					dI00, dI01, dI02, dI10, dI11, dI12, dI20, dI21, dI22, d_odata_dv,
					d_odata_dcom0, d_odata_dcom1, d_odata_dcom2, d_odata_dI00, d_odata_dI01,
					d_odata_dI02, d_odata_dI10, d_odata_dI11, d_odata_dI12, d_odata_dI20,
					d_odata_dI21, d_odata_dI22);
			checkErrorAfterKernelLaunch("reducedI<float> in compute_dv_dcom_dI_reduction");
		}
		else {
		reduce<float>(s, threads, blocks, whichKernel, dv, d_odata_dv);
		checkErrorAfterKernelLaunch("reduce<float> in compute_dv_dcom_dI_reduction");
		reduce<float>(s, threads, blocks, whichKernel, dcom0, d_odata_dcom0);
		checkErrorAfterKernelLaunch("reduce<float> in compute_dv_dcom_dI_reduction");
		reduce<float>(s, threads, blocks, whichKernel, dcom1, d_odata_dcom1);
		checkErrorAfterKernelLaunch("reduce<float> in compute_dv_dcom_dI_reduction");
		reduce<float>(s, threads, blocks, whichKernel, dcom2, d_odata_dcom2);
		checkErrorAfterKernelLaunch("reduce<float> in compute_dv_dcom_dI_reduction");
		reduce<float>(s, threads, blocks, whichKernel, dI00,  d_odata_dI00);
		checkErrorAfterKernelLaunch("reduce<float> in compute_dv_dcom_dI_reduction");
		reduce<float>(s, threads, blocks, whichKernel, dI01,  d_odata_dI01);
		checkErrorAfterKernelLaunch("reduce<float> in compute_dv_dcom_dI_reduction");
		reduce<float>(s, threads, blocks, whichKernel, dI02,  d_odata_dI02);
		checkErrorAfterKernelLaunch("reduce<float> in compute_dv_dcom_dI_reduction");
		reduce<float>(s, threads, blocks, whichKernel, dI10,  d_odata_dI10);
		checkErrorAfterKernelLaunch("reduce<float> in compute_dv_dcom_dI_reduction");
		reduce<float>(s, threads, blocks, whichKernel, dI11,  d_odata_dI11);
		checkErrorAfterKernelLaunch("reduce<float> in compute_dv_dcom_dI_reduction");
		reduce<float>(s, threads, blocks, whichKernel, dI12,  d_odata_dI12);
		checkErrorAfterKernelLaunch("reduce<float> in compute_dv_dcom_dI_reduction");
		reduce<float>(s, threads, blocks, whichKernel, dI20,  d_odata_dI20);
		checkErrorAfterKernelLaunch("reduce<float> in compute_dv_dcom_dI_reduction");
		reduce<float>(s, threads, blocks, whichKernel, dI21,  d_odata_dI21);
		checkErrorAfterKernelLaunch("reduce<float> in compute_dv_dcom_dI_reduction");
		reduce<float>(s, threads, blocks, whichKernel, dI22,  d_odata_dI22);
		checkErrorAfterKernelLaunch("reduce<float> in compute_dv_dcom_dI_reduction");
		}

		if (whichKernel < 3)
			s = (s + threads - 1) / threads;
		else
			s = (s + (threads*2-1)) / (threads*2);
		if (s > 1)
			printf("s is bigger than one");
	}

	/* Copy and assign */
	set_dv_dcom_di_krnl<<<1,1>>>(d_odata_dv, d_odata_dcom0, d_odata_dcom1,
			d_odata_dcom2, d_odata_dI00, d_odata_dI01, d_odata_dI02, d_odata_dI10,
			d_odata_dI11, d_odata_dI12, d_odata_dI20, d_odata_dI21, d_odata_dI22,
			dmod, c);

	/* Free up the temporary arrays */
	cudaFree(d_odata_dv);
	cudaFree(d_odata_dcom0);	cudaFree(d_odata_dcom1);	cudaFree(d_odata_dcom2);
	cudaFree(d_odata_dI00);		cudaFree(d_odata_dI01);		cudaFree(d_odata_dI02);
	cudaFree(d_odata_dI10);		cudaFree(d_odata_dI11);		cudaFree(d_odata_dI12);
	cudaFree(d_odata_dI20);		cudaFree(d_odata_dI21);		cudaFree(d_odata_dI22);
}

__global__ void dbg_pos_xlim_krnl(float *debug, float *d_odata_imax,
		float *d_odata_jmax, float *d_odata_imin, float *d_odata_jmin) {
	/* Single-threaded debug kernel */
	if (threadIdx.x == 0) {
		debug[0] = d_odata_imin[0];
		debug[1] = d_odata_imax[0];
		debug[2] = d_odata_jmin[0];
		debug[3] = d_odata_jmax[0];
		printf("\nimin: %g", debug[0]);
		printf("\nimax: %g", debug[1]);
		printf("\njmin: %g", debug[2]);
		printf("\njmax: %g", debug[3]);
	}
}
__global__ void dbg_print_device_array(float *in, int size) {
	/* Single-threaded debug kernel */
	int i;
	if (threadIdx.x == 0) {
		for (i=0; i<size; i++)
			printf("\ndev_array[%i]=%g", i, in[i]);
	}
}

__host__ void compute_xlim_ylim(struct dat_t *ddat, int size,
		int set, int frm, int src, float *iminflt, float *imaxflt, float *jminflt,
		float *jmaxflt, float *minmax_overall) {
	/* Function calculates the pos->xlim and pos->ylim values and also the
	 * imax_overall/imin_overall/jmax_overall/jmin_overall values and updates
	 * pos accordingly 	 */

	int s;						// array size
	int maxThreads = 256;		// max # of threads per block
	int maxBlocks = 256;		// max # of blocks per grid
	int whichKernel = 6;		// id of reduction kernel
	int numBlocks = 0;			// initialize numBlocks
	int numThreads = 0;			// initialize numThreads
	float *d_odata_imax;		// temp. float arrays for reduction output
	float *d_odata_imin;
	float *d_odata_jmax;
	float *d_odata_jmin;
	float2 xblock_ythread;		// used for return value of getNumBlocksAndThreads

	dim3 BLK,THD;
	BLK.x = floor((maxThreadsPerBlock - 1 + size)/maxThreadsPerBlock);
	THD.x = maxThreadsPerBlock; // Thread block dimensions

	/* Find number of blocks & threads needed for reduction call */
	xblock_ythread = getNumBlocksAndThreads(size, maxBlocks, maxThreads);
	numBlocks = xblock_ythread.x;
	numThreads = xblock_ythread.y;

	/* Create the streams (four) */
//	cudaStream_t stream1, stream2, stream3, stream4;
//	cudaStreamCreate(&stream1);
//	cudaStreamCreate(&stream2);
//	cudaStreamCreate(&stream3);
//	cudaStreamCreate(&stream4);

	/* Allocate memory for four device output data arrays d_odata_imax,
	 * d_odata_imin, d_odata_jmax, d_odata_jminwith enough elements to hold
	 * the reduction of each block during the first call */
	cudaCalloc((void**)&d_odata_imax,  sizeof(float), numBlocks);
	cudaCalloc((void**)&d_odata_imin,  sizeof(float), numBlocks);
	cudaCalloc((void**)&d_odata_jmax,  sizeof(float), numBlocks);
	cudaCalloc((void**)&d_odata_jmin,  sizeof(float), numBlocks);
//	gpuErrchk(cudaHostAlloc((void**)&d_odata_imax, sizeof(float)*4,
//			cudaHostAllocWriteCombined | cudaHostAllocMapped));
//	gpuErrchk(cudaHostAlloc((void**)&d_odata_imin, sizeof(float)*4,
//			cudaHostAllocWriteCombined | cudaHostAllocMapped));
//	gpuErrchk(cudaHostAlloc((void**)&d_odata_jmax, sizeof(float)*4,
//			cudaHostAllocWriteCombined | cudaHostAllocMapped));
//	gpuErrchk(cudaHostAlloc((void**)&d_odata_jmin, sizeof(float)*4,
//			cudaHostAllocWriteCombined | cudaHostAllocMapped));


	/* Call maxz for first time */
	maxz<float>(size, numThreads, numBlocks, whichKernel, imaxflt, d_odata_imax);
	checkErrorAfterKernelLaunch("maxz<float> for imaxflt in compute_xlim_ylim");

	maxz<float>(size, numThreads, numBlocks, whichKernel, jmaxflt, d_odata_jmax);
	checkErrorAfterKernelLaunch("maxz<float> for jmaxflt in compute_xlim_ylim");

	minz<float>(size, numThreads, numBlocks, whichKernel, iminflt, d_odata_imin);
	checkErrorAfterKernelLaunch("minz<float> for iminflt in compute_xlim_ylim");

	minz<float>(size, numThreads, numBlocks, whichKernel, jminflt, d_odata_jmin);
	checkErrorAfterKernelLaunch("minz<float> for jminflt in compute_xlim_ylim");
//	maxzexp<float>(size, numThreads, numBlocks, whichKernel, imaxflt, jmaxflt,
//			iminflt, jminflt, d_odata_imax, d_odata_jmax, d_odata_imin,
//			d_odata_jmin, &stream1, &stream2, &stream3, &stream4);

	/* Reset the original input arrays for later use as buffer */
	gpuErrchk(cudaMemset(imaxflt, 0, size*sizeof(float)));
	gpuErrchk(cudaMemset(iminflt, 0, size*sizeof(float)));
	gpuErrchk(cudaMemset(jmaxflt, 0, size*sizeof(float)));
	gpuErrchk(cudaMemset(jminflt, 0, size*sizeof(float)));

	/* Now sum partial block sums on GPU, using the maxz6<> kernel */
	s = numBlocks;

	while (s > 1)
	{
		int threads = 0, blocks = 0;
		xblock_ythread = getNumBlocksAndThreads(s, maxBlocks, maxThreads);
		blocks = xblock_ythread.x;
		threads = xblock_ythread.y;

		/* Copy the first d_odata back into d_idata2  */
		gpuErrchk(cudaMemcpy(imaxflt, d_odata_imax, s*sizeof(float), cudaMemcpyDeviceToHost));
		gpuErrchk(cudaMemcpy(iminflt, d_odata_imin, s*sizeof(float), cudaMemcpyDeviceToHost));
		gpuErrchk(cudaMemcpy(jmaxflt, d_odata_jmax, s*sizeof(float), cudaMemcpyDeviceToHost));
		gpuErrchk(cudaMemcpy(jminflt, d_odata_jmin, s*sizeof(float), cudaMemcpyDeviceToHost));

		maxz<float>(s, threads, blocks, whichKernel, imaxflt, d_odata_imax);
		checkErrorAfterKernelLaunch("maxz<float> for imaxflt in compute_xlim_ylim");
		maxz<float>(s, threads, blocks, whichKernel, jmaxflt, d_odata_jmax);
		checkErrorAfterKernelLaunch("maxz<float> for jmaxflt in compute_xlim_ylim");
		minz<float>(s, threads, blocks, whichKernel, iminflt, d_odata_imin);
		checkErrorAfterKernelLaunch("minz<float> for iminflt in compute_xlim_ylim");
		minz<float>(s, threads, blocks, whichKernel, jminflt, d_odata_jmin);
		checkErrorAfterKernelLaunch("minz<float> for jminflt in compute_xlim_ylim");
//		maxzexp<float>(s, threads, blocks, whichKernel, imaxflt, jmaxflt, iminflt,
//				jminflt, d_odata_imax, d_odata_jmax, d_odata_imin, d_odata_jmin,
//				&stream1, &stream2, &stream3, &stream4);

		if (whichKernel < 3)
			s = (s + threads - 1) / threads;
		else
			s = (s + (threads*2-1)) / (threads*2);
		if (s > 1)
			printf("s is bigger than one");
	}

	/* Sync streams */
//	gpuErrchk(cudaStreamSynchronize(stream1));
//	gpuErrchk(cudaStreamSynchronize(stream2));
//	gpuErrchk(cudaStreamSynchronize(stream3));
//	gpuErrchk(cudaStreamSynchronize(stream4));

	/* Calculate the min/max overall values (regardless of POS frame limits) */
	set_xlim_ylim_krnl<<<1,1>>>(ddat, set, frm, src, d_odata_imax, d_odata_imin,
			d_odata_jmax, d_odata_jmin, minmax_overall);
	checkErrorAfterKernelLaunch("set_xlim_ylim_krnl in compute_xlim_ylim");


	/* Nuke streams */
//	gpuErrchk(cudaStreamDestroy(stream1));
//	gpuErrchk(cudaStreamDestroy(stream2));
//	gpuErrchk(cudaStreamDestroy(stream3));
//	gpuErrchk(cudaStreamDestroy(stream4));
	cudaFree(d_odata_imax);
	cudaFree(d_odata_imin);
	cudaFree(d_odata_jmax);
	cudaFree(d_odata_jmin);
}

__host__ void c2af_deldop_add_o2_m2(
		float **temp_o2,
		float **temp_m2,
		float **temp_om,
		int size,
		int nframes) {
	/* Function reduces the input arrays for nframes-frames, once per frame.
	 * The input array is structured like this:  input[nframes][size]  */

	int maxThreads = maxThreadsPerBlock;
	int maxBlocks = 2048;
	int numBlocks = 0;			// initialize numBlocks
	int numThreads = 0;			// initialize numThreads
	float *d_odata_o2, *d_odata_m2, *d_odata_om;
	float *d_idata_o2, *d_idata_m2, *d_idata_om;
	float2 xblock_ythread;

	dim3 BLK,THD;
	BLK.x = floor((maxThreadsPerBlock - 1 + size)/maxThreadsPerBlock);
	THD.x = maxThreadsPerBlock; // Thread block dimensions

	/* Find number of blocks & threads needed to reduce ONE FRAME ONLY */
	xblock_ythread = getNumBlocksAndThreads(size, maxBlocks, maxThreads);
	numBlocks = xblock_ythread.x;
	numThreads = xblock_ythread.y;
	dim3 dimBlock(numThreads, 1, 1);
	dim3 dimGrid(numBlocks, 1, 1);

	/* Allocate memory for d_idata and d_odata */
	cudaCalloc((void**)&d_idata_o2, sizeof(float), size);
	cudaCalloc((void**)&d_odata_o2,  sizeof(float), numBlocks);
	cudaCalloc((void**)&d_idata_m2, sizeof(float), size);
	cudaCalloc((void**)&d_odata_m2,  sizeof(float), numBlocks);
	cudaCalloc((void**)&d_idata_om, sizeof(float), size);
	cudaCalloc((void**)&d_odata_om,  sizeof(float), numBlocks);

	for (int frm=0; frm<nframes; frm++) {
		c2af_set_data_krnl<<<BLK,THD>>>(temp_o2, d_idata_o2, frm, size);
		checkErrorAfterKernelLaunch("c2af_set_data_krnl in reduction.cu");

		c2af_set_data_krnl<<<BLK,THD>>>(temp_m2, d_idata_m2, frm, size);
		checkErrorAfterKernelLaunch("c2af_set_data_krnl in reduction.cu");

		c2af_set_data_krnl<<<BLK,THD>>>(temp_om, d_idata_om, frm, size);
		checkErrorAfterKernelLaunch("c2af_set_data_krnl in reduction.cu");

		/* Call reduction  */
		device_reduce_block_atomic_kernel<<< dimGrid, dimBlock>>>(d_idata_o2,
				d_odata_o2, size);
		checkErrorAfterKernelLaunch("device_reduce_block_atomic_kernel");

		device_reduce_block_atomic_kernel<<< dimGrid, dimBlock>>>(d_idata_m2,
				d_odata_m2, size);
		checkErrorAfterKernelLaunch("device_reduce_block_atomic_kernel");

		device_reduce_block_atomic_kernel<<< dimGrid, dimBlock>>>(d_idata_om,
				d_odata_om, size);
		checkErrorAfterKernelLaunch("device_reduce_block_atomic_kernel");

		deviceSyncAfterKernelLaunch("device_reduce_block_atomic_kernel");

		temp_o2[frm][0] = d_odata_o2[0];
		temp_m2[frm][0] = d_odata_m2[0];
		temp_om[frm][0] = d_odata_om[0];

		gpuErrchk(cudaMemset(d_odata_o2, 0, numBlocks*sizeof(float)));
		gpuErrchk(cudaMemset(d_odata_m2, 0, numBlocks*sizeof(float)));
		gpuErrchk(cudaMemset(d_odata_om, 0, numBlocks*sizeof(float)));

	}
	/* Output sum for each frame is in first entry for each frame in the
	 * input array	 */

	cudaFree(d_odata_o2);
	cudaFree(d_idata_o2);
	cudaFree(d_odata_m2);
	cudaFree(d_idata_m2);
	cudaFree(d_odata_om);
	cudaFree(d_idata_om);

}
//#endif // #ifndef _REDUCE_KERNEL_H_

