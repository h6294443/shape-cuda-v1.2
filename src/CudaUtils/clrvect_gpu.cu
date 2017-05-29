extern "C" {
#include "../shape/head.h"
}
__global__ void clrvect_krnl(struct dat_t *ddat, int size, int s, int f) {
	/* Multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;

	if (offset < size) {
		if (ddat->set[s].type == DELAY)
			ddat->set[s].desc.deldop.frame[f].fit_s[offset] = 0.0;
		if (ddat->set[s].type == DOPPLER)
			ddat->set[s].desc.doppler.frame[f].fit_s[offset] = 0.0;
	}
}

