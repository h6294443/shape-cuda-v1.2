extern "C" {
#include "../shape/head.h"
}

__global__ void clrvect_krnl(struct dat_t *ddat, int size, int s, int f, int dblflg) {
	/* Multi-threaded kernel */
	int idel, idop, offset = blockIdx.x * blockDim.x + threadIdx.x;

	if (offset < size) {
		if (ddat->set[s].type == DELAY) {
			if (dblflg) {
				idel = offset % ddat->set[s].desc.deldop.frame[f].ndel + 1;
				idop = offset / ddat->set[s].desc.deldop.frame[f].ndel + 1;
				ddat->set[s].desc.deldop.frame[f].fit[idel][idop] = 0.0;
			}
			else
				ddat->set[s].desc.deldop.frame[f].fit_s[offset] = 0.0;

		}
		if (ddat->set[s].type == DOPPLER) {
			if (dblflg)
				ddat->set[s].desc.doppler.frame[f].fit[offset+1] = 0.0;
			else
				ddat->set[s].desc.doppler.frame[f].fit_s[offset+1] = 0.0;
		}
	}
}
