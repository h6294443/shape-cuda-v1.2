
extern "C" {
#include "../shape/head.h"
}
__device__ int dbg_ndop1, dbg_ndel1, dbg_xlim0, dbg_xlim1, dbg_ylim0, dbg_ylim1;

__global__ void dbg_print_fit_krnl1(struct dat_t *ddat, int s, int f){
	/* Single-threaded debug kernel */
	if (threadIdx.x == 0) {
		switch (ddat->set[s].type) {
		case DELAY:
			dbg_ndel1 = ddat->set[s].desc.deldop.frame[f].ndel;
			dbg_ndop1 = ddat->set[s].desc.deldop.frame[f].ndop;
			dbg_xlim0 = ddat->set[s].desc.deldop.frame[f].pos.xlim[0];
			dbg_xlim1 = ddat->set[s].desc.deldop.frame[f].pos.xlim[1];
			dbg_ylim0 = ddat->set[s].desc.deldop.frame[f].pos.ylim[0];
			dbg_ylim1 = ddat->set[s].desc.deldop.frame[f].pos.ylim[1];
			break;
		case DOPPLER:
			dbg_ndop1 = ddat->set[s].desc.doppler.frame[f].ndop;
			dbg_xlim0 = ddat->set[s].desc.doppler.frame[f].pos.xlim[0];
			dbg_xlim1 = ddat->set[s].desc.doppler.frame[f].pos.xlim[1];
			dbg_ylim0 = ddat->set[s].desc.doppler.frame[f].pos.ylim[0];
			dbg_ylim1 = ddat->set[s].desc.doppler.frame[f].pos.ylim[1];
			break;
		}
	}
}
__global__ void dbg_print_fit_krnl2(struct dat_t *ddat, float *fit, int s, int f) {
	/* ndop-threaded kernel */
	int idop = blockIdx.x * blockDim.x + threadIdx.x;// +1 ;

	if (idop < dbg_ndop1) {
		fit[idop] = ddat->set[s].desc.doppler.frame[f].fit_s[idop];
		printf("fit_s[%i]=%g\n", idop, ddat->set[s].desc.doppler.frame[f].fit_s[idop]);
	}
}
__global__ void dbg_print_lc_fit_krnl(struct dat_t *ddat, double *fit, int s, int n) {
	/* ndop-threaded kernel */
	int i = blockIdx.x * blockDim.x + threadIdx.x + 1;

	if (i <= n) {
		fit[i] = ddat->set[s].desc.lghtcrv.fit[i];
	}
}
__global__ void dbg_print_fit_deldop_krnl2(struct dat_t *ddat, float *fit, int s, int f){
	/* ndel*ndop-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;

	if (offset < (dbg_ndop1*dbg_ndel1))
		fit[offset] = ddat->set[s].desc.deldop.frame[f].fit_s[offset];
}

__host__ void dbg_print_fit(struct dat_t *ddat, int s, int f, const char
		*filename_fit, int gpuid) {
	/* Debug function that prints all Doppler frame fit values to csv */

	int idop, nThreads, ndop, xlim[2], ylim[2];
	FILE *fp_fit;
	float *fit, *host_fit;
	dim3 BLK,THD;
	//gpuid = GPU0;
	cudaSetDevice(gpuid);
	printf("\n %sfile created",filename_fit);


	/* Launch 1st debug kernel to get ndop and xlim/ylim	 */
	dbg_print_fit_krnl1<<<1,1>>>(ddat, s, f);
	checkErrorAfterKernelLaunch("dbg_print_fit_krnl1");
	deviceSyncAfterKernelLaunch("dbg_print_fit_krnl2");
	gpuErrchk(cudaMemcpyFromSymbol(&xlim[0], dbg_xlim0, sizeof(int),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&xlim[1], dbg_xlim1, sizeof(int),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&ylim[0], dbg_ylim0, sizeof(int),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&ylim[1], dbg_ylim1, sizeof(int),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&ndop, dbg_ndop1, sizeof(int),
			0, cudaMemcpyDeviceToHost));

	nThreads = (xlim[1] - xlim[0] + 1) * (ylim[1] - ylim[0] + 1);
	cudaCalloc((void**)&fit, sizeof(float), ndop);
//	host_fit = (float *)malloc(sizeof(float) * ndop);
	int maxThreads = 128;
	BLK.x = floor((maxThreads - 1 + ndop)/maxThreads);
	THD.x = maxThreads; // Thread block dimensions

	dbg_print_fit_krnl2<<<BLK,THD>>>(ddat, fit, s, f);
	checkErrorAfterKernelLaunch("dbg_print_fit_krnl_2");
	deviceSyncAfterKernelLaunch("dbg_print_fit_krnl_2");
//	gpuErrchk(cudaMemcpy(&host_fit, fit, sizeof(float)*ndop,
//			cudaMemcpyDeviceToHost));

//	fp_fit = fopen(filename_fit, "w+");
//	fprintf(fp_fit, "idop , ");
//	for (idop=0; idop<ndop; idop++)
//		fprintf(fp_fit,	"\n%i , %g", idop, fit[idop]);
//	fprintf(fp_fit, "\nxlim0 , %i", xlim[0]);
//	fprintf(fp_fit, "\nxlim1 , %i", xlim[1]);
//	fprintf(fp_fit, "\nylim0 , %i", ylim[0]);
//	fprintf(fp_fit, "\nylim1 , %i", ylim[1]);
//	fprintf(fp_fit, "\nthreads , %i", nThreads);
//	fclose(fp_fit);
	cudaFree(fit);
//	free(host_fit);
}

__host__ void dbg_print_fit_host(struct dat_t *ddat, int s, int f, char *filename_fit) {
	/* Debug function that prints all Doppler frame fit values to csv */

	int idop, nThreads;
	FILE *fp_fit;
	nThreads = (ddat->set[s].desc.doppler.frame[f].pos.xlim[1]-
			ddat->set[s].desc.doppler.frame[f].pos.xlim[0]+1)*
					(ddat->set[s].desc.doppler.frame[f].pos.ylim[1]-
							ddat->set[s].desc.doppler.frame[f].pos.ylim[0]+1);

	printf("\n %sfile created",filename_fit);
	fp_fit = fopen(filename_fit, "w+");

	fprintf(fp_fit, "idel , ");

	for (idop=1; idop<=ddat->set[s].desc.doppler.frame[f].ndop; idop++)
		fprintf(fp_fit,	"\n%i , %g", idop, ddat->set[s].desc.doppler.frame[f].fit[idop]);

	fprintf(fp_fit, "\nxlim0 , %i", ddat->set[s].desc.doppler.frame[f].pos.xlim[0]);
	fprintf(fp_fit, "\nxlim1 , %i", ddat->set[s].desc.doppler.frame[f].pos.xlim[1]);
	fprintf(fp_fit, "\nylim0 , %i", ddat->set[s].desc.doppler.frame[f].pos.ylim[0]);
	fprintf(fp_fit, "\nylim1 , %i", ddat->set[s].desc.doppler.frame[f].pos.ylim[1]);
	fprintf(fp_fit, "\nthreads , %i", nThreads);
	fclose(fp_fit);
}

__host__ void dbg_print_deldop_fit(struct dat_t *ddat, int s, int f, char *filename_fit) {
	/* Debug function that prints all Doppler frame fit values to csv */

	int idop, ndop, idel, ndel, nbins, nThreads, offset, xlim[2], ylim[2];
	FILE *fp_fit;
	float *fit_dd;
	dim3 BLK,THD;
	printf("\n %sfile created",filename_fit);
	float *host_fit;

	/* Check for dual-GPU mode and set device accordingly (GPU0 for even
	 * frames and GPU1 for odd frames) 	 */
	if ((f%2)==0)		/* Even frames */
		gpuErrchk(cudaSetDevice(GPU0));
	else if ((f%2)==1)	/* Odd frames */
		gpuErrchk(cudaSetDevice(GPU1));

	/* Launch 1st debug kernel to get ndop and xlim/ylim	 */
	dbg_print_fit_krnl1<<<1,1>>>(ddat, s, f);
	checkErrorAfterKernelLaunch("dbg_print_fit_krnl1");
	deviceSyncAfterKernelLaunch("dbg_print_fit_krnl2");
	gpuErrchk(cudaMemcpyFromSymbol(&xlim[0], dbg_xlim0, sizeof(int),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&xlim[1], dbg_xlim1, sizeof(int),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&ylim[0], dbg_ylim0, sizeof(int),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&ylim[1], dbg_ylim1, sizeof(int),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&ndop, dbg_ndop1, sizeof(int),
			0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&ndel, dbg_ndel1, sizeof(int),
				0, cudaMemcpyDeviceToHost));

	nThreads = (xlim[1] - xlim[0] + 1) * (ylim[1] - ylim[0] + 1);
	nbins = ndop * ndel;
	cudaCalloc((void**)&fit_dd, sizeof(float), nbins);
	host_fit = (float *)malloc(sizeof(float) * nbins);

	BLK.x = floor((maxThreadsPerBlock - 1 + nbins)/maxThreadsPerBlock);
	THD.x = maxThreadsPerBlock; // Thread block dimensions

	dbg_print_fit_deldop_krnl2<<<BLK,THD>>>(ddat, fit_dd, s, f);
	checkErrorAfterKernelLaunch("dbg_print_fit_deldop_krnl_2");
	cudaMemcpy(host_fit, fit_dd, sizeof(float)*nbins, cudaMemcpyDeviceToHost);

	cudaSetDevice(GPU0);
	cudaDeviceSynchronize();
	cudaSetDevice(GPU1);
	cudaDeviceSynchronize();

	if ((f%2)==0) 	 	/* Even frames */
		gpuErrchk(cudaSetDevice(GPU0));
	else if ((f%2)==1) 	/* Odd frames */
		gpuErrchk(cudaSetDevice(GPU1));


	fp_fit = fopen(filename_fit, "w+");

	/* Print top corner idop/idel label */
	fprintf(fp_fit, "idop/idel , ");

	/* Print top row idel values */
	for (idel=1; idel<=ndel; idel++)
		fprintf(fp_fit, "%i , ", idel);

	/* Print first entry in every row (except 1st): idop */
	for (idop=1; idop<=ndop; idop++) {
		fprintf(fp_fit,	"\n%i , ", idop);

		/* Write the rest of the row values: fit[idel][idop] */
		for (idel=1; idel<=ndel; idel++) {
			offset = (idop-1)*ndel + (idel-1);
			fprintf(fp_fit, " %g , ", host_fit[offset]);//fit_dd[offset]);
		}
	}
	fprintf(fp_fit, "\nxlim0 , %i", xlim[0]);
	fprintf(fp_fit, "\nxlim1 , %i", xlim[1]);
	fprintf(fp_fit, "\nylim0 , %i", ylim[0]);
	fprintf(fp_fit, "\nylim1 , %i", ylim[1]);
	fprintf(fp_fit, "\nthreads , %i", nThreads);
	fclose(fp_fit);

	/* Back to GPU0 before returning */
	gpuErrchk(cudaSetDevice(GPU0));
}

__host__ void dbg_print_deldop_fit_host(struct dat_t *ddat, int s, int f, char *filename_fit) {
	/* Debug function that prints all Delay-Doppler frame fit values to csv */

	int idop, ndop, idel, ndel, nThreads, xlim[2], ylim[2];
	FILE *fp_fit;
	printf("\n %sfile created",filename_fit);

	for (idop=0;idop<2;idop++){
		xlim[idop] = ddat->set[s].desc.deldop.frame[f].pos.xlim[idop];
		ylim[idop] = ddat->set[s].desc.deldop.frame[f].pos.ylim[idop];}

	ndel = ddat->set[s].desc.deldop.frame[f].ndel;
	ndop = ddat->set[s].desc.deldop.frame[f].ndop;
	nThreads = (xlim[1] - xlim[0] + 1) * (ylim[1] - ylim[0] + 1);
	fp_fit = fopen(filename_fit, "w+");

	/* Print top corner idop/idel label */
	fprintf(fp_fit, "idop/idel , ");

	/* Print top row idel values */
	for (idel=1; idel<=ndel; idel++)
		fprintf(fp_fit, "%i , ", idel);

	/* Print first entry in every row (except 1st): idop */
	for (idop=1; idop<=ndop; idop++) {
		fprintf(fp_fit,	"\n%i , ", idop);

		/* Write the rest of the row values: fit[idel][idop] */
		for (idel=1; idel<=ndel; idel++)
			fprintf(fp_fit, " %g , ", ddat->set[s].desc.deldop.frame[f].fit[idel][idop]);
	}
	fprintf(fp_fit, "\nxlim0 , %i", xlim[0]);
	fprintf(fp_fit, "\nxlim1 , %i", xlim[1]);
	fprintf(fp_fit, "\nylim0 , %i", ylim[0]);
	fprintf(fp_fit, "\nylim1 , %i", ylim[1]);
	fprintf(fp_fit, "\nthreads , %i", nThreads);
	fclose(fp_fit);
}

__host__ void dbg_print_lc_fit(struct dat_t *ddat, int s, char *filename_fit, int n) {
	/* Debug function that prints lightcurve fit values */

	int i;
	FILE *fp_fit;
	double *fit;
	dim3 BLK,THD;

	cudaCalloc((void**)&fit, sizeof(double), n);
	fit -= 1;
	THD.x = maxThreadsPerBlock;
	BLK.x = floor((THD.x - 1 + n)/THD.x);

	dbg_print_lc_fit_krnl<<<BLK,THD>>>(ddat, fit, s, n);
	checkErrorAfterKernelLaunch("dbg_print_lc_fit_krnl");
	deviceSyncAfterKernelLaunch("dbg_print_lc_fit_krnl");

	fp_fit = fopen(filename_fit, "w+");
	fprintf(fp_fit, "i , ");
	for (i=1; i<=n; i++)
		fprintf(fp_fit,	"\n%i , %g", i, fit[i]);
	fclose(fp_fit);
	//cudaFree(fit);
}

__host__ void dbg_print_lc_fit_host(struct lghtcrv_t *lghtcrv, char *filename_fit, int n) {
	/* Debug function that prints light curve fit values (host version) */

	int i;
	FILE *fp_fit;
	dim3 BLK,THD;
	fp_fit = fopen(filename_fit, "w+");
	fprintf(fp_fit, "i , ");
	for (i=1; i<=n; i++)
		fprintf(fp_fit,	"\n%i , %g", i, lghtcrv->fit[i]);
	fclose(fp_fit);
}
