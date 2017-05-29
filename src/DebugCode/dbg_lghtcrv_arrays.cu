
extern "C" {
#include "../shape/head.h"
}

__global__ void dbg_copy_lghtcrv_arrays_krnl(struct dat_t *ddat, int set, int n,
		double *fit, double *obs, double *oneovervar) {
	/* n-threaded kernel */
	int i = blockIdx.x * blockDim.x + threadIdx.x + 1;

	if ((i>0) && (i<=n)) {
		fit[i] = ddat->set[set].desc.lghtcrv.fit[i];
		obs[i] = ddat->set[set].desc.lghtcrv.obs[i];
		oneovervar[i] = ddat->set[set].desc.lghtcrv.oneovervar[i];
	}
}
__global__ void dbg_copy_lghtcrv_xyy2_krnl(struct dat_t *ddat, int set, int ncalc,
		double *x, double *y, double *y2) {
	/* n-threaded kernel */
	int i = blockIdx.x * blockDim.x + threadIdx.x + 1;

	if ((i>0) && (i<=ncalc)) {
		x[i] = ddat->set[set].desc.lghtcrv.x[i];
		y[i] = ddat->set[set].desc.lghtcrv.y[i];
		y2[i] = ddat->set[set].desc.lghtcrv.y2[i];
	}
}
__global__ void dbg_copy_lghtcrv_pos_arrays_krnl(struct dat_t *ddat, int set,
		int npixels, float *b, float *cosi, float *cose, int i) {
	/* npixels-threaded debug kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;

	if (offset < npixels) {
		b[offset] = ddat->set[set].desc.lghtcrv.rend[i].pos.b_s[offset];
		cosi[offset] = ddat->set[set].desc.lghtcrv.rend[i].pos.cosi_s[offset];
		cose[offset] = ddat->set[set].desc.lghtcrv.rend[i].pos.cose_s[offset];
	}
}
__global__ void dbg_copy_lghtcrv_pos_arrays2_krnl(struct pos_t **pos, int f,
		int npixels, float *b, float *cosi, float *cose) {
	/* npixels-threaded debug kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;

	if (offset < npixels) {
		b[offset] = pos[f]->b_s[offset];
		cosi[offset] = pos[f]->cosi_s[offset];
		cose[offset] = pos[f]->cose_s[offset];
	}
}
__global__ void dbg_copy_lghtcrv_pos_bd_krnl(struct pos_t **pos, int f,
		int npixels, double *bd) {
	/* npixels-threaded debug kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;

	if (offset < npixels) {
		bd[offset] = pos[f]->b_d[offset];
	}
}
__global__ void dbg_copy_lghtcrv_pos_arrays_full_krnl(struct pos_t **pos, int f,
		int npixels, float *b, float *cosi, float *cose, float *zz) {
	/* npixels-threaded debug kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;

	if (offset < npixels) {
		b[offset] = pos[f]->b_s[offset];
		cosi[offset] = pos[f]->cosi_s[offset];
		cose[offset] = pos[f]->cose_s[offset];
		zz[offset] = pos[f]->z_s[offset];
	}
}

__host__ void dbg_print_lghtcrv_arrays(struct dat_t *ddat, int set, int n,
		char *filename) {
	/* Debug function that will print lghtcrv->fit, lghtcrv->obs, and
	 * lghtcrv->oneovervar for any given lightcurve in a specified dataset and
	 * frame	 */
	/* Each array starts at 1 and ends on n */

	double *fit, *obs, *oneovervar;
	FILE *fp;
	dim3 BLK,THD;
	int i, thd = 64;

	cudaCalloc((void**)&fit, sizeof(double), n);
	cudaCalloc((void**)&obs, sizeof(double), n);
	cudaCalloc((void**)&oneovervar, sizeof(double), n);
	fit -= 1;
	obs -= 1;
	oneovervar -= 1;

	BLK.x = floor((thd-1+n)/thd);
	THD.x = thd;
	dbg_copy_lghtcrv_arrays_krnl<<<BLK,THD>>>(ddat, set, n, fit, obs,
			oneovervar);
	checkErrorAfterKernelLaunch("dbg_copy_lghtcrv_arrays_krnl");
	deviceSyncAfterKernelLaunch("dbg_copy_lghtcrv_arrays_krnl");

	printf("\n\nFilename: %s",filename);
	fp = fopen(filename, "w+");

	/* Print top corner idop/idel label */
	fprintf(fp, "set %i , GPU fit, GPU obs, GPU oneovervar\n", set);

	/* Print row */
	for (i=1; i<=n; i++)
		fprintf(fp, "%i , %g, %g, %g\n", i,	fit[i],	obs[i],	oneovervar[i]);

	fclose(fp);
}
__host__ void dbg_print_lghtcrv_arrays_host(struct lghtcrv_t *lghtcrv, int set, int n, char *filename) {
	/* Debug function that prints the three lightcurve arrays fit, obs,
	 * and oneovervar to a single csv file*/

	int i;
	FILE *fp;

	printf("\n\nFilename: %s",filename);
	fp = fopen(filename, "w+");

	/* Print top corner idop/idel label */
	fprintf(fp, "set %i , CPU fit, CPU obs, CPU oneovervar\n", set);

	/* Print row */
	for (i=1; i<=n; i++)
		fprintf(fp, "%i , %g, %g, %g\n", i,
				lghtcrv->fit[i],
				lghtcrv->obs[i],
				lghtcrv->oneovervar[i]);

	fclose(fp);
}
__host__ void dbg_print_lghtcrv_xyy2(struct dat_t *ddat, int set, int ncalc,
		char *filename) {
	/* This debug function prints to a csv file 'filename' these arrays:
	 * lghtcrv->x[lghtcrv->ncalc]
	 * lghtcrv->y[lghtcrv->ncalc]
	 * lghtcrv->y2[lghtcrv->ncalc]
	 * lghtcrv->fit[lghtcrv->n]
	 *
	 * where the lghtcrv is specified by 'set' in 'ddat'	 */
	double *x, *y, *y2;
	FILE *fp;
	dim3 BLK,THD;
	int i, thd = 64;

	cudaCalloc((void**)&x, sizeof(double), ncalc);
	cudaCalloc((void**)&y, sizeof(double), ncalc);
	cudaCalloc((void**)&y2, sizeof(double), ncalc);
	x -= 1;	y -= 1;	y2 -= 1;

	BLK.x = floor((thd-1+ncalc)/thd);
	THD.x = thd;
	dbg_copy_lghtcrv_xyy2_krnl<<<BLK,THD>>>(ddat, set, ncalc, x, y, y2);
	checkErrorAfterKernelLaunch("dbg_copy_lghtcrv_xyy2_krnl");
	deviceSyncAfterKernelLaunch("dbg_copy_lghtcrv_xyy2_krnl");
	printf("\n\nFilename: %s",filename);
	fp = fopen(filename, "w+");

	/* Print top corner idop/idel label */
	fprintf(fp, "set %i , GPU x, GPU y, GPU y2\n", set);

	/* Print row */
	for (i=1; i<=ncalc; i++)
		fprintf(fp, "%i , %g, %g, %g\n", i,	x[i],	y[i],	y2[i]);

	fclose(fp);
}
__host__ void dbg_print_lghtcrv_xyy2_host(struct lghtcrv_t *lghtcrv, int set, int ncalc, char *filename) {
	/* This debug function prints to a csv file 'filename' these arrays:
		 * lghtcrv->x[lghtcrv->ncalc]
		 * lghtcrv->y[lghtcrv->ncalc]
		 * lghtcrv->y2[lghtcrv->ncalc]
		 *
		 * where the lghtcrv is specified by 'set' in 'ddat'	 */

	int i;
	FILE *fp;

	printf("\n\nFilename: %s",filename);
	fp = fopen(filename, "w+");

	/* Print top corner idop/idel label */
	fprintf(fp, "set %i , CPU x, CPU y, CPU y2\n", set);

	/* Print row */
	for (i=1; i<=ncalc; i++)
		fprintf(fp, "%i , %g, %g, %g\n", i,
				lghtcrv->x[i],
				lghtcrv->y[i],
				lghtcrv->y2[i]);

	fclose(fp);
}
__host__ void dbg_print_lghtcrv_pos_arrays(struct dat_t *ddat, int set, int f, int npixels, int n) {
	/* This debug function prints the GPU arrays:
	 *  - pos->cosi_s
	 *  - pos->cose_s
	 *  - pos->b_s,
	 *
	 *  all of length nPixels in the lghtcrv specified by 'set' in 'ddat' */
	float *b, *cosi, *cose;
	FILE *fp_b, *fp_cosi, *fp_cose;
	dim3 BLK,THD;
	char *filename;
	int i, j, pxa, thd = 256;

	cudaCalloc((void**)&b, sizeof(float), npixels);
	cudaCalloc((void**)&cosi, sizeof(float), npixels);
	cudaCalloc((void**)&cose, sizeof(float), npixels);

	BLK.x = floor((thd-1+npixels)/thd);	THD.x = thd;
	dbg_copy_lghtcrv_pos_arrays_krnl<<<BLK,THD>>>(ddat, set, npixels, b, cosi, cose, f);
	checkErrorAfterKernelLaunch("dbg_copy_lghtcrv_pos_arrays_krnl");
	deviceSyncAfterKernelLaunch("dbg_copy_lghtcrv_pos_arrays_krnl");

	filename = "dbg_lghtcrv_pos_b.csv";
	fp_b = fopen(filename, "w+");
	filename = "dbg_lghtcrv_pos_cosi_s.csv";
	fp_cosi = fopen(filename, "w+");
	filename = "dbg_lghtcrv_pos_cose_s.csv";
	fp_cose = fopen(filename, "w+");

	/* Print top corner set label */
	fprintf(fp_cosi,"s%i i%i , ", set,f);
	fprintf(fp_cose,"s%i i%i , ", set,f);
	fprintf(fp_b, 	"s%i i%i , ", set,f);

	/* Print i values along top of table */
	for (i=-n; i<=n; i++) {
		fprintf(fp_cosi, "%i, ", i);
		fprintf(fp_cose, "%i, ", i);
		fprintf(fp_b, "%i, ", i);
	}
	fprintf(fp_cosi, "\n");
	fprintf(fp_cose, "\n");
	fprintf(fp_b, "\n");

	for (j=-n; j<=n; j++) {
		fprintf(fp_b, "%i, ", j);	/* j-entry on far left */
		fprintf(fp_cosi, "%i, ", j);
		fprintf(fp_cose, "%i, ", j);

		for (i=-n; i<=n; i++) {
			pxa = (j+n)*(2*n+1) + (i+n);
			fprintf(fp_b, "%g, ", b[pxa]);
			fprintf(fp_cosi, "%g, ", cosi[pxa]);
			fprintf(fp_cose, "%g, ", cose[pxa]);
		}
		fprintf(fp_b, "\n");
		fprintf(fp_cosi, "\n");
		fprintf(fp_cose, "\n");
	}

	fclose(fp_b);
	fclose(fp_cosi);
	fclose(fp_cose);
}
__host__ void dbg_print_pos_bd(struct pos_t **pos, int f, int npixels, int n) {
	/* This debug function prints the GPU arrays:
	 *  - pos->b_d, currently used only experimentally in light curves
	 *
	 *  all of length nPixels in the lghtcrv specified by 'set' in 'ddat' */
	double *bd;
	FILE *fp_b;
	dim3 BLK,THD;
	char *filename;
	int i, j, pxa;
	THD.x = maxThreadsPerBlock;

	cudaCalloc((void**)&bd, sizeof(double), npixels);

	BLK.x = floor((THD.x - 1 + npixels ) / THD.x);
	dbg_copy_lghtcrv_pos_bd_krnl<<<BLK,THD>>>(pos, f, npixels, bd);
	checkErrorAfterKernelLaunch("dbg_copy_lghtcrv_pos_bd_krnl");
	deviceSyncAfterKernelLaunch("dbg_copy_lghtcrv_pos_bd_krnl");

	filename = "dbg_pos_bd.csv";
	fp_b = fopen(filename, "w+");

	fprintf(fp_b, 	"i%i , ", f);

	/* Print i values along top of table */
	for (i=-n; i<=n; i++)
		fprintf(fp_b, "%i, ", i);
	fprintf(fp_b, "\n");

	for (j=-n; j<=n; j++) {
		fprintf(fp_b, "%i, ", j);	/* j-entry on far left */

		for (i=-n; i<=n; i++) {
			pxa = (j+n)*(2*n+1) + (i+n);
			fprintf(fp_b, "%g, ", bd[pxa]);
		}
		fprintf(fp_b, "\n");
	}
	fclose(fp_b);
	cudaFree(bd);
}
__host__ void dbg_print_lghtcrv_pos_arrays_host(struct lghtcrv_t *lghtcrv,
		int f, int set) {
	/* This debug function prints the CPU arrays:
	 *  - pos->cosi
	 *  - pos->cose
	 *  - pos->b,
	 *
	 *  all of length nPixels in the lghtcrv specified by 'set' in 'ddat' */
	int i, j, n;
	FILE *fp_cosi, *fp_cose, *fp_b;
	struct pos_t *pos;
	char *fn;
	pos = &lghtcrv->rend[f].pos;
	n = pos->n;

	fn = "dbg_lghtcrv_pos_cosi_CPU.csv";
	fp_cosi = fopen(fn, "w+");
	fn = "dbg_lghtcrv_pos_cose_CPU.csv";
	fp_cose = fopen(fn, "w+");
	fn = "dbg_lghtcrv_pos_b_CPU.csv";
	fp_b = fopen(fn, "w+");

	/* Print top corner set label */
	fprintf(fp_cosi, "set %i , ", set);
	fprintf(fp_cose, "set %i , ", set);
	fprintf(fp_b, "set %i , ", set);

	/* Print i values along top of table */
	for (i=-n; i<=n; i++) {
		fprintf(fp_cosi, "%i, ", i);
		fprintf(fp_cose, "%i, ", i);
		fprintf(fp_b, "%i, ", i);
	}
	fprintf(fp_cosi, "\n");
	fprintf(fp_cose, "\n");
	fprintf(fp_b, "\n");

	for (j=-n; j<=n; j++) {
		fprintf(fp_b, "%i, ", j);	/* j-entry on far left */
		fprintf(fp_cosi, "%i, ", j);
		fprintf(fp_cose, "%i, ", j);

		for (i=-n; i<=n; i++) {
			fprintf(fp_b, "%g, ", pos->b[i][j]);
			fprintf(fp_cosi, "%g, ", pos->cosi[i][j]);
			fprintf(fp_cose, "%g, ", pos->cose[i][j]);
		}
		fprintf(fp_b, "\n");
		fprintf(fp_cosi, "\n");
		fprintf(fp_cose, "\n");
	}

	fclose(fp_b);
	fclose(fp_cosi);
	fclose(fp_cose);
}
__host__ void dbg_print_pos_arrays2(struct pos_t **pos, int f, int npixels, int n) {
	/* This debug function prints the GPU arrays:
	 *  - pos->cosi_s
	 *  - pos->cose_s
	 *  - pos->b_s,
	 *
	 *  all of length nPixels in the lghtcrv specified by 'set' in 'ddat' */
	float *b, *cosi, *cose;
	FILE *fp_b, *fp_cosi, *fp_cose;
	dim3 BLK,THD;
	char *filename;
	int i, j, pxa, thd = 256;

	cudaCalloc((void**)&b, sizeof(float), npixels);
	cudaCalloc((void**)&cosi, sizeof(float), npixels);
	cudaCalloc((void**)&cose, sizeof(float), npixels);

	BLK.x = floor((thd-1+npixels)/thd);	THD.x = thd;
	dbg_copy_lghtcrv_pos_arrays2_krnl<<<BLK,THD>>>(pos, f, npixels, b, cosi, cose);
	checkErrorAfterKernelLaunch("dbg_copy_lghtcrv_pos_arrays_krnl");
	deviceSyncAfterKernelLaunch("dbg_copy_lghtcrv_pos_arrays_krnl");

	filename = "dbg_pos_b.csv";
	fp_b = fopen(filename, "w+");
	filename = "dbg_pos_cosi_s.csv";
	fp_cosi = fopen(filename, "w+");
	filename = "dbg_pos_cose_s.csv";
	fp_cose = fopen(filename, "w+");

	/* Print top corner set label */
	fprintf(fp_cosi,"i%i , ", f);
	fprintf(fp_cose,"i%i , ", f);
	fprintf(fp_b, 	"i%i , ", f);

	/* Print i values along top of table */
	for (i=-n; i<=n; i++) {
		fprintf(fp_cosi, "%i, ", i);
		fprintf(fp_cose, "%i, ", i);
		fprintf(fp_b, "%i, ", i);
	}
	fprintf(fp_cosi, "\n");
	fprintf(fp_cose, "\n");
	fprintf(fp_b, "\n");

	for (j=-n; j<=n; j++) {
		fprintf(fp_b, "%i, ", j);	/* j-entry on far left */
		fprintf(fp_cosi, "%i, ", j);
		fprintf(fp_cose, "%i, ", j);

		for (i=-n; i<=n; i++) {
			pxa = (j+n)*(2*n+1) + (i+n);
			fprintf(fp_b, "%g, ", b[pxa]);
			fprintf(fp_cosi, "%g, ", cosi[pxa]);
			fprintf(fp_cose, "%g, ", cose[pxa]);
		}
		fprintf(fp_b, "\n");
		fprintf(fp_cosi, "\n");
		fprintf(fp_cose, "\n");
	}

	fclose(fp_b);
	fclose(fp_cosi);
	fclose(fp_cose);
}
