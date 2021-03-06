
extern "C" {
#include "../shape/head.h"
}

__global__ void dbg_copy_pos_arrays_full_krnl32(struct pos_t **pos, int f,
		int npixels, double *b, float *cosi, float *cose, float *zz) {
	/* npixels-threaded debug kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int n = pos[f]->n;
	int i = offset % (2*n+1) - n;
	int j = offset / (2*n+1) - n;

	if (offset < npixels) {
		b[offset] = pos[f]->b[i][j];
		cosi[offset] = pos[f]->cosi_s[offset];
		cose[offset] = pos[f]->cose_s[offset];
		zz[offset] = pos[f]->z_s[offset];
	}
}
__global__ void dbg_copy_pos_arrays_full_krnl64(struct pos_t **pos, int f,
		int npixels, double *b, double *cosi, double *cose, double *zz, int *fac) {
	/* npixels-threaded debug kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int n = pos[f]->n;
	int i = offset % (2*n+1) - n;
	int j = offset / (2*n+1) - n;
	if (offset < npixels) {
//		b[offset] = pos[f]->b[i][j];
		cosi[offset] = pos[f]->cosi[i][j];
		cose[offset] = pos[f]->cose[i][j];
		zz[offset] = pos[f]->z[i][j];
//		fac[offset] = pos[f]->f[i][j];
	}
}
__global__ void dbg_copy_lc_pos_arrays_full_krnl64(struct pos_t **pos, int f,
		int npixels, double *b, double *cosi, double *cose, double *zz,
		double *zill, double *cosill) {
	/* npixels-threaded debug kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int n = pos[f]->n;
	int i = offset % (2*n+1) - n;
	int j = offset / (2*n+1) - n;
	if (offset < npixels) {
		b[offset] = pos[f]->b[i][j];
		cosi[offset] = pos[f]->cosi[i][j];
		cose[offset] = pos[f]->cose[i][j];
		zz[offset] = pos[f]->z[i][j];
		zill[offset] = pos[f]->zill[i][j];
		cosill[offset] = pos[f]->cosill[i][j];
	}
}

__host__ void dbg_print_pos_arrays2_host(struct pos_t *pos) {
	/* This debug function prints the CPU arrays:
	 *  - pos->cosi
	 *  - pos->cose
	 *  - pos->b,
	 *
	 *  all of length nPixels in the lghtcrv specified by 'set' in 'ddat' */
	int i, j, n;
	FILE *fp_cosi, *fp_cose, *fp_b;
	const char *fn;
	n = pos->n;

	fn = "dbg_pos_cosi_CPU.csv";
	fp_cosi = fopen(fn, "w+");
	fn = "dbg_pos_cose_CPU.csv";
	fp_cose = fopen(fn, "w+");
	fn = "dbg_pos_b_CPU.csv";
	fp_b = fopen(fn, "w+");

	/* Print top corner set label */
	fprintf(fp_cosi, "set , ");
	fprintf(fp_cose, "set , ");
	fprintf(fp_b, "set , ");

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

__host__ void dbg_print_pos_z_host(struct pos_t *pos, const char *fn) {
	/* This debug function prints the CPU arrays:
	 *	 - pos->z
	 *
	 *  all of length nPixels in the lghtcrv specified by 'set' in 'ddat' */
	int i, j, n;
	FILE *fp_z;
	n = pos->n;
	fp_z = fopen(fn, "w+");

	/* Print top corner set label */
	fprintf(fp_z, "s?f?, ");

	/* Print i values along top of table */
	for (i=-n; i<=n; i++)
		fprintf(fp_z, "%i, ", i);
	fprintf(fp_z, "\n");

	for (j=-n; j<=n; j++) {
		fprintf(fp_z, "%i, ", j);	/* j-entry on far left */
		for (i=-n; i<=n; i++)
			fprintf(fp_z, "%g, ", pos->z[i][j]);
		fprintf(fp_z, "\n");
	}
	fclose(fp_z);
}

__host__ void dbg_print_pos_arrays_full32(struct pos_t **pos, int f,
		int npixels, int n) {
	/* This debug function prints the GPU arrays:
	 *  - pos->cosi_s
	 *  - pos->cose_s
	 *  - pos->b_s,
	 *
	 *  all of length nPixels in the lghtcrv specified by 'set' in 'ddat' */
	float *cosi, *cose, *zz;
	FILE *fp_b, *fp_cosi, *fp_cose, *fp_z;
	dim3 BLK,THD;
	const char *filename;
	int i, j, pxa, thd = 256;
	double *b;

	cudaCalloc1((void**)&b, sizeof(double), npixels);
	cudaCalloc1((void**)&cosi, sizeof(float), npixels);
	cudaCalloc1((void**)&cose, sizeof(float), npixels);
	cudaCalloc1((void**)&zz, sizeof(float), npixels);

	BLK.x = floor((thd-1+npixels)/thd);	THD.x = thd;
	dbg_copy_pos_arrays_full_krnl32<<<BLK,THD>>>(pos, f, npixels, b, cosi, cose, zz);
	checkErrorAfterKernelLaunch("dbg_copy_pos_arrays_full_krnl");
	deviceSyncAfterKernelLaunch("dbg_copy_pos_arrays_full_krnl");

	filename = "1080Ti_b_32.csv";
	fp_b = fopen(filename, "w+");
	filename = "1080Ti_cosis32.csv";
	fp_cosi = fopen(filename, "w+");
	filename = "1080Ti_coses32.csv";
	fp_cose = fopen(filename, "w+");
	filename = "1080Ti_zs32.csv";
	fp_z = fopen(filename, "w+");

	/* Print top corner set label */
	fprintf(fp_cosi,"i%i , ", f);
	fprintf(fp_cose,"i%i , ", f);
	fprintf(fp_b, 	"i%i , ", f);
	fprintf(fp_z,   "i%i , ", f);

	/* Print i values along top of table */
	for (i=-n; i<=n; i++) {
		fprintf(fp_cosi, "%i, ", i);
		fprintf(fp_cose, "%i, ", i);
		fprintf(fp_b, "%i, ", i);
		fprintf(fp_z, "%i, ", i);
	}
	fprintf(fp_cosi, "\n");
	fprintf(fp_cose, "\n");
	fprintf(fp_b, "\n");
	fprintf(fp_z, "\n");

	for (j=-n; j<=n; j++) {
		fprintf(fp_b, "%i, ", j);	/* j-entry on far left */
		fprintf(fp_cosi, "%i, ", j);
		fprintf(fp_cose, "%i, ", j);
		fprintf(fp_z, "%i, ", j);

		for (i=-n; i<=n; i++) {
			pxa = (j+n)*(2*n+1) + (i+n);
			fprintf(fp_b, "%g, ", b[pxa]);
			fprintf(fp_cosi, "%g, ", cosi[pxa]);
			fprintf(fp_cose, "%g, ", cose[pxa]);
			fprintf(fp_z, "%g, ", zz[pxa]);
		}
		fprintf(fp_b, "\n");
		fprintf(fp_cosi, "\n");
		fprintf(fp_cose, "\n");
		fprintf(fp_z, "\n");
	}

	fclose(fp_b);
	fclose(fp_cosi);
	fclose(fp_cose);
	fclose(fp_z);
	cudaFree(zz);
	cudaFree(cosi);
	cudaFree(cose);
	cudaFree(b);
}

__host__ void dbg_print_pos_arrays_full64(struct pos_t **pos, int f,
		int npixels, int n) {
	/* This debug function prints the GPU arrays:
	 *  - pos->cosi_s
	 *  - pos->cose_s
	 *  - pos->b_s,
	 *
	 *  all of length nPixels in the lghtcrv specified by 'set' in 'ddat' */
	double *b, *cosi, *cose, *zz;
	FILE *fp_b, *fp_cosi, *fp_cose, *fp_z, *fp_f;
	dim3 BLK,THD;
	const char *filename;
	int i, j, pxa, thd = 256, *fac;

//	cudaCalloc1((void**)&b, sizeof(double), npixels);
	cudaCalloc1((void**)&cosi, sizeof(double), npixels);
	cudaCalloc1((void**)&cose, sizeof(double), npixels);
	cudaCalloc1((void**)&zz, sizeof(double), npixels);
//	cudaCalloc1((void**)&fac, sizeof(int), npixels);

	BLK.x = floor((thd-1+npixels)/thd);	THD.x = thd;
	dbg_copy_pos_arrays_full_krnl64<<<BLK,THD>>>(pos, f, npixels, b, cosi, cose, zz, fac);
	checkErrorAfterKernelLaunch("dbg_copy_pos_arrays_full_krnl");
	deviceSyncAfterKernelLaunch("dbg_copy_pos_arrays_full_krnl");

//	filename = "1080Ti_b_64.csv";
//	fp_b = fopen(filename, "w+");
	filename = "1080Ti_cosi64.csv";
	fp_cosi = fopen(filename, "w+");
	filename = "1080Ti_cose64.csv";
	fp_cose = fopen(filename, "w+");
	filename = "1080Ti_z64.csv";
	fp_z = fopen(filename, "w+");
//	filename = "1080Ti_fac64.csv";
//	fp_f = fopen(filename, "w+");

	/* Print top corner set label */
	fprintf(fp_cosi,"i%i , ", f);
	fprintf(fp_cose,"i%i , ", f);
//	fprintf(fp_b, 	"i%i , ", f);
	fprintf(fp_z,   "i%i , ", f);
//	fprintf(fp_f,   "i%i , ", f);

	/* Print i values along top of table */
	for (i=-n; i<=n; i++) {
		fprintf(fp_cosi, "%i, ", i);
		fprintf(fp_cose, "%i, ", i);
//		fprintf(fp_b, "%i, ", i);
		fprintf(fp_z, "%i, ", i);
//		fprintf(fp_f, "%i, ", i);
	}
	fprintf(fp_cosi, "\n");
	fprintf(fp_cose, "\n");
//	fprintf(fp_b, "\n");
	fprintf(fp_z, "\n");
//	fprintf(fp_f, "\n");

	for (j=n; j>=(-n); j--) {
//		fprintf(fp_b, "%i, ", j);	/* j-entry on far left */
		fprintf(fp_cosi, "%i, ", j);
		fprintf(fp_cose, "%i, ", j);
		fprintf(fp_z, "%i, ", j);
//		fprintf(fp_f, "%i, ", j);

		for (i=-n; i<=n; i++) {
			pxa = (j+n)*(2*n+1) + (i+n);
//			fprintf(fp_b, "%g, ", b[pxa]);
			fprintf(fp_cosi, "%g, ", cosi[pxa]);
			fprintf(fp_cose, "%g, ", cose[pxa]);
			fprintf(fp_z, "%g, ", zz[pxa]);
//			fprintf(fp_f, "%i, ", fac[pxa]);
		}
//		fprintf(fp_b, "\n");
		fprintf(fp_cosi, "\n");
		fprintf(fp_cose, "\n");
		fprintf(fp_z, "\n");
//		fprintf(fp_f, "\n");
	}

//	fclose(fp_b);
	fclose(fp_cosi);
	fclose(fp_cose);
	fclose(fp_z);
//	fclose(fp_f);
	cudaFree(zz);
	cudaFree(cosi);
	cudaFree(cose);
	cudaFree(b);
	cudaFree(fac);
}

__host__ void dbg_print_lc_pos_arrays_full64(struct pos_t **pos, int f,
		int npixels, int n) {
	/* This debug function prints the GPU arrays:
	 *  - pos->cosi_s
	 *  - pos->cose_s
	 *  - pos->b_s,
	 *
	 *  all of length nPixels in the lghtcrv specified by 'set' in 'ddat' */
	double *b, *cosi, *cose, *zz, *zill, *cosill;
	FILE *fp_b, *fp_cosi, *fp_cose, *fp_z, *fp_zill, *fp_cosill;
	dim3 BLK,THD;
	const char *filename;
	int i, j, pxa, thd = 256, *fac;

	cudaCalloc1((void**)&b, sizeof(double), npixels);
	cudaCalloc1((void**)&cosi, sizeof(double), npixels);
	cudaCalloc1((void**)&cose, sizeof(double), npixels);
	cudaCalloc1((void**)&zz, sizeof(double), npixels);
	cudaCalloc1((void**)&cosill, sizeof(double), npixels);
	cudaCalloc1((void**)&zill, sizeof(double), npixels);
//	cudaCalloc1((void**)&fac, sizeof(int), npixels);

	BLK.x = floor((thd-1+npixels)/thd);	THD.x = thd;
	dbg_copy_lc_pos_arrays_full_krnl64<<<BLK,THD>>>(pos, f, npixels, b, cosi,
			cose, zz, zill, cosill);
	checkErrorAfterKernelLaunch("dbg_copy_lc_pos_arrays_full_krnl");
	deviceSyncAfterKernelLaunch("dbg_copy_lc_pos_arrays_full_krnl");

	filename = "1080Ti_b_64.csv";
	fp_b = fopen(filename, "w+");
	filename = "1080Ti_cosi64.csv";
	fp_cosi = fopen(filename, "w+");
	filename = "1080Ti_cose64.csv";
	fp_cose = fopen(filename, "w+");
	filename = "1080Ti_z64.csv";
	fp_z = fopen(filename, "w+");
	filename = "1080Ti_zill64.csv";
	fp_zill = fopen(filename, "w+");
	filename = "1080Ti_cosill64.csv";
	fp_cosill = fopen(filename, "w+");

	/* Print top corner set label */
	fprintf(fp_cosi,"i%i , ", f);
	fprintf(fp_cose,"i%i , ", f);
	fprintf(fp_b, 	"i%i , ", f);
	fprintf(fp_z,   "i%i , ", f);
	fprintf(fp_zill,"i%i , ", f);
	fprintf(fp_cosill,"i%i , ", f);

	/* Print i values along top of table */
	for (i=-n; i<=n; i++) {
		fprintf(fp_cosi, "%i, ", i);
		fprintf(fp_cose, "%i, ", i);
		fprintf(fp_b, "%i, ", i);
		fprintf(fp_z, "%i, ", i);
		fprintf(fp_zill, "%i, ", i);
		fprintf(fp_cosill, "%i, ", i);
	}
	fprintf(fp_cosi, "\n");
	fprintf(fp_cose, "\n");
	fprintf(fp_b, "\n");
	fprintf(fp_z, "\n");
	fprintf(fp_zill, "\n");
	fprintf(fp_cosill, "\n");

	for (j=n; j>=(-n); j--) {
		fprintf(fp_b, "%i, ", j);	/* j-entry on far left */
		fprintf(fp_cosi, "%i, ", j);
		fprintf(fp_cose, "%i, ", j);
		fprintf(fp_z, "%i, ", j);
		fprintf(fp_zill, "%i, ", j);
		fprintf(fp_cosill, "%i, ", j);

		for (i=-n; i<=n; i++) {
			pxa = (j+n)*(2*n+1) + (i+n);
			fprintf(fp_b, "%g, ", b[pxa]);
			fprintf(fp_cosi, "%g, ", cosi[pxa]);
			fprintf(fp_cose, "%g, ", cose[pxa]);
			fprintf(fp_z, "%g, ", zz[pxa]);
			fprintf(fp_zill, "%g, ", zill[pxa]);
			fprintf(fp_cosill, "%g, ", cosill[pxa]);
		}
		fprintf(fp_b, "\n");
		fprintf(fp_cosi, "\n");
		fprintf(fp_cose, "\n");
		fprintf(fp_z, "\n");
		fprintf(fp_zill, "\n");
		fprintf(fp_cosill, "\n");
	}

	fclose(fp_b);
	fclose(fp_cosi);
	fclose(fp_cose);
	fclose(fp_z);
	fclose(fp_zill);
	fclose(fp_cosill);
	cudaFree(zz);
	cudaFree(cosi);
	cudaFree(cose);
	cudaFree(b);
	cudaFree(zill);
	cudaFree(cosill);
}

__host__ void dbg_print_pos_arrays_full_host(struct pos_t *pos) {
	/* This debug function prints the CPU arrays:
	 *  - pos->cosi
	 *  - pos->cose
	 *  - pos->b,
	 *
	 *  all of length nPixels in the lghtcrv specified by 'set' in 'ddat' */
	int i, j, n;
	FILE *fp_cosi, *fp_cose, *fp_b, *fp_z, *fp_f;
	const char *fn;
	n = pos->n;

	fn = "CPU_pos-cosi.csv";
	fp_cosi = fopen(fn, "w+");
	fn = "CPU_pos-cose.csv";
	fp_cose = fopen(fn, "w+");
	fn = "CPU_pos-b.csv";
	fp_b = fopen(fn, "w+");
	fn = "CPU_pos-z.csv";
	fp_z = fopen(fn, "w+");
	fn = "CPU_pos-fac.csv";
	fp_f = fopen(fn, "w+");

	/* Print top corner set label */
	fprintf(fp_cosi, "set , ");
	fprintf(fp_cose, "set , ");
	fprintf(fp_b, "set , ");
	fprintf(fp_z, "set , ");
	fprintf(fp_f, "set , ");

	/* Print i values along top of table */
	for (i=-n; i<=n; i++) {
		fprintf(fp_cosi, "%i, ", i);
		fprintf(fp_cose, "%i, ", i);
		fprintf(fp_b, "%i, ", i);
		fprintf(fp_z, "%i, ", i);
		fprintf(fp_f, "%i, ", i);
	}
	fprintf(fp_cosi, "\n");
	fprintf(fp_cose, "\n");
	fprintf(fp_b, "\n");
	fprintf(fp_z, "\n");
	fprintf(fp_f, "\n");

	for (j=n; j>=(-n); j--) {
		fprintf(fp_b, "%i, ", j);	/* j-entry on far left */
		fprintf(fp_cosi, "%i, ", j);
		fprintf(fp_cose, "%i, ", j);
		fprintf(fp_z, "%i, ", j);
		fprintf(fp_f, "%i, ", j);

		for (i=-n; i<=n; i++) {
			fprintf(fp_b, "%g, ", pos->b[i][j]);
			fprintf(fp_z, "%g, ", pos->z[i][j]);
			fprintf(fp_cosi, "%g, ", pos->cosi[i][j]);
			fprintf(fp_cose, "%g, ", pos->cose[i][j]);
			fprintf(fp_f, "%i, ", pos->f[i][j]);
		}
		fprintf(fp_b, "\n");
		fprintf(fp_z, "\n");
		fprintf(fp_cosi, "\n");
		fprintf(fp_cose, "\n");
		fprintf(fp_f, "\n");
	}

	fclose(fp_b);
	fclose(fp_z);
	fclose(fp_cosi);
	fclose(fp_cose);
	fclose(fp_f);
}

__host__ void dbg_print_lc_pos_arrays_full_host(struct pos_t *pos) {
	/* This debug function prints the CPU arrays:
	 *  - pos->cosi
	 *  - pos->cose
	 *  - pos->b,
	 *
	 *  all of length nPixels in the lghtcrv specified by 'set' in 'ddat' */
	int i, j, n;
	FILE *fp_cosi, *fp_cose, *fp_cosill, *fp_b, *fp_z, *fp_zill, *fp_f;
	const char *fn;
	n = pos->n;

	fn = "CPU_pos-cosi.csv";
	fp_cosi = fopen(fn, "w+");
	fn = "CPU_pos-cosill.csv";
	fp_cosill = fopen(fn, "w+");
	fn = "CPU_pos-cose.csv";
	fp_cose = fopen(fn, "w+");
	fn = "CPU_pos-b.csv";
	fp_b = fopen(fn, "w+");
	fn = "CPU_pos-z.csv";
	fp_z = fopen(fn, "w+");
	fn = "CPU_pos-zill.csv";
	fp_zill = fopen(fn, "w+");
	fn = "CPU_pos-fac.csv";
	fp_f = fopen(fn, "w+");

	/* Print top corner set label */
	fprintf(fp_cosi, "set , ");
	fprintf(fp_cosill, "set , ");
	fprintf(fp_cose, "set , ");
	fprintf(fp_b, "set , ");
	fprintf(fp_z, "set , ");
	fprintf(fp_zill, "set , ");
	fprintf(fp_f, "set , ");

	/* Print i values along top of table */
	for (i=-n; i<=n; i++) {
		fprintf(fp_cosi, "%i, ", i);
		fprintf(fp_cosill, "%i, ", i);
		fprintf(fp_cose, "%i, ", i);
		fprintf(fp_b, "%i, ", i);
		fprintf(fp_z, "%i, ", i);
		fprintf(fp_zill, "%i, ", i);
		fprintf(fp_f, "%i, ", i);
	}
	fprintf(fp_cosi, "\n");
	fprintf(fp_cosill, "\n");
	fprintf(fp_cose, "\n");
	fprintf(fp_b, "\n");
	fprintf(fp_z, "\n");
	fprintf(fp_zill, "\n");
	fprintf(fp_f, "\n");

	for (j=n; j>=(-n); j--) {
		fprintf(fp_b, "%i, ", j);	/* j-entry on far left */
		fprintf(fp_cosi, "%i, ", j);
		fprintf(fp_cosill, "%i, ", j);
		fprintf(fp_cose, "%i, ", j);
		fprintf(fp_z, "%i, ", j);
		fprintf(fp_zill, "%i, ", j);
		fprintf(fp_f, "%i, ", j);

		for (i=-n; i<=n; i++) {
			fprintf(fp_b, "%g, ", pos->b[i][j]);
			fprintf(fp_z, "%g, ", pos->z[i][j]);
			fprintf(fp_zill, "%g, ", pos->zill[i][j]);
			fprintf(fp_cosi, "%g, ", pos->cosi[i][j]);
			fprintf(fp_cosill, "%g, ", pos->cosill[i][j]);
			fprintf(fp_cose, "%g, ", pos->cose[i][j]);
			fprintf(fp_f, "%i, ", pos->f[i][j]);
		}
		fprintf(fp_b, "\n");
		fprintf(fp_z, "\n");
		fprintf(fp_zill, "\n");
		fprintf(fp_cosi, "\n");
		fprintf(fp_cosill, "\n");
		fprintf(fp_cose, "\n");
		fprintf(fp_f, "\n");
	}

	fclose(fp_b);
	fclose(fp_z);
	fclose(fp_zill);
	fclose(fp_cosi);
	fclose(fp_cosill);
	fclose(fp_cose);
	fclose(fp_f);
}
