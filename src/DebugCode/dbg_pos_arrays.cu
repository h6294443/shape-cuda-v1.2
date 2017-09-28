
extern "C" {
#include "../shape/head.h"
}

__global__ void dbg_copy_pos_arrays_full_krnl(struct pos_t **pos, int f,
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
__host__ void dbg_print_pos_arrays_full(struct pos_t **pos, int f,
		int npixels, int n) {
	/* This debug function prints the GPU arrays:
	 *  - pos->cosi_s
	 *  - pos->cose_s
	 *  - pos->b_s,
	 *
	 *  all of length nPixels in the lghtcrv specified by 'set' in 'ddat' */
	float *b, *cosi, *cose, *zz;
	FILE *fp_b, *fp_cosi, *fp_cose, *fp_z;
	dim3 BLK,THD;
	const char *filename;
	int i, j, pxa, thd = 256;

	cudaCalloc1((void**)&b, sizeof(float), npixels);
	cudaCalloc1((void**)&cosi, sizeof(float), npixels);
	cudaCalloc1((void**)&cose, sizeof(float), npixels);
	cudaCalloc1((void**)&zz, sizeof(float), npixels);

	BLK.x = floor((thd-1+npixels)/thd);	THD.x = thd;
	dbg_copy_pos_arrays_full_krnl<<<BLK,THD>>>(pos, f, npixels, b, cosi, cose, zz);
	checkErrorAfterKernelLaunch("dbg_copy_pos_arrays_full_krnl");
	deviceSyncAfterKernelLaunch("dbg_copy_pos_arrays_full_krnl");

	filename = "GPU_pos-b_s.csv";
	fp_b = fopen(filename, "w+");
	filename = "GPU_pos-cosi_s.csv";
	fp_cosi = fopen(filename, "w+");
	filename = "GPU_pos-cose_s.csv";
	fp_cose = fopen(filename, "w+");
	filename = "GPU_pos-z_s.csv";
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
__host__ void dbg_print_pos_arrays_full_host(struct pos_t *pos) {
	/* This debug function prints the CPU arrays:
	 *  - pos->cosi
	 *  - pos->cose
	 *  - pos->b,
	 *
	 *  all of length nPixels in the lghtcrv specified by 'set' in 'ddat' */
	int i, j, n;
	FILE *fp_cosi, *fp_cose, *fp_b, *fp_z;
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

	/* Print top corner set label */
	fprintf(fp_cosi, "set , ");
	fprintf(fp_cose, "set , ");
	fprintf(fp_b, "set , ");
	fprintf(fp_z, "set , ");

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
			fprintf(fp_b, "%g, ", pos->b[i][j]);
			fprintf(fp_z, "%g, ", pos->z[i][j]);
			fprintf(fp_cosi, "%g, ", pos->cosi[i][j]);
			fprintf(fp_cose, "%g, ", pos->cose[i][j]);
		}
		fprintf(fp_b, "\n");
		fprintf(fp_z, "\n");
		fprintf(fp_cosi, "\n");
		fprintf(fp_cose, "\n");
	}

	fclose(fp_b);
	fclose(fp_z);
	fclose(fp_cosi);
	fclose(fp_cose);
}
