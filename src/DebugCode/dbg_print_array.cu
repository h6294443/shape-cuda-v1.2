
extern "C" {
#include "../shape/head.h"
}

__host__ void dbg_print_array(float *data, int x, int y) {
	/* Debug function that prints all elements in data to a csv in x col * y rows */

	int n, i, j;
	float *host;
	FILE *fp_fit;
	char *filename_fit;
	double *fit;
	dim3 BLK,THD;

	filename_fit = "dbg_array_cuda.csv";
	printf("\n %sfile created",filename_fit);
	printf("\n\nFilename: %s",filename_fit);

	n = x*y;
	host = (float *) malloc(n*sizeof(float));
	gpuErrchk(cudaMemcpy(host, data, n*sizeof(float), cudaMemcpyDeviceToHost));

	fp_fit = fopen(filename_fit, "w+");
	/* Print top corner idop/idel label */
	fprintf(fp_fit, "i/j , ");

	/* Print top row idel values */
	for (i=0; i<x; i++)
		fprintf(fp_fit, "%i , ", i);

	/* Print first entry in every row (except 1st): idop */
	for (j=1; j<y; j++) {
		fprintf(fp_fit,	"\n%i , ", j);

		/* Write the rest of the row values: fit[idel][idop] */
		for (i=0; i<x; i++)
			fprintf(fp_fit, " %g , ", host[j*x + i]);
	}
	fclose(fp_fit);
	free(host);
}
__host__ void dbg_print_array1D(float *data, int size) {
	/* Debug function that prints all elements in data to a csv */

	int i;
	float *host;
	FILE *fp_fit;
	char *filename_fit;
	double *fit;
	dim3 BLK,THD;

	filename_fit = "dbg_array1D_cuda.csv";
	printf("\n %sfile created",filename_fit);
	printf("\n\nFilename: %s",filename_fit);

	host = (float *) malloc(size*sizeof(float));
	gpuErrchk(cudaMemcpy(host, data, size*sizeof(float), cudaMemcpyDeviceToHost));

	fp_fit = fopen(filename_fit, "w+");
	/* Print top corner idop/idel label */
	fprintf(fp_fit, "i , \n");

	/* Print top row idel values */
	for (i=0; i<size; i++)
		fprintf(fp_fit, "%i , ", i);

	/* Go to second row */
	fprintf(fp_fit, "\n , ");

	/* Write the rest of the row values: fit[idel][idop] */
	for (i=0; i<size; i++)
		fprintf(fp_fit, " %g , ", host[i]);

	fclose(fp_fit);
	free(host);
}
__host__ void dbg_print_array1D_dbl(double *data, int size, int offset,
		char *filename) {
	/* Debug function that prints all elements in data to a csv */

	int i;
	FILE *fp_fit;

	printf("\n\nFilename: %s",filename);
	fp_fit = fopen(filename, "w+");

	/* Print top corner idop/idel label */
	fprintf(fp_fit, "i , \n");

	/* Print row */
	for (i=offset; i<=size; i++)
		fprintf(fp_fit, "%i , %g, \n", i, data[i]);

	fclose(fp_fit);
}
__host__ void dbg_print_array1(float *in, int size) {
	/* This debug function prints each array value */
	int i;

	for (i=0; i<size; i++) {
		printf("\narray[%i]=%g", i, in[i]);
	}
}
