
extern "C" {
#include "../shape/head.h"
}

__global__ void dbg_copy_facet_normals_krnl(struct mod_t *dmod, int nf, float3 *dnormals)
{
	/* nf-threaded kernel */
	int f = blockIdx.x * blockDim.x + threadIdx.x;

	if (f<nf) {
		dnormals[f].x = __double2float_rn(dmod->shape.comp[0].real.f[f].n[0]);
		dnormals[f].y = __double2float_rn(dmod->shape.comp[0].real.f[f].n[1]);
		dnormals[f].z = __double2float_rn(dmod->shape.comp[0].real.f[f].n[2]);
	}
}
__host__ void dbg_print_facet_normals_host(struct mod_t *mod, const char *fn) {
	/* This debug function prints all facet normals in a given model */
	int nf;
	FILE *fp_n;
	nf = mod->shape.comp[0].real.nf;
	fp_n = fopen(fn, "w+");

	/* Print top row */
	fprintf(fp_n, ", value, \n");

	for (int f=0; f<nf; f++) {
		fprintf(fp_n, "%i, %g, \n", f, mod->shape.comp[0].real.f[f].n[0]);
		fprintf(fp_n, "%i, %g, \n", f, mod->shape.comp[0].real.f[f].n[1]);
		fprintf(fp_n, "%i, %g, \n", f, mod->shape.comp[0].real.f[f].n[2]);
	}
	fclose(fp_n);
}
__host__ void dbg_print_facet_normals(struct mod_t *dmod, int nf, const char *fn) {
	/* This debug function prints all facet normals in a given model */
	FILE *fp_n;
	float3 *dnormals, *hnormals;
	dim3 BLK,THD;
	fp_n = fopen(fn, "w+");

	/* Allocate memory */
	gpuErrchk(cudaMalloc((void**)&dnormals, sizeof(float3) * nf));
	hnormals = (float3 *) malloc(nf*sizeof(float3));

	THD.x = maxThreadsPerBlock;
	BLK.x = floor((THD.x - 1 + nf)/THD.x);
	dbg_copy_facet_normals_krnl<<<BLK,THD>>>(dmod, nf, dnormals);
	checkErrorAfterKernelLaunch("copy_facet_normals_krnl");
	gpuErrchk(cudaMemcpy(hnormals, dnormals, sizeof(float3)*nf, cudaMemcpyDeviceToHost));

	/* Print top row */
	fprintf(fp_n, ", value, \n");

	for (int f=0; f<nf; f++) {
		fprintf(fp_n, "%i, %g, \n", f, hnormals[f].x);
		fprintf(fp_n, "%i, %g, \n", f, hnormals[f].y);
		fprintf(fp_n, "%i, %g, \n", f, hnormals[f].z);
	}
	fclose(fp_n);
}
