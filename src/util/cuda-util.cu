extern "C" {
#include "../shape/head.h"
}
__global__ void clrvect_krnl(struct dat_t *ddat, int s, int f, int nThreads) {
	/* multi-threaded kernel */
	int offset = blockIdx.x * blockDim.x + threadIdx.x;

	if (offset < nThreads){
		switch (ddat->set[s].type) {
		case DELAY:
			ddat->set[s].desc.deldop.frame[f].fit_s[offset] = 0.0;
			break;
		case DOPPLER:
			ddat->set[s].desc.doppler.frame[f].fit_s[offset] = 0.0;
			break;
		}
	}
}
__global__ void clrvect_af_krnl(struct dat_t *ddat, int s, int nframes,
		int nThreads, int frame_size) {
	/* multi-threaded kernel for all frames in a set */
	int total_offset = blockIdx.x * blockDim.x + threadIdx.x;
	int frm = total_offset / frame_size;
	int offset = offset % frame_size;

	if ((offset < nThreads) && (frm < nframes)) {
		switch (ddat->set[s].type) {
		case DELAY:
			ddat->set[s].desc.deldop.frame[frm].fit_s[offset] = 0.0;
			break;
		case DOPPLER:
			ddat->set[s].desc.doppler.frame[frm].fit_s[offset] = 0.0;
			break;
		}
	}
}

void cotrans_cuda(double y[3], double a[3][3], double x[3], int dir)
{
	double t[3];
	int i, j;

	if (dir==1)
		for (i=0;i<=2;i++) {
			t[i] = 0.0;
			for (j=0;j<=2;j++)
				t[i] += a[i][j]*x[j];
		}
	if (dir==(-1))
		for (i=0;i<=2;i++) {
			t[i] = 0.0;
			for (j=0;j<=2;j++)
				t[i] += a[j][i]*x[j];
		}
	for (i=0;i<=2;i++)
		y[i] = t[i];
}

void mmmul_cuda( double *x, double y[3][3], double *z)
{
	double t[3][3];
	int i, j, k;

	for (i=0;i<=2;i++)
		for (j=0;j<=2;j++) {
			t[i][j] = 0.0;
			for (k=0;k<=2;k++)
				t[i][j] += y[i][k]*z[k*3+j];
		}
	for (i=0;i<=2;i++)
		for (j=0;j<=2;j++)
			x[i*3+j] = t[i][j];
}

__device__ void dev_mmmul( double x[3][3], double y[3][3], double z[3][3])
{
  double t[3][3];
  int i, j, k;

  for (i=0;i<=2;i++)
	for (j=0;j<=2;j++) {
	  t[i][j] = 0.0;
	  for (k=0;k<=2;k++)
		t[i][j] += y[i][k]*z[k][j];
	}
  for (i=0;i<=2;i++)
	for (j=0;j<=2;j++)
	  x[i][j] = t[i][j];
}


void mtrnsps_cuda( double *a, double b[3][3])
{
	double t[3][3];
	int i, j;

	for (i=0;i<=2;i++)
		for (j=0;j<=2;j++)
			t[i][j] = b[j][i];
	for (i=0;i<=2;i++)
		for (j=0;j<=2;j++)
			a[i*3+j] = t[i][j];
	//a[i][j] = t[i][j];
}

void checkErrorAfterKernelLaunch(char *location) {
	cudaError_t cudaStatus;
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Kernel launch failed in %s: %s\n", location, cudaGetErrorString(cudaStatus));
	}
}

void deviceSyncAfterKernelLaunch(char *location) {
	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaError_t cudaStatus;
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess)
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching the kernel in %s.\n", cudaStatus, location);
}

__device__ double dev_dot( double x[3], double y[3])
{
	return x[0]*y[0]+x[1]*y[1]+x[2]*y[2];
}

__device__ double dev_normalize(double *u)
{
	int i;
	double norm;

	norm = 0.0;
	for (i=0; i<=2; i++)
		norm += u[i]*u[i];
	norm = sqrt(norm);
	if (norm != 0.0) {
		for (i=0; i<=2; i++)
			u[i] /= norm;
	}
	return norm;
}

__device__ void dev_cotrans1( double y[3], double *a, double x[3], int dir)
{
	double t[3];
	int i, j;

	if (dir==1)
		for (i=0;i<=2;i++) {
			t[i] = 0.0;
			for (j=0;j<=2;j++)
				t[i] += a[3*j+i]*x[j];
		}
	if (dir==(-1))
		for (i=0;i<=2;i++) {
			t[i] = 0.0;
			for (j=0;j<=2;j++)
				t[i] += a[3*i+j]*x[j];
		}
	for (i=0;i<=2;i++)
		y[i] = t[i];
}
__device__ void dev_cotrans2( double y[3], double a[3][3], double x[3], int dir)
{
	double t[3];
	int i, j;

	if (dir==1)
		for (i=0;i<=2;i++) {
			t[i] = 0.0;
			for (j=0;j<=2;j++)
				t[i] += a[i][j]*x[j];
		}
	if (dir==(-1))
		for (i=0;i<=2;i++) {
			t[i] = 0.0;
			for (j=0;j<=2;j++)
				t[i] += a[j][i]*x[j];
		}
	for (i=0;i<=2;i++)
		y[i] = t[i];
}
__device__ void dev_cotrans4(float3 *y, double a[3][3], double x[3], int dir, int f)
{
	double t[3];
	int i, j;

	if (dir==1)
		for (i=0;i<=2;i++) {
			t[i] = 0.0;
			for (j=0;j<=2;j++)
				t[i] += a[i][j]*x[j];
		}
	if (dir==(-1))
		for (i=0;i<=2;i++) {
			t[i] = 0.0;
			for (j=0;j<=2;j++)
				t[i] += a[j][i]*x[j];
		}
	y[f].x = t[0];
	y[f].y = t[1];
	y[f].z = t[2];
}
__device__ void dev_mtrnsps( double a[3][3], double b[3][3])
{
  double t[3][3];
  int i, j;

  for (i=0;i<=2;i++)
	for (j=0;j<=2;j++)
	  t[i][j] = b[j][i];
  for (i=0;i<=2;i++)
	for (j=0;j<=2;j++)
	  a[i][j] = t[i][j];
}
__device__ double radlaw_cuda(union radscat_t *radar, unsigned char *radtype,
		int ilaw, double cosinc, int c, int f)
{
	int irdl;
	double tan2inc, sin2inc, hagforsarg, incidence, rho1, rho2, angres, angratio;
	double diffxsec_proj = -9.99;  /* dummy value */

	switch (radtype[ilaw]) {
	case COSINELAW_DIFF:
		diffxsec_proj = radar[ilaw].RC.R.val*(radar[ilaw].RC.C.val + 1)
		* pow( cosinc, 2*radar[ilaw].RC.C.val - 1);
		break;
	case TABULARLAW:
		if (f < 0) {
			diffxsec_proj = 0.0;   /* blank sky */
		} else {
			incidence = acos(cosinc);
			angres = (PIE/2) / (radar[ilaw].tabular.n - 1);
			angratio = incidence/angres;
			irdl = (int) floor(angratio);
			rho1 = radar[ilaw].tabular.rho[irdl].val;
			rho2 = radar[ilaw].tabular.rho[irdl+1].val;
			diffxsec_proj = (rho1 + (angratio - irdl)*(rho2 - rho1)) / cosinc;
		}
		break;
	case GAUSSIANLAW:
		if (cosinc >= radar[ilaw].quasispec.cos_cutoff) {
			tan2inc = 1/(cosinc*cosinc) - 1;
			diffxsec_proj = radar[ilaw].quasispec.R.val*radar[ilaw].quasispec.C.val
					* exp( -radar[ilaw].quasispec.C.val * tan2inc)
			/ (cosinc*cosinc*cosinc*cosinc*cosinc);
		} else {
			diffxsec_proj = 0.0;
		}
		break;
	case HAGFORSLAW:
		if (cosinc >= radar[ilaw].quasispec.cos_cutoff) {
			sin2inc = 1 - cosinc*cosinc;
			hagforsarg = cosinc*cosinc*cosinc*cosinc + radar[ilaw].quasispec.C.val*sin2inc;
			diffxsec_proj = 0.5*radar[ilaw].quasispec.R.val*radar[ilaw].quasispec.C.val
					* pow( hagforsarg, -1.5) / cosinc;
		} else {
			diffxsec_proj = 0.0;
		}
		break;
	case COSINELAW_QS:
		if (cosinc >= radar[ilaw].quasispec.cos_cutoff)
			diffxsec_proj = radar[ilaw].quasispec.R.val*(radar[ilaw].quasispec.C.val + 1)
			* pow( cosinc, 2*radar[ilaw].quasispec.C.val - 1);
		else
			diffxsec_proj = 0.0;
		break;
	case GAUSSIAN_COSINE:
		diffxsec_proj = radar[ilaw].hybrid.diff.R.val*(radar[ilaw].hybrid.diff.C.val + 1)
		* pow( cosinc, 2*radar[ilaw].hybrid.diff.C.val - 1);
		if (cosinc >= radar[ilaw].hybrid.qs.cos_cutoff) {
			tan2inc = 1/(cosinc*cosinc) - 1;
			diffxsec_proj += radar[ilaw].hybrid.qs.R.val*radar[ilaw].hybrid.qs.C.val
					* exp( -radar[ilaw].hybrid.qs.C.val * tan2inc)
			/ (cosinc*cosinc*cosinc*cosinc*cosinc);
		}
		break;
	case HAGFORS_COSINE:
		diffxsec_proj = radar[ilaw].hybrid.diff.R.val*(radar[ilaw].hybrid.diff.C.val + 1)
		* pow( cosinc, 2*radar[ilaw].hybrid.diff.C.val - 1);
		if (cosinc >= radar[ilaw].hybrid.qs.cos_cutoff) {
			sin2inc = 1 - cosinc*cosinc;
			hagforsarg = cosinc*cosinc*cosinc*cosinc + radar[ilaw].hybrid.qs.C.val*sin2inc;
			diffxsec_proj += 0.5*radar[ilaw].hybrid.qs.R.val*radar[ilaw].hybrid.qs.C.val
					* pow( hagforsarg, -1.5) / cosinc;
		}
		break;
	case COSINE_COSINE:
		diffxsec_proj = radar[ilaw].hybrid.diff.R.val*(radar[ilaw].hybrid.diff.C.val + 1)
		* pow( cosinc, 2*radar[ilaw].hybrid.diff.C.val - 1);
		if (cosinc >= radar[ilaw].hybrid.qs.cos_cutoff)
			diffxsec_proj += radar[ilaw].hybrid.qs.R.val*(radar[ilaw].hybrid.qs.C.val + 1)
			* pow( cosinc, 2*radar[ilaw].hybrid.qs.C.val - 1);
		break;
	case HARMCOSINE_DIFF:
		if (f < 0) {
			diffxsec_proj = 0.0;   /* blank sky */
		} else {
			diffxsec_proj = radar[ilaw].harmcosine.local[c][f].R.val
					* (radar[ilaw].harmcosine.local[c][f].C.val + 1)
					* pow( cosinc, 2*radar[ilaw].harmcosine.local[c][f].C.val - 1);
		}
		break;
	case INHOCOSINE_DIFF:
		if (f < 0) {
			diffxsec_proj = 0.0;   /* blank sky */
		} else {
			diffxsec_proj = radar[ilaw].inhocosine.local[c][f].R.val
					* (radar[ilaw].inhocosine.local[c][f].C.val + 1)
					* pow( cosinc, 2*radar[ilaw].inhocosine.local[c][f].C.val - 1);
		}
		break;
	case NOLAW:
		printf("\n\npos2doppler-cuda.cu: can't set radar scattering law.\n\n");
	default:
		printf("\n\npos2doppler-cuda.cu: Unspecified error.\n\n");
	}
	return diffxsec_proj;
}

__device__ double dev_radlaw( struct photo_t *photo, int ilaw, double cosinc, int c, int f)
{
  int i;
  double tan2inc, sin2inc, hagforsarg, incidence, rho1, rho2, angres, angratio;
  double diffxsec_proj = -9.99;  /* dummy value */

  switch (photo->radtype[ilaw]) {
  case COSINELAW_DIFF:
      diffxsec_proj = photo->radar[ilaw].RC.R.val*(photo->radar[ilaw].RC.C.val + 1)
                      * pow( cosinc, 2*photo->radar[ilaw].RC.C.val - 1);
      break;
  case TABULARLAW:
      if (f < 0) {
          diffxsec_proj = 0.0;   /* blank sky */
      } else {
          incidence = acos(cosinc);
          angres = (PIE/2) / (photo->radar[ilaw].tabular.n - 1);
          angratio = incidence/angres;
          i = (int) floor(angratio);
          rho1 = photo->radar[ilaw].tabular.rho[i].val;
          rho2 = photo->radar[ilaw].tabular.rho[i+1].val;
          diffxsec_proj = (rho1 + (angratio - i)*(rho2 - rho1)) / cosinc;
      }
      break;
  case GAUSSIANLAW:
      if (cosinc >= photo->radar[ilaw].quasispec.cos_cutoff) {
          tan2inc = 1/(cosinc*cosinc) - 1;
          diffxsec_proj = photo->radar[ilaw].quasispec.R.val*photo->radar[ilaw].quasispec.C.val
                          * exp( -photo->radar[ilaw].quasispec.C.val * tan2inc)
                          / (cosinc*cosinc*cosinc*cosinc*cosinc);
      } else {
          diffxsec_proj = 0.0;
      }
      break;
  case HAGFORSLAW:
      if (cosinc >= photo->radar[ilaw].quasispec.cos_cutoff) {
          sin2inc = 1 - cosinc*cosinc;
          hagforsarg = cosinc*cosinc*cosinc*cosinc + photo->radar[ilaw].quasispec.C.val*sin2inc;
          diffxsec_proj = 0.5*photo->radar[ilaw].quasispec.R.val*photo->radar[ilaw].quasispec.C.val
                          * pow( hagforsarg, -1.5) / cosinc;
      } else {
          diffxsec_proj = 0.0;
      }
      break;
  case COSINELAW_QS:
      if (cosinc >= photo->radar[ilaw].quasispec.cos_cutoff)
        diffxsec_proj = photo->radar[ilaw].quasispec.R.val*(photo->radar[ilaw].quasispec.C.val + 1)
                        * pow( cosinc, 2*photo->radar[ilaw].quasispec.C.val - 1);
      else
        diffxsec_proj = 0.0;
      break;
  case GAUSSIAN_COSINE:
      diffxsec_proj = photo->radar[ilaw].hybrid.diff.R.val*(photo->radar[ilaw].hybrid.diff.C.val + 1)
                      * pow( cosinc, 2*photo->radar[ilaw].hybrid.diff.C.val - 1);
      if (cosinc >= photo->radar[ilaw].hybrid.qs.cos_cutoff) {
        tan2inc = 1/(cosinc*cosinc) - 1;
        diffxsec_proj += photo->radar[ilaw].hybrid.qs.R.val*photo->radar[ilaw].hybrid.qs.C.val
                         * exp( -photo->radar[ilaw].hybrid.qs.C.val * tan2inc)
                         / (cosinc*cosinc*cosinc*cosinc*cosinc);
      }
      break;
  case HAGFORS_COSINE:
      diffxsec_proj = photo->radar[ilaw].hybrid.diff.R.val*(photo->radar[ilaw].hybrid.diff.C.val + 1)
                      * pow( cosinc, 2*photo->radar[ilaw].hybrid.diff.C.val - 1);
      if (cosinc >= photo->radar[ilaw].hybrid.qs.cos_cutoff) {
        sin2inc = 1 - cosinc*cosinc;
        hagforsarg = cosinc*cosinc*cosinc*cosinc + photo->radar[ilaw].hybrid.qs.C.val*sin2inc;
        diffxsec_proj += 0.5*photo->radar[ilaw].hybrid.qs.R.val*photo->radar[ilaw].hybrid.qs.C.val
                         * pow( hagforsarg, -1.5) / cosinc;
      }
      break;
  case COSINE_COSINE:
      diffxsec_proj = photo->radar[ilaw].hybrid.diff.R.val*(photo->radar[ilaw].hybrid.diff.C.val + 1)
                      * pow( cosinc, 2*photo->radar[ilaw].hybrid.diff.C.val - 1);
      if (cosinc >= photo->radar[ilaw].hybrid.qs.cos_cutoff)
        diffxsec_proj += photo->radar[ilaw].hybrid.qs.R.val*(photo->radar[ilaw].hybrid.qs.C.val + 1)
                         * pow( cosinc, 2*photo->radar[ilaw].hybrid.qs.C.val - 1);
      break;
  case HARMCOSINE_DIFF:
      if (f < 0) {
          diffxsec_proj = 0.0;   /* blank sky */
      } else {
          diffxsec_proj = photo->radar[ilaw].harmcosine.local[c][f].R.val
                          * (photo->radar[ilaw].harmcosine.local[c][f].C.val + 1)
                          * pow( cosinc, 2*photo->radar[ilaw].harmcosine.local[c][f].C.val - 1);
      }
      break;
  case INHOCOSINE_DIFF:
      if (f < 0) {
          diffxsec_proj = 0.0;   /* blank sky */
      } else {
          diffxsec_proj = photo->radar[ilaw].inhocosine.local[c][f].R.val
                          * (photo->radar[ilaw].inhocosine.local[c][f].C.val + 1)
                          * pow( cosinc, 2*photo->radar[ilaw].inhocosine.local[c][f].C.val - 1);
      }
      break;
  case NOLAW:
      printf("radlaw.c: can't set radar scattering law = \"none\" when radar data are used\n");
      break;
  default:
      printf("radlaw.c: can't handle that radar scattering law yet\n");
  }

  return diffxsec_proj;
}

