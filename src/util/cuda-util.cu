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
	int offset = total_offset % frame_size;

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
__device__ void dev_mmmul2(double3 *x, double y[3][3], double3 *z, int frm)
{ /* This version turns the original double x[3][3] and double z[3][3] into
   * double3 pointers with nframes entries.  Selection is made via f  */
	double t[3][3];
	int i, f;
	f = 3*frm;

	for (i=0; i<=2; i++) {

		t[i][0] = 0.0;
		t[i][0] += y[i][0] * z[f+0].x;
		t[i][0] += y[i][1] * z[f+1].x;
		t[i][0] += y[i][2] * z[f+2].x;

		t[i][1] = 0.0;
		t[i][1] += y[i][0] * z[f+0].y;
		t[i][1] += y[i][1] * z[f+1].y;
		t[i][1] += y[i][2] * z[f+2].y;

		t[i][2] = 0.0;
		t[i][2] += y[i][0] * z[f+0].z;
		t[i][2] += y[i][1] * z[f+1].z;
		t[i][2] += y[i][2] * z[f+2].z;
	}

	for (i=0; i<=2; i++) {
		x[f+i].x = t[i][0];
		x[f+i].y = t[i][1];
		x[f+i].z = t[i][2];
	}
}
__device__ void dev_mmmul3(float3 *x, double y[3][3], float3 *z, int frm)
{ /* This version turns the original double x[3][3] and double z[3][3] into
   * double3 pointers with nframes entries.  Selection is made via f  */
	double t[3][3];
	int i, f;
	f = 3*frm;

	for (i=0; i<=2; i++) {

		t[i][0] = 0.0;
		t[i][0] += y[i][0] * z[f+0].x;
		t[i][0] += y[i][1] * z[f+1].x;
		t[i][0] += y[i][2] * z[f+2].x;

		t[i][1] = 0.0;
		t[i][1] += y[i][0] * z[f+0].y;
		t[i][1] += y[i][1] * z[f+1].y;
		t[i][1] += y[i][2] * z[f+2].y;

		t[i][2] = 0.0;
		t[i][2] += y[i][0] * z[f+0].z;
		t[i][2] += y[i][1] * z[f+1].z;
		t[i][2] += y[i][2] * z[f+2].z;
	}

	for (i=0; i<=2; i++) {
		x[f+i].x = t[i][0];
		x[f+i].y = t[i][1];
		x[f+i].z = t[i][2];
	}
}
__device__ int dev_vp_iround(double x)
{
  if (x < 0.0)
    return ((int)(x - 0.5));
  else
    return ((int)(x + 0.5));
}
__device__ int dev_vp_iroundf(float x)
{
  if (x < 0.0)
    return ((int)(x - 0.5));
  else
    return ((int)(x + 0.5));
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
__device__ double dev_dot2( double x[3], double3 y)
{
	/* This version replaces double y[3] with a double3 *y */
	return x[0]*y.x + x[1]*y.y + x[2]*y.z;
}
__device__ double dev_dot3( float3 x, double3 y)
{
	/* This version replaces double x[3] with a float3 *x */
	return x.x*y.x + x.y*y.y + x.z*y.z;
}
__device__ float dev_dot4(float3 x, float3 y)
{
	/* This version just uses two float3's and returns a float */
	return x.x*y.x + x.y*y.y + x.z*y.z;
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
__device__ float dev_normalize2(float3 u)
{
	int i;
	float norm;

	norm = 0.0;
	norm += u.x*u.x;
	norm += u.y*u.y;
	norm += u.z*u.z;

	norm = sqrt(norm);

	if (norm != 0.0) {
		u.x /= norm;
		u.y /= norm;
		u.z /= norm;
	}
	return norm;
}
__device__ void dev_cotrans1(double y[3], double *a, double x[3], int dir)
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
__device__ void dev_cotrans2(double y[3], double a[3][3], double x[3], int dir)
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
__device__ void dev_cotrans3(double y[3], double a[3][3], double x[3],
		int dir) {
	double t[3];
	int i, j;

	if (dir == 1)
		for (i = 0; i <= 2; i++) {
			t[i] = 0.0;
			for (j = 0; j <= 2; j++)
				t[i] += a[i][j] * x[j];
		}
	if (dir == (-1))
		for (i = 0; i <= 2; i++) {
			t[i] = 0.0;
			for (j = 0; j <= 2; j++)
				t[i] += a[j][i] * x[j];
		}
	for (i = 0; i <= 2; i++)
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
__device__ void dev_cotrans5(double3 *y, double a[3][3], double3 x,
		int dir) {
	/* This version replaces double y[3] and double x[3] with double3 y and double3 x */
	double t[3];
	int i;

	if (dir == 1)
		for (i = 0; i <= 2; i++) {
			t[i] = 0.0;
			t[i] += a[i][0] * x.x;
			t[i] += a[i][1] * x.y;
			t[i] += a[i][2] * x.z;
		}

	if (dir == (-1))
		for (i = 0; i <= 2; i++) {
			t[i] = 0.0;
			t[i] += a[0][i] * x.x;
			t[i] += a[1][i] * x.y;
			t[i] += a[2][i] * x.z;
		}

	y->x = t[0];
	y->y = t[1];
	y->z = t[2];
}
__device__ void dev_cotrans6(double y[3], double3 *a, double x[3], int dir, int frm) {
	/* This version replaces double a[3][3] with a double3 pointers of lenght
	 * nframes, selected with 'f'	 */
	double t[3];
	int i, j, f;
	f = frm*3;

	if (dir == 1)
		for (i = 0; i <= 2; i++) {
			t[i] = 0.0;
			t[i] += a[f+i].x * x[0];
			t[i] += a[f+i].y * x[1];
			t[i] += a[f+i].z * x[2];
		}

	if (dir == (-1)) {

		t[0] = 0.0;
		for (j=0; j<=2; j++) {
			t[0] += a[f+j].x * x[j];
			t[0] += a[f+j].x * x[j];
			t[0] += a[f+j].x * x[j];
		}
		t[1] = 0.0;
		for (j=0; j<=2; j++) {
			t[1] += a[f+j].y * x[j];
			t[1] += a[f+j].y * x[j];
			t[1] += a[f+j].y * x[j];
		}
		t[2] = 0.0;
		for (j=0; j<=2; j++) {
			t[2] += a[f+j].z * x[j];
			t[2] += a[f+j].z * x[j];
			t[2] += a[f+j].z * x[j];
		}
	}

	for (i = 0; i <= 2; i++)
		y[i] = t[i];
}
__device__ void dev_cotrans7(float3 *y, double3 *a, float3 x, int dir, int frm) {
	/* dev_cotrans7 is an iteration of dev_cotrans6.  It replaces the double y[3]
	 * with a float3.  It also replaces the double x[3] with a float3. It also
	 * replaces double3 *a with a float3 *a
	 * dev_cotrans6 - This version replaces double a[3][3] with a double3 pointers of lenght
	 * nframes, selected with 'f'	 */
	double3 t;
	int i, j, f;
	f = frm*3;

	if (dir == 1){
		t.x = 0.0;
		t.x += a[f+0].x * x.x;
		t.x += a[f+0].y * x.y;
		t.x += a[f+0].z * x.z;
		t.y = 0.0;
		t.y += a[f+1].x * x.x;
		t.y += a[f+1].y * x.y;
		t.y += a[f+1].z * x.z;
		t.z = 0.0;
		t.z += a[f+2].x * x.x;
		t.z += a[f+2].y * x.y;
		t.z += a[f+2].z * x.z;
	}

	if (dir == (-1)) {
		t.x = 0.0;
		t.x += a[f+0].x * x.x;
		t.x += a[f+1].x * x.y;
		t.x += a[f+2].x * x.z;
		t.y = 0.0;
		t.y += a[f+0].y * x.x;
		t.y += a[f+1].y * x.y;
		t.y += a[f+2].y * x.z;
		t.z = 0.0;
		t.z += a[f+0].z * x.x;
		t.z += a[f+1].z * x.y;
		t.z += a[f+2].z * x.z;
	}

	y->x = t.x;
	y->y = t.y;
	y->z = t.z;
}
__device__ void dev_cotrans8(float3 *y, float3 *a, float3 x, int dir, int frm)
{
	float3 t;
	int i, j, f=3*frm;

	if (dir==1) {
		t.x = 0.0;
		t.x += a[f+0].x * x.x;
		t.x += a[f+0].y * x.y;
		t.x += a[f+0].z * x.z;
		t.y = 0.0;
		t.y += a[f+1].x * x.x;
		t.y += a[f+1].y * x.y;
		t.y += a[f+1].z * x.z;
		t.z = 0.0;
		t.z += a[f+2].x * x.x;
		t.z += a[f+2].y * x.y;
		t.z += a[f+2].z * x.z;
	}
	if (dir==(-1)) {
		t.x = 0.0;
		t.x += a[f+0].x * x.x;
		t.x += a[f+1].x * x.y;
		t.x += a[f+2].x * x.z;
		t.y = 0.0;
		t.y += a[f+0].y * x.x;
		t.y += a[f+1].y * x.y;
		t.y += a[f+2].y * x.z;
		t.z = 0.0;
		t.z += a[f+0].z * x.x;
		t.z += a[f+1].z * x.y;
		t.z += a[f+2].z * x.z;
	}

	y->x = t.x;
	y->y = t.y;
	y->z = t.z;
}
__device__ void dev_cotrans9(float3 *y, double a[3][3], float3 x, int dir) {
	/* This version replaces double y[3] and double x[3] with double3 y and double3 x */
	double t[3];
	int i;

	if (dir == 1)
		for (i = 0; i <= 2; i++) {
			t[i] = 0.0;
			t[i] += a[i][0] * x.x;
			t[i] += a[i][1] * x.y;
			t[i] += a[i][2] * x.z;
		}

	if (dir == (-1))
		for (i = 0; i <= 2; i++) {
			t[i] = 0.0;
			t[i] += a[0][i] * x.x;
			t[i] += a[1][i] * x.y;
			t[i] += a[2][i] * x.z;
		}

	y->x = t[0];
	y->y = t[1];
	y->z = t[2];
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
__device__ void dev_mtrnsps2(double3 *a, double b[3][3], int frm)
{	/* This version splits the double a[3][3] of the original function into
	 * three separate double3 vector variables. b[3][3] remains unchanged. */
  double t[3][3];
  int i, j, f;
  f = frm *3;

  for (i=0;i<=2;i++)
	  for (j=0;j<=2;j++)
		  t[i][j] = b[j][i];

  a[f+0].x = t[0][0];
  a[f+0].y = t[0][1];
  a[f+0].z = t[0][2];
  a[f+1].x = t[1][0];
  a[f+1].y = t[1][1];
  a[f+1].z = t[1][2];
  a[f+2].x = t[2][0];
  a[f+2].y = t[2][1];
  a[f+2].z = t[2][2];
}
__device__ void dev_mtrnsps3(float3 *a, double b[3][3], int frm)
{	/* This version splits the double a[3][3] of the original function into
	 * three separate double3 vector variables. b[3][3] remains unchanged. */
  double t[3][3];
  int i, j, f;
  f = frm *3;

  for (i=0;i<=2;i++)
	  for (j=0;j<=2;j++)
		  t[i][j] = b[j][i];

  a[f+0].x = t[0][0];
  a[f+0].y = t[0][1];
  a[f+0].z = t[0][2];
  a[f+1].x = t[1][0];
  a[f+1].y = t[1][1];
  a[f+1].z = t[1][2];
  a[f+2].x = t[2][0];
  a[f+2].y = t[2][1];
  a[f+2].z = t[2][2];



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
__device__ double dev_radlaw_f( struct photo_t *photo, int ilaw, float cosinc, int c, int f)
{
  int i;
  float tan2inc, sin2inc, hagforsarg, incidence, rho1, rho2, angres, angratio;
  float diffxsec_proj = -9.99;  /* dummy value */

  switch (photo->radtype[ilaw]) {
  case COSINELAW_DIFF:
      diffxsec_proj = __double2float_rn(photo->radar[ilaw].RC.R.val)*
      	  	  	     (__double2float_rn(photo->radar[ilaw].RC.C.val) + 1)
                      * powf(cosinc, 2*__double2float_rn(photo->radar[ilaw].RC.C.val) - 1);
      break;
  case TABULARLAW:
      if (f < 0) {
          diffxsec_proj = 0.0;   /* blank sky */
      } else {
          incidence = acos(cosinc);
          angres = (PIE/2) / (photo->radar[ilaw].tabular.n - 1);
          angratio = incidence/angres;
          i = (int) floor(angratio);
          rho1 = __double2float_rn(photo->radar[ilaw].tabular.rho[i].val);
          rho2 = __double2float_rn(photo->radar[ilaw].tabular.rho[i+1].val);
          diffxsec_proj = (rho1 + (angratio - i)*(rho2 - rho1)) / cosinc;
      }
      break;
  case GAUSSIANLAW:
      if (cosinc >= photo->radar[ilaw].quasispec.cos_cutoff) {
          tan2inc = 1/(cosinc*cosinc) - 1;
          diffxsec_proj = __double2float_rn(photo->radar[ilaw].quasispec.R.val) *
        		  	  	  __double2float_rn(photo->radar[ilaw].quasispec.C.val) *
        		  	  	  exp(__double2float_rn(-photo->radar[ilaw].quasispec.C.val) *
        		  	  			  tan2inc) / (powf(cosinc,5));
      } else {
          diffxsec_proj = 0.0;
      }
      break;
  case HAGFORSLAW:
      if (cosinc >= photo->radar[ilaw].quasispec.cos_cutoff) {
          sin2inc = 1 - cosinc*cosinc;
          hagforsarg = powf(cosinc,4) + __double2float_rn(photo->radar[ilaw].quasispec.C.val)*sin2inc;
          diffxsec_proj = 0.5 * __double2float_rn(photo->radar[ilaw].quasispec.R.val) *
        		  __double2float_rn(photo->radar[ilaw].quasispec.C.val)
                          * powf( hagforsarg, -1.5) / cosinc;
      } else {
          diffxsec_proj = 0.0;
      }
      break;
  case COSINELAW_QS:
      if (cosinc >= photo->radar[ilaw].quasispec.cos_cutoff)
        diffxsec_proj = __double2float_rn(photo->radar[ilaw].quasispec.R.val) *
        (__double2float_rn(photo->radar[ilaw].quasispec.C.val) + 1)
        * pow( cosinc, 2*__double2float_rn(photo->radar[ilaw].quasispec.C.val) - 1);
      else
        diffxsec_proj = 0.0;
      break;
  case GAUSSIAN_COSINE:
      diffxsec_proj = __double2float_rn(photo->radar[ilaw].hybrid.diff.R.val) *
      	  (__double2float_rn(photo->radar[ilaw].hybrid.diff.C.val) + 1)
          * powf( cosinc, 2*__double2float_rn(photo->radar[ilaw].hybrid.diff.C.val) - 1);
      if (cosinc >= photo->radar[ilaw].hybrid.qs.cos_cutoff) {
        tan2inc = 1/(cosinc*cosinc) - 1;
        diffxsec_proj += __double2float_rn(photo->radar[ilaw].hybrid.qs.R.val) *
        		__double2float_rn(photo->radar[ilaw].hybrid.qs.C.val)
                         * exp(__double2float_rn(-photo->radar[ilaw].hybrid.qs.C.val) *
                        		 tan2inc)  / (pow(cosinc,5));
      }
      break;
  case HAGFORS_COSINE:
      diffxsec_proj = __double2float_rn(photo->radar[ilaw].hybrid.diff.R.val) *
      	  (__double2float_rn(photo->radar[ilaw].hybrid.diff.C.val) + 1) *
      	  	  pow( cosinc, 2*__double2float_rn(photo->radar[ilaw].hybrid.diff.C.val) - 1);
      if (cosinc >= photo->radar[ilaw].hybrid.qs.cos_cutoff) {
        sin2inc = 1 - cosinc*cosinc;
        hagforsarg = powf(cosinc,4) + __double2float_rn(photo->radar[ilaw].hybrid.qs.C.val)*sin2inc;
        diffxsec_proj += 0.5*__double2float_rn(photo->radar[ilaw].hybrid.qs.R.val) *
        		__double2float_rn(photo->radar[ilaw].hybrid.qs.C.val)
                         * powf( hagforsarg, -1.5) / cosinc;
      }
      break;
  case COSINE_COSINE:
      diffxsec_proj = __double2float_rn(photo->radar[ilaw].hybrid.diff.R.val) *
      	  (__double2float_rn(photo->radar[ilaw].hybrid.diff.C.val) + 1)
      	  * powf( cosinc, 2*__double2float_rn(photo->radar[ilaw].hybrid.diff.C.val) - 1);
      if (cosinc >= photo->radar[ilaw].hybrid.qs.cos_cutoff)
        diffxsec_proj += __double2float_rn(photo->radar[ilaw].hybrid.qs.R.val) *
        	(__double2float_rn(photo->radar[ilaw].hybrid.qs.C.val) + 1)
        	* powf( cosinc, 2*__double2float_rn(photo->radar[ilaw].hybrid.qs.C.val) - 1);
      break;
  case HARMCOSINE_DIFF:
      if (f < 0) {
          diffxsec_proj = 0.0;   /* blank sky */
      } else {
          diffxsec_proj = __double2float_rn(photo->radar[ilaw].harmcosine.local[c][f].R.val)
             * (__double2float_rn(photo->radar[ilaw].harmcosine.local[c][f].C.val) + 1)
             * powf( cosinc, 2*__double2float_rn(photo->radar[ilaw].harmcosine.local[c][f].C.val) - 1);
      }
      break;
  case INHOCOSINE_DIFF:
      if (f < 0) {
          diffxsec_proj = 0.0;   /* blank sky */
      } else {
          diffxsec_proj = __double2float_rn(photo->radar[ilaw].inhocosine.local[c][f].R.val)
              * (__double2float_rn(photo->radar[ilaw].inhocosine.local[c][f].C.val) + 1)
              * powf( cosinc, 2*__double2float_rn(photo->radar[ilaw].inhocosine.local[c][f].C.val) - 1);
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
