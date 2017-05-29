extern "C" {
#include "../shape/head.h"
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
__device__ void dev_cotrans3(float y[3], float3 *a, float x[3], int dir, int frm) {
	/* This version replaces double a[3][3] with a double3 pointers of lenght
	 * nframes, selected with 'f'	 */
	float t[3];
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
