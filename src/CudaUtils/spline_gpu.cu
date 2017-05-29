extern "C" {
#include "../shape/head.h"
}

__global__ void lghtcrv_spline_streams_krnl(struct dat_t *ddat, int set, double
		yp1, double ypn, double *u, int ncalc) {
	/*(double *x  - lghtcrv->x
	 * double *y  - lghtcrv->y
	 * int n      - calc
	 * double yp1 - 2.0e30
	 * double ypn - 2.0e30
	 * double *y2 - lghtcrv->y2)*/

	/* Multi-threaded kernel */
	int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
	int n = ncalc;
	int k = n - 1 - i;
	double *x = ddat->set[set].desc.lghtcrv.x;
	double *y = ddat->set[set].desc.lghtcrv.y;
	double *y2 = ddat->set[set].desc.lghtcrv.y2;
	double p, qn, sig, un;

	if (i == 1) {
		if (yp1 > 0.99e30)
			ddat->set[set].desc.lghtcrv.y2[1] = u[1] = 0.0;
		else {
			ddat->set[set].desc.lghtcrv.y2[1] = -0.5;
			u[1] = (3.0/(ddat->set[set].desc.lghtcrv.x[2]-ddat->set[set].desc.lghtcrv.x[1])) *
					((ddat->set[set].desc.lghtcrv.y[2]-ddat->set[set].desc.lghtcrv.y[1]) /
					 (ddat->set[set].desc.lghtcrv.x[2]-ddat->set[set].desc.lghtcrv.x[1])-yp1);
		}
	}
	__syncthreads();

	if ((i>=2) && (i<=(n-1))) {
		sig = (ddat->set[set].desc.lghtcrv.x[i] - ddat->set[set].desc.lghtcrv.x[i-1]) /
				(ddat->set[set].desc.lghtcrv.x[i+1] - ddat->set[set].desc.lghtcrv.x[i-1]);

		p = sig * ddat->set[set].desc.lghtcrv.y2[i-1] + 2.0;

		ddat->set[set].desc.lghtcrv.y2[i] = (sig-1.0)/p;

		u[i] = (ddat->set[set].desc.lghtcrv.y[i+1] - ddat->set[set].desc.lghtcrv.y[i]) /
				(ddat->set[set].desc.lghtcrv.x[i+1] - ddat->set[set].desc.lghtcrv.x[i]) -
				(ddat->set[set].desc.lghtcrv.y[i] - ddat->set[set].desc.lghtcrv.y[i-1]) /
				(ddat->set[set].desc.lghtcrv.x[i] - ddat->set[set].desc.lghtcrv.x[i-1]);

		u[i] = (6.0*u[i]/(ddat->set[set].desc.lghtcrv.x[i+1] - ddat->set[set].desc.lghtcrv.x[i-1]) -
				sig * u[i-1])/p;
	}
	__syncthreads();

	if (i == n) {
		if (ypn > 0.99e30)
			qn = un = 0.0;
		else {
			qn = 0.5;
			un = (3.0/(ddat->set[set].desc.lghtcrv.x[n] - ddat->set[set].desc.lghtcrv.x[n-1])) *
				(ypn-(ddat->set[set].desc.lghtcrv.y[n] - ddat->set[set].desc.lghtcrv.y[n-1]) /
				(ddat->set[set].desc.lghtcrv.x[n] - ddat->set[set].desc.lghtcrv.x[n-1]));
		}

		ddat->set[set].desc.lghtcrv.y2[n] = (un - qn * u[n-1]) /
				(qn * ddat->set[set].desc.lghtcrv.y2[n-1] + 1.0);

		for (k=n-1;k>=1;k--)
			ddat->set[set].desc.lghtcrv.y2[k] = ddat->set[set].desc.lghtcrv.y2[k] *
			ddat->set[set].desc.lghtcrv.y2[k+1] + u[k];
	}
}
__global__ void lghtcrv_spline_streams_f_krnl(struct dat_t *ddat, int set, float
		yp1, float ypn, float *u, int ncalc) {
	/*(double *x  - lghtcrv->x
	 * double *y  - lghtcrv->y
	 * int n      - calc
	 * double yp1 - 2.0e30
	 * double ypn - 2.0e30
	 * double *y2 - lghtcrv->y2)*/

	/* Multi-threaded kernel */
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int n = ncalc;
	int k = n - 1 - i;

	/* The following will be turned into floats */
	float *x = ddat->set[set].desc.lghtcrv.x_s;
	float *y = ddat->set[set].desc.lghtcrv.y_s;
	float *y2 = ddat->set[set].desc.lghtcrv.y2_s;
	float p,qn,sig,un;

	/* Perform single-thread task */
	if (i == 0) {
		if (yp1 > 0.99e30)
			y2[1]=u[1]=0.0;
		else {
			y2[1] = -0.5;
			u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
		}
	}
	__syncthreads();

	if (i > 1 && i < n-1) {

		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	__syncthreads();

	/* Perform another single-thread task */
	if (i == 1) {
		if (ypn > 0.99e30)
			qn=un=0.0;
		else {
			qn=0.5;
			un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
		}
		y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
	}
	__syncthreads();

	if (k <= (n-1) && k >= 1)
		y2[k]=y2[k]*y2[k+1]+u[k];

	__syncthreads();
}
__global__ void lghtcrv_spline_streams_test_krnl(struct dat_t *ddat, int set, double
		yp1, double ypn, double *u, int ncalc) {
	/*(double *x  - lghtcrv->x
	 * double *y  - lghtcrv->y
	 * int n      - calc
	 * double yp1 - 2.0e30
	 * double ypn - 2.0e30
	 * double *y2 - lghtcrv->y2)*/

	int i, k, n=ncalc;
	double p, qn, sig, un;
	double *x = ddat->set[set].desc.lghtcrv.x;
	double *y = ddat->set[set].desc.lghtcrv.y;
	double *y2 = ddat->set[set].desc.lghtcrv.y2;

	/* single-threaded kernel */

	if (threadIdx.x == 0) {
		u[0]=0.0;
		if (yp1 > 0.99e30)
			y2[1]=u[1]=0.0;
		else {
			y2[1] = -0.5;
			u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
		}

		for (i=2;i<=n-1;i++) {
			sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
			p=sig*y2[i-1]+2.0;
			y2[i]=(sig-1.0)/p;
			u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);

			u[i] = (6.0 * u[i] / (x[i+1]-x[i-1]) - sig*u[i-1] ) / p;
		}

		if (ypn > 0.99e30)
			qn=un=0.0;
		else {
			qn=0.5;
			un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
		}
		y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);

		for (k=n-1;k>=1;k--)
			y2[k]=y2[k]*y2[k+1]+u[k];
	}
}
