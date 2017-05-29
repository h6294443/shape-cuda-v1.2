extern "C" {
#include "../shape/head.h"
}

#define IMAX 11
#define NUSE 7
#define SHRINK 0.95
#define GROW 1.2
#define nv 12	// new to the CUDA version, this is to allow the array declaration at compile time

__device__ double dd[nv][NUSE],devx[IMAX];	/* defining declaration */

__device__ void  dev_bsstep(double *y, double *dydx, int nv1, double *xx, double htry,
             double eps, double *yscal, double *hdid, double *hnext, 
             void (*dev_derivs)(double,double *,double *))
{
	int i,j;

	double xsav,xest,h,errmax,temp;
	double ysav[nv],dysav[nv],yseq[nv],yerr[nv];
	static int nseq[IMAX+1]={0,2,4,6,8,12,16,24,32,48,64,96};

	h=htry;
	xsav=(*xx);
	for (i=1;i<=nv;i++) {
		ysav[i]=y[i];
		dysav[i]=dydx[i];
	}
	for (;;) {
		for (i=1;i<=IMAX;i++) {
			dev_mmid(ysav,dysav,nv,xsav,h,nseq[i],yseq,dev_derivs);
			xest=(temp=h/nseq[i],temp*temp);
			dev_rzextr(i,xest,yseq,y,yerr);
			errmax=0.0;
			for (j=1;j<=nv;j++)
				if (errmax < fabs(yerr[j]/yscal[j]))
					errmax=fabs(yerr[j]/yscal[j]);
			errmax /= eps;
			if (errmax < 1.0) {
				*xx += h;
				*hdid=h;
				*hnext = i==NUSE? h*SHRINK : i==NUSE-1?
					h*GROW : (h*nseq[NUSE-1])/nseq[i];
				return;
			}
		}
		h *= 0.25;
		for (i=1;i<=(IMAX-NUSE)/2;i++) h /= 2.0;
		if ((*xx+h) == (*xx)) printf("Step size underflow in BSSTEP");
	}
}

__device__ void dev_rzextr( int iest, double xest, double *yest, double *yz, double *dy)
{
	int m1,k,j;
	double yy,v,ddy,c,b1,b,fx[NUSE];

	devx[iest]=xest;
	if (iest == 1)
		for (j=1;j<=nv;j++) {
			yz[j]=yest[j];
			dd[j][1]=yest[j];
			dy[j]=yest[j];
		}
	else {
		m1=(iest < NUSE ? iest : NUSE);
		for (k=1;k<=m1-1;k++)
			fx[k+1]=devx[iest-k]/xest;
		for (j=1;j<=nv;j++) {
			yy=yest[j];
			v=dd[j][1];
			c=yy;
			dd[j][1]=yy;
			for (k=2;k<=m1;k++) {
				b1=fx[k]*v;
				b=b1-c;
				if (b) {
					b=(c-v)/b;
					ddy=c*b;
					c=b1*b;
				} else
					ddy=v;
				if (k != m1) v=dd[j][k];
				dd[j][k]=ddy;
				yy += ddy;
			}
			dy[j]=ddy;
			yz[j]=yy;
		}
	}
}



#undef IMAX
#undef NUSE
#undef SHRINK
#undef GROW
