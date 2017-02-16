/* Modified 2011 August 12 by CM: increased MAXSTP to 1000000; set yscal[i] = 1.0 for all i
       -- that is, treat eps as an absolute tolerance rather than fractional (see discussion
       in Numerical Recipes in C) on the assumption that all variables y are of order unity  */
/* Modified 2011 August 10 by CM: increased MAXSTP to 100000  */
extern "C" {

#include "head.h"
}

#define MAXSTP 1000000
#define TINY 1.0e-30
#define nvar 12 //new

__device__ int dkmax=0,dkount=0;  /* defining declaration */

//** FIX THESE POINTER LIMITS **//
__device__ double *dxp=0,*dyp[12],ddxsav=0;  /* defining declaration */
//** FIX THESE POINTER LIMITS **//

/* Note that the second argument - originally nvar, now nvar1, has been fixed to 12 here */

__device__ void dev_odeint( double ystart[12], int nvar1, double x1, double x2, double eps,
			double h1, double hmin, int *nok, int *nbad,
			void (*dev_derivs)(double,double *,double *),
			void (*drkqc)(double *,double *,int,double *,double,double,double *,double *,
					double *,void (*)(double,double *,double *)))
{ 	
  int nstp,i;
  double xsav,x,hnext,hdid,h; 	
  double yscal[nvar], y[nvar], dydx[nvar];

  x = x1;
  h = (x2 > x1) ? fabs(h1) : -fabs(h1);
  *nok = (*nbad) = dkount = 0;

  for (i=1;i<=nvar;i++)
	  y[i]=ystart[i];

  if (dkmax > 0)
	  xsav=x-ddxsav*2.0;

  for (nstp=1;nstp<=MAXSTP;nstp++) {
	(*dev_derivs)(x,y,dydx);
	for (i=1;i<=nvar;i++)
	  /* yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY; */
	  yscal[i]=1.0;

	if (dkmax > 0) {
	  if (fabs(x-xsav) > fabs(ddxsav)) {
		if (dkount < dkmax-1) {
		  dxp[++dkount]=x;
		  for (i=1;i<=nvar;i++)
			  dyp[i][dkount]=y[i];
		  xsav=x;
		}
	  }
	}
	if ((x+h-x2)*(x+h-x1) > 0.0)
		h=x2-x;

	(*drkqc)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,dev_derivs);
	if (hdid == h)
		++(*nok);
	else
		++(*nbad);

	if ((x-x2)*(x2-x1) >= 0.0) {
	  for (i=1;i<=nvar;i++) ystart[i]=y[i];
	  if (dkmax) {
		dxp[++dkount]=x;
		for (i=1;i<=nvar;i++)
			dyp[i][dkount]=y[i];
	  }
	  return;
	}
	if (fabs(hnext) <= hmin) printf("Step size too small in ODEINT");
	h=hnext;
  }
  printf("Too many steps in routine ODEINT");
}

#undef MAXSTP
#undef TINY
#undef nvar
