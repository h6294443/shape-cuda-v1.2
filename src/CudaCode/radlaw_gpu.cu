extern "C" {
#include "../shape/head.h"
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
