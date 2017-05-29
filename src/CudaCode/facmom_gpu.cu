/***************************************************************************
                                                                   facmom.c

Compute a triangular facet's contribution to a model component's
0, 1, and 2-order moments, under the assumption of uniform density.

There are four input parameters:

    1-3: The body coordinates of the three vertices which define the facet,
         fv0, fv1, and fv2

    4:   The unit vector (in body coordinates) normal to the facet, fn

         This normal is defined in the right-hand sense for traveling
         from vertex 0 to 1 to 2.

The three output parameters are this facet's contribution to the 
volume, to the first-order moment (volume * COM displacement), and to
the inertia tensor.  Of course it is not the 2-D facet itself which
contributes, but the 3-D region which the facet "sweeps out" as it is
projected parallel to the (body-fixed) z-axis onto the xy plane.

facmom permutes the three components of the four input vectors so that
all three components of the first-order moment and all six independent
elements of the inertia tensor can be computed.  (The other three
elements are then filled in by symmetry.)

Since the input coordinates are in km, the output is in km^3, km^4,
and km^5 for volume, first-order moment, and inertia tensor elements,
respectively.  To convert inertia tensor elements to kg-m^2, multiply
by (10^18 * assumed density in g/cm^3).

I've verified that the (sometimes ugly) expressions in the code below
are exact integrals rather than approximations.

Commented 2003 May 12 by CM
***************************************************************************/
extern "C" {
#include "../shape/head.h"
}

#define X 0
#define Y 1
#define Z 2

#define PERMUTE( a) tmp=a[Z];a[Z]=a[Y];a[Y]=a[X];a[X]=tmp;


__device__ void dev_facmom( double fv0[3], double fv1[3], double fv2[3], double fn[3],
             double *dv, double dvr[3], double dI[3][3])
{
  double jac, v0[3], v1[3], v2[3], n[3], r, s, t, tmp;
  int i, sign;
  
  for (i=0; i<=2; i++) {
    v0[i] = fv0[i];
    v1[i] = fv1[i];
    v2[i] = fv2[i];
    n[i] = fn[i];
  }

  /* Flip the sign of moment contributions from facets which face away from us;
   * note that this does NOT necessarily make these contributions negative */
  if (n[Z] < 0.0)   sign = -1;
  else			    sign = 1;

  /*  jac is double the area of the facet's projection onto the xy plane  */
  jac = fabs((v1[X] - v0[X])*(v2[Y] - v1[Y])- (v2[X] - v1[X])*(v1[Y] - v0[Y]));

  /* dv is the volume of the region between the facet and the xy plane
   * (projecting parallel to the z-axis)                                */
  *dv = sign*jac*(v0[Z] + v1[Z] + v2[Z])/6.0;

  for (i=0; i<=2; i++) {
    /* Evaluate three integrals over volume dv. r is the integral of z (i.e.,
     * the z-component of the first-order moment)  */
    r = sign*jac*( v0[Z]*(v0[Z] + v1[Z] + v2[Z]) +
                   v1[Z]*(v1[Z] + v2[Z]) + v2[Z]*v2[Z] )/24.0;

    /* s is the integral of -y*z (off-diagonal inertia tensor element I_yz)  */
    s = -(sign*jac/120.0) 
        * ( (v0[Y] + v1[Y] + v2[Y])*(v0[Z]*v0[Z] + v1[Z]*v1[Z] + v2[Z]*v2[Z])
             + 2.0*(v0[Y]*v0[Z] + v1[Y]*v1[Z] + v2[Y]*v2[Z])*(v0[Z] + v1[Z] + v2[Z])
             + (v0[Y]*v1[Z]*v2[Z] + v1[Y]*v0[Z]*v2[Z] + v2[Y]*v0[Z]*v1[Z]) );

    /* t is integral of z^2. Contributes to diagonal inertia tensor elements I_xx and I_yy */
    t = (sign*jac/60.0)
        * ( (v0[Z]*v0[Z] + v1[Z]*v1[Z] + v2[Z]*v2[Z])*(v0[Z] + v1[Z] + v2[Z]) 
            + v0[Z]*v1[Z]*v2[Z] );

    /*  Assign values, then permute coordinates (x,y,z) so that
        we can build up all vector components and tensor elements  */

    switch (i) {
    case 0:
        dvr[Z] = r;
        dI[Y][Z] = s;
        dI[X][X] = t;
        dI[Y][Y] = t;
        break;
    case 1:
        dvr[Y] = r;
        dI[X][Y] = s;
        dI[X][X] += t;
        dI[Z][Z] = t;
        break;
    case 2:
        dvr[X] = r;
        dI[X][Z] = s;
        dI[Y][Y] += t;
        dI[Z][Z] += t;
    }
    if (i < 2) {
      PERMUTE( v0);
      PERMUTE( v1);
      PERMUTE( v2);
      PERMUTE( n);
      jac = fabs((v1[X] - v0[X])*(v2[Y] - v1[Y])
                               - (v2[X] - v1[X])*(v1[Y] - v0[Y]) );
      if (n[Z] < 0.0)
        sign = -1;
      else
        sign = 1;
    }
  }

  /*  Use symmetry to complete the inertia tensor  */

  dI[Y][X] = dI[X][Y];
  dI[Z][X] = dI[X][Z];
  dI[Z][Y] = dI[Y][Z];
}
__device__ void dev_facmom_f( float3 fv0, float3 fv1, float3 fv2, float3 fn,
             float *dv, float3 *dvr, float3 dI[3])
{
  float jac, r, s, t;
  float3 v0, v1, v2, tmp, n;
  int i, sign;

  v0.x = fv0.x;
  v1.x = fv1.x;
  v2.x = fv2.x;
  n.x = fn.x;
  v0.y = fv0.y;
  v1.y = fv1.y;
  v2.y = fv2.y;
  n.y = fn.y;
  v0.z = fv0.z;
  v1.z = fv1.z;
  v2.z = fv2.z;
  n.z = fn.z;

  /* Flip the sign of moment contributions from facets which face away from us;
   * note that this does NOT necessarily make these contributions negative */
  if (n.z < 0.0)    sign = -1;
  else			    sign = 1;

  /*  jac is double the area of the facet's projection onto the xy plane  */
  jac = fabs((v1.x - v0.x)*(v2.y - v1.y) - (v2.x - v1.x)*(v1.y - v0.y));

  /* dv is the volume of the region between the facet and the xy plane
   * (projecting parallel to the z-axis)                                */
  *dv = sign*jac*(v0.z + v1.z + v2.z)/6.0;

  for (i=0; i<=2; i++) {
    /* Evaluate three integrals over volume dv. r is the integral of z (i.e.,
     * the z-component of the first-order moment)  */
    r = sign*jac*( v0.z*(v0.z + v1.z + v2.z) +
                   v1.z*(v1.z + v2.z) + v2.z*v2.z )/24.0;

    /* s is the integral of -y*z (off-diagonal inertia tensor element I_yz)  */
    s = -(sign*jac/120.0)
        * ( (v0.y + v1.y + v2.y)*(v0.z*v0.z + v1.z*v1.z + v2.z*v2.z)
             + 2.0*(v0.y*v0.z + v1.y*v1.z + v2.y*v2.z)*(v0.z + v1.z + v2.z)
             + (v0.y*v1.z*v2.z + v1.y*v0.z*v2.z + v2.y*v0.z*v1.z) );

    /* t is integral of z^2. Contributes to diagonal inertia tensor elements I_xx and I_yy */
    t = (sign*jac/60.0) * ( (v0.z*v0.z + v1.z*v1.z + v2.z*v2.z)*(v0.z + v1.z + v2.z)
            + v0.z*v1.z*v2.z );

    /*  Assign values, then permute coordinates (x,y,z) so that
        we can build up all vector components and tensor elements  */

    switch (i) {
    case 0:
        dvr->z = r;
        dI[1].z = s;
        dI[0].x = t;
        dI[1].y = t;
        break;
    case 1:
        dvr->y = r;
        dI[1].y = s;
        dI[1].x += t;
        dI[2].z = t;
        break;
    case 2:
        dvr->x = r;
        dI[0].z = s;
        dI[1].y += t;
        dI[2].z += t;
    }
    if (i < 2) {
      tmp = v0;
      v0.x = tmp.z;
      v0.y = tmp.x;
      v0.z = tmp.y;
      tmp = v1;
      v1.x = tmp.z;
      v1.y = tmp.x;
      v1.z = tmp.y;
      tmp = v2;
      v2.x = tmp.z;
      v2.y = tmp.x;
      v2.z = tmp.y;
      tmp = n;
      n.x = tmp.z;
      n.y = tmp.x;
      n.z = tmp.y;

      jac = fabs((v1.x - v0.x)*(v2.y - v1.y) - (v2.x - v1.x)*(v1.y - v0.y) );
      if (n.z < 0.0)
        sign = -1;
      else
        sign = 1;
    }
  }

  /*  Use symmetry to complete the inertia tensor  */

  dI[1].x = dI[0].y;
  dI[2].x = dI[0].z;
  dI[2].y = dI[1].z;
}

#undef X
#undef Y
#undef Z
#undef PERMUTE
