/*****************************************************************************************
                                                                            realize_mod.c

Takes a struct mod_t model and "realizes" its components as polyhedral solids made up of
triangular facets.

Modified 2016 July 9 by Matthias Engels:
	Adapted for use with shape-cuda.
------------------------------------------------------------------------------------------
Modified 2014 April 26 by CM:
    Increase the minimum permitted value of the highest-order coefficient in the cubic
        equation that locates an ovoid vertex: if the coefficient is smaller than this
        minimum, treat it as if it's zero and solve a quadratic equation instead

Modified 2014 March 22 by CM:
    Relax the tolerance for finding a valid ovoid vertex position

Modified 2014 March 10 by CM:
    Guard against roundoff problems when computing vertex positions for ovoid components
        with very small |k|

Modified 2014 February 10 by CM:
    Implement multiple radar and optical scattering laws

Modified 2013 August 28 by CM:
    Set the bad diameter flag for harmonic components with tiny or negative vertex
        displacements, and for harmonic and vertex components with tiny or negative
        "scale factor" values

Modified 2013 June 2 by CM:
    In the cubic_realroot routine, initialize nrealroots to avoid compilation warning
    Fix a comment

Modified 2013 May 20 by CM:
    Implement ovoid shape components

Modified 2012 July 4 by CM:
    Add test in "realize_coordinates" routine to avoid compilation warning

Modified 2011 September 2 by CM:
    Bug fix: the "check_surface" routine makes use of facet normals when identifying
        active vs. inactive vertices and facets, but facet normals weren't being computed
        until *after* check_surface was called
    Make the code more modular (and address the above bug) by introducing the
        "realize_coordinates" and "compute_moments" routines, as per the version of
        realize_mod in the SHERMAN package
    Store the area and the centroid coordinates of each facet
    Add "harmlambert" optical scattering law (compute facet angular coordinates)

Modified 2010 September 1 by CM:
    Add "facetnorm" argument to the rayfacint routine

Modified 2010 June 1 by CM:
    Change "scalefactor" parameter from a scalar to a 3-component vector

Modified 2010 March 19 by CM:
    Implement '=' state for vertex deviations

Modified 2009 November 15 by CM:
    In the "check_surface" routine, eliminate an unused variable and fix
        a couple of ambiguous nested if-then-else statements

Modified 2009 August 3 by CM:
    For the "harmlommel" "harmhapke" "harmkaas" and "harmcosine_diff"
        inhomogeneous scattering laws, compute the spherical coordinates
        (theta and phi) of each facet after each component's rotational and
        translational offsets have been applied rather than before, so that
        these laws can be used for multiple-component models
    For multiple-component models, use a more careful method (already used
        for facets) to determine which vertices are on the model's surface;
        also, for both vertices and facets, allow for a bit of roundoff
        error in this determination by adding a tolerance argument to the
        "rayfacint" routine
    For multiple-component models, determine the new "act" (active) flag
        for each model side
    For multiple-component models, fix a bug in computing the center of mass
        for individual components

Modified 2009 July 5 by CM:
    Turn each component's rotational offsets into a rotation matrix here
        rather than in the "read_mod" routine, in case the offsets are
        being allowed to float

Modified 2009 July 1 by CM:
    Add "check_surface" routine that determines which facets of a
        multiple-component model lie on the model's surface rather than
        interior to the model
    For multiple-component models, when computing the area and the moments
        of the overall model, ignore facets that lie interior to the model

Modified 2009 April 3 by CM:
    Fix slight bug in defining function a[i] = 1/radius^2 when a/b or b/c
        is tiny or negative for ellipsoid components
    Initialize the "baddiam_logfactor" parameter and set its value when
        2a, a/b, or b/c is tiny or negative for ellipsoid components

Modified 2007 August 10 by CM:
    Eliminate unused variable

Modified 2007 January 8 by CM:
    Define "scalefactor" state for vertex realizations of ellipsoid and
        harmonic components, not just its value

Modified 2006 October 1 by CM:
    Add "scalefactor" to harmonic and vertex shape structures
    Replace ellipsoid diameters D with two_a, a_over_b, b_over_c

Modified 2005 September 6 by CM:
    Add computation of facet angular coordinates for use with harmonic
        scattering laws

Modified 2005 August 17 by CM:
    Move computation of spherical harmonic functions afactor and bfactor
        from here to read_mod.c, so that it can be done just once per fit

Modified 2005 February 28 by CM:
    Initialize the "baddiam" parameter (flag indicating tiny or negative
        ellipsoid diameters) to 0 here rather than in bestfit.c so that
        it can be used for actions other than "fit"

Modified 2004 August 23 by CM:
    Eliminated newtheta and oldcostheta variables and THETATOL constant,
        since they weren't actually being used (i.e., the test in which
        they were included was always true)

Modified 2003 April 17 by CM:
    Added computation of component and model moments; this used to
        be done in function penalties (but wasn't always being done)
    Added code to cope with tiny or negative ellipsoid diameters;
        as a result, must now pass the model's parameter structure
        as an argument to realize_mod
    Added surface area computation for components and for the full model
 *****************************************************************************************/
extern "C" {
#include "../shape/head.h"
}
#define HAIRWIDTH 1.0e-7
#define SMALLRATIO 0.01
#define SMALLOVOIDK1 0.01
#define SMALLOVOIDK2 1.0e-6
#define OVOIDTOL 1.0e-6
#define MAXEDGE 100
#define EDGETOL 1.0e-14
#define RTOL 1000*EDGETOL
#define SMALLCOEFF3 1.0e-5

/* These 2 device variables are to get nf and nv from the GPU-located dmod file */
__device__ int dnv, dnf, dns;
__device__ double d_a[3];
__device__ float3 d_af;
__device__ float af_radius, a_over_b_f, b_over_c_f, k_asym_f, x0term_f, numer_f, denom_f, x0_f;
__device__ double a_radius, a_over_b, b_over_c, k_asym, x0term, numer, denom, x0;
__device__ int harmonic_scatlaw;
__device__ float rm_area=0.0, rm_ifarea=0.0, rm_vol=0.0, rm_ifvol=0.0,
		rm_dcom[3], rm_ifdcom[3], rm_dI[3][3], rm_ifdI[3][3];
static int nv, nf, ns;
static dim3 nvBLK,nvTHD,nfBLK,nfTHD,nsBLK,nsTHD;
__host__ void realize_coordinates_cuda(struct par_t *dpar, struct mod_t *dmod, unsigned char type);
__host__ void realize_coordinates_f_cuda( struct par_t *dpar, struct mod_t *dmod, unsigned char type);
__host__ void check_surface_cuda(struct mod_t *dmod);
__host__ void compute_moments_cuda(struct mod_t *dmod);

__global__ void set_diam_krnl(struct par_t *dpar, struct mod_t *dmod){
	/* This is a single-thread kernel */
	if (threadIdx.x == 0) {
		dpar->baddiam = 0;
		dpar->baddiam_logfactor = 0;
		dnv = dmod->shape.comp[0].real.nv;
		dnf = dmod->shape.comp[0].real.nf;
		dns = dmod->shape.comp[0].real.ns;
	}
	__syncthreads();
}
__global__ void ellipse_diameter_krnl(struct par_t *dpar, struct mod_t *dmod) {
	/* This is a single-thread kernel */
	float diam, diamratio;

	if (threadIdx.x == 0) {
		diam = (float)dmod->shape.comp[0].desc.ell.two_a.val;
		if (diam > HAIRWIDTH) {
			d_a[0] = 2.0/diam; /* 1/radii */
		} else {
			d_a[0] = (2.0/HAIRWIDTH) * (1 + HAIRWIDTH - diam);
			dpar->baddiam = 1;
			dpar->baddiam_logfactor += log(1 + HAIRWIDTH - diam);
		}
		diam = (2.0/d_a[0]);
		diamratio = dmod->shape.comp[0].desc.ell.a_over_b.val;
		if (diamratio > SMALLRATIO) {
			d_a[1] = 2.0/(diam/diamratio);
		} else {
			d_a[1] = (2.0/(diam/SMALLRATIO)) / (1 + SMALLRATIO - diamratio);
			dpar->baddiam = 1;
			dpar->baddiam_logfactor += log(1 + SMALLRATIO - diamratio);
		}
		diam = (2.0/d_a[1]);
		diamratio = dmod->shape.comp[0].desc.ell.b_over_c.val;
		if (diamratio > SMALLRATIO) {
			d_a[2] = 2.0/(diam/diamratio);
		} else {
			d_a[2] = (2.0/(diam/SMALLRATIO)) / (1 + SMALLRATIO - diamratio);
			dpar->baddiam = 1;
			dpar->baddiam_logfactor += log(1 + SMALLRATIO - diamratio);
		}
		d_a[0] *= d_a[0];
		d_a[1] *= d_a[1];
		d_a[2] *= d_a[2];
	}
}
__global__ void ellipse_diameter_f_krnl(struct par_t *dpar, struct mod_t *dmod) {
	/* This is a single-thread kernel */
	float diam, diamratio;

	if (threadIdx.x == 0) {
		diam = dmod->shape.comp[0].desc.ell.two_a.val;
		if (diam > HAIRWIDTH) {
			d_af.x = 2.0/diam; /* 1/radii */
		} else {
			d_af.x = (2.0/HAIRWIDTH) * (1 + HAIRWIDTH - diam);
			dpar->baddiam = 1;
			dpar->baddiam_logfactor += log(1 + HAIRWIDTH - diam);
		}
		diam = (2.0/d_af.x);
		diamratio = (float)dmod->shape.comp[0].desc.ell.a_over_b.val;
		if (diamratio > SMALLRATIO) {
			d_af.y = 2.0/(diam/diamratio);
		} else {
			d_af.y = (2.0/(diam/SMALLRATIO)) / (1 + SMALLRATIO - diamratio);
			dpar->baddiam = 1;
			dpar->baddiam_logfactor += log(1 + SMALLRATIO - diamratio);
		}
		diam = (2.0/d_af.y);
		diamratio = (float)dmod->shape.comp[0].desc.ell.b_over_c.val;
		if (diamratio > SMALLRATIO) {
			d_af.z = 2.0/(diam/diamratio);
		} else {
			d_af.z = (2.0/(diam/SMALLRATIO)) / (1 + SMALLRATIO - diamratio);
			dpar->baddiam = 1;
			dpar->baddiam_logfactor += log(1 + SMALLRATIO - diamratio);
		}
		d_af.x *= d_af.x;
		d_af.y *= d_af.y;
		d_af.z *= d_af.z;
	}
}
__global__ void ellipse_distance_krnl(struct par_t *dpar, struct mod_t *dmod)
{
	int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int j;
	double den;

	if (offset < dmod->shape.comp[0].real.nv) {

		/* Routine setuprealver (called by setupreal, which was called by
		 * read_mod) already created as many ellipsoid vertices as were needed
		 * for specified value of theta_steps, and initialized direction
		 * cosines u[j] for each vertex to be
		 * 		sin(theta)cos(phi), sin(theta)sin(phi), and cos(theta) for
		 * 		j=0, 1, and 2, respectively.
		 *
		 * These values are x/r, y/r, and z/r, where r is distance from origin
		 * to ellipsoid surface along direction (theta, phi) for given vertex.
		 * Since an ellipsoid has (x/a)^2 + (y/b)^2 + (z/c)^2 = 1, quantity
		 * "den" in code below is equal to 1/(r^2) for vertex i.
		 *
		 * Note that setuprealver initialized all vertex "base points" a[j] to
		 * be zero for ellipsoid components; hence "deviation" r is in fact the
		 * entire thing.		 */
		den = 0.0;
		for (j=0; j<=2; j++)
			den += d_a[j]*( dmod->shape.comp[0].real.v[offset].u[j]
			                                                     * dmod->shape.comp[0].real.v[offset].u[j] );
		dmod->shape.comp[0].real.v[offset].r.val = 1/sqrt(den);
	}
}
__global__ void ellipse_scalefactor_krnl(struct mod_t *dmod) {
	/* Single-threaded kernel */
	int j;
	if (threadIdx.x == 0) {
		dmod->shape.comp[0].real.scalefactor[0].state = dmod->shape.comp[0].desc.ell.two_a.state;
		dmod->shape.comp[0].real.scalefactor[1].state = dmod->shape.comp[0].desc.ell.a_over_b.state;
		dmod->shape.comp[0].real.scalefactor[2].state = dmod->shape.comp[0].desc.ell.b_over_c.state;
		for (j=0; j<=2; j++)
			dmod->shape.comp[0].real.scalefactor[j].val = 1.0;
	}
}
__global__ void set_ovoid_parameters_krnl(struct par_t *dpar, struct mod_t *dmod) {
	//, double a_radius, double a_over_b, double b_over_c, double
	//  k_asym, double x0term, double numer, double denom, double x0) {

	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		/*  Determine all shape parameters, making sure that none are out of bounds  */
		a_radius = dmod->shape.comp[0].desc.ovoid.two_a.val / 2;
		if (a_radius <= HAIRWIDTH/2) {
			dpar->baddiam = 1;
			dpar->baddiam_logfactor += log(1 + HAIRWIDTH - 2*a_radius);
			a_radius = (HAIRWIDTH/2) / (1 + HAIRWIDTH - 2*a_radius);
		}
		a_over_b = dmod->shape.comp[0].desc.ovoid.a_over_b.val;
		if (a_over_b <= SMALLRATIO) {
			dpar->baddiam = 1;
			dpar->baddiam_logfactor += log(1 + SMALLRATIO - a_over_b);
			a_over_b = SMALLRATIO / (1 + SMALLRATIO - a_over_b);
		}
		b_over_c = dmod->shape.comp[0].desc.ovoid.b_over_c.val;
		if (b_over_c <= SMALLRATIO) {
			dpar->baddiam = 1;
			dpar->baddiam_logfactor += log(1 + SMALLRATIO - b_over_c);
			b_over_c = SMALLRATIO / (1 + SMALLRATIO - b_over_c);
		}
		k_asym = dmod->shape.comp[0].desc.ovoid.k.val;
		if (fabs(k_asym) > 1 - SMALLVAL) {
			dpar->baddiam = 1;
			dpar->baddiam_logfactor += log(fabs(k_asym) + SMALLVAL);
			if (k_asym > 0.0)
				k_asym = 1 - SMALLVAL*(1 - SMALLVAL)/k_asym;
			else
				k_asym = -1 - SMALLVAL*(1 - SMALLVAL)/k_asym;
		}

		/* Compute x0, the x-offset that places the ovoid's center of mass at the
		 * origin; for small |k|, use an analytical approximation to avoid
		 * roundoff problems       */

		if (fabs(k_asym) > SMALLOVOIDK1) {
			x0term = 3*(1 - k_asym*k_asym)*log((1 + k_asym)/(1 - k_asym));
			numer = 2.0*k_asym*(3 - 2*k_asym*k_asym) - x0term;
			denom = 2.0*k_asym*(3 -   k_asym*k_asym) - x0term;
			x0 = (a_radius/k_asym)*(numer/denom);
		} else {
			x0 = 0.4*k_asym*a_radius;
		}
	}
}
__global__ void set_ovoid_parameters_f_krnl(struct par_t *dpar, struct mod_t *dmod) {
	//, double a_radius, double a_over_b, double b_over_c, double
	//  k_asym, double x0term, double numer, double denom, double x0) {

	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		/*  Determine all shape parameters, making sure that none are out of bounds  */
		af_radius = (float)(dmod->shape.comp[0].desc.ovoid.two_a.val / 2);
		if (af_radius <= HAIRWIDTH/2) {
			dpar->baddiam = 1;
			dpar->baddiam_logfactor += log(1 + HAIRWIDTH - 2*af_radius);
			af_radius = (HAIRWIDTH/2) / (1 + HAIRWIDTH - 2*af_radius);
		}
		a_over_b_f = (float)dmod->shape.comp[0].desc.ovoid.a_over_b.val;
		if (a_over_b_f <= SMALLRATIO) {
			dpar->baddiam = 1;
			dpar->baddiam_logfactor += log(1 + SMALLRATIO - a_over_b_f);
			a_over_b_f = SMALLRATIO / (1 + SMALLRATIO - a_over_b_f);
		}
		b_over_c_f = (float)dmod->shape.comp[0].desc.ovoid.b_over_c.val;
		if (b_over_c_f <= SMALLRATIO) {
			dpar->baddiam = 1;
			dpar->baddiam_logfactor += log(1 + SMALLRATIO - b_over_c_f);
			b_over_c_f = SMALLRATIO / (1 + SMALLRATIO - b_over_c_f);
		}
		k_asym_f = (float)dmod->shape.comp[0].desc.ovoid.k.val;
		if (fabs(k_asym_f) > 1 - SMALLVAL) {
			dpar->baddiam = 1;
			dpar->baddiam_logfactor += log(fabs(k_asym_f) + SMALLVAL);
			if (k_asym_f > 0.0)
				k_asym_f = 1 - SMALLVAL*(1 - SMALLVAL)/k_asym_f;
			else
				k_asym_f = -1 - SMALLVAL*(1 - SMALLVAL)/k_asym_f;
		}

		/* Compute x0, the x-offset that places the ovoid's center of mass at the
		 * origin; for small |k|, use an analytical approximation to avoid
		 * roundoff problems       */

		if (fabs(k_asym_f) > SMALLOVOIDK1) {
			x0term_f = 3*(1 - k_asym_f*k_asym_f)*log((1 + k_asym_f)/(1 - k_asym_f));
			numer_f = 2.0*k_asym_f * (3 - 2*k_asym_f*k_asym_f) - x0term_f;
			denom_f = 2.0*k_asym_f * (3 - k_asym_f*k_asym_f) - x0term_f;
			x0_f = (af_radius/k_asym_f) * (numer_f/denom_f);
		} else {
			x0_f = 0.4 * k_asym_f * af_radius;
		}
	}
}
__global__ void ovoid_distance_krnl(struct par_t *dpar, struct mod_t *dmod)
//double d_a[3], double a_radius, double a_over_b, double b_over_c, double
//k_asym, double x0term, double numer, double denom, double x0)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j, k, nrealroots;
	double a_over_c, h, alpha0, u_x, coeff[4], goodroot, realroot[3], x_over_a;

	if (i < dmod->shape.comp[0].real.nv) {

		a_over_c = a_over_b*b_over_c;
		h = a_over_b*a_over_b*dmod->shape.comp[0].real.v[i].u[1]
		                                                      *dmod->shape.comp[0].real.v[i].u[1] + a_over_c*a_over_c
		                                                      *dmod->shape.comp[0].real.v[i].u[2]*dmod->shape.comp[0].real.v[i].u[2];

		alpha0 = x0/a_radius;
		u_x = dmod->shape.comp[0].real.v[i].u[0];
		coeff[3] = (h - u_x*u_x)*k_asym*u_x;
		coeff[2] = (1 + 3*k_asym*alpha0)*u_x*u_x + h*(1 - k_asym*alpha0);
		coeff[1] = (k_asym - (2 + 3*k_asym*alpha0)*alpha0)*u_x;
		coeff[0] = -(1 - alpha0*alpha0)*(1 + k_asym*alpha0);
		if (fabs(k_asym) <= SMALLOVOIDK2) {

			/* |k| is very small, so guard against roundoff error by
			 * computing the vertex position for an ellipsoid (k = 0) and then
			 * applying a first-order correction for nonzero k  */
			goodroot = 1/sqrt(u_x*u_x + h);
			goodroot -= (coeff[3]*goodroot*goodroot*goodroot + coeff[1]*goodroot)
    			    								/ (3*coeff[3]*goodroot*goodroot + 2*coeff[2]*goodroot + coeff[1]);
		} else {

			/* |k| isn't very small, so solve the cubic equation  */
			nrealroots = cubic_realroots_cuda( coeff, realroot);
			goodroot = -HUGENUMBER;
			for (k=0; k<nrealroots; k++)
				if (realroot[k] >= 0.0) {
					x_over_a = realroot[k]*u_x;
					if (fabs(x_over_a - alpha0) - 1 < OVOIDTOL)
						goodroot = MAX( goodroot, realroot[k]);
				}
		}
		if (goodroot < 0.0)
			printf("Can't compute vertex displacement for ovoid vertex %d\n", i);

		dmod->shape.comp[0].real.v[i].r.val = goodroot*a_radius;

		/* Assign scalefactor values  */
		dmod->shape.comp[0].real.scalefactor[0].state = dmod->shape.comp[0].desc.ovoid.two_a.state;
		dmod->shape.comp[0].real.scalefactor[1].state = dmod->shape.comp[0].desc.ovoid.a_over_b.state;
		dmod->shape.comp[0].real.scalefactor[2].state = dmod->shape.comp[0].desc.ovoid.b_over_c.state;
		for (j=0; j<=2; j++)
			dmod->shape.comp[0].real.scalefactor[j].val = 1.0;
	}
}
__global__ void ovoid_distance_f_krnl(struct par_t *dpar, struct mod_t *dmod)
//double d_a[3], double a_radius, double a_over_b, double b_over_c, double
//k_asym, double x0term, double numer, double denom, double x0)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j, nrealroots;
	float a_over_c, h, alpha0, u_x, goodroot, x_over_a;
	float4 coeff;
	float3 realroot;

	if (i < dnv) {

		a_over_c = a_over_b_f * b_over_c_f;
		h = a_over_b_f * a_over_b_f * dmod->shape.comp[0].real.v[i].u[1]
		                                                              *dmod->shape.comp[0].real.v[i].u[1] + a_over_c *a_over_c
		                                                              *dmod->shape.comp[0].real.v[i].u[2] * dmod->shape.comp[0].real.v[i].u[2];

		alpha0 = x0_f/af_radius;
		u_x = (float)dmod->shape.comp[0].real.v[i].u[0];

		coeff.z = (h - u_x*u_x)*k_asym_f*u_x;
		coeff.y = (1 + 3*k_asym_f * alpha0) * u_x*u_x + h*(1 - k_asym_f*alpha0);
		coeff.x = (k_asym_f - (2 + 3*k_asym_f * alpha0) * alpha0) * u_x;
		coeff.w = -(1 - alpha0*alpha0)*(1 + k_asym_f*alpha0);
		if (fabs(k_asym_f) <= SMALLOVOIDK2) {

			/* |k| is very small, so guard against roundoff error by
			 * computing the vertex position for an ellipsoid (k = 0) and then
			 * applying a first-order correction for nonzero k  */
			goodroot = 1/sqrt(u_x*u_x + h);
			goodroot -= (coeff.z*goodroot*goodroot*goodroot + coeff.x*goodroot)
    			   			/ (3*coeff.z*goodroot*goodroot + 2*coeff.y*goodroot + coeff.x);
		} else {

			/* |k| isn't very small, so solve the cubic equation  */
			nrealroots = cubic_realroots_f_cuda( coeff, &realroot);
			goodroot = -HUGENUMBER;

			if (realroot.x >= 0.0) {
				x_over_a = realroot.x * u_x;
				if (fabs(x_over_a - alpha0) - 1 < OVOIDTOL)
					goodroot = MAX(goodroot, realroot.x);
			}
			if (nrealroots > 1 && realroot.y >= 0.0) {
				x_over_a = realroot.y * u_x;
				if (fabs(x_over_a - alpha0) - 1 < OVOIDTOL)
					goodroot = MAX(goodroot, realroot.y);
			}
			if (nrealroots == 3 && realroot.z >= 0.0) {
				x_over_a = realroot.z * u_x;
				if (fabs(x_over_a - alpha0) - 1 < OVOIDTOL)
					goodroot = MAX(goodroot, realroot.z);
			}
		}
		if (goodroot < 0.0)
			printf("Can't compute vertex displacement for ovoid vertex %d\n", i);

		dmod->shape.comp[0].real.v[i].r.val = (double)goodroot*af_radius;

		/* Assign scalefactor values  */
		dmod->shape.comp[0].real.scalefactor[0].state = dmod->shape.comp[0].desc.ovoid.two_a.state;
		dmod->shape.comp[0].real.scalefactor[1].state = dmod->shape.comp[0].desc.ovoid.a_over_b.state;
		dmod->shape.comp[0].real.scalefactor[2].state = dmod->shape.comp[0].desc.ovoid.b_over_c.state;
		for (j=0; j<=2; j++)
			dmod->shape.comp[0].real.scalefactor[j].val = 1.0;
	}
}
__global__ void harmonic_krnl(struct par_t *dpar, struct mod_t *dmod)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int L, l, m;
	double r;

	if (i < dnv) {

		L = dmod->shape.comp[0].desc.har.nhar;
		r = 0.0;

		for (l=0; l<=L; l++) {
			r += dmod->shape.comp[0].desc.har.a[l][0].val
					* dmod->shape.comp[0].real.v[i].afactor[l][0];
			for (m=1; m<=l; m++)
				r += dmod->shape.comp[0].desc.har.a[l][m].val
				* dmod->shape.comp[0].real.v[i].afactor[l][m]
				                                           + dmod->shape.comp[0].desc.har.b[l][m].val
				                                           * dmod->shape.comp[0].real.v[i].bfactor[l][m];
		}
		if (r > HAIRWIDTH/2) {
			dmod->shape.comp[0].real.v[i].r.val = r;
		} else {
			dpar->baddiam = 1;
			dpar->baddiam_logfactor += log(1 + HAIRWIDTH - 2*r) / ((L+1)*(L+1));
			dmod->shape.comp[0].real.v[i].r.val = (HAIRWIDTH/2) / (1 + HAIRWIDTH - 2*r);
		}
	}
}
__global__ void harmonic_f_krnl(struct par_t *dpar, struct mod_t *dmod)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int L, l, m;
	float r;

	if (i < dnv) {

		L = dmod->shape.comp[0].desc.har.nhar;
		r = 0.0;

		for (l=0; l<=L; l++) {
			r += __double2float_rn(dmod->shape.comp[0].desc.har.a[l][0].val)
							* __double2float_rn(dmod->shape.comp[0].real.v[i].afactor[l][0]);
			for (m=1; m<=l; m++)
				r += dmod->shape.comp[0].desc.har.a[l][m].val
				* dmod->shape.comp[0].real.v[i].afactor[l][m]
				                                           + dmod->shape.comp[0].desc.har.b[l][m].val
				                                           * dmod->shape.comp[0].real.v[i].bfactor[l][m];
		}
		if (r > HAIRWIDTH/2) {
			dmod->shape.comp[0].real.v[i].r.val = r;
		} else {
			dpar->baddiam = 1;
			dpar->baddiam_logfactor += log(1 + HAIRWIDTH - 2*r) / ((L+1)*(L+1));
			dmod->shape.comp[0].real.v[i].r.val = (HAIRWIDTH/2) / (1 + HAIRWIDTH - 2*r);
		}
	}
}
__global__ void harmonic_scalefactor_krnl(struct par_t *dpar, struct mod_t *dmod)
{
	// This is a 3-thread single thread kernel
	int j = threadIdx.x;

	if (j < 3){
		if (j > 0 && dmod->shape.comp[0].desc.har.scalefactor[j].state == '=')
			dmod->shape.comp[0].desc.har.scalefactor[j].val
			= dmod->shape.comp[0].desc.har.scalefactor[j-1].val;
		dmod->shape.comp[0].real.scalefactor[j].state = dmod->shape.comp[0].desc.har.scalefactor[j].state;
		dmod->shape.comp[0].real.scalefactor[j].val = dmod->shape.comp[0].desc.har.scalefactor[j].val;
		if (dmod->shape.comp[0].real.scalefactor[j].val <= SMALLRATIO) {
			dpar->baddiam = 1;
			dpar->baddiam_logfactor += log(1 + SMALLRATIO - dmod->shape.comp[0].real.scalefactor[j].val);
			dmod->shape.comp[0].real.scalefactor[j].val = SMALLRATIO
					/ (1 + SMALLRATIO - dmod->shape.comp[0].real.scalefactor[j].val);
		}
	}
}
__global__ void vertex_update_dev_krnl(struct par_t *dpar, struct mod_t *dmod)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int v_mirror;
	if (i < dmod->shape.comp[0].real.nv) {
		if (dmod->shape.comp[0].real.v[i].r.state == '=') {
			v_mirror = dmod->shape.comp[0].real.v[i].v_mirror;
			dmod->shape.comp[0].real.v[i].r.val =
					dmod->shape.comp[0].real.v[v_mirror].r.val;
		}
	}
}
__global__ void vertex_scalefactor_krnl(struct par_t *dpar, struct mod_t *dmod)
{
	// This is a 3-thread single thread kernel
	//int j = threadIdx.x;

	if (threadIdx.x == 0) {
		for (int j=0; j<=2; j++) {
			if (j > 0 && dmod->shape.comp[0].desc.ver.scalefactor[j].state == '=')
				dmod->shape.comp[0].desc.ver.scalefactor[j].val
				= dmod->shape.comp[0].desc.ver.scalefactor[j-1].val;

			dmod->shape.comp[0].real.scalefactor[j].val = dmod->shape.comp[0].desc.ver.scalefactor[j].val;

			if (dmod->shape.comp[0].real.scalefactor[j].val <= SMALLRATIO) {
				dpar->baddiam = 1;
				dpar->baddiam_logfactor += log(1 + SMALLRATIO - dmod->shape.comp[0].real.scalefactor[j].val);
				dmod->shape.comp[0].real.scalefactor[j].val = SMALLRATIO
						/ (1 + SMALLRATIO - dmod->shape.comp[0].real.scalefactor[j].val);
			}
		}
	}
	__syncthreads();
}
__global__ void calc_vertex_co_krnl(struct par_t *dpar, struct mod_t *dmod)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j;
	if (i < dnv){
		for (j=0; j<=2; j++)
			dmod->shape.comp[0].real.v[i].x[j] = dmod->shape.comp[0].real.scalefactor[j].val
			* (dmod->shape.comp[0].real.v[i].u[j] * dmod->shape.comp[0].real.v[i].r.val
					+ dmod->shape.comp[0].real.v[i].a[j]);
	}
}
__global__ void perform_rotation_krnl_old(struct par_t *dpar, struct mod_t *dmod)
{
	/* Single-threaded kernel */
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	double x[3];
	int j, k;

	if (threadIdx.x == 0) {

		if (!(dmod->shape.comp[0].rot[0].val == 0 &&
				dmod->shape.comp[0].rot[1].val == 0 &&
				dmod->shape.comp[0].rot[2].val == 0    )) {
			if (i <dmod->shape.comp[0].real.nv){

				for (j=0; j<=2; j++) {
					x[j] = 0.0;
					for (k=0; k<=2; k++)
						x[j] += dmod->shape.comp[0].m[j][k] * dmod->shape.comp[0].real.v[i].x[k];
				}
				for (j=0; j<=2; j++)
					dmod->shape.comp[0].real.v[i].x[j] = x[j];
			}
		}
	}
}
__global__ void perform_rotation_krnl(struct par_t *dpar, struct mod_t *dmod)
{
	/* nv-threaded kernel */
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	double x[3];
	int j, k, test;

	if (i == 0) {

		if (!(dmod->shape.comp[0].rot[0].val == 0 &&
				dmod->shape.comp[0].rot[1].val == 0 &&
				dmod->shape.comp[0].rot[2].val == 0    ))
			test = 1;
	}
	__syncthreads();

	if (i < dnv && test == 1) {

		for (j=0; j<=2; j++) {
			x[j] = 0.0;
			for (k=0; k<=2; k++)
				x[j] += dmod->shape.comp[0].m[j][k] * dmod->shape.comp[0].real.v[i].x[k];
		}
		for (j=0; j<=2; j++)
			dmod->shape.comp[0].real.v[i].x[j] = x[j];
	}
}
__global__ void perform_rotation_f_krnl(struct par_t *dpar, struct mod_t *dmod)
{
	/* nv-threaded kernel */
	int test, i = blockIdx.x * blockDim.x + threadIdx.x;
	float3 x, vx;

	__shared__ float m[3][3];

	if (i == 0) {
		if (!(dmod->shape.comp[0].rot[0].val == 0 &&
				dmod->shape.comp[0].rot[1].val == 0 &&
				dmod->shape.comp[0].rot[2].val == 0    ))
			test = 1;
	}
	__syncthreads();

	/* Now assign the shared m array for each warp */
	if (threadIdx.x == 0 && test == 1) {
		m[0][0] = __double2float_rn(dmod->shape.comp[0].m[0][0]);
		m[0][1] = __double2float_rn(dmod->shape.comp[0].m[0][1]);
		m[0][2] = __double2float_rn(dmod->shape.comp[0].m[0][2]);
		m[1][0] = __double2float_rn(dmod->shape.comp[0].m[1][0]);
		m[1][1] = __double2float_rn(dmod->shape.comp[0].m[1][1]);
		m[1][2] = __double2float_rn(dmod->shape.comp[0].m[1][2]);
		m[2][0] = __double2float_rn(dmod->shape.comp[0].m[2][0]);
		m[2][1] = __double2float_rn(dmod->shape.comp[0].m[2][1]);
		m[2][2] = __double2float_rn(dmod->shape.comp[0].m[2][2]);
	}
	__syncthreads();

	if (i < dnv && test == 1) {
		vx.x = (float)dmod->shape.comp[0].real.v[i].x[0];
		vx.y = (float)dmod->shape.comp[0].real.v[i].x[1];
		vx.z = (float)dmod->shape.comp[0].real.v[i].x[2];

		x.x = 0;
		x.x += m[0][0] * vx.x;
		x.x += m[0][1] * vx.y;
		x.x += m[0][2] * vx.z;
		x.y = 0;
		x.y += m[1][0] * vx.x;
		x.y += m[1][1] * vx.y;
		x.y += m[1][2] * vx.z;
		x.z = 0;
		x.z += m[2][0] * vx.x;
		x.z += m[2][1] * vx.y;
		x.z += m[2][2] * vx.z;

		dmod->shape.comp[0].real.v[i].x[0] = x.x;
		dmod->shape.comp[0].real.v[i].x[1] = x.y;
		dmod->shape.comp[0].real.v[i].x[2] = x.z;
	}
}
__global__ void perform_translation_old_krnl(struct par_t *dpar, struct mod_t *dmod)
{
	/* Single-threaded kernel */
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j;
	if (threadIdx.x == 0) {
		if (!(dmod->shape.comp[0].off[0].val == 0.0 &&
				dmod->shape.comp[0].off[1].val == 0.0 &&
				dmod->shape.comp[0].off[2].val == 0.0    )) {
			if (i <dmod->shape.comp[0].real.nv){

				for (j=0; j<=2; j++)
					dmod->shape.comp[0].real.v[i].x[j] += dmod->shape.comp[0].off[j].val;
			}
		}
	}
}
__global__ void perform_translation_krnl(struct par_t *dpar, struct mod_t *dmod)
{
	/* nv-threaded kernel */
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j, test;
	if (i == 0) {

		if (!(dmod->shape.comp[0].off[0].val == 0.0 &&
				dmod->shape.comp[0].off[1].val == 0.0 &&
				dmod->shape.comp[0].off[2].val == 0.0    ))
			test = 1;
	}
	__syncthreads();

	if (i < dnv && test == 1){
		for (j=0; j<=2; j++)
			dmod->shape.comp[0].real.v[i].x[j] += dmod->shape.comp[0].off[j].val;
	}
}
__global__ void set_optical_params_krnl(struct par_t *dpar, struct mod_t *dmod)
{
	/* Single-thread kernel */
	int ilaw;
	harmonic_scatlaw = 0;
	if (threadIdx.x == 0) {
		for (ilaw=0; ilaw<dmod->photo.noptlaws; ilaw++)
			if (dmod->photo.opttype[ilaw] == HARMLAMBERT || dmod->photo.opttype[ilaw] == HARMLOMMEL
					|| dmod->photo.opttype[ilaw] == HARMHAPKE
					|| dmod->photo.opttype[ilaw] == HARMKAAS)
				harmonic_scatlaw = 1;
		for (ilaw=0; ilaw<dmod->photo.nradlaws; ilaw++)
			if (dmod->photo.radtype[ilaw] == HARMCOSINE_DIFF)
				harmonic_scatlaw = 1;
	}
}
__global__ void calc_vertex_nrmls_krnl(struct mod_t *dmod)
{
	/* nv-threaded kernel */
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	double n[3];
	int j, k, naf, f;
	if (i < dmod->shape.comp[0].real.nv){
		n[0] = n[1] = n[2] = 0.0;
		naf = dmod->shape.comp[0].real.v[i].naf;
		for (j=0; j<naf; j++) {
			f = dmod->shape.comp[0].real.v[i].af[j];
			n[0] += dmod->shape.comp[0].real.f[f].n[0];
			n[1] += dmod->shape.comp[0].real.f[f].n[1];
			n[2] += dmod->shape.comp[0].real.f[f].n[2];
		}
		dev_normalize( n);
		for (k=0; k<=2; k++)
			dmod->shape.comp[0].real.v[i].n[k] = n[k];
	}
}
__global__ void calc_vertex_nrmls_f_krnl(struct mod_t *dmod)
{
	/* nv-threaded kernel */
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	float3 n;
	int j, naf, f;

	if (i < dnv){
		n.x = n.y = n.z = 0.0;
		naf = dmod->shape.comp[0].real.v[i].naf;
		for (j=0; j<naf; j++) {
			f = dmod->shape.comp[0].real.v[i].af[j];
			n.x += __double2float_rn(dmod->shape.comp[0].real.f[f].n[0]);
			n.y += __double2float_rn(dmod->shape.comp[0].real.f[f].n[1]);
			n.z += __double2float_rn(dmod->shape.comp[0].real.f[f].n[2]);
		}
		dev_normalize2( n);
		dmod->shape.comp[0].real.v[i].n[0] = n.x;
		dmod->shape.comp[0].real.v[i].n[1] = n.y;
		dmod->shape.comp[0].real.v[i].n[2] = n.z;
	}
}
__global__ void facet_krnl(struct par_t *dpar, struct mod_t *dmod)
{
	/* For each facet of this component, compute the outward unit normal,
	 * the area, the mean coordinates of the three corner vertices, and
	 * the corresponding angular coordinates (for some scattering laws)    */
	/* nf-threaded kernel */

	int f = blockIdx.x * blockDim.x + threadIdx.x;
	int j;

	if (f < dnf) {

		dmod->shape.comp[0].real.f[f].area = dev_facnrm(dmod->shape.comp[0].real, f);

		for (j=0; j<=2; j++)
			dmod->shape.comp[0].real.f[f].x[j] = (dmod->shape.comp[0].real.v[ dmod->shape.comp[0].real.f[f].v[0] ].x[j] +
					dmod->shape.comp[0].real.v[ dmod->shape.comp[0].real.f[f].v[1] ].x[j] +
					dmod->shape.comp[0].real.v[ dmod->shape.comp[0].real.f[f].v[2] ].x[j]   )/3;

		if (harmonic_scatlaw) {
			dmod->shape.comp[0].real.f[f].theta = atan2( sqrt(dmod->shape.comp[0].real.f[f].x[0]*dmod->shape.comp[0].real.f[f].x[0] +
					dmod->shape.comp[0].real.f[f].x[1]*dmod->shape.comp[0].real.f[f].x[1]   ),
					dmod->shape.comp[0].real.f[f].x[2]);
			dmod->shape.comp[0].real.f[f].phi = atan2( dmod->shape.comp[0].real.f[f].x[1], dmod->shape.comp[0].real.f[f].x[0]);
		}
	}
}
__global__ void facet_f_krnl(struct par_t *dpar, struct mod_t *dmod)
{
	/* For each facet of this component, compute the outward unit normal,
	 * the area, the mean coordinates of the three corner vertices, and
	 * the corresponding angular coordinates (for some scattering laws)    */
	/* nf-threaded kernel */

	int f = blockIdx.x * blockDim.x + threadIdx.x;
	int j;

	if (f < dnf) {

		dmod->shape.comp[0].real.f[f].area = dev_facnrm(dmod->shape.comp[0].real, f);

		for (j=0; j<=2; j++)
			dmod->shape.comp[0].real.f[f].x[j] = (dmod->shape.comp[0].real.v[ dmod->shape.comp[0].real.f[f].v[0] ].x[j] +
					dmod->shape.comp[0].real.v[ dmod->shape.comp[0].real.f[f].v[1] ].x[j] +
					dmod->shape.comp[0].real.v[ dmod->shape.comp[0].real.f[f].v[2] ].x[j]   )/3;

		if (harmonic_scatlaw) {
			dmod->shape.comp[0].real.f[f].theta = atan2( sqrt(dmod->shape.comp[0].real.f[f].x[0]*dmod->shape.comp[0].real.f[f].x[0] +
					dmod->shape.comp[0].real.f[f].x[1]*dmod->shape.comp[0].real.f[f].x[1]   ),
					dmod->shape.comp[0].real.f[f].x[2]);
			dmod->shape.comp[0].real.f[f].phi = atan2( dmod->shape.comp[0].real.f[f].x[1], dmod->shape.comp[0].real.f[f].x[0]);
		}
	}
}
__global__ void set_real_active_vert_krnl(struct mod_t *dmod)
{
	/* nv-threaded kernel */
	int v = blockIdx.x * blockDim.x + threadIdx.x;

	if (v < dnv) //dmod->shape.comp[0].real.nv)
		dmod->shape.comp[0].real.v[v].act = 1;
}
__global__ void set_real_active_facet_krnl(struct mod_t *dmod)
{
	/* nf-threaded kernel */
	int f = blockIdx.x * blockDim.x + threadIdx.x;

	if (f < dmod->shape.comp[0].real.nf)
		dmod->shape.comp[0].real.f[f].act = 1;
}
__global__ void set_real_active_side_krnl(struct mod_t *dmod)
{
	/* ns-threaded kernel */
	int k = blockIdx.x * blockDim.x + threadIdx.x;

	if (k < dmod->shape.comp[0].real.ns)
		dmod->shape.comp[0].real.s[k].act = 1;
}

__host__ void realize_mod_cuda( struct par_t *dpar, struct mod_t *dmod,
		unsigned char type)
{

	/*  We need to realize each model component as a polyhedral solid with
      triangular facets.  The first step is to call realize_coordinates,
      which computes the displacement of each vertex in this realization,
      represented as a base displacement plus a vertex deviation (either
      positive or negative) along a specified set of direction cosines.
      Additionally, for each facet it computes the outward unit normal,
      the area, the mean coordinates of the corner vertices, and (for
      some scattering laws) the corresponding angular coordinates.        */

	if (FLOAT)
		realize_coordinates_cuda(dpar, dmod, type);
	else
		realize_coordinates_cuda(dpar, dmod, type);

	/*  For multiple-component models, figure out which facets lie on
      the model's surface and which fall within some other component;
      such facets will have their "act" (active) flag reset to zero.   */

	check_surface_cuda(dmod);

	/*  Compute the area and moments (volume, center of mass, and
      inertia tensor) of each component and of the overall model  */

	compute_moments_cuda(dmod);
}

/*  Compute the vertex coordinates and (if necessary) facet angular coordinates
    for each component of the model's vertex realization                         */
__host__ void realize_coordinates_cuda( struct par_t *dpar, struct mod_t *dmod, unsigned char type)
{
	dim3 BLK, THD;
	THD.x = maxThreadsPerBlock;
	/* Loop over all model components, realizing each one as a polyhedral solid
	 * with triangular facets. Compute displacement of each vertex in this
	 * realization, represented as a base displacement plus a vertex deviation
	 * (positive or negative) along a specified set of direction cosines*/

	/*  Call Kernel to initialize flag for tiny/negative ellipsoid diameters  */
	set_diam_krnl<<<1,1>>>(dpar, dmod);//, dnv, dnf);
	checkErrorAfterKernelLaunch("set_diam_krnl");

	/* Note:  The CUDA-code assumes a single-component model for now.  */
	/* Loop over all model components, realizing each one as a polyhedral solid
	 * with triangular facets. Compute the displacement of each vertex in this
	 * realization, represented as a base displacement plus a vertex deviation
	 * (positive or negative) along a specified set of direction cosines.  */

	/* Copy nf and nv back from device copies dnf and dnv; used as launch
	 * parameters below */
	gpuErrchk(cudaMemcpyFromSymbol(&nv, dnv, sizeof(nv), 0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&nf, dnf, sizeof(nv), 0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&ns, dns, sizeof(nv), 0, cudaMemcpyDeviceToHost));
	/* Calculate launch parameters for kernels going over vertices (nvBLK) and
	 * facets (nfBLK) */
	nvTHD.x = nfTHD.x = THD.x = maxThreadsPerBlock;
	nvBLK.x = floor((THD.x - 1 + nv) / THD.x);
	nfBLK.x = floor((THD.x - 1 + nf) / THD.x);

	/* Check component type & create corresponding vertex realization.  */
	switch (type) {
	case ELLIPSE:
		/* To avoid negative diameters/very small positive diameters,
		 * adjust the function a[i] = 1/radius[i]^2 so it monotonically
		 * increases as diameter[i] decreases through zero and beyond,
		 * rather than being symmetric about zero diameter. Also set flag
		 * "baddiam" when any diameter is very small or negative, so that
		 * extra penalties can later be applied to this model. */

		/* Launch ellipse diameter kernel */
		ellipse_diameter_krnl<<<1,1>>>(dpar, dmod);
		checkErrorAfterKernelLaunch("ellipse_diameter_krnl");

		/* Kernel finds distance of each vertex to ellipsoid's center     */
		ellipse_distance_krnl<<<nvBLK,nvTHD>>>(dpar, dmod);
		checkErrorAfterKernelLaunch("ellipse_distance_krnl");

		/* Launch kernel to set real->scalefactor */
		ellipse_scalefactor_krnl<<<1,1>>>(dmod);
		checkErrorAfterKernelLaunch("ellipse_scalefactor_krnl");
		break;
	case OVOID:
		/*  Determine all shape parameters, making sure that none are out of bounds  */
		set_ovoid_parameters_krnl<<<1,1>>>(dpar, dmod);
		checkErrorAfterKernelLaunch("set_ovoid_parameters_krnl");

		/* Kernel finds distance of each vertex to ovoid's center     */
		ovoid_distance_krnl<<<nvBLK,nvTHD>>>(dpar, dmod);
		checkErrorAfterKernelLaunch("ovoid_distance_krnl");
		break;
	case HARMONIC:
		/* Kernel sets parameters associated with harmonic model     */
		harmonic_krnl<<<nvBLK,nvTHD>>>(dpar, dmod);
		checkErrorAfterKernelLaunch("harmonic_krnl");

		THD.x = 3;
		harmonic_scalefactor_krnl<<<1,THD>>>(dpar, dmod);
		checkErrorAfterKernelLaunch("harmonic_scalefactor_krnl");
		break;
	case VERTEX:
		/* The vertex type is its own realization, but we still need to update
		 * the values of the "scale factor" parameters and update any vertex
		 * deviations that have the '=' state    */
		vertex_update_dev_krnl<<<nvBLK,nvTHD>>>(dpar, dmod);
		checkErrorAfterKernelLaunch("vertex_update_dev_kernel");

		THD.x = 3;
		vertex_scalefactor_krnl<<<1,THD>>>(dpar, dmod);
		checkErrorAfterKernelLaunch("vertex_scalefactor_krnl");
		break;
	default:
		printf("realize_mod.c: don't know that component type\n");
	}      /* end of switch statement for component type */

	/*  Calculate vertex coordinates for this component  */
	calc_vertex_co_krnl<<<nvBLK,nvTHD>>>(dpar, dmod);
	checkErrorAfterKernelLaunch("calc_vertex_co_krnl");

	/* Use this component's rotational offset angles to create comp[c].m, the
	 * rotation matrix that will be applied to the vertex coordinates  */
	euler2mat_realize_mod_krnl<<<1,1>>>(dmod);
	checkErrorAfterKernelLaunch("dev_euler2mat");

	/* If needed, perform rotation on this component  */
	perform_rotation_krnl<<<nvBLK,nvTHD>>>(dpar, dmod);
	checkErrorAfterKernelLaunch("perform_rotation_krnl");

	/*  If needed, perform translation on this component  */
	perform_translation_krnl<<<nvBLK,nvTHD>>>(dpar, dmod);
	checkErrorAfterKernelLaunch("perform_translation_krnl, line 651");

	/* 	Figure out if optical/radar harmonic scattering laws are in use     *
	 *  and set the flag harmonic_scatlaw accordingly				        */
	set_optical_params_krnl<<<1,1>>>(dpar, dmod);
	checkErrorAfterKernelLaunch("set_optical_params_krnl, line 656");

	/* For each facet of this component, compute outward unit normal, area,
	 * mean coordinates of the three corner vertices, and corresponding angular
	 * coordinates (for some scattering laws)    */
	facet_krnl<<<nfBLK,nfTHD>>>(dpar, dmod);
	checkErrorAfterKernelLaunch("facet_krnl, line 662");

	/* Calculate vertex normals for this component as normalized sums of the
	 * facet normals for all facets attached to each vertex     */
	calc_vertex_nrmls_krnl<<<nvBLK,nvTHD>>>(dmod);
	checkErrorAfterKernelLaunch("calc_vertex_nrmls, line 667");
}
__host__ void realize_coordinates_f_cuda( struct par_t *dpar, struct mod_t *dmod, unsigned char type)
{
	dim3 BLK, THD;
	THD.x = maxThreadsPerBlock;
	/* Loop over all model components, realizing each one as a polyhedral solid
	 * with triangular facets. Compute displacement of each vertex in this
	 * realization, represented as a base displacement plus a vertex deviation
	 * (positive or negative) along a specified set of direction cosines*/

	/*  Call Kernel to initialize flag for tiny/negative ellipsoid diameters  */
	set_diam_krnl<<<1,1>>>(dpar, dmod);//, dnv, dnf);
	checkErrorAfterKernelLaunch("set_diam_krnl");

	/* Note:  The CUDA-code assumes a single-component model for now.  */
	/* Loop over all model components, realizing each one as a polyhedral solid
	 * with triangular facets. Compute the displacement of each vertex in this
	 * realization, represented as a base displacement plus a vertex deviation
	 * (positive or negative) along a specified set of direction cosines.  */

	/* Copy nf and nv back from device copies dnf and dnv; used as launch
	 * parameters below */
	gpuErrchk(cudaMemcpyFromSymbol(&nv, dnv, sizeof(nv), 0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&nf, dnf, sizeof(nv), 0, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpyFromSymbol(&ns, dns, sizeof(nv), 0, cudaMemcpyDeviceToHost));
	/* Calculate launch parameters for kernels going over vertices (nvBLK) and
	 * facets (nfBLK) */
	nvTHD.x = nfTHD.x = THD.x = maxThreadsPerBlock;
	nvBLK.x = floor((THD.x - 1 + nv) / THD.x);
	nfBLK.x = floor((THD.x - 1 + nf) / THD.x);

	/* Check component type & create corresponding vertex realization.  */
	switch (type) {
	case ELLIPSE:
		/* To avoid negative diameters/very small positive diameters,
		 * adjust the function a[i] = 1/radius[i]^2 so it monotonically
		 * increases as diameter[i] decreases through zero and beyond,
		 * rather than being symmetric about zero diameter. Also set flag
		 * "baddiam" when any diameter is very small or negative, so that
		 * extra penalties can later be applied to this model. */

		/* Launch ellipse diameter kernel */
		ellipse_diameter_krnl<<<1,1>>>(dpar, dmod);
		checkErrorAfterKernelLaunch("ellipse_diameter_krnl");

		/* Kernel finds distance of each vertex to ellipsoid's center     */
		ellipse_distance_krnl<<<nvBLK,nvTHD>>>(dpar, dmod);
		checkErrorAfterKernelLaunch("ellipse_distance_krnl");

		/* Launch kernel to set real->scalefactor */
		ellipse_scalefactor_krnl<<<1,1>>>(dmod);
		checkErrorAfterKernelLaunch("ellipse_scalefactor_krnl");
		break;
	case OVOID:
		/*  Determine all shape parameters, making sure that none are out of bounds  */
		set_ovoid_parameters_krnl<<<1,1>>>(dpar, dmod);
		checkErrorAfterKernelLaunch("set_ovoid_parameters_krnl");

		/* Kernel finds distance of each vertex to ovoid's center     */
		ovoid_distance_krnl<<<nvBLK,nvTHD>>>(dpar, dmod);
		checkErrorAfterKernelLaunch("ovoid_distance_krnl");
		break;
	case HARMONIC:
		/* Kernel sets parameters associated with harmonic model     */
		harmonic_krnl<<<nvBLK,nvTHD>>>(dpar, dmod);
		checkErrorAfterKernelLaunch("harmonic_krnl");

		THD.x = 3;
		harmonic_scalefactor_krnl<<<1,THD>>>(dpar, dmod);
		checkErrorAfterKernelLaunch("harmonic_scalefactor_krnl");
		break;
	case VERTEX:
		/* The vertex type is its own realization, but we still need to update
		 * the values of the "scale factor" parameters and update any vertex
		 * deviations that have the '=' state    */
		vertex_update_dev_krnl<<<nvBLK,nvTHD>>>(dpar, dmod);
		checkErrorAfterKernelLaunch("vertex_update_dev_kernel");

		BLK.x = 1;	THD.x = 3;
		vertex_scalefactor_krnl<<<1,1>>>(dpar, dmod);
		checkErrorAfterKernelLaunch("vertex_scalefactor_krnl");
		break;
	default:
		printf("realize_mod.c: don't know that component type\n");
	}      /* end of switch statement for component type */

	/*  Calculate vertex coordinates for this component  */
	calc_vertex_co_krnl<<<nvBLK,nvTHD>>>(dpar, dmod);
	checkErrorAfterKernelLaunch("calc_vertex_co_krnl");

	/* Use this component's rotational offset angles to create comp[c].m, the
	 * rotation matrix that will be applied to the vertex coordinates  */
	euler2mat_realize_mod_krnl<<<1,1>>>(dmod);
	checkErrorAfterKernelLaunch("dev_euler2mat");

	/* If needed, perform rotation on this component  */
	/* NOTE:  Still needs to be debugged, doesn't run as is */
	perform_rotation_f_krnl<<<nvBLK,nvTHD>>>(dpar, dmod);
	checkErrorAfterKernelLaunch("perform_rotation_krnl");

	/*  If needed, perform translation on this component  */
	perform_translation_krnl<<<nvBLK,nvTHD>>>(dpar, dmod);
	checkErrorAfterKernelLaunch("perform_translation_krnl, line 651");

	/* 	Figure out if optical/radar harmonic scattering laws are in use     *
	 *  and set the flag harmonic_scatlaw accordingly				        */
	set_optical_params_krnl<<<1,1>>>(dpar, dmod);
	checkErrorAfterKernelLaunch("set_optical_params_krnl, line 656");

	/* For each facet of this component, compute outward unit normal, area,
	 * mean coordinates of the three corner vertices, and corresponding angular
	 * coordinates (for some scattering laws)    */
	facet_krnl<<<nfBLK,nfTHD>>>(dpar, dmod);
	checkErrorAfterKernelLaunch("facet_krnl, line 662");

	/* Calculate vertex normals for this component as normalized sums of the
	 * facet normals for all facets attached to each vertex     */
	calc_vertex_nrmls_f_krnl<<<nvBLK,nvTHD>>>(dmod);
	checkErrorAfterKernelLaunch("calc_vertex_nrmls, line 667");
}
/*.....................................................................................*/

/*  Determine which vertices, facets, and sides of a multiple-component
    model lie interior to the model rather than on the model's surface,
    and reset their "act" (active) flags to zero                         */

__host__ void check_surface_cuda(struct mod_t *dmod)
{
	dim3 THD; THD.x = maxThreadsPerBlock;

	/* Calculate launch parameters */
	nvBLK.x = floor((maxThreadsPerBlock - 1 + nv) / maxThreadsPerBlock);
	nfBLK.x = floor((maxThreadsPerBlock - 1 + nf) / maxThreadsPerBlock);
	nsBLK.x = floor((maxThreadsPerBlock - 1 + ns) / maxThreadsPerBlock);

	/* 1-component model: flag all vertices and facets as active, then return  */
	set_real_active_vert_krnl<<<nvBLK,THD>>>(dmod);
	checkErrorAfterKernelLaunch("set_real_active_vert_krnl");

	set_real_active_facet_krnl<<<nfBLK,THD>>>(dmod);
	checkErrorAfterKernelLaunch("set_real_active_vert_krnl");

	set_real_active_side_krnl<<<nsBLK,THD>>>(dmod);
	checkErrorAfterKernelLaunch("set_real_active_side_krnl");

	return;
}

__global__ void comp_moments_1stinit_krnl(struct mod_t *dmod, int c) {
	/* Single-thread kernel */
	int j, k;
	if (threadIdx.x == 0) {
		dmod->shape.area = 0.0;
		dmod->shape.volume = 0.0;
		for (k=0; k<=2; k++) {
			dmod->shape.com[k] = 0.0;
			for (j=0; j<=2; j++)
				dmod->shape.inertia[k][j] = 0.0;
		}
	}
}
__global__ void comp_moments_2ndinit_krnl(struct mod_t *dmod, float area1,
		float area2, int c) {
	/* Single-threaded kernel - meant to initialize the individual component
	 * com and inertia arrays */
	if (threadIdx.x == 0) {
		int j, k;
		dmod->shape.comp[c].area = area1;
		dmod->shape.area = area2;
		dmod->shape.comp[0].volume = 0.0;
		for (k=0; k<=2; k++) {
			dmod->shape.comp[0].com[k] = 0.0;
			for (j=0; j<=2; j++)
				dmod->shape.comp[0].inertia[k][j] = 0.0;
		}
		dmod->shape.comp[0].area = 0.0; // actually 1st step in calculating surface area
	}
}
__global__ void comp_moments_facet_krnl(struct mod_t *dmod, int c, float *dvarr,
		float *dcom0, float *dcom1, float *dcom2, float *dI00, float *dI01,
		float *dI02, float *dI10, float *dI11, float *dI12, float *dI20,
		float *dI21, float *dI22)
{
	/* nf-threaded kernel */
	int f = blockIdx.x * blockDim.x + threadIdx.x;
	double dI[3][3], dcom[3], dv;
	if (f < dmod->shape.comp[0].real.nf)
	{
		dev_facmom(dmod->shape.comp[c].real.v[ dmod->shape.comp[0].real.f[f].v[0] ].x,
				dmod->shape.comp[c].real.v[ dmod->shape.comp[0].real.f[f].v[1] ].x,
				dmod->shape.comp[c].real.v[ dmod->shape.comp[0].real.f[f].v[2] ].x,
				dmod->shape.comp[c].real.f[f].n,	&dv, dcom, dI);

		/* Assign calculated dv, dcom, dI to each facet for later parallel reduction */
		dvarr[f]  	= (float) dv;
		dcom0[f]	= (float)dcom[0];
		dcom1[f]	= (float)dcom[1];
		dcom2[f]	= (float)dcom[2];
		dI00[f] 	= (float)dI[0][0];
		dI01[f] 	= (float)dI[0][1];
		dI02[f] 	= (float)dI[0][2];
		dI10[f] 	= (float)dI[1][0];
		dI11[f] 	= (float)dI[1][1];
		dI12[f] 	= (float)dI[1][2];
		dI20[f] 	= (float)dI[2][0];
		dI21[f] 	= (float)dI[2][1];
		dI22[f] 	= (float)dI[2][2];
	}
}
__global__ void comp_moments_facet_f_krnl(struct mod_t *dmod, int c, float *dvarr,
		float *dcom0, float *dcom1, float *dcom2, float *dI00, float *dI01,
		float *dI02, float *dI10, float *dI11, float *dI12, float *dI20,
		float *dI21, float *dI22)
{
	/* nf-threaded kernel */
	int f = blockIdx.x * blockDim.x + threadIdx.x;
	float3 dI[3], dcom, fv0, fv1, fv2, fn;
	float dv;
	if (f < dnf)
	{
		fv0.x = dmod->shape.comp[c].real.v[dmod->shape.comp[c].real.f[f].v[0]].x[0];
		fv0.y = dmod->shape.comp[c].real.v[dmod->shape.comp[c].real.f[f].v[0]].x[1];
		fv0.z = dmod->shape.comp[c].real.v[dmod->shape.comp[c].real.f[f].v[0]].x[2];

		fv1.x = dmod->shape.comp[c].real.v[dmod->shape.comp[c].real.f[f].v[1]].x[0];
		fv1.y = dmod->shape.comp[c].real.v[dmod->shape.comp[c].real.f[f].v[1]].x[1];
		fv1.z = dmod->shape.comp[c].real.v[dmod->shape.comp[c].real.f[f].v[1]].x[2];

		fv2.x = dmod->shape.comp[c].real.v[dmod->shape.comp[c].real.f[f].v[2]].x[0];
		fv2.y = dmod->shape.comp[c].real.v[dmod->shape.comp[c].real.f[f].v[2]].x[1];
		fv2.z = dmod->shape.comp[c].real.v[dmod->shape.comp[c].real.f[f].v[2]].x[2];

		fn.x = dmod->shape.comp[c].real.f[f].n[0];
		fn.y = dmod->shape.comp[c].real.f[f].n[1];
		fn.z = dmod->shape.comp[c].real.f[f].n[2];

		dev_facmom_f(fv0, fv1, fv2, fn, &dv, &dcom, dI);

		/* Assign calculated dv, dcom, dI to each facet for later parallel reduction */
		dvarr[f]  	= dv;
		dcom0[f]	= dcom.x;
		dcom1[f]	= dcom.y;
		dcom2[f]	= dcom.z;
		dI00[f] 	= dI[0].x;
		dI01[f] 	= dI[0].y;
		dI02[f] 	= dI[0].z;
		dI10[f] 	= dI[1].x;
		dI11[f] 	= dI[1].y;
		dI12[f] 	= dI[1].z;
		dI20[f] 	= dI[2].x;
		dI21[f] 	= dI[2].y;
		dI22[f] 	= dI[2].z;
	}
}
__global__ void comp_moments_facets_old_krnl(struct mod_t *dmod)
{
	/* nf-threaded kernel */
	int f = blockIdx.x * blockDim.x + threadIdx.x;
	int j, k;
	double dI[3][3], dcom[3], dv;
	if (f < dmod->shape.comp[0].real.nf)
	{
		/* Calculate surface area for this component; for active facets, also add
		 * the contributions to the area of the overall model    */

		dmod->shape.comp[0].area += dmod->shape.comp[0].real.f[f].area;
		if (dmod->shape.comp[0].real.f[f].act)
			dmod->shape.area += dmod->shape.comp[0].real.f[f].area;

		dev_facmom( dmod->shape.comp[0].real.v[ dmod->shape.comp[0].real.f[f].v[0] ].x,
				dmod->shape.comp[0].real.v[ dmod->shape.comp[0].real.f[f].v[1] ].x,
				dmod->shape.comp[0].real.v[ dmod->shape.comp[0].real.f[f].v[2] ].x,
				dmod->shape.comp[0].real.f[f].n,	&dv, dcom, dI);
		dmod->shape.comp[0].volume += dv;
		for (j=0; j<=2; j++) {
			dmod->shape.comp[0].com[j] += dcom[j];
			for (k=0; k<=2; k++)
				dmod->shape.comp[0].inertia[j][k] += dI[j][k];
		}

		if (dmod->shape.comp[0].real.f[f].act) {
			dmod->shape.volume += dv;
			for (j=0; j<=2; j++) {
				dmod->shape.com[j] += dcom[j];
				for (k=0; k<=2; k++)
					dmod->shape.inertia[j][k] += dI[j][k];
			}
		}
	}
}
__global__ void comp_moments_facets_atomics_krnl(struct mod_t *dmod)
{
	/* nf-threaded kernel */
	int f = blockIdx.x * blockDim.x + threadIdx.x;
	int j, k;
	double dI[3][3], dcom[3], dv;
	if (f < dmod->shape.comp[0].real.nf)
	{
		/* Calculate surface area for this component; for active facets, also add
		 * the contributions to the area of the overall model    */
		atomicAdd(&rm_area, (float)dmod->shape.comp[0].real.f[f].area);

		if (dmod->shape.comp[0].real.f[f].act)
			atomicAdd(&rm_ifarea, (float)dmod->shape.comp[0].real.f[f].area);

		dev_facmom( dmod->shape.comp[0].real.v[ dmod->shape.comp[0].real.f[f].v[0] ].x,
				dmod->shape.comp[0].real.v[ dmod->shape.comp[0].real.f[f].v[1] ].x,
				dmod->shape.comp[0].real.v[ dmod->shape.comp[0].real.f[f].v[2] ].x,
				dmod->shape.comp[0].real.f[f].n,	&dv, dcom, dI);

		atomicAdd(&rm_vol, (float)dv);

		for (j=0; j<=2; j++) {
			atomicAdd(&rm_dcom[j], (float)dcom[j]);
			for (k=0; k<=2; k++)
				atomicAdd(&rm_dI[j][k], (float)dI[j][k]);
			dmod->shape.comp[0].inertia[j][k] += dI[j][k];
		}

		if (dmod->shape.comp[0].real.f[f].act) {
			atomicAdd(&rm_vol, (float)dv);
			for (j=0; j<=2; j++) {
				atomicAdd(&rm_ifdcom[j], dcom[j]);
				for (k=0; k<=2; k++)
					atomicAdd(&rm_ifdI[j][k], (float)dI[j][k]);
			}
		}
	}
}
__global__ void comp_moments_facets_at2_krnl(struct mod_t *dmod) {
	/* Single-threaded kernel */
	if (threadIdx.x ==0) {
		int i, j;
		dmod->shape.comp[0].area = (double)rm_area;
		dmod->shape.area = (double)rm_ifarea;
		dmod->shape.comp[0].volume = (double)rm_vol;
		dmod->shape.volume = (double)rm_ifvol;

		for (i=0; i<3; i++) {
			dmod->shape.comp[0].com[i] = rm_dcom[i];
			dmod->shape.com[i] = rm_ifdcom[i];
			for (j=0; j<3; j++) {
				dmod->shape.comp[0].inertia[i][j] = rm_dI[i][j];
				dmod->shape.inertia[i][j] = rm_ifdI[i][j];
			}
		}
	}
}
__global__ void comp_moments_com_krnl(struct mod_t *dmod)
{
	/* Single-thread kernel */
	if (threadIdx.x == 0) {
		int j;

		for (j=0; j<=2; j++) {
			dmod->shape.comp[0].com[j] /= dmod->shape.comp[0].volume;
			dmod->shape.com[j] /= dmod->shape.volume;
		}
		j = 2;
	}
}

/*  Compute the area and the 0,1,2-order moments (volume, center of mass, and
    inertia tensor) of each component and of the overall model, assuming uniform
    density and ignoring interior facets' contributions to the overall model      */

__host__ void compute_moments_cuda( struct mod_t *dmod)
{
	float area1=0.0, area2=0.0, *dv, *dcom0, *dcom1, *dcom2, *dI00, *dI01, *dI02,
			*dI10, *dI11, *dI12, *dI20, *dI21, *dI22;
	int c=0, size;
	/*  Initialize the model's surface area, volume, center-of-mass (COM)
	 * displacement, and inertia tensor  */
	comp_moments_1stinit_krnl<<<1,1>>>(dmod, c);
	checkErrorAfterKernelLaunch("comp_moments_init_krnl, line 945");

	gpuErrchk(cudaMemcpyFromSymbol(&size, dnf, sizeof(int), 0, cudaMemcpyDeviceToHost));

	/* CUDA note:  Only single-component models for now.
	 * Loop over all model components, computing areas and moments (volume,
	 * center of mass, and inertia tensor); COM and inertia tensor are computed
	 * assuming uniform density. For multiple-component models, when computing
	 * the area and the moments for overall model, ignore facets interior to
	 * the model (i.e., that are inside some other component).         */
	/* Note that area2 (area of active facets summed up) is not currently
	 * implemented. A single-component model is assumed, in which case every
	 * facet is active and area1=area2 */
	//	for (c=0; c<dmod->shape.ncomp; c++) {
	area1 = compute_model_area(dmod, c, size);
	/*area2 = compute_model_area2(dmod, c, size);*/
	area2 = area1;

	/* Allocate temporary dv, dcom, dI pointers */
	cudaCalloc1((void**)&dv, sizeof(float), size);
	cudaCalloc1((void**)&dcom0, sizeof(float), size);
	cudaCalloc1((void**)&dcom1, sizeof(float), size);
	cudaCalloc1((void**)&dcom2, sizeof(float), size);
	cudaCalloc1((void**)&dI00, sizeof(float), size);
	cudaCalloc1((void**)&dI01, sizeof(float), size);
	cudaCalloc1((void**)&dI02, sizeof(float), size);
	cudaCalloc1((void**)&dI10, sizeof(float), size);
	cudaCalloc1((void**)&dI11, sizeof(float), size);
	cudaCalloc1((void**)&dI12, sizeof(float), size);
	cudaCalloc1((void**)&dI20, sizeof(float), size);
	cudaCalloc1((void**)&dI21, sizeof(float), size);
	cudaCalloc1((void**)&dI22, sizeof(float), size);

	/* Set area and initialize per-component COM and Inertia arrays */
	comp_moments_2ndinit_krnl<<<1,1>>>(dmod, area1, area2, c);
	checkErrorAfterKernelLaunch("comp_moments_2ndinit_krnl (compute_moments_cuda)");

	/* Load the temporary arrays with data */
//	if (FLOAT)
//		comp_moments_facet_f_krnl<<<nfBLK,nfTHD>>>(dmod, c, dv, dcom0, dcom1, dcom2,
//				dI00, dI01, dI02, dI10, dI11, dI12, dI20, dI21, dI22);
//	else
		comp_moments_facet_krnl<<<nfBLK,nfTHD>>>(dmod, c, dv, dcom0, dcom1, dcom2,
				dI00, dI01, dI02, dI10, dI11, dI12, dI20, dI21, dI22);
	checkErrorAfterKernelLaunch("comp_moments_facets_krnl (compute_moments_cuda)");

	/* Calculate surface area for this component; for active facets, also add
	 * the contributions to the area of the overall model    */

	dvdI_reduce_single(dmod, dv, dcom0,	dcom1, dcom2, dI00, dI01, dI02,
			dI10, dI11, dI12, dI20, dI21, dI22, size, c);

	/* This kernel computes the overall COM vector */
	comp_moments_com_krnl<<<1,1>>>(dmod);
	checkErrorAfterKernelLaunch("comp_moments_facets_krnl, line 963");

	/* Free up the temporary arrays */
	cudaFree(dv);
	cudaFree(dcom0);	cudaFree(dcom1);	cudaFree(dcom2);
	cudaFree(dI00);		cudaFree(dI01);		cudaFree(dI02);
	cudaFree(dI10);		cudaFree(dI11);		cudaFree(dI12);
	cudaFree(dI20);		cudaFree(dI21);		cudaFree(dI22);
}

/*  Find all real roots of a cubic equation, using methods given in section 5.6 of
    Numerical Recipes in C.  Element 3 of the input coeff vector is the cubic
    coefficient while element 0 is the constant term.  Up to three real roots are
    stored in the output realroot vector, with any unused elements set to a large
    negative dummy value.  The return value is the number of real roots found.
    The routine includes several tests for coefficients that are equal to zero;
    those tests assume that nonzero coefficients are of order unity.                */
__device__ int cubic_realroots_cuda( double *coeff, double *realroot)
{
	int nrealroots, bsign;
	double a, b, c, discriminant, q, qsqrt, r, r2minusq3, rsign, s, t, theta;
	nrealroots = 0;
	realroot[0] = realroot[1] = realroot[2] = -HUGENUMBER;

	if (fabs(coeff[3]) < SMALLCOEFF3) {
		/*  cubic term is zero  */
		a = coeff[2];
		b = coeff[1];
		c = coeff[0];

		if (fabs(a) < SMALLVAL) {

			if (fabs(b) < SMALLVAL) {
				/*  Error: the cubic, quadratic, and linear terms are zero  */
				if (fabs(c) < SMALLVAL)
					printf("cubic_realroots in realize_mod.c: all four coefficients are zero\n");
				else
					printf("cubic_realroots in realize_mod.c: only the constant term is nonzero\n");

			} else {
				/*  linear equation  */
				realroot[0] = -c/b;
				nrealroots = 1;
			}

		} else {
			/*  quadratic equation  */
			discriminant = b*b - 4*a*c;
			if (discriminant < 0.0)
				printf("cubic_realroots in realize_mod.c: quadratic equation has no real roots\n");
			if (fabs(b) < SMALLVAL) {
				realroot[0] = sqrt(discriminant)/(2*a);
				realroot[1] = -realroot[0];
			} else {
				bsign = (b < 0.0) ? -1 : 1;
				q = -0.5*(b + bsign*sqrt(discriminant));
				realroot[0] = q/a;
				realroot[1] = c/q;
			}
			nrealroots = 2;
		}
	} else {
		/*  cubic term is nonzero: scale to standard form x^3 + ax^2 + b^x + c = 0  */
		a = coeff[2]/coeff[3];
		b = coeff[1]/coeff[3];
		c = coeff[0]/coeff[3];

		/* Check if there is one real root or three. Write out test quantity
		 * r^2 - q^3 explicitly in terms of coefficients a, b, and c in order
		 * to cancel high-order terms and thus reduce the likelihood of
		 * roundoff problems           */

		q = (a*a - 3*b)/9;
		r = (2*a*a*a - 9*a*b + 27*c)/54;

		r2minusq3 = (4*a*a*a*c - a*a*b*b - 18*a*b*c + 27*c*c + 4*b*b*b)/108;
		if (r2minusq3 >= 0.0) {
			/*  one real root  */
			rsign = (r < 0.0) ? -1 : 1;
			s = -rsign*pow( fabs(r) + sqrt(r2minusq3), 1.0/3);
			t = (fabs(s) >= SMALLVAL) ? q/s : 0.0;
			realroot[0] = s + t - a/3;
			nrealroots = 1;
		} else {
			/*  three real roots  */
			qsqrt = sqrt(q);
			theta = acos(r/(q*qsqrt));
			realroot[0] = -2*qsqrt*cos(theta/3) - a/3;
			realroot[1] = -2*qsqrt*cos((theta + 2*PIE)/3) - a/3;
			realroot[2] = -2*qsqrt*cos((theta - 2*PIE)/3) - a/3;
			nrealroots = 3;
		}
	}
	return nrealroots;
}
__device__ int cubic_realroots_f_cuda( float4 coeff, float3 *realroot)
{
	int nrealroots, bsign;
	float a, b, c, discriminant, q, qsqrt, r, r2minusq3, rsign, s, t, theta;
	nrealroots = 0;
	realroot->x = realroot->y = realroot->z = -HUGENUMBER;

	if (fabs(coeff.z) < SMALLCOEFF3) {
		/*  cubic term is zero  */
		a = coeff.y;
		b = coeff.x;
		c = coeff.w;

		if (fabs(a) < SMALLVAL) {

			if (fabs(b) < SMALLVAL) {
				/*  Error: the cubic, quadratic, and linear terms are zero  */
				if (fabs(c) < SMALLVAL)
					printf("cubic_realroots in realize_mod.c: all four coefficients are zero\n");
				else
					printf("cubic_realroots in realize_mod.c: only the constant term is nonzero\n");

			} else {
				/*  linear equation  */
				realroot->x = -c/b;
				nrealroots = 1;
			}

		} else {
			/*  quadratic equation  */
			discriminant = b*b - 4*a*c;
			if (discriminant < 0.0)
				printf("cubic_realroots in realize_mod.c: quadratic equation has no real roots\n");
			if (fabs(b) < SMALLVAL) {
				realroot->x = sqrt(discriminant)/(2*a);
				realroot->y = -realroot->x;
			} else {
				bsign = (b < 0.0) ? -1 : 1;
				q = -0.5*(b + bsign*sqrt(discriminant));
				realroot->x = q/a;
				realroot->y = c/q;
			}
			nrealroots = 2;
		}
	} else {
		/*  cubic term is nonzero: scale to standard form x^3 + ax^2 + b^x + c = 0  */
		a = coeff.y/coeff.z;
		b = coeff.x/coeff.z;
		c = coeff.w/coeff.z;

		/* Check if there is one real root or three. Write out test quantity
		 * r^2 - q^3 explicitly in terms of coefficients a, b, and c in order
		 * to cancel high-order terms and thus reduce the likelihood of
		 * roundoff problems           */
		q = (a*a - 3*b)/9;
		r = (2*a*a*a - 9*a*b + 27*c)/54;

		r2minusq3 = (4*a*a*a*c - a*a*b*b - 18*a*b*c + 27*c*c + 4*b*b*b)/108;
		if (r2minusq3 >= 0.0) {
			/*  one real root  */
			rsign = (r < 0.0) ? -1 : 1;
			s = -rsign*powf( fabs(r) + sqrt(r2minusq3), 1.0/3);
			t = (fabs(s) >= SMALLVAL) ? q/s : 0.0;
			realroot->x = s + t - a/3;
			nrealroots = 1;
		} else {
			/*  three real roots  */
			qsqrt = sqrt(q);
			theta = acos(r/(q*qsqrt));
			realroot->x = -2*qsqrt*cos(theta/3) - a/3;
			realroot->y = -2*qsqrt*cos((theta + 2*PIE)/3) - a/3;
			realroot->z = -2*qsqrt*cos((theta - 2*PIE)/3) - a/3;
			nrealroots = 3;
		}
	}
	return nrealroots;
}

#undef HAIRWIDTH
#undef SMALLRATIO
#undef SMALLOVOIDK1
#undef SMALLOVOIDK2
#undef OVOIDTOL
#undef MAXEDGE
#undef EDGETOL
#undef RTOL
#undef SMALLCOEFF3

__device__ double dev_facnrm( struct vertices_t verts, int fi)
{
	int i;
	double a[3], b[3], area;

	for (i=0; i<=2; i++) {
		a[i] = verts.v[verts.f[fi].v[1]].x[i] - verts.v[verts.f[fi].v[0]].x[i];
		b[i] = verts.v[verts.f[fi].v[2]].x[i] - verts.v[verts.f[fi].v[1]].x[i];
	}
	area = 0.5*dev_cross( verts.f[fi].n, a, b);
	dev_normalize( verts.f[fi].n);
	return area;
}
__device__ float dev_facnrm_f( struct vertices_t verts, int fi)
{
	float3 a, b;
	float area;
	int3 idx;
	float3 xv0, xv1, xv2;
	idx.x = verts.f[fi].v[0];
	idx.y = verts.f[fi].v[1];
	idx.z = verts.f[fi].v[2];
	xv0.x = __double2float_rn(verts.v[idx.x].x[0]);
	xv0.y = __double2float_rn(verts.v[idx.x].x[1]);
	xv0.z = __double2float_rn(verts.v[idx.x].x[2]);
	xv1.x = __double2float_rn(verts.v[idx.y].x[0]);
	xv1.y = __double2float_rn(verts.v[idx.y].x[1]);
	xv1.z = __double2float_rn(verts.v[idx.y].x[2]);
	xv2.x = __double2float_rn(verts.v[idx.z].x[0]);
	xv2.y = __double2float_rn(verts.v[idx.z].x[1]);
	xv2.z = __double2float_rn(verts.v[idx.z].x[2]);

	a.x =  verts.v[verts.f[fi].v[1]].x[0] - verts.v[verts.f[fi].v[0]].x[0];
	a.y =  verts.v[verts.f[fi].v[1]].x[1] - verts.v[verts.f[fi].v[0]].x[1];
	a.z =  verts.v[verts.f[fi].v[1]].x[2] - verts.v[verts.f[fi].v[0]].x[2];

	b.x  = verts.v[verts.f[fi].v[2]].x[0] - verts.v[verts.f[fi].v[1]].x[0];
	b.y =  verts.v[verts.f[fi].v[2]].x[1] - verts.v[verts.f[fi].v[1]].x[1];
	b.z =  verts.v[verts.f[fi].v[2]].x[2] - verts.v[verts.f[fi].v[1]].x[2];

	area = 0.5*dev_cross_f( verts.f[fi].n, a, b);
	dev_normalize( verts.f[fi].n);
	return area;
}
__device__ double dev_cross( double z[3], double x[3], double y[3])
{
	double zz[3];

	zz[0] = x[1]*y[2]-x[2]*y[1];
	zz[1] = x[2]*y[0]-x[0]*y[2];
	zz[2] = x[0]*y[1]-x[1]*y[0];
	z[0] = zz[0];
	z[1] = zz[1];
	z[2] = zz[2];
	return sqrt(z[0]*z[0]+z[1]*z[1]+z[2]*z[2]);
}
__device__ float dev_cross_f( double z[3], float3 x, float3 y)
{
	float3 zz;
	float area;
	zz.x = x.y*y.z-x.z*y.y;
	zz.y = x.z*y.x-x.x*y.z;
	zz.z = x.x*y.y-x.y*y.x;
	z[0] = zz.x;
	z[1] = zz.y;
	z[2] = zz.z;

	area = sqrt(z[0]*z[0] + z[1]*z[1] + z[2]*z[2]);
	return area;
}
