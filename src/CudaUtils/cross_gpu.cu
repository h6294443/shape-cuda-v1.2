extern "C" {
#include "../shape/head.h"
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
