extern "C" {
#include "../shape/head.h"
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
