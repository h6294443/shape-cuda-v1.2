extern "C" {
#include "../shape/head.h"
}
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
