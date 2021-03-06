extern "C" {
#include "../shape/head.h"
}

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
//#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

__host__ void mnbrak_gpu(
		double *ax,
		double *bx,
		double *cx,
		double *fa,
		double *fb,
		double *fc,
		double (*func)
		(double,
				struct vertices_t**,
				unsigned char*,
				unsigned char*,
				int*,
				int*,
				int*,
				int,
				int,
				cudaStream_t*),
				struct vertices_t **verts,
				unsigned char *htype,
				unsigned char *dtype,
				int *nframes,
				int *nviews,
				int *lc_n,
				int nsets,
				int nf,
				cudaStream_t *bf_stream )
{
	double ulim,u,r,q,fu,dum;

	*fa = (*func)(*ax, verts, htype, dtype, nframes, nviews, lc_n, nsets, nf, bf_stream);
//	printf("mnbrak ax, %g\n", *ax);
//	printf("mnbrak fa, %g\n", *fa);
	*fb = (*func)(*bx, verts, htype, dtype, nframes, nviews, lc_n, nsets, nf, bf_stream);
//	printf("mnbrak bx, %g\n", *bx);
//	printf("mnbrak fb, %g\n", *fb);
	if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}

	*cx = (*bx)+GOLD*(*bx-*ax);
	*fc = (*func)(*cx, verts, htype, dtype, nframes, nviews, lc_n, nsets, nf, bf_stream);

	while (*fb > *fc) {

		r = (*bx-*ax) * (*fb-*fc);
		q = (*bx-*cx) * (*fb-*fa);
		u = (*bx) - ((*bx-*cx)*q - (*bx-*ax)*r)/
				(2.0 * SIGN( MAX( fabs(q-r),TINY), q-r));
		ulim = (*bx) + GLIMIT * (*cx-*bx);

		if ((*bx-u)*(u-*cx) > 0.0) {
			fu=(*func)(u, verts, htype, dtype, nframes, nviews, lc_n, nsets, nf, bf_stream);

			if (fu < *fc) {
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			} else if (fu > *fb) {
				*cx = u;
				*fc = fu;
				return;
			}
			u = (*cx) + GOLD * (*cx-*bx);
			fu = (*func)(u, verts, htype, dtype, nframes, nviews, lc_n, nsets, nf, bf_stream);

		} else if ((*cx-u) * (u-ulim) > 0.0) {
			fu = (*func)(u, verts, htype, dtype, nframes, nviews, lc_n, nsets, nf, bf_stream);

			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
				SHFT(*fb,*fc,fu,(*func)(u, verts, htype, dtype, nframes, nviews, lc_n, nsets, nf, bf_stream))
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			u = ulim;
			fu = (*func)(u, verts, htype, dtype, nframes, nviews, lc_n, nsets, nf, bf_stream);
		} else {
			u = (*cx)+GOLD*(*cx-*bx);
			fu = (*func)(u, verts, htype, dtype, nframes, nviews, lc_n, nsets, nf, bf_stream);
		}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
	}
}

__host__ void mnbrak_gpu_dbg(
		double *ax,
		double *bx,
		double *cx,
		double *fa,
		double *fb,
		double *fc,
		double (*func)
		(double,
				struct vertices_t**,
				unsigned char*,
				unsigned char*,
				int*,
				int*,
				int*,
				int,
				int,
				cudaStream_t*),
				struct vertices_t **verts,
				unsigned char *htype,
				unsigned char *dtype,
				int *nframes,
				int *nviews,
				int *lc_n,
				int nsets,
				int nf,
				cudaStream_t *bf_stream )
{
	double ulim,u,r,q,fu,dum;

	int dbg=0;
	printf("mnbrak start\n");

	*fa = (*func)(*ax, verts, htype, dtype, nframes, nviews, lc_n, nsets, nf, bf_stream);
	*fb = (*func)(*bx, verts, htype, dtype, nframes, nviews, lc_n, nsets, nf, bf_stream);
	if (*fb > *fa) {
		printf("if (*fb > *fa), 1\n");
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	else	printf("if (*fb > *fa), 0\n");

	*cx = (*bx)+GOLD*(*bx-*ax);
	*fc = (*func)(*cx, verts, htype, dtype, nframes, nviews, lc_n, nsets, nf, bf_stream);

	printf("prior to while loop\n");
	printf("ax, %3.8g\n", *ax);
	printf("bx, %3.8g\n", *bx);
	printf("cx, %3.8g\n", *cx);

	while (*fb > *fc) {
		dbg++;
		printf("start of loop #, %i\n", dbg);

		r = (*bx-*ax) * (*fb-*fc);
		q = (*bx-*cx) * (*fb-*fa);
		u = (*bx) - ((*bx-*cx)*q - (*bx-*ax)*r)/
				(2.0 * SIGN( MAX( fabs(q-r),TINY), q-r));
		ulim = (*bx) + GLIMIT * (*cx-*bx);

		printf("r, %3.8g\n", r);
		printf("q, %3.8g\n", q);
		printf("u, %3.8g\n", u);
		printf("ulim, %3.8g\n", ulim);

		if ((*bx-u)*(u-*cx) > 0.0) {

			printf("if ((*bx-u)*(u-*cx) > 0.0), 1\n");

			fu=(*func)(u, verts, htype, dtype, nframes, nviews, lc_n, nsets, nf, bf_stream);

			printf("fu, %3.8g\n", fu);

			if (fu < *fc) {

				printf("if (fu < *fc), 1\n");

				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;

				printf("*ax, %3.8g\n", *ax);
				printf("*bx, %3.8g\n", *bx);
				printf("*fa, %3.8g\n", *fa);
				printf("*fb, %3.8g\n", *fb);

				return;
			} else if (fu > *fb) {

				printf("if ((*bx-u)*(u-*cx) > 0.0), 0\n");
				printf("if (fu > *fb), 1\n");

				*cx = u;
				*fc = fu;

				printf("*cx, %3.8g\n", *cx);
				printf("*fc, %3.8g\n", *fc);

				return;
			}

			u = (*cx) + GOLD * (*cx-*bx);
			fu = (*func)(u, verts, htype, dtype, nframes, nviews, lc_n, nsets, nf, bf_stream);

			printf("u, %3.8g\n", u);
			printf("fu, %3.8g\n", fu);

		} else if ((*cx-u) * (u-ulim) > 0.0) {
			printf("if ((*bx-u)*(u-*cx) > 0.0), 0\n");
			printf("else if ((*cx-u)*(u-ulim) > 0.0), 1\n");

			fu = (*func)(u, verts, htype, dtype, nframes, nviews, lc_n, nsets, nf, bf_stream);

			printf("fu, %3.8g\n", fu);

			if (fu < *fc) {
				printf("if (fu < *fc), 1\n");

				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
				SHFT(*fb,*fc,fu,(*func)(u, verts, htype, dtype, nframes, nviews, lc_n, nsets, nf, bf_stream))
			}
			else printf("if (fu < *fc), 0\n");

		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {

			printf("if ((*bx-u)*(u-*cx) > 0.0), 0\n");
			printf("if ((u-ulim)*(ulim-*cx) > 0.0), 1\n");

			u = ulim;
			fu = (*func)(u, verts, htype, dtype, nframes, nviews, lc_n, nsets, nf, bf_stream);

			printf("u, %3.8g\n", u);
			printf("fu, %3.8g\n", fu);

		} else {
			printf("if ((*bx-u)*(u-*cx) > 0.0), 0\n");

			u = (*cx)+GOLD*(*cx-*bx);
			fu = (*func)(u, verts, htype, dtype, nframes, nviews, lc_n, nsets, nf, bf_stream);

			printf("u, %3.8g\n", u);
			printf("fu, %3.8g\n", fu);
		}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
	}
	printf("loops in mnbrak_gpu, %i\n", dbg);
}

__host__ void mnbrak_MFS_gpu(
		double *ax,
		double *bx,
		double *cx,
		double *fa,
		double *fb,
		double *fc,
		double (*func)
		(double,
				struct vertices_t**,
				int*,
				int,
				int,
				cudaStream_t*),
				struct vertices_t **verts,
				int *nviews,
				int nsets,
				int nf,
				cudaStream_t *bf_stream )
{
	double ulim,u,r,q,fu,dum;

	*fa = (*func)(*ax, verts, nviews, nsets, nf, bf_stream);
	*fb = (*func)(*bx, verts, nviews, nsets, nf, bf_stream);
	if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}

	*cx = (*bx)+GOLD*(*bx-*ax);
	*fc = (*func)(*cx, verts, nviews, nsets, nf, bf_stream);

	while (*fb > *fc) {

		r = (*bx-*ax) * (*fb-*fc);
		q = (*bx-*cx) * (*fb-*fa);
		u = (*bx) - ((*bx-*cx)*q - (*bx-*ax)*r)/
				(2.0 * SIGN( MAX( fabs(q-r),TINY), q-r));
		ulim = (*bx) + GLIMIT * (*cx-*bx);

		if ((*bx-u)*(u-*cx) > 0.0) {
			fu=(*func)(u, verts, nviews, nsets, nf, bf_stream);

			if (fu < *fc) {
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			} else if (fu > *fb) {
				*cx = u;
				*fc = fu;
				return;
			}
			u = (*cx) + GOLD * (*cx-*bx);
			fu = (*func)(u, verts, nviews, nsets, nf, bf_stream);

		} else if ((*cx-u) * (u-ulim) > 0.0) {
			fu = (*func)(u, verts, nviews, nsets, nf, bf_stream);

			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
				SHFT(*fb,*fc,fu,(*func)(u, verts, nviews, nsets, nf, bf_stream))
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			u = ulim;
			fu = (*func)(u, verts, nviews, nsets, nf, bf_stream);
		} else {
			u = (*cx)+GOLD*(*cx-*bx);
			fu = (*func)(u, verts, nviews, nsets, nf, bf_stream);
		}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
	}
}

#undef GOLD
#undef GLIMIT
#undef TINY
#undef MAX
#undef SIGN
#undef SHFT
