extern "C" {
#include "../shape/head.h"
}

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
//#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);


__host__ void mnbrak_pthreads(
		/* Arguments to mnbrak */
		double *ax,
		double *bx,
		double *cx,
		double *fa,
		double *fb,
	    double *fc,
	    /* Objective function with data types */
	    double (*func)(
	    		double,
	    		struct vertices_t**,
	    		struct vertices_t**,
	    		unsigned char*,
	    		unsigned char*,
	    		unsigned char*,
	    		int*,
	    		int*,
	    		int*,
	    		int*,
	    		int,
	    		int,
	    		int,
	    		pthread_t,
	    		pthread_t,
	    		cudaStream_t*,
	    		cudaStream_t*),

	    		/* And now the arguments to objective function */
	    		struct vertices_t **verts0,
	    		struct vertices_t **verts1,
	    		unsigned char *htype,
	    		unsigned char *dtype0,
	    		unsigned char *dtype1,
	    		int *nframes,
	    		int *nviews,
	    		int *lc_n,
	    		int *GPUID,
	    		int nsets,
	    		int nf,
	    		int max_frames,
	    		pthread_t thread1,
	    		pthread_t thread2,
	    		cudaStream_t *gpu0_stream,
	    		cudaStream_t *gpu1_stream)
{
	double ulim,u,r,q,fu,dum;
	*fa=(*func)(*ax, verts0, verts1, htype, dtype0, dtype1, nframes, nviews,
			lc_n, GPUID, nsets, nf, max_frames, thread1, thread2, gpu0_stream,
			gpu1_stream);
	*fb=(*func)(*bx, verts0, verts1, htype, dtype0, dtype1, nframes, nviews,
			lc_n, GPUID, nsets, nf, max_frames, thread1, thread2, gpu0_stream,
			gpu1_stream);
	if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx=(*bx)+GOLD*(*bx-*ax);
 	*fc=(*func)(*cx, verts0, verts1, htype, dtype0, dtype1, nframes, nviews,
 			lc_n, GPUID, nsets, nf, max_frames, thread1, thread2, gpu0_stream,
 			gpu1_stream);
	while (*fb > *fc) {
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
			(2.0*SIGN(MAX(fabs(q-r),TINY),q-r));
		ulim=(*bx)+GLIMIT*(*cx-*bx);
		if ((*bx-u)*(u-*cx) > 0.0) {
			fu=(*func)(u, verts0, verts1, htype, dtype0, dtype1, nframes,
					nviews, lc_n, GPUID, nsets, nf, max_frames, thread1,
					thread2, gpu0_stream, gpu1_stream);
			if (fu < *fc) {
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			} else if (fu > *fb) {
				*cx=u;
				*fc=fu;
				return;
			}
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u, verts0, verts1, htype, dtype0, dtype1, nframes,
					nviews, lc_n, GPUID, nsets, nf, max_frames, thread1,
					thread2, gpu0_stream, gpu1_stream);
		} else if ((*cx-u)*(u-ulim) > 0.0) {
			fu=(*func)(u, verts0, verts1, htype, dtype0, dtype1, nframes,
					nviews, lc_n, GPUID, nsets, nf, max_frames, thread1,
					thread2, gpu0_stream, gpu1_stream);
			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
				SHFT(*fb,*fc,fu,(*func)(u, verts0, verts1, htype, dtype0,
						dtype1, nframes, nviews, lc_n, GPUID, nsets, nf,
						max_frames, thread1, thread2, gpu0_stream, gpu1_stream))
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			u=ulim;
			fu=(*func)(u, verts0, verts1, htype, dtype0, dtype1, nframes, nviews, lc_n, GPUID,
					nsets, nf, max_frames, thread1, thread2, gpu0_stream, gpu1_stream);
		} else {
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u, verts0, verts1, htype, dtype0, dtype1, nframes, nviews, lc_n, GPUID,
					nsets, nf, max_frames, thread1, thread2, gpu0_stream, gpu1_stream);
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
