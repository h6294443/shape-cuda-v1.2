/* Flags */
extern int HMT;			/* Use Host multi-threading							*/
extern int HMT_threads;	/* Number of host threads 	*/

/* Functions for host multi-threading : */
double bestfit_hmt(struct par_t *par, struct mod_t *mod, struct dat_t *dat);
double brent_abs_hmt(double ax, double bx, double cx, double (*f) (double,
		pthread_t*), double tol, double abstol, double *xmin, pthread_t
		*hmt_thread);
void calc_fits_hmt( struct par_t *par, struct mod_t *mod, struct dat_t *dat,
		pthread_t *hmt_thread);
void mnbrak_hmt(double *ax, double *bx, double *cx, double *fa, double *fb,
	    double *fc, double (*func)(double, pthread_t*), pthread_t *hmt_stream);
void vary_params_hmt(struct par_t *par, struct mod_t *mod, struct dat_t *dat,
		double *deldop_zmax, double *rad_xsec, double *opt_brightness, double
		*cos_subradarlat, pthread_t *hmt_thread);
