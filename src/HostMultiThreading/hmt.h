/* Flags */
extern int HMT;			/* Use Host multi-threading							*/
extern int HMT_threads;	/* Number of host threads 	*/

/* Functions for host multi-threading : */
double bestfit_hmt(struct par_t *par, struct mod_t *mod, struct dat_t *dat);
void vary_params_hmt(struct par_t *par, struct mod_t *mod, struct dat_t *dat,
		double *deldop_zmax, double *rad_xsec, double *opt_brightness, double
		*cos_subradarlat, pthread_t *hmt_thread);
