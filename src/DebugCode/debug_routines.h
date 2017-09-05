/* Cuda (and host!) debug functions */
__host__ void dbg_print_deldop_fit_host(struct dat_t *ddat, int s, int f, char *filename);
__host__ void dbg_print_fit_host(struct dat_t *ddat, int s, int f, char *filename);
__host__ void dbg_print_fit(struct dat_t *ddat, int s, int f, const char *filename, int gpuid);
__host__ void dbg_print_deldop_fit(struct dat_t *ddat, int s, int , char *filename);
__host__ void dbg_print_array(float *data, int x, int y);
__host__ void dbg_print_array1D(float *data, int size);
__host__ void dbg_print_array1(float *in, int size);
__host__ void dbg_print_array1D_dbl(double *data, int size, int offset, char *filename);
__host__ void dbg_print_lghtcrv_arrays(struct dat_t *ddat, int set, int n, char *filename);
__host__ void dbg_print_lghtcrv_arrays_host(struct lghtcrv_t *lghtcrv, int set, int n, char *filename);
__host__ void dbg_print_lghtcrv_xyy2(struct dat_t *ddat, int set, int ncalc, char *filename);
__host__ void dbg_print_lghtcrv_xyy2_host(struct lghtcrv_t *lghtcrv, int set, int ncalc, char *filename);
__host__ void dbg_print_lghtcrv_pos_arrays(struct dat_t *ddat, int set, int f, int npixels, int n);
__host__ void dbg_print_lghtcrv_pos_arrays_host(struct lghtcrv_t *lghtcrv, int f, int set);
__host__ void dbg_print_pos_arrays2(struct pos_t **pos, int f, int npixels, int n);
__host__ void dbg_print_pos_arrays2_host(struct pos_t *pos);
__host__ void dbg_print_pos_arrays_full(struct pos_t **pos, int f, int f_real, int npixels, int n);
__host__ void dbg_print_pos_arrays_full_host(struct pos_t *pos);
__host__ void dbg_print_facet_normals_host(struct mod_t *mod, char *fn);
__host__ void dbg_print_facet_normals(struct mod_t *dmod, int nf, char *fn);
__host__ void dbg_print_pos_bd(struct pos_t **pos, int f, int npixels, int n);
__host__ void dbg_print_lc_fit(struct dat_t *ddat, int s, char *filename_fit, int n);
__host__ void dbg_print_lc_fit_host(struct lghtcrv_t *lghtcrv, char *filename_fit, int n);
__host__ void dbg_print_pos_z_host(struct pos_t *pos, const char *fn);
