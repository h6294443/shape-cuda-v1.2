/*****************************************************************************************
                                                                                bestfit.c

Iterate over all floating parameters, at each step adjusting just one parameter x in order
to minimize objective(x), the objective function (reduced chi-square plus penalties).
Continue until the fractional reduction in objective(x) due to a full pass through the
parameter list is less than term_prec.  Return the final value of the objective function.
__________________________________________________________________________________________
Modified 2016 July 7 by Matt Engels:
	Adapted for use in shape-cuda.
------------------------------------------------------------------------------------------

 *****************************************************************************************/
extern "C" {
#include "../shape/head.h"
}
static __device__ double *hotparam_mgpu;
static __device__ int dmax_frames_mgpu;
static struct par_t *spar, *sdev_par;
static struct mod_t *smod, *sdev_mod;
static struct dat_t *sdat, *sdev_dat;

static int newsize_mgpu, newshape_mgpu, newspin_mgpu, newphoto_mgpu,
	newdelcor_mgpu, newdopscale_mgpu, newxyoff_mgpu, showvals_mgpu=0,
	vary_delcor0_size_mgpu, vary_delcor0_shapespin_mgpu,
	vary_dopscale_spin_mgpu, vary_dopscale_sizeshape_mgpu, vary_alb_size_mgpu,
	vary_alb_shapespin_mgpu, vary_hapke_mgpu, call_vary_params_mgpu,
	check_posbnd_mgpu, check_badposet_mgpu, check_badradar_mgpu;
static double deldop_zmax_mgpu, deldop_zmax_mgpu_save, cos_subradarlat_mgpu,
	cos_subradarlat_mgpu_save, rad_xsec_mgpu, rad_xsec_mgpu_save,
	opt_brightness_mgpu, opt_brightness_mgpu_save, baddiam_factor_mgpu,
	badphoto_factor_mgpu, posbnd_factor_mgpu, badposet_factor_mgpu,
	badradar_factor_mgpu, baddopscale_factor_mgpu;
static unsigned char type;

static double hotparamval;

__host__ double objective_mgpu(double x, unsigned char *htype, unsigned char
		*dtype, int *nframes, int *nviews, int *lc_n, int nsets, int nf,
		cudaStream_t *gpu0_stream, cudaStream_t *gpu1_stream);
__device__ double mgpu_hotparamval, mgpu_dummyval=0.0;
__device__ int mgpu_partype;

__global__ void mgpu_get_flags_krnl(struct par_t *dpar, unsigned char *flags) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		flags[0] = dpar->baddiam;
		flags[1] = dpar->badphoto;
		flags[2] = dpar->posbnd;
		flags[3] = dpar->badposet;
		flags[4] = dpar->badradar;
		flags[5] = dpar->baddopscale;
	}
}
__global__ void mgpu_ocs_get_flags_krnl(struct par_t *dpar, unsigned char *flags,
		double *dlogfactors) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		flags[0] = dpar->baddiam;
		flags[1] = dpar->badphoto;
		flags[2] = dpar->posbnd;
		flags[3] = dpar->badposet;
		flags[4] = dpar->badradar;
		flags[5] = dpar->baddopscale;

		dlogfactors[0] = dpar->bad_objfactor;
		dlogfactors[1] = dpar->baddiam_logfactor;
		dlogfactors[2] = dpar->badphoto_logfactor;
		dlogfactors[3] = dpar->posbnd_logfactor;
		dlogfactors[4] = dpar->badposet_logfactor;
		dlogfactors[5] = dpar->badradar_logfactor;
		dlogfactors[6] = dpar->baddopscale_logfactor;
	}
}
__global__ void mgpu_set_hotparam_initial_krnl(int max_frames) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		hotparam_mgpu = &mgpu_dummyval;
		dmax_frames_mgpu = max_frames;
	}
}
__global__ void mgpu_set_hotparam_pntr_krnl(double **fpntr,
		int *fpartype, int p) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		hotparam_mgpu = fpntr[p];	/* This is pointing at a device variable */
		mgpu_partype = fpartype[p];  /* parameter type */
	}
}
__global__ void mgpu_get_hotparam_val_krnl() {
	/* Single threaded kernel */
	if (threadIdx.x == 0)
		mgpu_hotparamval = *hotparam_mgpu;
}
__global__ void mgpu_mult_hotparam_val_krnl(double factor) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0)
		*hotparam_mgpu *= factor;
}
__global__ void mgpu_set_hotparam_val_krnl(double newvalue) {
	/* Single-threaded kernel */
	if (threadIdx.x == 0) {
		*hotparam_mgpu = newvalue;
		mgpu_hotparamval = newvalue;
	}
}

__host__ double bestfit_mgpu(struct par_t *dpar, struct mod_t *dmod,
		struct dat_t *ddat, struct par_t *par, struct mod_t *mod,
		struct dat_t *dat)
{
	char hostname[MAXLEN], dofstring[MAXLEN];
	int i, iter=0, p, cntr, first_fitpar, partype, keep_iterating=1, ilaw, nf, term_maxiter;
	long pid_long;
	pid_t pid;
	double beginerr, enderr, ax, bx, cx, obja, objb, objc, xmin, final_chi2,
		final_redchi2, dummyval2, dummyval3, dummyval4, delta_delcor0,
		dopscale_factor, radalb_factor, optalb_factor, *hfparstep, *hfpartol,
		*hfparabstol, objfunc_start, term_prec;
	unsigned char *flags, *hflags, *htype, *dtype, action, avoid_badpos, term_badmodel;
	int nsets, *nframes, *lc_n, *nviews, nfpar, *hfpartype, npar_update, max_frames=0,
			max_streams=0;
	struct vertices_t **verts;
	dim3 THD, BLK;

	/* This section collects parameters used for CUDA kernel launches throughout
	 * the program.  The cudaStreams created here are used/re-used for the
	 * lifetime of one program run */
	nsets = dat->nsets;
	nfpar = par->nfpar;
	nf = mod->shape.comp[0].real.nf;
	action = par->action;
	npar_update = par->npar_update;
	avoid_badpos = par->avoid_badpos;
	objfunc_start = par->objfunc_start;
	term_prec = par->term_prec;
	term_badmodel = par->term_badmodel;
	type = mod->shape.comp[0].type;
	htype 	= (unsigned char *) malloc(nsets*sizeof(unsigned char));
	nframes = (int *) malloc(nsets*sizeof(int));
	lc_n	= (int *) malloc(nsets*sizeof(int));
	nviews 	= (int *) malloc(nsets*sizeof(int));
	gpuErrchk(cudaMalloc((void**)&dtype, sizeof(unsigned char)*nsets));
	gpuErrchk(cudaMalloc((void**)&verts, sizeof(struct vertices_t*)*2));

	for (int s=0; s<nsets; s++) {
		htype[s] = dat->set[s].type;
		switch (htype[s]) {
		case DELAY:
			nframes[s] = dat->set[s].desc.deldop.nframes;
			nviews[s]  = dat->set[s].desc.deldop.nviews;
			lc_n[s]    = 0;
			break;
		case DOPPLER:
			nframes[s] = dat->set[s].desc.doppler.nframes;
			nviews[s]  = dat->set[s].desc.doppler.nviews;
			lc_n[s]    = 0;
			break;
		case POS:
			nframes[s] = dat->set[s].desc.poset.nframes;
			nviews[s]  = dat->set[s].desc.poset.nviews;
			lc_n[s]    = 0;
			break;
		case LGHTCRV:
			nframes[s] = dat->set[s].desc.lghtcrv.ncalc;
			nviews[s]  = dat->set[s].desc.lghtcrv.nviews;
			lc_n[s]    = dat->set[s].desc.lghtcrv.n;
			break;
		}
		if (nframes[s]>max_frames)	max_frames = nframes[s];
	}
	gpuErrchk(cudaMemcpy(dtype, htype, sizeof(unsigned char)*nsets,
			cudaMemcpyHostToDevice));

	/* The following check is necessary to ensure the dVdIdCOM reduction has
	 * enough streams to operate properly	 */
	if (max_frames < 13) max_streams = 13;
	else	max_streams = max_frames;

	/* Create streams for gpu0 (the only gpu in single-GPU mode) */
	gpuErrchk(cudaSetDevice(GPU0));
	cudaStream_t gpu0_stream[max_streams];
	for (int f=0; f<max_streams; f++)
		gpuErrchk(cudaStreamCreate(&gpu0_stream[f]));
	cudaStream_t gpu1_stream[max_streams];
	/* Check for multi-GPU flag, then change device and create gpu1 streams */
	gpuErrchk(cudaSetDevice(GPU1));
	for (int f=0; f<max_streams; f++)
		gpuErrchk(cudaStreamCreate(&gpu1_stream[f]));
	gpuErrchk(cudaSetDevice(GPU0));	/* Back to our default device */

	/*..........................End section..................................*/

	/* Get the hostname of host machine and the PID */
	(void) gethostname(hostname, MAXLEN-1);
	pid = getpid();
	pid_long = (long) pid;  /* Assumes pid_t fits in a long */
	printf("#\n# CUDA fit (pid %ld on %s)\n", pid_long, hostname);
	fflush(stdout);

	/* Allocate memory for pointers, steps, and tolerances on bothhost and
	 * device. fpntr remains a cudaMallocManaged allocation because it is a
	 * double pointer.  */
	gpuErrchk(cudaMalloc((void**)&sdev_par, sizeof(struct par_t)));
	gpuErrchk(cudaMemcpy(sdev_par, &par, sizeof(struct par_t), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMalloc((void**)&sdev_mod, sizeof(struct mod_t)));
	gpuErrchk(cudaMemcpy(sdev_mod, &mod, sizeof(struct mod_t), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMalloc((void**)&sdev_dat, sizeof(struct dat_t)));
	gpuErrchk(cudaMemcpy(sdev_dat, &dat, sizeof(struct dat_t), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMalloc((void**)&flags, sizeof(unsigned char) * 7));
	gpuErrchk(cudaMalloc((void**)&fparstep,   sizeof(double)  * nfpar));
	gpuErrchk(cudaMalloc((void**)&fpartol,    sizeof(double)  * nfpar));
	gpuErrchk(cudaMalloc((void**)&fparabstol, sizeof(double)  * nfpar));
	gpuErrchk(cudaMalloc((void**)&fpartype,   sizeof(int) 	  * nfpar));
	cudaCalloc1((void**)&fpntr,	  sizeof(double*), nfpar);
	hfparstep 	 = (double *) malloc(nfpar*sizeof(double));
	hfpartol	 = (double *) malloc(nfpar*sizeof(double));
	hfparabstol  = (double *) malloc(nfpar*sizeof(double));
	hfpartype 	 = (int *) 	  malloc(nfpar*sizeof(int));
	hflags 		 = (unsigned char *) malloc(7*sizeof(unsigned char));

	for (i=0; i<nfpar; i++)
		gpuErrchk(cudaMalloc((void**)&fpntr[i], sizeof(double) * 1));

	/* Set vertices shortcut and also set max_frames (the maximum number of
	 * frames for any one set) to device so that objective_mgpu
	 * can retrieve it later */
//	gpuErrchk(cudaDeviceSynchronize());

//	set_verts_shortcut_krnl<<<1,1>>>(dmod, verts, max_frames);
//	checkErrorAfterKernelLaunch("set_verts_shortcut_krnl");

	/* Initialize static global pointers used by objective(x) below
      to be compatible with "Numerical Recipes in C" routines       */
	spar = par;			smod = mod;			sdat = dat;
	sdev_par = dpar;	sdev_mod = dmod;	sdev_dat = ddat;

	/*  Initialize static global parameters  */
	newsize_mgpu = newshape_mgpu = newspin_mgpu = newphoto_mgpu = newdelcor_mgpu = newdopscale_mgpu = newxyoff_mgpu = 1;
	deldop_zmax_mgpu = deldop_zmax_mgpu_save = 0.0;
	cos_subradarlat_mgpu = cos_subradarlat_mgpu_save = 0.0;
	rad_xsec_mgpu = rad_xsec_mgpu_save = 0.0;
	opt_brightness_mgpu = opt_brightness_mgpu_save = 0.0;
	vary_delcor0_size_mgpu = (par->vary_delcor0 != VARY_NONE);
	vary_delcor0_shapespin_mgpu = (par->vary_delcor0 == VARY_ALL);
	vary_dopscale_spin_mgpu = (par->vary_dopscale != VARY_NONE);
	vary_dopscale_sizeshape_mgpu = (par->vary_dopscale == VARY_ALL);
	vary_alb_size_mgpu = (par->vary_radalb != VARY_NONE || par->vary_optalb != VARY_NONE);
	vary_alb_shapespin_mgpu = (par->vary_radalb == VARY_ALL || par->vary_optalb == VARY_ALL);
	vary_hapke_mgpu = 0;
	if (par->vary_optalb != VARY_NONE)
		for (ilaw=0; ilaw<mod->photo.noptlaws; ilaw++)
			if (mod->photo.opttype[ilaw] == HAPKE || mod->photo.opttype[ilaw] == HARMHAPKE
					|| mod->photo.opttype[ilaw] == INHOHAPKE)
				vary_hapke_mgpu = 1;
	call_vary_params_mgpu = (par->vary_delcor0 != VARY_NONE || par->vary_dopscale != VARY_NONE
			|| par->vary_radalb != VARY_NONE
			|| par->vary_optalb != VARY_NONE);

	/*  Initialize local parameters  */
	delta_delcor0 = 0.0;
	dopscale_factor = radalb_factor = optalb_factor = 1.0;

//	/* Do-nothing kernel */
//	mgpu_set_hotparam_val_krnl<<<1,1>>>(0);
	/* The following call sets up the parameter lists allocated above and copy
	 * the device contents to host copies */
	mkparlist_gpu(dpar, dmod,	ddat, fparstep, fpartol, fparabstol, fpartype,
			fpntr, nfpar, nsets);
	gpuErrchk(cudaMemcpy(hfparstep, fparstep, sizeof(double)*nfpar, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(hfpartol, fpartol, sizeof(double)*nfpar, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(hfparabstol, fparabstol, sizeof(double)*nfpar, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(hfpartype,	fpartype, sizeof(int)*nfpar, cudaMemcpyDeviceToHost));

	/* Compute deldop_zmax_mgpu_save, cos_subradarlat_mgpu_save, rad_xsec_mgpu_save, and
	 * opt_brightness_mgpu_save for the initial model  */
//	call_vary_params_mgpu=1;
	if (call_vary_params_mgpu)
	{
		realize_mod_gpu(dpar, dmod, type, nf, gpu0_stream);

		realize_spin_mgpu(dpar, dmod, ddat, htype, nframes, nviews,
					nsets, gpu0_stream, gpu1_stream);

		realize_photo_gpu(dpar, dmod, 1.0, 1.0, 0, nf);  /* set R_save to R */

		vary_params_mgpu(dpar, dmod, ddat, action, &deldop_zmax_mgpu_save,
					&rad_xsec_mgpu_save, &opt_brightness_mgpu_save, &cos_subradarlat_mgpu_save,
					nframes, lc_n, nviews, htype, dtype, nf, nsets,
					gpu0_stream, gpu1_stream, max_frames);
	}
	printf("rad_xsec_mgpu: %f\n", rad_xsec_mgpu_save);
	printf("deldop_zmax_mgpu: %f\n", (float)deldop_zmax_mgpu_save);

	/* Point hotparam to a dummy variable (dummyval) rather than to a model pa-
	 * rameter; then call objective(0.0) to set dummy variable = 0.0, realize
	 * the initial model, calculate the fits, return initial model's objective
	 * function as enderr.                          */
	mgpu_set_hotparam_initial_krnl<<<1,1>>>(max_frames);
	checkErrorAfterKernelLaunch("mgpu_set_hotparam_initial_krnl");

	enderr = objective_mgpu(0.0, htype, dtype, nframes,
				nviews, lc_n, nsets, nf, gpu0_stream, gpu1_stream);

	printf("#\n# searching for best fit ...\n");
	printf("%4d %8.6f to begin", 0, enderr);

//	/* Launch single-thread kernel to retrieve flags in dev_par */
//	/*		flags[0] = dpar->baddiam;
//			flags[1] = dpar->badphoto;
//			flags[2] = dpar->posbnd;
//			flags[3] = dpar->badposet;
//			flags[4] = dpar->badradar;
//			flags[5] = dpar->baddopscale;*/
//
//	mgpu_get_flags_krnl<<<1,1>>>(dpar, flags);
//	checkErrorAfterKernelLaunch("mgpu_get_flags_krnl");
//	gpuErrchk(cudaMemcpy(hflags, flags, sizeof(unsigned char)*7,
//			cudaMemcpyDeviceToHost));
//
//	/* Now act on the flags just retrieved from dev_par */
//	if (hflags[0])		printf("  (BAD DIAMS)");
//	if (hflags[1])		printf("  (BAD PHOTO)");
//	if (hflags[2])		printf("  (BAD POS)");
//	if (hflags[3])		printf("  (BAD POSET)");
//	if (hflags[4])		printf("  (BAD RADAR)");
//	if (hflags[5])		printf("  (BAD DOPSCALE)");		printf("\n");
//	fflush(stdout);
//
//	/* Display the region within each delay-Doppler or Doppler frame that, ac-
//	 * cording to initial model, has nonzero power. A warning is displayed if
//	 * any region extends beyond the data limits: the vignetting is too tight,
//	 * or else some model parameter (such as a delay correction polynomial co-
//	 * efficient) is seriously in error.   */
//	show_deldoplim_gpu(ddat, htype, nsets, nframes, max_frames);
//
//	/* Set the starting fit parameter for the first iteration only  */
//	first_fitpar = par->first_fitpar;
//	term_maxiter = par->term_maxiter;
//	if (first_fitpar < 0 || first_fitpar >= nfpar) {
//		printf("ERROR: need 0 <= first_fitpar < nparams (%d)\n", nfpar);
//		bailout("bestfit.c\n");
//	}
//
//	/* Iteratively adjust model; for each iteration, step through all free pa-
//	 * rameters, adjusting one parameter at a time so as to minimize the objec-
//	 * tive function at each step. Stop when fractional decrease in the objec-
//	 * tive function from one iteration to the next is less than term_prec.   */
//
//	do {
//		showvals_mgpu = 1;        /* show reduced chi-square and penalties at beginning */
//		beginerr = enderr;
//		printf("# iteration %d %f", ++iter, beginerr);
//
//		/* Launch single-thread kernel to retrieve flags in dev_par */
//		mgpu_get_flags_krnl<<<1,1>>>(dpar, flags);
//		checkErrorAfterKernelLaunch("mgpu_get_flags_krnl");
//		gpuErrchk(cudaMemcpy(hflags, flags, sizeof(unsigned char)*7,
//				cudaMemcpyDeviceToHost));
//
//		/* Now act on the flags just retrieved from dev_par */
//		if (hflags[0])		printf("  (BAD DIAMS)");
//		if (hflags[1])		printf("  (BAD PHOTO)");
//		if (hflags[2])		printf("  (BAD POS)");
//		if (hflags[3])		printf("  (BAD POSET)");
//		if (hflags[4])		printf("  (BAD RADAR)");
//		if (hflags[5])		printf("  (BAD DOPSCALE)");		printf("\n");
//		fflush(stdout);
//
//		/* Show breakdown of chi-square by data type    */
//		chi2_gpu(dpar, ddat, htype, dtype, nframes, lc_n, 1,
//				nsets, gpu0_stream, max_frames);
//
//		/*  Loop through the free parameters  */
//		cntr = first_fitpar % npar_update;
//		//p = first_fitpar = 1;
//		for (p=first_fitpar; p<nfpar; p++) {
//
////		p = first_fitpar;
//			/*  Adjust only parameter p on this try  */
//			mgpu_set_hotparam_pntr_krnl<<<1,1>>>(fpntr, fpartype, p);
//			checkErrorAfterKernelLaunch("mgpu_set_hotparam_pntr_krnl");
//			gpuErrchk(cudaMemcpyFromSymbol(&partype, mgpu_partype, sizeof(int),
//					0, cudaMemcpyDeviceToHost));
//
//			newsize_mgpu = newshape_mgpu = newspin_mgpu = newphoto_mgpu = newdelcor_mgpu = newdopscale_mgpu
//					= newxyoff_mgpu = 0;
//			if 		(partype == SIZEPAR)		newsize_mgpu	 	= 1;
//			else if (partype == SHAPEPAR)		newshape_mgpu 	= 1;
//			else if (partype == SPINPAR)		newspin_mgpu 	= 1;
//			else if (partype == PHOTOPAR)		newphoto_mgpu 	= 1;
//			else if (partype == DELCORPAR)		newdelcor_mgpu 	= 1;
//			else if (partype == DOPSCALEPAR)	newdopscale_mgpu	= 1;
//			else if (partype == XYOFFPAR)		newxyoff_mgpu 	= 1;
//
//			/* If this is a size parameter AND model extends beyond POS frame
//			 * AND the "avoid_badpos" parameter is turned on, shrink model by
//			 * 5% at a time until it fits within the POS frame.
//			 * We must start with the redundant model evaluation for the un-
//			 * changed value of the size parameter, in case the first call to
//			 * objective displays reduced chi-square and the penalty functions.  */
//			if (avoid_badpos && partype == SIZEPAR) {
//				mgpu_get_flags_krnl<<<1,1>>>(dpar, flags);
//				checkErrorAfterKernelLaunch("mgpu_get_flags_krnl");
//				gpuErrchk(cudaMemcpy(hflags, flags, sizeof(unsigned char)*7,
//						cudaMemcpyDeviceToHost));
//
//				/* Get value of (*hotparam_mgpu) */
//				mgpu_get_hotparam_val_krnl<<<1,1>>>();
//				checkErrorAfterKernelLaunch("mgpu_get_hotparam_val_krnl");
//				gpuErrchk(cudaMemcpyFromSymbol(&hotparamval, mgpu_hotparamval,
//						sizeof(double),	0, cudaMemcpyDeviceToHost));
//
//				while (hflags[2]) {
//					objective_mgpu(hotparamval, verts, htype, dtype,
//							nframes, nviews, lc_n, nsets, nf, gpu0_stream);
//
//					mgpu_get_flags_krnl<<<1,1>>>(dpar, flags);
//					checkErrorAfterKernelLaunch("mgpu_get_flags_krnl");
//					gpuErrchk(cudaMemcpy(hflags, flags, sizeof(unsigned char)*7,
//							cudaMemcpyDeviceToHost));
//
//					if (hflags[2]) {
//						/* Set the value pointed to by hotparam_mgpu to 0.95 of its
//						 * previous value */
//						bf_mult_hotparam_mgpu_val_krnl<<<1,1>>>(0.95);
//						checkErrorAfterKernelLaunch("mgpu_mult_hotparam_val_krnl");
//					}
//				}
//			}
//
//			/* Get value of (*hotparam_mgpu) so that mnbrak can use it*/
//			mgpu_get_hotparam_val_krnl<<<1,1>>>();
//			checkErrorAfterKernelLaunch("mgpu_get_hotparam_val_krnl");
//			gpuErrchk(cudaMemcpyFromSymbol(&hotparamval, mgpu_hotparamval,
//					sizeof(double),	0, cudaMemcpyDeviceToHost));
//
//			/* Use Numerical Recipes routine mnbrak to bracket a minimum in the
//			 * objective function (reduced chi-square plus penalties) objec-
//			 * tive(x), where x is the value of parameter p.  As initial trial
//			 * parameter values, use ax (unadjusted value) and bx, that value
//			 * incremented by the appropriate step size (length_step,spin_step,
//			 * etc.). mnbrak returns 3 parameter values, with bx between ax
//			 * and cx; note that ax and bx are changed from their input values.
//			 * It also returns the 3 corresponding objective(x) values, where
//			 * objb is less than obja and objc.  Hence there is at least one
//			 * local minimum (but not necessarily *any* global minimum)
//			 * somewhere between ax and cx.          */
//			ax = hotparamval;
//			bx = ax + hfparstep[p]; /* par usage us fine here */
//
//			mnbrak_gpu(&ax, &bx, &cx, &obja, &objb, &objc,
//					objective_mgpu, verts, htype, dtype, nframes,
//					nviews, lc_n, nsets, nf, gpu0_stream);
//
//			/* Before homing in on local minimum, initialize flags that will
//			 * tell us if model extended beyond POS frame (sky rendering) for
//			 * any trial parameter value(s), if it extended beyond any POS ima-
//			 * ges, and if it was too wide in delay-Doppler space         */
//			check_posbnd_mgpu = 0;
//			check_badposet_mgpu = 0;
//			check_badradar_mgpu = 0;
//
//			/* Now use Numerical Recipes function brent to find local minimum -
//			 * that is, to find xmin, the best value of x, to within the
//			 * *fractional* tolerance specified for parameter p (length_tol,
//			 * spin_tol, etc.). brent's return value is the minimized objective
//			 * function, objective(xmin). If more than one local minimum bet-
//			 * ween ax and cx, brent might not find the best one. brent_abs is
//			 * a modified version of brent that has an absolute fitting tole-
//			 * rance as one of its arguments, in addition to the existing
//			 * fractional tolerance.                                      */
//			enderr = brent_abs_gpu(ax, bx, cx, objective_mgpu, hfpartol[p],
//					hfparabstol[p], &xmin, verts, htype, dtype, nframes, nviews, lc_n,
//					nsets, nf, gpu0_stream);
//
//			/* Realize whichever part(s) of the model has changed.
//			 *
//			 * The code here is somewhat opaque because more than one part of
//			 * the model may have changed - if the "vary_delcor0" "vary_radalb"
//			 * and/or "vary_optalb" parameter is being used to permit joint pa-
//			 * rameter adjustments. Before calling the vary_params routine, the
//			 * size/shape and spin states must be realized (realize_mod and
//			 * realize_spin); if albedos are being varied jointly with other
//			 * parameters, the photometric state must also be realized
//			 * (realize_photo); and in either case the 0th-order delay correc-
//			 * tion polynomial coefficients must be reset to their saved
//			 * values via the appropriate call to realize_delcor.          */
//			/* Set the value pointed to by hotparam_mgpu to 0.95 of its
//			 * previous value (*hotparam_mgpu) = xmin; */
//			mgpu_set_hotparam_val_krnl<<<1,1>>>(xmin);
//			checkErrorAfterKernelLaunch("mgpu_set_hotparam_val_krnl");
//			gpuErrchk(cudaMemcpyFromSymbol(&hotparamval, mgpu_hotparamval,
//					sizeof(double),	0, cudaMemcpyDeviceToHost));
//
//			if (newsize_mgpu || newshape_mgpu)
//				realize_mod_gpu(dpar, dmod, type, nf, gpu0_stream);
//			if (newspin_mgpu) {
//				realize_spin_gpu(dpar, dmod, ddat, htype, nframes,
//						nviews, nsets, gpu0_stream);
//			}
//			if ((newsize_mgpu && vary_alb_size_mgpu) || ((newshape_mgpu ||
//					newspin_mgpu) && vary_alb_shapespin_mgpu))
//				realize_photo_gpu(dpar, dmod, 1.0, 1.0, 1);  /* set R to R_save */
//			if ((newsize_mgpu && vary_delcor0_size_mgpu) || ((newshape_mgpu || newspin_mgpu)
//					&& vary_delcor0_shapespin_mgpu)) {
//				realize_delcor_gpu(ddat, 0.0, 1, nsets, nframes);  /* set delcor0 to delcor0_save */
//			}
//			if ((newspin_mgpu && vary_dopscale_spin_mgpu) || ((newsize_mgpu || newshape_mgpu)
//					&& vary_dopscale_sizeshape_mgpu))
//				realize_dopscale_gpu(dpar, ddat, 1.0, 1, nsets, dtype);  /* set dopscale to dopscale_save */
//			if (call_vary_params_mgpu) {
//				/* Call vary_params to get the adjustments to 0th-order delay
//				 * correction polynomial coefficients, to Doppler scaling fac-
//				 * tors, and to radar and optical albedos                  */
//
//				vary_params_gpu(dpar,dmod,ddat,11,&deldop_zmax_mgpu,
//						&rad_xsec_mgpu, &opt_brightness_mgpu, &cos_subradarlat_mgpu,
//						nframes, lc_n, nviews, verts, htype, dtype, nf, nsets,
//						gpu0_stream, max_frames);
//
//				delta_delcor0 = (deldop_zmax_mgpu - deldop_zmax_mgpu_save)*KM2US;
//				if (cos_subradarlat_mgpu != 0.0)
//					dopscale_factor = cos_subradarlat_mgpu_save/cos_subradarlat_mgpu;
//				if (rad_xsec_mgpu != 0.0)
//					radalb_factor = rad_xsec_mgpu_save/rad_xsec_mgpu;
//				if (opt_brightness_mgpu != 0.0)
//					optalb_factor = opt_brightness_mgpu_save/opt_brightness_mgpu;
//			}
//			if ((newsize_mgpu && vary_alb_size_mgpu) || ((newshape_mgpu || newspin_mgpu) &&
//					vary_alb_shapespin_mgpu)) {
//				realize_photo_gpu(dpar, dmod, radalb_factor, optalb_factor, 2);  /* reset R, then R_save */
//
//				/* Must update opt_brightness_mgpu_save for Hapke optical scattering
//				 * law, since single-scattering albedo w isn't just an overall
//				 * scaling factor  */
//				if (vary_hapke_mgpu) {
//					vary_params_gpu(dpar,dmod,ddat,12,&dummyval2,
//							&dummyval3,&opt_brightness_mgpu,&dummyval4,
//							nframes, lc_n, nviews, verts, htype, dtype, nf, nsets,
//							gpu0_stream, max_frames);
//				}
//			} else if (newphoto_mgpu) {
//				rad_xsec_mgpu_save = rad_xsec_mgpu;
//				opt_brightness_mgpu_save = opt_brightness_mgpu;
//				realize_photo_gpu(dpar, dmod, 1.0, 1.0, 0);  /* set R_save to R */
//			}
//			if ((newsize_mgpu && vary_delcor0_size_mgpu) || ((newshape_mgpu || newspin_mgpu) &&
//					vary_delcor0_shapespin_mgpu)) {
//				deldop_zmax_mgpu_save = deldop_zmax_mgpu;
//				realize_delcor_gpu(ddat, delta_delcor0, 2, nsets, nframes);  /* reset delcor0, then delcor0_save */
//			} else if (newdelcor_mgpu)
//				realize_delcor_gpu(ddat, 0.0, 0, nsets, nframes);  /* set delcor0_save to delcor0 */
//
//			if ((newspin_mgpu && vary_dopscale_spin_mgpu) || ((newsize_mgpu || newshape_mgpu) &&
//					vary_dopscale_sizeshape_mgpu)) {
//				cos_subradarlat_mgpu_save = cos_subradarlat_mgpu;
//				realize_dopscale_gpu(dpar, ddat, dopscale_factor, 2, nsets, dtype);  /* reset dopscale, then dopscale_save */
//			} else if (newdopscale_mgpu) {
//				realize_dopscale_gpu(dpar, ddat, 1.0, 0, nsets, dtype);  /* set dopscale_save to dopscale */
//			}
//			if (newxyoff_mgpu)
//				realize_xyoff_gpu(ddat, nsets, dtype);
//
//			/* If the model extended beyond POS frame (sky rendering) for any
//			 * trial parameter value(s), if it extended beyond any plane-of-
//			 * sky fit frames, or if it was too wide in delay-Doppler space,
//			 * evaluate model for best-fit parameter value to check if these
//			 * problems persist - that is, to update "posbnd" "badposet" and
//			 * "badradar" parameters for updated model.
//			 * (This needn't be done for "baddiam" "badphoto" flags: if we've
//			 * just finished adjusting an ellipsoid dimension or photometric
//			 * parameter, realize_mod or realize_photo was called in code block
//			 * above in order to realize the changed portion of model, and that
//			 * call updated corresponding flag. Also we needn't worry about the
//			 * "baddopscale" flag, since realize_dopscale was called above if
//			 * Doppler scaling factors were changed.) The call to objective
//			 * (*hotparam_mgpu) first sets *hotparam_mgpu (the parameter that we just
//			 * adjusted) equal to itself (i.e., no change) and then calls
//			 * calc_fits to evaluate the model for all datasets.          */
//			if (check_posbnd_mgpu || check_badposet_mgpu || check_badradar_mgpu)
//				objective_mgpu(hotparamval, verts, htype, dtype,
//						nframes, nviews, lc_n, nsets, nf, gpu0_stream);
//
//			/* Launch single-thread kernel to retrieve flags in dev_par */
//			mgpu_get_flags_krnl<<<1,1>>>(dpar, flags);
//			checkErrorAfterKernelLaunch("mgpu_get_flags_krnl");
//			gpuErrchk(cudaMemcpy(hflags, flags, sizeof(unsigned char)*7,
//					cudaMemcpyDeviceToHost));
//			/* Display the objective function after each parameter adjustment.  */
//			printf("%4d %8.6f %d", p, enderr, iround(par->fpartype[p]));
//			if (hflags[0])		printf("  (BAD DIAMS)");
//			if (hflags[1])		printf("  (BAD PHOTO)");
//			if (hflags[2])		printf("  (BAD POS)");
//			if (hflags[3])		printf("  (BAD POSET)");
//			if (hflags[4])		printf("  (BAD RADAR)");
//			if (hflags[5])		printf("  (BAD DOPSCALE)");
//			printf("\n");
//			fflush(stdout);
//
//			/* Display reduced chi-square and individual penalty values after
//			 * every 20th parameter adjustment. Setting showvals_mgpu to 1 here
//			 * means that these things will be displayed next time objective(x)
//			 * is evaluated - at start of NEXT parameter adjustment.  Specifi-
//			 * cally, they will be displayed when routine mnbrak evaluates
//			 * objective(x) for *unadjusted* parameter value ax (see comment
//			 * above).
//			 * Also rewrite model and obs files after every 20th parameter
//			 * adjustment. Most of obs file doesn't change, but some floating
//			 * parameters (i.e. delay correction polynomial coefficients) do.  */
//			if (++cntr >= npar_update) {
//				cntr = 0;
//				showvals_mgpu = 1;
//				calc_fits_gpu(dpar, dmod, ddat, verts, nviews,
//						nframes, lc_n, htype, nsets, nf, gpu0_stream, max_frames);
//				chi2_gpu(dpar, ddat, htype, dtype, nframes,
//						lc_n, 0, nsets, gpu0_stream, max_frames);
//
//				//write_mod( par, mod);
//				//write_dat( par, dat);
//			}
//		}  // End fitpar loop
//
//		/* End of this iteration: Write model and data to disk, and display the
//		 * region within each delay-Doppler or Doppler frame for which model
//		 * power is nonzero.                                               */
//		if (cntr != 0) {
//			calc_fits_gpu(dpar, dmod, ddat, verts, nviews,
//					nframes, lc_n, htype, nsets, nf, gpu0_stream, max_frames);
//			chi2_gpu(dpar, ddat, htype, dtype, nframes,
//					lc_n, 0, nsets, gpu0_stream, max_frames);
//
//			//write_mod( par, mod);
//			//write_dat( par, dat);
//		}
//		show_deldoplim_gpu(ddat, htype, nsets, nframes, max_frames);
//
//		/* Check if we should start a new iteration  */
//		if (iter == term_maxiter) {
//			/* Just completed last iteration permitted by "term_maxiter" para-
//			 * meter, so stop iterating; note that since iter is 1-based, this
//			 * test is always false if "term_maxiter" = 0 (its default value)  */
//			keep_iterating = 0;
//
//		} else if (first_fitpar > 0) {
//			/* Just completed partial iteration (possible for iteration 1): if
//			 * "objfunc_start" parameter was given, check if fractional decrea-
//			 * se in objective function *relative to objfunc_start* during the
//			 * just-completed iteration was larger than term_prec, thus
//			 * justifying a new iteration; if it wasn't specified, definitely
//			 * proceed to a new iteration.                            */
//			if (objfunc_start > 0.0)
//				keep_iterating = ((objfunc_start - enderr)/enderr >= term_prec);
//			else
//				keep_iterating = 1;
//			first_fitpar = 0;     /* for all iterations after the first iteration */
//
//		} else if (term_badmodel && (hflags[0] || hflags[1] || hflags[2] ||
//				hflags[3] || hflags[4] || hflags[5]) ) {
//
//			/* Just completed a full iteration, stop iterating because "term_
//			 * badmodel" parameter is turned on and model has a fatal flaw: it
//			 * extends beyond POS frame OR it one or more illegal photometric
//			 * parameters OR it has one or more tiny or negative ellipsoid dia-
//			 * meters OR it has plane-of-sky fit frames too small to "contain"
//			 * model OR it is too wide in delay-Doppler space for (delay-)
//			 * Doppler fit frames to be correctly constructed OR it has out-of-
//			 * range values for one or more Doppler scaling factors    */
//			keep_iterating = 0;
//
//		} else {
//			/* Just completed a full iteration and the model has no fatal flaws
//			 * (or else the "term_badmodel" parameter is turned off): keep
//			 * iterating if fractional decrease objective function during the
//			 * just-completed iteration was greater than term_prec         */
//			keep_iterating = ((beginerr - enderr)/enderr >= term_prec);
//		}
//
//	} while (keep_iterating);
//
//		/* Show final values of reduced chi-square, individual penalty functions,
//		 * and the objective function  */
//		final_chi2 = chi2_gpu(dpar, ddat, htype, dtype, nframes,
//				lc_n, 1, nsets, gpu0_stream, max_frames);
//
//		final_redchi2 = final_chi2/dat->dof;
//		printf("# search completed\n");
//
//		/* Launch single-thread kernel to get these final flags from dev->par:
//		 * pen.n, baddiam, badphoto, posbnd, badposet, badradar, baddopscale */
//		/* Launch single-thread kernel to retrieve flags in dev_par */
//		mgpu_get_flags_krnl<<<1,1>>>(dpar, flags);
//		checkErrorAfterKernelLaunch("mgpu_get_flags_krnl");
//		gpuErrchk(cudaMemcpy(hflags, flags, sizeof(unsigned char)*7,
//				cudaMemcpyDeviceToHost));
//
//		if (par->pen.n > 0 || hflags[0] || hflags[1] || hflags[2]	|| hflags[3] ||
//				hflags[4] || hflags[5]) {
//			printf("#\n");
//			printf("# %15s %e\n", "reduced chi2", final_redchi2);
//			if (par->pen.n > 0) {
//				par->showstate = 1;
//			penalties_gpu(dpar, dmod, ddat);
//			par->showstate = 0;
//		}
//		if (hflags[0])
//			printf("# objective func multiplied by %.1f: illegal ellipsoid diameters\n",
//					baddiam_factor_mgpu);
//		if (hflags[1])
//			printf("# objective func multiplied by %.1f: illegal photometric parameters\n",
//					badphoto_factor_mgpu);
//		if (hflags[2])
//			printf("# objective func multiplied by %.1f: model extends beyond POS frame\n",
//					posbnd_factor_mgpu);
//		if (hflags[3])
//			printf("# objective func multiplied by %.1f: "
//					"model extends beyond plane-of-sky fit image\n",
//					badposet_factor_mgpu);
//		if (hflags[4])
//			printf("# objective func multiplied by %.1f: "
//					"model is too wide in delay-Doppler space to construct fit image\n",
//					badradar_factor_mgpu);
//		if (hflags[5])
//			printf("# objective func multiplied by %.1f: illegal Doppler scaling factors\n",
//					baddopscale_factor_mgpu);
//		printf("# ----------------------------\n");
//		printf("# %15s %e\n", "objective func", enderr);
//		printf("#\n");
//	}
//	intifpossible( dofstring, MAXLEN, dat->dof, SMALLVAL, "%f");
	printf("# final chi2 = %e for %s dof (reduced chi2 = %f)\n",
			final_chi2, dofstring, final_redchi2);
	printf("#\n");
	fflush(stdout);

	/* Destroy the streams */
	cudaSetDevice(GPU0);
	for (int f=0; f<max_frames; f++)
		cudaStreamDestroy(gpu0_stream[f]);



	free(hflags);
	free(htype);
	free(nframes);
	free(lc_n);
	free(nviews);
	free(hfparstep);
	free(hfpartol);
	free(hfparabstol);
//	free(fpartype);
	cudaFree(sdev_par);
	cudaFree(sdev_mod);
	cudaFree(sdev_dat);
	cudaFree(fparstep);
	cudaFree(fpartol);
	cudaFree(fparabstol);
	cudaFree(fpartype);
	cudaFree(fpntr);
	cudaFree(flags);
	cudaFree(dtype);
	cudaFree(verts);
	cudaDeviceReset();
	//cudaProfilerStop();
	return enderr;
}

/* objective_mgpu is a version of objective_cuda that takes an extra
 * argument - the cudaStreams created in bestfit_cuda2. The goal is to
 * reduce overhead from stream creation/destruction to a minimum by having
 * just one set number of streams per program run.  */
__host__ double objective_mgpu(
		double x,
		unsigned char *htype,
		unsigned char *dtype,
		int *nframes,
		int *nviews,
		int *lc_n,
		int nsets,
		int nf,
		cudaStream_t *gpu0_stream,
		cudaStream_t *gpu1_stream)
{
	double err, pens, delta_delcor0, dopscale_factor, radalb_factor,
		optalb_factor, *dlogfactors, *hlogfactors;
	unsigned char *dflags, *hflags;
	int max_frames;

	gpuErrchk(cudaSetDevice(GPU0));
	gpuErrchk(cudaMalloc((void**)&dflags, sizeof(unsigned char)*7));
	gpuErrchk(cudaMalloc((void**)&dlogfactors, sizeof(double)*7));
	hflags 	 	= (unsigned char *) malloc(7*sizeof(unsigned char));
	hlogfactors	= (double *) malloc(7*sizeof(double));

	/* Initialize local parameters  */
	delta_delcor0 = 0.0;
	dopscale_factor = radalb_factor = optalb_factor = 1.0;

	/* Assign new trial value to the model parameter being adjusted  */
	mgpu_set_hotparam_val_krnl<<<1,1>>>(x);	//(*hotparam_mgpu) = x;
	checkErrorAfterKernelLaunch("mgpu_set_hotparam_val_krnl (in objective_cuda)");
	gpuErrchk(cudaMemcpyFromSymbol(&max_frames, dmax_frames_mgpu,
			sizeof(int),	0, cudaMemcpyDeviceToHost));

	/* Realize whichever part(s) of the model have changed, then calculate root's
	 * contribution to chi-square.
	 * The code here is somewhat opaque because more than one part of the model
	 * may have changed - if the "vary_delcor0" "vary_dopscale" "vary_radalb" and
	 * /or "vary_optalb" parameter is being used to permit joint parameter ad-
	 * justments. Before calling the vary_params routine, the size/shape and spin
	 * states must be realized (realize_mod and realize_spin); if albedos are
	 * being varied jointly with other parameters, the photometric state must
	 * also be realized (realize_photo); and in either case the 0th-order delay
	 * correction polynomial coefficients and the Doppler scaling factors must be
	 * reset to their saved values via the appropriate calls to realize_delcor
	 * and realize_dopscale, respectively.*/

	if (newsize_mgpu || newshape_mgpu)
		realize_mod_gpu(sdev_par, sdev_mod, type, nf, gpu0_stream);
	if (newspin_mgpu)
		realize_spin_mgpu(sdev_par, sdev_mod, sdev_dat, htype, nframes,
				nviews, nsets, gpu0_stream, gpu1_stream);

	if ((newsize_mgpu && vary_alb_size_mgpu) || ((newshape_mgpu || newspin_mgpu) && vary_alb_shapespin_mgpu))
		realize_photo_gpu(sdev_par, sdev_mod, 1.0, 1.0, 1, nf);  /* set R to R_save */
	if ((newsize_mgpu && vary_delcor0_size_mgpu) || ((newshape_mgpu || newspin_mgpu) && vary_delcor0_shapespin_mgpu))
		realize_delcor_gpu(sdev_dat, 0.0, 1, nsets, nframes);  /* set delcor0 to delcor0_save */

	if ((newspin_mgpu && vary_dopscale_spin_mgpu) || ((newsize_mgpu || newshape_mgpu) && vary_dopscale_sizeshape_mgpu))
		realize_dopscale_gpu(sdev_par, sdev_dat, 1.0, 1, nsets, dtype);  /* set dopscale to dopscale_save */
	if (call_vary_params_mgpu) {
		/* Call vary_params to get the trial adjustments to 0th-order delay correc-
		 * tion polynomial coefficients, to Doppler scaling factors,and to radar
		 * and optical albedos, then send them to the branch nodes  */

		vary_params_mgpu(sdev_par, sdev_mod, sdev_dat, spar->action,
				&deldop_zmax_mgpu, &rad_xsec_mgpu, &opt_brightness_mgpu,
				&cos_subradarlat_mgpu, nframes, lc_n, nviews, htype, dtype, nf,
				nsets, gpu0_stream, gpu1_stream, max_frames);

		delta_delcor0 = (deldop_zmax_mgpu - deldop_zmax_mgpu_save)*KM2US;
		if (cos_subradarlat_mgpu != 0.0)
			dopscale_factor = cos_subradarlat_mgpu_save/cos_subradarlat_mgpu;
		if (rad_xsec_mgpu != 0.0)
			radalb_factor = rad_xsec_mgpu_save/rad_xsec_mgpu;
		if (opt_brightness_mgpu != 0.0)
			optalb_factor = opt_brightness_mgpu_save/opt_brightness_mgpu;
	}

	if ((newsize_mgpu && vary_alb_size_mgpu) || ((newshape_mgpu || newspin_mgpu) && vary_alb_shapespin_mgpu))
		realize_photo_gpu(sdev_par, sdev_mod, radalb_factor, optalb_factor, 1, nf);  /* adjust R */
	else if (newphoto_mgpu)
		realize_photo_gpu(sdev_par, sdev_mod, 1.0, 1.0, 0, nf);  /* set R_save to R */
	if ((newsize_mgpu && vary_delcor0_size_mgpu) || ((newshape_mgpu || newspin_mgpu) && vary_delcor0_shapespin_mgpu)) {
		realize_delcor_gpu(sdev_dat, delta_delcor0, 1, nsets, nframes);  /* adjust delcor0 */
	}
	else if (newdelcor_mgpu) {
		realize_delcor_gpu(sdev_dat, 0.0, 0, nsets, nframes);  /* set delcor0_save to delcor0 */
	}
	if ((newspin_mgpu && vary_dopscale_spin_mgpu) || ((newsize_mgpu || newshape_mgpu) && vary_dopscale_sizeshape_mgpu))
		realize_dopscale_gpu(sdev_par, sdev_dat, dopscale_factor, 1, nsets, dtype);  /* adjust dopscale */
	else if (newdopscale_mgpu)
		realize_dopscale_gpu(sdev_par, sdev_dat, 1.0, 0, nsets, dtype);  /* set dopscale_save to dopscale */
	if (newxyoff_mgpu)
		realize_xyoff_gpu(sdev_dat, nsets, dtype);

	calc_fits_mgpu(sdev_par, sdev_mod, sdev_dat, nviews, nframes, lc_n, htype,
			nsets, nf, gpu0_stream, gpu1_stream, max_frames);

//	calc_fits_gpu(sdev_par, sdev_mod, sdev_dat, verts, nviews, nframes, lc_n,
//			htype, nsets, nf, gpu0_stream, max_frames);
//	err = chi2_gpu(sdev_par, sdev_dat, htype, dtype, nframes, lc_n, 0, nsets,
//			gpu0_stream, max_frames);
//
//	/* Divide chi-square by DOF to get reduced chi-square.    */
//	err /= sdat->dof;
//
//	/* If bestfit has set showvals_mgpu = 1, display reduced chi-square. Then set
//	 * spar->showstate = 1, so that when function penalties is called later,
//	 * it "knows" that it should display the individual penalty values.
//	 * Reset showstate to 0 if showvals_mgpu = 0.  */
//	if (showvals_mgpu) {
//		printf("# %15s %e\n", "reduced chi2", err);
//		spar->showstate = 1;
//	}
//	else
//		spar->showstate = 0;
//
//	/* Compute penalties and add to reduced chi-square. Individual penalty values
//	 * will be displayed if we set spar->showstate = 1 a few lines back.        */
//	pens = penalties_gpu(sdev_par, sdev_mod, sdev_dat);
//	err += pens;
//
//	/* Double the objective function if there's an ellipsoid component with tiny
//	 * or negative diameter, if any optical photometric parameters have invalid
//	 * values, if any portion of the model lies outside specified POS window or
//	 * outside any plane-of-sky fit image, or if model is too wide in delay-Dopp-
//	 * ler space for any (delay-)Doppler fit image to be correctly constructed.
//	 * This effectively rules out any models with any of these flaws.         */
//	/* NOTE: TO-DO: baddiam may need to come from elsewhere other than spar.
//	 * However, bestfit gets called only once and spar/smod/sdat gets copied
//	 * only once.
//	 * flags[0] = dpar->baddiam;
//		flags[1] = dpar->badphoto;
//		flags[2] = dpar->posbnd;
//		flags[3] = dpar->badposet;
//		flags[4] = dpar->badradar;
//		flags[5] = dpar->baddopscale;
//
//		dlogfactors[0] = dpar->bad_objfactor;
//		dlogfactors[1] = dpar->baddiam_logfactor;
//		dlogfactors[2] = dpar->badphoto_logfactor;
//		dlogfactors[3] = dpar->posbnd_logfactor;
//		dlogfactors[4] = dpar->badposet_logfactor;
//		dlogfactors[5] = dpar->badradar_logfactor;
//		dlogfactors[6] = dpar->baddopscale_logfactor;
//	 */
//	mgpu_ocs_get_flags_krnl<<<1,1>>>(sdev_par, dflags, dlogfactors);
//	checkErrorAfterKernelLaunch("mgpu_get_flags_krnl");
//	gpuErrchk(cudaMemcpy(hflags, dflags, sizeof(unsigned char)*7,
//			cudaMemcpyDeviceToHost));
//	gpuErrchk(cudaMemcpy(hlogfactors, dlogfactors, sizeof(double)*6,
//			cudaMemcpyDeviceToHost));
//
//	if (hflags[0]) {
//		baddiam_factor_mgpu = hlogfactors[0] * exp(hlogfactors[1]);
//		err *= baddiam_factor_mgpu;
//		if (showvals_mgpu)
//			printf("# objective func multiplied by %.1f: illegal ellipsoid diameters\n",
//					baddiam_factor_mgpu);
//	}
//	if (hflags[1]) {
//		badphoto_factor_mgpu = hlogfactors[0] * exp(hlogfactors[2]);
//		err *= badphoto_factor_mgpu;
//		if (showvals_mgpu)
//			printf("# objective func multiplied by %.1f: illegal photometric parameters\n",
//					badphoto_factor_mgpu);
//	}
//	if (hflags[2]) {
//		check_posbnd_mgpu = 1;     /* tells bestfit about this problem */
//		posbnd_factor_mgpu = hlogfactors[0] * exp(hlogfactors[3]);
//		err *= posbnd_factor_mgpu;
//		if (showvals_mgpu)
//			printf("# objective func multiplied by %.1f: model extends beyond POS frame\n",
//					posbnd_factor_mgpu);
//	}
//	if (hflags[3]) {
//		check_badposet_mgpu = 1;     /* tells bestfit about this problem */
//		badposet_factor_mgpu = hlogfactors[0] * exp(hlogfactors[4]);
//		err *= badposet_factor_mgpu;
//		if (showvals_mgpu)
//			printf("# objective func multiplied by %.1f: plane-of-sky fit frame too small\n",
//					badposet_factor_mgpu);
//	}
//	if (hflags[4]) {
//		check_badradar_mgpu = 1;     /* tells bestfit about this problem */
//		badradar_factor_mgpu = hlogfactors[0] * exp(hlogfactors[5]);
//		err *= badradar_factor_mgpu;
//		if (showvals_mgpu)
//			printf("# objective func multiplied by %.1f: model too wide in delay-Doppler space\n",
//					badradar_factor_mgpu);
//	}
//	if (hflags[5]) {
//		baddopscale_factor_mgpu = hlogfactors[0] * exp(hlogfactors[6]);
//		err *= baddopscale_factor_mgpu;
//		if (showvals_mgpu)
//			printf("# objective func multiplied by %.1f: illegal Doppler scaling factors\n",
//					baddopscale_factor_mgpu);
//	}

	/* Reset showvals_mgpu to 0 if it had been 1 (i.e., turn off display of reduced
	 * chi-square and the individual penalty values).  */
	if (showvals_mgpu)
		fflush( stdout);
	showvals_mgpu = 0;

	free(hflags);
	free(hlogfactors);
	cudaFree(dflags);
	cudaFree(dlogfactors);
	return err;
}

