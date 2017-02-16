extern "C" {
#include "../shape/head.h"
}
//__host__ void copy_par_to_device(struct par_t *hpar)
//{
//	/* NOTE:  The double pointer dev_par_fpntr contains pointers that point to
//	 * 		  host memory.  This won't work.  Fix later.  Though kernel
//	 * 		  debugging shows dev_par_fpartype = 105.6...
//	 */
//	int size_int = sizeof(int)*hpar->nfpar;
//	int size_dbl = sizeof(double)*hpar->nfpar;
//	int size_dblpntr = sizeof(double*)*hpar->nfpar;
//	int size_par = sizeof(hpar);
//
//	gpuErrchk(cudaMalloc((void**)&dev_par, 				size_par));
//	gpuErrchk(cudaMalloc((void**)&dev_par_fparstep,		size_dbl));
//	gpuErrchk(cudaMalloc((void**)&dev_par_fpartol,		size_dbl));
//	gpuErrchk(cudaMalloc((void**)&dev_par_fparabstol,	size_dbl));
//	gpuErrchk(cudaMalloc((void**)&dev_par_fpntr,		size_dblpntr));
//	gpuErrchk(cudaMalloc((void**)&dev_par_fpartype,		size_int));
//
//	gpuErrchk(cudaMemcpy(dev_par, hpar, size_par, cudaMemcpyHostToDevice));
//	gpuErrchk(cudaMemcpy(dev_par_fparstep, hpar->fparstep, size_dbl,
//			cudaMemcpyHostToDevice));
//	gpuErrchk(cudaMemcpy(dev_par_fpartol, hpar->fpartol, size_dbl,
//			cudaMemcpyHostToDevice));
//	gpuErrchk(cudaMemcpy(dev_par_fparabstol, hpar->fparabstol, size_dbl,
//			cudaMemcpyHostToDevice));
//	gpuErrchk(cudaMemcpy(dev_par_fpntr, hpar->fpntr, size_dblpntr,
//			cudaMemcpyHostToDevice));
//	gpuErrchk(cudaMemcpy(dev_par_fpartype, hpar->fpartype, size_int,
//			cudaMemcpyHostToDevice));
//
//}
//
//void copy_CUDA_structs(struct par_t *hpar, struct mod_t *hmod, struct dat_t *hdat)
//{
//	gpuErrchk(cudaMalloc((void**)&dev_par, sizeof(struct par_t)*1));
//	gpuErrchk(cudaMalloc((void**)&dev_mod, sizeof(hmod)*1));
//	gpuErrchk(cudaMalloc((void**)&dev_dat, sizeof(hdat)*1));
//
//	gpuErrchk(cudaMemcpy(dev_par, &hpar, sizeof(struct par_t), cudaMemcpyHostToDevice));
//	gpuErrchk(cudaMemcpy(dev_mod, &hmod, sizeof(hmod), cudaMemcpyHostToDevice));
//	gpuErrchk(cudaMemcpy(dev_dat, &hdat, sizeof(hdat), cudaMemcpyHostToDevice));
//
//}

__host__ void gpuAssert(cudaError_t code, const char *file, int line)
{
	if (code != cudaSuccess)
	{
		fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		exit(code);
	}
}

/* To-Do:  finish allocating double and triple pointers. */
/* To-Do:  anything with structures inside structures at the endpoints (such as lots of param_t's at the end of long chain structures)
 * 		   may need to have those final param_t's declared, allocated, copied
 * NOTE:  Most of the commented out allocations are those of vectors declared with definite size, i.e. double x[3]	*/
//__host__ void copy_mod_to_device(struct mod_t *hmod) {
//	/* Assumes single component model */
//	/*.................................................................................................................*/
//
//	/* Allocate and copy the main parent structure first */
//	gpuErrchk(cudaMalloc((void**)&dev_mod, sizeof(hmod)*1));
//	gpuErrchk(cudaMemcpy(dev_mod, &hmod, sizeof(hmod), cudaMemcpyHostToDevice));
//
//	/* Allocate mod->spin memory */
//	int angsz1 = sizeof(struct param_t) * 3;
//	int dblsz = sizeof(double) * MAXIMP;
//	int angsz2 = sizeof(struct param_t) * MAXIMP * 3;
//
//
//	gpuErrchk(cudaMalloc((void**)&dev_mod_spin, sizeof(hmod->spin)));
//	//gpuErrchk(cudaMalloc((void**)&dev_mod_spin_angle, sizeof(hmod->spin.angle)));
//	//	gpuErrchk(cudaMalloc((void**)&dev_mod_spin_omega, sizeof(hmod->spin.omega)));
//	//	gpuErrchk(cudaMalloc((void**)&dev_mod_spin_omegadot, sizeof(hmod->spin.omegadot)));
//	//	gpuErrchk(cudaMalloc((void**)&dev_mod_spin_t_impulse, sizeof(hmod->spin.t_impulse)));
//	//	gpuErrchk(cudaMalloc((void**)&dev_mod_spin_impulse, sizeof(hmod->spin.impulse)));
//	//	gpuErrchk(cudaMalloc((void**)&dev_mod_spin_inertia, sizeof(hmod->spin.inertia)));
//
//	/* Copy mod->spin contents */
//	gpuErrchk(cudaMemcpy(dev_mod_spin, &hmod->spin, sizeof(hmod->spin), cudaMemcpyHostToDevice));
//	//gpuErrchk(cudaMemcpy(dev_mod_spin_angle, &hmod->spin.angle, angsz1, cudaMemcpyHostToDevice));
//	//gpuErrchk(cudaMemcpy(dev_mod_spin_omega, hmod->spin.omega, angsz1, cudaMemcpyHostToDevice));
//	//gpuErrchk(cudaMemcpy(dev_mod_spin_omegadot, &hmod->spin.omegadot, angsz1, cudaMemcpyHostToDevice));
//	//gpuErrchk(cudaMemcpy(dev_mod_spin_t_impulse, hmod->spin.t_impulse, dblsz, cudaMemcpyHostToDevice));
//	//gpuErrchk(cudaMemcpy(dev_mod_spin_impulse, &hmod->spin.impulse, angsz2, cudaMemcpyHostToDevice));
//	//gpuErrchk(cudaMemcpy(dev_mod_spin_inertia, &hmod->spin.inertia, angsz1, cudaMemcpyHostToDevice));
//
//	/*..................................................................................................................*/
//
//	/* Allocate mod->shape */
//	/* mod->shape.comp[0].real (vertices_t structure) */
//	int cmp_sz = sizeof(hmod->shape.comp[0]);
//	int shp_sz = sizeof(hmod->shape);
//	int inertia_sz = sizeof(double) * 9;
//	int off_sz = sizeof(struct param_t) * 3;
//	int ver_sz = sizeof(struct vertices_t);
//	int f_sz = sizeof(struct facet_t);
//	int s_sz = sizeof(struct side_t);
//	int v_sz = sizeof(struct vertex_t);
//	int ns = hmod->shape.comp[0].real.ns;
//	int nf = hmod->shape.comp[0].real.nf;
//	int nv = hmod->shape.comp[0].real.nv;
//	int afactor = sizeof(double*) * nv;
//	int int_sz = sizeof(int);
//	int pint_sz = sizeof(int*);
//	int dbl_sz = sizeof(double);
//	int par_sz = sizeof(struct param_t);
//
//	gpuErrchk(cudaMalloc((void**)&dev_mod_shape, shp_sz));
//	gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp, cmp_sz));
//	//	gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_inertia, inertia_sz));
//	//	gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_com, inertia_sz/3));
//	//	gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_m, inertia_sz));
//	//	gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_off, off_sz));
//	//	gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_rot, off_sz));
//	gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_real, ver_sz));
//	gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_real_f, f_sz*nf));
//	gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_real_v, v_sz*nv));
//	gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_real_s, s_sz*ns));
//	gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_real_s_v[2], ns*int_sz));	// dbl *[]
//	gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_real_s_f[2], ns*int_sz));	// dbl *[]
//	gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_real_v_a[3], ns*dbl_sz));	// dbl *[]
//	gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_real_v_u[3], ns*dbl_sz));	// dbl *[]
//	gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_real_v_x[3], ns*dbl_sz));	// dbl *[]
//	gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_real_v_n[3], ns*dbl_sz));	// dbl *[]
//	gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_real_v_af, nv*pint_sz));	// dbl **
//	gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_real_v_as, nv*pint_sz));	// dbl **
//	gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_real_v_afactor, afactor));	// dbl ***
//	gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_real_v_bfactor, afactor));	// dbl ***
//	gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_real_f_v[3], nf*int_sz));	// dbl *[]
//	gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_real_f_s[3], nf*int_sz));	// dbl *[]
//	gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_real_f_n[3], nf*dbl_sz));	// dbl *[]
//	gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_real_f_x[3], nf*dbl_sz));	// dbl *[]
//
//	/* Copy mod->shape contents */
//	gpuErrchk(cudaMemcpy(dev_mod_shape, &hmod->shape, shp_sz, cudaMemcpyHostToDevice));
//	gpuErrchk(cudaMemcpy(dev_mod_shape_comp, &hmod->shape.comp[0], cmp_sz, cudaMemcpyHostToDevice));
//	gpuErrchk(cudaMemcpy(dev_mod_shape_comp_inertia, &hmod->shape.comp[0].inertia, inertia_sz, cudaMemcpyHostToDevice));
//	gpuErrchk(cudaMemcpy(dev_mod_shape_comp_com, &hmod->shape.comp[0].com, inertia_sz/3, cudaMemcpyHostToDevice));
//	gpuErrchk(cudaMemcpy(dev_mod_shape_comp_m, &hmod->shape.comp[0].m, inertia_sz, cudaMemcpyHostToDevice));
//	gpuErrchk(cudaMemcpy(dev_mod_shape_comp_off, &hmod->shape.comp[0].off, off_sz, cudaMemcpyHostToDevice));
//	gpuErrchk(cudaMemcpy(dev_mod_shape_comp_rot, &hmod->shape.comp[0].rot, off_sz, cudaMemcpyHostToDevice));
//	gpuErrchk(cudaMemcpy(dev_mod_shape_comp_real,&hmod->shape.comp[0].real,ver_sz, cudaMemcpyHostToDevice));
//	gpuErrchk(cudaMemcpy(dev_mod_shape_comp_real_f, &hmod->shape.comp[0].real.f, f_sz*nf, cudaMemcpyHostToDevice));
//	gpuErrchk(cudaMemcpy(dev_mod_shape_comp_real_v, &hmod->shape.comp[0].real.v, v_sz*nv, cudaMemcpyHostToDevice));
//	gpuErrchk(cudaMemcpy(dev_mod_shape_comp_real_s, &hmod->shape.comp[0].real.s, s_sz*ns, cudaMemcpyHostToDevice));
//	gpuErrchk(cudaMemcpy(dev_mod_shape_comp_real_s, &hmod->shape.comp[0].real.s->v, 2*ns*int_sz, cudaMemcpyHostToDevice));
//	gpuErrchk(cudaMemcpy(dev_mod_shape_comp_real_f, &hmod->shape.comp[0].real.s->f, 2*ns*int_sz, cudaMemcpyHostToDevice));
//	gpuErrchk(cudaMemcpy(dev_mod_shape_comp_real_v_a, &hmod->shape.comp[0].real.v->a, 3*nv*dbl_sz, cudaMemcpyHostToDevice));
//	gpuErrchk(cudaMemcpy(dev_mod_shape_comp_real_v_u, &hmod->shape.comp[0].real.v->u, 3*nv*dbl_sz, cudaMemcpyHostToDevice));
//	gpuErrchk(cudaMemcpy(dev_mod_shape_comp_real_v_x, &hmod->shape.comp[0].real.v->x, 3*nv*dbl_sz, cudaMemcpyHostToDevice));
//	gpuErrchk(cudaMemcpy(dev_mod_shape_comp_real_v_n, &hmod->shape.comp[0].real.v->n, 3*nv*dbl_sz, cudaMemcpyHostToDevice));
//	gpuErrchk(cudaMemcpy(dev_mod_shape_comp_real_v_af, &hmod->shape.comp[0].real.v->af, nv*pint_sz, cudaMemcpyHostToDevice));
//	gpuErrchk(cudaMemcpy(dev_mod_shape_comp_real_v_as, &hmod->shape.comp[0].real.v->as, nv*pint_sz, cudaMemcpyHostToDevice));
//	gpuErrchk(cudaMemcpy(dev_mod_shape_comp_real_v_afactor, &hmod->shape.comp[0].real.v->afactor, afactor, cudaMemcpyHostToDevice));
//	gpuErrchk(cudaMemcpy(dev_mod_shape_comp_real_v_bfactor, &hmod->shape.comp[0].real.v->bfactor, afactor, cudaMemcpyHostToDevice));
//	gpuErrchk(cudaMemcpy(dev_mod_shape_comp_real_f_v, &hmod->shape.comp[0].real.f->v, 3*nf*int_sz, cudaMemcpyHostToDevice));
//	gpuErrchk(cudaMemcpy(dev_mod_shape_comp_real_f_s, &hmod->shape.comp[0].real.f->s, 3*nf*int_sz, cudaMemcpyHostToDevice));
//	gpuErrchk(cudaMemcpy(dev_mod_shape_comp_real_f_n, &hmod->shape.comp[0].real.f->n, 3*nf*dbl_sz, cudaMemcpyHostToDevice));
//	gpuErrchk(cudaMemcpy(dev_mod_shape_comp_real_f_x, &hmod->shape.comp[0].real.f->x, 3*nf*dbl_sz, cudaMemcpyHostToDevice));
//
//
//	if (hmod->shape.comp[0].type == ELLIPSE) {
//		int ell_sz = sizeof(struct ellipse_t);
//		gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_desc_ell, ell_sz));
//		gpuErrchk(cudaMemcpy(dev_mod_shape_comp_desc_ell, &hmod->shape.comp[0].desc.ell, ell_sz, cudaMemcpyHostToDevice));
//	}
//
//	if (hmod->shape.comp[0].type == OVOID) {
//		int ov_sz = sizeof(struct ovoid_t);
//		gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_desc_ovoid, ov_sz));
//		gpuErrchk(cudaMemcpy(dev_mod_shape_comp_desc_ovoid, &hmod->shape.comp[0].desc.ovoid,  ov_sz, cudaMemcpyHostToDevice));
//	}
//
//	/* To-Do:  These double pointer allocations and memcpy's need attention */
//	if (hmod->shape.comp[0].type == HARMONIC) {
//		int har_sz = sizeof(struct harmonic_t);
//		gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_desc_har, har_sz));
//		gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_desc_har_a, par_sz));			// dbl **
//		gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_desc_har_b, par_sz));			// dbl **
//		gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_desc_har_a_save, dbl_sz));		// dbl **
//		gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_desc_har_b_save, dbl_sz));		// dbl **
//
//		gpuErrchk(cudaMemcpy(dev_mod_shape_comp_desc_har, &hmod->shape.comp[0].desc.har, har_sz, cudaMemcpyHostToDevice));
//		gpuErrchk(cudaMemcpy(dev_mod_shape_comp_desc_har_a, &hmod->shape.comp[0].desc.har.a, par_sz, cudaMemcpyHostToDevice));
//		gpuErrchk(cudaMemcpy(dev_mod_shape_comp_desc_har_b, &hmod->shape.comp[0].desc.har.b, par_sz, cudaMemcpyHostToDevice));
//		gpuErrchk(cudaMemcpy(dev_mod_shape_comp_desc_har_a_save, &hmod->shape.comp[0].desc.har.a_save, dbl_sz, cudaMemcpyHostToDevice));
//		gpuErrchk(cudaMemcpy(dev_mod_shape_comp_desc_har_b_save, &hmod->shape.comp[0].desc.har.b_save, dbl_sz, cudaMemcpyHostToDevice));
//	}
//
//	if (hmod->shape.comp[0].type == VERTEX) {
//		int nv1 = hmod->shape.comp[0].desc.ver.nv;
//		int ns1 = hmod->shape.comp[0].desc.ver.ns;
//		int nf1 = hmod->shape.comp[0].desc.ver.nf;
//		int dpt_sz = sizeof(double*);
//		gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_desc_ver, ver_sz));
//		gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_desc_ver_s, s_sz*ns1));
//		//		gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_desc_ver_s_v[2], int_sz*ns1));
//		//		gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_desc_ver_s_f[2], int_sz*ns1));
//		gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_desc_ver_f, f_sz*nf1));
//		//		gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_desc_ver_f_v[3], int_sz*nf1));
//		//		gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_desc_ver_f_s[3], int_sz*nf1));
//		//		gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_desc_ver_f_n[3], dbl_sz*nf1));
//		//		gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_desc_ver_f_x[3], dbl_sz*nf1));
//		gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_desc_ver_v, v_sz*nv1));
//		//		gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_desc_ver_v_a[3], dbl_sz*nv1));
//		//		gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_desc_ver_v_u[3], dbl_sz*nv1));
//		//		gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_desc_ver_v_x[3], dbl_sz*nv1));
//		//		gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_desc_ver_v_n[3], dbl_sz*nv1));
//		gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_desc_ver_v_af, nv1*int_sz));	// dbl **
//		gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_desc_ver_v_as, nv1*int_sz));	// dbl **
//		gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_desc_ver_v_afactor, nv1*dbl_sz));	// dbl ***
//		gpuErrchk(cudaMalloc((void**)&dev_mod_shape_comp_desc_ver_v_bfactor, nv1*dbl_sz));	// dbl ***
//
//		gpuErrchk(cudaMemcpy(dev_mod_shape_comp_desc_ver, &hmod->shape.comp[0].desc.ver, ver_sz, cudaMemcpyHostToDevice));
//		gpuErrchk(cudaMemcpy(dev_mod_shape_comp_desc_ver_s, &hmod->shape.comp[0].desc.ver.s, s_sz*ns1, cudaMemcpyHostToDevice));
//		gpuErrchk(cudaMemcpy(dev_mod_shape_comp_desc_ver_s_v, &hmod->shape.comp[0].desc.ver.s->v, 2*int_sz*ns1, cudaMemcpyHostToDevice));
//		gpuErrchk(cudaMemcpy(dev_mod_shape_comp_desc_ver_s_f, &hmod->shape.comp[0].desc.ver.s->f, 2*int_sz*ns1, cudaMemcpyHostToDevice));
//		gpuErrchk(cudaMemcpy(dev_mod_shape_comp_desc_ver_f, &hmod->shape.comp[0].desc.ver.f, f_sz*nf1, cudaMemcpyHostToDevice));
//		gpuErrchk(cudaMemcpy(dev_mod_shape_comp_desc_ver_f_v, &hmod->shape.comp[0].desc.ver.f->v, 3*int_sz*nf1, cudaMemcpyHostToDevice));
//		gpuErrchk(cudaMemcpy(dev_mod_shape_comp_desc_ver_f_s, &hmod->shape.comp[0].desc.ver.f->s, 3*int_sz*nf1, cudaMemcpyHostToDevice));
//		gpuErrchk(cudaMemcpy(dev_mod_shape_comp_desc_ver_f_n, &hmod->shape.comp[0].desc.ver.f->n, 3*dbl_sz*nf1, cudaMemcpyHostToDevice));
//		gpuErrchk(cudaMemcpy(dev_mod_shape_comp_desc_ver_f_x, &hmod->shape.comp[0].desc.ver.f->x, 3*dbl_sz*nf1, cudaMemcpyHostToDevice));
//		gpuErrchk(cudaMemcpy(dev_mod_shape_comp_desc_ver_v, &hmod->shape.comp[0].desc.ver.v, v_sz*nv1, cudaMemcpyHostToDevice));
//		gpuErrchk(cudaMemcpy(dev_mod_shape_comp_desc_ver_v_a, &hmod->shape.comp[0].desc.ver.v->a, 3*int_sz*nf1, cudaMemcpyHostToDevice));
//		gpuErrchk(cudaMemcpy(dev_mod_shape_comp_desc_ver_v_u, &hmod->shape.comp[0].desc.ver.v->u, 3*int_sz*nf1, cudaMemcpyHostToDevice));
//		gpuErrchk(cudaMemcpy(dev_mod_shape_comp_desc_ver_v_x, &hmod->shape.comp[0].desc.ver.v->x, 3*dbl_sz*nf1, cudaMemcpyHostToDevice));
//		gpuErrchk(cudaMemcpy(dev_mod_shape_comp_desc_ver_v_n, &hmod->shape.comp[0].desc.ver.v->n, 3*dbl_sz*nf1, cudaMemcpyHostToDevice));
//		gpuErrchk(cudaMemcpy(dev_mod_shape_comp_desc_ver_v_af, &hmod->shape.comp[0].desc.ver.v->af, pint_sz, cudaMemcpyHostToDevice));
//		gpuErrchk(cudaMemcpy(dev_mod_shape_comp_desc_ver_v_as, &hmod->shape.comp[0].desc.ver.v->as, pint_sz, cudaMemcpyHostToDevice));
//		gpuErrchk(cudaMemcpy(dev_mod_shape_comp_desc_ver_v_afactor, &hmod->shape.comp[0].desc.ver.v->afactor, dpt_sz, cudaMemcpyHostToDevice));
//		gpuErrchk(cudaMemcpy(dev_mod_shape_comp_desc_ver_v_bfactor, &hmod->shape.comp[0].desc.ver.v->bfactor, dpt_sz, cudaMemcpyHostToDevice));
//	}
//
//	/*............................................................................................................*/
//	/* Allocate mod->photo device pointers */
//	/* NOTE:  the following allocation and copying of mod->photo contents assumes nrl = nol = 0 or 1, but no more */
//	/* First some helper variables for sizing */
//	int nrl = hmod->photo.nradlaws;
//	int nol = hmod->photo.noptlaws;
//	int u_sz = sizeof(unsigned char);
//	int pho_sz = sizeof(struct photo_t);
//	int tab_sz = sizeof(struct tabular_t);
//	int rc_sz = sizeof(struct RC_t);
//	int rcpt_sz = sizeof(struct RC_t*);
//	int parptr_sz = sizeof(struct param_t*) * hmod->photo.radar->tabular.n;
//	int dblptr_sz = sizeof(double*) * hmod->photo.radar->tabular.n;
//	int qspc_sz = sizeof(struct quasispec_t);
//	int hyr_sz = sizeof(struct hybridradar_t);
//	int hrcs_sz = sizeof(struct harmcosine_t);
//	int incs_sz = sizeof(struct inhocosine_t);
//	int r_sz = sizeof(struct R_t);
//	int rpt_sz = sizeof(struct R_t*);
//	int hmR_sz = sizeof(struct harmR_t);
//	int inho_sz = sizeof(struct inhoR_t);
//	int hpk_sz = sizeof(struct hapke_t);
//	int phpk_sz = sizeof(struct hapke_t*);
//	int hmhpk_sz = sizeof(struct harmhapke_t);
//	int inhpk_sz = sizeof(struct inhohapke_t);
//	int kas_sz = sizeof(struct kaas_t);
//	int pkas_sz = sizeof(struct kaas_t*);
//	int hmkas_sz = sizeof(struct harmkaas_t);
//	int inkas_sz = sizeof(struct inhokaas_t);
//
//	/* Allocate gpu copy of mod->photo and memcpy */
//	gpuErrchk(cudaMalloc((void**)&dev_mod_photo, sizeof(hmod->photo)));
//	gpuErrchk(cudaMalloc((void**)&dev_mod_photo_radtype, u_sz*nrl));
//	gpuErrchk(cudaMalloc((void**)&dev_mod_photo_opttype, u_sz*nol));
//	gpuErrchk(cudaMemcpy(dev_mod_photo, &hmod->photo, pho_sz, cudaMemcpyHostToDevice));
//	gpuErrchk(cudaMemcpy(dev_mod_photo_radtype, &hmod->photo.radtype, u_sz*nrl, cudaMemcpyHostToDevice));
//	gpuErrchk(cudaMemcpy(dev_mod_photo_opttype, &hmod->photo.opttype, u_sz*nol, cudaMemcpyHostToDevice));
//
//	/* Check for nrl and nol being 0 or 1 */
//	if (nrl >= 1) { /* This checks that nradlaws is at least 1 */
//
//		if (nrl > 1)
//			printf("\nShape-cuda V1.0 currently supports only one optical and one radar law maximum.");
//
//		/* Allocate gpu copy of mod->photo.radar and memcpy, depending on type */
//		if (hmod->photo.radtype[0] == COSINELAW_DIFF) {
//
//			gpuErrchk(cudaMalloc((void**)&dev_mod_photo_radar_RC, rc_sz));
//			gpuErrchk(cudaMemcpy(dev_mod_photo_radar_RC, &hmod->photo.radar->RC, rc_sz, cudaMemcpyHostToDevice));
//		}
//		if (hmod->photo.radtype[0] == TABULARLAW) {
//			gpuErrchk(cudaMalloc((void**)&dev_mod_photo_radar_tabular, tab_sz));
//			gpuErrchk(cudaMalloc((void**)&dev_mod_photo_radar_tabular_rho, parptr_sz));
//			gpuErrchk(cudaMalloc((void**)&dev_mod_photo_radar_tabular_rho_save, dblptr_sz));
//			gpuErrchk(cudaMemcpy(dev_mod_photo_radar_tabular, &hmod->photo.radar->tabular,  tab_sz, cudaMemcpyHostToDevice));
//			gpuErrchk(cudaMemcpy(dev_mod_photo_radar_tabular_rho, &hmod->photo.radar->tabular.rho, parptr_sz, cudaMemcpyHostToDevice));
//			gpuErrchk(cudaMemcpy(dev_mod_photo_radar_tabular_rho_save, &hmod->photo.radar->tabular.rho_save, dblptr_sz, cudaMemcpyHostToDevice));
//		}
//		if (hmod->photo.radtype[0] == GAUSSIANLAW || hmod->photo.radtype[0] == HAGFORSLAW || hmod->photo.radtype[0] == COSINELAW_QS) {
//			gpuErrchk(cudaMalloc((void**)&dev_mod_photo_radar_quasispec, qspc_sz));
//			gpuErrchk(cudaMemcpy(dev_mod_photo_radar_quasispec, &hmod->photo.radar->quasispec, qspc_sz, cudaMemcpyHostToDevice));
//		}
//		if (hmod->photo.radtype[0] == GAUSSIAN_COSINE || hmod->photo.radtype[0] == HAGFORS_COSINE || hmod->photo.radtype[0] == COSINE_COSINE) {
//			gpuErrchk(cudaMalloc((void**)&dev_mod_photo_radar_hybrid, hyr_sz));
//			gpuErrchk(cudaMemcpy(dev_mod_photo_radar_hybrid, &hmod->photo.radar->hybrid, hyr_sz, cudaMemcpyHostToDevice));
//		}
//		if (hmod->photo.radtype[0] == HARMCOSINE_DIFF) {
//			gpuErrchk(cudaMalloc((void**)&dev_mod_photo_radar_harmcosine, hrcs_sz));
//			gpuErrchk(cudaMalloc((void**)&dev_mod_photo_radar_harmcosine_local, rcpt_sz)); // so what is this - a [][] pointer of RC_t's.  Per pixel? How many?
//			gpuErrchk(cudaMemcpy(dev_mod_photo_radar_harmcosine, &hmod->photo.radar->harmcosine, hrcs_sz, cudaMemcpyHostToDevice));
//			gpuErrchk(cudaMemcpy(dev_mod_photo_radar_harmcosine_local, &hmod->photo.radar->harmcosine.local, rcpt_sz, cudaMemcpyHostToDevice)); //dbl pointer
//		}
//		if (hmod->photo.radtype[0] == INHOCOSINE_DIFF) {
//			gpuErrchk(cudaMalloc((void**)&dev_mod_photo_radar_inhocosine, incs_sz));
//			gpuErrchk(cudaMalloc((void**)&dev_mod_photo_radar_inhocosine_local, rcpt_sz)); // so what is this - a [][] pointer of RC_t's.  Per pixel? How many?
//			gpuErrchk(cudaMemcpy(dev_mod_photo_radar_inhocosine, &hmod->photo.radar->inhocosine, incs_sz, cudaMemcpyHostToDevice));
//			gpuErrchk(cudaMemcpy(dev_mod_photo_radar_inhocosine_local, &hmod->photo.radar->inhocosine.local, rcpt_sz, cudaMemcpyHostToDevice)); //dbl pointer
//		}
//	}
//
//	if (nol >= 1) {/* This checks that noptlaws is at least 1 */
//
//		if (nol > 1)
//			printf("\nShape-cuda V1.0 currently supports only one optical and one radar law maximum.");
//
//		/* Allocate gpu copy of mod->photo.optical and memcpy, depending on type */
//		if (hmod->photo.opttype[0] == GEOMETRICAL || hmod->photo.opttype[0] == LAMBERTLAW || hmod->photo.opttype[0] == LOMMEL) {
//			gpuErrchk(cudaMalloc((void**)&dev_mod_photo_optical_R, r_sz));
//			gpuErrchk(cudaMemcpy(dev_mod_photo_optical_R, &hmod->photo.optical->R, r_sz, cudaMemcpyHostToDevice));
//		}
//		if (hmod->photo.opttype[0] == HARMLAMBERT || hmod->photo.opttype[0] == HARMLOMMEL) {
//			gpuErrchk(cudaMalloc((void**)&dev_mod_photo_optical_harmR, hmR_sz));
//			gpuErrchk(cudaMalloc((void**)&dev_mod_photo_optical_harmR_local, rpt_sz));	// double pointer, but how many?
//			gpuErrchk(cudaMemcpy(dev_mod_photo_optical_harmR, &hmod->photo.optical->harmR, hmR_sz, cudaMemcpyHostToDevice));
//			gpuErrchk(cudaMemcpy(dev_mod_photo_optical_harmR_local, &hmod->photo.optical->harmR.local, rpt_sz, cudaMemcpyHostToDevice));
//		}
//		if (hmod->photo.opttype[0] == INHOLAMBERT || hmod->photo.opttype[0] == INHOLOMMEL) {
//			gpuErrchk(cudaMalloc((void**)&dev_mod_photo_optical_inhoR, inho_sz));
//			gpuErrchk(cudaMalloc((void**)&dev_mod_photo_optical_inhoR_local, rpt_sz));	// double pointer, but how many?
//			gpuErrchk(cudaMemcpy(dev_mod_photo_optical_inhoR, &hmod->photo.optical->inhoR, inho_sz, cudaMemcpyHostToDevice));
//			gpuErrchk(cudaMemcpy(dev_mod_photo_optical_inhoR_local, &hmod->photo.optical->inhoR.local, rpt_sz, cudaMemcpyHostToDevice));
//		}
//		if (hmod->photo.opttype[0] == HAPKE) {
//			gpuErrchk(cudaMalloc((void**)&dev_mod_photo_optical_hapke, hpk_sz));
//			gpuErrchk(cudaMemcpy(dev_mod_photo_optical_hapke, &hmod->photo.optical->hapke, hpk_sz, cudaMemcpyHostToDevice));
//		}
//		if (hmod->photo.opttype[0] == HARMHAPKE) {
//			gpuErrchk(cudaMalloc((void**)&dev_mod_photo_optical_harmhapke, hmhpk_sz));
//			gpuErrchk(cudaMalloc((void**)&dev_mod_photo_optical_harmhapke_local, phpk_sz));	// double pointer, but how many?
//			gpuErrchk(cudaMemcpy(dev_mod_photo_optical_harmhapke, &hmod->photo.optical->harmhapke, hmhpk_sz, cudaMemcpyHostToDevice));
//			gpuErrchk(cudaMemcpy(dev_mod_photo_optical_harmhapke_local, &hmod->photo.optical->harmhapke.local, phpk_sz, cudaMemcpyHostToDevice));
//		}
//		if (hmod->photo.opttype[0] == INHOHAPKE) {
//			gpuErrchk(cudaMalloc((void**)&dev_mod_photo_optical_inhohapke, inhpk_sz));
//			gpuErrchk(cudaMalloc((void**)&dev_mod_photo_optical_inhohapke_local, phpk_sz));	// double pointer, but how many?
//			gpuErrchk(cudaMemcpy(dev_mod_photo_optical_inhohapke, &hmod->photo.optical->inhohapke, inhpk_sz, cudaMemcpyHostToDevice));
//			gpuErrchk(cudaMemcpy(dev_mod_photo_optical_inhohapke_local, &hmod->photo.optical->inhohapke.local, phpk_sz, cudaMemcpyHostToDevice));
//		}
//		if (hmod->photo.opttype[0] == KAASALAINEN) {
//			gpuErrchk(cudaMalloc((void**)&dev_mod_photo_optical_kaas, kas_sz));
//			gpuErrchk(cudaMemcpy(dev_mod_photo_optical_kaas, &hmod->photo.optical->kaas, kas_sz, cudaMemcpyHostToDevice));
//		}
//		if (hmod->photo.opttype[0] == HARMKAAS) {
//			gpuErrchk(cudaMalloc((void**)&dev_mod_photo_optical_harmkaas, hmkas_sz));
//			gpuErrchk(cudaMalloc((void**)&dev_mod_photo_optical_harmkaas_local, pkas_sz));	// double pointer, but how many?
//			gpuErrchk(cudaMemcpy(dev_mod_photo_optical_harmkaas, &hmod->photo.optical->harmkaas, hmkas_sz, cudaMemcpyHostToDevice));
//			gpuErrchk(cudaMemcpy(dev_mod_photo_optical_harmkaas_local, &hmod->photo.optical->harmkaas.local, pkas_sz, cudaMemcpyHostToDevice));
//		}
//		if (hmod->photo.opttype[0] == INHOKAAS) {
//			gpuErrchk(cudaMalloc((void**)&dev_mod_photo_optical_inhokaas, inkas_sz));
//			gpuErrchk(cudaMalloc((void**)&dev_mod_photo_optical_inhokaas_local, pkas_sz));	// double pointer, but how many?
//			gpuErrchk(cudaMemcpy(dev_mod_photo_optical_inhokaas, &hmod->photo.optical->inhokaas, inkas_sz, cudaMemcpyHostToDevice));
//			gpuErrchk(cudaMemcpy(dev_mod_photo_optical_inhokaas_local, &hmod->photo.optical->inhokaas.local, pkas_sz, cudaMemcpyHostToDevice));
//		}
//	}
//}

