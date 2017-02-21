################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CU_SRCS += \
../src/shape/apply_photo_cuda.cu \
../src/shape/bestfit_CUDA.cu \
../src/shape/bsstep_cuda.cu \
../src/shape/calc_fits_cuda.cu \
../src/shape/calc_fits_cuda_af.cu \
../src/shape/cel_cuda.cu \
../src/shape/chi2_cuda.cu \
../src/shape/chi2_cuda_af.cu \
../src/shape/copy_to_device.cu \
../src/shape/debug_routines.cu \
../src/shape/deldopoffs_cuda.cu \
../src/shape/dopoffs_cuda.cu \
../src/shape/facmom_cuda.cu \
../src/shape/gamma_trans_cuda.cu \
../src/shape/inteuler_cuda.cu \
../src/shape/mkparlist_cuda.cu \
../src/shape/mmid_cuda.cu \
../src/shape/odeint_cuda.cu \
../src/shape/penalties_cuda.cu \
../src/shape/pos2deldop_cuda.cu \
../src/shape/pos2deldop_cuda_2.cu \
../src/shape/pos2deldop_cuda_af.cu \
../src/shape/pos2doppler-cuda.cu \
../src/shape/pos2doppler_cuda_2.cu \
../src/shape/pos2doppler_cuda_af.cu \
../src/shape/posvis-cuda.cu \
../src/shape/posvis-cuda_2.cu \
../src/shape/posvis_streams.cu \
../src/shape/realize_delcor_cuda.cu \
../src/shape/realize_dopscale_cuda.cu \
../src/shape/realize_mod_cuda.cu \
../src/shape/realize_photo_cuda.cu \
../src/shape/realize_spin_cuda.cu \
../src/shape/realize_spin_cuda_af.cu \
../src/shape/realize_xyoff_cuda.cu \
../src/shape/reduction.cu \
../src/shape/show_deldoplim_cuda.cu \
../src/shape/vary_params_cuda.cu \
../src/shape/vary_params_cuda_af.cu 

C_SRCS += \
../src/shape/apply_photo.c \
../src/shape/bailout.c \
../src/shape/bestfit.c \
../src/shape/calc_fits.c \
../src/shape/calc_orbit.c \
../src/shape/check_cuda.c \
../src/shape/chi2.c \
../src/shape/convex_hull.c \
../src/shape/covar.c \
../src/shape/delcorinit.c \
../src/shape/deldopoffs.c \
../src/shape/diag_inertia.c \
../src/shape/dopoffs.c \
../src/shape/ephem2mat.c \
../src/shape/facmom.c \
../src/shape/facnrm.c \
../src/shape/gamma_trans.c \
../src/shape/global_cuda.c \
../src/shape/init.c \
../src/shape/inteuler.c \
../src/shape/map_radar.c \
../src/shape/merge_comps.c \
../src/shape/mirror_dat.c \
../src/shape/mirror_mod.c \
../src/shape/mkparlist.c \
../src/shape/penalties.c \
../src/shape/photofacets.c \
../src/shape/photoharm.c \
../src/shape/pos2deldop.c \
../src/shape/pos2doppler.c \
../src/shape/posclr.c \
../src/shape/posmask.c \
../src/shape/posvis.c \
../src/shape/proj_area.c \
../src/shape/radlaw.c \
../src/shape/rayfacint.c \
../src/shape/read_dat.c \
../src/shape/read_ephem.c \
../src/shape/read_mod.c \
../src/shape/read_par.c \
../src/shape/readparam.c \
../src/shape/realize_delcor.c \
../src/shape/realize_dopscale.c \
../src/shape/realize_mod.c \
../src/shape/realize_photo.c \
../src/shape/realize_spin.c \
../src/shape/realize_xyoff.c \
../src/shape/ref_mod.c \
../src/shape/sample_mod.c \
../src/shape/shape-cuda-v1.0.c \
../src/shape/show_deldoplim.c \
../src/shape/show_moments.c \
../src/shape/slice.c \
../src/shape/split_mod.c \
../src/shape/vary_params.c \
../src/shape/view_mod.c \
../src/shape/write_dat.c \
../src/shape/write_mod.c \
../src/shape/write_pos.c \
../src/shape/write_wf.c 

OBJS += \
./src/shape/apply_photo.o \
./src/shape/apply_photo_cuda.o \
./src/shape/bailout.o \
./src/shape/bestfit.o \
./src/shape/bestfit_CUDA.o \
./src/shape/bsstep_cuda.o \
./src/shape/calc_fits.o \
./src/shape/calc_fits_cuda.o \
./src/shape/calc_fits_cuda_af.o \
./src/shape/calc_orbit.o \
./src/shape/cel_cuda.o \
./src/shape/check_cuda.o \
./src/shape/chi2.o \
./src/shape/chi2_cuda.o \
./src/shape/chi2_cuda_af.o \
./src/shape/convex_hull.o \
./src/shape/copy_to_device.o \
./src/shape/covar.o \
./src/shape/debug_routines.o \
./src/shape/delcorinit.o \
./src/shape/deldopoffs.o \
./src/shape/deldopoffs_cuda.o \
./src/shape/diag_inertia.o \
./src/shape/dopoffs.o \
./src/shape/dopoffs_cuda.o \
./src/shape/ephem2mat.o \
./src/shape/facmom.o \
./src/shape/facmom_cuda.o \
./src/shape/facnrm.o \
./src/shape/gamma_trans.o \
./src/shape/gamma_trans_cuda.o \
./src/shape/global_cuda.o \
./src/shape/init.o \
./src/shape/inteuler.o \
./src/shape/inteuler_cuda.o \
./src/shape/map_radar.o \
./src/shape/merge_comps.o \
./src/shape/mirror_dat.o \
./src/shape/mirror_mod.o \
./src/shape/mkparlist.o \
./src/shape/mkparlist_cuda.o \
./src/shape/mmid_cuda.o \
./src/shape/odeint_cuda.o \
./src/shape/penalties.o \
./src/shape/penalties_cuda.o \
./src/shape/photofacets.o \
./src/shape/photoharm.o \
./src/shape/pos2deldop.o \
./src/shape/pos2deldop_cuda.o \
./src/shape/pos2deldop_cuda_2.o \
./src/shape/pos2deldop_cuda_af.o \
./src/shape/pos2doppler-cuda.o \
./src/shape/pos2doppler.o \
./src/shape/pos2doppler_cuda_2.o \
./src/shape/pos2doppler_cuda_af.o \
./src/shape/posclr.o \
./src/shape/posmask.o \
./src/shape/posvis-cuda.o \
./src/shape/posvis-cuda_2.o \
./src/shape/posvis.o \
./src/shape/posvis_streams.o \
./src/shape/proj_area.o \
./src/shape/radlaw.o \
./src/shape/rayfacint.o \
./src/shape/read_dat.o \
./src/shape/read_ephem.o \
./src/shape/read_mod.o \
./src/shape/read_par.o \
./src/shape/readparam.o \
./src/shape/realize_delcor.o \
./src/shape/realize_delcor_cuda.o \
./src/shape/realize_dopscale.o \
./src/shape/realize_dopscale_cuda.o \
./src/shape/realize_mod.o \
./src/shape/realize_mod_cuda.o \
./src/shape/realize_photo.o \
./src/shape/realize_photo_cuda.o \
./src/shape/realize_spin.o \
./src/shape/realize_spin_cuda.o \
./src/shape/realize_spin_cuda_af.o \
./src/shape/realize_xyoff.o \
./src/shape/realize_xyoff_cuda.o \
./src/shape/reduction.o \
./src/shape/ref_mod.o \
./src/shape/sample_mod.o \
./src/shape/shape-cuda-v1.0.o \
./src/shape/show_deldoplim.o \
./src/shape/show_deldoplim_cuda.o \
./src/shape/show_moments.o \
./src/shape/slice.o \
./src/shape/split_mod.o \
./src/shape/vary_params.o \
./src/shape/vary_params_cuda.o \
./src/shape/vary_params_cuda_af.o \
./src/shape/view_mod.o \
./src/shape/write_dat.o \
./src/shape/write_mod.o \
./src/shape/write_pos.o \
./src/shape/write_wf.o 

CU_DEPS += \
./src/shape/apply_photo_cuda.d \
./src/shape/bestfit_CUDA.d \
./src/shape/bsstep_cuda.d \
./src/shape/calc_fits_cuda.d \
./src/shape/calc_fits_cuda_af.d \
./src/shape/cel_cuda.d \
./src/shape/chi2_cuda.d \
./src/shape/chi2_cuda_af.d \
./src/shape/copy_to_device.d \
./src/shape/debug_routines.d \
./src/shape/deldopoffs_cuda.d \
./src/shape/dopoffs_cuda.d \
./src/shape/facmom_cuda.d \
./src/shape/gamma_trans_cuda.d \
./src/shape/inteuler_cuda.d \
./src/shape/mkparlist_cuda.d \
./src/shape/mmid_cuda.d \
./src/shape/odeint_cuda.d \
./src/shape/penalties_cuda.d \
./src/shape/pos2deldop_cuda.d \
./src/shape/pos2deldop_cuda_2.d \
./src/shape/pos2deldop_cuda_af.d \
./src/shape/pos2doppler-cuda.d \
./src/shape/pos2doppler_cuda_2.d \
./src/shape/pos2doppler_cuda_af.d \
./src/shape/posvis-cuda.d \
./src/shape/posvis-cuda_2.d \
./src/shape/posvis_streams.d \
./src/shape/realize_delcor_cuda.d \
./src/shape/realize_dopscale_cuda.d \
./src/shape/realize_mod_cuda.d \
./src/shape/realize_photo_cuda.d \
./src/shape/realize_spin_cuda.d \
./src/shape/realize_spin_cuda_af.d \
./src/shape/realize_xyoff_cuda.d \
./src/shape/reduction.d \
./src/shape/show_deldoplim_cuda.d \
./src/shape/vary_params_cuda.d \
./src/shape/vary_params_cuda_af.d 

C_DEPS += \
./src/shape/apply_photo.d \
./src/shape/bailout.d \
./src/shape/bestfit.d \
./src/shape/calc_fits.d \
./src/shape/calc_orbit.d \
./src/shape/check_cuda.d \
./src/shape/chi2.d \
./src/shape/convex_hull.d \
./src/shape/covar.d \
./src/shape/delcorinit.d \
./src/shape/deldopoffs.d \
./src/shape/diag_inertia.d \
./src/shape/dopoffs.d \
./src/shape/ephem2mat.d \
./src/shape/facmom.d \
./src/shape/facnrm.d \
./src/shape/gamma_trans.d \
./src/shape/global_cuda.d \
./src/shape/init.d \
./src/shape/inteuler.d \
./src/shape/map_radar.d \
./src/shape/merge_comps.d \
./src/shape/mirror_dat.d \
./src/shape/mirror_mod.d \
./src/shape/mkparlist.d \
./src/shape/penalties.d \
./src/shape/photofacets.d \
./src/shape/photoharm.d \
./src/shape/pos2deldop.d \
./src/shape/pos2doppler.d \
./src/shape/posclr.d \
./src/shape/posmask.d \
./src/shape/posvis.d \
./src/shape/proj_area.d \
./src/shape/radlaw.d \
./src/shape/rayfacint.d \
./src/shape/read_dat.d \
./src/shape/read_ephem.d \
./src/shape/read_mod.d \
./src/shape/read_par.d \
./src/shape/readparam.d \
./src/shape/realize_delcor.d \
./src/shape/realize_dopscale.d \
./src/shape/realize_mod.d \
./src/shape/realize_photo.d \
./src/shape/realize_spin.d \
./src/shape/realize_xyoff.d \
./src/shape/ref_mod.d \
./src/shape/sample_mod.d \
./src/shape/shape-cuda-v1.0.d \
./src/shape/show_deldoplim.d \
./src/shape/show_moments.d \
./src/shape/slice.d \
./src/shape/split_mod.d \
./src/shape/vary_params.d \
./src/shape/view_mod.d \
./src/shape/write_dat.d \
./src/shape/write_mod.d \
./src/shape/write_pos.d \
./src/shape/write_wf.d 


# Each subdirectory must supply rules for building sources it contributes
src/shape/%.o: ../src/shape/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-8.0/bin/nvcc -G -g -lineinfo -pg -O0 -gencode arch=compute_35,code=sm_35 -m64 -odir "src/shape" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-8.0/bin/nvcc -G -g -lineinfo -pg -O0 --compile -m64  -x c -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/shape/%.o: ../src/shape/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-8.0/bin/nvcc -G -g -lineinfo -pg -O0 -gencode arch=compute_35,code=sm_35 -m64 -odir "src/shape" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-8.0/bin/nvcc -G -g -lineinfo -pg -O0 --compile --relocatable-device-code=true -gencode arch=compute_35,code=compute_35 -gencode arch=compute_35,code=sm_35 -m64  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


