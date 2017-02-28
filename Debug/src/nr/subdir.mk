################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CU_SRCS += \
../src/nr/gammln_cuda.cu \
../src/nr/plgndr_cuda.cu 

C_SRCS += \
../src/nr/bcucof.c \
../src/nr/bcuint.c \
../src/nr/brent_abs.c \
../src/nr/bsstep.c \
../src/nr/caldat.c \
../src/nr/cel.c \
../src/nr/complex.c \
../src/nr/ellf.c \
../src/nr/gammln.c \
../src/nr/hpsort.c \
../src/nr/indexx.c \
../src/nr/jacobi.c \
../src/nr/julday.c \
../src/nr/laguer.c \
../src/nr/lubksb.c \
../src/nr/ludcmp.c \
../src/nr/mmid.c \
../src/nr/mnbrak.c \
../src/nr/nrutil.c \
../src/nr/odeint.c \
../src/nr/piksrt.c \
../src/nr/plgndr.c \
../src/nr/qsimp.c \
../src/nr/rf.c \
../src/nr/rtbis.c \
../src/nr/rzextr.c \
../src/nr/sncndn.c \
../src/nr/sort3.c \
../src/nr/spline.c \
../src/nr/splint.c \
../src/nr/svbksb.c \
../src/nr/svdcmp.c \
../src/nr/trapzd.c \
../src/nr/zroots.c 

OBJS += \
./src/nr/bcucof.o \
./src/nr/bcuint.o \
./src/nr/brent_abs.o \
./src/nr/bsstep.o \
./src/nr/caldat.o \
./src/nr/cel.o \
./src/nr/complex.o \
./src/nr/ellf.o \
./src/nr/gammln.o \
./src/nr/gammln_cuda.o \
./src/nr/hpsort.o \
./src/nr/indexx.o \
./src/nr/jacobi.o \
./src/nr/julday.o \
./src/nr/laguer.o \
./src/nr/lubksb.o \
./src/nr/ludcmp.o \
./src/nr/mmid.o \
./src/nr/mnbrak.o \
./src/nr/nrutil.o \
./src/nr/odeint.o \
./src/nr/piksrt.o \
./src/nr/plgndr.o \
./src/nr/plgndr_cuda.o \
./src/nr/qsimp.o \
./src/nr/rf.o \
./src/nr/rtbis.o \
./src/nr/rzextr.o \
./src/nr/sncndn.o \
./src/nr/sort3.o \
./src/nr/spline.o \
./src/nr/splint.o \
./src/nr/svbksb.o \
./src/nr/svdcmp.o \
./src/nr/trapzd.o \
./src/nr/zroots.o 

CU_DEPS += \
./src/nr/gammln_cuda.d \
./src/nr/plgndr_cuda.d 

C_DEPS += \
./src/nr/bcucof.d \
./src/nr/bcuint.d \
./src/nr/brent_abs.d \
./src/nr/bsstep.d \
./src/nr/caldat.d \
./src/nr/cel.d \
./src/nr/complex.d \
./src/nr/ellf.d \
./src/nr/gammln.d \
./src/nr/hpsort.d \
./src/nr/indexx.d \
./src/nr/jacobi.d \
./src/nr/julday.d \
./src/nr/laguer.d \
./src/nr/lubksb.d \
./src/nr/ludcmp.d \
./src/nr/mmid.d \
./src/nr/mnbrak.d \
./src/nr/nrutil.d \
./src/nr/odeint.d \
./src/nr/piksrt.d \
./src/nr/plgndr.d \
./src/nr/qsimp.d \
./src/nr/rf.d \
./src/nr/rtbis.d \
./src/nr/rzextr.d \
./src/nr/sncndn.d \
./src/nr/sort3.d \
./src/nr/spline.d \
./src/nr/splint.d \
./src/nr/svbksb.d \
./src/nr/svdcmp.d \
./src/nr/trapzd.d \
./src/nr/zroots.d 


# Each subdirectory must supply rules for building sources it contributes
src/nr/%.o: ../src/nr/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-8.0/bin/nvcc -I/home/matt/git/cfitsio -G -g -lineinfo -pg -O0 -gencode arch=compute_35,code=sm_35 -m64 -odir "src/nr" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-8.0/bin/nvcc -I/home/matt/git/cfitsio -G -g -lineinfo -pg -O0 --compile -m64  -x c -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/nr/%.o: ../src/nr/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-8.0/bin/nvcc -I/home/matt/git/cfitsio -G -g -lineinfo -pg -O0 -gencode arch=compute_35,code=sm_35 -m64 -odir "src/nr" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-8.0/bin/nvcc -I/home/matt/git/cfitsio -G -g -lineinfo -pg -O0 --compile --relocatable-device-code=true -gencode arch=compute_35,code=compute_35 -gencode arch=compute_35,code=sm_35 -m64  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


