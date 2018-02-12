################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/nr/bcucof.c \
../src/nr/bcuint.c \
../src/nr/bessj.c \
../src/nr/bessj0.c \
../src/nr/bessj1.c \
../src/nr/bessk.c \
../src/nr/bessk0.c \
../src/nr/bessk1.c \
../src/nr/bessy.c \
../src/nr/bessy0.c \
../src/nr/bessy1.c \
../src/nr/brent.c \
../src/nr/brent_abs.c \
../src/nr/bsstep.c \
../src/nr/caldat.c \
../src/nr/cel.c \
../src/nr/complex.c \
../src/nr/ellf.c \
../src/nr/f1dim.c \
../src/nr/factrl.c \
../src/nr/gammln.c \
../src/nr/gasdev.c \
../src/nr/hpsort.c \
../src/nr/indexx.c \
../src/nr/jacobi.c \
../src/nr/julday.c \
../src/nr/laguer.c \
../src/nr/linmin.c \
../src/nr/lubksb.c \
../src/nr/ludcmp.c \
../src/nr/mmid.c \
../src/nr/mnbrak.c \
../src/nr/nrutil.c \
../src/nr/odeint.c \
../src/nr/piksrt.c \
../src/nr/plgndr.c \
../src/nr/poidev.c \
../src/nr/powell.c \
../src/nr/qsimp.c \
../src/nr/ran1.c \
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
../src/nr/zbrent.c \
../src/nr/zroots.c 

OBJS += \
./src/nr/bcucof.o \
./src/nr/bcuint.o \
./src/nr/bessj.o \
./src/nr/bessj0.o \
./src/nr/bessj1.o \
./src/nr/bessk.o \
./src/nr/bessk0.o \
./src/nr/bessk1.o \
./src/nr/bessy.o \
./src/nr/bessy0.o \
./src/nr/bessy1.o \
./src/nr/brent.o \
./src/nr/brent_abs.o \
./src/nr/bsstep.o \
./src/nr/caldat.o \
./src/nr/cel.o \
./src/nr/complex.o \
./src/nr/ellf.o \
./src/nr/f1dim.o \
./src/nr/factrl.o \
./src/nr/gammln.o \
./src/nr/gasdev.o \
./src/nr/hpsort.o \
./src/nr/indexx.o \
./src/nr/jacobi.o \
./src/nr/julday.o \
./src/nr/laguer.o \
./src/nr/linmin.o \
./src/nr/lubksb.o \
./src/nr/ludcmp.o \
./src/nr/mmid.o \
./src/nr/mnbrak.o \
./src/nr/nrutil.o \
./src/nr/odeint.o \
./src/nr/piksrt.o \
./src/nr/plgndr.o \
./src/nr/poidev.o \
./src/nr/powell.o \
./src/nr/qsimp.o \
./src/nr/ran1.o \
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
./src/nr/zbrent.o \
./src/nr/zroots.o 

C_DEPS += \
./src/nr/bcucof.d \
./src/nr/bcuint.d \
./src/nr/bessj.d \
./src/nr/bessj0.d \
./src/nr/bessj1.d \
./src/nr/bessk.d \
./src/nr/bessk0.d \
./src/nr/bessk1.d \
./src/nr/bessy.d \
./src/nr/bessy0.d \
./src/nr/bessy1.d \
./src/nr/brent.d \
./src/nr/brent_abs.d \
./src/nr/bsstep.d \
./src/nr/caldat.d \
./src/nr/cel.d \
./src/nr/complex.d \
./src/nr/ellf.d \
./src/nr/f1dim.d \
./src/nr/factrl.d \
./src/nr/gammln.d \
./src/nr/gasdev.d \
./src/nr/hpsort.d \
./src/nr/indexx.d \
./src/nr/jacobi.d \
./src/nr/julday.d \
./src/nr/laguer.d \
./src/nr/linmin.d \
./src/nr/lubksb.d \
./src/nr/ludcmp.d \
./src/nr/mmid.d \
./src/nr/mnbrak.d \
./src/nr/nrutil.d \
./src/nr/odeint.d \
./src/nr/piksrt.d \
./src/nr/plgndr.d \
./src/nr/poidev.d \
./src/nr/powell.d \
./src/nr/qsimp.d \
./src/nr/ran1.d \
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
./src/nr/zbrent.d \
./src/nr/zroots.d 


# Each subdirectory must supply rules for building sources it contributes
src/nr/%.o: ../src/nr/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-9.1/bin/nvcc -G -g -O0 -gencode arch=compute_61,code=sm_61 -m64 -odir "src/nr" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-9.1/bin/nvcc -G -g -O0 --compile -m64  -x c -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


