################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CU_SRCS += \
../src/astro/\ hapke_cuda.cu 

C_SRCS += \
../src/astro/ec2eq.c \
../src/astro/eq2ec.c \
../src/astro/hapke.c \
../src/astro/jdcal.c \
../src/astro/kepler.c 

OBJS += \
./src/astro/\ hapke_cuda.o \
./src/astro/ec2eq.o \
./src/astro/eq2ec.o \
./src/astro/hapke.o \
./src/astro/jdcal.o \
./src/astro/kepler.o 

CU_DEPS += \
./src/astro/\ hapke_cuda.d 

C_DEPS += \
./src/astro/ec2eq.d \
./src/astro/eq2ec.d \
./src/astro/hapke.d \
./src/astro/jdcal.d \
./src/astro/kepler.d 


# Each subdirectory must supply rules for building sources it contributes
src/astro/\ hapke_cuda.o: ../src/astro/\ hapke_cuda.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-8.0/bin/nvcc -I/home/matt/git/cfitsio -G -g -lineinfo -pg -O0 -gencode arch=compute_35,code=sm_35 -m64 -odir "src/astro" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-8.0/bin/nvcc -I/home/matt/git/cfitsio -G -g -lineinfo -pg -O0 --compile --relocatable-device-code=true -gencode arch=compute_35,code=compute_35 -gencode arch=compute_35,code=sm_35 -m64  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/astro/%.o: ../src/astro/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-8.0/bin/nvcc -I/home/matt/git/cfitsio -G -g -lineinfo -pg -O0 -gencode arch=compute_35,code=sm_35 -m64 -odir "src/astro" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-8.0/bin/nvcc -I/home/matt/git/cfitsio -G -g -lineinfo -pg -O0 --compile -m64  -x c -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


