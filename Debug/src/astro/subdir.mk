################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/astro/ec2eq.c \
../src/astro/eq2ec.c \
../src/astro/hapke.c \
../src/astro/jdcal.c \
../src/astro/kepler.c 

OBJS += \
./src/astro/ec2eq.o \
./src/astro/eq2ec.o \
./src/astro/hapke.o \
./src/astro/jdcal.o \
./src/astro/kepler.o 

C_DEPS += \
./src/astro/ec2eq.d \
./src/astro/eq2ec.d \
./src/astro/hapke.d \
./src/astro/jdcal.d \
./src/astro/kepler.d 


# Each subdirectory must supply rules for building sources it contributes
src/astro/%.o: ../src/astro/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-9.1/bin/nvcc -G -g -O0 -gencode arch=compute_61,code=sm_61 -m64 -odir "src/astro" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-9.1/bin/nvcc -G -g -O0 --compile -m64  -x c -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


