################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CU_SRCS += \
../src/util/cuda-util.cu \
../src/util/euler_cuda.cu 

C_SRCS += \
../src/util/addsuffix.c \
../src/util/allwhite.c \
../src/util/changepath.c \
../src/util/checkposet.c \
../src/util/clrmat.c \
../src/util/clrvect.c \
../src/util/cotrans.c \
../src/util/countdata.c \
../src/util/createdir.c \
../src/util/cross.c \
../src/util/distance.c \
../src/util/dot.c \
../src/util/euler.c \
../src/util/getdouble.c \
../src/util/getint.c \
../src/util/gettstr.c \
../src/util/intifpossible.c \
../src/util/iround.c \
../src/util/lowcase.c \
../src/util/matinv.c \
../src/util/mmmul.c \
../src/util/mtrnsps.c \
../src/util/normalize.c \
../src/util/readline.c \
../src/util/resampim.c \
../src/util/swapendian.c \
../src/util/tempname.c \
../src/util/timestring.c \
../src/util/ucvector.c \
../src/util/vecnorm.c 

OBJS += \
./src/util/addsuffix.o \
./src/util/allwhite.o \
./src/util/changepath.o \
./src/util/checkposet.o \
./src/util/clrmat.o \
./src/util/clrvect.o \
./src/util/cotrans.o \
./src/util/countdata.o \
./src/util/createdir.o \
./src/util/cross.o \
./src/util/cuda-util.o \
./src/util/distance.o \
./src/util/dot.o \
./src/util/euler.o \
./src/util/euler_cuda.o \
./src/util/getdouble.o \
./src/util/getint.o \
./src/util/gettstr.o \
./src/util/intifpossible.o \
./src/util/iround.o \
./src/util/lowcase.o \
./src/util/matinv.o \
./src/util/mmmul.o \
./src/util/mtrnsps.o \
./src/util/normalize.o \
./src/util/readline.o \
./src/util/resampim.o \
./src/util/swapendian.o \
./src/util/tempname.o \
./src/util/timestring.o \
./src/util/ucvector.o \
./src/util/vecnorm.o 

CU_DEPS += \
./src/util/cuda-util.d \
./src/util/euler_cuda.d 

C_DEPS += \
./src/util/addsuffix.d \
./src/util/allwhite.d \
./src/util/changepath.d \
./src/util/checkposet.d \
./src/util/clrmat.d \
./src/util/clrvect.d \
./src/util/cotrans.d \
./src/util/countdata.d \
./src/util/createdir.d \
./src/util/cross.d \
./src/util/distance.d \
./src/util/dot.d \
./src/util/euler.d \
./src/util/getdouble.d \
./src/util/getint.d \
./src/util/gettstr.d \
./src/util/intifpossible.d \
./src/util/iround.d \
./src/util/lowcase.d \
./src/util/matinv.d \
./src/util/mmmul.d \
./src/util/mtrnsps.d \
./src/util/normalize.d \
./src/util/readline.d \
./src/util/resampim.d \
./src/util/swapendian.d \
./src/util/tempname.d \
./src/util/timestring.d \
./src/util/ucvector.d \
./src/util/vecnorm.d 


# Each subdirectory must supply rules for building sources it contributes
src/util/%.o: ../src/util/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-8.0/bin/nvcc -G -g -lineinfo -pg -O0 -gencode arch=compute_35,code=sm_35 -m64 -odir "src/util" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-8.0/bin/nvcc -G -g -lineinfo -pg -O0 --compile -m64  -x c -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/util/%.o: ../src/util/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-8.0/bin/nvcc -G -g -lineinfo -pg -O0 -gencode arch=compute_35,code=sm_35 -m64 -odir "src/util" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-8.0/bin/nvcc -G -g -lineinfo -pg -O0 --compile --relocatable-device-code=true -gencode arch=compute_35,code=compute_35 -gencode arch=compute_35,code=sm_35 -m64  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


