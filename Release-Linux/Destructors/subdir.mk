################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Destructors/AtmosphDesctruct.cpp \
../Destructors/BasinDestruct.cpp \
../Destructors/ForestDestruct.cpp \
../Destructors/GroveDestruct.cpp \
../Destructors/TrackingDestruct.cpp 

OBJS += \
./Destructors/AtmosphDesctruct.o \
./Destructors/BasinDestruct.o \
./Destructors/ForestDestruct.o \
./Destructors/GroveDestruct.o \
./Destructors/TrackingDestruct.o

CPP_DEPS += \
./Destructors/AtmosphDesctruct.d \
./Destructors/BasinDestruct.d \
./Destructors/ForestDestruct.d \
./Destructors/GroveDestruct.d \
./Destructors/TrackingDestruct.d


# Each subdirectory must supply rules for building sources it contributes
Destructors/%.o: ../Destructors/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -ggdb -DCPU_LITTLE_ENDIAN -I"../includes" -O3 -ggdb -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


