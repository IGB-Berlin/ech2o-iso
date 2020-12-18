################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Atmosphere/AdjustPrecip.cpp \
../Atmosphere/AdvanceClimateMaps.cpp \
../Atmosphere/CountNumZones.cpp \
../Atmosphere/InitiateClimateMap.cpp \
../Atmosphere/UpdateClimateMap.cpp 

OBJS += \
./Atmosphere/AdjustPrecip.o \
./Atmosphere/AdvanceClimateMaps.o \
./Atmosphere/CountNumZones.o \
./Atmosphere/InitiateClimateMap.o \
./Atmosphere/UpdateClimateMap.o 

CPP_DEPS += \
./Atmosphere/AdjustPrecip.d \
./Atmosphere/AdvanceClimateMaps.d \
./Atmosphere/CountNumZones.d \
./Atmosphere/InitiateClimateMap.d \
./Atmosphere/UpdateClimateMap.d 


# Each subdirectory must supply rules for building sources it contributes
Atmosphere/%.o: ../Atmosphere/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -ggdb -DCPU_LITTLE_ENDIAN -I"../includes" -O3 -ggdb -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


