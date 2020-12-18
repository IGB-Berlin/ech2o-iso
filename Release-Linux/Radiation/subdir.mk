################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Radiation/CalcAerodynResist.cpp \
../Radiation/GrndHeat.cpp \
../Radiation/LatHeat.cpp \
../Radiation/MeltHeat.cpp \
../Radiation/NetRad.cpp \
../Radiation/NetRad_water.cpp \
../Radiation/RainHeat.cpp \
../Radiation/SensHeat.cpp \
../Radiation/SolveSurfaceFluxes.cpp \
../Radiation/snowheat.cpp 

OBJS += \
./Radiation/CalcAerodynResist.o \
./Radiation/GrndHeat.o \
./Radiation/LatHeat.o \
./Radiation/MeltHeat.o \
./Radiation/NetRad.o \
./Radiation/NetRad_water.o \
./Radiation/RainHeat.o \
./Radiation/SensHeat.o \
./Radiation/SolveSurfaceFluxes.o \
./Radiation/snowheat.o 

CPP_DEPS += \
./Radiation/CalcAerodynResist.d \
./Radiation/GrndHeat.d \
./Radiation/LatHeat.d \
./Radiation/MeltHeat.d \
./Radiation/NetRad.d \
./Radiation/NetRad_water.d \
./Radiation/RainHeat.d \
./Radiation/SensHeat.d \
./Radiation/SolveSurfaceFluxes.d \
./Radiation/snowheat.d 


# Each subdirectory must supply rules for building sources it contributes
Radiation/%.o: ../Radiation/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -ggdb -DCPU_LITTLE_ENDIAN -I"../includes" -O3 -ggdb -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


