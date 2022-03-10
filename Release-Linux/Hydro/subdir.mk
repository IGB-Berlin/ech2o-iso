################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Hydro/AdvanceBCMaps.cpp \
../Hydro/AdvanceLAIMaps.cpp \
../Hydro/CalcCatchArea.cpp \
../Hydro/CalcFracMobileWater.cpp \
../Hydro/CalcInitialStreamStorage.cpp \
../Hydro/CalcPropLayers.cpp \
../Hydro/CalcRootDistrib.cpp \
../Hydro/CalcSoilMoistureProfile.cpp \
../Hydro/CalcSoilResist.cpp \
../Hydro/CalcTPDMoisture.cpp \
../Hydro/CalculateForestGrowth.cpp \
../Hydro/CalculateSatArea.cpp \
../Hydro/CanopyInterception.cpp \
../Hydro/ChannelEvaporation.cpp \
../Hydro/GWrouting.cpp \
../Hydro/Infilt_GreenAmpt.cpp \
../Hydro/InitiateBCMap.cpp \
../Hydro/KinematicWave.cpp \
../Hydro/RichardsEquation.cpp \
../Hydro/SnowOutputPhase.cpp \
../Hydro/SoilEvapotranspiration.cpp \
../Hydro/SoilWaterRedistribution.cpp \
../Hydro/SolveCanopyFluxes.cpp \
../Hydro/SolveSurfaceEnergyBalance.cpp \
../Hydro/UpdateBCMap.cpp \
../Hydro/UpdateSnowPack.cpp 

OBJS += \
./Hydro/AdvanceBCMaps.o \
./Hydro/AdvanceLAIMaps.o \
./Hydro/CalcCatchArea.o \
./Hydro/CalcFracMobileWater.o \
./Hydro/CalcInitialStreamStorage.o \
./Hydro/CalcPropLayers.o \
./Hydro/CalcRootDistrib.o \
./Hydro/CalcSoilMoistureProfile.o \
./Hydro/CalcSoilResist.o \
./Hydro/CalcTPDMoisture.o \
./Hydro/CalculateForestGrowth.o \
./Hydro/CalculateSatArea.o \
./Hydro/CanopyInterception.o \
./Hydro/ChannelEvaporation.o \
./Hydro/GWrouting.o \
./Hydro/Infilt_GreenAmpt.o \
./Hydro/InitiateBCMap.o \
./Hydro/KinematicWave.o \
./Hydro/RichardsEquation.o \
./Hydro/SnowOutputPhase.o \
./Hydro/SoilEvapotranspiration.o \
./Hydro/SoilWaterRedistribution.o \
./Hydro/SolveCanopyFluxes.o \
./Hydro/SolveSurfaceEnergyBalance.o \
./Hydro/UpdateBCMap.o \
./Hydro/UpdateSnowPack.o 

CPP_DEPS += \
./Hydro/AdvanceBCMaps.d \
./Hydro/AdvanceLAIMaps.d \
./Hydro/CalcCatchArea.d \
./Hydro/CalcFracMobileWater.d \
./Hydro/CalcInitialStreamStorage.d \
./Hydro/CalcPropLayers.d \
./Hydro/CalcRootDistrib.d \
./Hydro/CalcSoilMoistureProfile.d \
./Hydro/CalcSoilResist.d \
./Hydro/CalcTPDMoisture.d \
./Hydro/CalculateForestGrowth.d \
./Hydro/CalculateSatArea.d \
./Hydro/CanopyInterception.d \
./Hydro/ChannelEvaporation.d \
./Hydro/GWrouting.d \
./Hydro/Infilt_GreenAmpt.d \
./Hydro/InitiateBCMap.d \
./Hydro/KinematicWave.d \
./Hydro/RichardsEquation.d \
./Hydro/SnowOutputPhase.d \
./Hydro/SoilEvapotranspiration.d \
./Hydro/SoilWaterRedistribution.d \
./Hydro/SolveCanopyFluxes.d \
./Hydro/SolveSurfaceEnergyBalance.d \
./Hydro/UpdateBCMap.d \
./Hydro/UpdateSnowPack.d 


# Each subdirectory must supply rules for building sources it contributes
Hydro/%.o: ../Hydro/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -ggdb -DCPU_LITTLE_ENDIAN -I"../includes" -O3 -ggdb -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


