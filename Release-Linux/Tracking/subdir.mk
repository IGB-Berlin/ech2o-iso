################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Tracking/AdvanceBCMaps_Iso.cpp \
../Tracking/BCupstreamMixing.cpp \
../Tracking/CalcInitTPD.cpp \
../Tracking/CalcTPDtoLayers.cpp \
../Tracking/CalcTrck_L1L2.cpp \
../Tracking/CalcTrck_SoilAv.cpp \
../Tracking/CheckMapsTrck.cpp \
../Tracking/FCdownstream.cpp \
../Tracking/Fractionation_Estorage.cpp \
../Tracking/IncrementAge.cpp \
../Tracking/InitiateBCMap_Iso.cpp \
../Tracking/MixingTPD_postET.cpp \
../Tracking/MixingV_canopy.cpp \
../Tracking/MixingV_down.cpp \
../Tracking/MixingV_evapS.cpp \
../Tracking/MixingV_evapW.cpp \
../Tracking/MixingV_latup.cpp \
../Tracking/MixingV_sapflow.cpp \
../Tracking/MixingV_snow.cpp \
../Tracking/MixingV_through.cpp \
../Tracking/OutletVals.cpp \
../Tracking/TracerMixing.cpp \
../Tracking/UpdateBCMap_Iso.cpp

OBJS += \
./Tracking/AdvanceBCMaps_Iso.o \
./Tracking/BCupstreamMixing.o \
./Tracking/CalcInitTPD.o \
./Tracking/CalcTPDtoLayers.o \
./Tracking/CalcTrck_L1L2.o \
./Tracking/CalcTrck_SoilAv.o \
./Tracking/CheckMapsTrck.o\
./Tracking/FCdownstream.o \
./Tracking/Fractionation_Estorage.o \
./Tracking/IncrementAge.o \
./Tracking/InitiateBCMap_Iso.o \
./Tracking/MixingTPD_postET.o \
./Tracking/MixingV_canopy.o \
./Tracking/MixingV_down.o \
./Tracking/MixingV_evapS.o \
./Tracking/MixingV_evapW.o \
./Tracking/MixingV_latup.o \
./Tracking/MixingV_sapflow.o \
./Tracking/MixingV_snow.o \
./Tracking/MixingV_through.o \
./Tracking/OutletVals.o \
./Tracking/TracerMixing.o \
./Tracking/UpdateBCMap_Iso.o

CPP_DEPS += \
./Tracking/AdvanceBCMaps_Iso.d \
./Tracking/BCupstreamMixing.d \
./Tracking/CalcInitTPD.d \
./Tracking/CalcTPDtoLayers.d \
./Tracking/CalcTrck_L1L2.d \
./Tracking/CalcTrck_SoilAv.d \
./Tracking/CheckMapsTrck.d \
./Tracking/FCdownstream.d \
./Tracking/Fractionation_Estorage.d \
./Tracking/IncrementAge.d \
./Tracking/InitiateBCMap_Iso.d \
./Tracking/MixingTPD_postET.d \
./Tracking/MixingV_canopy.d \
./Tracking/MixingV_down.d \
./Tracking/MixingV_evapS.d \
./Tracking/MixingV_evapW.d \
./Tracking/MixingV_latup.d \
./Tracking/MixingV_sapflow.d \
./Tracking/MixingV_snow.d \
./Tracking/MixingV_through.d \
./Tracking/OutletVals.d \
./Tracking/TracerMixing.d \
./Tracking/UpdateBCMap_Iso.d

# Each subdirectory must supply rules for building sources it contributes
Tracking/%.o: ../Tracking/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -ggdb -DCPU_LITTLE_ENDIAN -I"../includes" -O3 -ggdb -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


