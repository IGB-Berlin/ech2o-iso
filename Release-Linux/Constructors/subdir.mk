################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Constructors/AtmosphConstruct.cpp \
../Constructors/BasinConstruct.cpp \
../Constructors/CheckMaps.cpp \
../Constructors/ForestConstruct.cpp \
../Constructors/GroveConstruct.cpp \
../Constructors/ReportConstruct.cpp \
../Constructors/checkForestDatabase.cpp \
../Constructors/TrackingConstruct.cpp 

OBJS += \
./Constructors/AtmosphConstruct.o \
./Constructors/BasinConstruct.o \
./Constructors/CheckMaps.o \
./Constructors/ForestConstruct.o \
./Constructors/GroveConstruct.o \
./Constructors/ReportConstruct.o \
./Constructors/checkForestDatabase.o \
./Constructors/TrackingConstruct.o 

CPP_DEPS += \
./Constructors/AtmosphConstruct.d \
./Constructors/BasinConstruct.d \
./Constructors/CheckMaps.d \
./Constructors/ForestConstruct.d \
./Constructors/GroveConstruct.d \
./Constructors/ReportConstruct.d \
./Constructors/checkForestDatabase.d \
./Constructors/TrackingConstruct.d


# Each subdirectory must supply rules for building sources it contributes
Constructors/%.o: ../Constructors/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -ggdb -DCPU_LITTLE_ENDIAN -I"../includes" -O3 -ggdb -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


