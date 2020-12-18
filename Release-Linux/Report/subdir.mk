################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Report/RenameFile.cpp \
../Report/ReportTimeSeries.cpp \
../Report/ReportVectCells.cpp 

OBJS += \
./Report/RenameFile.o \
./Report/ReportTimeSeries.o \
./Report/ReportVectCells.o 

CPP_DEPS += \
./Report/RenameFile.d \
./Report/ReportTimeSeries.d \
./Report/ReportVectCells.d 


# Each subdirectory must supply rules for building sources it contributes
Report/%.o: ../Report/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -ggdb -DCPU_LITTLE_ENDIAN -I"../includes" -O3 -ggdb -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


