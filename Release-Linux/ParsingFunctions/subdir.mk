################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../ParsingFunctions/ConfigFile.cpp \
../ParsingFunctions/FileRoot.cpp \
../ParsingFunctions/ParseString.cpp 

OBJS += \
./ParsingFunctions/ConfigFile.o \
./ParsingFunctions/FileRoot.o \
./ParsingFunctions/ParseString.o 

CPP_DEPS += \
./ParsingFunctions/ConfigFile.d \
./ParsingFunctions/FileRoot.d \
./ParsingFunctions/ParseString.d 


# Each subdirectory must supply rules for building sources it contributes
ParsingFunctions/%.o: ../ParsingFunctions/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -ggdb -DCPU_LITTLE_ENDIAN -I"../includes" -O3 -ggdb -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


