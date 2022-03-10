################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Budgets/AccountFluxes.cpp \
../Budgets/AccountBCFluxes.cpp \
../Budgets/AccountRelArea.cpp \
../Budgets/AccountStorages.cpp \
../Budgets/AccountTrckBCFluxes.cpp \
../Budgets/AccountTrckStorages.cpp \
../Budgets/MassBalanceError.cpp \
../Budgets/TotalBoundaryInflow.cpp \
../Budgets/TotalDeepGrndFlow.cpp \
../Budgets/TotalDeepGWtoChn.cpp \
../Budgets/TotalDeepGW.cpp \
../Budgets/TotalEvaporation.cpp \
../Budgets/TotalEvaporationS.cpp \
../Budgets/TotalEvaporationC.cpp \
../Budgets/TotalEvaporationI.cpp \
../Budgets/TotalGWtoChn.cpp \
../Budgets/TotalTranspiration.cpp \
../Budgets/TotalLeakage.cpp \
../Budgets/TotalOvlndFlow.cpp \
../Budgets/TotalPrecipitation.cpp \
../Budgets/TotalRecharge.cpp \
../Budgets/TotalSaturationArea.cpp \
../Budgets/TotalSrftoChn.cpp \
../Budgets/TotalStorage.cpp \
../Budgets/TotalGrndFlow.cpp \
../Budgets/TrckBalanceError.cpp

OBJS += \
./Budgets/AccountFluxes.o \
./Budgets/AccountBCFluxes.o \
./Budgets/AccountRelArea.o \
./Budgets/AccountStorages.o \
./Budgets/AccountTrckBCFluxes.o \
./Budgets/AccountTrckStorages.o \
./Budgets/MassBalanceError.o \
./Budgets/TotalBoundaryInflow.o \
./Budgets/TotalDeepGrndFlow.o \
./Budgets/TotalDeepGWtoChn.o \
./Budgets/TotalDeepGW.o \
./Budgets/TotalEvaporation.o \
./Budgets/TotalEvaporationS.o \
./Budgets/TotalEvaporationC.o \
./Budgets/TotalEvaporationI.o \
./Budgets/TotalGrndFlow.o \
./Budgets/TotalGWtoChn.o \
./Budgets/TotalTranspiration.o \
./Budgets/TotalLeakage.o \
./Budgets/TotalOvlndFlow.o \
./Budgets/TotalPrecipitation.o \
./Budgets/TotalRecharge.o \
./Budgets/TotalSaturationArea.o \
./Budgets/TotalSrftoChn.o \
./Budgets/TotalStorage.o \
./Budgets/TrckBalanceError.o

CPP_DEPS += \
./Budgets/AccountFluxes.d \
./Budgets/AccountBCFluxes.d \
./Budgets/AccountRelArea.d \
./Budgets/AccountStorages.d \
./Budgets/AccountTrckBCFluxes.d \
./Budgets/AccountTrckStorages.d \
./Budgets/MassBalanceError.d \
./Budgets/TotalBoundaryInflow.d \
./Budgets/TotalDeepGrndFlow.d \
./Budgets/TotalDeepGWtoChn.d \
./Budgets/TotalDeepGW.d \
./Budgets/TotalEvaporation.d \
./Budgets/TotalEvaporationS.d \
./Budgets/TotalEvaporationC.d \
./Budgets/TotalEvaporationI.d \
./Budgets/TotalGrndFlow.d \
./Budgets/TotalGWtoChn.d \
./Budgets/TotalTranspiration.d \
./Budgets/TotalLeakage.d \
./Budgets/TotalOvlndFlow.d \
./Budgets/TotalPrecipitation.d \
./Budgets/TotalRecharge.d \
./Budgets/TotalSaturationArea.d \
./Budgets/TotalSrftoChn.d \
./Budgets/TotalStorage.d \
./Budgets/TrckBalanceError.d


# Each subdirectory must supply rules for building sources it contributes
Budgets/%.o: ../Budgets/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -ggdb -DCPU_LITTLE_ENDIAN -I"../includes" -O3 -ggdb -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


