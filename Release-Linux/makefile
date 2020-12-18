################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include wrappers/subdir.mk
-include Tracking/subdir.mk
-include Report/subdir.mk
-include Radiation/subdir.mk
-include ParsingFunctions/subdir.mk
-include IO/subdir.mk
-include Hydro/subdir.mk
-include Grid/subdir.mk
-include Functions/subdir.mk
-include Forest/subdir.mk
-include Destructors/subdir.mk
-include Control/subdir.mk
-include Constructors/subdir.mk
-include Budgets/subdir.mk
-include Atmosphere/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: ech2o_iso

# Tool invocations
ech2o_iso: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ -ggdb -L"../lib" -L"/usr/local/lib" -L"${DEST}/lib64" -L"${DEST}/lib" -L"${BOOSTDIR}/lib" -I"${DEST}/include" -I"${BOOSTDIR}/include/boost" -fopenmp -o "ech2o_iso" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '
	$(MAKE) --no-print-directory post-build

# Other Targets
clean:
	-$(RM) $(CC_DEPS)$(C++_DEPS)$(EXECUTABLES)$(C_UPPER_DEPS)$(CXX_DEPS)$(OBJS)$(CPP_DEPS)$(C_DEPS) ech2o_iso
	-@echo ' '

post-build:
#	-strip ech2o_iso
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY: post-build

-include ../makefile.targets 
