################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/getopt_pp/getopt_pp.cpp 

OBJS += \
./src/getopt_pp/getopt_pp.o 

CPP_DEPS += \
./src/getopt_pp/getopt_pp.d 


# Each subdirectory must supply rules for building sources it contributes
src/getopt_pp/%.o: ../src/getopt_pp/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -DLOG_MAX_LEVEL=0 -O3 -g3 -pedantic -pedantic-errors -Wall -Werror -c -fmessage-length=0 -fno-strict-aliasing -std=c++0x -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


