################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/em/EM.cpp \
../src/em/hoNullModel.cpp 

OBJS += \
./src/em/EM.o \
./src/em/hoNullModel.o 

CPP_DEPS += \
./src/em/EM.d \
./src/em/hoNullModel.d 


# Each subdirectory must supply rules for building sources it contributes
src/em/%.o: ../src/em/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -DLOG_MAX_LEVEL=0 -O3 -g3 -pedantic -pedantic-errors -Wall -Werror -c -fmessage-length=0 -fno-strict-aliasing -std=c++0x -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


