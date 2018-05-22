################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/nucleotides/conservationScores.cpp \
../src/nucleotides/motifRegion.cpp 

OBJS += \
./src/nucleotides/conservationScores.o \
./src/nucleotides/motifRegion.o 

CPP_DEPS += \
./src/nucleotides/conservationScores.d \
./src/nucleotides/motifRegion.d 


# Each subdirectory must supply rules for building sources it contributes
src/nucleotides/%.o: ../src/nucleotides/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -DLOG_MAX_LEVEL=0 -O3 -g3 -pedantic -pedantic-errors -Wall -Werror -c -fmessage-length=0 -fno-strict-aliasing -std=c++0x -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


