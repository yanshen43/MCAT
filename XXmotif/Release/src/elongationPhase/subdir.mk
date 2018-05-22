################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/elongationPhase/Extension.cpp \
../src/elongationPhase/ExtensionTable.cpp \
../src/elongationPhase/Kmer.cpp \
../src/elongationPhase/elongationCandidates.cpp 

OBJS += \
./src/elongationPhase/Extension.o \
./src/elongationPhase/ExtensionTable.o \
./src/elongationPhase/Kmer.o \
./src/elongationPhase/elongationCandidates.o 

CPP_DEPS += \
./src/elongationPhase/Extension.d \
./src/elongationPhase/ExtensionTable.d \
./src/elongationPhase/Kmer.d \
./src/elongationPhase/elongationCandidates.d 


# Each subdirectory must supply rules for building sources it contributes
src/elongationPhase/%.o: ../src/elongationPhase/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -DLOG_MAX_LEVEL=0 -O3 -g3 -pedantic -pedantic-errors -Wall -Werror -c -fmessage-length=0 -fno-strict-aliasing -std=c++0x -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


