################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/AbstractKmer.cpp \
../src/Globals.cpp \
../src/LogTable.cpp \
../src/NullModel.cpp \
../src/ProgressBar.cpp \
../src/SmallKmer.cpp \
../src/StartModels.cpp \
../src/ThresholdChecker.cpp \
../src/UngappedKmer.cpp \
../src/UniversalKmer.cpp \
../src/alphabet.cpp \
../src/backgroundDistribution.cpp \
../src/branch_and_bound.cpp \
../src/log.cpp \
../src/main.cpp \
../src/output.cpp \
../src/pValCalculation.cpp \
../src/seqset.cpp 

OBJS += \
./src/AbstractKmer.o \
./src/Globals.o \
./src/LogTable.o \
./src/NullModel.o \
./src/ProgressBar.o \
./src/SmallKmer.o \
./src/StartModels.o \
./src/ThresholdChecker.o \
./src/UngappedKmer.o \
./src/UniversalKmer.o \
./src/alphabet.o \
./src/backgroundDistribution.o \
./src/branch_and_bound.o \
./src/log.o \
./src/main.o \
./src/output.o \
./src/pValCalculation.o \
./src/seqset.o 

CPP_DEPS += \
./src/AbstractKmer.d \
./src/Globals.d \
./src/LogTable.d \
./src/NullModel.d \
./src/ProgressBar.d \
./src/SmallKmer.d \
./src/StartModels.d \
./src/ThresholdChecker.d \
./src/UngappedKmer.d \
./src/UniversalKmer.d \
./src/alphabet.d \
./src/backgroundDistribution.d \
./src/branch_and_bound.d \
./src/log.d \
./src/main.d \
./src/output.d \
./src/pValCalculation.d \
./src/seqset.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -DLOG_MAX_LEVEL=0 -O3 -g3 -pedantic -pedantic-errors -Wall -Werror -c -fmessage-length=0 -fno-strict-aliasing -std=c++0x -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


