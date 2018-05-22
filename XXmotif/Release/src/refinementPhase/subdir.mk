################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/refinementPhase/Compare_Motifs.cpp \
../src/refinementPhase/Iterate_Motif.cpp \
../src/refinementPhase/Motif.cpp \
../src/refinementPhase/MotifContainer.cpp \
../src/refinementPhase/Refine_Motifs.cpp \
../src/refinementPhase/StartPosUpdater.cpp \
../src/refinementPhase/empiricalPvalCalibration.cpp 

OBJS += \
./src/refinementPhase/Compare_Motifs.o \
./src/refinementPhase/Iterate_Motif.o \
./src/refinementPhase/Motif.o \
./src/refinementPhase/MotifContainer.o \
./src/refinementPhase/Refine_Motifs.o \
./src/refinementPhase/StartPosUpdater.o \
./src/refinementPhase/empiricalPvalCalibration.o 

CPP_DEPS += \
./src/refinementPhase/Compare_Motifs.d \
./src/refinementPhase/Iterate_Motif.d \
./src/refinementPhase/Motif.d \
./src/refinementPhase/MotifContainer.d \
./src/refinementPhase/Refine_Motifs.d \
./src/refinementPhase/StartPosUpdater.d \
./src/refinementPhase/empiricalPvalCalibration.d 


# Each subdirectory must supply rules for building sources it contributes
src/refinementPhase/%.o: ../src/refinementPhase/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -DLOG_MAX_LEVEL=0 -O3 -g3 -pedantic -pedantic-errors -Wall -Werror -c -fmessage-length=0 -fno-strict-aliasing -std=c++0x -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


