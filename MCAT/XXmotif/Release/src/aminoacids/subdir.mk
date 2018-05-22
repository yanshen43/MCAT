################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/aminoacids/DisoCons.cpp \
../src/aminoacids/MProGlobal.cpp \
../src/aminoacids/NNet.cpp \
../src/aminoacids/NNetSupplInfProvider.cpp \
../src/aminoacids/NNode.cpp \
../src/aminoacids/NullModelCompositional.cpp \
../src/aminoacids/StateLib.cpp \
../src/aminoacids/alphabet.cpp \
../src/aminoacids/columnState.cpp \
../src/aminoacids/madonaPro.cpp \
../src/aminoacids/utils.cpp 

OBJS += \
./src/aminoacids/DisoCons.o \
./src/aminoacids/MProGlobal.o \
./src/aminoacids/NNet.o \
./src/aminoacids/NNetSupplInfProvider.o \
./src/aminoacids/NNode.o \
./src/aminoacids/NullModelCompositional.o \
./src/aminoacids/StateLib.o \
./src/aminoacids/alphabet.o \
./src/aminoacids/columnState.o \
./src/aminoacids/madonaPro.o \
./src/aminoacids/utils.o 

CPP_DEPS += \
./src/aminoacids/DisoCons.d \
./src/aminoacids/MProGlobal.d \
./src/aminoacids/NNet.d \
./src/aminoacids/NNetSupplInfProvider.d \
./src/aminoacids/NNode.d \
./src/aminoacids/NullModelCompositional.d \
./src/aminoacids/StateLib.d \
./src/aminoacids/alphabet.d \
./src/aminoacids/columnState.d \
./src/aminoacids/madonaPro.d \
./src/aminoacids/utils.d 


# Each subdirectory must supply rules for building sources it contributes
src/aminoacids/%.o: ../src/aminoacids/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -DLOG_MAX_LEVEL=0 -O3 -g3 -pedantic -pedantic-errors -Wall -Werror -c -fmessage-length=0 -fno-strict-aliasing -std=c++0x -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


