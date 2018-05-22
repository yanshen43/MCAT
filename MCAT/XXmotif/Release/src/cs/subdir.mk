################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../src/cs/alphabet.cc \
../src/cs/amino_acid.cc \
../src/cs/blast_hits.cc \
../src/cs/blosum_matrix.cc \
../src/cs/gonnet_matrix.cc \
../src/cs/log.cc 

OBJS += \
./src/cs/alphabet.o \
./src/cs/amino_acid.o \
./src/cs/blast_hits.o \
./src/cs/blosum_matrix.o \
./src/cs/gonnet_matrix.o \
./src/cs/log.o 

CC_DEPS += \
./src/cs/alphabet.d \
./src/cs/amino_acid.d \
./src/cs/blast_hits.d \
./src/cs/blosum_matrix.d \
./src/cs/gonnet_matrix.d \
./src/cs/log.d 


# Each subdirectory must supply rules for building sources it contributes
src/cs/%.o: ../src/cs/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -DLOG_MAX_LEVEL=0 -O3 -g3 -pedantic -pedantic-errors -Wall -Werror -c -fmessage-length=0 -fno-strict-aliasing -std=c++0x -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


