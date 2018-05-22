################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/seqFormat/AlignStream.cpp \
../src/seqFormat/Alignment.cpp \
../src/seqFormat/Clustalw.cpp \
../src/seqFormat/CustomMultFasta.cpp \
../src/seqFormat/Fasta.cpp \
../src/seqFormat/SeqExceptions.cpp 

OBJS += \
./src/seqFormat/AlignStream.o \
./src/seqFormat/Alignment.o \
./src/seqFormat/Clustalw.o \
./src/seqFormat/CustomMultFasta.o \
./src/seqFormat/Fasta.o \
./src/seqFormat/SeqExceptions.o 

CPP_DEPS += \
./src/seqFormat/AlignStream.d \
./src/seqFormat/Alignment.d \
./src/seqFormat/Clustalw.d \
./src/seqFormat/CustomMultFasta.d \
./src/seqFormat/Fasta.d \
./src/seqFormat/SeqExceptions.d 


# Each subdirectory must supply rules for building sources it contributes
src/seqFormat/%.o: ../src/seqFormat/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -DLOG_MAX_LEVEL=0 -O3 -g3 -pedantic -pedantic-errors -Wall -Werror -c -fmessage-length=0 -fno-strict-aliasing -std=c++0x -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


