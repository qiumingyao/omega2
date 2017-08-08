################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Config.cpp \
../src/DataSet.cpp \
../src/Edge.cpp \
../src/OverlapGraph.cpp \
../src/Read.cpp \
../src/main.cpp 

OBJS += \
./src/Config.o \
./src/DataSet.o \
./src/Edge.o \
./src/OverlapGraph.o \
./src/Read.o \
./src/main.o 

CPP_DEPS += \
./src/Config.d \
./src/DataSet.d \
./src/Edge.d \
./src/OverlapGraph.d \
./src/Read.d \
./src/main.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -std=c++0x -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


