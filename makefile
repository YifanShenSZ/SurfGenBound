###########################################
#                                         #
#        Makefile for SurfGenBound        #
#                                         #
###########################################

# Flags for user to tune
compiler = ifort
flag = -mkl -m64 -xCORE-AVX2 -mtune=core-avx2 -no-prec-div -fp-model fast=2 -parallel -O3 -ipo

# User does not have to take care of following variables
libSrc = $(addprefix Fortran-Library_v1.0.0/, General.f90 Mathematics.f90 LinearAlgebra.f90 \
mkl_rci.f90 NonlinearOptimization.f90 Chemistry.f90 \
GeometryTransformation.f90)
src = $(addprefix source/, DiabaticHamiltonian.f90 \
Basic.f90 \
HdLeastSquareFit.f90 Analyzation.f90 ElectronicStructure.f90 \
NadVibSInterface.f90 \
Main.f90)
# Faster compilation for debugging
inc = ~/Library/Fortran-Library/include
lib = ~/Library/Fortran-Library/lib
obj = DiabaticHamiltonian.o \
Basic.o \
HdLeastSquareFit.o Analyzation.o ElectronicStructure.o \
NadVibSInterface.o \
Main.o

# release
SurfGenBound.exe: $(libSrc) $(src)
	$(compiler) -ipo $(flag) $^ -o $@

# debug
debug.exe: $(obj)
	$(compiler) $(flag) -I$(inc) -L$(lib) -lFL $^ -o $@

DiabaticHamiltonian.o: source/DiabaticHamiltonian.f90
	$(compiler) $(flag) -I$(inc) -c $^

Basic.o: source/Basic.f90
	$(compiler) $(flag) -I$(inc) -c $^

ElectronicStructure.o: source/ElectronicStructure.f90
	$(compiler) $(flag) -I$(inc) -c $^

HdLeastSquareFit.o: source/HdLeastSquareFit.f90
	$(compiler) $(flag) -I$(inc) -c $^

Analyzation.o: source/Analyzation.f90
	$(compiler) $(flag) -I$(inc) -c $^

NadVibSInterface.o: source/NadVibSInterface.f90
	$(compiler) $(flag) -I$(inc) -c $^

Main.o: source/Main.f90
	$(compiler) $(flag) -I$(inc) -c $^

.PHONY: clean
clean:
	rm *.mod *.o
