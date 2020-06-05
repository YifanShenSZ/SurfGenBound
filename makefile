###########################################
#                                         #
#        Makefile for SurfGenBound        #
#                                         #
###########################################

# Flags for user to tune
compiler = ifort
flag = -mkl -m64 -xCORE-AVX2 -mtune=core-avx2 -no-prec-div -fp-model fast=2 -parallel -O3

# User does not have to take care of following variables
libSrc = $(addprefix Fortran-Library_v1.0.0/, General.f90 Mathematics.f90 LinearAlgebra.f90 \
mkl_rci.f90 NonlinearOptimization.f90 Chemistry.f90 \
GeometryTransformation.f90)
src = $(addprefix source/, DiabaticHamiltonian.f90 \
Basic.f90 \
HdLeastSquareFit.f90 Analyzation.f90 ElectronicStructure.f90 \
NadVibSInterface.f90 \
Main.f90)

# release
SurfGenBound.exe: $(libSrc) $(src)
	$(compiler) $(flag) -ipo $^ -o $@

# debug
debug.exe: DiabaticHamiltonian.o Basic.o HdLeastSquareFit.o Analyzation.o ElectronicStructure.o NadVibSInterface.o Main.o
	$(compiler) $(flag) -lFL $^ -o $@

%.o: source/%.f90
	$(compiler) $(flag) -I~/Library/Fortran-Library/include -c $<

.PHONY: clean
clean:
	rm *.mod *.o
