###########################################
#                                         #
#        Makefile for SurfGenBound        #
#                                         #
###########################################

compiler = ifort
flag = -m64 -xCORE-AVX2 -mtune=core-avx2 -no-prec-div -fp-model fast=2 -mkl -parallel -O3
MyLibSrc = $(addprefix Fortran-Library_v1.0.0/, General.f90 Mathematics.f90 LinearAlgebra.f90 \
mkl_rci.f90 NonlinearOptimization.f90 Chemistry.f90 \
GeometryTransformation.f90)
src = $(addprefix source/, DiabaticHamiltonian.f90 \
Basic.f90 \
HdLeastSquareFit.f90 Analyzation.f90 ElectronicStructure.f90 \
NadVibSInterface.f90 \
Main.f90)
# Faster compilation for debugging
# Link to dynamic library version of Fortran-Library for fast debugging
MyInc = ~/Library/Fortran-Library/include
MyLib = ~/Library/Fortran-Library/lib
# Keep old object files
obj = DiabaticHamiltonian.o Basic.o ElectronicStructure.o HdLeastSquareFit.o Analyzation.o NadVibSInterface.o Main.o

SurfGenBound.exe: $(MyLibSrc) $(src)
	$(compiler) $(flag) -ipo $^ -o SurfGenBound.exe

debug: $(obj)
	$(compiler) $(flag) -I$(MyInc) -L$(MyLib) -lFL $^ -o debug

DiabaticHamiltonian.o: DiabaticHamiltonian.f90
	$(compiler) $(flag) -I$(MyInc) -c $^

Basic.o: Basic.f90
	$(compiler) $(flag) -I$(MyInc) -c $^

ElectronicStructure.o: ElectronicStructure.f90
	$(compiler) $(flag) -I$(MyInc) -c $^

HdLeastSquareFit.o: HdLeastSquareFit.f90
	$(compiler) $(flag) -I$(MyInc) -c $^

Analyzation.o: Analyzation.f90
	$(compiler) $(flag) -I$(MyInc) -c $^

NadVibSInterface.o: NadVibSInterface.f90
	$(compiler) $(flag) -I$(MyInc) -c $^

Main.o: Main.f90
	$(compiler) $(flag) -I$(MyInc) -c $^

.PHONY: clean
clean:
	rm -f *.mod *.o
