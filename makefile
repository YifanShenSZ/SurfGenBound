###########################################
#                                         #
#        Makefile for SurfGenBound        #
#                                         #
###########################################

compiler = ifort
flag = -m64 -xCORE-AVX2 -mtune=core-avx2 -mkl -O3 -no-prec-div -fp-model fast=2
MyLibSrcDir = ~/Library/Fortran-Library
MyLibSrc = $(MyLibSrcDir)/General.f90 $(MyLibSrcDir)/Mathematics.f90 $(MyLibSrcDir)/LinearAlgebra.f90 \
$(MyLibSrcDir)/mkl_rci.f90 $(MyLibSrcDir)/NonlinearOptimization.f90 $(MyLibSrcDir)/Chemistry.f90 \
$(MyLibSrcDir)/GeometryTransformation.f90
src = DiabaticHamiltonian.f90 Basic.f90 ElectronicStructure.f90 HdLeastSquareFit.f90 Analyzation.f90 NadVibSInterface.f90 Main.f90
# Faster compilation for debugging
# Link to dynamic library version of Fortran-Library for fast debugging
MyInc = /usr/local/include
MyLib = /usr/local/lib
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
