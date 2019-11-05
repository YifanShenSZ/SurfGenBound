###########################################
#                                         #
#        Makefile for SurfGenBound        #
#                                         #
###########################################

compiler = ifort
flag = -m64 -xavx -mtune=corei7-avx -mkl -O3 -no-prec-div -fp-model fast=2
MyLibSrcDir = ~/Library/Fortran-Library
MyLibSrc = $(MyLibSrcDir)/General.f90 $(MyLibSrcDir)/Mathematics.f90 $(MyLibSrcDir)/LinearAlgebra.f90 \
$(MyLibSrcDir)/mkl_rci.f90 $(MyLibSrcDir)/NonlinearOptimization.f90 $(MyLibSrcDir)/Chemistry.f90 \
$(MyLibSrcDir)/GeometryTransformation.f90
src = DiabaticHamiltonian.f90 Basic.f90 ElectronicStructure.f90 HdLeastSquareFit.f90 Analyzation.f90 NadVibSInterface.f90 Main.f90

SurfGenBound.exe: $(MyLibSrc) $(src)
	$(compiler) $(flag) -ipo $^ -o SurfGenBound.exe

.PHONY: clean
clean:
	rm -f *.mod
