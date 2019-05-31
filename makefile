###########################################
#                                         #
#        Makefile for SurfGenBound        #
#                                         #
###########################################

compiler = ifort
MyLibDir = /home-4/yshen57@jhu.edu/Library/Fortran-Library
MyLib = $(MyLibDir)/General.f90 $(MyLibDir)/Mathematics.f90 $(MyLibDir)/LinearAlgebra.f90 $(MyLibDir)/MKL_RCI.f90 $(MyLibDir)/NonlinearOptimization.f90 $(MyLibDir)/GeometryTransformation.f90 $(MyLibDir)/Chemistry.f90
src = DiabaticHamiltonian.f90 Basic.f90 ElectronicStructure.f90 HdLeastSquareFit.f90 Analyzation.f90 Main.f90
exe = SurfGenBound.exe
# -static causes strange bug in ifort 2018
flag = -u -mkl -O3 -ipo -no-prec-div -fp-model fast=2 -march=core-avx2 -mtune=core-avx2

$(exe): $(MyLib) $(src)
	$(compiler) $(flag) $^ -o $(exe)

clean:
	rm $(exe)
	rm *.mod