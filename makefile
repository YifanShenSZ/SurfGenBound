###########################################
#                                         #
#        Makefile for SurfGenBound        #
#                                         #
###########################################

compiler = ifort
MyLibDir = ~/Library/Fortran-Library
MyLib = $(MyLibDir)/General.f90 $(MyLibDir)/Mathematics.f90 $(MyLibDir)/LinearAlgebra.f90 \
$(MyLibDir)/mkl_rci.f90 $(MyLibDir)/NonlinearOptimization.f90 $(MyLibDir)/Chemistry.f90 \
$(MyLibDir)/GeometryTransformation.f90
src = DiabaticHamiltonian.f90 Basic.f90 ElectronicStructure.f90 HdLeastSquareFit.f90 Analyzation.f90 NadVibSInterface.f90 Main.f90
exe = SurfGenBound.exe
# -static causes strange bug in ifort 2018, so cannot use -fast
flag = -m64 -xCORE-AVX2 -mtune=core-avx2 -mkl -ipo -O3 -no-prec-div -fp-model fast=2

$(exe): $(MyLib) $(src)
	$(compiler) $(flag) $^ -o $(exe)

clean:
	rm *.mod