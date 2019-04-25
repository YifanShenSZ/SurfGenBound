###########################################
#                                         #
#        Makefile for SurfGenBound        #
#                                         #
###########################################

compiler = ifort
# This is my library directory on MARCC
MyLibDir = /home-4/yshen57@jhu.edu/Library
MyLib = $(MyLibDir)/General.f90 $(MyLibDir)/Mathematics.f90 $(MyLibDir)/LinearAlgebra.f90 $(MyLibDir)/MKL_RCI.f90 $(MyLibDir)/NonlinearOptimization.f90 $(MyLibDir)/GeometryTransformation.f90 $(MyLibDir)/Nonadiabatic.f90
src = DiabaticHamiltonian.f90 Basic.f90 ElectronicStructure.f90 HdLeastSquareFit.f90 Analyzation.f90 NadVibS.f90 Main.f90
exe = SurfGenBound.exe
flag = -u -mkl -fast -march=core-avx2

$(exe): $(MyLib) $(src)
	$(compiler) $(flag) $^ -o $(exe)

clean:
	rm $(exe)
	rm *.mod