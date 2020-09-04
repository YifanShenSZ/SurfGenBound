# SurfGenBound
Surface generation package for bounded system

SurfGenBound casts the elements of diabatic Hamiltonian into simple polynomials of the internal coordinate

A reference point is required, whose ground state energy becomes the energy zero point, internal coordinate becomes the origin

## Featured utilities
Construct diabatic Hamiltonian from least square fitting ab initio data:
* Energy, energy gradient, interstate coupling are adopted
* Various optimization technique is provided. The default choice is a good starting point, although you may want to try out other recipes to polish your surface

Analyze the fitted surface:
* Evaluate H & ▽H in diabatic & adiabatic & nondegenerate representation on specified geometries
* Search for minimum on specified adiabatic or diabatic state, then analyze vibration
* Search for mex between specified adiabatic states, then orthogonalize gh & prepare input along gh to evaluate & draw double cone
* Shift origin

Generate input file for NadVibS, a nonadiabatic vibrational spectrum simulation package

## Installation
1. Copy mkl_rci.f90 & mkl_dfti.f90 from MKL installation path to Fortran-Library_v1.0.0 (for gnu compiler, additionally modify DJACOBI in mkl_rci.f90 to specify the dimension of x, fjac, f explicitly: x(n), fjac(m,n), f(m))
2. `make`

## Source
To see what this package is capable of in detail, you may open certain source file then fold all: routines are categorized and folded into different sections (VS code is recommended: press ctrl+k+0)

Source code level from bottom to top:
1. Fortran-Library_v1.0.0 (see https://github.com/YifanShenSZ/Fortran-Library for more details)
2. DiabaticHamiltonian
3. Basic
4. HdLeastSquareFit, Analyzation, ElectronicStructure
5. NadVibSInterface
6. Main

## Cite this work
> 1. Y. Shen and D. R. Yarkony, J. Phys. Chem. A 2020, 124, 22, 4539–4548 https://doi.org/10.1021/acs.jpca.0c02763
> 2. Y. Shen and D. R. Yarkony, J. Phys. Chem. Lett. 2020, 11, 17, 7245–7252 https://doi.org/10.1021/acs.jpclett.0c02199