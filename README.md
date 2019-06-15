# SurfGenBound
Surface generation package for bounded system

SurfGenBound casts the elements of diabatic Hamiltonian into simple polynomials of the internal coordinate

This version requires a reference point, whose ground state energy becomes the energy zero point, internal coordinate becomes the origin

Utilities:
1. Construct diabatic Hamiltonian from least square fitting ab initio data
* Energy, energy gradient, interstate coupling are adopted
* Various optimization technique is provided. The default choice is a good starting point, although you may want to try out other recipes to polish your surface
2. Analyze the fitted surface
* Evaluate H & ▽H in diabatic & adiabatic & nondegenerate representation on specified geometries
* Search for minimum on specified adiabatic or diabatic state, then analyze vibration
* Search for mex between specified adiabatic states, then orthogonalize gh & prepare input along gh to evaluate & draw double cone
* Shift origin
3. Generate input file for NadVibS, a nonadiabatic vibrational spectrum simulation package

Dependency: my Fortran-Library, as written in MyLib of makefile