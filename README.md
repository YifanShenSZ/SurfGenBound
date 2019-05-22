# SurfGenBound
Surface generation package for bounded system

SurfGenBound casts the elements of diabatic Hamiltonian into simple polynomials of the internal coordinate

This version requires a reference point, whose ground state energy becomes the energy zero point, internal coordinate becomes the origin

Utilities:
1. Construct diabatic Hamiltonian from least square fitting ab initio data
* Energy, energy gradient, interstate coupling are adopted
* Various optimization technique is provided. The default choice is a good starting point, although you may want to try out other recipes to polish your surface
2. Analyze the fitted surface
* Evaluate the surface on several input geometries, including energy, energy gradient, nonadiabatic coupling, diabatic Hamiltonian matrix, etc
* Search for minimum on specified adiabatic surface
* Search for minimum energy crossing point between specified adiabatic states
* Shift origin
3. Generate input file for NadVibS, a nonadiabatic vibrational spectrum simulation package

Dependency: my Fortran-Library, as written in MyLib of makefile