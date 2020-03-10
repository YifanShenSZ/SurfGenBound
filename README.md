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
`make`

## Source
To see what this package is capable of in detail, you may open certain source file then fold all: routines are categorized and folded into different sections (VS code is recommended: press ctrl+k+0)

Source code level from bottom to top:
1. DiabaticHamiltonian
2. Basic
3. HdLeastSquareFit, Analyzation, ElectronicStructure
4. NadVibSInterface
5. Main