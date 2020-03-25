# SurfGenBound
H^d construction of cyclopentoxy

## Content
The contents are divided into 4 parts in corresponding directories:
1. Training set: MRCISD data
2. input: collection of all kinds of input files
3. fit: how this H^d is constructed, with all raw input and output
4. analyze: some analyzation of H^d, with all raw input and output

## Technical details
This fit is carried out by:
1. Get a quadric expansion:
* Start fitting on min
* Continue on Fewest-Anion
2. Move on to selected 4th + others 2nd order expansion:
* Based on the quadric expansion and far scan data, generate prior estimation
* Try out different regularization parameter on Fewest

At each step, trust region is applied first, then continue with other line search method to polish the fit if necessary

In practise, Dai-Yun conjugate gradient under Wolfe condition yields the lowest loss

Compiler and platform:
* The compiler adopted in this work is ifort 2018.3.222
* The operating systems used in this work are Ubuntu 18.04 and CentOS 7