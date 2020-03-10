# SurfGenBound
H^d construction of cyclopentoxy

The fit is carried out by:
1. Get a quadric expansion:
* Start fitting on min
* Continue on Fewest-Anion
2. Move on to selected 4th + others 2nd order expansion:
* Based on the quadric expansion and far scan data, generate prior estimation
* Try out different regularization parameter on Fewest

At each step, trust region is applied first, then continue with other line search method to polish the fit if necessary

In practise, Dai-Yun conjugate gradient under Wolfe condition yields the lowest loss

Compiler dependence:
* This work is done with ifort 2018.3.222
* ifort 2017.7.259 and ifort 2019.4 would yield different result