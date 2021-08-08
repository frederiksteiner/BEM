# BEM
Matlab code of a project done for the master course: Boundary element methods

For description of code see here:

_function mbp = add_fourier_shape_parametrisation(mbp, vs, no)_: implements parametrization in fourier domain

_function dcd = interpolate_cauchy_data(mbp, ns, u, dnu)_: interpolates the cauchy data (Neumann and Dirichlet data onto the meshes contained in mbp)

_function us = potential_cauchy(dcd, kappa, xs)_: Calculates values of solution of equation with the help of the cauchy data in the points xs

_function dcd = nystroem_cauchy(mbp, ns, g, h, kappa, btf)_: implements discretization of boundary integral operator to calculate the missing boundary data

Then there are similar functions, where brakhage werner was used. These are all denoted with _bw_ in their name.

Main 4-8 are all tasks we had to do for the assignement.

For the full project description in german, see also: http://cm.dmi.unibas.ch/teaching/bem/projekt.pdf
  
