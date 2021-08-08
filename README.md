# BEM
Matlab code of a project done for the master course: Boundary element methods

For description of code see here:

function mbp = add_fourier_shape_parametrisation(mbp, vs, no): implements parametrization in fourier domain

function dcd = interpolate_cauchy_data(mbp, ns, u, dnu): interpolates the cauchy data (Neumann and Dirichlet data onto the meshes contained in mbp)

function us = potential_cauchy(dcd, kappa, xs): Calculates values of solution of equation with the help of the cauchy data in the points xs

function dcd = nystroem_cauchy(mbp, ns, g, h, kappa, btf): implements discretization of boundary integral operator to calculate the missing boundary data

Then there are similar functions, where brakhage werner was used. These are all denoted with bw in their name.

Main 4-8 are all tasks we had to do for the assignement.

For the full project description in german, see also: http://cm.dmi.unibas.ch/teaching/bem/projekt.pdf
  
