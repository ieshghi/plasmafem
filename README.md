# plasmafem
Finite element method with finite elements for solveing Poisson equation and Grad-Shafranov equation, with 2nd order errors.

Currently can solve both of them on any boundary, with second order errors.

File list:

mgmres.f90 -sparse linear equation solver, downloaded from http://people.sc.fsu.edu/~jburkardt/f_src/mgmres/mgmres.html
mesh.f90 -contains a few utilities, like a linspace() routine, 3x3 matrix inverter, determinant calculator, and a mesh generator for a rectangle. Also contains the exact solution and right hand side of the equation we want to solve
meshgs.f90 -same as mesh.f90, but is then imported in the G-S version of the code, not the Poisson one.
grad.f90 -G-S finite element solver
pois.f90 -Poisson finite element solver
plot.py -Python routines for plotting, debugging, checking errors.

in files/:
ia.dat,ja.dat,arr.dat,fu.dat,p.dat,t.dat,x.dat : used by plot.py to analyse the latest finite element solution. Contain array info and solutions.
conv.dat : lists the maximum error for the latest few runs of the code. Each run of pois.f90 or grad.f90 appends the max error to the file, so the file has to be reinitialised for each attempt at testing convergence.
exact.dat : contains exact solution for the latest run of pois.f90 or grad.f90

infiles/ is a folder containing files used for the communication between distmesh and our Fortran programs.
