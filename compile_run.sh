#!/bin/bash

gfortran -c functions.f90
gfortran -c mesh.f90
gfortran -c mgmres.f90
gfortran -c curvestuff.f90
gfortran -o testing testing.f90 mesh.o curvestuff.o functions.o l2dacquadAllfiles/*.o mgmres.o -lfftw3 -llapack
./testing
