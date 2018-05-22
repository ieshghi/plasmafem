#!/bin/bash

gfortran -c functions.f90
gfortran -c mesh.f90
gfortran -c mgmres.f90
gfortran -c curvestuff.f90
gfortran -o testing testing.f90 mesh.o curvestuff.o functions.o l2dacquadAllfiles/*.o mgmres.o -lfftw3 -llapack

rm files/*.dat

echo Folder name: 
read fold
echo Converging on meshes in $fold

cd infiles/$fold

for d in */; do
    cp $d/* ../
    cd ../../
    ./testing #>/dev/null
    echo Mesh done
    cd infiles/$fold
done
echo Completed meshes

