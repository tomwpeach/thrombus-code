#!/bin/bash
./cleanAll.sh
decomposePar
mpirun -np 4 porousScalarFoam -parallel 
reconstructPar


