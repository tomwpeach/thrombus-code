#!/bin/bash
./cleanAll.sh
decomposePar
mpirun -np 16 porousScalarFoam -parallel 
reconstructPar


