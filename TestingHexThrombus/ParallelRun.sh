#!/bin/bash
decomposePar | tee run_output.txt
mpirun -np 16 porousScalarFoam -parallel | tee -a run_output.txt
reconstructPar | tee -a run_output.txt