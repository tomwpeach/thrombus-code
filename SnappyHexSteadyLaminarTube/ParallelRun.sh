#!/bin/bash
decomposePar | tee run_output.txt
mpirun -np 4 simpleFoam -parallel | tee -a run_output.txt
reconstructPar | tee -a run_output.txt