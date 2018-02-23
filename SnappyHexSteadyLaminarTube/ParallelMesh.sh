#!/bin/bash

blockMesh > mesh_output.txt
decomposePar | tee mesh_output.txt
mpirun -np 4 snappyHexMesh -overwrite -parallel | tee -a mesh_output.txt
reconstructParMesh -constant | tee -a mesh_output.txt
