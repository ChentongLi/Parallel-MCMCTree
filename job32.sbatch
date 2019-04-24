#!/bin/bash

#SBATCH --job-name=ParaTree-mpi
#SBATCH --output=mpintask32.out

#SBATCH --ntasks=32
#SBATCH --mem=10g
#SBATCH --partition=sandyb

module load openmpi

mpirun ./out
