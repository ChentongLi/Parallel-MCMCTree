module load openmpi
mpicc main.c -lm -lgsl -lgslcblas -o out
