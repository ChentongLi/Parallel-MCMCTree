#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>
#include "header/MCMC.h"

int main(int argc, char *argv[], char *envp[]){
    
    MPI_Init(&argc, &argv);
    MCMC();
    MPI_Finalize();
    
    return 0;

}
