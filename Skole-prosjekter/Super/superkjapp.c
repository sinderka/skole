#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "mpi.h"

int main(int argc, char **argv){

int my_id, num_procs;

int k = atof(argv[1]);

MPI_Init(&argc,&argv);
MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

long long int i;
double limi = pow(M_PI,2)/6, sum = 0,j, summer =0;


#pragma omp parallel for schedule (static) private(i) reduction(+:sum)
	for (i = my_id; i < ldexp(1,k)+1;i =i+num_procs){
		j = (double) (i+1);
		sum = sum + 1/(j*j);
	}
	
	MPI_Reduce(&sum,&summer,1, MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	
	
	if (my_id == 0){
	printf("Difference with n = 2^%d : %e\n",k, limi-summer);
	}
MPI_Finalize();
return 0;
}

