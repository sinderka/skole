#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char **argv){

int my_id, num_procs, mine[2];
MPI_Status status;


int k = atof(argv[1]);

MPI_Init(&argc,&argv);
MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

int i;
double limi = pow(M_PI,2)/6, sum = 0,j, summer =0, temp;

if(my_id == 0){
int sendes[2], ant = ldexp(1,k);
mine[0] = 0;
mine[1] = ant/num_procs;

for(i = 1; i < num_procs;i++){
sendes[0] = ant/num_procs*i;
sendes[1] = ant/num_procs*(i+1);
MPI_Send(sendes, 2, MPI_INT,i,101+i,MPI_COMM_WORLD);
}

}else {
MPI_Recv(mine, 2, MPI_INT,0,101+my_id,MPI_COMM_WORLD, &status);

}

	for (i = mine[0]; i < mine[1];i++){
	temp = (double) 1/((i+1)*(i+1));
		sum = sum + temp;
	}
	MPI_Reduce(&sum,&summer,1, MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	
	
	if (my_id == 0){
	printf("Difference with n = 2^%d : %e\n",k, limi-summer);
	}
MPI_Finalize();
return 0;
}

