#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

int main(int argc, char **argv){

int k = atof(argv[1]);

	int i;
	double limi = pow(M_PI,2)/6, sum=0, *v;
	v = (double*) malloc(ldexp(1,k) * sizeof(double));

#pragma omp parallel for schedule (static) private(i) reduction(+:sum)
	for (i = 0; i < ldexp(1,k) ; i++){
		v[i] = (double) 1/((i+1)*(i+1));
		sum = sum+ v[i];
		}

	printf("Difference with n = 2^%d : %e\n",k, limi-sum);

	free(v);
	return 0;
}


