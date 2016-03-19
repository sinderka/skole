#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char **argv){

	int i, k = 3;
	double limi = pow(M_PI,2)/6, sum = 0, *v;
	v = (double*) malloc(ldexp(1,14) * sizeof(double));

	for (i = 0; i < ldexp(1,14);i++){
		v[i] =(double) 1/((i+1)*(i+1));
		sum = sum + v[i];
		if (i >= ldexp(1,k)-1){
			printf("Difference with n = 2^%d : %e\n",k, limi-sum);
			k++;
		}
	}
	free(v);
	return 0;
}

