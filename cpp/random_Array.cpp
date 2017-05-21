#include <cstdlib> // for rand() and srand()
#include <ctime> // for time()
#include <limits> // for RAND_MAX

#include "tools.h"

// Skal returnere en double vektor av størrelse n
double* randomArray(long n, int seed){
//INPUT: n , length of desired array
// 		 seed , to ensure calls back to back does not get same random numbers
//OUTPUT: rand_vec , a random vector of length n

    srand((u_int)(time(0)) + seed);

    double *rand_vec = initArray(n);

    for (int ii = 0; ii < n ; ii++) {

        // random number between -1 and 1
        rand_vec[ii] = ((double)rand()/RAND_MAX)*2-1;
    }
    return rand_vec;
}

