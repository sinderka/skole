///////////////////////////////////
////Matrix-manipulation tools//////
///////////////////////////////////
///////////////////////////////////

#include "../../../lapack/CBLAS/include/cblas.h"    //for cblas
#include "../../../lapack/LAPACKE/include/lapacke.h"//for lapack

#include <iostream> // for cout and cin
#include <cstdlib>  // for rand() and srand()
#include <ctime>    // for time()
#include <limits>   // for RAND_MAX

#include "tools.h"

// allocate

double* Tools::initArray(long length){
// creates an array of length n with only zeros
//INPUT: length , length of array
//OUTPUT: vec , array with length doubles

    double *vec = new double[length] ();

    return vec;
}

double* randomArray(long length, int seed){
//creates a vector of length n with radom elements
//INPUT: length , length of desired array
// 		 seed , a number given to make sure calls back to back does not return identical random numbers
//OUTPUT: rand_vec , a random vector of length length

    srand((u_int)(time(0)) + seed);

    double *rand_vec = Tools::initArray(length);

    for (int ii = 0; ii < length ; ii++) {

        // random number between -1 and 1
        rand_vec[ii] = ((double)rand()/RAND_MAX)*2-1;
    }
    return rand_vec;
}

void Tools::print(double* vec, long n_row, long n_col){
//Prints array with n_row rows and n_col cols
//INPUT: vec , an array of length n_row*n_col
// 		 n_row , number of rows to print vec in
//		 n_col , number of cols to print vec in

    for (long ii = 0; ii < n_row; ii++) {

		for (long jj = 0; jj < n_col; jj++){

            printf("%.10f\t", vec[ii + jj*n_row] );
		}
		std::cout << "\n";
    }
    std::cout << "\n";
}
















