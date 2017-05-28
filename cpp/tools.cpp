///////////////////////////////////
////Matrix-manipulation tools//////
///////////////////////////////////
///////////////////////////////////

#include "../../../lapack/CBLAS/include/cblas.h"    //for cblas
#include "../../../lapack/LAPACKE/include/lapacke.h"//for lapack

#include <cmath>    // for pow
#include <iostream> // for cout and cin
#include <cstdlib>  // for rand() and srand()
#include <ctime>    // for time()
#include <limits>   // for RAND_MAX

#include "tools.h"

// allocate

double* Tools::initArray(long n){
// creates an array of length n with only zeros
//INPUT: n , size of array
//OUTPUT: vec , array with n doubles

    double *vec = new double[n] ();

    return vec;
}

double* randomArray(long n, int seed){
//creates a vector of length n with radom elements
//INPUT: n , length of desired array
// 		 seed , a number given to make sure calls back to back does not return identical random numbers
//OUTPUT: rand_vec , a random vector of length n

    srand((u_int)(time(0)) + seed);

    double *rand_vec = Tools::initArray(n);

    for (int ii = 0; ii < n ; ii++) {

        // random number between -1 and 1
        rand_vec[ii] = ((double)rand()/RAND_MAX)*2-1;
    }
    return rand_vec;
}

// norm

double Tools::norm(double * vec, long n, int p) {
//Returns the p-norm of the vector
//INPUT: vec, array of n doubles
//		 n , length of array
//		 p , the norm
//OUTPUT: sum , returns the p norm of vec, max norm assumed when p > 99

    double value = 0;

    if (p > 99) { 
    //returns maximum number in vec
        
        for (long ii = 0; ii < n; ii++) {

            value = (value < fabs(vec[ii])) ? fabs(vec[ii]) : value;
        }

		return value;

    } else if (p == 0) { 
    // returns number of elements in vec

        return n;

    } else { 
    // calculates p norm

        for (long ii = 0; ii < n; ii++) {

            value += pow(fabs(vec[ii]),p);
        }

		return pow(value,1.0/p);
    }
    // Something very wrong happend
    return -1;
}

// assign

void Tools::arrayToArray(double *A, long start_col_A, long start_row_A, long height_A,
                  double *B, long start_col_B, long start_row_B, long height_B,
                  long n_cols, long n_rows) {
// copies elements from B[start_row+:n_rows][start_col+:n_col] to A[start_row+:n_rows][start_col+:n_col] 
// matrices are collumnwise
//INPUT: A , A matrix with height height_A
//		 B , A matrix with height height_B
//		 n_cols , the number of columns to be copied from B to A
// 		 n_rows , the number of rows to be coiped from B to A

if ( start_row_A + n_rows > height_A || start_row_B + n_rows > height_B ) {
//Dimensions cannot be correct
    std::cout << ("Dimension-error in tools.cpp: arrayToArray");
}

    for (long ii = 0; ii < n_cols; ii++) {

        for (long jj = 0; jj < n_rows; jj ++) {

            A[(start_col_A + ii)*height_A + jj + start_row_A] = B[(start_col_B + ii)*height_B + jj + start_row_B];

        }
    }


/*
    int x, y,i,j;

    for (long ii = 0; ii < n_cols*n_rows; ii++) {
        j = (ii/n_rows - (ii%n_cols));
        i = (ii % n_cols);
        x = (start_col_A + i )*height_A + j + start_row_A; 
        y = (start_col_B + i )*height_B + j + start_row_B;

    A[x] = B[y];

    }
*/
}

// print

void Tools::print(auto str) {
//Prints str to consoll
//INPUT: str , any type

    std::cout << str << "\n";
}



void Tools::print(std::string str) {
//Prints str to consoll
//INPUT: str , a string

    std::cout << str << "\n";
}



void Tools::print(char* str) {
//Prints str to consoll
//INPUT: str , a string

    std::cout << str << "\n";
}


void Tools::print(double* vec, long n, long m){
//Prints array with n rows and m cols
//INPUT: vec , an array of length n*m
// 		 n , number of rows to print vec in
//		 m , number of cols to print vec in

    std::string str;

    for (long ii = 0; ii < n; ii++) {

		for (long jj = 0; jj < m; jj++){

            printf("%.6f\t", vec[ii + jj*n] );
            std::cout << str;

		}

		std::cout << "\n";
    }
    std::cout << std::endl;
}
















