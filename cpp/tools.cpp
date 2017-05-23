///////////////////////////////////
////Matrix-manipulation tools//////
///////////////////////////////////
///////////////////////////////////

#include "../../../lapack/CBLAS/include/cblas.h"
#include "../../../lapack/LAPACKE/include/lapacke.h"

#include <cmath> // for pow

#include "io.h"
#include <iostream>

// allocate

double* initArray(long n){
//INPUT: n , size of array
//OUTPUT: vec , array with n doubles


    double *vec = new double[n] ();

    return vec;
}

// norm

double norm(double * vec, long n, int p) {
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

void arrayToArray(double *A, long start_col_A, long start_row_A, long height_A,
                  double *B, long start_col_B, long start_row_B, long height_B,
                  long n_cols, long n_rows) {
// copies elements from B[start_row+:n_rows][start_col+:n_col] to A[start_row+:n_rows][start_col+:n_col] 
// matrices are collumnwise
//INPUT: A , 
//		 b , A matrix with hieght height_B
//		 n , number of rows in A, and length of b
// 		 m , the number of cols in A
// 	 	 col , the col to be assigned value to A from b
//OUTPUT: void

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

void subtractDiag(double *A,int n,double value) {

    for (int ii = 0; ii < n*n; ii += n+1) {
        A[ii] -= value;
    }

}



////BURDE TA INN ET HINT OM HVORDAN A SER UT!
double* explisitIntegration(double *A,int n,double *F,int t_n,double t_s) {
// Solves the equation du/dt = A*u + F
//INPUT: A, an m x m array
//       m, integer, and height of A and F
//       F, an m*t_n array
//       t_n number of steps in time
//       t_s length of time-steps
//OUTPUT: U

double *mat = initArray(n*n);

double *U = initArray(n*t_n); // vil helst klare meg uten denne!

//mat = eye(n) - t_s * A
cblas_daxpy(n,t_s,A,1,mat,1);
subtractDiag(mat,n,1.0);

//U = F
cblas_daxpy(n,t_s,F,1,U,1);

//mat = t_s * mat
//cblas_dscal(n,t_s,mat,1);

//mat = eye(m) - mat;



int INFO, IPIV;

//U[:][2] = mat\( ht * F[:][2] ) );
LAPACKE_dgelsy(n,1,A,n,IPIV,U[n*2],n,INFO) // arg, denne gir ut A som en lu-faktorisering av A i A
//LAPACKE_dgelsy(n,1,mat,n,mat,);




for (int ii = 3; ii < t_n-1 ; ii += 2 ) {

    //U[:][i] = U[:][i-1] + ht*(A*U[:][i-1] + F[:][i-1] );
    //                        daxpy(cblas_dgemv(A*U) +F[:][i-1]);

    //U[:][i+1] = mat\( U[:][i] + ht*F[:][i+1] );
    
}

return U;

}











