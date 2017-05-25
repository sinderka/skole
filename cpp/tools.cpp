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




void integrateArnoldi(double *A,int n,double *F,int t_n,double t_s,double eps) {
// Solves the equation du/dt = A*u + eps*e_1*F[-1][:] with u(0) = 0 
//INPUT: A, an m x m array
//       m, integer, and height of A and F
//       F, an m*t_n array
//       t_n number of steps in time
//       t_s length of time-steps
//OUTPUT: U

    int inc = 1;
    double *mat = initArray(n*n);

    //mat = inv(eye(n) - t_s * mat)
    int *IPIV = new int[n] ();
    double *WORK;
    int INFO;
    cblas_daxpy(n*n,t_s,A,inc,mat,inc);
    subtractDiag(mat,n,1.0);

//lapack_int LAPACKE_dgetrf( int matrix_layout, lapack_int m, lapack_int n,
//                           double* a, lapack_int lda, lapack_int* ipiv );

//    LAPACKE_dgetrf(LAPACK_COL_MAJOR,n,n,mat,n,IPIV);
    //LAPACKE_dgetri(LAPACK_COL_MAJOR,n,mat,n,IPIV);

    //A = eps*A
    cblas_dscal(n*n,eps,A,inc);

    // Pointers to F[:][i]
    double *F_pntr_old, *F_pntr_current, F_pntr_next;

    //remember old values
    double F_old = eps*F[n-1];
    double F_current = eps*F[n*2-1];
    double F_next = eps*F[n*3-1];

    F_pntr_current = &F[n];

    cblas_dscal(n,F_current,mat,n); // Er dette rett, eller skal det være noe med eps

    for (int ii = 2; ii < t_n-1 ; ii += 2 ) {

        //Saving the values needed from F
        double F_old = F_current;                   //F[:][ii-1]
        double F_current = F_next;                  //F[:][ii]
        double F_next = eps*F[n*(ii+1) + n-1];      //eps*F[:][ii+1]

        //Making pointers to F
        double *F_pntr_old = & F[n*(ii-1) + n-1];   //F[-1][ii-1]
        double *F_pntr_current = & F[n*ii + n-1];   //F[-1][ii]
        double *F_pntr_next = & F[n*(ii+1) + n-1];  //F[-1][ii+1]


        // U[:][ii] = U[:][ii-1] + A*U[:][ii-1] + F_old*e_1; 

        // F[:][ii] = A*U[:][ii-1]
        cblas_dgemv(CblasColMajor,CblasNoTrans,n,n,1.0,A,n,F_pntr_old,inc,0.0,F_pntr_current, inc);

        // F[1][ii] += F_old
        F[n*ii + n - 1] += F_old;

        //F[:][ii] += F[:][ii-1]
        cblas_daxpy(n,1.0,F_pntr_old,inc,F_pntr_current,inc);
        

        // U[:][ii+1] = mat*( U[:][ii] + F_next*e_1 )
        F[n*ii + n - 1] += F_next;

        cblas_dgemv(CblasColMajor,CblasNoTrans,n,n,1.0,mat,n,F_pntr_current,inc,0.0,F_pntr_next,inc);
        
        F[n*ii + n - 1] -= F_next;
    }
}











