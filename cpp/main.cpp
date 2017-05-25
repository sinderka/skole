//#include "cblas/include/cblas.h" // for BLAS


#include <iostream>
#include <string>


//#include "../../include/cblas.h"
/*
// For QT
#include "mainwindow.h"
#include <QApplication>
#include <QPushButton>
*/
#include "../../../lapack/CBLAS/include/cblas.h"
#include "../../../lapack/LAPACKE/include/lapacke.h"

// My project files
#include "io.h"
#include "random_Array.h"
#include "tools.h"
#include "krylov.h"

#include <stdio.h>

//void print();


//#include "io.h"
//#include "test.cpp"


int main(int argc, char ** argv ) {


/*
    long n = 5;
    long m = 5;
    int k = 3;
    int inc = 1;
    double *vec1 = new double[n] {1,2,3,4,5};

//    double *vec2 = new double[n] {1,2,3,4,5};

    double *mat = new double[n*m] {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25};


    
    OrthogonalSet set;

    set = arnoldi(mat,vec1,n,0.001,k);  

    printArray(set.V,n,set.k);

    printArray(set.H,set.k,set.k);

    printArray(set.v,n,1);

    std::cout << set.eps << "\n";

    std::cout << set.k << "\n";


    n = 6;
    m = 6;
    k = 3;
    double *vec2 = new double[n] {1,1,1,1,1,1};

    double *mat2 = new double[n*m] ();
    mat2[0] = 1; mat2[7] = 2; mat2[14] = 3; mat2[21] = 4; mat2[28] = 5; mat2[35] = 6; 



    OrthogonalSet set1;

    set1 = arnoldi(mat2,vec2,n,0.001,k);  

    printArray(set1.V,n,set1.k);

    printArray(set1.H,set1.k,set1.k);

    printArray(set1.v,n,1);

    std::cout << set1.eps << "\n";

    std::cout << set1.k << "\n";
*/


    int n = 3;
    int t_n = 10;
    double sim_time = 2;
    double t_s = sim_time/t_n;
    double eps = -1;
    double *H = new double[n*n] {5,5,5,
                                 5,5,5,
                                 5,5,5};

    double *F = new double[n*t_n] {1,1,1,1,1,1,1,1,1,1,
                                   2,2,2,2,2,2,2,2,2,2,
                                   3,3,3,3,3,3,3,3,3,3};
    cblas_dscal(n*n,1.0,H,n);

    integrateArnoldi(H,n,F,t_n,t_s,eps);

    printArray(F,n,t_n);

   double A[5][3] = {{1,2,3},{4,5,1},{3,5,2},{4,1,4},{2,5,3}};
   double b[5][2] = {{-10,12},{14,16},{18,-3},{14,12},{16,16}};
   lapack_int info,m,lda,ldb,nrhs;
   lapack_int nn;
   m = 5;
   nn = 3;
   nrhs = 2;
   lda = 5;
   ldb = 5;

    info = LAPACKE_dgels(LAPACK_COL_MAJOR,'N',m,nn,nrhs,*A,lda,*b,ldb);

    //printArray(b,2,5);



    return 0;
}




