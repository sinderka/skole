//#include "cblas/include/cblas.h" // for BLAS


#include <iostream>
#include <string>
#include <cblas.h>

//#include "../../include/cblas.h"
/*
// For QT
#include "mainwindow.h"
#include <QApplication>
#include <QPushButton>
*/
#include "../../../lapack-3.7.0/CBLAS/include/cblas.h"
#include "../../../lapack-3.7.0/LAPACKE/include/lapacke.h"

// My project files
#include "io.h"
#include "random_Array.h"
#include "tools.h"
#include "krylov.h"

#include <stdio.h>


int main(int argc, char ** argv ) {


    long n = 5, k = 3;
    int t_n = 10;
    double sim_time = 2;
    double t_s = sim_time/t_n;
    double tol = 0.01;
    double *A = new double[n*n] {1,2,3,4,5,
                                 6,7,8,9,10,
                                 11,12,13,14,15,
                                 16,17,18,19,20,
                                 21,22,23,24,25};

    double *v = new double[n]   {1,2,3,4,5};

    Krylov krylov_obj(A,v,tol,n,k);

    krylov_obj.printK();

    krylov_obj.arnoldi();

    krylov_obj.printK();

/*
    double *F = new double[n*t_n] {1,1,1,1,1,1,1,1,1,1,
                                   2,2,2,2,2,2,2,2,2,2,
                                   3,3,3,3,3,3,3,3,3,3};
    cblas_dscal(n*n,1.0,H,n);
*/

//    b og A skal ikke endre seg i Arnoldi!



    return 0;
}




