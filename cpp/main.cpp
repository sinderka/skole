


#include <iostream>
#include <string>
#include <stdio.h>

// For QT
/*
#include "mainwindow.h"
#include <QApplication>
#include <QPushButton>
*/

#include "../../../lapack-3.7.0/CBLAS/include/cblas.h"
#include "../../../lapack-3.7.0/LAPACKE/include/lapacke.h"

// My project files
#include "tools.h"
#include "krylov.h"




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

    double *A2 = new double[6*6] ();
    A2[0] = 1; A2[7] = 2; A2[14] = 3; A2[21] = 4; A2[28] = 5; A2[35] = 6; 
    double *v2 = new double[6] {1,1,1,1,1,1};

    Krylov krylov_obj(A2,v2,tol,6,k);

    krylov_obj.print();

    krylov_obj.arnoldi();

    krylov_obj.print();

    

/*
    double *F = new double[n*t_n] {1,1,1,1,1,1,1,1,1,1,
                                   2,2,2,2,2,2,2,2,2,2,
                                   3,3,3,3,3,3,3,3,3,3};
    cblas_dscal(n*n,1.0,H,n);
*/

//    b og A skal ikke endre seg i Arnoldi!



    return 0;
}




