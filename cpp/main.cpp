


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
#include "problems.h"




int main(int argc, char ** argv ) {


    long nnn = 4, mmm = 4;
    long dim_f = 2;

    double * A = new double[nnn*nnn] {};
    A[0] = 1; A[5] = 2; A[10] = 3; A[15] = 4;
    Tools::print(A,nnn,nnn);

    double * b = new double[nnn] {1,2,3,4};

    long max_restarts = 20;
    long depth = 0;
    long max_depth = 2;

    double tol = 0.001;

    long (*fs)(long,long) =  reduceSize;
    int (*sol_met)(double*,double*,long) = solve;



    double * y = project(A,b,nnn,dim_f,max_restarts,depth,max_depth,tol,fs,sol_met);
    //int INFO = solve(A,b,nnn);
    Tools::print(y,nnn,1);

    //std::cout << INFO <<"\n";

    return 0;
}




