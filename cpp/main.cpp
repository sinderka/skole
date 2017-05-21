#include "cblas/include/cblas.h" // for BLAS


#include <iostream>
#include <string>

//#include "../../include/lapacke.h"
//#include "../../include/cblas.h"
/*
// For QT
#include "mainwindow.h"
#include <QApplication>
#include <QPushButton>
*/


// My project files
#include "io.h"
#include "random_Array.h"
#include "tools.h"
#include "krylov.h"

//pragma once
#include <stdio.h>

//void print();


//#include "io.h"
//#include "test.cpp"


int main(int argc, char ** argv ) {

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

    return 0;
}




