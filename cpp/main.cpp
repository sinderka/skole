


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


    long nnn = 8000, mmm = 16;
    long dim_f = 2;

    std::string structure = "diagonal";
    std::string type = "nice";

    Problem prob1_obj(nnn,structure,type);

/*
    std::cout << "A = \n";
    Tools::print(prob1_obj.getM_A(),nnn,nnn);
    std::cout << "b = \n";
    Tools::print(prob1_obj.getM_b(),nnn,1);
    std::cout << "x = \n";
    Tools::print(prob1_obj.getM_x(),nnn,1);
*/
/*

    double * A = new double[nnn*nnn] {};
    A[0] = 1; A[5] = 2; A[10] = 3; A[15] = 4;
    Tools::print(A,nnn,nnn);

    double * b = new double[nnn] {1,2,3,4};
*/
    long max_restarts = 4000;
    long depth = 0;
    long max_depth = 5;

    double tol = 0.000000000001;

    long (*fs)(long,long) =  reduceSize;
    int (*sol_met)(double*,double*,long) = solve;



    double * y = project(prob1_obj.getM_A(),prob1_obj.getM_b(),prob1_obj.getM_n(),dim_f,max_restarts,depth,max_depth,tol,fs,sol_met);
    //int INFO = solve(A,b,nnn);
/*
    Tools::print(y,nnn,1);
*/  
    std::cout << "max difference(rrKPM) = "<< prob1_obj.compareSolution(y) << "\n";


    solve(prob1_obj.getM_A(),prob1_obj.getM_b(),prob1_obj.getM_n());

    std::cout << "max difference(lapacke) = "<< prob1_obj.compareSolution(prob1_obj.getM_b()) << "\n";


    //std::cout << INFO <<"\n";

    return 0;
}




