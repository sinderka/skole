


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

/*
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
    long n2 = 6; long k2 = 3;
    double *A2 = new double[6*6] ();
    A2[0] = 1; A2[7] = 2; A2[14] = 3; A2[21] = 4; A2[28] = 5; A2[35] = 6; 
    double *v2 = new double[6] {1,1,1,1,1,1};

    Krylov krylov_obj(A2,v2,tol,n2,k2);

    //krylov_obj.print();

    //projMet (int max_restarts,int t_n, double t_s)
//    double *F = krylov_obj.project(4000,t_n,t_s);

    //krylov_obj.print();

    std::cout << "sadada\n" ;

    //Tools::print(F,6,t_s);


    std::cout << "INTE\n";

    double *F = new double[n*t_n] {1,1,1,1,1,1,1,1,1,1,
                                   2,2,2,2,2,2,2,2,2,2,
                                   3,3,3,3,3,3,3,3,3,3};

    krylov_obj.integrate(F,t_n,t_s,tol);
*/

//    Problem p;
    long n3 = 5, t_n3 = 6; 
    Test1 p1();
    Test1 p2(n3,t_n3); //

    Tools::print(p2.getM_A(),n3,n3);
    Tools::print(p2.getM_v(),n3,1);
    Tools::print(p2.getM_f(),1,n3);

    return 0;
}




