


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


    long n = 5, k = 3;
    int t_n = 10;
    double sim_time = 2;
    double t_s = sim_time/t_n;
    double tol = 0.01;

    long n2 = 6; long k2 = 3;
    double *A2 = new double[3*3] {1,4,7,2,5,8,3,6,9};
    //A2[0] = 1; A2[4] = 2; A2[8] = 3;// A2[21] = 4; A2[28] = 5; A2[35] = 6; 
    double *v2 = new double[6] {1,1,1};

    double *F = new double[n*t_n] {0,0,0,
                                   0,0,2,
                                   0,0,3,
                                   0,0,4,
                                   0,0,5,
                                   0,0,6};

    Tools::print(F,3,6);

    double *A3 = Tools::initArray(3);

    Krylov krylov_obj(A2,v2,0.5,3,3);

    krylov_obj.print();

    //double *F = krylov_obj.project(1200,t_n,t_s);
    krylov_obj.integrate(F,6,(double)1/6,5);    

    krylov_obj.print();

    Tools::print(F,3,6);

    std::cout << "sadada\n" ;

    //Tools::print(F,6,t_s);

/*
    std::cout << "INTE\n";

    double *F = new double[n*t_n] {1,1,1,1,1,1,1,1,1,1,
                                   2,2,2,2,2,2,2,2,2,2,
                                   3,3,3,3,3,3,3,3,3,3};

    krylov_obj.integrate(F,t_n,t_s,tol);


//    Problem p;
    /*
    double *A = new double[25] {1,2,3,4,5,
                                 6,7,8,9,10,
                                 11,12,13,14,15,
                                 16,17,18,19,20,
                                 21,22,23,24,25};

    double *v = new double[5]   {1,2,3,4,5};

    long n3 = 5, t_n3 = 6; 
    long k3 = 5, max_restarts = 4;
    double tol3 = 0.001;
//    Test1 p1();
    Test1 p2(n3,t_n3); //

    Tools::print(p2.getM_A(),n3,n3);
    Tools::print(p2.getM_v(),n3,1);
    Tools::print(p2.getM_f(),1,t_n3);
    std::cout << p2.getM_n() << "\n" ;

    //double *A3 = p2.getM_A();
    //double *v3 = p2.getM_v();
    

    //Krylov krylov_obj3(A,v,tol3,n3,k3);

    Krylov krylov_obj3(p2.getM_A(),p2.getM_v(),tol,p2.getM_n(),k3);

    krylov_obj3.print();

    krylov_obj3.integrate();

    krylov_obj3.print();

    //krylov_obj3.project(max_restarts, p2.getM_t_n(),p2.getM_t_s() );

    krylov_obj3.print();
    */
    return 0;
}




