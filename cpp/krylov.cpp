#include "../../../lapack/CBLAS/include/cblas.h"            //http://www.netlib.org/blas/
#include "../../../lapack-3.7.0/LAPACKE/include/lapacke.h"  //http://www.netlib.org/lapack/explore-html/index.html


#include "krylov.h"
#include "tools.h"

#include <stdio.h>      // for printf
#include <string.h>     // for memmove
#include <algorithm>    // for std::copy

#include <iostream>



Krylov::Krylov(double *A, double *v, double tol,long n, long k) {

    m_k = ( n < k ) ? n : k;

    m_n = n;

    m_A = A; 

    m_V = Tools::initArray(m_n*m_k);
    m_H = Tools::initArray(m_k*m_k);

    m_v = Tools::initArray(n); // Dårlig løsning!
    memmove(m_v,v,m_n*sizeof(double));

    //Tools::arrayToArray(m_v,0,0,m_n,v,0,0,m_n,1,m_n);

    m_b = v;

    m_eps = tol+1;
    m_tol = tol;

}

void Krylov::print() {

    Tools::print((char*) ("A: "));
    Tools::print(m_A,m_n,m_n);
    Tools::print((char*) ("V: "));
    Tools::print(m_V,m_n,m_k);
    Tools::print((char*) ("H: "));
    Tools::print(m_H,m_k,m_k);
    Tools::print((char*) ("v: "));
    Tools::print(m_v,m_n,1);
    printf("tol: %f,\teps: %f,\tn: %li,\tk: %li\n",m_tol,m_eps,m_n,m_k);
    
}


void Krylov::arnoldi() {
//Performs m_k Arnoldi-iterations or until m_eps < m_tol


    double *v_pntr;
    m_eps = Tools::norm(m_v,m_n,2);
    long ii = 0;
    int inc = 1;

    //v *= 1/eps
    cblas_dscal(m_n,1/m_eps, m_v,inc);

    for (ii = 0; ii < m_k; ii++) {

        v_pntr = &m_V[ii*m_n];        
        //V[:][ii] = b
        // burde bruke memmove
        memmove(v_pntr,m_v,m_n*sizeof(double));
        
        //v_ii = V[:][ii]
        v_pntr = &m_V[ii*m_n];

        //b = A * v_ii
        cblas_dgemv(CblasColMajor, CblasNoTrans,  m_n, m_n, 1.0, m_A,m_n,v_pntr,inc,  0.0 ,m_v, inc);

        for (long jj = 0; jj <= ii; jj++) {

            // v_jj = V[:][jj]
            v_pntr = &m_V[jj*m_n];

            //H[ii][jj] = v_jj^T * v;
            m_H[ii*m_k + jj] = cblas_ddot(m_n,v_pntr,inc,m_v,inc);

            //b = b - H(ii,jj) * V[:][i]
            cblas_daxpy(m_n,-m_H[ii*m_k + jj],v_pntr,inc,m_v,inc);
        }

        m_eps = Tools::norm(m_v,m_n,2);

        //v = 1/eps
        cblas_dscal(m_n,1/m_eps, m_v,inc);


        if (m_eps < m_tol) {

            //H = H[1:ii+1][1:ii+1]
            double *H_pntr_new;
            double *H_pntr_old;

            for (long kk = 1; kk < ii+1; kk++ ) {

                H_pntr_new = &m_H[kk*(ii+1)];
                H_pntr_old = &m_H[m_k*kk];

                memmove(H_pntr_new, H_pntr_old, (ii+1)*sizeof(double));

            }

            memmove(m_V, m_V, m_n*(ii+1)*sizeof(double));
            m_v = Tools::initArray(m_n);

            m_k = ii+1; 

            return;

        } else if (ii < m_k-1) {
            
            m_H[ii*m_n + ii+1] = m_eps;
        } else {

            return;
        }
 
    }
}


double* Krylov::project (int max_restarts,int t_n, double t_s) {

    double *F = Tools::initArray(m_n*t_n);
    double *G = Tools::initArray(m_k*t_n);

    double eps = Tools::norm(m_v,m_n,2);
    int itr = 0;
    print();


    while (m_tol > m_eps || max_restarts > itr ) {

        arnoldi();

        integrate(G,t_n,t_s,eps);

        eps = m_eps;

        //F = m_V*G + F;
        cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,m_n,m_n,m_n,1.0,m_V,m_n,G,m_n,1.0,F,m_n);

        itr += 1;
    }

    return F;
}

void Krylov::subtractDiag(double *A,int n,double value) {

    for (int ii = 0; ii < n*n; ii += n+1) {
        A[ii] -= value;
    }

}

void Krylov::integrate(double *F,int t_n,double t_s,double eps) {
// Solves the equation du/dt = A*u + eps*e_1*F[-1][:] with u(0) = 0 
//INPUT: A, an m x m array
//       m, integer, and height of A and F
//       F, an m*t_n array
//       t_n number of steps in time
//       t_s length of time-steps
//OUTPUT: U

    int inc = 1;
    double *mat = Tools::initArray(m_n*m_n);

    //mat = inv(eye(n) - t_s * mat)
    int *IPIV = new int[m_n] ();
    double *WORK;
    int INFO;
    cblas_daxpy(m_n*m_n,t_s,m_H,inc,mat,inc);
    subtractDiag(mat,m_n,1.0);


    LAPACKE_dgetrf(LAPACK_COL_MAJOR,m_n,m_n,mat,m_n,IPIV);
    LAPACKE_dgetri(LAPACK_COL_MAJOR,m_n,mat,m_n,IPIV);

    //A = eps*A
    cblas_dscal(m_n*m_n,eps,m_H,inc);

    // Pointers to F[:][i]
    double *F_pntr_old, *F_pntr_current, F_pntr_next;

    //remember old values
    double F_old = eps*F[m_n-1];
    double F_current = eps*F[m_n*2-1];
    double F_next = eps*F[m_n*3-1];

    F_pntr_current = &F[m_n];

    cblas_dscal(m_n,F_current,mat,m_n);

    for (int ii = 2; ii < t_n-1 ; ii += 2 ) {

        //Saving the values needed from F
        double F_old = F_current;                       //F[:][ii-1]
        double F_current = F_next;                      //F[:][ii]
        double F_next = eps*F[m_n*(ii+1) + m_n-1];      //eps*F[:][ii+1]

        //Making pointers to F
        double *F_pntr_old = & F[m_n*(ii-1) + m_n-1];   //F[-1][ii-1]
        double *F_pntr_current = & F[m_n*ii + m_n-1];   //F[-1][ii]
        double *F_pntr_next = & F[m_n*(ii+1) + m_n-1];  //F[-1][ii+1]


        // U[:][ii] = U[:][ii-1] + A*U[:][ii-1] + F_old*e_1; 

        // F[:][ii] = A*U[:][ii-1]
        cblas_dgemv(CblasColMajor,CblasNoTrans,m_n,m_n,1.0,m_H,m_n,F_pntr_old,inc,0.0,F_pntr_current, inc);

        // F[1][ii] += F_old
        F[m_n*ii + m_n - 1] += F_old;

        //F[:][ii] += F[:][ii-1]
        cblas_daxpy(m_n,1.0,F_pntr_old,inc,F_pntr_current,inc);
        

        // U[:][ii+1] = mat*( U[:][ii] + F_next*e_1 )
        F[m_n*ii + m_n - 1] += F_next;

        cblas_dgemv(CblasColMajor,CblasNoTrans,m_n,m_n,1.0,mat,m_n,F_pntr_current,inc,0.0,F_pntr_next,inc);
        
        F[m_n*ii + m_n - 1] -= F_next;
    }

    Tools::print(F,m_k,t_n);
}





























