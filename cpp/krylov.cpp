#include "../../../lapack/CBLAS/include/cblas.h"            //http://www.netlib.org/blas/
#include "../../../lapack-3.7.0/LAPACKE/include/lapacke.h"  //http://www.netlib.org/lapack/explore-html/index.html


#include "krylov.h"
#include "tools.h"
#include "linear_System.h"

#include <stdio.h>      // for printf
#include <string.h>     // for memmove
#include <cstring>      // for memset

#include <iostream>


double arnoldi(double *A, double *b, long nnn, long &mmm, double tol,double *H, double *V) {
// Arnoldi(A,b,nnn,mmm,tol,H,V,v);
//n -> nnn
//k -> mmm
//Performs Arnoldi-iterations while until m_eps < m_tol and ii < k

    int inc = 1;

    double *v_pntr;
    double eps = cblas_dnrm2(nnn,b,inc);


    //b = 1/eps*b
    cblas_dscal(nnn,1/eps, b,inc);

    for (int ii = 0; ii < mmm; ii++) {

        v_pntr = &V[ii*nnn];
        
        //V[:][ii] = b
        memmove(v_pntr,b,nnn*sizeof(double));
        
        //v_ii = V[:][ii]
        v_pntr = &V[ii*nnn];

        //v = A * v_ii
        cblas_dgemv(CblasColMajor, CblasNoTrans,  nnn, nnn, 1.0, A,nnn,v_pntr,inc,  0.0 ,b, inc);

        for (long jj = 0; jj <= ii; jj++) {

            // v_jj = V[:][jj]
            v_pntr = &V[jj*nnn];

            //H[ii][jj] = v_jj^T * v;
            H[ii*mmm + jj] = cblas_ddot(nnn,v_pntr,inc,b,inc);

            //b = b - H(ii,jj) * V[:][i]
            cblas_daxpy(nnn,-H[ii*mmm + jj],v_pntr,inc,b,inc);
        }

        //eps = Tools::norm(m_v,m_n,2);
        eps = cblas_dnrm2(nnn,b,inc);

        //b = 1/eps*b
        cblas_dscal(nnn,1/eps, b,inc);


        if (eps < tol) {

            //H = H[1:ii+1][1:ii+1]
            double *H_pntr_new;
            double *H_pntr_old;

            for (long kk = 1; kk < ii+1; kk++ ) {

                H_pntr_new = &H[kk*(ii+1)];
                H_pntr_old = &H[mmm*kk];

                memmove(H_pntr_new, H_pntr_old, (ii+1)*sizeof(double));

            }

            memmove(V, V, nnn*(ii+1)*sizeof(double));
            std::memset(b, 0, nnn*sizeof(double));

            mmm = ii+1; 

            return eps;

        } else if (ii < mmm-1) {
            H[ii*mmm + ii+1] = eps;
        } else {
            return eps;
        }
 
    }
    return eps;

}

long reduceSize(long nnn, long nnn_f) {

    long size = (long) pow(nnn,0.5);

    if ( size < nnn_f) {

        return nnn_f;
    }
    return size;
}

int solve(double *A, double *b,long nnn) {

    int *IPIV = new int[nnn] ();
    int n_col_b = 1;

    int INFO = LAPACKE_dgesv(LAPACK_COL_MAJOR,nnn,n_col_b,A,nnn,IPIV,b,nnn);

    delete [] IPIV;

    return INFO;
}


double * project(double *A, double *b,  long nnn, long nnn_f, long max_restarts, 
                 long depth, long max_depth, double tol,long fs(long,long),int sol_met(double*,double*,long) ) {

    std::cout << "depth: " << depth << "\n";

    long mmm = fs(nnn,nnn_f);

    double eps = cblas_dnrm2(nnn,b,1);
    double new_eps;
    long itr = 0;
    int inc = 1;
    double *ptr;

    double * H = Tools::initArray(mmm*mmm);
    double * V = Tools::initArray(nnn*mmm);
    double * y = Tools::initArray(nnn);
    double * x = Tools::initArray(mmm);


    std::cout << "mmm = " << mmm << "\n";



    new_eps = arnoldi(A,b,nnn,mmm,tol,H,V);

    x[0] = eps;

    if (mmm == nnn_f || depth >= max_depth) {
        int INFO = sol_met(H,x,mmm);
    } else {
        x = project(H, x, mmm, nnn_f, max_restarts, depth+1, max_depth, tol, fs, sol_met );
    }
    std::cout << "x = \n";
    Tools::print(x,mmm,1);

    
    cblas_dgemv(CblasColMajor,CblasNoTrans,nnn,mmm,1.0,V,nnn,x,inc,1.0,y, inc);
    std::cout << "y = \n";
    Tools::print(y,nnn,1);

    std::cout <<"eps = " << eps << ", tol = " << tol << ", itr = " << itr << ", max_restarts = " << max_restarts << "\n";
    std::cout << "=============================================================================\n";
    eps = new_eps;
    do {

        new_eps = arnoldi(A,b,nnn,mmm,tol,H,V);
        

        //x = eps*x[-1]
        x[0] = eps*x[mmm-1];
        ptr = &x[1];
        memset(ptr, 0, (mmm-1)*sizeof(double));

        std::cout << "x = \n";
        Tools::print(x,mmm,1);

        std::cout <<"eps = " << eps << ", tol = " << tol << ", itr = " << itr << ", max_restarts = " << max_restarts << "\n";
        std::cout << "=============================================================================\n";

        //x = A\x
        int INFO = sol_met(H,x,mmm);

        std::cout << "(after)x = \n";
        Tools::print(x,mmm,1);

        //y += V*x
        cblas_dgemv(CblasColMajor,CblasNoTrans,nnn,mmm,-1.0,V,nnn,x,inc,1.0,y, inc);

        eps = new_eps;
        itr += 1;
    } while (eps*x[mmm-1] > tol && itr < max_restarts);
    return y;
}
/*
void Krylov::addDiag(double *A,int n,double value) {

    for (int ii = 0; ii < n*n; ii += n+1) {
        A[ii] += value;
    }

}

void Krylov::integrate(double *F,int t_n,double eps) {
// Solves the equation du/dt = A*u + eps*e_1*F[-1][:] with u(0) = 0 
//INPUT: A, an m x m array
//       m, integer, and height of A and F
//       F, an m*t_n array
//       t_n number of steps in time
//OUTPUT: U

    int inc = 1;

    //mat = inv(eye(n)/ht - m_H)
    double *mat = Tools::initArray(m_k*m_k);
    cblas_daxpy(m_k*m_k,-1.0,m_H,inc,mat,inc);
    int *IPIV = new int[m_k] ();
    addDiag(mat,m_k,(double)1/eps);
    LAPACKE_dgetrf(LAPACK_COL_MAJOR,m_k,m_k,mat,m_k,IPIV);
    LAPACKE_dgetri(LAPACK_COL_MAJOR,m_k,mat,m_k,IPIV);

    //A = eps*A
    cblas_dscal(m_k*m_k,eps,m_H,inc);

    // Pointers to F[:][i]
    double *F_pntr_old, *F_pntr_current, *F_pntr_next;
    F_pntr_current = &F[m_k];

    //remember old values
    double F_old = eps*F[m_k - 1];
    double F_current = eps*F[m_k*2 - 1];
    double F_next = eps*F[m_k*2 - 1];
    double F_far = eps*F[m_k*3 - 1];

    //F[:][1] = F[-1][1]*mat[:][0]
    std::memset(F_pntr_current, 0, m_k*sizeof(double));
    cblas_daxpy(m_k,F_current/eps,mat,inc,F_pntr_current,inc);

    for (int ii = 2; ii < t_n-1 ; ii += 2 ) {

        //Saving the values needed from F
        F_old = F_next;                             //F[:][ii-1]
        F_current = F_far;                          //F[:][ii]
        F_next = eps*F[m_k*(ii+1) + m_k - 1 ];      //F[:][ii+1]
        F_far = eps*F[m_k*(ii+2) + m_k - 1 ];       //F[:][ii+2]

        //Making pointers to F
        F_pntr_old = & F[m_k*(ii-1)];   //F[-1][ii-1]
        F_pntr_current = & F[m_k*ii];   //F[-1][ii]
        F_pntr_next = & F[m_k*(ii+1)];  //F[-1][ii+1]


        // U[:][ii] = U[:][ii-1] + A*U[:][ii-1] + F_old*e_1; 

        // F[:][ii] = A*U[:][ii-1]
        cblas_dgemv(CblasColMajor,CblasNoTrans,m_k,m_k,1.0,m_H,m_k,F_pntr_old,inc,0.0,F_pntr_current, inc);

        // F[0][ii] += F_old
        F[m_k*ii] += F_old;

        //F[:][ii] += F[:][ii-1]
        cblas_daxpy(m_k,1.0,F_pntr_old,inc,F_pntr_current,inc);

        // F[:][ii+1] = mat*( U[:][ii] + F_next*e_1 )
        F[m_k*ii] += F_next;
        cblas_dgemv(CblasColMajor,CblasNoTrans,m_k,m_k,1/eps,mat,m_k,F_pntr_current,inc,0.0,F_pntr_next,inc);
        F[m_k*ii] -= F_next;

    }

    delete[] IPIV;
    delete[] mat;
    
}








*/




















