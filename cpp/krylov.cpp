#include "cblas/include/cblas.h" // for BLAS

#include "krylov.h"
#include "io.h"
#include "random_Array.h"
#include "tools.h"

#include <iostream>

#include <stdio.h>
#include <string.h>

//http://www.netlib.org/blas/
// definition of BLAS-functions:
// DGEMV    (matrix-vector multiplication) args: T, m,n,a,A,lda,x,in,b,y,inc; y = a * A * x + b * y
// DDOT     (vector-vector multiplication) args: n,x,inc,y,inc;               a = x' * y
// DAXPY    (scale and add vector)         args: n,a,x,inc,y,inc;             y = a * x + y
// DSCAL    (scale vector)                 args: n,a,x,inc;                   x = a * x
// DGEMM    (matrix-matrix multiplication) args: T_A,T_B, n_A, m_A,m_B,a,A,lda,B,ldb,b,C,ldc

OrthogonalSet arnoldi(double *A, double *b, long n, double e, long k ) {
//INPUT: A , an array acting as a n*n matrix
//		 b , an vector of length n
// 	 	 n , number of rows in A and b
// 		 e , a convergence qrterion
// 		 k , a maximum number of krylov iterations
//OUTPUT: an OrthogonalSet of the arnoldi iterations


    OrthogonalSet set;

    k = ( n < k ) ? n : k;

    double *V = initArray(n*k);
    double *H = initArray(k*k);

    double *v_pntr;
    double eps = norm(b,n,2);
    long ii = 0;
	int inc = 1;

    //b = 1/eps
    cblas_dscal(n,1/eps, b,inc);

    for (ii = 0; ii < k; ii++) {

        //V[:][ii] = b
        arrayToArray(V,ii,0,n,b,0,0,n,1,n);
        
        //v_ii = V[:][ii]
        v_pntr = &V[ii*n];

        //b = A * v_ii
        cblas_dgemv(CblasColMajor, CblasNoTrans,  n, n, 1.0, A,n,v_pntr,inc,  0.0 ,b, inc);

        for (long jj = 0; jj <= ii; jj++) {

            // v_jj = V[:][jj]
            v_pntr = &V[jj*n];

            //H[ii][jj] = v_jj^T * b;
            H[ii*k + jj] = cblas_ddot(n,v_pntr,inc,b,inc);

            //b = b - H(ii,jj) * V[:][i]
            cblas_daxpy(n,-H[ii*k + jj],v_pntr,inc,b,inc);
        }

        eps = norm(b,n,2);

        //b = 1/eps
        cblas_dscal(n,1/eps, b,inc);


        if (eps < e) {
			
            // Burde gjøres på en annen måte
            
            set.V = initArray(n*(ii+1)); // memcpy?
            set.H = initArray((ii+1)*(ii+1)); //memcpy?
            set.v = initArray(n); // no residual

            arrayToArray(set.V,0,0,n, V,0,0,n,ii+1 ,n);
            arrayToArray(set.H,0,0,ii+1,H,0,0,k,ii+1,ii+1);



            delete[] V;
			delete[] H;
			delete[] b;
            set.eps = 0;
            set.k = ii+1; 

            return set;
            
        } else if (ii < n-1) {

            //H[ii][ii+1] = eps
            H[ii*n + ii+1] = eps;
        }
	}

    set.V = V;
    set.H = H;

    set.v = b;

    set.eps = eps;
    set.k = ii;

    return set;
}





int projMet (double *A,double *b,int n,int k,double e,int max_restarts,t_n) {

OrthogonalSet set;

double *F = initArray(n*t_n); // burde få denne arrayen fra en annen funksjon

double set.eps = 10;
double eps = norm(b,n,2);
int itr = 0;



while (e > set.eps || max_restarts < itr ) {

set = arnoldi(A,b,n,e,k); // set.b og b er samme element?

intMet(G,n,inc,H,k,hn,ht); // G er "in place" endret

hn = set.eps; // Endret kriteria for å angi hvor mange iterasjoner er nok, gidder ikke sjekke hvor stort det største elementet er.


//cblas, kan løses ved DGEMM
//next_step = V[:][1:k] * G; matrix matrix ting
// F += next_step;  // gjerne med plusser og!

//v = V[:][set.k] peker hit!

//
itr += 1;
}

    return 0;
}















