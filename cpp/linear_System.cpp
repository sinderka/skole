#include "../../../lapack-3.7.0/LAPACKE/include/lapacke.h"  //http://www.netlib.org/lapack/explore-html/index.html
#include "tools.h"                                          // for Tools::print(array,n_row,n_col)

#include <stdio.h>                                          // for printf

void LinearSystem::print() {

    Tools::print((char*) ("A: "));
    Tools::print(m_A,m_dim,m_dim);
    Tools::print((char*) ("b: "));
    Tools::print(m_b,m_dim,1);
    printf("n: %li\n",m_dim);
    
}

void LinearSystem::solve() {

    int IPIV = new int[m_dim] ();
    int INFO;
    int n_col_b = 1;

    int INFO = LAPACKE_dgesv(LAPACK_COL_MAJOR,m_dim,n_col_b,m_A,m_dim,IPIV,m_b,m_dim);

    delete [] IPIV;

    return INFO
}

void LinearSystem::expm() {

    return -1;


}
