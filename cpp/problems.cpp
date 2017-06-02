#include "../../../lapack-3.7.0/CBLAS/include/cblas.h"

#include <math.h> // for fabs

#include <cstdlib>  // for rand() and srand()
#include <ctime>    // for time()
#include <limits>   // for RAND_MAX

#include "problems.h"
#include "tools.h"


double randomNumber(long seed) {

    srand((u_int)(time(0)) + seed);

    return ((double)rand()/RAND_MAX)*2-1;

}

double niceNumber(long ii) {

    static double counter = 0;

    counter += 1;
    
    return counter;
}
    


///////////////Problem///////////////
Problem::Problem(long n,std::string structure="diagonal", std::string type = "nice") {

    m_A = Tools::initArray(n*n);
    m_x = Tools::initArray(n);
    m_b = Tools::initArray(n);
    m_n = n;

    double (*f_pntr)(long);

    if ( type.compare("nice") ) {

        f_pntr = randomNumber;

    } else if (type.compare("random")) {

        f_pntr = niceNumber;

    }

    if ( !structure.compare("diagonal") ) {

        for (long ii = 0; ii < m_n; ii++ ) {

            m_A[ii + m_n*ii] = f_pntr(ii);
            m_x[ii] = f_pntr(ii);

        }

    } else if (!structure.compare("dense")) {
    
        for( long ii = 0; ii < m_n*m_n; ii++) {

            m_A[ii] = f_pntr(ii); 

        }

        for( long ii = 0; ii < m_n; ii++) {

            m_x[ii] = f_pntr(ii); 

        }

    } 

    int inc = 1;
    cblas_dgemv(CblasColMajor,CblasNoTrans,m_n,m_n,1.0,m_A,m_n,m_x,inc,0.0,m_b, inc);


}


double Problem::compareSolution(double *F) {

    double max_diff = fabs(F[0] - m_x[0] );
    double diff;

    for (int ii = 0; ii < m_n; ii ++) {

            diff = fabs(F[ii] - m_x[ii] );
            max_diff = ( max_diff < diff ) ? diff : max_diff ;
    }
    return max_diff;
}

