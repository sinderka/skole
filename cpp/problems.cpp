#include <math.h> // for fabs
#include <limits.h>

#include "problems.h"
#include "tools.h"


///////////////Problem///////////////

void Problem::setProblem(long n, long t_n, double t_s, double *solution) {

    m_n = n;
    m_t_n = t_n;

    m_t_s = t_s;
    
    m_solution = solution;

}

double Problem::compareSolution(double *F) {

    double max_diff = fabs(F[0] - m_solution[0] );
    double diff;

    for (int ii = 0; ii < m_n; ii ++) {

        for (int jj = 0; jj < m_t_n; jj ++) {

            diff = fabs(F[jj*m_n + ii] - m_solution[jj*m_n + ii] );
            max_diff = ( max_diff < diff ) ? diff : max_diff ;
        }
    }
    return max_diff;
}

////////////////Test1////////////////

Test1::Test1(long n,long t_n) /*: Problem()*/ {

    double value = (double) 1/(n);

    // set parent attr
    double *solution = Tools::initArray(n*t_n); // u = v*f*x

    for (int ii = 0; ii < n; ii++) {

        for (int jj = 0; jj < t_n; jj++) {

            solution[jj*n + ii] = value*ii;
        }
    }
    double t_s = (double)1/t_n;
    setProblem( n, t_n, t_s, solution);

    // instanciate variables

    m_v = Tools::initArray(m_n);               // v = ones(n,1);
    m_A = Tools::initArray(m_n*m_n);           // A[ii][ii] = -1/m_t_s, A[ii][ii+1] = 1/m_t_s

    for (int ii = 0; ii < m_n-1; ii++) {

        m_v[ii] = 1;

        m_A[ii*m_n +ii] = - value;
        m_A[ii*m_n +ii+m_n] = value;
    }

    m_v[m_n-1] = 1;
    m_A[m_n*m_n-1] = -value;    
    
    m_f = Tools::initArray(m_t_n);             // f = ones(n,1);

    for (int ii = 0; ii < m_t_n; ii++) {
    
        m_f[ii] = 1;                
    } 

    


}
