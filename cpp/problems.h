#ifndef PROBLEM_H
#define PROBLEM_H

#include <string>

class Problem {

    protected: 

        double *m_A;
        double *m_b;
        
        long m_n;

        double * m_x;
        
    public:

        Problem(long n,std::string structure, std::string type);

        double compareSolution(double *x);

        double* getM_A() {return m_A;}
        double* getM_b() {return m_b;}
        long getM_n() {return m_n;}

        double * getM_x() {return m_x;}
};



#endif
