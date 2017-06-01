#ifndef PROBLEM_H
#define PROBLEM_H

class Problem {

    protected: 

        long m_n;
        long m_t_n;

        double m_t_s;

        double * m_solution;
        
    public:

        Problem() {}

        double compareSolution(double *F);
        void setProblem(long n, long t_n, double* solution);

        long getM_n() {return m_n;}
        long getM_t_n() {return m_t_n;} 

        double * getM_solution() {return m_solution;}
};

class Test1 : public Problem {

// du/dx = v*f from 0 to 1 transformed to the equation:
// A = v*f, where A[i][i] = -1/h, A[i][i+1] = 1/h, and otherwise zero, v = ones(n,1), f = ones(1,t_n)

    private:

        double *m_A;
        double *m_v;
        double *m_f;

    public:
        
        Test1() {}
        Test1(long n,long t_n);

         double * getM_A() {return m_A;}
         double * getM_v() {return m_v;}
         double * getM_f() {return m_f;}


};





#endif
