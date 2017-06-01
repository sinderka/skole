#ifndef KRYLOV_H
#define KRYLOV_H


/*
class Krylov {

    private:

        double *m_A;        //Column major n*n matrix

        double *m_b;        //n vector

        long m_n;           //Height of m_b, and size of m_n

    public:
    
        Krylov(double *A, double *v, long n) : m_A(A), m_b(v), m_n(n) {}

        void print();

        void arnoldi(double *A, );

        double* project(int max_restarts);

//        void addDiag(double *A,int n,double value);

//        void integrate(double *F, int t_n, double eps);

        double * getM_A() {return m_A;}
        double * getM_b() {return m_b;}

//        double * getM_V() {return m_V;}
//        double * getM_H() {return m_H;}

//        double * getM_v() {return m_v;}

//        double getM_eps() {return m_eps;}

        long getM_n() {return m_n;}
//        long getM_k() {return m_k;}
};


*/

double arnoldi(double *A, double *b, long dim, long &dim_n, double tol,double *H, double *V);
double * project(double *A, double *b,  long dim, long dim_f, long max_restarts, 
                 long depth, long max_depth, double tol,long fs(long,long),int sol_met(double*,double*,long) );
long reduceSize(long dim, long dim_f);
int solve(double *A, double *b, long n);


#endif
