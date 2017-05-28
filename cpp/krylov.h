#ifndef KRYLOV_H
#define KRYLOV_H


class Krylov {

    private:

        double *m_A;        //Column major n*n matrix
        double *m_b;        //Initial vector

        double *m_V;        //Column major n*k matrix
        double *m_H;        //Column major k*k matrix

        double *m_v;        //Vector with height n

        double m_eps;       //Residual
        double m_tol;       //Residual

        long m_n;           // Height of V and v
        long m_k;           // Size of H, width of V

    public:
    
        Krylov(double *A, double *v, double tol,long n, long k);

        void print();

        void arnoldi();

        double* projMet(int max_restarts, int t_n, double t_s);

        void subtractDiag(double *A,int n,double value);

        void integrate(double *F, int t_n, double t_s, double eps);

        double * getM_A() {return m_A;}
        double * getM_b() {return m_b;}

        double * getM_V() {return m_V;}
        double * getM_H() {return m_H;}

        double * getM_v() {return m_v;}

        double getM_eps() {return m_eps;}

        long getM_k() {return m_k;}
};



#endif
