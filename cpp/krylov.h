#ifndef KRYLOV_H
#define KRYLOV_H


class Krylov {

    private:

        double *m_A;        //Column major n*n matrix;

        double *m_V;        //Column major n*k matrix
        double *m_H;        //Column major k*k matrix

        double *m_v;        //Vector with height n

        double m_eps;       //Residual
        double m_tol;       //Residual

        long m_n;           // Height of V and v
        long m_k;           // Size of H, width of V

    public:
    
        Krylov(double *A, double *v, double tol,long n, long k);

        void printK();

        void arnoldi();

        double* projMet(int max_restarts, int t_n, double t_s);

        void subtractDiag(double *A,int n,double value);

        void integrate(double *F, int t_n, double t_s, double eps);

};



#endif
