#ifndef KRYLOV_H
#define KRYLOV_H

class LinearSystem {

    public:

        double *m_A;        //Column major n*n matrix

        double *m_b;        //n vector

        long m_dim;           //Height of m_b, and size of m_n

        LinearSystem(double *A, double *v, long dim) : m_A(A), m_b(v), m_dim(dim) {}

        void print();

        int solve();
        int expm(); 


        double * getM_A() {return m_A;}
        double * getM_b() {return m_b;}
        long getM_dim() {return m_dim;}
};

#endif
