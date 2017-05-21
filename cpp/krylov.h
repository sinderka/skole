#ifndef KPM_H
#define KPM_H

struct OrthogonalSet {

    double *V;
    double *H;

    double *v;

    long k;
    double eps; 
    

};


OrthogonalSet arnoldi(double *A, double *b, long n, double e, long k );

#endif
