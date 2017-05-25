#ifndef KPM_H
#define KPM_H

struct OrthogonalSet {

    double *V;
    double *H;

    double *v;

    long k;
    double eps; 
    

};


void arnoldi(OrthogonalSet &set,double *A, long n, long k,double e );

#endif
