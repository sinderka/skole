#ifndef TOOLS_H
#define TOOLS_H


// allocate
double* initArray(long n);

// norm
double norm(double * vec, long n, int p);

// assign
void arrayToArray(double *A, long start_col_A, long start_row_A, long height_A,
                  double *B, long start_col_B, long start_row_B, long height_B,
                  long n_cols, long n_rows);

void integrateArnoldi(double *A,int n,double *F,int t_n,double t_s,double eps);


#endif
