#ifndef TEST_H
#define TEST_H
//////////tools//////////

#include "krylov.h"

struct arnoldiAns {
	int k;
	double V;
	int V_c;
	double h;
	int h_c;
	double eps;
	double v;
	int v_c;

};

void tryNorm(double *vec,int n,int p,int ans);
void normTest();

void tryArrayToArray(double *A, int start_col_A, int start_row_A, int n_A, int m_A,
					 double *B, int start_col_B, int start_row_B, int n_B, int m_B,
					 int n_cols, int n_rows, double *ans );
void arrayToArrayTest();

void toolsTest();

//////////krylov//////////
void tryArnoldi(OrthogonalSet set, double *mat,int n,double e,int k,arnoldiAns ans);
void arnoldiTest();
void krylovTest();

#endif
