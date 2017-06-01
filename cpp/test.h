#ifndef TEST_H
#define TEST_H

#include "tools.h"
#include "krylov.h"

//tools

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


//arnoldi

void tryArnoldi(Krylov krylov_obj,arnoldiAns ans);
void arnoldiTest();

void tryIntegrate(Krylov krylov_obj,double *F, int t_n,double *ans);
void integrateTest();

void tryProject(Krylov krylov_obj, double * F, int t_n, double * ans, int max_restarts);
void projectTest();


void krylovTest();

#endif
