//#include <cmath>    
#include <iostream> // for cout
#include <cmath>    // for fabs

#include "krylov.h"
#include "tools.h"



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


void tryArnoldi(Krylov krylov_obj, arnoldiAns ans) {

	bool result = false; 
	
	try {
        krylov_obj.arnoldi();

        double *A = krylov_obj.getM_A();
        double *b = krylov_obj.getM_b();

        double *H = krylov_obj.getM_H();
        double *V = krylov_obj.getM_V();

        double *v = krylov_obj.getM_v();

        double eps = krylov_obj.getM_eps();
        long k = krylov_obj.getM_k();

        result =  ((eps - ans.eps) < 0.0001);
        result *= ((H[ans.h_c] - ans.h) < 0.0001);

		result *= ((V[ans.V_c] - ans.V) < 0.0001);
		result *= (k == ans.k);

		result *= ((v[ans.v_c] - ans.v) < 0.0001);

	} catch (int e) {
		std::cout << "Error while running arnoldi in test\n";
		result = false;
        return;
	}

	if ( !result ) {

		std::cout << "Error in arnoldi, got\n";
        krylov_obj.print();

        std::cout << "expected: \n";
        std::cout << "V["<< ans.V_c << "] = " << ans.V << "\n";
        std::cout << "H["<< ans.h_c << "] = " << ans.h << "\n";
        std::cout << "v["<< ans.v_c << "] = " << ans.v << "\n";
        std::cout << "k = " << ans.k << "\n";
        std::cout << "eps = " << ans.eps << "\n";

	}
}


void arnoldiTest() {


    // Creating test object
    long n1 = 5, k1 = 3;
    int t_n1 = 10;
    double sim_time1 = 2;
    double t_s1 = sim_time1/t_n1;
    double tol1 = 0.01;
    double *A1 = new double[n1*n1] {1,2,3,4,5,
                                 6,7,8,9,10,
                                 11,12,13,14,15,
                                 16,17,18,19,20,
                                 21,22,23,24,25};

    double *v1 = new double[n1]   {1,2,3,4,5};

    Krylov krylov_obj1(A1,v1,tol1,n1,k1);
    
    // Creating 
	arnoldiAns ans1;
	ans1.k = 2;
	ans1.V = -0.3814; ans1.V_c = 9;
	ans1.h = -4.5455; ans1.h_c = 3;
	ans1.v = 0.0; ans1.v_c = 4;
    ans1.eps = 0;

    tryArnoldi(krylov_obj1,ans1);

    
    long n2 = 6, k2 = 3;
    int t_n2 = 10;
    double sim_time2 = 2;
    double t_s2 = sim_time2/t_n2;
    double tol2 = 0.01;
    double *A2 = new double[n2*n2] ();
    A2[0] = 1; A2[7] = 2; A2[14] = 3; A2[21] = 4; A2[28] = 5; A2[35] = 6; 
    double *v2 = new double[n2] {1,1,1,1,1,1};

    Krylov krylov_obj2(A2,v2,tol2,n2,k2);

	arnoldiAns ans2;
	ans2.k = 3;
	ans2.V = 0.545545; ans2.V_c = 17;
	ans2.h = 3.5; ans2.h_c = 8;
	ans2.v = 0.372678; ans2.v_c = 5;
    ans2.eps = 1.31747;
    
    tryArnoldi(krylov_obj2,ans2);

}

void tryIntegrate(Krylov krylov_obj,double *F, int t_n,double *ans) {

    bool result = true;

    try {

        krylov_obj.arnoldi();

        krylov_obj.integrate(F,6,krylov_obj.getM_eps());

        for (int ii = 0; ii < t_n; ii++) {

            result *= (fabs(F[ii] - ans[ii]) < 0.001);
        }


    } catch (int e) {

		std::cout << "Error while running krylov integration in test\n";
		result = false;
        return;
    }

    if (!result) {

        std::cout << "Error in integration with Arnoldi, got\n";

        Tools::print(F,krylov_obj.getM_k(),t_n);
    
        std::cout << "expected: \n";

        Tools::print(ans,krylov_obj.getM_k(),t_n);    
    }
}


void integrateTest() {

    long n = 3;
    int k = 1;
    int t_n = 6;
    double tol = 0.001;
    double *F = new double[t_n] {0,2,3,4,5,6};
    double *A = new double[n*n] {1,4,7,2,5,8,3,6,9};
    double *v = new double[t_n] {1,1,1};

    double *ans = new double[t_n] {0.0,-0.1345540323,-0.2691080646,-0.2666443151,-0.2641805655,-0.4012434598 }; 

    Krylov krylov_obj(A,v,tol,n,k);

    tryIntegrate(krylov_obj,F,t_n,ans);

}

void tryProject(Krylov krylov_obj, double * F, int t_n, double * ans, int max_restarts) {

    bool result = true;

    double *U;

    try {

        U = krylov_obj.project(F,max_restarts,t_n);

        for (int ii = 0; ii < t_n*krylov_obj.getM_n(); ii++) {

            result *= (fabs(U[ii] - ans[ii]) < 0.001);
        }


    } catch (int e) {

		std::cout << "Error while running projection in test\n";
		result = false;
        return;

    }

    if (!result) {

        std::cout << "Error in project, got\n";

        Tools::print(U,krylov_obj.getM_n(),t_n);
    
        std::cout << "expected: \n";

        Tools::print(ans,krylov_obj.getM_n(),t_n);  
    }
}

void projectTest() {


    long n = 3;
    int k = 2, t_n = 6;
    double tol = 0.001; 
    double *v = new double[n] {1,1,1};
    double *F = new double[k*t_n] {0,0,0,2,0,3,0,4,0,5,0,6};
    double *A = new double[n*n] {1,4,7,2,5,8,3,6,9};

    double * ans = new double[n*t_n] {0.0,0.0,0.0,0.0481125224,-0.0962250449,-0.2405626122,0.0962250449,-0.1924500897,
                                        -0.4811252243,0.3127313958,-0.0481125224,-0.4089564407,0.5292377468,0.0962250449,
                                        -0.3367876570,0.5653221386,-0.2646188734,-1.0945598853 };
/*
    double * ans = new double[n*t_n] {0.0,0.0481125224,0.0962250449,0.3127313958,0.5292377468,0.5653221386,
                                    0.0,-0.0962250449,-0.1924500897,-0.0481125224,0.0962250449,-0.2646188734,
                                    0.0,-0.2405626122,-0.4811252243,-0.4089564407,-0.3367876570,-1.0945598853};
*/

    Krylov krylov_obj(A,v,tol,n,k);

    tryProject(krylov_obj, F,t_n,ans,2);

}



void krylovTest() {

    arnoldiTest();
    integrateTest();
    projectTest();

}






