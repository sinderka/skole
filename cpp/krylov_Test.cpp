//#include <cmath>    
#include <iostream> // for cout

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

void krylovTest() {

    arnoldiTest();

}






