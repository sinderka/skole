#include <cmath>
#include <iostream>

#include "krylov.h"
#include "io.h"


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


void tryArnoldi(double *mat,double *vec,int n,double e,int k,arnoldiAns ans) {

    OrthogonalSet set;

	bool result = false; 
	
	try {
		set = arnoldi(mat,vec,n,e,k);

        result =  (set.eps - ans.eps < 0.0001);
        result *= (set.H[ans.h_c] - ans.h < 0.0001);
		result *= (set.V[ans.V_c] - ans.V < 0.0001);
		result *= (set.k == ans.k);
		result *= (set.v[ans.v_c] - ans.v < 0.0001);

	} catch (int e) {
		std::cout << "Error while running arrayToArray in test\n";
		result = false;
	}

	if ( !result ) {
		std::cout << "Error in arnoldi, got\n";
        std::cout << "eps: \t" << set.eps << ", expected " << ans.eps << "\n";
        std::cout << "H[end]: " << set.H[ans.h_c] << ", expected " << ans.h << "\n";
        std::cout << "V[end]: " << set.V[ans.V_c] << ", expected " << ans.V << "\n";
        std::cout << "k: \t" << set.k << ", expected " << ans.k << "\n";
        std::cout << "v: \t" << set.v[ans.v_c] << ", expected " << ans.v << "\n";
//		std::cout << "\texpected\n";
//        std::cout << "eps: " << ans << "\n";
	}
}


void arnoldiTest() {


    int n = 5;
    double *vec1 = new double[n] {1,2,3,4,5};
    double *mat1 = new double[n*n] {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25};

	arnoldiAns ans1;
	ans1.k = 2;
	ans1.V = -0.3814; ans1.V_c = 9;
	ans1.h = -4.5455; ans1.h_c = 3;
	ans1.v = 0.0; ans1.v_c = 4;
	ans1.eps = 0;
    tryArnoldi(mat1,vec1,n,0.01,3,ans1);


    n = 6;
    double *vec2 = new double[n] {1,1,1,1,1,1};
    double *mat2 = new double[n*n] ();
    mat2[0] = 1; mat2[7] = 2; mat2[14] = 3; mat2[21] = 4; mat2[28] = 5; mat2[35] = 6; 

	arnoldiAns ans2;
	ans2.k = 3;
	ans2.V = 0.545545; ans2.V_c = 17;
	ans2.h = 3.5; ans2.h_c = 8;
	ans2.v = 0.372678; ans2.v_c = 5;
	ans2.eps = 1.31747;
    tryArnoldi(mat2,vec2,n,0.01,3,ans2);

}

void krylovTest() {

    arnoldiTest();

}






