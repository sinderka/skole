#include <cmath>
#include <iostream>

#include "tools.h"



void tryNorm(double *vec,int n,int p,int ans) {

    double result;

    try {
        result = Tools::norm(vec,n,p);
    } catch (int e) {
        std::cout << "Error while running norm in test\n";    
    }
    if (result != ans) {
        std::cout << "Error in norm with p = " << p << ", got " << result << ", expected " << ans << "\n";
    }

}


void normTest() {
    double *vec = new double[4]; //{1,4,2,2};
    vec[0] = 1;
    vec[1] = 4;
    vec[2] = 2;
    vec[3] = 2;

    tryNorm(vec,4,0,4);
    tryNorm(vec,4,1,9);
    tryNorm(vec,4,2,5);
    tryNorm(vec,4,100,4);
    
}


void tryArrayToArray(double *A,int start_col_A,int start_row_A, int n_A, int m_A,
					 double *B,int start_col_B,int start_row_B, int n_B, int m_B,
					 int n_cols,int n_rows,double *ans ) {

	int result(0); 
	
	try {
		Tools::arrayToArray(A,start_col_A,start_row_A,n_A,B,start_col_B, start_row_B,n_B, n_cols, n_rows );

		for (int ii = 0; ii < n_A*m_A; ii++) {
			result += (ans[ii] != A[ii]);
		}

	} catch (int e) {
		std::cout << "Error while running arrayToArray in test\n";
		result = 1;
	}


	if ( result > 0 ) {
		std::cout << "Error in arrayToArray, got\n";
		Tools::print(A,n_A,m_A);
		std::cout << "\texpected\n";
		Tools::print(ans,n_A,m_A);
	}
}


void arrayToArrayTest() {

	double *mat11 = new double[15] ();
	double *mat12 = new double[15] ();
	double *mat13 = new double[15] ();
	double *mat14 = new double[15] ();
	double *mat2 = new double[3] {2,3,4};
	double *mat3 = new double[4] {5,7,6,8};
	double *mat4 = new double[5] {10,11,12,13,14};
	
	int n = 3;
	int m = 5;

	double *ans1 = new double[15] {0,0,0,2,3,4,0,0,0,0,0,0,0,0,0};
	double *ans2 = new double[15] {0,0,0,0,5,7,0,6,8,0,0,0,0,0,0};
	double *ans3 = new double[15] {0,0,0,0,8,0,0,0,0,0,0,0,0,0,0};
	double *ans4 = new double[15] {0,10,0,0,11,0,0,12,0,0,13,0,0,14,0};

	tryArrayToArray(mat11,1,0,n,m,mat2,0,0,3,1,1,3,ans1);


	tryArrayToArray(mat12,1,1,n,m,mat3,0,0,2,2,2,2,ans2);
	tryArrayToArray(mat13,1,1,n,m,mat3,1,1,2,2,1,1,ans3);
	tryArrayToArray(mat14,0,1,n,m,mat4,0,0,1,5,5,1,ans4);

}

void toolsTest() {

    normTest();
	arrayToArrayTest();

}






