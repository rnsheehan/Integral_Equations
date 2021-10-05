#ifndef TEST_ROUTINES_H
#define TEST_ROUTINES_H

// Declaration of the testing namespace
// Functions for testing methods go here
// R. Sheehan 4 - 10 - 2021

namespace testing {

	void gaussj_test(); 

	void ludcmp_test();

	void numerical_integration(); 

	void fredholm_solve(); 

	void volterra_solve(); 

	double gfe(double t); 

	double akfe(double t, double s); 

	double gve(int k, double t); 

	double akve(int k, int l, double t, double s); 

	double g1(double x); 
	
	double g2(double x); 

	double g3(double x); 

	double g4(double x); 

	double g5(double x); 
}

#endif


