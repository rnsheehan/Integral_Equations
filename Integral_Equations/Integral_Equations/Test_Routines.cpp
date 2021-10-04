#ifndef ATTACH_H
#include "Attach.h"
#endif

void testing::gaussj_test()
{
	// Ensure the Gauss elimination scheme is working correctly
	// R. Sheehan 22 - 6 - 2018

	int nrows = 4, ncols = 4, ncolsb = 1;
	std::vector<std::vector<double>> A;
	std::vector<std::vector<double>> b;

	A = vecut::array_2D(nrows, ncols);
	b = vecut::array_2D(nrows, ncolsb);

	// Fill the matrix A
	A[0][0] = 1; A[0][1] = -1; A[0][2] = 2; A[0][3] = -1;
	A[1][0] = 2; A[1][1] = -2; A[1][2] = 3; A[1][3] = -3;
	A[2][0] = 1; A[2][1] = 1; A[2][2] = 1; A[2][3] = 0;
	A[3][0] = 1; A[3][1] = -1; A[3][2] = 4; A[3][3] = 3;

	// Fill the vector B
	b[0][0] = -8;
	b[1][0] = -20;
	b[2][0] = -2;
	b[3][0] = 4;

	std::cout << "Matrix A is\n";
	vecut::print_mat(A); 
	
	std::cout << "\nVector b is\n";
	vecut::print_mat(b); 	

	lin_alg::gaussj(A, nrows, b, ncolsb);

	std::cout << "\nInverse of A is\n";
	vecut::print_mat(A);
	
	std::cout << "\nSolution of system of equations is\n";
	vecut::print_mat(b);	
	std::cout << "\n";

	A.clear(); b.clear();
}

void testing::ludcmp_test()
{
	// Ensure the LU Decomposition scheme is working correctly
	// R. Sheehan 4 - 10 - 2021

	int nrows = 4, ncols = 4, ncolsb = 1; 
	double d = 0.0; 
	std::vector<std::vector<double>> A;
	std::vector<std::vector<double>> Ainv;
	std::vector<double> b(nrows);
	std::vector<int> indx(nrows); 

	A = vecut::array_2D(nrows, ncols);
	Ainv = vecut::array_2D(nrows, ncols);

	// Fill the matrix A
	A[0][0] = 1; A[0][1] = -1; A[0][2] = 2; A[0][3] = -1;
	A[1][0] = 2; A[1][1] = -2; A[1][2] = 3; A[1][3] = -3;
	A[2][0] = 1; A[2][1] = 1; A[2][2] = 1; A[2][3] = 0;
	A[3][0] = 1; A[3][1] = -1; A[3][2] = 4; A[3][3] = 3;

	// Fill the vector B
	b[0] = -8;
	b[1] = -20;
	b[2] = -2;
	b[3] = 4;

	std::cout << "Matrix A is\n";
	vecut::print_mat(A);

	std::cout << "\nVector b is\n";
	vecut::print_vec(b);

	lin_alg::ludcmp(A, nrows, indx, d); 
	lin_alg::lubksb(A, nrows, indx, b); 

	std::cout << "\nLU decomposition of A is\n";
	vecut::print_mat(A);

	std::cout << "\nSolution of system of equations is\n";
	vecut::print_vec(b);
	
	lin_alg::ludet(A, d, nrows); 

	std::cout << "Determinant of the Matrix A is "<<d<<"\n";

	lin_alg::luinv(A, Ainv, indx, nrows); 

	std::cout << "\nInverse of the Matrix A is\n"; 
	vecut::print_mat(Ainv);
	std::cout << "\n";

	A.clear(); b.clear(); Ainv.clear(); indx.clear(); 
}

void testing::numerical_integration()
{
	// testing of the numerical integration routines
	// R. Sheehan 4 - 10 - 2021

	std::cout << "Integration by Gaussian Quadrature\n"; 
	std::cout << "Integral g1 on [-1, 1] = " << quad::qgaus(g1, -1, 1) << "\n"; 
	std::cout << "Integral g2 on [1, 10] = " << quad::qgaus(g2, 1, 10) << "\n"; 
	std::cout << "Integral g3 on [3, 6] = " << quad::qgaus(g3, 3, 6) << "\n\n"; 

	int n = 10;
	double x1 = 1, x2 = 10;
	std::vector<double> X(n);
	std::vector<double> W(n);

	quad::gauleg(x1, x2, X, W, n);

	std::cout << "Gauss-Legendre Quadrature Nodes and Weights\n";
	std::cout << "N = " << n << "\n";
	for (int i = 0; i < n; i++) {
		std::cout << i << " , " << X[i] << " , " << W[i] << "\n";
	}
	std::cout << "\n";

	double sum = 0.0; 
	double exact = g5(x2) - g5(x1); 
	for (int i = 0; i < n; i++) {
		sum += W[i] * g4(X[i]); 
	}
	std::cout << "Gauss-Legendre Integral of g4 on [1, 10] = " << sum << "\n"; 
	std::cout << "Exact value of Integral of g4 on [1, 10] = " << exact << "\n"; 
	std::cout << "Error in Numerical Integral of g4 on [1, 10] = " << fabs(exact - sum) << "\n\n"; 
}

double testing::g1(double x)
{
	return (exp(-(x * x)));
	// exact = sqrt(pi/4) Erf(x)
	// exact = sqrt(pi) Erf(1) ~ 1.49365
	//error ~ 10^{-30} over [-1,1]
}

double testing::g2(double x)
{
	return (x * x * log(x));
	// exact = (x^{3}/3) Log(x) - (x^{3}/9) + C
	// exact = -111 + (1000/9)Log[1000] ~ 656.528
	//error ~ 10^{-12} over [1,10]
}

double testing::g3(double x)
{
	return ((2.0 * x) / (x * x - 4.0));
	// exact = Log[x^{2}-4] + C
	// exact = Log[32/5] ~ 1.8563
	//error ~ 10^{-15} over[3, 6]
}

double testing::g4(double x)
{
	return (x*exp(-x));
	// exact = -(1 + x) e^{-x} + C
	// exact = (2 e^{9} - 11)/e^{10} ~ 0.735259
	//error ~ 10^{-15} over[1, 10]
}

double testing::g5(double x)
{
	return (-1.0 * (1 + x) * exp(-x));
}