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