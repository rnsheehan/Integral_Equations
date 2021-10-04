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
	/*for (int i = 0; i < nrows; i++) {
		for (int j = 0; j < ncols; j++)
			std::cout << A[i][j] << " ";
		std::cout << "\n";
	}*/
	std::cout << "\nVector b is\n";
	vecut::print_mat(b); 
	/*for (int i = 0; i < nrows; i++) {
		std::cout << b[i][0] << "\n";
	}*/

	lin_alg::gaussj(A, nrows, b, ncolsb);

	std::cout << "\nInverse of A is\n";
	vecut::print_mat(A);
	/*for (int i = 0; i < nrows; i++) {
		for (int j = 0; j < ncols; j++)
			std::cout << A[i][j] << " ";
		std::cout << "\n";
	}*/
	std::cout << "\nSolution of system of equations is\n";
	vecut::print_mat(b);
	/*for (int i = 0; i < nrows; i++) {
		std::cout << b[i][0] << "\n";
	}*/

	A.clear(); b.clear();
}