#ifndef ATTACH_H
#include "Attach.h"
#endif

// Definitions of the methods in testing namespace
// R. Sheehan 4 - 10 - 2021

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

	std::cout << "\nDeterminant of the Matrix A is "<<d<<"\n";

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

void testing::fredholm_solve()
{
	// Solve the Fredholm equation of the second kind for f(t)
	// f(t) = \int_{0}^{\pi/2} akfe(t, s) f(t) ds + gfe(t)
	// exact solution for original choice of gfe and akfe is sqrt(t)
	// R. Sheehan 5 - 10 - 2021

	int N = 10; 
	double a = 0, b = PI_2, exact, err, pos = PI_4, ans; 
	std::vector<double> t(N); 
	std::vector<double> w(N); 
	std::vector<double> f(N);

	ieqn::fred2(N, a, b, t, f, w, gfe, akfe); 

	std::cout << "Numerical Solution of Fredholm Equation\n"; 
	std::cout << "i , t[i] , f[i], sqrt(t[i]) , error\n"; 
	for (int i = 0; i < N; i++) {
		exact = sqrt(t[i]); 
		err = fabs(f[i] - exact); 
		std::cout << i << " , " << t[i] << " , " << f[i] << " , " << sqrt(t[i]) << " , " << err << "\n";
	}
	
	// Perform interpolation of the solution using Nystrom method
	ans = ieqn::fredin(pos, N, a, b, t, f, w, gfe, akfe);

	std::cout << "\nNystrom Interpolation of Solution\n"; 
	std::cout << "pos , appr , exact\n"; 
	std::cout << pos << " , " << ans << " , " << sqrt(pos) << "\n\n";	

	t.clear(); w.clear(); f.clear(); 
}

void testing::volterra_solve()
{
	// Solve the pair of coupled Volterra equations
	// f_{1}(t) = \int_{0}^{t} akve0 f_{1}(s) ds + \int_{0}^{t} akve1 f_{2}(s) ds + gve0
	// f_{2}(t) = \int_{0}^{t} akve1 f_{1}(s) ds + \int_{0}^{t} akve1 f_{2}(s) ds + gve1
	// exact solution is f_{1}(t) = e^{-t}, f_{2}(t) = 2 sin(t)
	// R. Sheehan 5 - 10 - 2021

	int n = 25, m = 2; 
	double h = 0.01, t0 = 0.0; 
	std::vector<double> t(n); 
	std::vector<std::vector<double>> f; 

	f = vecut::array_2D(m, n); 

	// solve the pair of volterra equations
	ieqn::voltra(n, m, t0, h, t, f, gve, akve); 

	std::cout << "pos , f1-appr, f1-exact, f2-appr, f2-exact\n"; 
	for (int i = 0; i < n; i++) {
		std::cout << t0 << " , " << f[0][i] << " , " << exp(-t0) << " , " << f[1][i] << " , " << 2.0 * sin(t0) << "\n"; 
		t0 += h; 
	}
	std::cout << "\n"; 
}

double testing::gfe(double t)
{
	// rhs function for Fredholm equation example
	double c1 = 2.25;
	double c2 = pow(PI_2, c1)/c1; // (pi/2)^2.25 / 2.25
	return (sqrt(t) - c2 * pow(t, 0.75)); 
}

double testing::akfe(double t, double s)
{
	// kernel function for Fredholm equation example
	return pow(t * s, 0.75); 
}

double testing::gve(int k, double t)
{
	// rhs function for Volterra equation example

	return (k == 0 ? cosh(t) + t * sin(t) : 2.0 * sin(t) + t * (template_funcs::DSQR(sin(t)) + exp(t))); 
}

double testing::akve(int k, int l, double t, double s)
{
	// kernel function for Volterra equation example

	return ((k == 0) ? (l == 0 ? -exp(t - s) : -cos(t - s)) : (l == 0 ? -exp(t + s) : -t * cos(s))); 
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