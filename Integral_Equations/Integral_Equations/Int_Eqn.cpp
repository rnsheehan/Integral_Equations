#ifndef ATTACH_H
#include "Attach.h"
#endif

// Definitions of the functions needed to numerically compute the solutions to integral equations
// R. Sheehan 5 - 10 - 2021

void ieqn::fred2(int n, double a, double b, std::vector<double>& t, std::vector<double>& f, std::vector<double>& w, double (*g)(double), double (*ak)(double, double))
{
	// Solves a linear Fredholm equation of the second kind. On input, a and b are the limits of integration, and n is the number of points to use in the 
	// Gaussian quadrature. g and ak are user-supplied external functions that respectively return g(t) and lK(t, s). The routine returns
	// arrays t[0..n-1] and f[0..n-1] containing the abscissas ti of the Gaussian quadrature and the solution f at these abscissas.Also returned is the 
	// array w[0..n-1] of Gaussian weights for use with the Nystrom interpolation routine fredin.

	try {
		bool c1 = n > 0 ? true : false; 
		bool c2 = b > a ? true : false; 
		bool c3 = t.size() > 0 ? true : false;
		bool c4 = t.size() == f.size() ? true : false;
		bool c5 = t.size() == w.size() ? true : false;
		bool c10 = c1 && c2 && c3 && c4 && c5; 

		if (c10) {
			int i, j;
			double d = 0.0;

			// initialise the arrays to be used
			std::vector<int> indx(n); 

			std::vector<std::vector<double>> omk; 

			omk = vecut::array_2D(n, n); 

			quad::gauleg(a, b, t, w, n); // generate the nodes and weights to be used for the integration

			// form the array 1 - l K
			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++)
					omk[i][j] = (double)(i == j) - (*ak)(t[i], t[j]) * w[j];
				f[i] = (*g)(t[i]);
			}

			// solve the set of linear equations
			lin_alg::ludcmp(omk, n, indx, d);

			lin_alg::lubksb(omk, n, indx, f);
			
			indx.clear(); omk.clear(); 
		}
		else {
			std::string reason = "Error: void fred2()\n";
			if(!c1) reason += "No. nodes is not correct\n";
			if(!c2) reason += "Integration range is not correct\n";
			if (!c3 || !c4 || !c5) reason += "Storage arrays are not correctly sized\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double ieqn::fredin(double x, int n, double a, double b, std::vector<double>& t, std::vector<double>& f, std::vector<double>& w, double (*g)(double), double (*ak)(double, double))
{
	// Given arrays t[0..n-1] and w[0..n-1] containing the abscissas and weights of the Gaussian quadrature, and given the solution array f[0..n-1] 
	// from fred2, this function returns the value of f at x using the Nystrom interpolation formula. On input, aand b are the limits of integration,
	// and n is the number of points used in the Gaussian quadrature. g and ak are user-supplied external functions that respectively return g(t) and l K(t, s).

	try{
		bool c1 = n > 0 ? true : false;
		bool c2 = b > a ? true : false;
		bool c3 = t.size() > 0 ? true : false;
		bool c4 = t.size() == f.size() ? true : false;
		bool c5 = t.size() == w.size() ? true : false;
		bool c10 = c1 && c2 && c3 && c4 && c5;

		if (c10) {
			int i;
			double sum = 0.0;

			for (i = 0; i < n; i++) sum += (*ak)(x, t[i]) * w[i] * f[i];

			return (*g)(x) + sum;
		}
		else {
			return 0.0; 
			std::string reason = "Error: double fredin()\n";
			if (!c1) reason += "No. nodes is not correct\n";
			if (!c2) reason += "Integration range is not correct\n";
			if (!c3 || !c4 || !c5) reason += "Storage arrays are not correctly sized\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void ieqn::voltra(int n, int m, double t0, double h, std::vector<double>& t, std::vector< std::vector<double> >& f, double (*g)(int, double), double (*ak)(int, int, double, double))
{
	// Solves a set of m linear Volterra equations of the second kind using the extended trapezoidal rule. On input, t0 is the starting point of the 
	// integration and n - 1 is the number of steps of size h to be taken. g(k, t) is a user-supplied external function that returns gk(t),
	// while ak(k, l, t, s)	is another user-supplied external function that returns the (k, l) element of the matrix K(t, s).
	// The solution is returned in f[0..m-1][0..n-1], with the corresponding abscissas in t[0..n-1].

	try {
		bool c1 = n > 0 ? true : false;
		bool c2 = m > 0 ? true : false;
		bool c3 = t.size() > 0 ? true : false;
		bool c4 = f.size() == m ? true : false;
		bool c5 = f[0].size() == n ? true : false;
		bool c10 = c1 && c2 && c3 && c4 && c5;

		if (c10) {
			int i, j, k, l;
			double d = 0.0, sum;

			// initialise arrays to be used
			std::vector<int> indx(m); 
			std::vector<double> b(m); 
			std::vector<std::vector<double>> a; 
			
			a = vecut::array_2D(m, m); 
			
			// accumulate rhs of linear equations in sum
			// lhs goes in matrix a
			t[0] = t0;
			for (k = 0; k < m; k++) f[k][0] = (*g)(k, t[0]);
			for (i = 1; i < n; i++) {
				t[i] = t[i - 1] + h;
				for (k = 0; k < m; k++) {
					sum = (*g)(k, t[i]);
					for (l = 0; l < m; l++) {
						sum += 0.5 * h * (*ak)(k, l, t[i], t[0]) * f[l][0];
						for (j = 1; j < i; j++)
							sum += h * (*ak)(k, l, t[i], t[j]) * f[l][j];
						a[k][l] = (k == l) - 0.5 * h * (*ak)(k, l, t[i], t[i]);
						/*if(k==l)
							a[k][l] = 1.0 - 0.5 * h * (*ak)(k, l, t[i], t[i]);
						else
							a[k][l] = -0.5 * h * (*ak)(k, l, t[i], t[i]);*/
					}
					b[k] = sum;
				}
				
				lin_alg::ludcmp(a, m, indx, d);

				lin_alg::lubksb(a, m, indx, b);
				
				for (k = 0; k < m; k++) f[k][i] = b[k];
			}

			b.clear(); indx.clear(); a.clear(); 
		}
		else {
			std::string reason = "Error: void voltra()\n";
			if (!c1) reason += "No. nodes is not correct\n";
			if (!c2) reason += "No. equations is not correct\n";
			if (!c3 || !c4 || !c5) reason += "Storage arrays are not correctly sized\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}