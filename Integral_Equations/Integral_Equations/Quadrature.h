#ifndef QUADRATURE_H
#define QUADRATURE_H

// Namespace containing functions for performing numerical integration 
// R. Sheehan 4 - 10 - 2021

namespace quad {
	double qgaus(double (*f)(double), double a, double b); // Gaussian Quadrature

	void gauleg(double x1, double x2, std::vector<double>& x, std::vector<double>& w, int &n);
}

#endif

