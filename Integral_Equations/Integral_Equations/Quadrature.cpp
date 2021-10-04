#ifndef ATTACH_H
#include "Attach.h"
#endif 

// Implementation of the routines in the quadrature namespace
// R. Sheehan 4 - 10 - 2021

double quad::qgaus(double (*f)(double), double a, double b)
{
	// The algorithm qgaus is taken from Numerical Recipes in C, section 4.5
	// This function will implement Gauss-Legendre Quadrature when the function is known
	// It will return the value of the integral of the function f by 10/20 point Gauss-Legendre Integration
	// 20-point Gauss-Legendre quadrature is accurate to 14 decimal places
	// R. Sheehan 30 - 5 - 2007

	try {
		if (b > a) {

			int j;
			double xr, xm, dx, s;

			//10 point Gauss-Legendre Quadrature
			int npts = 6;
			static double x[] = { 0.0, 0.1488743390, 0.4333953941, 0.6794095682, 0.8650633667, 0.9739065285 }; // Nodes taken from CRC tables 
			static double w[] = { 0.0, 0.2955242247, 0.2692667193, 0.2190863625, 0.1494513492, 0.0666713443 }; // Weights taken from CRC tables

			//20 point Gauss-Legendre Quadrature
			//Nodes and Weights taken from Abramowitz and Stegun 
			//int npts = 11;
			/*static double x[] = {0.0, 0.076526521133497333755, 0.227785851141645078080, 0.373706088715419560673, 0.510867001950827098004,
			0.636053680726515025453, 0.746331906460150792614, 0.839116971822218823395, 0.912234428251325905868,
			0.963971927277913791268, 0.993128599185094924786 };*/
			/*static double w[] = {0.0,0.152753387130725850698,0.149172986472603746788,0.142096109318382051329,0.131688638449176626898,
			0.118194531961518417312,0.101930119817240435037,0.083276741576704748725,0.062672048334109063570,
			0.040601429800386941331,0.017614007139152118312 };*/

			xm = 0.5 * (b + a);
			xr = 0.5 * (b - a);
			s = 0.0;
			for (j = 1; j < npts; j++) {
				dx = xr * x[j];
				s += w[j] * (f(xm + dx) + f(xm - dx));
			}

			return	s *= xr; // Scale the answer back to the range of integration
		}
		else {
			return 0.0; 
			std::string reason = "Error: double quad::qgaus(double (*f)(double), double a, double b)\n";
			reason += "Integration range is not correct\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}	
}

void quad::gauleg(double x1, double x2, std::vector<double>& x, std::vector<double>& w, int& n)
{
	// Given the lower and upper limits of integration x1 and x2, this routine returns arrays x[0..n-1]
	// and w[0..n - 1] of length n, containing the abscissasand weights of the Gauss - Legendre n - point
	// quadrature formula.
	// R. Sheehan 4 - 10 - 2021

	try {
		bool c1 = x2 > x1 ? true : false; 
		bool c2 = n > 0 ? true : false; 
		bool c3 = x.size() > 0 ? true : false; 
		bool c4 = x.size() == w.size() ? true : false; 
		bool c10 = c1 && c2 && c3 && c4; 

		if (c10) {
			int m, j, i;
			double z1, z, xm, xl, pp, p3, p2, p1;

			m = (n + 1) / 2;
			xm = 0.5 * (x2 + x1);
			xl = 0.5 * (x2 - x1);
			for (i = 0; i < m; i++) {
				z = cos(PI * (i + 0.75) / (n + 0.5));
				do {
					p1 = 1.0;
					p2 = 0.0;
					for (j = 0; j < n; j++) {
						p3 = p2;
						p2 = p1;
						p1 = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j+1);
					}
					pp = n * (z * p1 - p2) / (z * z - 1.0);
					z1 = z;
					z = z1 - p1 / pp;
				} while (fabs(z - z1) > EPS);
				x[i] = xm - xl * z;
				x[n - 1 - i] = xm + xl * z;
				w[i] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
				w[n - 1 - i] = w[i];
			}
		}
		else {
			std::string reason = "Error: void gauleg(double x1, double x2, std::vector<double>& x, std::vector<double>& w, int& n)\n";
			if(!c1) reason += "Integration range is not correct\n";
			if(!c2) reason += "Integration range is not correct\n";
			if(!c3 || !c4) reason += "Integration range is not correct\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}