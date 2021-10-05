#ifndef INT_EQN_H
#define INT_EQN_H

// Namespace that contains the methods used to numerically compute the solutions of Fredholm and Volterra equations of the second kind
// Fredholm equations of the 1st kind are not treated as they are numerically too difficult, at least according to NRinC
// This would require a study of techniques used to solve inverse problems, which is more than I want to do at this point
// Also, homogeneous equations i.e. those with g = 0, are not treated
// As I understand it homogeneous equations reduce to a solution of a system of linear equations or a matrix eigenvalue problem
// The code for treating problems in which the Kernel is singular is also not included. I'm not interested enough to pursue it. 
// R. Sheehan 5 - 10 - 2021

namespace ieqn {
	void fred2(int n, double a, double b, std::vector<double> &t, std::vector<double>& f, std::vector<double>& w, double (*g)(double), double (*ak)(double, double));

	void voltra(int n, int m, double t0, double h, std::vector<double>& t, std::vector< std::vector<double> >& f, double (*g)(int, double), double (*ak)(int, int, double, double));

	double fredin(double x, int n, double a, double b, std::vector<double>& t, std::vector<double>& f, std::vector<double>& w, double (*g)(double), double (*ak)(double, double));
}

#endif
