#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

// Implementation of functions required for linear algebra calculations
// R. Sheehan 22 - 6 - 2018

namespace lin_alg {
	
	double inf_norm(std::vector<double>& vec);

	// Implementation of Gaussian elimination with partial pivoting 
	void gaussj(std::vector< std::vector< double > > &a, int &n, std::vector< std::vector< double > > &b, int &m);
}

#endif