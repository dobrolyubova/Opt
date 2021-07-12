#pragma once
#include "vect_alg.h"

// Linear algenbra vector operations in R and C

namespace rlinalg
{
	// dot product
	double rlinalg::operator*( const std::vector< double >& x, const std::vector< double >& y)
	{
		size_t sz = x.size();
		double res = 0.;
//#pragma omp parallel for reduction(+:res) num_threads(THREADS) 
		for (int i = 0; i < sz; ++i)
			res += x[i] * y[i];

		return res;
	}

	std::vector<double> operator*(const double alpha,  const std::vector<double>& x)
	{
		size_t sz = x.size();
		std::vector<double> res( sz, 0.);
			
		for (int i = 0; i < sz; ++i)
			res[i] = alpha*x[i];

		return res;
	}

	std::vector<double> operator+(const std::vector<double> &x, const std::vector<double>& y)
	{
		size_t sz = x.size();
		std::vector<double> res(sz, 0.);

		for (int i = 0; i < sz; ++i)
			res[i] = x[i] + y[i];

		return res;
	}

	std::vector<double> operator-(const std::vector<double> &x, const std::vector<double>& y)
	{
		size_t sz = x.size();
		std::vector<double> res(sz, 0.);

		for (int i = 0; i < sz; ++i)
			res[i] = x[i] - y[i];

		return res;
	}

	std::vector<double> mult(const matrix& M, const std::vector<double>& x)
	{
		size_t sz = x.size();
		std::vector<double> res(sz, 0.);
		for (size_t i = 0; i < sz; ++i)
			res[i] = M[i] * x;
		return res;
	}
}
