#pragma once
#include <vector>
typedef std::vector<std::vector<double>> matrix;
namespace rlinalg
{
	// dot product in R 
	double operator*(const std::vector< double >&x, const std::vector< double > &y);
	std::vector<double> operator*(const double alpha,  const std::vector< double>& x);
	std::vector<double> operator+ (const std::vector<double>& x, const std::vector<double>& y);
	std::vector<double> operator- (const std::vector<double>& x, const std::vector<double>& y);
	
	std::vector<double> mult(const matrix& M, const std::vector<double>& x);
}
