#pragma once
#include <vector>
#include <functional>
#include "vect_alg.h"
namespace optimize
{
	double backtracking(std::function<double(std::vector<double>& x)> f, 
		std::function<std::vector<double>(std::vector<double>& x)> grad_f,
		std::vector<double>& xk, std::vector<double> &sk, const std::vector<double>& pk, const std::vector<double>& gradf_xk, double f_xk,
		int maxIt = 100, double shrink_step = 0.5, double c1 = 0.0001, double c2 = 0.9);

	void update_Hessian(matrix& H, const std::vector<double>& yk, const std::vector<double> sk);

	bool BFGS(std::function<double(std::vector<double>& x)> f, std::function<std::vector<double>(std::vector<double>& x)> grad_f,
		std::vector<double>& xk, int maxIt = 1000, double tol = 1.e-8);
}
