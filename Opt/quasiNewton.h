#pragma once
#include <vector>
#include <functional>
#include "vect_alg.h"
namespace optimize
{
	std::pair<bool, double> backtracking(std::function<double(std::vector<double>& x)> f, 
		std::function<std::vector<double>(std::vector<double>& x)> grad_f,
		std::vector<double>& xk, std::vector<double> &sk, const std::vector<double>& pk, const std::vector<double>& gradf_xk, double f_xk,
		const int maxIt = 100, const double shrink_step = 0.5, const double c1 = 0.0001, const double c2 = 0.9);

	// More-Thuente line search
	std::pair<bool, double> MT_linesearch(std::function<double(std::vector<double>& x)> f,
		std::function<std::vector<double>(std::vector<double>& x)> grad_f,
		std::vector<double>& xk, std::vector<double>& sk, const std::vector<double>& pk, const std::vector<double>& gradf_xk, double f_xk,
		const int maxIt = 1000, const double delta = 1.1, const double c1 = 0.0001, const double c2 = 0.9,
		const double alpha_max = 2.0);

	double MT_step_selection(
		const double& a_l, const double& a_u, const double& a_k,
		const double& f_l, const double& f_u, const double& f_k,
		const double& g_l, const double& g_u, const double& g_k, const double delta = 0.66
	);

	// Mininum of a quadratic function that interpolates fa, ga, and fb
	double MT_quadratic_interp(const double& a, const double& b, const double& fa, const double& ga, const double& fb);

	double MT_secant(const double& a, const double& b, const double& ga, const double& gb);

	// Mininum of a cubic function that interpolates fa, ga, fb and gb
	double MT_cubic_interp(const double& a, const double& b, const double& fa, const double& fb, const double& ga, const double& gb);

	void update_Hessian(matrix& H, const std::vector<double>& yk, const std::vector<double> sk);

	bool BFGS(std::function<double(std::vector<double>& x)> f, std::function<std::vector<double>(std::vector<double>& x)> grad_f,
		std::vector<double>& xk, const int maxIt = 1000, const double tol = 1.e-8);
}
