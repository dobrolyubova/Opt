#pragma once
#include <iostream>
#include <iomanip>
#include "quasiNewton.h"
#include <math.h>

using namespace rlinalg;

double optimize::MT_step_selection(
	const double& a_l, const double& a_u, const double& a_k,
	const double& f_l, const double& f_u, const double& f_k,
	const double& g_l, const double& g_u, const double& g_k, const double delta
)
{
	if (fabs(a_l - a_u) < std::numeric_limits<double>::epsilon())
		return a_l;

	// ac is a minimizer of the cubic interpolation of f_l, f_k, g_l, g_k
	// aq is a minimizer of the quadratic interpolation of f_l, g_l, f_k
	double ac = MT_cubic_interp(a_l, a_k, f_l, f_k, g_l, g_k);
	double aq = MT_quadratic_interp(a_l, a_k, f_l, g_l, f_k);
	// Case 1: f_k > f_l
	if (f_k > f_l)
		return (fabs(ac - a_l) < fabs(aq - a_l)) ? 
		ac :
		(0.5*(aq + ac));

	// as is a secant step
	double as = MT_secant(a_l, a_k, g_l, g_k);
	// Case 2: f_k <= f_l, g_k * g_l < 0
	if (g_k * g_l < 0.)
		return (fabs(ac - a_k) >= fabs(as - a_k)) ? ac : as;

	// Case 3: f_k <= f_l, g_k * g_l >= 0, |g_k| < |g_l|
	if (fabs(g_k) < fabs(g_l))
	{
		double ak1 = (fabs(ac - a_k) < fabs(as - a_k)) ? ac : as;
		return (a_k > a_l) ?
			std::min(a_k + delta * (a_u - a_k), ak1) :
			std::max(a_k + delta * (a_u - a_k), ak1);
	}

	// ac: is a minimizer of the cubic interpolation of f_k, f_u, g_k, g_u
	ac = MT_cubic_interp(a_k, a_u, f_k, f_u, g_k, g_u);
	// Case 4: f_k <= f_l, g_k * g_l >= 0, |g_k| > |g_l|
	return (a_k > a_l) ?
		std::min(a_k + delta * (a_u - a_k), ac) :
		std::max(a_k + delta * (a_u - a_k), ac);
}


std::pair<bool, double> optimize::MT_linesearch(std::function<double(std::vector<double>& x)> f, 
	std::function<std::vector<double>(std::vector<double>& x)> grad_f, 
	std::vector<double>& xk, std::vector<double>& sk, const std::vector<double>& pk, 
	const std::vector<double>& gradf_xk, double f_xk, const int maxIt, const double delta,
	const double c1, const double c2, const double alpha_max)
{
	std::cout << "MT search\n";

	double alpha_min = 0.,
		alpha_l = 0., alpha_u = std::numeric_limits<double>::infinity();

	std::vector<double> xk1(xk.size(), 0.), x_tmp(xk.size(), 0.);

	double alpha = 1.;
	sk = alpha * pk;

	xk1 = xk + sk;
	double f_xk1 = f(xk1);
	double desc_k1 = grad_f(xk1) * pk;
	double desc_k = gradf_xk * pk;
	if (desc_k > 0)
	{
		std::cout << "Direction is increasing objective function value!\n";
		return std::pair<bool, double>(false, alpha);
	}

	x_tmp = xk + alpha_l * pk;
	double f_l = f(x_tmp);
	double desc_l = grad_f(x_tmp) * pk;

	x_tmp = xk + alpha_u * pk;
	double f_u = f(x_tmp);
	double desc_u = grad_f(x_tmp) * pk;

	// strong Wolfe conditions
	if (f_xk1 <= f_xk + c1 * desc_k && fabs(desc_k1) <= c2 * fabs(desc_k))
	{
		xk = xk1;
		return std::pair<bool, double>(true, alpha);
	}
	int iter = 0;
	for (iter; iter < maxIt; ++iter)
	{
		// auxillary function, More&Thuente
		double	psi_k1 = f_xk1 - f_xk - c1 * alpha * desc_k,
				dpsi_k1 = desc_k1 - c1 * desc_k;
		double	psi_l = f_l - f_xk - c1 * alpha * desc_k,
				dpsi_l = desc_l - c1 * alpha * desc_k;
		double	psi_u = f_u - f_xk - c1 * alpha * desc_k,
				dpsi_u = desc_u - c1 * alpha * desc_k;

		double alpha_k1;
		if (f_xk1 > f_l)
		{
			// Case U1: f(alpha_k) > f(alpha_l)
			alpha_k1 = MT_step_selection(alpha_l, alpha_u, alpha,
				f_l, f_u, f_xk1,
				psi_l, psi_u, psi_k1);
			alpha_u = alpha;
			f_u = f_xk1;
			psi_u = psi_k1;
		}
		else if (psi_k1 * (alpha_l - alpha) > 0) 
		{
			// Case U2: f(alpha_k) <= f(alpha_l), dpsi(alpha_k) * (alpha_l - alpha_k) > 0
			alpha_k1 = std::min(alpha_max, alpha + delta * (alpha - alpha_l));

			alpha_l = alpha;
			f_l = f_xk1;
			psi_l = psi_k1;
		}
		else {
			// Case U3: f(alpha_k) <= f(alpha_l), dpsi(alpha_k) * (alpha_l - alpha_k) <= 0
			alpha_k1 = MT_step_selection(alpha_l, alpha_u, alpha,
				f_l, f_u, f_xk1,
				psi_l, psi_u, psi_k1);
			alpha_u = alpha_l;
			f_u = f_l;
			psi_u = psi_l;

			alpha_l = alpha;
			f_l = f_xk1;
			psi_l = psi_k1;
		}

		
		if (alpha_k1 > alpha_max)
		{
			std::cout << "Step size too large!\n";
		//	xk = xk1;
			return std::pair<bool, double>(false, alpha);
		}
		if (alpha_k1 < alpha_min)
		{
			std::cout << "Step size too small!\n";
		//	xk = xk1;
			return std::pair<bool, double>(false, alpha);
		}

		alpha = alpha_k1;
		x_tmp = xk + alpha_l * pk;
		f_l = f(x_tmp);

		if (fabs(alpha - alpha_max) < std::numeric_limits<double>::epsilon())
		{
			std::cout << "Max step size reached!\n";
			xk = xk1;
			return std::pair<bool, double>(true, alpha);
		}
	
		xk1 = xk + alpha * pk;
		f_xk1 = f(xk1);
		desc_k1 = grad_f(xk1) * pk;

		// strong Wolfe conditions
		if (f_xk1 <= f_xk + c1 * desc_k && fabs(desc_k1) <= c2 * fabs(desc_k))
		{
			xk = xk1;
			return std::pair<bool, double>(true, alpha);
		}
	}
}


std::pair<bool, double> optimize::backtracking(std::function<double(std::vector<double>& x)> f,
	std::function<std::vector<double>(std::vector<double>& x)> grad_f,
	std::vector<double>& xk, std::vector<double>& sk, const std::vector<double>& pk,
	const std::vector<double>& gradf_xk, double f_xk, const int maxIt, const double shrink_step,
	const double c1, const double c2)
{
	double inc_step = 2.1;

	bool down = false;
	std::vector<double> xk1(xk.size(), 0.);
	double alpha = 1.;
	double desc_k = gradf_xk * pk;

	if (desc_k > 0)
	{
		std::cout << "Direction is increasing objective function value!\n";
		return std::pair<bool, double>(false, alpha);
	}
	int i = 0;
	for (i; i < maxIt; ++i)
	{
		sk = alpha * pk;
		xk1 = xk + sk;
		double f_xk1 = f(xk1);
		double desc_k1 = grad_f(xk1) * pk;

		if (f_xk1 <= f_xk + c1 * alpha * desc_k) // First Wolfe condition
		{
			if ((desc_k1 >= c2 * desc_k)) // Second Wolfe condition
			{
				xk = xk1;
				break;
			}
			else
				alpha *= inc_step;
		}
		else
		{
			alpha *= shrink_step;
			down = true;
		}
	}
	if (i < maxIt)
		return std::pair<bool, double>(true, alpha);
	return std::pair<bool, double>(false, alpha);
}

void optimize::update_Hessian(matrix& H, const std::vector<double>& yk, const std::vector<double> sk)
{
	double skyk = sk * yk;
	double yTHy = yk * mult(H, yk);

	double scalar1 = (skyk + yTHy) / (skyk * skyk);
	double scalar2 = 1. / skyk;

	size_t sz = yk.size();
	matrix H_old = H;

	// Sherman–Morrison formula 
	for(size_t i = 0; i < sz; ++i)
		for (size_t j = 0; j < sz; ++j)
		{
			double M1_ij = 0., M2_ij = 0.;
			for (size_t k = 0; k < sz; ++k)
			{
				M1_ij += H_old[i][k] * yk[k] * sk[j];
				M2_ij += H_old[k][j] * yk[k] * sk[i];
			}
			H[i][j] += scalar1 * sk[i] * sk[j] - scalar2 * (M1_ij + M2_ij);
			
		}
}

bool optimize::BFGS(std::function<double(std::vector<double>& x)> f, 
	std::function<std::vector<double>(std::vector<double>& x)> grad_f, 
	std::vector<double>& xk, const int maxIt, const double tol)
{
	size_t sz = xk.size();
	matrix invHes(sz, std::vector<double>(sz, 0.));
	// taking identity matrix as initial inverse Hessian
	for (int i = 0; i < sz; ++i)
		invHes[i][i] = 1.;
	
	std::vector<double> xk1(sz, 0.), df_k(sz, 0.), df_k1(sz, 0.0), pk(sz, 0.), sk(sz, 0.), y(sz, 0.);
	double fx_k = f(xk), fx_k1;
	double criterion = 1.;
	size_t iter = 0;
	df_k = grad_f(xk);
	double del = fx_k;
	for (iter; iter < maxIt && criterion > tol; ++iter)
	{
		// direction
		pk = -1.0 * mult(invHes, df_k);
		xk1 = xk;
		// choose step, xk, sk are updated inside the funtion
		auto alpha = backtracking(f, grad_f, xk1, sk, pk, df_k, fx_k);
		if(!alpha.first)
			alpha = MT_linesearch(f, grad_f, xk1, sk, pk, df_k, fx_k);
		
		// preferably, use some logger object here
		if (!alpha.first)
		{
			std::cout << "No good step found! Terminating...\n";

			std::cout << "Tolerance = " << std::setprecision(15) << tol << std::endl;
			std::cout << "Iterations: " << iter << std::endl;
			std::cout << "Optimization success: false " << std::endl;
			std::cout << "Stopping criterion: " << criterion << std::endl;
			std::cout << "f(x*) = " << f(xk) << std::endl;
			return false;
		}
		
		df_k1 = grad_f(xk1);
		y = df_k1 - df_k;

		update_Hessian(invHes, y, sk);
		fx_k1 = f(xk1);
		
		criterion = std::max(fabs(fx_k1), fabs(fx_k1 - fx_k));
		std::cout << "res = " << std::setprecision(15) << criterion << std::endl;
		xk = xk1;
		fx_k = fx_k1;
		df_k = df_k1;
	}
	std::cout << "Tolerance = " << std::setprecision(15) << tol << std::endl;
	std::cout << "Iterations: " << iter << std::endl;
	std::cout << "Optimization success: " << (iter < maxIt) << std::endl;
	std::cout << "Stopping criterion: " << criterion << std::endl;
	std::cout << "f(x*) = " << f(xk) << std::endl;
	if (iter < maxIt)
		return true;
	else
		return false;
}

