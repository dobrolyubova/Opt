#pragma once
#include <iostream>
#include <iomanip>
#include "quasiNewton.h"


using namespace rlinalg;

double optimize::backtracking(std::function<double(std::vector<double>& x)> f, 
	std::function<std::vector<double>(std::vector<double>& x)> grad_f, 
	std::vector<double>& xk, std::vector<double>& sk, const std::vector<double>& pk, 
	const std::vector<double>& gradf_xk, double f_xk, int maxIt, double shrink_step, double c1, double c2)
{
	bool down = false;
	std::vector<double> xk1(xk.size(), 0.);
	double alpha = 1.;
	double desc_k = gradf_xk * pk;

	for (int i = 0; i < maxIt; ++i)
	{
		sk = alpha * pk;
		xk1 = xk + sk;
		double f_xk1 = f(xk1);
		double desc_k1 = grad_f(xk1) * pk;

		if (f_xk1 <= f_xk + c1 * alpha * desc_k) // First Wolfe condition
		{
			if ((desc_k1 >= c2 * desc_k) || down) // Second Wolfe condition
			{
				xk = xk1;
				break;
			}
			else
				alpha /= shrink_step; // shrink_step < 1

		}
		else
		{
			alpha *= shrink_step;
			down = true;
		}
	}
	return alpha;
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

bool optimize::BFGS(std::function<double(std::vector<double>& x)> f, std::function<std::vector<double>(std::vector<double>& x)> grad_f, std::vector<double>& xk, int maxIt, double tol)
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
	for (iter; iter < maxIt && criterion > tol; ++iter)
	{
		// direction
		pk = -1.0 * mult(invHes, df_k);
		xk1 = xk;
		// choose step, xk, sk are updated inside the funtion
		double alpha = backtracking(f, grad_f, xk1, sk, pk, df_k, fx_k);
		
		df_k1 = grad_f(xk1);
		y = df_k1 - df_k;

		update_Hessian(invHes, y, sk);
		fx_k1 = f(xk1);
		
		criterion = fabs(fx_k1 - fx_k);
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

