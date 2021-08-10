#include <iostream>
#include <conio.h>
#include <time.h>
#include "quasiNewton.h"

double cond(std::vector<double>& x)
{
	size_t n = x.size();
	double f_sum = 0., x_sum = 0., x_prod = 1.0;
	double b = n + 1.0;
	for (size_t i = 0; i < n; ++i)
	{
		x_sum += x[i];
		x_prod *= x[i];
	}
	return x_prod - 1.;
}

double line_eq(std::vector<double>& x, int i)
{
	size_t n = x.size();
	double f_sum = 0., x_sum = 0., x_prod = 1.0;
	double b = n + 1.0;
	for (size_t i = 0; i < n; ++i)
	{
		x_sum += x[i];
		x_prod *= x[i];
	}

	return (i < n - 1) ? (x[i] + x_sum - b) : (x_prod - 1.0);
}

double func(std::vector<double>& x)
{
	size_t n = x.size();
	double f_sum = 0., x_sum = 0., x_prod = 1.0;
	double b = n + 1.0;
	for (size_t i = 0; i < n; ++i)
	{
		x_sum += x[i];
		x_prod *= x[i];
	}

	for (size_t i = 0; i < n - 1; ++i)
		f_sum += (x[i] + x_sum - b) * (x[i] + x_sum - b);
	f_sum += (x_prod - 1.0) * (x_prod - 1.0);
	
	return f_sum;
}
std::vector<double> dfunc(std::vector<double>& x)
{
	size_t n = x.size();
	std::vector<double> der(n, 0.);
	double f_sum = 0., x_sum = 0., x_prod = 1.0;
	double b = n + 1.0;
	for (size_t i = 0; i < n; ++i)
	{
		x_sum += x[i];
		x_prod *= x[i];
	}

	for (size_t i = 0; i < n; ++i)
	{
		for (size_t j = 0; j < n - 1; ++j)
		{
			if (i == j)
				der[i] += 4.* (x[j] + x_sum - b);
			else
				der[i] += 2.* (x[j] + x_sum - b);
		}
		double prod_trunc = 1.;
		for (size_t k = 0; k < n; ++k)
			if (i != k)
				prod_trunc *= x[k];
		der[i] += 2. * (x_prod - 1.) * prod_trunc;
		
	}
	return der;
}
int main()
{
	std::vector<double> xk(10, -6.75);
	
	clock_t start = clock();
	bool success = optimize::BFGS(func, dfunc, xk, 1000, 1.e-14);
	clock_t end = clock();
	double time = (double)(end - start) / CLOCKS_PER_SEC;
	std::cout.precision(15);
	std::cout << "x* = (";
	for (auto& x : xk)
		std::cout << x << ",\n";
	std::cout << ")";
	std::cout << "Time: " << time << "\n";
	std::cout << func(xk) << "\n";
	for (size_t i = 0; i < xk.size(); ++i)
		std::cout << line_eq(xk, i) << "\n";
	_getch();
}
