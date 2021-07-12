#include <iostream>
#include <conio.h>
#include <time.h>
#include "quasiNewton.h"

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
			der[i] = (x[j] + x_sum - b);
			if (i == j)
				der[i] *= 4.;
			else
				der[i] *= 2.;

			double prod_trunc = 1.;
			for (size_t k = 0; k < n; ++k)
				if (i != k)
					prod_trunc *= x[k];
			der[i] += 2. * (x_prod - 1.) * prod_trunc;
		}
	}
	return der;
}
int main()
{
	std::vector<double> xk(50, 0.);
	clock_t start = clock();
	bool success = optimize::BFGS(func, dfunc, xk);
	clock_t end = clock();
	double time = (double)(end - start) / CLOCKS_PER_SEC;
	std::cout.precision(15);
	std::cout << "x* = (";
	for (auto& x : xk)
		std::cout << x << ", ";
	std::cout << ")";
	std::cout << "Time: " << time;
	_getch();
}
