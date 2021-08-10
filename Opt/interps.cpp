#pragma once
#include "quasiNewton.h"
// by https://github.com/yixuan/

double optimize::MT_quadratic_interp(const double& a, const double& b, const double& fa, const double& ga, const double& fb)
{
	double ba = b - a;
	return a + 0.5 * ba * ba * ga / (fa - fb + ba * ga);
}

double optimize::MT_secant(const double& a, const double& b, const double& ga, const double& gb)
{
	return b - (b - a) * gb / (gb - ga);
}

double optimize::MT_cubic_interp(const double& a, const double& b, const double& fa, const double& fb, const double& ga, const double& gb)
{
	if (fabs(a - b) < std::numeric_limits<double>::epsilon())
		return a;

	const double ba = b - a;
	const double ba2 = ba * ba;
	const double ba3 = ba2 * ba;
	const double fba = fb - fa;
	const double z = (ga + gb) * ba - 2. * fba;
	const double w = fba * ba - ga * ba2;

	// If c3 = z/(b-a)^3 == 0, reduce to quadratic problem
	const double endmin = (fa < fb) ? a : b;
	if (fabs(z) < std::numeric_limits<double>::epsilon())
	{
		const double c2 = fba / ba2 - ga / ba;
		const double c1 = fba / ba - (a + b) * c2;
		// Global minimum, can be infinity
		const double globmin = -c1 / (2. * c2);
		// If c2 <= 0, or globmin is outside [a, b], then the minimum is achieved at one end point
		return (c2 > 0. && globmin >= a && globmin <= b) ? globmin : endmin;
	}

	// v = c1 / c2
	const double v = (-2. * a * w + ga * ba3 + a * (a + 2. * b) * z) /
		(w - (2. * a + b) * z);
	// u = c2 / (3 * c3), may be very large if c3 ~= 0
	const double u = (w / z - (2. * a + b)) / 3.;
	// q'(x) = c1 + 2 * c2 * x + 3 * c3 * x^2 = 0
	// x1 = -u * (1 + sqrt(1 - v/u))
	// x2 = -u * (1 - sqrt(1 - v/u)) = -v / (1 + sqrt(1 - v/u))

	// If q'(x) = 0 has no solution in [a, b], q(x) is monotone in [a, b]
	// Case I: no solution globally, 1 - v/u <= 0
	if (v / u >= 1.)
		return endmin;
	// Case II: no solution in [a, b]
	const double vu = 1. + sqrt(1. - v / u);
	const double sol1 = -u * vu;
	const double sol2 = -v / vu;
	if ((sol1 - a) * (sol1 - b) >= 0. && (sol2 - a) * (sol2 - b) >= 0.)
		return endmin;

	// Now at least one solution is in (a, b)
	// Check the second derivative
	// q''(x) = 2 * c2 + 6 * c3 * x;
	const double c3 = z / ba3;
	const double c2 = 3. * c3 * u;
	const double qpp1 = 2. * c2 + 6. * c3 * sol1;
	const double sol = (qpp1 > 0.) ? sol1 : sol2;
	// If the local minimum is not in [a, b], return one of the end points
	if ((sol - a) * (sol - b) >= 0.)
		return endmin;

	// Compare the local minimum to the end points
	const double c1 = v * c2;
	const double fsol = fa + c1 * (sol - a) + c2 * (sol * sol - a * a) +
		c3 * (sol * sol * sol - a * a * a);
	return (fsol < std::min(fa, fb)) ? sol : endmin;
}
