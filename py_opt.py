import numpy as np
from scipy.optimize import minimize

def f(x):
	n = len(x)
	f_sum = 0.0
	b = n + 1
	for i in range (0, n - 1):
		f_sum += (x[i] + sum(x) - b)**2
	f_sum += (np.prod(x) - 1.0)**2
	return f_sum
    
def jac(x):
	n = len(x)
	b = n + 1
	der = np.zeros(n)
	for i in range (0, n):
		for j in range(0, n-1):
			if i == j:
				der[i] +=4.0*(x[j] + sum(x) - b)
			else:
				der[i] += 2.0*(x[j] + sum(x) - b)
		prod = 1.0
		for j in range(0, n):
			if i != j:
				prod *= x[j]
		der[i] += 2.0*(np.prod(x) - 1)*prod
	return der

x0 = np.zeros(10)
for i in range(0, len(x0)):
	x0[i] = 0.12

x0[2] = 1.7
x0[6] = 3.25

res = minimize(f, x0, method='BFGS', jac=jac,
               options={'disp': True}, tol=1.e-14)
print(res)

for i in range(0,  len(x0)):
	print(str(res.x[i] ) )