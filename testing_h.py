import numpy as np
import quadrature as quad
import matplotlib.pyplot as plt

def f(x):
    return 1/(1e-2+x**2)

n = [1, 2, 5, 7] #n in the intervall self
N = 50 #intervalls
a = -1
b = 1

int_exact = 20*np.arctan(10)
sum_appr_oNC = []#quad.calculate_quadrature_sum_oNC(f, N, n, a, b)
sum_appr_cNC = []#quad.calculate_quadrature_sum_oNC(f, N, n, a, b)
sum_appr_GL = []#quad.calculate_quadrature_sum_oNC(f, N, n, a, b)
for i in n:
    sum_appr_oNC.append(quad.calculate_sum_quadrature_oNC(f, N, i, a, b))
    sum_appr_cNC.append(quad.calculate_sum_quadrature_cNC(f, N, i, a, b))
    sum_appr_GL.append(quad.calculate_sum_quadrature_GL(f, N, i, a, b))

err_oNC = np.abs(np.subtract(int_exact, sum_appr_oNC))
err_cNC = np.abs(np.subtract(int_exact, sum_appr_cNC))
err_GL = np.abs(np.subtract(int_exact, sum_appr_GL))

plt.figure(1)
plt.semilogy(n, err_oNC, label="Sum oNC")
plt.semilogy(n, err_cNC, label="Sum cNC")
plt.semilogy(n, err_GL, label="Sum GL")
plt.legend()
plt.grid()
plt.show()

