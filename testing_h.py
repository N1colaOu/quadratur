import numpy as np
import quadrature as quad
import matplotlib.pyplot as plt

def f(x):
    return 1/(1e-2+x**2)

n = [1, 2, 5, 7] #n in the intervall self
N = 75 #intervalls
a = -1
b = 1

int_exact = 20*np.arctan(10)
sum_appr_oNC = []
sum_appr_cNC = []
sum_appr_GL = []
for i in n:
    sum_appr_oNC.append(quad.calculate_sum_quadrature(f, N, i, a, b, method=quad.get_oNC)) # the sum for open newton, depending on n
    sum_appr_cNC.append(quad.calculate_sum_quadrature(f, N, i, a, b, method=quad.get_cNC))
    sum_appr_GL.append(quad.calculate_sum_quadrature(f, N, i, a, b, method=quad.get_GL))

err_oNC = np.abs(np.subtract(int_exact, sum_appr_oNC)) # the error for each n of the sums
err_cNC = np.abs(np.subtract(int_exact, sum_appr_cNC))
err_GL = np.abs(np.subtract(int_exact, sum_appr_GL))

plt.figure(1)
plt.semilogy(n, err_oNC, label="Sum oNC") # log plots
plt.semilogy(n, err_cNC, label="Sum cNC")
plt.semilogy(n, err_GL, label="Sum GL")
plt.legend()
plt.xlabel("n")
plt.ylabel("Error")
plt.title("Semi-Log Error")
plt.grid()
plt.show()

