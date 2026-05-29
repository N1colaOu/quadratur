import numpy as np
import matplotlib.pyplot as plt
import quadrature as quad

def f1(x):
    return np.exp(x)
def f2(x):
    return np.sqrt(1-x**2)
def f3(x):
    return 1/(1e-2+x**2)

a1, b1 = [0, 2]
a2, b2 = [-1, 1]
a3, b3 = [-1, 1]
n = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] # all n-s to calculate

int1_appr_oNC = [] # each function and each quadrature
int1_appr_cNC = []
int1_appr_GL = []

int2_appr_oNC = []
int2_appr_cNC = []
int2_appr_GL = []

int3_appr_oNC = []
int3_appr_cNC = []
int3_appr_GL = []
for i in n:
    x_oNC, w_oNC = quad.get_oNC(i)
    x_cNC, w_cNC = quad.get_cNC(i)
    x_GL, w_GL = quad.get_GL(i)
    int1_appr_oNC.append(quad.calculate_quadrature(f1, x_oNC, w_oNC, a1, b1))
    int1_appr_cNC.append(quad.calculate_quadrature(f1, x_oNC, w_oNC, a1, b1))
    int1_appr_GL.append(quad.calculate_quadrature(f1, x_oNC, w_oNC, a1, b1))
    
    int2_appr_oNC.append(quad.calculate_quadrature(f2, x_oNC, w_oNC, a2, b2))
    int2_appr_cNC.append(quad.calculate_quadrature(f2, x_oNC, w_oNC, a2, b2))
    int2_appr_GL.append(quad.calculate_quadrature(f2, x_oNC, w_oNC, a2, b2))
    
    int3_appr_oNC.append(quad.calculate_quadrature(f3, x_oNC, w_oNC, a3, b3))
    int3_appr_cNC.append(quad.calculate_quadrature(f3, x_oNC, w_oNC, a3, b3))
    int3_appr_GL.append(quad.calculate_quadrature(f3, x_oNC, w_oNC, a3, b3))

int1_exact = np.exp(2)-1
int2_exact = np.pi/2
int3_exact = 20*np.arctan(10)

err1_oNC = np.abs(np.subtract(int1_exact, int1_appr_oNC))
err1_cNC = np.abs(np.subtract(int1_exact, int1_appr_cNC))
err1_GL = np.abs(np.subtract(int1_exact, int1_appr_GL))

err2_oNC = np.abs(np.subtract(int2_exact, int2_appr_oNC))
err2_cNC = np.abs(np.subtract(int2_exact, int2_appr_cNC))
err2_GL = np.abs(np.subtract(int2_exact, int2_appr_GL))

err3_oNC = np.abs(np.subtract(int3_exact, int3_appr_oNC))
err3_cNC = np.abs(np.subtract(int3_exact, int3_appr_cNC))
err3_GL = np.abs(np.subtract(int3_exact, int3_appr_GL))

plt.figure(1)
plt.semilogy(n, err1_oNC, label="oNC-1")
plt.semilogy(n, err1_cNC, label="cNC-1")
plt.semilogy(n, err1_GL, label="GL-1")
plt.legend()
plt.grid()
plt.title("Error Estimate for 1. Function")

plt.figure(2)
plt.semilogy(n, err2_oNC, label = "oNC-2")
plt.semilogy(n, err2_cNC, label="cNC-2")
plt.semilogy(n, err2_GL, label="GL-2")
plt.legend()
plt.grid()
plt.title("Error Estimate for 2. Function")

plt.figure(3)
plt.semilogy(n, err3_oNC, label="oNC-3")
plt.semilogy(n, err3_cNC, label="cNC-3")
plt.semilogy(n, err3_GL, label="GL-3")
plt.legend()
plt.grid()
plt.title("Error Estimate for 3. Function")

plt.show()