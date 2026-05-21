import numpy as np
import quadrature as quad
import error as err
import matplotlib.pyplot as plt

def naive_inter(x, y, xi):#copied from last HW
    n = len(x)
    assert n == len(y) #make sure both inputs are the same size

    A = np.vander(x, increasing=True) #create vandermonde matrix
    coeffs = np.linalg.solve(A, y) #solve lin sys
    m = len(xi)
    
    yi = np.zeros(m)
    for i in range(n):
        yi += coeffs[i] * pow(xi, i) #evaluate polynom at xi
        
    return yi


def pol_4(x):
    return (x**4)*1 + (x**3)*2 + (x**2)*1 + (x)*1 + 1 

n = 2
n_plot = 1000
x_plot = np.linspace(-1, 1, n_plot)

x_LG, w_LG = quad.get_LG(n)
x_oNC, w_oNC = quad.get_oNC(n)
x_cNC, w_cNC = quad.get_cNC(n)

y_oNC = pol_4(x_oNC)
y_inter_oNC = naive_inter(x_oNC, y_oNC, x_plot)
y_cNC = pol_4(x_cNC)
y_inter_cNC = naive_inter(x_cNC, y_cNC, x_plot)
y_LG = pol_4(x_LG)
y_inter_LG = naive_inter(x_LG, y_LG, x_plot)

plt.grid()
plt.plot(y_inter_oNC, label="Open Newton")
plt.scatter(x_oNC, y_oNC)
plt.plot(y_inter_cNC, label="Closed Newton")
plt.scatter(x_cNC, y_cNC)
plt.plot(y_inter_LG, label="Legendre Gauss")
plt.scatter(x_LG, y_LG)
plt.legend()
plt.show()
