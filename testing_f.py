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


n = 2
a0 = 1
a1 = 2
a2 = 3
a3 = 4
a4 = 5

pol_4 = lambda x : (x**4)*a4 + (x**3)*a3 + (x**2)*a2 + (x)*a1 + a0 # our polynom

n_plot = 1000 # used to visualize our polynoms
x_plot = np.linspace(-1, 1, n_plot)
x_GL, w_GL = quad.get_GL(n)
x_oNC, w_oNC = quad.get_oNC(n)
x_cNC, w_cNC = quad.get_cNC(n)
y_oNC = pol_4(x_oNC)
y_inter_oNC = naive_inter(x_oNC, y_oNC, x_plot) # we interpolate each set of points to see what we actually integrate in the end
y_cNC = pol_4(x_cNC)
y_inter_cNC = naive_inter(x_cNC, y_cNC, x_plot)
y_GL = pol_4(x_GL)
y_inter_GL = naive_inter(x_GL, y_GL, x_plot)

plt.grid()
plt.plot(x_plot, y_inter_oNC, label="Open Newton")
plt.scatter(x_oNC, y_oNC)
plt.plot(x_plot, y_inter_cNC, label="Closed Newton")
plt.scatter(x_cNC, y_cNC)
plt.plot(x_plot, y_inter_GL, label="Legendre Gauss")
plt.scatter(x_GL, y_GL)
plt.legend()
plt.show()

int_exact = a4/5 + a2/3 + a0 # exact integration of our polynom
int_exact *= 2
int_appr_oNC = quad.calculate_quadrature(pol_4, x_oNC, w_oNC, -1, 1)
int_appr_cNC = quad.calculate_quadrature(pol_4, x_cNC, w_cNC, -1, 1)
int_appr_GL = quad.calculate_quadrature(pol_4, x_GL, w_GL, -1, 1)
err_oNC = np.abs(int_exact - int_appr_oNC)
err_cNC = np.abs(int_exact - int_appr_cNC)
err_GL = np.abs(int_exact - int_appr_GL)
print("Error Open Newton Cotes: ", err_oNC)
print("Error Closed Newton Cotes: ", err_cNC)
print("Error Legendre-Gauss: ", err_GL)



