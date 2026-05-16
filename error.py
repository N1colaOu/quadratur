import quadrature as quad
import numpy.polynomial.polynomial as poly
import fractions as fr
import numpy as np
import matplotlib.pyplot as plt

def get_converg_rate(n, a = 0, b = 1, tol = 1e-10, open = True):
    if open:
        x = quad.get_oNC(n)
    else:
        x = quad.get_cNC(n)

    weights = quad.get_weights(n, x)
    power = fr.Fraction(1)
    rate = fr.Fraction()
    diff = fr.Fraction()

    while diff < tol: # if the difference is too big, then we assume they are different
        def f(x):
            return x**rate
        pol_int = fr.Fraction(b**power/power) - fr.Fraction(a**power/power) # integrated numerically
        appr = quad.calculate_quadrature(f, n, a, b, open=open)
        diff = np.abs(appr-pol_int)
        power += 1
        rate += 1

    return rate-1


    

def plot_converg(n_s = 1, n_e = 2, open = True, tol = 1e-10):
    n = []
    r_anal= []
    r_appr = []
    for i in range(n_s, n_e+1):
        n.append(i)
        r_anal.append(i+1)
        r_appr.append(get_converg_rate(i, open=False, tol = tol))
 
    plt.plot(n, r_anal, label="Analytical Convergence")
    plt.plot(n, r_appr, label="Calculated Convergence")
    plt.title("Convergnece Rate Analysis")
    plt.xlabel("Amount of Intervalls [n]")
    plt.ylabel("Convergence Rate [r]")
    plt.grid()
    plt.legend()
    plt.show()


