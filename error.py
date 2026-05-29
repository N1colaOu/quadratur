import quadrature as quad
import fractions as fr
import numpy as np
import matplotlib.pyplot as plt

def get_converg_rate(t, w, a = 0, b = 1, tol = 1e-10):

    power = int(1)
    rate = int(0)
    diff = fr.Fraction()

    while diff < tol: # if the difference is too big, then we assume they are different
        f = lambda x : x**rate
        pol_int = fr.Fraction(b**power/power) - fr.Fraction(a**power/power) # integrated exactly
        appr = quad.calculate_quadrature(f, t, w, a, b)
        diff = fr.Fraction(np.abs(appr-pol_int))
        power += 1
        rate += 1

    return rate-2

def plot_converg(n, method, tol = 1e-10):#one convergence func, dependant on the method
    if method == quad.get_oNC:
        title = "Open Newton Cotes"
        conv = lambda x: x+1
    elif method == quad.get_cNC:
        title = "Closed Newton Cotes"
        conv = lambda x: x+1
    elif method == quad.get_GL:
        title = "Gauss-Legendre"
        conv = lambda x: 2*x+1
    else:
        pass

    n_arr = []
    r_theory= []
    r_appr = []
    for i in range(1, n+1):
        t, w = method(i)
        n_arr.append(i)
        r_theory.append(conv(i))
        r_appr.append(get_converg_rate(t, w, tol = tol))
 
    plt.figure()
    plt.plot(n_arr, r_theory, label="Analytical Convergence")
    plt.plot(n_arr, r_appr, label="Calculated Convergence")
    plt.title("Convergnece Rate Analysis of " + title)
    plt.xlabel("Amount of Intervalls [n]")
    plt.ylabel("Convergence Rate [r]")
    plt.grid()
    plt.legend()



