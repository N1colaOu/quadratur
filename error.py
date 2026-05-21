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


    

def plot_converg_NC(n_s = 1, n_e = 2, tol = 1e-10):
    n = []
    r_anal= []
    r_appr_oNC = []
    r_appr_cNC = []
    for i in range(n_s, n_e+1):
        t_oNC, w_oNC = quad.get_oNC(i)
        t_cNC, w_cNC = quad.get_cNC(i)
        n.append(i)
        r_anal.append(i+1)
        r_appr_oNC.append(get_converg_rate(t_oNC, w_oNC, tol = tol))
        r_appr_cNC.append(get_converg_rate(t_cNC, w_cNC, tol = tol))
 
    plt.plot(n, r_anal, label="Analytical Convergence")
    plt.plot(n, r_appr_oNC, label="Calculated Convergence Open")
    plt.plot(n, r_appr_cNC, label="Calculated Convergence Closed")
    plt.title("Convergnece Rate Analysis")
    plt.xlabel("Amount of Intervalls [n]")
    plt.ylabel("Convergence Rate [r]")
    plt.grid()
    plt.legend()
    plt.show()

def plot_converg_LG(n_s = 1, n_e = 2, tol = 1e-10):
    n = []
    r_anal= []
    r_appr_LG = []
    for i in range(n_s, n_e+1):
        t_LG, w_LG = quad.get_LG(i)
        n.append(i)
        r_anal.append(2*i+1)
        r_appr_LG.append(get_converg_rate(t_LG, w_LG, tol = tol))
 
    plt.plot(n, r_anal, label="Analytical Convergence")
    plt.plot(n, r_appr_LG, label="Calculated Convergence Legendre")
    plt.title("Convergnece Rate Analysis")
    plt.xlabel("Amount of Intervalls [n]")
    plt.ylabel("Convergence Rate [r]")
    plt.grid()
    plt.legend()
    plt.show()


