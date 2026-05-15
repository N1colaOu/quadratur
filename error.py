import quadrature as quad
import numpy.polynomial.polynomial as poly
import numpy as np
import matplotlib.pyplot as plt

def get_converg_rate(n, a = 0, b = 1, tol = 1e-10, open = True):
    if open:
        x = quad.get_oNC(n)
    else:
        x = quad.get_cNC(n)

    weights = quad.get_weights(n, x)
    power = 1
    rate = 0
    diff = 0

    while diff < tol: # if the difference is too big, then we assume they are different
        def f(x):
            return x**rate
        pol_int = b**power/power - a**power/power # integrated numerically
        appr = quad.calculate_quadrature(f, n, a, b, open=open)
        diff = np.abs(appr-pol_int)
        power += 1
        rate += 1

    return rate-1


    

def plot_converg(n, r, open = True):
    pass


