import numpy as np
import fractions as fr
import numpy.polynomial.polynomial as poly

def get_cNC(n):
    #n is the amount of intervalls, x_0 is start pos (def -1) and x_n is end (def 1)
    x = np.linspace(-1, 1, n+1, dtype=float)
    return x

def get_oNC(n):
    #analog to the closed version we divide the interval, but in more teilintervalls and we take the vals without both ends
    x = np.linspace(-1, 1, n+3, dtype=float)
    return x[1:n+2]

def get_lagrange_pol_i(n, x, i): # get ith lagrange polynom
    lagr_pol = 1
    for j in range(n + 1):
        if i  != j:
            x_coeff = fr.Fraction(-x[j]/(x[i] - x[j]))
            free_coeff = fr.Fraction(1/(x[i] - x[j]))
            to_mult =  poly.Polynomial([x_coeff, free_coeff])
            lagr_pol = poly.polymul(to_mult, lagr_pol) # we keep updating the product with the next term needed to be multiplied
    return poly.Polynomial(lagr_pol)

def get_lagrange_base(n, x): # get the whole base of polynomials
    lagr_base = []
    for i in range(n+1): 
        lagr_base.append(get_lagrange_pol_i(n, x, i))
    return lagr_base

def get_weights_i(n, li): # we get the ith weight by integrating the ith base pol
    wi = 0
    coeffs = li.convert().coef
    for i in range(1, n+2, 2): # we go over only the odd indices, even ones are 0 (n+1 inclusive)
        to_add = fr.Fraction(coeffs[i-1]/i)
        wi += to_add
    wi *= 2
    return wi

def get_weights(n, x): # we get all the weights
    w = []
    l = get_lagrange_base(n, x)
    for i in range(n+1):
        w.append(get_weights_i(n, l[i]))
    return w

def print_base(base): # help method for pretty printing
    for b in base:
        print(b)

def calculate_quadrature(f, n, a, b, open = True):
    if open:
        t = get_oNC(n)
    else:
        t = get_cNC(n)
    t_wrapper = (b-a)/2 * t + (a+b)/2
    omegas = get_weights(n, t)
    res = 0
    for i in range(n+1):
        res += omegas[i] * f(t_wrapper[i])
    res *= (b-a)/2
    return res