import numpy as np
import fractions as fr

def get_cNC(n):
    #n is the amount of intervalls, x_0 is start pos (def -1) and x_n is end (def 1)
    x = np.linspace(-1, 1, n+1, dtype=float)
    return x

def get_oNC(n):
    #analog to the closed version we divide the interval, but in more teilintervalls and we take the vals without both ends
    x = np.linspace(-1, 1, n+3, dtype=float)
    return x[1:n+2]

def get_LGC(n):
    pass

def get_lagrange_pol_i(x, i): # get ith lagrange polynom
    n = len(x) - 1
    lagr_pol = [fr.Fraction(1)]
    for j in range(n + 1):
        if i  != j:
            x_coeff = fr.Fraction(-x[j]/(x[i] - x[j]))
            free_coeff = fr.Fraction(1/(x[i] - x[j]))
            to_mult = [x_coeff, free_coeff]
            lagr_pol = polmult(lagr_pol, to_mult)
    return lagr_pol

def polmult(coeffs1, coeffs2): # we multiply both polynomials
    m = len(coeffs1)
    n = len(coeffs2)
    prod = [fr.Fraction()] * (m + n - 1)
    for i in range(m):
        for j in range(n):
            prod[i + j] += coeffs1[i] * coeffs2[j]
    return prod

def get_lagrange_base(x): # get the whole base of polynomials
    n = len(x) - 1
    lagr_base = []
    for i in range(n+1): 
        lagr_base.append(get_lagrange_pol_i(x, i))
    return lagr_base

def get_weights_i(li): # we get the ith weight by integrating the ith base pol
    n = len(li) - 1
    wi = 0
    for i in range(1, n+2, 2): # we go over only the odd indices, even ones are 0 (n+1 inclusive)
        to_add = fr.Fraction(li[i-1]/i)
        wi += to_add
    wi *= 2
    return wi

def get_weights(x): # we get all the weights
    n = len(x) - 1
    w = []
    l = get_lagrange_base(x)
    for i in range(n+1):
        w.append(get_weights_i(l[i]))
    return w

def print_base(base): # help method for pretty printing
    for b in base:
        print(b)

def calculate_quadrature(f, t, a, b):
    n = len(t) - 1
    t_wrapper = fr.Fraction((b-a)/2)*t + fr.Fraction((a+b)/2) #we bring the interval from [a,b] to [-1,1]
    omegas = get_weights(t)
    res = fr.Fraction()
    for i in range(n+1):
        res += omegas[i] * f(t_wrapper[i])
    res *= fr.Fraction((b-a)/2) # we multiplty to keep consistent with the intervall
    return res # we return the result of the quadrature formula