import fractions as fr
import numpy.polynomial.polynomial as poly
import numbers as num

def get_lagrange_pol_i(n, x, i):
    lagr_pol = 1
    for j in range(n + 1):
        if i  != j:
            x_coeff =  -x[j] / (x[i] - x[j])
            """fr.Fraction(-x[j], (x[i] - x[j]))"""
            free_coeff = 1 / (x[i] - x[j])
            """fr.Fraction(1, (x[i] - x[j]))""" 
            to_mult =  poly.Polynomial([x_coeff, free_coeff])
            lagr_pol = poly.polymul(to_mult, lagr_pol)
    return poly.Polynomial(lagr_pol)

def get_lagrange_base(n, x):
    lagr_base = []
    for i in range(n+1):
        lagr_base.append(get_lagrange_pol_i(n, x, i))
    return lagr_base

def get_weights_i(n, li):
    wi = fr.Fraction()
    for i in range(1, n+1, 2):
        to_add = fr.Fraction((li[i] - 1),i)
        wi += to_add
    wi *= 2
    return wi

def get_weights(n, l):
    w = []
    for i in range(n+1):
        w.append(get_weights_i(n, l[i]))
    return w

def print_base(base):
    for b in base:
        print(b)