import numpy as np
import fractions as fr

def get_lagrange_pol_i(x, i): # get ith lagrange polynom
    n = len(x) - 1
    lagr_pol = [fr.Fraction(1)]
    for j in range(n + 1):
        if i  != j:
            free_coeff = fr.Fraction(-x[j]/(x[i] - x[j]))
            x_coeff = fr.Fraction(1/(x[i] - x[j]))
            to_mult = [free_coeff, x_coeff]
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

def get_lagrange_weights_i(li): # we get the ith weight by integrating the ith base pol
    n = len(li) - 1
    wi = 0
    for i in range(1, n+2, 2): # we go over only the odd indices, even ones are 0 (n+1 inclusive)
        to_add = fr.Fraction(li[i-1]/i)
        wi += to_add
    wi *= 2
    return float(wi)

def get_lagrange_weights(x): # we get all the weights
    n = len(x) - 1
    w = []
    l = get_lagrange_base(x)
    for i in range(n+1):
        w.append(get_lagrange_weights_i(l[i]))
    return w

def print_base(base): # help method for pretty printing
    for b in base:
        print(b)

#############################################################################################
def get_cNC(n, a=-1, b=1):
    #n is the amount of intervalls, x_0 is start pos (def -1) and x_n is end (def 1)
    x = np.linspace(a, b, n+1, dtype=float)
    omegas = get_lagrange_weights(x)
    return x, omegas

def get_oNC(n, a=-1, b=1):
    #analog to the closed version we divide the interval, but in more teilintervalls and we take the vals without both ends
    x = np.linspace(a, b, n+3, dtype=float)
    omegas = get_lagrange_weights(x[1:n+2])
    return x[1:n+2], omegas

def get_GL(n): #n as defined in the lectures is the last index of the points so: x_0 ... x_n, therfore n+1 points in total
    
    jacobian = np.zeros([n+1, n+1])    
    beta = [i/np.sqrt(4*i**2-1) for i in range(1, n+1)]

    jacobian += np.diag(beta, 1)
    jacobian += np.diag(beta, -1)
    
    ew, ev = np.linalg.eig(jacobian)
    v = ev[0]
    weights = 2*v.T*v
    return ew, weights

def calculate_quadrature(f, t, w, a, b):
    n = len(t) - 1
    t_wrapper = fr.Fraction((b-a)/2)*t + fr.Fraction((a+b)/2) #we bring the interval from [a,b] to [-1,1]
    res = fr.Fraction()
    for i in range(n+1):
        res += w[i] * f(t_wrapper[i])
    res *= fr.Fraction((b-a)/2) # we multiplty to keep consistent with the intervall
    return res # we return the result of the quadrature formula

def calculate_sum_quadrature(f, N, n, a, b, method):#gets the quadrature, by summing all N intervalls, in which [ab] is divided
    points = np.linspace(a, b, N)
    res = fr.Fraction()
    for i in range(N-1):
        a_new = points[i]
        b_new = points[i+1]
        t, w = method(n)
        res += calculate_quadrature(f, t, w, a_new, b_new)
    return res