import numpy as np
import quadrature as quad

def f(x):
    return 1/(1e-2+x**2)

n = [1, 2, 5] #n in the intervall self
N = 50 #intervalls

sum_appr_oNC = quad.calculate_quadrature_sum_oNC(f, N, n, a, b)
sum_appr_cNC = quad.calculate_quadrature_sum_oNC(f, N, n, a, b)
sum_appr_GL = quad.calculate_quadrature_sum_oNC(f, N, n, a, b)
