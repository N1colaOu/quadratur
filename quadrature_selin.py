import numpy as np
import fractions as fr

# Task a)
def get_lagrange_weights(x):
    """
    Computes the Newton-Cotes weights for a given set of nodes.
    The indices of the list represent the exponents of x (e.g., index 0 is x^0).
    """
    n = len(x) - 1
    weights = []
    
    # Outer loop: compute each weight w_0, w_1, ..., w_n
    for i in range(n + 1):
        x_i = x[i]
        
        # Step 1: Compute the numeric denominator D_i
        denominator = fr.Fraction(1)
        for j in range(n + 1):
            if j != i:
                denominator *= (x_i - x[j])
                
        # Step 2: Build the numerator polynomial, starting with P(x) = 1
        numerator_poly = [fr.Fraction(1)]
        for j in range(n + 1):
            if j != i:
                # Inline polynomial multiplication with linear factor (x - x_j)
                # Shifting indices simulates multiplication by 'x'
                result = [fr.Fraction(0)] * (len(numerator_poly) + 1)
                for k in range(len(numerator_poly)):
                    result[k + 1] += numerator_poly[k]      # Multiply by x
                    result[k] += numerator_poly[k] * (-x[j]) # Multiply by -x_j
                numerator_poly = result
                
        # Step 3: Combine numerator and denominator to get Lagrange polynomial L_i(x)
        lagrange_poly = [c / denominator for c in numerator_poly]
        
        # Step 4: Integrate L_i(x) over [-1, 1] using symmetry
        # Odd exponents integrate to 0 and are skipped via 'k % 2 == 0'
        integral_value = fr.Fraction(0)
        for k in range(len(lagrange_poly)):
            if k % 2 == 0:
                integral_value += lagrange_poly[k] * fr.Fraction(2, k + 1)
                
        # Convert to float for compatibility with the team's system
        weights.append(float(integral_value))
        
    return weights

def get_cNC(n):
    """Generates closed Newton-Cotes nodes and weights on [-1, 1]"""
    x = np.linspace(-1, 1, n + 1, dtype=float)
    omegas = get_lagrange_weights(x)
    return x, omegas

def get_oNC(n):
    """Generates open Newton-Cotes nodes and weights on [-1, 1]"""
    x = np.linspace(-1, 1, n + 3, dtype=float)
    omegas = get_lagrange_weights(x[1:n + 2])
    return x[1:n + 2], omegas

def get_LG(n):
    """Task d): Fetches optimal Gauss-Legendre nodes and weights on [-1, 1] using NumPy"""
    nodes, weights = np.polynomial.legendre.leggauss(n)
    return nodes, weights

def calculate_quadrature(f, t, w, a, b):
    """
    Task b): Computes the integral over [a, b] using the universal interval transformation.
    Maps nodes from [-1, 1] to [a, b] and scales the weights accordingly.
    """
    factor = float(b - a) / 2.0
    integral_sum = 0.0
    
    for i in range(len(t)):
        # Linear transformation formula: mapping nodes and scaling weights
        x_transformed = float(a) + factor * (float(t[i]) + 1.0)
        w_transformed = float(w[i]) * factor
        
        # Basic quadrature formula: sum of (w_i * f(x_i))
        integral_sum += w_transformed * f(x_transformed)
        
    return integral_sum
