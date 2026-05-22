import quadrature_selin as quad
import numpy as np

def get_converg_rate(t, w, a=0, b=1, tol=1e-12):
    """
    Task c): Determines the degree of precision by testing monomials x^d over [a, b].
    Utilizes an infinite loop that breaks via 'return' once the numerical error exceeds tolerance.
    """
    d = 0
    while True:
        # Dynamically define the monomial function f(x) = x^d
        monomial = lambda x: x**d
        
        # Calculate numerical integral using the custom quadrature module
        num_integral = quad.calculate_quadrature(monomial, t, w, a, b)
        
        # Calculate exact analytical integral value: (b^(d+1) - a^(d+1)) / (d+1)
        exact_integral = (b**(d + 1) - a**(d + 1)) / (d + 1)
        error = abs(num_integral - exact_integral)
        
        # If the numerical result is no longer exact, the breaking point is reached
        if error > tol:
            # Returns the last successfully integrated exponent
            return d - 1
            
        # If successful, increment exponent to test the next higher monomial degree
        d += 1
