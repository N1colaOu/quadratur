import quadrature_selin as quad
import error_selin as err
import matplotlib.pyplot as plt
import numpy as np

# Test for tasks b) & d)
print("--- Integration Test (f(x) = x^2 from 0 to 2) ---")
f_test = lambda x: x**2

# Test closed Newton-Cotes (n=2, Simpson's rule)
t_c, w_c = quad.get_cNC(2)
print("NC Closed (n=2):", quad.calculate_quadrature(f_test, t_c, w_c, 0, 2))

# Test Gauss-Legendre (n=2, integrates degree 2 exactly with only 2 nodes)
t_g, w_g = quad.get_LG(2)
print("Gauss-Legendre (n=2):", quad.calculate_quadrature(f_test, t_g, w_g, 0, 2))


# Plotting for c)
n_values = list(range(1, 9))
closed_precisions = []
open_precisions = []

# Collect the algebraic degree of precision for each n
for n in n_values:
    t_o, w_o = quad.get_oNC(n)
    t_c, w_c = quad.get_cNC(n)
    
    closed_precisions.append(err.get_converg_rate(t_c, w_c))
    open_precisions.append(err.get_converg_rate(t_o, w_o))

# Generate the matplotlib diagram
plt.figure(figsize=(10, 6))
plt.plot(n_values, closed_precisions, 'o-', label='Closed Newton-Cotes', color='blue', linewidth=2)
plt.plot(n_values, open_precisions, 's--', label='Open Newton-Cotes', color='red', linewidth=2)

plt.title('Degree of Precision for Newton-Cotes Formulas', fontsize=14)
plt.xlabel('Degree of the formula (n)', fontsize=12)
plt.ylabel('Degree of Precision (d)', fontsize=12)
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend(fontsize=11)

# Save the output visualization as a PNG image
plt.savefig('degree_of_precision.png')
print("\nPlot successfully saved as 'degree_of_precision.png'!")


# Plotting for task e)
print ("\nGenerating convergence plot for task e) ")

# 1. Define a non-polynomial function and its exact analytical integral value
f_conv = lambda x: np.exp(x)
exact_val = np.exp(1.0) - 1.0

# 2. Set the range of n (number of nodes) we want to evaluate
n_range = list(range(1, 9))
error_cNC = []
error_LG = []

# 3. Calculate the absolute error for each n
for n in n_range:
    # Fetch nodes and weights for both methods
    t_c, w_c = quad.get_cNC(n)
    t_g, w_g = quad.get_LG(n)
    
    # Compute the numerical integrals using your task b) function
    res_cNC = quad.calculate_quadrature(f_conv, t_c, w_c, 0, 1)
    res_LG = quad.calculate_quadrature(f_conv, t_g, w_g, 0, 1)
    
    # Calculate absolute differences (errors) and append them to our lists
    error_cNC.append(abs(res_cNC - exact_val))
    error_LG.append(abs(res_LG - exact_val))

# 4. Create a log-scaled visualization plot
plt.figure(figsize=(10, 6))

# Use semilogy because errors drop exponentially; a normal scale would look flat
plt.semilogy(n_range, error_cNC, 'o-', label='Closed Newton-Cotes', color='blue', linewidth=2)
plt.semilogy(n_range, error_LG, 's--', label='Gauss-Legendre', color='green', linewidth=2)

# 5. Add professional scientific chart configurations
plt.title('Convergence Analysis: Newton-Cotes vs. Gauss-Legendre ($f(x)=e^x$)', fontsize=14)
plt.xlabel('Number of intervals / parameters (n)', fontsize=12)
plt.ylabel('Absolute Quadrature Error (Log Scale)', fontsize=12)
plt.grid(True, which="both", linestyle=':', alpha=0.6)
plt.legend(fontsize=11)

# Save the second chart to your repository
plt.savefig('convergence_comparison.png')
print("Plot successfully saved as 'convergence_comparison.png'!")


# Task f) - case study using a polynomial of degree 4

print("\n--- Running study for task f) ---")

# Define the polynomial coefficients globally so they are easy to change during the Testat
# p(x) = c4*x^4 + c3*x^3 + c2*x^2 + c1*x + c0
COEFFS = [1.0, 1.0, 1.0, 2.0, 1.0] # Order: [c0, c1, c2, c3, c4]

def p_favorite(x):
    """Evaluates our favorite degree 4 polynomial using the flexible COEFFS list."""
    return COEFFS[4]*x**4 + COEFFS[3]*x**3 + COEFFS[2]*x**2 + COEFFS[1]*x + COEFFS[0]

def exact_p_integral():
    """Computes the exact analytical integral of p(x) over [-1, 1] using calculus."""
    # Integral of x^k over [-1, 1] is 2/(k+1) for even k, and 0 for odd k
    val = COEFFS[0]*2.0 + COEFFS[2]*(2.0/3.0) + COEFFS[4]*(2.0/5.0)
    return val

def naive_interpolation(x_nodes, y_nodes, x_fine):
    """Reconstructs and evaluates the interpolation polynomial using a Vandermonde matrix."""
    A = np.vander(x_nodes, increasing=True)
    poly_coeffs = np.linalg.solve(A, y_nodes)
    
    y_fine = np.zeros(len(x_fine))
    for i in range(len(poly_coeffs)):
        y_fine += poly_coeffs[i] * (x_fine ** i)
    return y_fine

# Set up parameters required by the task
n_f = 2
x_plot = np.linspace(-1, 1, 1000)
y_exact_plot = p_favorite(x_plot)

# 1. Fetch standard nodes and weights on [-1, 1]
t_c, w_c = quad.get_cNC(n_f)
t_o, w_o = quad.get_oNC(n_f)
t_g, w_g = quad.get_LG(n_f + 1) # Using 3 nodes to compare fairly with n+1 nodes

# 2. Compute numeric quadrature values over [-1, 1]
quad_cNC = quad.calculate_quadrature(p_favorite, t_c, w_c, -1, 1)
quad_oNC = quad.calculate_quadrature(p_favorite, t_o, w_o, -1, 1)
quad_LG  = quad.calculate_quadrature(p_favorite, t_g, w_g, -1, 1)
true_int = exact_p_integral()

# Print out the results for the verbal comparison
print(f"Exact Analytical Integral:  {true_int}")
print(f"Closed Newton-Cotes Result: {quad_cNC}  (Error: {abs(quad_cNC - true_int)})")
print(f"Open Newton-Cotes Result:   {quad_oNC}  (Error: {abs(quad_oNC - true_int)})")
print(f"Gauss-Legendre Result:      {quad_LG}  (Error: {abs(quad_LG - true_int)})")

# 3. Generate the polynomial interpolation visual trajectories
y_inter_cNC = naive_interpolation(t_c, p_favorite(t_c), x_plot)
y_inter_oNC = naive_interpolation(t_o, p_favorite(t_o), x_plot)
y_inter_LG  = naive_interpolation(t_g, p_favorite(t_g), x_plot)

# 4. Plotting
plt.figure(figsize=(11, 7))
plt.plot(x_plot, y_exact_plot, 'k-', label='Original $p(x)$ (Degree 4)', linewidth=3)

plt.plot(x_plot, y_inter_cNC, 'b-', label='Closed NC Interpolation', alpha=0.8)
plt.scatter(t_c, p_favorite(t_c), color='blue', s=60, zorder=5)

plt.plot(x_plot, y_inter_oNC, 'r--', label='Open NC Interpolation', alpha=0.8)
plt.scatter(t_o, p_favorite(t_o), color='red', s=60, zorder=5)

plt.plot(x_plot, y_inter_LG, 'g-.', label='Gauss-Legendre Interpolation', alpha=0.8)
plt.scatter(t_g, p_favorite(t_g), color='green', s=60, zorder=5)

plt.title('Task f): Polynomial Interpolation Trajectories vs Original $p(x)$', fontsize=14)
plt.xlabel('x', fontsize=12)
plt.ylabel('y', fontsize=12)
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend(fontsize=11)

plt.savefig('polynomial_interpolation.png')
print("\nInterpolation plot successfully saved as 'polynomial_interpolation.png'!")


