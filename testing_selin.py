import quadrature_selin as quad
import error_selin as err
import matplotlib.pyplot as plt


# Test for tasks b) & d)
print("--- SELIN'S INTEGRATION TEST (f(x) = x^2 from 0 to 2) ---")
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
