import numpy as np
import fractions as fr
import matplotlib.pyplot as plt

#d)
def gauss_legendre(n, a, b):

    m = n + 1 #wegen p_(n+1)
    #Nebendiagonale
    beta = np.array([ i / np.sqrt(4*i**2 - 1) for i in range(1, m) ])

    #Jacobi-Matrix
    A = np.diag(beta, 1) + np.diag(beta, -1)

    #Eigenwerte und Eigenvektoren
    eigenvalues, eigenvectors = np.linalg.eig(A)

    #Stützstellen
    nodes = eigenvalues

    #Gewichte
    weights = 2 * (eigenvectors[0, :] ** 2)

    #Mapping der Stützstellen von [-1, 1] auf [a, b]
    x_m = 0.5 * (b - a) * nodes + 0.5 * (a + b)
    w_m = 0.5 * (b - a) * weights

    return x_m, w_m


#Test der Gewichte und Stützstellen für n=2
x, w = gauss_legendre(2, -1, 1)
print('Stützstellen:', x)
print('Gewichte:', w)

def quadrature(f, n, a, b):
    x ,w = gauss_legendre(n, a, b)
    return np.sum(w * f(x))

#e)
#plotten um Genauigkeitsgrad 2n+1 zu verifizieren
#Fehler Matrix
def exact_integral(k):
    if k % 2 == 1:
        return 0
    else:
        return 2 / (k + 1)
#Funktion zum berechnen eines exakten Integrals für eine beliebige Funktion    

Nmax = 10
Kmax = 25

error = np.zeros((Nmax, Kmax))

for n in range(1, Nmax + 1):
    for k in range(Kmax):
        Q = quadrature(lambda x: x**k, n, -1, 1)
        I = exact_integral(k)
        error[n-1, k] = abs(I - Q)


#Plot
plt.figure(figsize=(8,5))
plt.imshow(
    np.log10(error + 1e-16),
    aspect='auto',
    origin='lower',
    cmap='viridis',
)
plt.colorbar(label='log10 error')
plt.xlabel('k (degree of monomial)')
plt.ylabel('n (number of nodes)')
plt.title('Error of Gauss-Legendre Quadrature')
plt.show()

#f)Visualierung von Gauß-Legendre für beliebeige Funktion n = 2
#Polynom 4. Grades
def f(x):
    return 0.3*x**4 - x**3 + 0.5*x**2 + x + 1 #Beispielpolynom

#Plot der Funktion
x_plot = np.linspace(-1, 1, 500)

plt.figure(figsize=(8,5))
#Gauß-Legendre Polynom
plt.plot(x_plot, f(x_plot), label='Gauß-Legendre', color='blue')

#Gauß-Legendre
x_nodes, w_nodes = gauss_legendre(2, -1, 1)

plt.scatter(x_nodes, f(x_nodes), color='red', label='Stützstellen n = 2')
plt.title('Gauß-Legendre Quadratur mit n=2')
plt.grid()
plt.legend()
plt.show()

#g)Untersuchen der Fehler für verschiedene Funktionen
#transformation der Funktion um das Intervall [-1,1] zu nutzen
def f1(x):
    return np.exp(x)

def f2(x):
    return np.sqrt(1-x**2)

def f3(x):
    return 1/((10)**(-2)+x**2)
I1 = np.exp(2) - 1
I2 = np.pi / 2
I3 = 20 * np.arctan(10)


n_vals = np.arange(1, 11)
err1 = []
err2 = []
err3 = []

for n in n_vals:
    err1.append(abs(I1 - quadrature(f1, n, 0, 2)))
    err2.append(abs(I2 - quadrature(f2, n, -1, 1)))
    err3.append(abs(I3 - quadrature(f3, n, -1, 1)))

#Plotten der Fehler
plt.figure(figsize=(8,5))
plt.semilogy(n_vals, err1, marker='o', label='f1(x) = exp(x)')
plt.semilogy(n_vals, err2, marker = 's', label = 'f2(x) = sqrt(1-x^2)')
plt.semilogy(n_vals, err3, marker = 'd', label = 'f3(x) = 1/((10)^(-2)+x^2))')

plt.xlabel('n (number of nodes)')
plt.ylabel('Absolute Error')
plt.title('Error of Gauss-Legendre Quadrature for Different Functions')
plt.grid(True, which='both')
plt.legend()
plt.show()

#summierte Quadraturregel
def composite_quadrature(f, a, b, N, n):
     #Teilintervalle
     intervals = np.linspace(a, b, N + 1)
     total = 0

     for i in range(N):
         
         ai = intervals[i]
         bi = intervals[i + 1]

         #Gauß-Legendre auf Teilintervall [ai, bi]
         x, w = gauss_legendre(n, ai, bi)
         total += np.sum(w * f(x))
         return total
     
I_exact = 20 * np.arctan(10)
N_vals = np.arange(1, 101)

err_n1 = []
err_n2 = []
err_n5 = []

for N in N_vals:
    Q1 = composite_quadrature(f3, -1, 1, N, 1)
    Q2 = composite_quadrature(f3, -1, 1, N, 2)
    Q5 = composite_quadrature(f3, -1, 1, N, 5)

    err_n1.append(abs(I_exact - Q1))
    err_n2.append(abs(I_exact - Q2))
    err_n5.append(abs(I_exact - Q5))

#Plotten der Fehler
plt.figure
plt.semilogy(N_vals, err_n1, marker='o', label='n=1')
plt.semilogy(N_vals, err_n2, marker='s', label='n=2')
plt.semilogy(N_vals, err_n5, marker='d', label='n=5')
plt.xlabel('N (number of subintervals)')
plt.ylabel('Absolute Error')
plt.title('Error of Composite Gauss-Legendre Quadrature for f3(x)')
plt.grid(True, which='both')
plt.legend()
plt.show()