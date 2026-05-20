import numpy as np
import fractions as fr

#d)
def gauss_legendre(n):

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

    return nodes, weights


#Test der Gewichte und Stützstellen für n=2
x, w = gauss_legendre(2)
print('Stützstellen:', x)
print('Gewichte:', w)

def integrate(f, n):
    x ,w = gauss_legendre(n)
    return np.sum(w * f(x))

#Test
f = lambda x: x**2
result = integrate(f, 2)
print('Integralwert:', result)
