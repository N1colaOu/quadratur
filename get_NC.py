import numpy as np

def get_cNC(n, x_0 = -1, x_n = 1):
    #n is the amount of intervalls, x_0 is start pos (def -1) and x_n is end (def 1)
    x = np.linspace(x_0, x_n, n+1)
    return x

def get_oNC(n, x_0 = -1, x_n = 1):
    #analog to the closed version we divide the interval, but in more teilintervalls and we take the vals without both ends
    x = np.linspace(x_0, x_n, n+3)
    return x[1:n+2]