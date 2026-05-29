import error as err
import quadrature as quad
import matplotlib.pyplot as plt

n = 10

err.plot_converg(n, method = quad.get_oNC) # we can see they are equal on the even n
err.plot_converg(n, method = quad.get_cNC) # we can see they are equal on the even n
err.plot_converg(n, method = quad.get_GL, tol=1e-12) # they converge, with lowe tolearnce for higher n
plt.show()