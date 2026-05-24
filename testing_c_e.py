import error as err

err.plot_converg_NC(n_e=10) # we can see they are equal on the even n
err.plot_converg_GL(n_e=10, tol=1e-12) # they converge, with lowe tolearnce for higher n