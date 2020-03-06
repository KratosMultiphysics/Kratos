"""This module contains the Newton-Raphson algorithm.

Author: Armin Geiser
"""

import numpy.linalg as la
import scipy.sparse.linalg as sla

def newton_raphson_solve(calculate_system, x_initial, max_iterations=100, tolerance=1e-7):
    """Solves the nonlinear system defined by the `calculate_system` callback.

    The array with the initial solution is updated during the solve and contains
    the solution at convergence.

    Parameters
    ----------
    calculate_system : function
        This function is called several times to evaluate the function (rhs) and the
        functions derivatives (lhs) with a given state (x)

        It should look like this:

        def calculate_system(x):
            ...
            return lhs, rhs
    x_initial : ndarray
        Initial guess of the solution
    max_iterations : int
        Maximum number of iterations
    tolerance : float
        Convergence tolerance value for the residual norm

    Raises
    ----------
    RuntimeError
        If the algorithm does not converge within `max_iterations`
    """
    x = x_initial
    residual_norm = None

    for i in range(1, max_iterations + 1):


        # calculate left and right hand side
        lhs, rhs = calculate_system(x)

        # calculate residual
        residual_norm = la.norm(rhs)

        # check convergence
        if residual_norm < tolerance:
            print('  Newthon-Raphson converged in step {}.'.format(i))
            print('  Residual norm: {}.'.format(residual_norm))
            return x, i

        # compute delta_x
        delta_x = sla.spsolve(lhs, rhs)

        # update x
        x -= delta_x

    raise RuntimeError('Newthon-Raphson did not converge after {} steps. Residual norm: {}'
                       .format(max_iterations, residual_norm))
