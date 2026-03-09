import numpy as np
from collections.abc import Iterable
from ode_solve import ode_solve
from RationalApproximation import compute_RationalApproximation_AAA
from RationalApproximation import compute_RationalApproximation_AAA_new
from time import time
import multiprocessing as mp
import os

class fde_solve:
    '''
    solves (1+\mathcal{L})^{\alpha}\mathcal{L}^{\beta}\psi = f in [0,H]
                                                      \psi = 0 at z=0
                                                  d\psi/dz = 0 at z=H
    '''
    def __init__(self, dof, alpha, coef, domain_height=1, t=0, Robin_const=None, beta=None, z_grid=None):
        self.dof = dof
        self.alpha = alpha
        self.beta = beta
        self.coef = coef
        self.domain_height = domain_height
        self.t = t
        self.Robin_const=Robin_const
        self.z_grid = z_grid
        self.reset_ode(coef)
        self.c, self.d = compute_RationalApproximation_AAA_new(alpha, beta, nPoints=1000)
        self.nModes = self.c.size
        self.grid = self.ode_solve.grid[:]
        self.h = self.ode_solve.grid.h
        if self.beta is None:
            self.kappa_alpha = self.coef(self.grid)**self.alpha
        else:
            self.kappa_alpha = self.coef(self.grid)**(self.alpha+self.beta)

    def reset_ode(self, coef):
        if isinstance(coef, Iterable):
            self.ode_solve = [ ode_solve(self.dof, coefj, domain_height=self.domain_height, grid=self.z_grid) for coefj in coef ]
            self.anisotrop = True
            self.apply_sqrtMass = self.ode_solve[0].apply_sqrtMass
            self.apply_Mass = self.ode_solve[0].apply_Mass
        else:
            self.ode_solve = ode_solve(self.dof, coef, domain_height=self.domain_height, grid=self.z_grid)
            self.anisotrop = False
            self.apply_sqrtMass = self.ode_solve.apply_sqrtMass
            self.apply_Mass = self.ode_solve.apply_Mass
        

    def reset_kappa(self, coef):
        try:
            det = 1
            for coefj in coef:
                det = det * coefj(self.grid)
            det = det**(1/len(coef))
        except:
            det = coef(self.grid)
        if self.beta is None:
            self.kappa_alpha = det**self.alpha
        else:
            self.kappa_alpha = det**(self.alpha+self.beta)

    def reset_parameters(self, coef=None, t=None, Robin_const=None):
        self.Robin_const=Robin_const
        if coef is not None:
            self.coef = coef
            self.reset_ode(coef)
            self.reset_kappa(coef)

    def reset_jac(self, grad_coef):
        self.grad_term1 = [ self.alpha * der_coef(self.grid)/self.coef(self.grid) for der_coef in grad_coef ]
        self.gradA = [ ode_solve(self.dof, der_coef, domain_height=self.domain_height) for der_coef in grad_coef ]

        self.e1 = np.zeros(self.dof)
        self.e1[0] = 1

    def __call__(self, f, k1, k2, **kwargs):
        self.t = kwargs.get('t', 0)
        self.Robin_const = kwargs.get('Robin_const', None)
        self.adjoint = kwargs.get('adjoint', False)
        jac = kwargs.get('jac', False)
        grad_coef = kwargs.get('grad_coef', False)
        self.component = kwargs.get('component', 0)
        self.kwargs = kwargs
        self.kwargs.__delitem__('Robin_const')

        if f is not None:
            if np.linalg.norm(f) == 0:
                return np.zeros(self.z_grid.size)

        self.rhs = 1*f  ### variant of copying vector (important for not modifying the input)
        self.k = (k1, k2)
        if not self.adjoint:
            self.rhs *= self.kappa_alpha

        args = list(enumerate((self,)*self.nModes))
        self.modes = np.array(list(map(func, args)))
        psi_n_approx = self.modes.T @ self.c #np.sum(*, axis=0)

        return psi_n_approx


def func(args):
    i, self = args
    if self.anisotrop:
        return self.ode_solve[self.component](  self.d[i], self.rhs, self.k[0], self.k[1], Robin_const=self.Robin_const, **self.kwargs)
    else:
        return self.ode_solve(  self.d[i], self.rhs, self.k[0], self.k[1], Robin_const=self.Robin_const, **self.kwargs)

def func_jac(args):
    i, self = args
    return self.ode_solve(  self.d[i], self.Bf_conj, self.k[0], self.k[1],
                            t=self.t, Robin_const=self.Robin_const, adjoint=False)  

def func_bc(args):
    i, self = args
    return -self.modes_jac[i][0]*self.modes[i][0]/self.h

def func_assemble(args):
    i, self = args
    gradA_mode = self.gradA[self.comp].apply_matvec(-1, self.modes[i], self.k[0], self.k[1])
    return np.sum(self.modes_jac[i]*gradA_mode)

if __name__ == "__main__":

    import matplotlib.pyplot as plt

    '''TEST by method of manufactured solutions'''

    ## TODO