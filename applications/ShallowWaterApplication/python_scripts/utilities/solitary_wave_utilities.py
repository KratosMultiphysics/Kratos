import numpy as np
from mpmath import sech


class SolitaryWaveSolution:
    """Base class for analytical solutions of a solitary wave."""

    def __init__(self, depth, amplitude):
        self.depth = depth
        self.amplitude = amplitude
        self.gravity = 9.81
        self.amplitude1 = self.amplitude
        self.amplitude2 = 0

    def eta(self, x, t):
        phase = self.k * (self.c * t - x)
        return self.amplitude1 * sech(phase)**2 + self.amplitude2 * sech(phase)**4

    @property
    def k(self):
        raise Exception("SolitaryWaveSolution. The wavenumber 'k' is not defined in the base class.")

    @property
    def c(self):
        return np.sqrt(self.gravity * (self.amplitude + self.depth))


class GoringSolution(SolitaryWaveSolution):
    """Goring analytical solution.

    K. Guizien and E. Barthelemy, Accuracy of solitary wave generation by a piston wave maker.
    Journal of Hydraulic Research, February 2010
    """

    @property
    def k(self):
        return np.sqrt(3 / 4 * self.amplitude / self.depth**3)


class RayleighSolution(SolitaryWaveSolution):
    """Goring analytical solution.

    K. Guizien and E. Barthelemy, Accuracy of solitary wave generation by a piston wave maker.
    Journal of Hydraulic Research, February 2010
    """

    @property
    def k(self):
        return np.sqrt(3 * self.amplitude / 4 / self.depth**2 / (self.depth + self.amplitude))


class WeiSolution(SolitaryWaveSolution):
    """Analytical solution for a solitary wave with the modified Boussinesq equations.

    G. Wei and J. T. Kirby, Time-dependent numerical Code for extended Boussinesq equations.
    Journal of Waterway, Port, Coastal and Ocean Engineering, September 1995

    O. Nwogu, Alternative form of Boussinesq for nearshore wave propagation.
    Journal of Waterway, Port, Coastal and Ocean Engineering, 1993
    """

    def __init__(self, depth, amplitude):
        super().__init__(depth, amplitude)
        self.beta = -0.531
        self.alpha = 0.5 * self.beta**2 + self.beta
        self.delta = self.amplitude / self.depth

        coefficients = [
            2*self.alpha,                                     # C^6
            -3*self.alpha - 1 / 3 - 2*self.alpha*self.delta,  # C^4
            2*self.delta*(self.alpha + 1 / 3),                # C^2
            self.alpha + 1 / 3]                               # 1
        roots = np.roots(coefficients)
        self.c_dimless = np.sqrt(roots[0])

        gh = self.gravity * self.depth
        c2 = gh * self.c_dimless**2
        gha3 = gh * (self.alpha + 1 / 3)
        self.amplitude1 = (c2 - gh) / 3 / (gha3 - self.alpha * c2) * self.depth
        self.amplitude2 = -(c2 - gh)**2 / 2 / gh / c2 * (gha3 + 2 * self.alpha * c2) / (gha3 - self.alpha * c2) * self.depth

    @property
    def c(self):
        return np.sqrt(self.gravity * self.depth) * self.c_dimless

    @property
    def k(self):
        gh = self.gravity * self.depth
        c2 = self.c**2
        gha3 = gh * (self.alpha + 1 / 3)
        return np.sqrt((c2 - gh) / 4 / (gha3 - self.alpha * c2)) / self.depth
