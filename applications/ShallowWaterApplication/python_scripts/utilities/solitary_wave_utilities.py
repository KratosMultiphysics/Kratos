from numpy import roots, cosh, errstate
from math import pi, sqrt, tanh

@errstate(over='ignore')
def sech(x):
    return 1/cosh(x)

class SolitaryWaveSolution:
    """Base class for analytical solutions of a solitary wave."""

    def __init__(self, depth, gravity=9.81, *, amplitude):
        self.depth = depth
        self.gravity = gravity
        self.amplitude = amplitude
        self.amplitude1 = amplitude
        self.amplitude2 = 0

    def eta(self, x, t):
        phase = self.wavenumber * (self.phase_speed * t - x)
        return self.amplitude1 * sech(phase)**2 + self.amplitude2 * sech(phase)**4

    def u(self, x, t):
        eta = self.eta(x, t)
        return self.phase_speed * eta / (self.depth + eta)

    @property
    def wavenumber(self):
        raise Exception("SolitaryWaveSolution. The wavenumber is not defined in the base class.")

    @property
    def phase_speed(self):
        return sqrt(self.gravity * (self.amplitude + self.depth))

    @property
    def frequency(self):
        return self.wavenumber * self.phase_speed

    @property
    def wavelength(self):
        return 2 * pi / self.wavenumber

    @property
    def period(self):
        return 2 * pi / self.frequency


class GoringSolution(SolitaryWaveSolution):
    """Goring analytical solution.

    K. Guizien and E. Barthelemy, Accuracy of solitary wave generation by a piston wave maker.
    Journal of Hydraulic Research, February 2010
    """

    @property
    def wavenumber(self):
        return sqrt(3 / 4 * self.amplitude / self.depth**3)


class RayleighSolution(SolitaryWaveSolution):
    """Goring analytical solution.

    K. Guizien and E. Barthelemy, Accuracy of solitary wave generation by a piston wave maker.
    Journal of Hydraulic Research, February 2010
    """

    @property
    def wavenumber(self):
        return sqrt(3 * self.amplitude / 4 / self.depth**2 / (self.depth + self.amplitude))


class BoussinesqSolution(SolitaryWaveSolution):
    """Analytical solution for a solitary wave with the modified Boussinesq equations.

    G. Wei and J. T. Kirby, Time-dependent numerical Code for extended Boussinesq equations.
    Journal of Waterway, Port, Coastal and Ocean Engineering, September 1995

    O. Nwogu, Alternative form of Boussinesq for nearshore wave propagation.
    Journal of Waterway, Port, Coastal and Ocean Engineering, 1993
    """

    def __init__(self, depth, gravity=9.81, *, amplitude):
        super().__init__(depth, gravity, amplitude=amplitude)
        self.beta = -0.531
        self.alpha = 0.5 * self.beta**2 + self.beta
        self.delta = self.amplitude / self.depth

        coefficients = [
            2*self.alpha,                                     # C^6
            -3*self.alpha - 1 / 3 - 2*self.alpha*self.delta,  # C^4
            2*self.delta*(self.alpha + 1 / 3),                # C^2
            self.alpha + 1 / 3]                               # 1
        c_roots = roots(coefficients)
        self.c_dimless = sqrt(c_roots[0])

        gh = self.gravity * self.depth
        c2 = gh * self.c_dimless**2
        gha3 = gh * (self.alpha + 1 / 3)
        self.amplitude1 = (c2 - gh) / 3 / (gha3 - self.alpha * c2) * self.depth
        self.amplitude2 = -(c2 - gh)**2 / 2 / gh / c2 * (gha3 + 2 * self.alpha * c2) / (gha3 - self.alpha * c2) * self.depth

    @property
    def phase_speed(self):
        return sqrt(self.gravity * self.depth) * self.c_dimless

    @property
    def wavenumber(self):
        gh = self.gravity * self.depth
        c2 = self.phase_speed**2
        gha3 = gh * (self.alpha + 1 / 3)
        return sqrt((c2 - gh) / 4 / (gha3 - self.alpha * c2)) / self.depth

    def u(self, x, t):
        gh = self.gravity * self.depth
        c2 = self.phase_speed**2
        horizontal_velocity = (c2 - gh) / self.phase_speed
        phase = self.wavenumber * (self.phase_speed * t - x)
        return horizontal_velocity * sech(phase)**2

    def a(self, x, t):
        gh = self.gravity * self.depth
        c2 = self.phase_speed**2
        horizontal_velocity = (c2 - gh) / self.phase_speed
        phase = self.wavenumber * (self.phase_speed * t - x)
        return -2 * self.frequency * horizontal_velocity * tanh(phase) * sech(phase)**2
