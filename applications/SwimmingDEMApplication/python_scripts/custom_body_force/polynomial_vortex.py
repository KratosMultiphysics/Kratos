import KratosMultiphysics
import numpy as np

## Import base class file
from custom_body_force.manufactured_solution import ManufacturedSolution

def CreateManufacturedSolution(custom_settings):
    return PolynomialVortex(custom_settings)

class PolynomialVortex(ManufacturedSolution):
    def __init__(self, settings):

        default_settings = KratosMultiphysics.Parameters("""
            {
                "velocity"    : 1.0,
                "length"      : 1.0,
                "viscosity"   : 1.0e-2,
                "density"     : 1.0,
                "frequency"   : 1.0,
                "damping"     : 1.0
            }
            """
            )

        settings.ValidateAndAssignDefaults(default_settings)

        self.a = 16*np.sqrt(3**3)
        self.U = settings["velocity"].GetDouble()
        self.sq_au = np.sqrt(self.a * self.U)
        self.L = settings["length"].GetDouble()
        self.rho = settings["density"].GetDouble()
        self.nu = settings["viscosity"].GetDouble() / self.rho
        self.T = self.L / self.U
        self.omega = settings["frequency"].GetDouble() * 2 * np.pi
        self.damp = settings["damping"].GetDouble()

    def f(self, x):
        return self.sq_au * x**2 * (self.L - x)**2

    def g(self, t):
        return np.cos(self.omega * t) * np.e**(-self.damp * t)

    def df(self, x):
        return 2 * self.sq_au * x * (self.L - x)**2 - 2 * self.sq_au * x**2 * (self.L - x)

    def dg(self, t):
        return -self.omega * np.sin(self.omega * t) * np.e**(-self.damp * t) -self.damp * np.cos(self.omega * t) * np.e**(-self.damp * t)

    def ddf(self, x):
        return 2 * self.sq_au * (self.L**2 - 6*self.L*x + 6*x**2)

    def dddf(self, x):
        return 12 * self.sq_au * (2*x - self.L)

    def u1(self, x1, x2, t):
        return self.f(x1) * self.df(x2) * self.g(t)

    def u2(self, x1, x2, t):
        return -self.df(x1) * self.f(x2) * self.g(t)

    def du1dt(self, x1, x2, t):
        return self.f(x1) * self.df(x2) * self.dg(t)

    def du2dt(self, x1, x2, t):
        return -self.df(x1) * self.f(x2) * self.dg(t)

    def du11(self, x1, x2, t):
        return self.df(x1) * self.df(x2) * self.g(t)

    def du12(self, x1, x2, t):
        return self.f(x1) * self.ddf(x2) * self.g(t)

    def du21(self, x1, x2, t):
        return -self.ddf(x1) * self.f(x2) * self.g(t)

    def du22(self, x1, x2, t):
        return -self.df(x1) * self.df(x2) * self.g(t)

    def du111(self, x1, x2, t):
        return self.ddf(x1) * self.df(x2) * self.g(t)

    def du122(self, x1, x2, t):
        return self.f(x1) * self.dddf(x2) * self.g(t)

    def du211(self, x1, x2, t):
        return -self.dddf(x1) * self.f(x2) * self.g(t)

    def du222(self, x1, x2, t):
        return -self.df(x1) * self.ddf(x2) * self.g(t)
