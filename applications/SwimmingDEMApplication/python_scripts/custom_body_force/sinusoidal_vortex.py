import KratosMultiphysics
import numpy as np

## Import base class file
from custom_body_force.manufactured_solution import ManufacturedSolution

def CreateManufacturedSolution(custom_settings):
    return ManufacturedVortex(custom_settings)

class ManufacturedVortex(ManufacturedSolution):
    def __init__(self, settings):

        default_settings = KratosMultiphysics.Parameters("""
            {
                "velocity"    : 1.0,
                "length"      : 1.0,
                "viscosity"   : 1.0e-2,
                "density"     : 1.0,
                "time_factor" : 1.0
            }
            """
            )

        settings.ValidateAndAssignDefaults(default_settings)

        self.U = settings["velocity"].GetDouble()
        self.L = settings["length"].GetDouble()
        self.rho = settings["density"].GetDouble()
        self.nu = settings["viscosity"].GetDouble() / self.rho
        self.T = self.L / self.U
        self.omega = settings["time_factor"].GetDouble() * 2 * np.pi / self.T
        self.k = np.pi / self.L

    # Auxiliary functions

    def a(self, t):
        return np.cos(self.omega * t) * np.e**(-t)
        # return self.U * np.sin(self.omega * t)

    def f(self, x):
        return np.sin(self.k * x)

    def g(self, x):
        return np.sin(2 * self.k * x)

    # First derivatives

    def da(self, t):
        return -self.omega * np.sin(self.omega * t) * np.e**(-t) - np.cos(self.omega * t) * np.e**(-t)
        # return self.omega * self.U * np.cos(self.omega * t)

    def df(self, x):
        return self.k * np.cos(self.k * x)

    def dg(self, x):
        return 2 * self.k * np.cos(2 * self.k * x)

    # Second derivatives

    def ddf(self, x):
        return - self.k**2 * np.sin(self.k * x)

    def ddg(self, x):
        return -4 * self.k**2 * np.sin(2 * self.k * x)

    # Velocity

    def u1(self, x1, x2, t):
        return self.a(t) * self.f(x1) * self.g(x2)

    def u2(self, x1, x2, t):
        return -self.a(t) * self.g(x1) * self.f(x2)

    # Velocity derivatives

    def du11(self, x1, x2, t):
        return self.a(t) * self.df(x1) * self.g(x2)

    def du12(self, x1, x2, t):
        return self.a(t) * self.f(x1) * self.dg(x2)

    def du21(self, x1, x2, t):
        return -self.a(t) * self.dg(x1) * self.f(x2)

    def du22(self, x1, x2, t):
        return -self.a(t) * self.g(x1) * self.df(x2)

    # Velocity second derivatives

    def du111(self, x1, x2, t):
        return self.a(t) * self.ddf(x1) * self.g(x2)

    def du122(self, x1, x2, t):
        return self.a(t) * self.f(x1) * self.ddg(x2)

    def du211(self, x1, x2, t):
        return -self.a(t) * self.ddg(x1) * self.f(x2)

    def du222(self, x1, x2, t):
        return -self.a(t) * self.g(x1) * self.ddf(x2)

    # Accelerations

    def du1dt(self, x1, x2, t):
        return self.da(t) * self.f(x1) * self.g(x2)

    def du2dt(self, x1, x2, t):
        return -self.da(t) * self.g(x1) * self.f(x2)
