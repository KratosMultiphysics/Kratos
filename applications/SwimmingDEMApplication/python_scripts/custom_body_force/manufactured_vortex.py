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
        self.nu = settings["viscosity"].GetDouble()
        self.rho = settings["density"].GetDouble()
        self.T = self.L / self.U
        self.omega = settings["time_factor"].GetDouble() * 2 * np.pi / self.T
        self.k = np.pi / self.L

    def amplitude(self, t):
        return self.U * np.cos(self.omega * t)

    def ustatic1(self, x1, x2):
        return np.sin(2 * self.k * x2) * np.sin(self.k * x1)

    def ustatic2(self, x1, x2):
        return - np.sin(2 * self.k * x1) * np.sin(self.k * x2)

    def u1(self, x1, x2, t):
        return self.amplitude(t) * self.ustatic1(x1, x2)

    def u2(self, x1, x2, t):
        return self.amplitude(t) * self.ustatic2(x1, x2)

    def damplitudedt(self, t):
        return self.omega * self.U * np.cos(self.omega * t)

    def dustatic11(self, x1, x2):
        return self.k * np.sin(2 * self.k * x2) * np.cos(self.k * x1)

    def dustatic12(self, x1, x2):
        return 2 * self.k * np.cos(2 * self.k * x2) * np.sin(self.k * x1)

    def dustatic21(self, x1, x2):
        return - 2 * self.k * np.sin(2 * self.k * x1) * np.cos(self.k * x2)

    def dustatic22(self, x1, x2):
        return - self.k * np.cos(2 * self.k * x1) * np.sin(self.k * x2)

    def du11(self, x1, x2, t):
        return self.amplitude(t) * self.dustatic11(x1, x2)

    def du12(self, x1, x2, t):
        return self.amplitude(t) * self.dustatic12(x1, x2)

    def du21(self, x1, x2, t):
        return self.amplitude(t) *self. dustatic21(x1, x2)

    def du22(self, x1, x2, t):
        return self.amplitude(t) * self.dustatic22(x1, x2)

    def dustatic111(self, x1, x2):
        return - self.k**2 * np.sin(2 * self.k * x2) * np.sin(self.k * x1)

    def dustatic122(self, x1, x2):
        return - 4 * self.k**2 * np.sin(2 * self.k * x2) * np.sin(self.k * x1)

    def dustatic211(self, x1, x2):
        return - 4 * self.k**2 * np.cos(2 * self.k * x2) * np.cos(self.k * x1)

    def dustatic222(self, x1, x2):
        return - self.k**2 * np.cos(2 * self.k * x1) * np.cos(self.k * x2)

    def du111(self, x1, x2, t):
        return self.amplitude(t) * self.dustatic111(x1, x2)

    def du122(self, x1, x2, t):
        return self.amplitude(t) * self.dustatic122(x1, x2)

    def du211(self, x1, x2, t):
        return self.amplitude(t) * self.dustatic211(x1, x2)

    def du222(self, x1, x2, t):
        return self.amplitude(t) * self.dustatic222(x1, x2)

    def du1dt(self, x1, x2, t):
        return self.damplitudedt(t) * self.ustatic1(x1, x2)

    def du2dt(self, x1, x2, t):
        return self.damplitudedt(t) * self.ustatic2(x1, x2)
