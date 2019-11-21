import numpy as np
from scipy.linalg import qr, solve_triangular

from KratosMultiphysics.CoSimulationApplication.co_simulation_component import CoSimulationComponent
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return ConvergenceAcceleratorAITKEN(parameters)


class ConvergenceAcceleratorAITKEN(CoSimulationComponent):
    def __init__(self, parameters):
        super().__init__()

        self.settings = parameters["settings"]
        self.omega_max = self.settings["omega_max"].GetDouble()

        self.omega = self.omega_max
        self.added = False
        self.rcurr = None

    def Predict(self, r_in):
        r = r_in.GetNumpyArray()
        # Calculate return value if sufficient data available
        if not self.added:
            raise RuntimeError("No information to predict")
        dx = self.omega * r
        dx_out = r_in.deepcopy()
        dx_out.SetNumpyArray(dx)
        return dx_out

    def Update(self, x_in, xt_in):
        x = x_in.GetNumpyArray().reshape(-1, 1)
        xt = xt_in.GetNumpyArray().reshape(-1, 1)
        r = xt - x
        rprev = self.rcurr
        self.rcurr = r
        if self.added:
            # Aitken Relaxation
            # Update omega
            self.omega *= -float(rprev.T @ (r - rprev) / np.linalg.norm(r - rprev, 2) ** 2)
        else:
            # Set first value of omega in a timestep
            self.omega = np.sign(self.omega) * min(abs(self.omega), self.omega_max)
            self.added = True

    def IsReady(self):
        return self.added

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        self.added = False
        self.rcurr = None
