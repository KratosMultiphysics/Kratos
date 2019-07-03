import numpy as np

from KratosMultiphysics.CoSimulationApplication.co_simulation_component import CoSimulationComponent
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return ConvergenceCriterionAbsoluteNorm(parameters)


class ConvergenceCriterionAbsoluteNorm(CoSimulationComponent):
    def __init__(self, parameters):
        super().__init__()

        settings = parameters["settings"]
        self.tolerance = settings["tolerance"].GetDouble()
        self.order = settings["order"].GetInt()

        self.last_norm = 0.0
        self.is_updated = False

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        self.last_norm = 0.0
        self.is_updated = False

    def Update(self, r):
        self.last_norm = np.linalg.norm(r.GetNumpyArray(), self.order)
        self.is_updated = True

    def IsSatisfied(self):
        cs_tools.PrintInfo("Norm: " + str(self.last_norm))
        if not self.is_updated:
            return False
        else:
            return self.last_norm < self.tolerance
