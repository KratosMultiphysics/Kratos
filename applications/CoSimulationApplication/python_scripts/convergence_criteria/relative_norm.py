import numpy as np

from KratosMultiphysics.CoSimulationApplication.co_simulation_component import CoSimulationComponent
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return ConvergenceCriterionRelativeNorm(parameters)


class ConvergenceCriterionRelativeNorm(CoSimulationComponent):
    def __init__(self, parameters):
        super().__init__()

        settings = parameters["settings"]
        self.tolerance = settings["tolerance"].GetDouble()
        self.order = settings["order"].GetInt()

        self.initial_norm = 0
        self.last_norm = 0
        self.first_iteration = True

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        self.initial_norm = 0
        self.last_norm = 0
        self.first_iteration = True

    def Update(self, r):
        if self.first_iteration:
            self.initial_norm = np.linalg.norm(r, self.order)
        else:
            self.last_norm = np.linalg.norm(r, self.order)

    def IsSatisfied(self):
        if self.first_iteration or self.initial_norm < np.finfo(type(self.initial_norm)).eps:
            return False
        else:
            return self.last_norm/self.initial_norm < self.tolerance
