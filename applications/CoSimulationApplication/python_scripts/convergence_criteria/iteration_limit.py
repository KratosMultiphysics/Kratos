from KratosMultiphysics.CoSimulationApplication.co_simulation_component import CoSimulationComponent
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return ConvergenceCriterionIterationLimit(parameters)


class ConvergenceCriterionIterationLimit(CoSimulationComponent):
    def __init__(self, parameters):
        super().__init__()

        settings = parameters["settings"]
        self.maximum = settings["maximum"].GetInt()

        self.iteration = 0

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        self.iteration = 0

    def Update(self, _unused):
        self.iteration += 1

    def IsSatisfied(self):
        cs_tools.PrintInfo("Iteration: " + str(self.iteration))
        return self.iteration >= self.maximum
