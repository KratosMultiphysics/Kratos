from KratosMultiphysics.CoSimulationApplication.co_simulation_component import CoSimulationComponent
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return ConvergenceAcceleratorIQNI(parameters)


class ConvergenceAcceleratorIQNI(CoSimulationComponent):
    def __init__(self, parameters):
        super().__init__()
        self.parameters = parameters

    def Predict(self, r):
        return 1

    def Update(self, x, xt):
        pass
