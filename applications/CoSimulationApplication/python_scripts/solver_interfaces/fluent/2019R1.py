from KratosMultiphysics.CoSimulationApplication.co_simulation_component import CoSimulationComponent
from KratosMultiphysics.CoSimulationApplication.co_simulation_interface import CoSimulationInterface
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return SolverInterfaceFluent2019R1(parameters)


class SolverInterfaceFluent2019R1(CoSimulationComponent):
    def __init__(self, parameters):
        super().__init__()

        self.settings = parameters["settings"]
