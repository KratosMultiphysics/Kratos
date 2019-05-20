from KratosMultiphysics.CoSimulationApplication.co_simulation_component import CoSimulationComponent
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return SolverInterfacePipeFlow(parameters)


class SolverInterfacePipeFlow(CoSimulationComponent):
    def __init__(self, parameters):
        self.parameters = parameters
        self.settings = parameters["settings"]
