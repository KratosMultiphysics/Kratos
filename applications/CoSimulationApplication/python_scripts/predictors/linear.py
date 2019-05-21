from KratosMultiphysics.CoSimulationApplication.co_simulation_component import CoSimulationComponent
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return PredictorLinear(parameters)


class PredictorLinear(CoSimulationComponent):
    def __init__(self, _unused):
        super().__init__()
        self.data_prev_iter = []
        self.data_current_iter = []
