from KratosMultiphysics.CoSimulationApplication.co_simulation_component import CoSimulationComponent
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return PredictorLinear(parameters)


# Class PredictorLinear: Linear extrapolation based on the last two time steps, assuming constant time step size.
class PredictorLinear(CoSimulationComponent):
    def __init__(self, _unused):
        super().__init__()

        self.data_prev = []
        self.data_last = []

    def Initialize(self, x):
        super().Initialize()

        self.data_last = x
        self.data_prev = x

    def Predict(self):
        return self.data_last*2.0 - self.data_prev

    def Update(self, x):
        self.data_prev = self.data_last
        self.data_last = x
