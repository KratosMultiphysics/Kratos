from KratosMultiphysics.CoSimulationApplication.co_simulation_component import CoSimulationComponent
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return PredictorLinear(parameters)


# Class PredictorLinear: Linear extrapolation based on the last two time steps, assuming constant time step size.
class PredictorLinear(CoSimulationComponent):
    def __init__(self, _unused):
        super().__init__()

        self.updated = False
        self.data_prev = []
        self.data_last = []

    def Initialize(self, x):
        super().Initialize()

        self.data_last = x.GetNumpyArray()
        self.data_prev = self.data_last

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        self.updated = False

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        if not self.updated:
            raise Exception("Not updated")

    def Predict(self, x):
        if not self.updated:
            y = self.data_last * 2.0 - self.data_prev
            x.SetNumpyArray(y)
            return x
        else:
            raise Exception("Already updated")

    def Update(self, x):
        if not self.updated:
            self.data_prev = self.data_last
            self.data_last = x.GetNumpyArray()
            self.updated = True
        else:
            raise Exception("Already updated")
