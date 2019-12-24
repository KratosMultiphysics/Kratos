from KratosMultiphysics.CoSimulationApplication.predictors.predictor import Predictor
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return PredictorLinear(parameters)


# Class PredictorLinear: Linear extrapolation based on the last two time steps, assuming constant time step size.
class PredictorLinear(Predictor):
    def __init__(self, _unused):
        super().__init__(_unused)

        self.order = 1

    def Predict(self, x):
        return self.Linear(x)
