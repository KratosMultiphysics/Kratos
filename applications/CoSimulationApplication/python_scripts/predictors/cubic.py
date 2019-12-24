from KratosMultiphysics.CoSimulationApplication.predictors.predictor import Predictor
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return PredictorCubic(parameters)


# Class PredictorCubic: Cubic extrapolation based on the last four time steps, assuming constant time step size.
class PredictorCubic(Predictor):
    def __init__(self, _unused):
        super().__init__(_unused)

        self.order = 3

    def Predict(self, x):
        if len(self.dataprev) < 3:
            return self.Linear(x)
        elif len(self.dataprev) == 3:
            return self.Quadratic(x)
        else:
            return self.Cubic(x)
