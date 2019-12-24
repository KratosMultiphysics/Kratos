from KratosMultiphysics.CoSimulationApplication.predictors.predictor import Predictor
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return PredictorQuadratic(parameters)


# Class PredictorQuadratic: Quadratic extrapolation based on the last three time steps, assuming constant time step size.
class PredictorQuadratic(Predictor):
    def __init__(self, _unused):
        super().__init__(_unused)

        self.order = 2

    def Predict(self, x):
        if len(self.dataprev) < 3:
            return self.Linear(x)
        else:
            return self.Quadratic(x)
