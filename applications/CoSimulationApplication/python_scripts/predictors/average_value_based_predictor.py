from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_base_predictor import CosimulationBasePredictor

# Other imports
import numpy as np
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

# Predictor implemented according to:
# "A new staggered scheme for fluid-structure interaction"; W.G. Dettmer and D. Peric
# Numerical Methods in Engineering 2013; 93; 1-22

def Create(predictor_settings, solvers):
    return AverageValuePredictor(predictor_settings, solvers)

class AverageValuePredictor(CosimulationBasePredictor):
    # @param beta factor for weighting last and current value of the predicted values. Can be set in interval: [0, 1.0]
    def __init__(self, settings, solvers):
        super(AverageValuePredictor, self).__init__(settings, solvers)
        if "beta" in self.settings:
            self.beta = self.settings["beta"]
            if self.beta > 1.0 or self.beta < 0:
                raise Exception("Wrong value for beta. Admissible interval [0.0, 1.0]")
        else:
            self.beta = 0.5
        # TODO check buffer_size

        # TODO add comment why we do this
        self.data_array_prediction = np.array([])
        self.data_array_aux        = np.array([])

    def Predict(self):
        solver = self.solvers[self.settings["solver"]]
        data_name = self.settings["data_name"]
        cs_tools.ImportArrayFromSolver(solver, data_name, self.data_array_prediction, 0)
        cs_tools.ImportArrayFromSolver(solver, data_name, self.data_array_aux, 1)

        self.data_array_prediction = 2*self.data_array_prediction - self.data_array_aux

        self._UpdateData(self.data_array_prediction)

    def FinalizeSolutionStep(self):
        solver = self.solvers[self.settings["solver"]]
        data_name = self.settings["data_name"]
        cs_tools.ImportArrayFromSolver(solver, data_name, self.data_array_aux, 0)

        self.data_array_prediction = self.beta * self.data_array_aux + (1-self.beta) * self.data_array_prediction

        self._UpdateData(self.data_array_prediction)


    def _Name(self):
        return self.__class__.__name__

