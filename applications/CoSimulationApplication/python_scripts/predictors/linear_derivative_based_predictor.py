from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_base_predictor import CosimulationBasePredictor

# Other imports
import numpy as np
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

def Create(predictor_settings, solver):
    return LinearDerivativeBasedPredictor(predictor_settings, solver)

class LinearDerivativeBasedPredictor(CosimulationBasePredictor):
    def __init__(self, settings, solver):
        super(LinearDerivativeBasedPredictor, self).__init__(settings, solver)
        # TODO add comment why we do this
        self.data_array            = np.array([])
        self.derivative_data_array = np.array([])

        # TODO check buffer size!

    def Predict(self):
        data_name = self.settings["data_name"].GetString()
        deriv_data_name = self.settings["derivative_data_name"].GetString()
        cs_tools.ImportArrayFromSolver(self.solver, deriv_data_name, self.derivative_data_array, 1)
        cs_tools.ImportArrayFromSolver(self.solver, data_name, self.data_array, 1)

        # TODO get the delta-time from the ModelPart-ProcessInfo in the future
        # => this will be accessible through the data-object (=> has the model)
        self.data_array += self.solver.GetDeltaTime() * self.derivative_data_array

        self._UpdateData(self.data_array)

    def _Name(self):
        return self.__class__.__name__

