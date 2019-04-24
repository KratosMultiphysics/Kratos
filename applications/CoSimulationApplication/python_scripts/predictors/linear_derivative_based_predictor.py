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

        self.interface_derivative_data = self.solver.GetInterfaceData(self.settings["derivative_data_name"].GetString())
        # TODO check buffer size!

    def Predict(self):
        data  = self.interface_data.GetNumpyArray(1)
        derivative_data  = self.interface_derivative_data.GetNumpyArray(1)

        # TODO get the delta-time from the ModelPart-ProcessInfo in the future
        # => this will be accessible through the data-object (=> has the model)
        data += self.solver.GetDeltaTime() * derivative_data

        self._UpdateData(data)

    def _Name(self):
        return self.__class__.__name__
