from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_base_predictor import CosimulationBasePredictor

# Other imports
import numpy as np
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import classprint

def Create(predictor_settings, solver):
    return StandardLinearPredictor(predictor_settings, solver)

class StandardLinearPredictor(CosimulationBasePredictor):
    def __init__(self, settings, solver):
        super(StandardLinearPredictor, self).__init__(settings, solver)

        # TODO check buffer size!

    def Predict(self):
        current_data  = self.interface_data.GetNumpyArray(0)
        previous_data  = self.interface_data.GetNumpyArray(1)

        predicted_data = 2*current_data - previous_data

        self._UpdateData(predicted_data)

    def _Name(self):
        return self.__class__.__name__
