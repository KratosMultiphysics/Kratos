from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_predictor import CosimulationPredictor

# Other imports
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import classprint

def Create(settings, solver):
    cs_tools.SettingsTypeCheck(settings)
    return StandardLinearPredictor(settings, solver)

class StandardLinearPredictor(CosimulationPredictor):
    def Predict(self):
        current_data  = self.interface_data.GetNumpyArray(0)
        previous_data  = self.interface_data.GetNumpyArray(1)

        predicted_data = 2*current_data - previous_data

        self._UpdateData(predicted_data)
