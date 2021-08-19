# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_predictor import CoSimulationPredictor

# Other imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

def Create(settings, solver_wrapper):
    cs_tools.SettingsTypeCheck(settings)
    return LinearPredictor(settings, solver_wrapper)

class LinearPredictor(CoSimulationPredictor):
    def Predict(self):
        current_data  = self.interface_data.GetData(0)
        previous_data  = self.interface_data.GetData(1)

        predicted_data = 2*current_data - previous_data

        self._UpdateData(predicted_data)
