# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_predictor import CoSimulationPredictor

# Other imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

def Create(settings, solver_wrapper):
    cs_tools.SettingsTypeCheck(settings)
    return LinearPredictor(settings, solver_wrapper)

class LinearPredictor(CoSimulationPredictor):
    def Predict(self):
        if not self.historical_data_accessor.TimeBufferIsInitialized():
            if self.echo_level > 0:
                cs_tools.cs_print_info(self._ClassName(), "Skipped prediction because time buffer is not yet filled")
            return

        current_data  = self.historical_data_accessor.GetData(0)
        previous_data  = self.historical_data_accessor.GetData(1)

        predicted_data = 2*current_data - previous_data

        self._UpdateData(predicted_data)
