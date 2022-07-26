# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_predictor import CoSimulationPredictor

# Other imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

def Create(settings, solver_wrapper):
    cs_tools.SettingsTypeCheck(settings)
    return LinearDerivativeBasedPredictor(settings, solver_wrapper)

class LinearDerivativeBasedPredictor(CoSimulationPredictor):
    def __init__(self, settings, solver_wrapper):
        super().__init__(settings, solver_wrapper)
        self.interface_derivative_data = solver_wrapper.GetInterfaceData(self.settings["derivative_data_name"].GetString())

    def Predict(self):
        if not self.interface_data.IsDefinedOnThisRank(): return

        data  = self.interface_data.GetData(1)
        derivative_data  = self.interface_derivative_data.GetData(1)

        delta_time = self.interface_data.GetModelPart().ProcessInfo[KM.DELTA_TIME]
        if abs(delta_time) < 1E-15:
            raise Exception("delta-time is almost zero!")

        data += delta_time * derivative_data

        self._UpdateData(data)

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "derivative_data_name" : "UNSPECIFIED"
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults
