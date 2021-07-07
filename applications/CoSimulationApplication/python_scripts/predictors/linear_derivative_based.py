# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_predictor import CoSimulationPredictor, HistoricalDataAccessor

# Other imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

def Create(settings, solver_wrapper):
    cs_tools.SettingsTypeCheck(settings)
    return LinearDerivativeBasedPredictor(settings, solver_wrapper)

class LinearDerivativeBasedPredictor(CoSimulationPredictor):
    def __init__(self, settings, solver_wrapper):
        super().__init__(settings, solver_wrapper)
        self.historical_derivative_data_accessor = HistoricalDataAccessor(solver_wrapper.GetInterfaceData(self.settings["derivative_data_name"].GetString()), self._GetMinimumBufferSize())

    def Predict(self):
        data  = self.historical_data_accessor.GetData(1)
        derivative_data = self.historical_derivative_data_accessor.GetData(1)

        # TODO this is a hack and currently is only supported by Kratos solvers!
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
