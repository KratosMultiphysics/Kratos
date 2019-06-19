from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_predictor import CoSimulationPredictor

# Other imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

def Create(settings, solver_wrapper):
    cs_tools.SettingsTypeCheck(settings)
    return LinearDerivativeBasedPredictor(settings, solver_wrapper)

class LinearDerivativeBasedPredictor(CoSimulationPredictor):
    def __init__(self, settings, solver_wrapper):
        super(LinearDerivativeBasedPredictor, self).__init__(settings, solver_wrapper)
        self.interface_derivative_data = solver_wrapper.GetInterfaceData(self.settings["derivative_data_name"].GetString())
        self.solver_wrapper = solver_wrapper

    def Predict(self):
        data  = self.interface_data.GetNumpyArray(1)
        derivative_data  = self.interface_derivative_data.GetNumpyArray(1)

        delta_time = self.interface_data.GetModelPart().ProcessInfo[cs_tools.cs_data_structure.DELTA_TIME]

        data += delta_time * derivative_data

        self._UpdateData(data)

    @classmethod
    def _GetDefaultSettings(cls):
        this_defaults = cs_tools.cs_data_structure.Parameters("""{
            "derivative_data_name" : "UNSPECIFIED"
        }""")
        this_defaults.AddMissingParameters(super(LinearDerivativeBasedPredictor, cls)._GetDefaultSettings())
        return this_defaults
