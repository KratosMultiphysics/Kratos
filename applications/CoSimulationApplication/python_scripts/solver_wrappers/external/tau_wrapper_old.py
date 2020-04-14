from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper

# Other imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

def Create(settings, solver_name):
    return TAUWrapper(settings, solver_name)

class TAUWrapper(CoSimulationSolverWrapper):
    """This class serves as wrapper for the CFD solver TAU
    """
    def __init__(self, settings, solver_name):
        super(TAUWrapper, self).__init__(settings, solver_name)

        # settings_defaults = KM.Parameters("""{
        #     "model_parts_read" : { },
        #     "model_parts_send" : { },
        #     "model_parts_recv" : { }
        # }""")

        # self.settings["solver_wrapper_settings"].ValidateAndAssignDefaults(settings_defaults)

        # cs_tools.CreateMainModelPartsFromCouplingData(self.data_dict.values(), self.model, self.name)
        # cs_tools.AllocateHistoricalVariablesFromCouplingData(self.data_dict.values(), self.model, self.name)

    def Initialize(self):
        super(TAUWrapper, self).Initialize()
        to_io_data_config = {
            "type" : "control_signal",
            "signal" : "Initialize"
        }
        self.ExportData(to_io_data_config)

    def SolveSolutionStep(self):
        super(TAUWrapper, self).SolveSolutionStep()
        to_io_data_config = {
            "type" : "control_signal",
            "signal" : "SolveSolutionStep"
        }
        self.ExportData(to_io_data_config)

    def Finalize(self):
        super(TAUWrapper, self).Finalize()
        to_io_data_config = {
            "type" : "control_signal",
            "signal" : "Finalize"
        }
        self.ExportData(to_io_data_config)

    def AdvanceInTime(self, current_time):
        return 100.0 # TODO find a better solution here... maybe get time from solver through IO

    def PrintInfo(self):
        cs_tools.cs_print_info(self._ClassName(), "printing info...")
        ## TODO print additional stuff with higher echo-level

    def _GetIOType(self):
        return self.settings["io_settings"]["type"].GetString()