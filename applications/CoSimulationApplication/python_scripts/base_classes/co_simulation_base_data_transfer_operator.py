from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Other imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

class CoSimulationBaseDataTransferOperator(object):
    def __init__(self, settings):
        self.settings = settings
        self.settings.ValidateAndAssignDefaults(self._GetDefaultSettings())
        self.echo_level = self.settings["echo_level"].GetInt()

    def TransferData(self, solver, data):
        pass

    def PrintInfo(self):
        pass

    def SetEchoLevel(self, level):
        self.echo_level = level

    def _Name(self):
        return self.__class__.__name__

    def Check(self):
        pass

    @classmethod
    def _GetDefaultSettings(cls):
        return cs_tools.cs_data_structure.Parameters("""{
            "echo_level" : 0
        }""")
