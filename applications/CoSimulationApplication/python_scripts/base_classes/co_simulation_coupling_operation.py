from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

class CoSimulationCouplingOperation(object):
    def __init__(self, settings):
        self.settings = settings
        self.settings.ValidateAndAssignDefaults(self._GetDefaultSettings())
        self.echo_level = self.settings["echo_level"].GetInt()

    def Initialize(self):
        pass

    def Finalize(self):
        pass


    def InitializeSolutionStep(self):
        pass

    def FinalizeSolutionStep(self):
        pass


    def InitializeCouplingIteration(self):
        pass

    def FinalizeCouplingIteration(self):
        pass


    def Execute(self):
        pass


    def PrintInfo(self):
        pass

    def _Name(self):
        return self.__class__.__name__

    def Check(self):
        pass

    @classmethod
    def _GetDefaultSettings(cls):
        return cs_tools.cs_data_structure.Parameters("""{
            "echo_level" : 0
        }""")
