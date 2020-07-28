# Importing the Kratos Library
import KratosMultiphysics as KM

# CoSimulation imports
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import SettingsTypeCheck

class CoSimulationCouplingOperation(object):
    """Baseclass for the coupling operations used for CoSimulation
    This class can be used to customize the behavior of the CoSimulation,
    by providing a large interface and access to the solvers/models
    """
    def __init__(self, settings, parent_coupled_solver_process_info):
        SettingsTypeCheck(settings)
        self.settings = settings
        self.settings.ValidateAndAssignDefaults(self._GetDefaultSettings())
        self.process_info = parent_coupled_solver_process_info
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
        raise NotImplementedError('"Execute" is not implemented for {}!'.format(self._ClassName))


    def PrintInfo(self):
        pass

    def Check(self):
        pass

    @classmethod
    def _ClassName(cls):
        return cls.__name__

    @classmethod
    def _GetDefaultSettings(cls):
        return KM.Parameters("""{
            "type"       : "UNSPECIFIED",
            "echo_level" : 0
        }""")
