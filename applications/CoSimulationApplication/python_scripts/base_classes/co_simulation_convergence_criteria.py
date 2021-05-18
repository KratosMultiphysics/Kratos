# Importing the Kratos Library
import KratosMultiphysics as KM

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
import KratosMultiphysics.CoSimulationApplication.colors as colors

class CoSimulationConvergenceCriteria:
    """Baseclass for the convergence criteria used for CoSimulation
    Checks if convergence was achieved in a (strongly) coupled simulation
    """
    def __init__(self, settings):
        self.settings = settings
        self.settings.RecursivelyValidateAndAssignDefaults(self._GetDefaultParameters())

        self.echo_level = self.settings["echo_level"].GetInt()

    def Initialize(self):
        pass

    def Finalize(self):
        pass

    def InitializeSolutionStep(self):
        pass

    def FinalizeSolutionStep(self):
        pass

    def InitializeNonLinearIteration(self):
        pass

    def FinalizeNonLinearIteration(self):
        pass

    def IsConverged(self, residual, current_data):
        raise NotImplementedError('"IsConverged" has to be implemented in the derived class!')

    def PrintInfo(self):
        cs_tools.cs_print_info("Convergence Criteria", colors.bold(self._ClassName()))

    def Check(self):
        cs_tools.cs_print_warning("Convergence Criteria", colors.bold(self._ClassName()), 'does not implement "Check"')

    @classmethod
    def _ClassName(cls):
        return cls.__name__

    @classmethod
    def _GetDefaultParameters(cls):
        return KM.Parameters("""{
            "type"       : "UNSPECIFIED",
            "solver"     : "UNSPECIFIED",
            "data_name"  : "UNSPECIFIED",
            "echo_level" : 0
        }""")
