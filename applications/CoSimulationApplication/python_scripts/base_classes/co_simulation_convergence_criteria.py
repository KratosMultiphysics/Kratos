# Importing the Kratos Library
import KratosMultiphysics as KM

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
import KratosMultiphysics.CoSimulationApplication.colors as colors

class CoSimulationConvergenceCriteria(object):
    """Baseclass for the convergence criteria used for CoSimulation
    Checks if convergence was achieved in a (strongly) coupled simulation
    """
    def __init__(self, settings):
        self.settings = settings
        self.settings.RecursivelyValidateAndAssignDefaults(self._GetDefaultSettings())

        if self.settings["criteria_composition"].GetString() == "energy_conjugate":
            if self.settings["conjugate_data_name"].GetString() == "UNSPECIFIED":
                self.__RaiseException('Energy conjugate criteria composition requires energy conjugate variables to be specified in "data_name" and "conjugate_data_name".')

        if "domain_difference" in self.settings["criteria_options"].GetStringArray():
            if self.settings["solver_domain_two"].GetString() == "UNSPECIFIED":
                self.__RaiseException('Domain difference requires "solver_domain_two" to be set to the second domain.')

        self.echo_level = self.settings["echo_level"].GetInt()
        self.ignore_first_convergence = self.settings["ignore_first_convergence"].GetBool()

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
    def _GetDefaultSettings(cls):
        return KM.Parameters("""{
            "type"       : "UNSPECIFIED",
            "solver"     : "UNSPECIFIED",
            "data_name"  : "UNSPECIFIED",
            "criteria_composition" : "primal",
            "criteria_options" : [],
            "conjugate_data_name" : "UNSPECIFIED",
            "solver_domain_two" : "UNSPECIFIED",
            "ignore_first_convergence" : false,
            "echo_level" : 0
        }""")
