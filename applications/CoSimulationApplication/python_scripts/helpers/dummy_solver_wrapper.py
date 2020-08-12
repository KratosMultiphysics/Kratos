# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper

# Other imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

def Create(settings, model, solver_name):
    return DummySolverWrapper(settings, model, solver_name)

class DummySolverWrapper(CoSimulationSolverWrapper):
    """This class serves as dummy for testing, it does not solve anything
    It only imports a mesh that can be used in the testing workflow

    Note that this is only an example, other configurations are of course also possible
    """
    def __init__(self, settings, model, solver_name):
        super().__init__(settings, model, solver_name)

        self.time_step = self.settings["solver_wrapper_settings"]["time_step"].GetDouble()
        self.model_part = self.model.CreateModelPart(self.settings["solver_wrapper_settings"]["main_model_part_name"].GetString())

        self.model_part.ProcessInfo[KM.DOMAIN_SIZE] = self.settings["solver_wrapper_settings"]["domain_size"].GetInt()

        cs_tools.AllocateHistoricalVariablesFromCouplingData(self.data_dict.values(), self.model, self.name)

    def Initialize(self):
        severity = KM.Logger.GetDefaultOutput().GetSeverity()
        KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING) # mute MP-IO
        model_part_io = KM.ModelPartIO(self.settings["solver_wrapper_settings"]["mdpa_file_name"].GetString())
        model_part_io.ReadModelPart(self.model_part)
        KM.Logger.GetDefaultOutput().SetSeverity(severity)

        super().Initialize()

    def AdvanceInTime(self, current_time):
        return current_time + self.time_step

    def PrintInfo(self):
        cs_tools.cs_print_info("DummySolver", self._ClassName())
        ## TODO print additional stuff with higher echo-level
