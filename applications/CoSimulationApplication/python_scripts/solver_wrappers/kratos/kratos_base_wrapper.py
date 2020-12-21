# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper

# Other imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
from importlib import import_module


class ThreadManager:
    """Class for setting and ressting the number of threads a context should use."""
    def __init__(self, num_threads=None):
        self.num_threads = num_threads
        if self.num_threads:
            self.num_threads_orig = KM.OpenMPUtils().GetNumThreads()

    def __enter__(self):
        if self.num_threads:
            KM.OpenMPUtils().SetNumThreads(min(self.num_threads, self.num_threads_orig))

    def __exit__(self, exc_type, exc_value, traceback):
        if self.num_threads:
            KM.OpenMPUtils().SetNumThreads(self.num_threads_orig)


def Create(settings, model, solver_name):
    return KratosBaseWrapper(settings, model, solver_name)

class KratosBaseWrapper(CoSimulationSolverWrapper):
    """This class serves as basis for the kratos-wrappers
    It uses the AnalysisStage as black-box interface to Kratos
    """
    def __init__(self, settings, model, solver_name):
        super().__init__(settings, model, solver_name)

        input_file_name = self.settings["solver_wrapper_settings"]["input_file"].GetString()
        if not input_file_name.endswith(".json"):
            input_file_name += ".json"

        with open(input_file_name,'r') as parameter_file:
            self.project_parameters = KM.Parameters(parameter_file.read())

        # this creates the AnalysisStage, creates the MainModelParts and allocates the historial Variables on the MainModelParts:
        self._analysis_stage = self.__GetAnalysisStage()

        if self.settings["solver_wrapper_settings"].Has("num_threads"):
            omp_num_threads = self.settings["solver_wrapper_settings"]["num_threads"].GetInt()
            self.thread_manager = ThreadManager(omp_num_threads)
        else:
            self.thread_manager = ThreadManager()


    def Initialize(self):
        with self.thread_manager:
            self._analysis_stage.Initialize() # this reades the Meshes
        super().Initialize()

    def Finalize(self):
        super().Finalize()
        with self.thread_manager:
            self._analysis_stage.Finalize()

    def AdvanceInTime(self, current_time):
        with self.thread_manager:
            new_time = self._analysis_stage._GetSolver().AdvanceInTime(current_time)
        self._analysis_stage.time = new_time # only needed to print the time correctly
        return new_time

    def InitializeSolutionStep(self):
        with self.thread_manager:
            self._analysis_stage.InitializeSolutionStep()

    def Predict(self):
        with self.thread_manager:
            self._analysis_stage._GetSolver().Predict()

    def SolveSolutionStep(self):
        with self.thread_manager:
            self._analysis_stage._GetSolver().SolveSolutionStep()
        super().SolveSolutionStep()

    def FinalizeSolutionStep(self):
        with self.thread_manager:
            self._analysis_stage.FinalizeSolutionStep()

    def OutputSolutionStep(self):
        with self.thread_manager:
            self._analysis_stage.OutputSolutionStep()

    def _CreateAnalysisStage(self):
        raise Exception('The "KratosBaseWrapper" can only be used when specifying "analysis_stage_module", otherwise the creation of the AnalysisStage must be implemented in the derived class!')

    def __GetAnalysisStage(self):
        if self.settings["solver_wrapper_settings"].Has("analysis_stage_module"):
            analysis_stage_module = import_module(self.settings["solver_wrapper_settings"]["analysis_stage_module"].GetString())
            return analysis_stage_module.Create(self.model, self.project_parameters)
        else:
            return self._CreateAnalysisStage()

    def PrintInfo(self):
        cs_tools.cs_print_info("KratosSolver", self._ClassName())
        cs_tools.cs_print_info("KratosSolver", 'Using AnalysisStage "{}", defined in module "{}'.format(self._analysis_stage.__class__.__name__, self._analysis_stage.__class__.__module__))
        ## TODO print additional stuff with higher echo-level
