# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper

# Other imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

# Import importlib to be able to load analysis stages from a string
from importlib import import_module

# Import string for capwords function
import string

class ThreadManager:
    """Class for setting and ressting the number of threads a context should use."""
    def __init__(self, num_threads=None):
        self.num_threads = num_threads
        if self.num_threads:
            self.num_threads_orig = KM.ParallelUtilities.GetNumThreads()

    def __enter__(self):
        if self.num_threads:
            KM.ParallelUtilities.SetNumThreads(min(self.num_threads, self.num_threads_orig))

    def __exit__(self, exc_type, exc_value, traceback):
        if self.num_threads:
            KM.ParallelUtilities.SetNumThreads(self.num_threads_orig)


def Create(settings, model, solver_name):
    return KratosBaseWrapper(settings, model, solver_name)

class KratosBaseWrapper(CoSimulationSolverWrapper):
    """This class serves as basis for the kratos-wrappers
    It uses the AnalysisStage as black-box interface to Kratos
    """
    def __init__(self, settings, model, solver_name):
        # We try to read the input file
        if settings["solver_wrapper_settings"].Has("input_file"):
            input_file_name = settings["solver_wrapper_settings"]["input_file"].GetString()
            if not input_file_name.endswith(".json"):
                input_file_name += ".json"

            with open(input_file_name,'r') as parameter_file:
                self.project_parameters = KM.Parameters(parameter_file.read())
        else: # The settings are in the root Parameters
            self.project_parameters = settings["solver_wrapper_settings"]

        super().__init__(settings, model, solver_name)

        if self.settings["solver_wrapper_settings"].Has("num_threads"):
            omp_num_threads = self.settings["solver_wrapper_settings"]["num_threads"].GetInt()
            self.thread_manager = ThreadManager(omp_num_threads)
        else:
            self.thread_manager = ThreadManager()

        # this creates the AnalysisStage, creates the MainModelParts and allocates the historial Variables on the MainModelParts:
        with self.thread_manager:
            self._analysis_stage = self.__GetAnalysisStage()

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
            new_time = self._analysis_stage._AdvanceTime(current_time)
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
            module_name = self.settings["solver_wrapper_settings"]["analysis_stage_module"].GetString()
            analysis_stage_module = import_module(module_name)
            if hasattr(analysis_stage_module, "Create"):
                return analysis_stage_module.Create(self.model, self.project_parameters)
            else:
                KM.Logger.PrintWarning("KratosBaseWrapper", f'The analysis_stage_module "{module_name}" does not have a "Create" function, trying to create the AnalysisStage directly...')
                # We assume that the name of the AnalysisStage is the same as the name of the module in PascalCase instead of snake_case
                file_name = module_name.split(".")[-1]
                # Convert Snake case to Pascal case
                analysis_stage_name = string.capwords(file_name.replace("_", " ")).replace(" ", "")

                # Getting the analysis class
                if hasattr(analysis_stage_module, analysis_stage_name):
                    analysis = getattr(analysis_stage_module, analysis_stage_name)
                else:
                    KM.Logger.PrintWarning("KratosBaseWrapper", f'The analysis_stage_module "{module_name}" does not follow the standard way to define the analysis stage name "{analysis_stage_name}" . Trying to retrieve from a custom definition')
                    if self.settings["solver_wrapper_settings"].Has("analysis_name"):
                        analysis_stage_name = self.settings["solver_wrapper_settings"]["analysis_stage_name"].GetString()
                        if hasattr(analysis_stage_module, analysis_stage_name):
                            analysis = getattr(analysis_stage_module, analysis_stage_name)
                        else:
                            raise Exception(f'"{module_name}" does not have a "{analysis_stage_name}" class!')
                    else:
                        raise Exception(f'"{module_name}" does not have a "{analysis_stage_name}" class! Please provide a custom "analysis_stage_name" in your settings')

                return analysis(self.model, self.project_parameters)
        else:
            return self._CreateAnalysisStage()

    def PrintInfo(self):
        super().PrintInfo()
        cs_tools.cs_print_info("KratosSolver", self._ClassName())
        cs_tools.cs_print_info("KratosSolver", 'Using AnalysisStage "{}", defined in module "{}'.format(self._analysis_stage.__class__.__name__, self._analysis_stage.__class__.__module__))

    def _CheckDataCommunicatorIsConsistentlyDefined(self, import_settings, mpi_settings):
        """
        Checking if the data-comm used for the solver (specified in the import-settings,
        see "distributed_import_model_part_utility") is consistent with the one that should be
        created by the solver-wrapper
        """
        solver_uses_custom_data_comm = import_settings.Has("data_communicator_name")
        creating_new_data_comm       = mpi_settings.Has("data_communicator_name")

        if not creating_new_data_comm:
            # nothing to check if no new data-comm is created
            return

        if creating_new_data_comm and not solver_uses_custom_data_comm:
            import_settings.AddValue("data_communicator_name", mpi_settings["data_communicator_name"])
            if self.echo_level > 0:
                cs_tools.cs_print_info("KratosSolver", self._ClassName(), 'Using data commnicator with name "{}"'.format(mpi_settings["data_communicator_name"].GetString()))
            return

        if not solver_uses_custom_data_comm and not creating_new_data_comm:
            # using all ranks aka the default data comm hence nothing to do here
            return

        # check if both settings use the same DataCommunicator
        solver_data_comm_name = import_settings["data_communicator_name"].GetString()
        data_comm_creation_name = mpi_settings["data_communicator_name"].GetString()
        if solver_data_comm_name != data_comm_creation_name:
            err_msg  = 'Names of data communicators do not match!\n'
            err_msg += '    Name specified in "model_import_settings: {}\n'.format(solver_data_comm_name)
            err_msg += '    Name specified in "mpi_settings": {}'.format(data_comm_creation_name)
            raise Exception(err_msg)
