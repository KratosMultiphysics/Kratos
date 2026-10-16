# Importing the Kratos Library
import KratosMultiphysics as KM
from KratosMultiphysics.process_factory import KratosProcessFactory
from KratosMultiphysics.analysis_stage import AnalysisStage

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.factories.solver_wrapper_factory as solver_wrapper_factory
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
import KratosMultiphysics.CoSimulationApplication.colors as colors
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_process import CoSimulationProcess

# Other imports
import sys

class CoSimulationAnalysis(AnalysisStage):
    """AnalysisStage of the CoSimulationApplication.
    It does NOT override the "RunSolutionLoop" method!
    """
    def __init__(self, cosim_settings, models=None):
        # Note: deliberately NOT calling the base-class constructor since arguments are different

        if not isinstance(cosim_settings, KM.Parameters):
            raise Exception("Input is expected to be provided as a Kratos Parameters object")

        self.cosim_settings = cosim_settings
        self.models = models
        self.model = KM.Model() # Dummy model, that will be used in the processes factory

        # this contains only the optional parameters, not the ones that have to be specified
        problem_data_defaults = KM.Parameters("""{
            "problem_name" : "default_co_simulation",
            "print_colors" : false,
            "echo_level"   : 1
        }""")

        problem_data = cosim_settings["problem_data"]

        problem_data.AddMissingParameters(problem_data_defaults)

        colors.PRINT_COLORS = problem_data["print_colors"].GetBool()
        self.echo_level = problem_data["echo_level"].GetInt()

        self.parallel_type = problem_data["parallel_type"].GetString()
        is_distributed_run = KM.IsDistributedRun()
        if self.parallel_type == "OpenMP":
            if is_distributed_run:
                cs_tools.cs_print_warning("Parallel Type", 'Specified "OpenMP" as "parallel_type", but Kratos is running in "MPI", please check your setup!')
        elif self.parallel_type == "MPI":
            if not is_distributed_run:
                cs_tools.cs_print_warning("Parallel Type", 'Specified "MPI" as "parallel_type", but Kratos is running in "OpenMP", please check your setup!')
        else:
            raise Exception('The "parallel_type" can be either "OpenMP" or "MPI"')

        if problem_data.Has("flush_stdout"):
            self.flush_stdout = problem_data["flush_stdout"].GetBool()
        else:
            # flush by default only in OpenMP, can decrease performance in MPI
            self.flush_stdout = (self.parallel_type == "OpenMP")

        self._GetSolver() # this creates the solver

    def Initialize(self):
        ## Here we initialize user-provided processes
        self.__CreateListOfProcesses()

        # Initialize
        for process in self._GetListOfProcesses():
            process.ExecuteBeforeCoSimulationInitialize()

        self._GetSolver().Initialize()
        # TODO: Separate the BeforeSolutionLoop() of the analysis from the Initialize() of the analysis and add the corresponding processes

        for process in self._GetListOfProcesses():
            process.ExecuteAfterCoSimulationInitialize()

        # Doing checks
        self._GetSolver().Check()
        for process in self._GetListOfProcesses():
            process.Check()

        if self.echo_level > 0:
            self._GetSolver().PrintInfo()

        ## Stepping and time settings
        self.end_time = self.cosim_settings["problem_data"]["end_time"].GetDouble()
        self.time = self.cosim_settings["problem_data"]["start_time"].GetDouble()
        self.step = 0

        if self.flush_stdout:
            CoSimulationAnalysis.Flush()

    def Finalize(self):
        # Finalize
        for process in self._GetListOfProcesses():
            process.ExecuteBeforeCoSimulationExecuteFinalize()

        self._GetSolver().Finalize()

        for process in self._GetListOfProcesses():
            process.ExecuteAfterCoSimulationExecuteFinalize()

    def InitializeSolutionStep(self):
        self.step += 1
        cs_tools.cs_print_info(colors.bold("\ntime={0:.12g}".format(self.time)+ " | step="+ str(self.step)))

        # InitializeSolutionStep
        for process in self._GetListOfProcesses():
            process.ExecuteBeforeCoSimulationInitializeSolutionStep()

        self._GetSolver().InitializeSolutionStep()

        for process in self._GetListOfProcesses():
            process.ExecuteAfterCoSimulationInitializeSolutionStep()

    def FinalizeSolutionStep(self):
        # FinalizeSolutionStep
        for process in self._GetListOfProcesses():
            process.ExecuteBeforeCoSimulationFinalizeSolutionStep()

        self._GetSolver().FinalizeSolutionStep()

        for process in self._GetListOfProcesses():
            process.ExecuteAfterCoSimulationFinalizeSolutionStep()

    def OutputSolutionStep(self):
        # FinalizeSolutionStep
        for process in self._GetListOfProcesses():
            process.ExecuteBeforeCoSimulationOutputSolutionStep()

        self._GetSolver().OutputSolutionStep()

        for process in self._GetListOfProcesses():
            process.ExecuteAfterCoSimulationOutputSolutionStep()

        if self.flush_stdout:
            CoSimulationAnalysis.Flush()

    def _GetSolver(self, solver_name=""):
        solver = super()._GetSolver()
        if solver_name == "":
            return solver
        else:
            return solver._GetSolver(solver_name)

    @staticmethod
    def Flush():
        sys.stdout.flush()
        KM.Logger.Flush()

    def _CreateSolver(self):
        """Create the solver
        """
        problem_name = self.cosim_settings["problem_data"]["problem_name"].GetString()
        return solver_wrapper_factory.CreateSolverWrapper(self.cosim_settings["solver_settings"], self.models, problem_name)

    def _GetSimulationName(self):
        """Returns the name of the Simulation
        """
        simulation_name = type(self).__name__ + ": " + self._GetSolver()._ClassName() + " coupling "
        list_solvers = self.cosim_settings["solver_settings"]["solvers"].keys()
        for solver_name in list_solvers:
            simulation_name += self._GetSolver(solver_name)._ClassName() + " - "
        simulation_name = simulation_name[0:-3] # Removing last " - "
        return simulation_name

    def __CreateListOfProcesses(self):
        """This function creates the processes and the output-processes
        """
        order_processes_initialization = self._GetOrderOfProcessesInitialization()
        self._list_of_processes        = self._CreateProcesses("processes", order_processes_initialization)

    def _CreateProcesses(self, parameter_name, initialization_order):
        """Create a list of processes
        Format:
        "processes" : {
            initial_processes : [
                { proces_specific_params },
                { proces_specific_params }
            ],
            boundary_processes : [
                { proces_specific_params },
                { proces_specific_params }
            ]
        }
        The order of intialization can be specified by setting it in "initialization_order"
        if e.g. the "boundary_processes" should be constructed before the "initial_processes", then
        initialization_order should be a list containing ["boundary_processes", "initial_processes"]
        see the functions _GetOrderOfProcessesInitialization and _GetOrderOfOutputProcessesInitialization
        """

        # First creating raw list
        list_of_processes_raw = []
        factory = KratosProcessFactory(self.model)
        if self.cosim_settings.Has(parameter_name):
            processes_params = self.cosim_settings[parameter_name]

            # first initialize the processes that depend on the order
            for processes_names in initialization_order:
                if processes_params.Has(processes_names):
                    list_of_processes_raw += factory.ConstructListOfProcesses(processes_params[processes_names])

            # then initialize the processes that don't depend on the order
            for name, value in processes_params.items():
                if not name in initialization_order:
                    list_of_processes_raw += factory.ConstructListOfProcesses(value) # Does this work? or should it be processes[name]

        # Filter only CoSimulationProcesses
        list_of_processes = []
        for process in list_of_processes_raw:
            if issubclass(process, CoSimulationProcess):
                list_of_processes.append(process)
            else:
                cs_tools.cs_print_info("Process", type(process).__name__ + " not added to the list of processes because is not derived from CoSimulationProcess")

        return list_of_processes

if __name__ == '__main__':
    if len(sys.argv) != 2:
        err_msg  = 'Wrong number of input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '    "python co_simulation_analysis.py <cosim-parameter-file>.json"\n'
        raise Exception(err_msg)

    parameter_file_name = sys.argv[1]

    with open(parameter_file_name,'r') as parameter_file:
        parameters = KM.Parameters(parameter_file.read())

    simulation = CoSimulationAnalysis(parameters)
    simulation.Run()
