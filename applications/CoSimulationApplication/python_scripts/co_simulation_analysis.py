# Importing the Kratos Library
import KratosMultiphysics as KM
from KratosMultiphysics.analysis_stage import AnalysisStage

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.factories.solver_wrapper_factory as solver_wrapper_factory
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
import KratosMultiphysics.CoSimulationApplication.colors as colors

# Other imports
import sys

class CoSimulationAnalysis(AnalysisStage):
    """AnalysisStage of the CoSimulationApplication.
    It contains all necessary modifications to make CoSimulation work both with Kratos and pyKratos
    It does NOT override the "RunSolutionLoop" method!
    """
    def __init__(self, cosim_settings, models=None):
        # Note: deliberately NOT calling the base-class constructor since arguments are different

        if not isinstance(cosim_settings, KM.Parameters):
            raise Exception("Input is expected to be provided as a Kratos Parameters object")

        self.cosim_settings = cosim_settings
        self.models = models

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
        self._GetSolver().Initialize()
        self._GetSolver().InitializeCouplingInterfaceData()
        self._GetSolver().Check()

        if self.echo_level > 0:
            self._GetSolver().PrintInfo()

        ## Stepping and time settings
        self.end_time = self.cosim_settings["problem_data"]["end_time"].GetDouble()
        self.time = self.cosim_settings["problem_data"]["start_time"].GetDouble()
        self.step = 0

        if self.flush_stdout:
            CoSimulationAnalysis.Flush()

    def Finalize(self):
        self._GetSolver().Finalize()

    def InitializeSolutionStep(self):
        self.step += 1
        cs_tools.cs_print_info(colors.bold("\ntime={0:.12g}".format(self.time)+ " | step="+ str(self.step)))

        self._GetSolver().InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        self._GetSolver().FinalizeSolutionStep()

    def OutputSolutionStep(self):
        self._GetSolver().OutputSolutionStep()

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
