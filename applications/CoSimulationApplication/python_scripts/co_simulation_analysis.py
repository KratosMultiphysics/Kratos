from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM
from KratosMultiphysics.analysis_stage import AnalysisStage

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.factories.solver_wrapper_factory as solvers_wrapper_factory
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
import KratosMultiphysics.CoSimulationApplication.colors as colors

# Other imports
import sys

class CoSimulationAnalysis(AnalysisStage):
    """AnalysisStage of the CoSimulationApplication.
    It contains all necessary modifications to make CoSimulation work both with Kratos and pyKratos
    It does NOT override the "RunSolutionLoop" method!
    """
    def __init__(self, cosim_settings):
        # Note: deliberately NOT calling the base-class constructor,
        # since this would break the python-only version due to the type-checks

        self.cosim_settings = cosim_settings

        problem_data_defaults = KM.Parameters("""{
            "problem_name" : "default_co_simulation",
            "print_colors" : false,
            "flush_stdout" : false,
            "echo_level"   : 0
        }""")

        problem_data = cosim_settings["problem_data"]

        problem_data.AddMissingParameters(problem_data_defaults)

        colors.PRINT_COLORS = problem_data["print_colors"].GetBool()
        self.flush_stdout = problem_data["flush_stdout"].GetBool()
        self.echo_level = problem_data["echo_level"].GetInt()

    def Initialize(self):
        self._GetSolver().Initialize()
        self._GetSolver().Check()

        if self.echo_level > 0:
            self._GetSolver().PrintInfo()

        ## Stepping and time settings
        self.end_time = self.cosim_settings["problem_data"]["end_time"].GetDouble()
        self.time = self.cosim_settings["problem_data"]["start_time"].GetDouble()
        self.step = 0

        if self.flush_stdout:
            sys.stdout.flush()

    def Finalize(self):
        self._GetSolver().Finalize()

    def InitializeSolutionStep(self):
        self.step += 1
        cs_tools.cs_print_info(colors.bold("time={0:.12g}".format(self.time)+ " | step="+ str(self.step)))

        self._GetSolver().InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        self._GetSolver().FinalizeSolutionStep()

    def OutputSolutionStep(self):
        self._GetSolver().OutputSolutionStep()

        if self.flush_stdout:
            sys.stdout.flush()

    def _CreateSolver(self):
        """Create the solver
        """
        problem_name = self.cosim_settings["problem_data"]["problem_name"].GetString()
        return solvers_wrapper_factory.CreateSolverWrapper(self.cosim_settings["solver_settings"], problem_name)

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
