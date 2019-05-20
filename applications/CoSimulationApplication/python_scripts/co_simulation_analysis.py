from __future__ import print_function, absolute_import, division  # makes backward compatible with python 2.6 and 2.7

import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import ImportDataStructure


class CoSimulationAnalysis(object):
    """
    The base class for the CoSimulation-AnalysisStage
    """
    def __init__(self, settings):
        self.settings = settings
        self.flush_stdout = False

        self.echo_level = 0
        if "echo_level" in self.settings["problem_data"]:
            self.echo_level = self.settings["problem_data"]["echo_level"].GetInt()

        self.end_time = self.settings["problem_data"]["end_time"].GetDouble()
        self.time = self.settings["problem_data"]["start_time"].GetDouble()
        self.step = 0

    def Run(self):
        self.Initialize()
        self.RunSolutionLoop()
        self.Finalize()

    def Initialize(self):
        self._GetSolver().Initialize()
        self._GetSolver().Check()
        self._GetSolver().PrintInfo()

    def RunSolutionLoop(self):
        while self.time < self.end_time:
            self.InitializeSolutionStep()
            self.Predict()
            self.SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

    def Finalize(self):
        self._GetSolver().Finalize()

    def InitializeSolutionStep(self):
        self.step += 1
        self.time = self._GetSolver().AdvanceInTime(self.time)
        self._GetSolver().InitializeSolutionStep()

    def Predict(self):
        self._GetSolver().Predict()

    def SolveSolutionStep(self):
        self._GetSolver().SolveSolutionStep()

    def FinalizeSolutionStep(self):
        self._GetSolver().FinalizeSolutionStep()

    def OutputSolutionStep(self):
        self._GetSolver().OutputSolutionStep()

    def PrintInfo(self):
        self._GetSolver().PrintInfo()

    def _GetSolver(self):
        if not hasattr(self, '_solver'):
            self._solver = self._CreateSolver()
        return self._solver

    def _CreateSolver(self):
        if "coupled_solver_settings" in self.settings.keys():
            import KratosMultiphysics.CoSimulationApplication.custom_coupled_solvers.coupled_solver_factory as coupled_solver_factory
            self._solver = coupled_solver_factory.CreateCoupledSolver(self.settings)
            return self._solver
        elif "solvers" in self.settings.keys():
            num_solvers = len(self.settings["solvers"])
            if num_solvers > 1 or num_solvers == 0:
                raise Exception("More than one or no solvers defined without coupled solver!")
            else:
                import KratosMultiphysics.CoSimulationApplication.custom_solver_interfaces.co_simulation_solver_factory as solver_factory
                self._solver = solver_factory.CreateSolverInterface("Solver", self.settings["solvers"][0])
                return self._solver


if __name__ == '__main__':
    from sys import argv

    # Check number of command line arguments
    if len(argv) != 2:
        err_msg = 'Wrong number of input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '    "python co_simulation_analysis.py <cosim-parameter-file>.json"\n'
        raise Exception(err_msg)

    # Import data structure
    parameter_file_name = argv[1]
    cs_data_structure = ImportDataStructure(parameter_file_name)

    # Import parameters using the data structure
    with open(parameter_file_name, 'r') as parameter_file:
        parameters = cs_data_structure.Parameters(parameter_file.read())

    simulation = CoSimulationAnalysis(parameters)
    simulation.Run()