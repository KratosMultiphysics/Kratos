from __future__ import print_function, absolute_import, division  # makes backward compatible with python 2.6 and 2.7

from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import ImportDataStructure
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools


class CoSimulationAnalysis(object):
    def __init__(self, parameters):
        self.parameters = parameters
        self.settings = parameters["settings"]

        self.echo_level = self.settings["echo_level"].GetInt()

        self.start_step = self.settings["start_step"].GetInt()
        self.stop_step = self.settings["stop_step"].GetInt()
        self.step = self.start_step

        self.coupled_solver = cs_tools.CreateInstance(self.parameters["coupled_solver"])

    def Run(self):
        self.Initialize()
        self.RunSolutionLoop()
        self.Finalize()

    def Initialize(self):
        self.coupled_solver.Initialize()
        self.coupled_solver.Check()
        self.coupled_solver.PrintInfo()

    def RunSolutionLoop(self):
        for self.step in range(self.start_step, self.stop_step):
            self.coupled_solver.InitializeSolutionStep()
            self.coupled_solver.SolveSolutionStep()
            self.coupled_solver.FinalizeSolutionStep()
            self.coupled_solver.OutputSolutionStep()

    def Finalize(self):
        self.coupled_solver.Finalize()


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
