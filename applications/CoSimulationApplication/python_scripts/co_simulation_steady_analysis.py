from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM
from KratosMultiphysics.analysis_stage import AnalysisStage

# CoSimulation imports
from KratosMultiphysics.CoSimulationApplication.co_simulation_analysis import CoSimulationAnalysis
# Other imports
import sys

class CoSimulationSteadyAnalysis(CoSimulationAnalysis):
    """AnalysisStage of the CoSimulationApplication.
    It contains all necessary modifications to make CoSimulation work both with Kratos and pyKratos
    It does NOT override the "RunSolutionLoop" method!
    """
    def RunSolutionLoop(self):
        self.InitializeSolutionStep()
        self._GetSolver().SolveSolutionStep()
        self.FinalizeSolutionStep()
        self.OutputSolutionStep()

        if self.flush_stdout:
            sys.stdout.flush()


if __name__ == '__main__':
    if len(sys.argv) != 2:
        err_msg  = 'Wrong number of input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '    "python co_simulation_analysis.py <cosim-parameter-file>.json"\n'
        raise Exception(err_msg)

    parameter_file_name = sys.argv[1]

    with open(parameter_file_name,'r') as parameter_file:
        parameters = KM.Parameters(parameter_file.read())

    simulation = CoSimulationSteadyAnalysis(parameters)
    simulation.Run()
