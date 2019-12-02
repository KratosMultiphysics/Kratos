from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from KratosMultiphysics.EmpireApplication.co_simulation_analysis import CoSimulationAnalysis

import sys

class CoSimulationSteadyAnalysis(CoSimulationAnalysis):

    def RunSolutionLoop(self):
        self.time = self._GetSolver().AdvanceInTime(self.time)
        self.InitializeSolutionStep()
        self._GetSolver().SolveSolutionStep()
        self.FinalizeSolutionStep()
        self.OutputSolutionStep()

        if self.flush_stdout:
            sys.stdout.flush()