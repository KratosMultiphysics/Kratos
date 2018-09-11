from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import sys
import co_simulation_tools as cs_tools

class CoSimulationAnalysis(object):
    """
    The base class for the CoSimulation-AnalysisStage
    """
    def __init__(self, cosim_settings):
        if (type(cosim_settings) != dict):
            raise Exception("Input is expected to be provided as a python dictionary")

        self.cosim_settings = cosim_settings
        self.flush_stdout = False

        if "parallel_type" in self.cosim_settings["problem_data"]:
            parallel_type = self.cosim_settings["problem_data"]["parallel_type"]
            if parallel_type == "OpenMP":
                self.flush_stdout = True
                cs_tools.PRINTING_RANK = True
            elif parallel_type == "MPI":
                self.flush_stdout = False
                cs_tools.COSIM_SPACE = CoSimulationMPISpace()
                cs_tools.PRINTING_RANK = (cs_tools.COSIM_SPACE.Rank() == 0)
            else:
                raise Exception('"parallel_type" can only be "OpenMP" or "MPI"!')

        if "flush_terminal" in self.cosim_settings["problem_data"]:
            self.flush_stdout = self.cosim_settings["problem_data"]["parallel_type"]

        self.echo_level = 0
        if "echo_level" in self.cosim_settings["problem_data"]:
            self.echo_level = self.cosim_settings["problem_data"]["echo_level"]

    def Run(self):
        self.Initialize()
        self.RunSolutionLoop()
        self.Finalize()

    def RunSolutionLoop(self):
        print("")
        while self.time < self.end_time:
            print("")
            self.step += 1
            self.time = self._GetSolver().AdvanceInTime(self.time)
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            self._GetSolver().SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

            if self.flush_stdout:
                sys.stdout.flush()

    def Initialize(self):
        self._GetSolver().Initialize()
        self._GetSolver().Check()

        if self.echo_level > 0:
            self._GetSolver().PrintInfo()

        ## Stepping and time settings
        self.end_time = self.cosim_settings["problem_data"]["end_time"]
        self.time = self.cosim_settings["problem_data"]["start_time"]
        self.step = 0

        if self.flush_stdout:
            sys.stdout.flush()

    def Finalize(self):
        self._GetSolver().Finalize()

    def InitializeSolutionStep(self):
        csprint(0, bold("time={0:.12g}".format(self.time)+ " | step="+ str(self.step)))

        self._GetSolver().InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        self._GetSolver().FinalizeSolutionStep()

    def OutputSolutionStep(self):
        self._GetSolver().OutputSolutionStep()

    def _GetSolver(self):
        if not hasattr(self, '_solver'):
            self._solver = self._CreateSolver()
        return self._solver

    def _CreateSolver(self):
        """Create the solver
        """
        import co_simulation_solvers.python_solvers_wrapper_co_simulation as solvers_wrapper
        return solvers_wrapper.CreateSolver(self.cosim_settings["solver_settings"], level=0)

if __name__ == '__main__':
    from sys import argv
    import json

    if len(argv) != 2:
        err_msg =  'Wrong number of input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '    "python co_simulation_analysis.py <cosim-parameter-file>.json"\n'
        raise Exception(err_msg)

    parameter_file_name = argv[1]

    with open(parameter_file_name,'r') as parameter_file:
        parameters = json.load(parameter_file)

    simulation = CoSimulationAnalysis(parameters)
    simulation.Run()
