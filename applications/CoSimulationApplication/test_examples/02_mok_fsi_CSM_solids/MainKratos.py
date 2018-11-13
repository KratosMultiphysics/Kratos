from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import KratosMultiphysics.CoSimulationApplication as CoSimulationApplication
import co_simulation_tools as cs_tools
from CoSimulationApplication import *
import sys

class CoSimulationAnalysis(object):
    """
    The base class for the CoSimulation-AnalysisStage
    """
    def __init__(self, cosim_settings):
        self.cosim_settings = cosim_settings
        self.flush_stdout = False

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
            sys.stdout.flush()

    def Initialize(self):
        self._GetSolver().Initialize()
        self._GetSolver().Check()

        if self.echo_level > 0:
            self._GetSolver().PrintInfo()

        ## Stepping and time settings
        self.end_time = self.cosim_settings["problem_data"]["end_time"].GetDouble()
        self.time = self.cosim_settings["problem_data"]["start_time"].GetDouble()
        self.step = 0

        sys.stdout.flush()

    def Finalize(self):
        self._GetSolver().Finalize()

    def PrintInfo(self):
        self._GetSolver().PrintInfo()

    def InitializeSolutionStep(self):
        print( cs_tools.bcolors.GREEN + cs_tools.bcolors.BOLD +"Time = {0:.10f}".format(self.time) + " | Step = " + str(self.step) + cs_tools.bcolors.ENDC )
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
        if("coupled_solver_settings" in self.cosim_settings.keys()):
            import custom_co_simulation_coupled_solvers.co_simulation_coupled_solver_factory as coupled_solver_factory
            self._solver = coupled_solver_factory.CreateCoupledSolver(self.cosim_settings)
            return self._solver
        elif ("solvers" in self.cosim_settings.keys()):
            num_solvers = len(self.cosim_settings["solvers"])
            if(num_solvers > 1 or num_solvers == 0):
                Exception("More than one or no solvers defined with out coupled solver !")
            else:
                import co_simulation_solver_factory as solver_factory
                self._solver = solver_factory.CreateSolverInterface(self.cosim_settings["solvers"][0])
                return self._solver



if __name__ == '__main__':
    from sys import argv

    if len(argv) != 2:
        err_msg =  'Wrong number of input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '    "python co_simulation_analysis.py <cosim-parameter-file>.json"\n'
        raise Exception(err_msg)

    parameter_file_name = argv[1]
    global cs_data_structure
    cs_data_structure = cs_tools.ImportDataStructure(parameter_file_name) 

    # Now we import actual parameters from the cs_data_structure
    with open(parameter_file_name,'r') as parameter_file:
        parameters = cs_data_structure.Parameters(parameter_file.read())

    simulation = CoSimulationAnalysis(parameters)
    simulation.Run()

