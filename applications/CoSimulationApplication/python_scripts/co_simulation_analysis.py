from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

from KratosMultiphysics.analysis_stage import AnalysisStage
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import csprint, bold, CheckCoSimulationSettingsAndAssignDefaults

import sys

class CoSimulationAnalysis(AnalysisStage):
    def __init__(self, model, cosim_settings):
        # Note: deliberately NOT calling the base-class constructor, since this would
        # break the python-only version due to the type-checks

        self.model = model

        # CheckCoSimulationSettingsAndAssignDefaults(cosim_settings)

        problem_data = cosim_settings["problem_data"]
        self.cosim_settings = cosim_settings

        cs_tools.PRINT_COLORS = problem_data["print_colors"].GetBool()

        self.flush_stdout = False#problem_data["flush_stdout"].GetBool()
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
        csprint(0, bold("time={0:.12g}".format(self.time)+ " | step="+ str(self.step)))

        self._GetSolver().InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        self._GetSolver().FinalizeSolutionStep()

    def OutputSolutionStep(self):
        self._GetSolver().OutputSolutionStep()

    def _CreateSolver(self):
        """Create the solver
        """
        problem_name = "default"
        if "problem_name" in self.cosim_settings["problem_data"]:
            problem_name = self.cosim_settings["problem_data"]["problem_name"].GetString()
        import KratosMultiphysics.CoSimulationApplication.coupled_solvers.co_simulation_coupled_solver_factory as solvers_wrapper
        return solvers_wrapper.CreateCoupledSolver(self.model, self.cosim_settings["solver_settings"], problem_name)

if __name__ == '__main__':
    from sys import argv

    if len(argv) != 2:
        err_msg  = 'Wrong number of input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '    "python co_simulation_analysis.py <cosim-parameter-file>.json"\n'
        raise Exception(err_msg)

    parameter_file_name = argv[1]
    global cs_data_structure
    cs_data_structure = cs_tools.ImportDataStructure(parameter_file_name)

    # Now we import actual parameters from the cs_data_structure
    with open(parameter_file_name,'r') as parameter_file:
        # import json
        # parameters = json.load(parameter_file)
        parameters = cs_data_structure.Parameters(parameter_file.read())

    print(type(parameters))
    model = cs_data_structure.Model()

    simulation = CoSimulationAnalysis(model, parameters)
    simulation.Run()



class BaseSolverInterface(object):
    pass


class BaseCompositeSolver(BaseSolverInterface):
    pass