from __future__ import print_function, absolute_import, division
import KratosMultiphysics as KM

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

from KratosMultiphysics.CoSimulationApplication.co_simulation_analysis import CoSimulationAnalysis

from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import UsingPyKratos
using_pykratos = UsingPyKratos()

import os, subprocess

class TestFLOWerCoupling(KratosUnittest.TestCase):
    '''TODO add description
    '''

    def test_dummy_flower_solver(self):
        self._createTest("FLOWer_coupling", "cosim", "dummy_flower_solver")
        self._runTest()

    def _createTest(self, problem_dir_name, parameter_file_name, ext_parameter_file_name):
        self.problem_dir_name = problem_dir_name

        full_parameter_file_name = os.path.join(problem_dir_name, parameter_file_name + '_parameters.json')
        self.ext_parameter_file_name = os.path.join(problem_dir_name, ext_parameter_file_name + '_parameters.json')

        with open(full_parameter_file_name, 'r') as parameter_file:
            self.cosim_parameters = KM.Parameters(parameter_file.read())

        if not using_pykratos:
            # To avoid many prints
            echo_level = self.cosim_parameters["problem_data"]["echo_level"].GetInt()
            if (echo_level == 0):
                KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
            else:
                KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.INFO)

    def _runTest(self):
        p = subprocess.Popen(
            ["python3", "../python_scripts/helpers/dummy_flower_solver.py", self.ext_parameter_file_name],
            cwd=os.path.dirname(os.path.abspath(__file__)))

        CoSimulationAnalysis(self.cosim_parameters).Run()

        p.communicate()

        kratos_utils.DeleteTimeFiles(self.problem_dir_name)

    # called only once for this class, opposed of tearDown()
    @classmethod
    def tearDownClass(cls):
        # Cleaning
        kratos_utils.DeleteDirectoryIfExisting("__pycache__")


if __name__ == '__main__':
    KratosUnittest.main()
