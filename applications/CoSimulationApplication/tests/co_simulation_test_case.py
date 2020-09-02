import KratosMultiphysics as KM

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

from KratosMultiphysics.CoSimulationApplication.co_simulation_analysis import CoSimulationAnalysis

from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import UsingPyKratos
using_pykratos = UsingPyKratos()

import os, subprocess

class CoSimulationTestCase(KratosUnittest.TestCase):
    '''This class is the basis for the testing the framework
    It can be used to test complete cases with the "CoSimulation-Analysis"
    '''

    def _createTest(self, problem_dir_name, parameter_file_name):
        self.problem_dir_name = problem_dir_name

        full_parameter_file_name = os.path.join(problem_dir_name, parameter_file_name + '_parameters.json')

        with open(full_parameter_file_name, 'r') as parameter_file:
            self.cosim_parameters = KM.Parameters(parameter_file.read())

        if not using_pykratos:
            # To avoid many prints
            echo_level = self.cosim_parameters["problem_data"]["echo_level"].GetInt()
            if (echo_level == 0):
                KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
            else:
                KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.INFO)

        self.addCleanup(kratos_utils.DeleteTimeFiles, self.problem_dir_name)

    def _runTest(self):
        CoSimulationAnalysis(self.cosim_parameters).Run()

    def _runTestWithExternal(self, subprocess_args_list):
        self.addCleanup(self.__TerminateSubProcessIfRunning)

        self.p = subprocess.Popen(subprocess_args_list)

        CoSimulationAnalysis(self.cosim_parameters).Run()

        self.p.communicate()

        self.assertEqual(self.p.returncode, 0)

    def __TerminateSubProcessIfRunning(self):
        poll = self.p.poll()
        if poll == None:
            self.p.terminate()


    # called only once for this class, opposed of tearDown()
    @classmethod
    def tearDownClass(cls):
        # Cleaning
        kratos_utils.DeleteDirectoryIfExisting("__pycache__")