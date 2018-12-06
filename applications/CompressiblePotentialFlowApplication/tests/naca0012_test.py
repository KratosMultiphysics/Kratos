from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as UnitTest

# Other imports
import test_PotentialAnalysis
from test_PotentialAnalysis import PotentialFlowAnalysis
import os

class controlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

class WorkFolderScope:
    def __init__(self, work_folder):
        self.currentPath = os.getcwd()
        self.scope = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),work_folder))

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, exc_type, exc_value, traceback):
        os.chdir(self.currentPath)

class PotentialFlowTestFactory(UnitTest.TestCase):

    def setUp(self):
        print('hello world')
        # Within this location context:
        #with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):

        '''
        # Reading the ProjectParameters
        with open(self.file_name + "_parameters.json",'r') as parameter_file:
            print('self.file_name = ', self.file_name)
            ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())

        # To avoid many prints
        #if (ProjectParameters["problem_data"]["echo_level"].GetInt() == 0):
        #    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

        # Creating the test
        model = KratosMultiphysics.Model()
        self.test = test_PotentialAnalysis.PotentialAnalysis(model, ProjectParameters)
        self.test.Initialize()
        '''

    def test_execution(self):

        with WorkFolderScope(self.work_folder):
            self._run_test()
        

    def _run_test(self):
        model = KratosMultiphysics.Model()
        with open(self.file_name + "_parameters.json",'r') as parameter_file:
            print('self.file_name = ', self.file_name)
            ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())

        potential_flow_analysis = PotentialFlowAnalysis(model,ProjectParameters)
        potential_flow_analysis.Run()


class Naca0012Test(PotentialFlowTestFactory):
    file_name = "naca0012_Case_5"
    work_folder = "naca0012_tests"


if __name__ == '__main__':
    UnitTest.main()
