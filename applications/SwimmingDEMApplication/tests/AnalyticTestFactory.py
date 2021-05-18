import os

# Importing the Kratos Library
import KratosMultiphysics

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
import tests_python_scripts.analytic_tests_scripts.analytic_test_analysis as analytic_test_analysis
MultipleGhostsTestAnalysis = analytic_test_analysis.MultipleGhostsTestAnalysis

# This utility will control the execution scope
class controlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

# General test factory
class TestFactory(KratosUnittest.TestCase):

    def setUp(self):
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # Setting parameters

            with open(self.file_parameters,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())

            # Create Model
            model = KratosMultiphysics.Model()

            self.test = MultipleGhostsTestAnalysis(model, parameters)

    def test_execution(self):
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.test.Run()

    def tearDown(self):
        pass