import os

# Importing the Kratos Library
from KratosMultiphysics import Model, Parameters, Logger

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
from tests_python_scripts.recovery_scripts.linear_standard_test import LinearStandardTestAnalysis

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
                parameters = Parameters(parameter_file.read())

            # Create Model
            model = Model()

            # To avoid too many prints
            Logger.GetDefaultOutput().SetSeverity(Logger.Severity.DETAIL)

            self.test = LinearStandardTestAnalysis(model, parameters)

    def test_execution(self):
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.test.Run()

    def tearDown(self):
        pass
