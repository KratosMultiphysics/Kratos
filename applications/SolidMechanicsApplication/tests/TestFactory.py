import os

# Import Kratos
import KratosMultiphysics

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
import MainSolid

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
            self.model = KratosMultiphysics.Model()
            if( self.file_parameters == None ):
                self.file_parameters = self.file_name + "_parameters.json"
            self.test = MainSolid.Solution(self.model, self.file_parameters, self.file_name)

    def test_execution(self):
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.test.Run()

    def tearDown(self):
        pass
