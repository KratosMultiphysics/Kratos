import os

# Import Kratos
from KratosMultiphysics import *

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosExecuteIGATest as ExecuteTest
#import KratosExecuteConvergenceAcceleratorTest as ExecuteConvergenceAcceleratorTest

# This utiltiy will control the execution scope in case we need to acces files or we depend
# on specific relative locations of the files.

# TODO: Should we move this to KratosUnittest?
class controlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)
		
class StructuralMechanichsTestFactory(KratosUnittest.TestCase):
    def setUp(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # Initialize GiD  I/O
            parameter_file = open(self.file_name + "_parameters.json", 'r')
            ProjectParameters = Parameters(parameter_file.read())

            # Creating the model part
            self.test = ExecuteTest.KratosExecuteIGATest(ProjectParameters)

    def test_execution(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.test.Solve()

    def tearDown(self):
        pass		
		
class Penalty_Test(StructuralMechanichsTestFactory):
	file_name = "Shell_Test/Penalty_Test"

class Lagrange_Test(StructuralMechanichsTestFactory):
	file_name = "Shell_Test/Lagrange_Test"	