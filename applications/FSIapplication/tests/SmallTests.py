import os

# Import Kratos
from KratosMultiphysics import *

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosExecuteMapperTest as ExecuteMapperTest
import KratosExecuteConvergenceAcceleratorTest as ExecuteConvergenceAcceleratorTest
import KratosExecuteFSIProblemEmulatorTest as ExecuteFSIProblemEmulatorTest

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


class MapperTestFactory(KratosUnittest.TestCase):

    def setUp(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # Get the ProjectParameters file
            parameter_file = open(self.file_name + "_parameters.json", 'r')
            ProjectParameters = Parameters(parameter_file.read())

            # Create the test
            self.test = ExecuteMapperTest.KratosExecuteMapperTest(ProjectParameters)

    def test_execution(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.test.Solve()

    def tearDown(self):
        pass


class FSIProblemEmulatorTestFactory(KratosUnittest.TestCase):

    def setUp(self):
        self.test_list = []
        
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # Iterate in the convergence accelerators test list
            for parameter_file_name in self.file_name_list:
                # Get the ProjectParameters file
                parameter_file = open(parameter_file_name + "_parameters.json", 'r')
                ProjectParameters = Parameters(parameter_file.read())

                # Create the test
                self.test_list.append(ExecuteFSIProblemEmulatorTest.KratosExecuteFSIProblemEmulatorTest(ProjectParameters))

    def test_execution(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # Iterate in the convergence accelerators test list
            for test in self.test_list:
                test.Solve()

    def tearDown(self):
        pass
        

class NonConformantOneSideMap2D_test1(MapperTestFactory):
    file_name = "NonConformantOneSideMap2D_test1/NonConformantOneSideMap2D_test1"


class NonConformantOneSideMap2D_test2(MapperTestFactory):
    file_name = "NonConformantOneSideMap2D_test2/NonConformantOneSideMap2D_test2"
                      
                
class FSIProblemEmulatorTest(FSIProblemEmulatorTestFactory):
    file_name_1 = "FSIProblemEmulatorTest/FSIProblemEmulatorTest_Aitken"
    file_name_2 = "FSIProblemEmulatorTest/FSIProblemEmulatorTest_MVQN"
    file_name_3 = "FSIProblemEmulatorTest/FSIProblemEmulatorTest_MVQN_recursive"
    
    file_name_list = [file_name_1, file_name_2, file_name_3]
