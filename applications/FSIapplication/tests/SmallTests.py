import os

# Import Kratos
from KratosMultiphysics import *

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosExecuteConvergenceAcceleratorTest as ExecuteConvergenceAcceleratorTest
import KratosExecuteFSIProblemEmulatorTest as ExecuteFSIProblemEmulatorTest

try:
    from KratosMultiphysics.StructuralMechanicsApplication import *
    from KratosMultiphysics.FluidDynamicsApplication import *
    missing_external_dependencies = False
    missing_application = ''
except ImportError as e:
    missing_external_dependencies = True
    # extract name of the missing application from the error message
    import re
    missing_application = re.search(r'''.*'KratosMultiphysics\.(.*)'.*''','{0}'.format(e)).group(1)

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


@KratosUnittest.skipIf(missing_external_dependencies, "Missing required application: {0}".format(missing_application))
class FSIProblemEmulatorTest(FSIProblemEmulatorTestFactory):
    file_name_1 = "FSIProblemEmulatorTest/FSIProblemEmulatorTest_Aitken"
    file_name_2 = "FSIProblemEmulatorTest/FSIProblemEmulatorTest_MVQN"
    file_name_3 = "FSIProblemEmulatorTest/FSIProblemEmulatorTest_MVQN_recursive"

    file_name_list = [file_name_1, file_name_2, file_name_3]
