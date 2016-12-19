import os

# Import Kratos
from KratosMultiphysics import *

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosExecuteEmbeddedTest as ExecuteEmbeddedTest

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


class EmbeddedTestFactory(KratosUnittest.TestCase):

    def setUp(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # Get the ProjectParameters file
            parameter_file = open(self.file_name + "_parameters.json", 'r')
            ProjectParameters = Parameters(parameter_file.read())

            # Create the test
            self.test = ExecuteEmbeddedTest.KratosExecuteEmbeddedTest(ProjectParameters)

    def test_execution(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.test.Solve()

    def tearDown(self):
        pass
#
#
# class EmbeddedReservoirTestFactory(KratosUnittest.TestCase):
#
#     def setUp(self):
#         # Within this location context:
#         with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
#             # Get the ProjectParameters file
#             parameter_file = open(self.file_name + "_parameters.json", 'r')
#             ProjectParameters = Parameters(parameter_file.read())
#
#             # Create the test
#             self.test = ExecuteEmbeddedTest.KratosExecuteEmbeddedTest(ProjectParameters)
#
#     def test_execution(self):
#         # Within this location context:
#         with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
#             self.test.Solve()
#
#     def tearDown(self):
#         pass


class EmbeddedArtificialCompressibilityTest(EmbeddedTestFactory):
    file_name = "EmbeddedArtificialCompressibilityTest/EmbeddedArtificialCompressibilityTest"


class EmbeddedCouetteTest(EmbeddedTestFactory):
    file_name = "EmbeddedCouetteTest/EmbeddedCouetteTest"


class EmbeddedCouetteImposedTest(EmbeddedTestFactory):
    file_name = "EmbeddedCouetteImposedTest/EmbeddedCouetteImposedTest"


class EmbeddedReservoirTest(EmbeddedTestFactory):
    file_name = "EmbeddedReservoirTest/EmbeddedReservoirTest"
