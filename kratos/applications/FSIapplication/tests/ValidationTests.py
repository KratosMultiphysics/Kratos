import os

# Import Kratos
from KratosMultiphysics import *

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosExecuteMokBenchmark as ExecuteMokBenchmark

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


class MokBenchmarkFactory(KratosUnittest.TestCase):

    def setUp(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # Get the ProjectParameters file
            parameter_file = open(self.file_name, 'r')
            ProjectParameters = Parameters(parameter_file.read())

            # Create the test
            self.test = ExecuteMokBenchmark.KratosExecuteMokBenchmark(ProjectParameters)

    def test_execution(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.test.Solve()

    def tearDown(self):
        pass


class MokBenchmarkTest(MokBenchmarkFactory):
    file_name = "MokBenchmark/ProjectParameters.json"
