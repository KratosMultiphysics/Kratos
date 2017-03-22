import os

# Import Kratos
from KratosMultiphysics import *

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosExecuteMeshTest as ExecuteMeshTest

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


class MeshTestFactory(KratosUnittest.TestCase):

    def setUp(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # Get the ProjectParameters file
            parameter_file = open(self.file_name + "_parameters.json", 'r')
            ProjectParameters = Parameters(parameter_file.read())

            # Create the test
            self.test = ExecuteMeshTest.KratosExecuteMeshTest(ProjectParameters)

    def test_execution(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.test.Solve()

    def tearDown(self):
        pass


class MeshTest2D3NLaplacianTest(MeshTestFactory):
    file_name = "MeshTest2D3N/MeshTest2D3NLaplacianTest"


class MeshTest2D4NLaplacianTest(MeshTestFactory):
    file_name = "MeshTest2D4N/MeshTest2D4NLaplacianTest"


class MeshTest3D4NLaplacianTest(MeshTestFactory):
    file_name = "MeshTest3D4N/MeshTest3D4NLaplacianTest"


class MeshTest2D3NStructuralSimilarityTest(MeshTestFactory):
    file_name = "MeshTest2D3N/MeshTest2D3NStructuralSimilarityTest"


class MeshTest2D4NStructuralSimilarityTest(MeshTestFactory):
    file_name = "MeshTest2D4N/MeshTest2D4NStructuralSimilarityTest"


class MeshTest3D4NStructuralSimilarityTest(MeshTestFactory):
    file_name = "MeshTest3D4N/MeshTest3D4NStructuralSimilarityTest"
