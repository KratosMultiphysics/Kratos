import os

# Import Kratos
from KratosMultiphysics import *

try:
    import KratosMultiphysics.ExternalSolversApplication
    have_external_solvers = True
except ImportError as e:
    have_external_solvers = False

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosExecuteEmbeddedTest as ExecuteEmbeddedTest
import KratosExecuteManufacturedSolutionTest as ExecuteManufacturedSolutionTest

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

class ManufacturedSolutionTestFactory(KratosUnittest.TestCase):

    def setUp(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # Get the ProjectParameters file
            parameter_file = open(self.file_name + "_parameters.json", 'r')
            ProjectParameters = Parameters(parameter_file.read())

            # Create the test
            self.test = ExecuteManufacturedSolutionTest.KratosExecuteManufacturedSolutionTest(ProjectParameters)

    def test_execution(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.test.Solve()

    def tearDown(self):
        pass


class EmbeddedArtificialCompressibilityTest(EmbeddedTestFactory):
    file_name = "EmbeddedArtificialCompressibilityTest/EmbeddedArtificialCompressibilityTest"


class EmbeddedCouette2DTest(EmbeddedTestFactory):
    file_name = "EmbeddedCouette2DTest/EmbeddedCouette2DTest"


class EmbeddedCouette3DTest(EmbeddedTestFactory):
    file_name = "EmbeddedCouette3DTest/EmbeddedCouette3DTest"


class EmbeddedCouette2DImposedTest(EmbeddedTestFactory):
    file_name = "EmbeddedCouetteImposed2DTest/EmbeddedCouetteImposed2DTest"


class EmbeddedCouette3DImposedTest(EmbeddedTestFactory):
    file_name = "EmbeddedCouetteImposed3DTest/EmbeddedCouetteImposed3DTest"


class EmbeddedReservoirTest(EmbeddedTestFactory):
    file_name = "EmbeddedReservoirTest/EmbeddedReservoirTest"


class EmbeddedSlipBoundaryConditionTest(EmbeddedTestFactory):
    file_name = "EmbeddedSlipBoundaryConditionTest/EmbeddedSlipBoundaryConditionTest"


class EmbeddedSlipReservoirTest(EmbeddedTestFactory):
    file_name = "EmbeddedSlipReservoirTest/EmbeddedSlipReservoirTest"

    
@KratosUnittest.skipUnless(have_external_solvers, "Missing required application: ExternalSolversApplication")
class ManufacturedSolutionTest(ManufacturedSolutionTestFactory):
    file_name = "ManufacturedSolutionTest/ManufacturedSolutionTest"


class NavierStokesWallConditionTest(EmbeddedTestFactory):
    file_name = "NavierStokesWallConditionTest/NavierSokesWallConditionTest"
