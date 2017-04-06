import os

# Import Kratos
from KratosMultiphysics import *

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
import Kratos_Execute_Solid_Test as Execute_Test

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
            self.test = Execute_Test.Kratos_Execute_Test(ProjectParameters)

    def test_execution(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.test.Solve()

    def tearDown(self):
        pass


class SimpleMeshMovingTest(StructuralMechanichsTestFactory):
    file_name = "mesh_moving_test/simple_mesh_moving_test"


class DynamicBossakTests(StructuralMechanichsTestFactory):
    file_name = "dynamic_test/dynamic_bossak_test"


class DynamicNewmarkTests(StructuralMechanichsTestFactory):
    file_name = "dynamic_test/dynamic_newmark_test"


class SprismMembranePatchTests(StructuralMechanichsTestFactory):
    file_name = "sprism_test/patch_membrane_test"


class SprismBendingPatchTests(StructuralMechanichsTestFactory):
    file_name = "sprism_test/patch_bending_test"


class ShellQ4ThickBendingRollUpTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_Thick__BendingRollUp_test"


class ShellQ4ThickDrillingRollUpTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_Thick__DrillingRollUp_test"


class ShellT3ThinBendingRollUpTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_Thin__BendingRollUp_test"


class ShellT3ThinDrillingRollUpTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_Thin__DrillingRollUp_test"


class EigenQ4Thick2x2PlateTests(StructuralMechanichsTestFactory):
    file_name = "eigen_test/Eigen_Q4_Thick_2x2_Plate_test"


class EigenTL3D8NCubeTests(StructuralMechanichsTestFactory):
    file_name = "eigen_test/Eigen_TL_3D8N_Cube_test"

class Simple3D2NTrussTest(StructuralMechanichsTestFactory):
    file_name = "3D2N_TRUSS_test/nonlinear_3D2NTruss_test"

class Simple3D2NBeamCrTest(StructuralMechanichsTestFactory):
    file_name = "3D2N_BEAMCR_test/nonlinear_3D2NBeamCr_test"
     