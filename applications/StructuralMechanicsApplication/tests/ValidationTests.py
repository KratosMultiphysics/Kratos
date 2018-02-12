import os

# Import Kratos
from KratosMultiphysics import *

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
import Kratos_Execute_Structural_Test as Execute_Test

# This utility will control the execution scope in case we need to access files or we depend
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
    
class SprismPanTests(StructuralMechanichsTestFactory):
    file_name = "sprism_test/pan_test"
    
class PendulusTLTest(StructuralMechanichsTestFactory):
    file_name = "pendulus_test/pendulus_TL_test"
    
class PendulusULTest(StructuralMechanichsTestFactory):
    file_name = "pendulus_test/pendulus_UL_test"

class ShellT3AndQ4LinearStaticUnstructScordelisLoRoofTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_linear_static_unstruct_scordelis_lo_roof"

class ShellT3AndQ4LinearStaticUnstructUnstructPinchedCylinderTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_linear_static_unstruct_pinched_cylinder"

class ShellT3AndQ4LinearStaticUnstructPinchedHemisphereTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_linear_static_unstruct_pinched_hemisphere"

class ShellT3AndQ4LinearStaticUnstructClampedCylinderOrthotropicTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_linear_static_unstruct_clamped_cylinder_orthotropic"

class ShellT3AndQ4NonLinearStaticUnstructHingedCylRoofSnapthroughTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_nonlinear_static_unstruct_hinged_cyl_roof_snapthrough"

class ShellT3AndQ4NonLinearStaticUnstructHingedCylRoofSnapthroughOrthotropicTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_nonlinear_static_unstruct_hinged_cyl_roof_snapthrough_orthotropic"

class ShellT3AndQ4NonLinearDynamicUnstructOscillatingPlateTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_nonlinear_dynamic_unstruct_oscillating_plate"

class ShellT3AndQ4NonLinearDynamicUnstructOscillatingPlateLumpedTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_nonlinear_dynamic_unstruct_oscillating_plate_lumped"

class ShellT3AndQ4NonLinearDynamicStructPendulusTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_nonlinear_dynamic_struct_pendulus"

class ShellT3AndQ4NonLinearDynamicStructPendulusLumpedTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_nonlinear_dynamic_struct_pendulus_lumped"

class ShellT3AndQ4NonLinearDynamicUnstructPendulusTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_nonlinear_dynamic_unstruct_pendulus"

class ShellT3AndQ4NonLinearDynamicUnstructPendulusLumpedTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_nonlinear_dynamic_unstruct_pendulus_lumped"
