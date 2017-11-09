import os

# Import Kratos
from KratosMultiphysics import *

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
import Kratos_Execute_Structural_Test as Execute_Test

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
    
class IsotropicDamageSimoJuPSTest(StructuralMechanichsTestFactory):
    file_name = "cl_test/IsotropicDamageSimoJu/PlaneStress_FourPointShear_test"

class ShellQ4ThickBendingRollUpTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_Thick__BendingRollUp_test"

class ShellQ4ThickDrillingRollUpTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_Thick__DrillingRollUp_test"

class ShellT3ThinBendingRollUpTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_Thin__BendingRollUp_test"

class ShellT3ThinDrillingRollUpTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_Thin__DrillingRollUp_test"

class ShellT3IsotropicScordelisTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_Isotropic_Scordelis_test"

class ShellT3ThickLinearStaticTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_Thick_linear_static_test"

class ShellT3ThickNonLinearStaticTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_Thick_nonlinear_static_test"

class ShellT3ThickLinearDynamicTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_Thick_linear_dynamic_test"

class ShellT3ThickNonLinearDynamicTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_Thick_nonlinear_dynamic_test"

class ShellQ4ThinLinearStaticTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_Thin_linear_static_test"

class ShellQ4ThinNonLinearStaticTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_Thin_nonlinear_static_test"

class ShellQ4ThinLinearDynamicTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_Thin_linear_dynamic_test"

class ShellQ4ThinNonLinearDynamicTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_Thin_nonlinear_dynamic_test"
