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


class StructrualMechanichsTestFactory(KratosUnittest.TestCase):

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


class DynamicBossakTests(StructrualMechanichsTestFactory):
    file_name = "dynamic_test/dynamic_bossak_test"


class DynamicNewmarkTests(StructrualMechanichsTestFactory):
    file_name = "dynamic_test/dynamic_newmark_test"


class SprismMembranePatchTests(StructrualMechanichsTestFactory):
    file_name = "sprism_test/patch_membrane_test"


class SprismBendingPatchTests(StructrualMechanichsTestFactory):
    file_name = "sprism_test/patch_bending_test"


class ShellQ4ThickBendingRollUpTests(StructrualMechanichsTestFactory):
    file_name = "Shell_Q4_Thick__BendingRollUp_test/Shell_Q4_Thick__BendingRollUp_test"


class ShellQ4ThickDrillingRollUpTests(StructrualMechanichsTestFactory):
    file_name = "Shell_Q4_Thick__DrillingRollUp_test/Shell_Q4_Thick__DrillingRollUp_test"


class ShellT3IsotropicScordelisTests(StructrualMechanichsTestFactory):
    file_name = "Shell_T3_Isotropic_Scordelis_test/Shell_T3_Isotropic_Scordelis_test"


class ShellT3ThinBendingRollUpTests(StructrualMechanichsTestFactory):
    file_name = "Shell_T3_Thin__BendingRollUp_test/Shell_T3_Thin__BendingRollUp_test"


class ShellT3ThinDrillingRollUpTests(StructrualMechanichsTestFactory):
    file_name = "Shell_T3_Thin__DrillingRollUp_test/Shell_T3_Thin__DrillingRollUp_test"
