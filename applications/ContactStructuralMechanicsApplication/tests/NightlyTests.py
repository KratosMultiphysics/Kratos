import os

# Import Kratos
from KratosMultiphysics import *

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
import Kratos_Execute_Contact_Structural_Test as Execute_Test

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

class ALMMeshMovingMatchingTestContact(StructuralMechanichsTestFactory):
    file_name = "ALM_frictionless_contact_test_2D/mesh_moving_matching_test"
    
class ALMMeshMovingNotMatchingTestContact(StructuralMechanichsTestFactory):
    file_name = "ALM_frictionless_contact_test_2D/mesh_moving_notmatching_test"
    
class ALMTaylorPatchTestContact(StructuralMechanichsTestFactory):
    file_name = "ALM_frictionless_contact_test_2D/taylor_patch_test"

class ALMTaylorPatchDynamicTestContact(StructuralMechanichsTestFactory):
    file_name = "ALM_frictionless_contact_test_2D/taylor_patch_dynamic_test"
    
class ALMHertzSimpleSphereTestContact(StructuralMechanichsTestFactory):
    file_name = "ALM_frictionless_contact_test_2D/simple_hertz_sphere_plate_test"
    
class ALMHertzSimpleTestContact(StructuralMechanichsTestFactory):
    file_name = "ALM_frictionless_contact_test_2D/hertz_simple_test"
    
class ALMHertzSphereTestContact(StructuralMechanichsTestFactory):
    file_name = "ALM_frictionless_contact_test_2D/hertz_sphere_plate_test"
    
class ALMHertzCompleteTestContact(StructuralMechanichsTestFactory):
    file_name = "ALM_frictionless_contact_test_2D/hertz_complete_test"
