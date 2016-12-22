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
    
class BasicCATest(StructuralMechanichsTestFactory):
    file_name = "CA_test/basic_CA_test"
    
class SolidCATest(StructuralMechanichsTestFactory):
    file_name = "CA_test/solid_CA_test"
    
class HyperSimplePatchTestContact(StructuralMechanichsTestFactory):
    file_name = "contact_test_2D/hyper_simple_patch_test"
    
class SimplePatchTestContact(StructuralMechanichsTestFactory):
    file_name = "contact_test_2D/simple_patch_test"
    
class SimpleSlopePatchTestContact(StructuralMechanichsTestFactory):
    file_name = "contact_test_2D/simple_slope_patch_test"
    
class SimplePatchNotMatchingATestContact(StructuralMechanichsTestFactory):
    file_name = "contact_test_2D/simple_patch_notmatching_a_test"
    
class SimplePatchNotMatchingBTestContact(StructuralMechanichsTestFactory):
    file_name = "contact_test_2D/simple_patch_notmatching_b_test"
    
class TaylorPatchTestContact(StructuralMechanichsTestFactory):
    file_name = "contact_test_2D/taylor_patch_test"

class TaylorPatchDynamicTestContact(StructuralMechanichsTestFactory):
    file_name = "contact_test_2D/taylor_patch_dynamic_test"
    
class HertzSimpleSphereTestContact(StructuralMechanichsTestFactory):
    file_name = "contact_test_2D/simple_hertz_sphere_plate_test"
    
class HertzSimpleTestContact(StructuralMechanichsTestFactory):
    file_name = "contact_test_2D/hertz_simple_test"
    
class HertzSphereTestContact(StructuralMechanichsTestFactory):
    file_name = "contact_test_2D/hertz_sphere_plate_test"
    
class HertzCompleteTestContact(StructuralMechanichsTestFactory):
    file_name = "contact_test_2D/hertz_complete_test"

class ThreeDSimplestPatchMatchingTestContact(StructuralMechanichsTestFactory):
    file_name = "contact_test_3D/3D_contact_simplest_patch_matching_test"
    
class ThreeDSimplestTrianglePatchMatchingTestContact(StructuralMechanichsTestFactory):
    file_name = "contact_test_3D/3D_contact_simplest_triangle_patch_matching_test"
    
class ThreeDPatchMatchingTestContact(StructuralMechanichsTestFactory):
    file_name = "contact_test_3D/3D_contact_patch_matching_test"

class ThreeDPatchNotMatchingTestContact(StructuralMechanichsTestFactory):
    file_name = "contact_test_3D/3D_contact_patch_nonmatching_test"
    
