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
    
class SimplePatchTestTwoDMeshTying(StructuralMechanichsTestFactory):
    file_name = "mesh_tying_test/simple_patch_test_2D"
    
class SimpleSlopePatchTestTwoDMeshTying(StructuralMechanichsTestFactory):
    file_name = "mesh_tying_test/hyper_simple_slope_patch_test_2D"
    
class SimplestPatchTestThreeDMeshTying(StructuralMechanichsTestFactory):
    file_name = "mesh_tying_test/3D_contact_simplest_patch_matching_test"
    
class SimplePatchTestThreeDMeshTying(StructuralMechanichsTestFactory):
    file_name = "mesh_tying_test/simple_patch_test_3D"
    
class ALMHyperSimplePatchTestContact(StructuralMechanichsTestFactory):
    file_name = "ALM_frictionless_contact_test_2D/hyper_simple_patch_test"
    
class ALMHyperSimpleSlopePatchTestContact(StructuralMechanichsTestFactory):
    file_name = "ALM_frictionless_contact_test_2D/hyper_simple_slope_patch_test"
    
class ALMTwoDPatchComplexGeomTestContact(StructuralMechanichsTestFactory):
    file_name = "ALM_frictionless_contact_test_2D/contact_patch_complex_geom_test"
    
class ALMTwoDPatchComplexGeomSlopeTestContact(StructuralMechanichsTestFactory):
    file_name = "ALM_frictionless_contact_test_2D/contact_patch_complex_geom_slope_test"
    
class ALMSimplePatchTestContact(StructuralMechanichsTestFactory):
    file_name = "ALM_frictionless_contact_test_2D/simple_patch_test"
    
class ALMSimpleSlopePatchTestContact(StructuralMechanichsTestFactory):
    file_name = "ALM_frictionless_contact_test_2D/hyper_simple_slope_patch_test"
    
class ALMSimplePatchNotMatchingATestContact(StructuralMechanichsTestFactory):
    file_name = "ALM_frictionless_contact_test_2D/simple_patch_notmatching_a_test"
    
class ALMSimplePatchNotMatchingBTestContact(StructuralMechanichsTestFactory):
    file_name = "ALM_frictionless_contact_test_2D/simple_patch_notmatching_b_test"

class ALMThreeDSimplestPatchMatchingTestContact(StructuralMechanichsTestFactory):
    file_name = "ALM_frictionless_contact_test_3D/3D_contact_simplest_patch_matching_test"
    
class ALMThreeDSimplestPatchMatchingAdaptativeTestContact(StructuralMechanichsTestFactory):
    file_name = "ALM_frictionless_contact_test_3D/3D_contact_simplest_patch_matching_adaptative_test"
    
class ALMThreeDPatchComplexGeomTestContact(StructuralMechanichsTestFactory):
    file_name = "ALM_frictionless_contact_test_3D/3D_contact_patch_complex_geom_test"
    
class ALMThreeDSimplestPatchMatchingSlopeTestContact(StructuralMechanichsTestFactory):
    file_name = "ALM_frictionless_contact_test_3D/3D_contact_simplest_patch_matching_slope_test"
    
class ALMThreeDPatchMatchingTestContact(StructuralMechanichsTestFactory):
    file_name = "ALM_frictionless_contact_test_3D/3D_contact_patch_matching_test"

class ALMThreeDPatchNotMatchingTestContact(StructuralMechanichsTestFactory):
    file_name = "ALM_frictionless_contact_test_3D/3D_contact_patch_nonmatching_test"
