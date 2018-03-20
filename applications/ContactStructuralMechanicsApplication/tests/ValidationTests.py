import os

# Import Kratos
from KratosMultiphysics import *

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
import Kratos_Execute_Contact_Structural_Test as Execute_Test

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

            # Checking if frictionless_by_components is defined
            try:
                self.frictionless_by_components
            except AttributeError:
                self.frictionless_by_components = False

            # Creating the model part
            self.test = Execute_Test.Kratos_Execute_Test(ProjectParameters, self.frictionless_by_components)

    def test_execution(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.test.Solve()

    def tearDown(self):
        pass

class ALMTaylorPatchDynamicTestContact(StructuralMechanichsTestFactory):
    file_name = "ALM_frictionless_contact_test_2D/taylor_patch_dynamic_test"

class ALMMeshMovingMatchingTestContact(StructuralMechanichsTestFactory):
    file_name = "ALM_frictionless_contact_test_2D/mesh_moving_matching_test"

class ALMMeshMovingNotMatchingTestContact(StructuralMechanichsTestFactory):
    file_name = "ALM_frictionless_contact_test_2D/mesh_moving_notmatching_test"

class ALMIroningTestContact(StructuralMechanichsTestFactory):
    file_name = "ALM_frictionless_contact_test_2D/ironing_test"

class ALMIroningDieTestContact(StructuralMechanichsTestFactory):
    file_name = "ALM_frictionless_contact_test_2D/ironing_die_test"

class LargeDisplacementPatchTestHexa(StructuralMechanichsTestFactory):
    file_name = "mesh_tying_test/3D_contact_patch_test_large_disp_hexa"

class ALMLargeDisplacementPatchTestTetra(StructuralMechanichsTestFactory):
    file_name = "ALM_frictionless_contact_test_3D/3D_contact_patch_test_large_disp_tetra"

class ALMLargeDisplacementPatchTestHexa(StructuralMechanichsTestFactory):
    file_name = "ALM_frictionless_contact_test_3D/3D_contact_patch_test_large_disp_hexa"
    
class ComponentsALMTaylorPatchDynamicTestContact(StructuralMechanichsTestFactory):
    frictionless_by_components = True
    file_name = "ALM_frictionless_contact_test_2D/taylor_patch_dynamic_test"

class ComponentsALMMeshMovingMatchingTestContact(StructuralMechanichsTestFactory):
    frictionless_by_components = True
    file_name = "ALM_frictionless_contact_test_2D/mesh_moving_matching_test"

class ComponentsALMMeshMovingNotMatchingTestContact(StructuralMechanichsTestFactory):
    frictionless_by_components = True
    file_name = "ALM_frictionless_contact_test_2D/mesh_moving_notmatching_test"

class ComponentsALMLargeDisplacementPatchTestTetra(StructuralMechanichsTestFactory):
    frictionless_by_components = True
    file_name = "ALM_frictionless_contact_test_3D/3D_contact_patch_test_large_disp_tetra"

class ComponentsALMLargeDisplacementPatchTestHexa(StructuralMechanichsTestFactory):
    frictionless_by_components = True
    file_name = "ALM_frictionless_contact_test_3D/3D_contact_patch_test_large_disp_hexa"
