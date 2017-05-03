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

class SolidMechanichsTestFactory(KratosUnittest.TestCase):

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

class SDTwoDShearQuaPatchTest(SolidMechanichsTestFactory):
    file_name = "patch_test/small_disp/patch_test_2D_shear_qua"
    
class SDTwoDShearTriPatchTest(SolidMechanichsTestFactory):
    file_name = "patch_test/small_disp/patch_test_2D_shear_tri"
    
class SDTwoDTensionQuaPatchTest(SolidMechanichsTestFactory):
    file_name = "patch_test/small_disp/patch_test_2D_tension_qua"
    
class SDTwoDTensionTriPatchTest(SolidMechanichsTestFactory):
    file_name = "patch_test/small_disp/patch_test_2D_tension_tri"
    
class SDThreeDShearHexaPatchTest(SolidMechanichsTestFactory):
    file_name = "patch_test/small_disp/patch_test_3D_shear_hexa"
    
class SDThreeDShearTetraPatchTest(SolidMechanichsTestFactory):
    file_name = "patch_test/small_disp/patch_test_3D_shear_tetra"
    
class SDThreeDTensionHexaPatchTest(SolidMechanichsTestFactory):
    file_name = "patch_test/small_disp/patch_test_3D_tension_hexa"
    
class SDThreeDTensionTetraPatchTest(SolidMechanichsTestFactory):
    file_name = "patch_test/small_disp/patch_test_3D_tension_tetra"
    
class TLTwoDShearQuaPatchTest(SolidMechanichsTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_2D_shear_qua"
    
class TLTwoDShearTriPatchTest(SolidMechanichsTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_2D_shear_tri"
    
class TLTwoDTensionQuaPatchTest(SolidMechanichsTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_2D_tension_qua"
    
class TLTwoDTensionTriPatchTest(SolidMechanichsTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_2D_tension_tri"
    
class TLThreeDShearHexaPatchTest(SolidMechanichsTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_3D_shear_hexa"
    
class TLThreeDShearTetraPatchTest(SolidMechanichsTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_3D_shear_tetra"
    
class TLThreeDTensionHexaPatchTest(SolidMechanichsTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_3D_tension_hexa"
    
class TLThreeDTensionTetraPatchTest(SolidMechanichsTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_3D_tension_tetra"
    
class ULTwoDShearQuaPatchTest(SolidMechanichsTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_2D_shear_qua"
    
class ULTwoDShearTriPatchTest(SolidMechanichsTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_2D_shear_tri"
    
class ULTwoDTensionQuaPatchTest(SolidMechanichsTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_2D_tension_qua"
    
class ULTwoDTensionTriPatchTest(SolidMechanichsTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_2D_tension_tri"
    
class ULThreeDShearHexaPatchTest(SolidMechanichsTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_3D_shear_hexa"
    
class ULThreeDShearTetraPatchTest(SolidMechanichsTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_3D_shear_tetra"
    
class ULThreeDTensionHexaPatchTest(SolidMechanichsTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_3D_tension_hexa"
    
class ULThreeDTensionTetraPatchTest(SolidMechanichsTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_3D_tension_tetra"
