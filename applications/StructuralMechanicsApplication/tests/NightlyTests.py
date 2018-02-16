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

class SDTwoDShearQuaPatchTest(StructuralMechanichsTestFactory):
    file_name = "patch_test/small_disp/patch_test_2D_shear_qua"
    
class SDTwoDShearTriPatchTest(StructuralMechanichsTestFactory):
    file_name = "patch_test/small_disp/patch_test_2D_shear_tri"
    
class SDTwoDTensionQuaPatchTest(StructuralMechanichsTestFactory):
    file_name = "patch_test/small_disp/patch_test_2D_tension_qua"
    
class SDTwoDTensionTriPatchTest(StructuralMechanichsTestFactory):
    file_name = "patch_test/small_disp/patch_test_2D_tension_tri"
    
class SDThreeDShearHexaPatchTest(StructuralMechanichsTestFactory):
    file_name = "patch_test/small_disp/patch_test_3D_shear_hexa"
    
class SDThreeDShearTetraPatchTest(StructuralMechanichsTestFactory):
    file_name = "patch_test/small_disp/patch_test_3D_shear_tetra"
    
class SDThreeDTensionHexaPatchTest(StructuralMechanichsTestFactory):
    file_name = "patch_test/small_disp/patch_test_3D_tension_hexa"
    
class SDThreeDTensionTetraPatchTest(StructuralMechanichsTestFactory):
    file_name = "patch_test/small_disp/patch_test_3D_tension_tetra"
    
class TLTwoDShearQuaPatchTest(StructuralMechanichsTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_2D_shear_qua"
    
class TLTwoDShearTriPatchTest(StructuralMechanichsTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_2D_shear_tri"
    
class TLTwoDTensionQuaPatchTest(StructuralMechanichsTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_2D_tension_qua"
    
class TLTwoDTensionTriPatchTest(StructuralMechanichsTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_2D_tension_tri"
    
class TLThreeDShearHexaPatchTest(StructuralMechanichsTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_3D_shear_hexa"
    
class TLThreeDShearTetraPatchTest(StructuralMechanichsTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_3D_shear_tetra"
    
class TLThreeDTensionHexaPatchTest(StructuralMechanichsTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_3D_tension_hexa"
    
class TLThreeDTensionTetraPatchTest(StructuralMechanichsTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_3D_tension_tetra"
    
class ULTwoDShearQuaPatchTest(StructuralMechanichsTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_2D_shear_qua"
    
class ULTwoDShearTriPatchTest(StructuralMechanichsTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_2D_shear_tri"
    
class ULTwoDTensionQuaPatchTest(StructuralMechanichsTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_2D_tension_qua"
    
class ULTwoDTensionTriPatchTest(StructuralMechanichsTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_2D_tension_tri"
    
class ULThreeDShearHexaPatchTest(StructuralMechanichsTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_3D_shear_hexa"
    
class ULThreeDShearTetraPatchTest(StructuralMechanichsTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_3D_shear_tetra"
    
class ULThreeDTensionHexaPatchTest(StructuralMechanichsTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_3D_tension_hexa"
    
class ULThreeDTensionTetraPatchTest(StructuralMechanichsTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_3D_tension_tetra"

class SprismMembranePatchTests(StructuralMechanichsTestFactory):
    file_name = "sprism_test/patch_membrane_test"

class SprismBendingPatchTests(StructuralMechanichsTestFactory):
    file_name = "sprism_test/patch_bending_test"

class EigenQ4Thick2x2PlateTests(StructuralMechanichsTestFactory):
    file_name = "eigen_test/Eigen_Q4_Thick_2x2_Plate_test"

class EigenTL3D8NCubeTests(StructuralMechanichsTestFactory):
    file_name = "eigen_test/Eigen_TL_3D8N_Cube_test"
    
class Eigen3D3NThinCircleTests(StructuralMechanichsTestFactory):
    file_name = "eigen_test/Eigen_3D3N_Thin_Circle_test"
    
class Fofi4PointTentnoCableTests(StructuralMechanichsTestFactory):
    file_name = "formfinding_test/Fofi_4Point_Tent_noCable_test"
	
class Fofi4PointTentCableTests(StructuralMechanichsTestFactory):
    file_name = "formfinding_test/Fofi_4Point_Tent_Cable_test"	
    
class MembraneQ4PointLoadTests(StructuralMechanichsTestFactory):
    file_name = "membrane_test/Membrane_Q4_PointLoad_test"
    
class MembraneQ4TrussPointLoadTests(StructuralMechanichsTestFactory):
    file_name = "membrane_test/Membrane_Q4_Truss_PointLoad_test"  
  
class Simple3D2NTrussTest(StructuralMechanichsTestFactory):
    file_name = "truss_test/nonlinear_3D2NTruss_test"
    
class Simple3D2NTrussLinearTest(StructuralMechanichsTestFactory):
    file_name = "truss_test/linear_3D2NTruss_test"

class Simple3D2NTrussDynamicTest(StructuralMechanichsTestFactory):
    file_name = "truss_test/dynamic_3D2NTruss_test"

class Simple3D2NBeamCrTest(StructuralMechanichsTestFactory):
    file_name = "beam_test/nonlinear_3D2NBeamCr_test"
     
class Simple3D2NBeamCrLinearTest(StructuralMechanichsTestFactory):
    file_name = "beam_test/linear_3D2NBeamCr_test"
    
class Simple3D2NBeamCrDynamicTest(StructuralMechanichsTestFactory):
    file_name = "beam_test/dynamic_3D2NBeamCr_test"  

class Simple2D2NBeamCrTest(StructuralMechanichsTestFactory):
    file_name = "beam_test/nonlinear_2D2NBeamCr_test"
    
class IsotropicDamageSimoJuPSTest(StructuralMechanichsTestFactory):
    file_name = "cl_test/IsotropicDamageSimoJu/PlaneStress_FourPointShear_test"

class ShellT3IsotropicLinearStaticStructScordelisLoRoofTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_isotropic_linear_static_struct_scordelis_lo_roof"

class ShellT3AndQ4LinearStaticStructScordelisLoRoofTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_linear_static_struct_scordelis_lo_roof"

class ShellT3AndQ4LinearStaticStructPinchedCylinderTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_linear_static_struct_pinched_cylinder"

class ShellT3AndQ4LinearStaticStructPinchedHemisphereTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_linear_static_struct_pinched_hemisphere"

class ShellT3AndQ4LinearStaticStructClampedCylinderOrthotropicTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_linear_static_struct_clamped_cylinder_orthotropic"

class ShellT3AndQ4NonLinearStaticStructHingedCylRoofSnapthroughTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_nonlinear_static_struct_hinged_cyl_roof_snapthrough"
    
class ShellT3AndQ4NonLinearStaticStructHingedCylRoofSnapthroughOrthotropicTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_nonlinear_static_struct_hinged_cyl_roof_snapthrough_orthotropic"

class ShellT3AndQ4NonLinearDynamicStructOscillatingPlateTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_nonlinear_dynamic_struct_oscillating_plate"

class ShellT3AndQ4NonLinearDynamicStructOscillatingPlateLumpedTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_nonlinear_dynamic_struct_oscillating_plate_lumped"

### OLD Tests Start, will be removed soon, Philipp Bucher, 31.01.2018 |---
class ShellQ4ThickBendingRollUpTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_Thick__BendingRollUp_test"
class ShellQ4ThickDrillingRollUpTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_Thick__DrillingRollUp_test"
class ShellQ4ThickOrthotropicLaminateLinearStaticTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_Thick_orthotropic_laminate_linear_static_test"
class ShellT3ThinBendingRollUpTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_Thin__BendingRollUp_test"
class ShellT3ThinDrillingRollUpTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_Thin__DrillingRollUp_test"
class ShellT3ThinOrthotropicLaminateLinearStaticTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_Thin_orthotropic_laminate_linear_static_test"
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
class ShellT3ThickOrthotropicLaminateLinearStaticTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_Thick_orthotropic_laminate_linear_static_test"
class ShellQ4ThinLinearStaticTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_Thin_linear_static_test"
class ShellQ4ThinNonLinearStaticTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_Thin_nonlinear_static_test"
class ShellQ4ThinLinearDynamicTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_Thin_linear_dynamic_test"
class ShellQ4ThinNonLinearDynamicTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_Thin_nonlinear_dynamic_test"
class ShellQ4ThinOrthotropicLaminateLinearStaticTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_Thin_orthotropic_laminate_linear_static_test"
### ---| OLD Tests End
