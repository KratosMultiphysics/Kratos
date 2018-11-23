from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
import structural_mechanics_analysis

# Other imports
import os

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


class StructuralMechanicsTestFactory(KratosUnittest.TestCase):
    def setUp(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):

            # Reading the ProjectParameters
            with open(self.file_name + "_parameters.json",'r') as parameter_file:
                ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())

            self.modify_parameters(ProjectParameters)

            # To avoid many prints
            if (ProjectParameters["problem_data"]["echo_level"].GetInt() == 0):
                KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

            # Creating the test
            model = KratosMultiphysics.Model()
            self.test = structural_mechanics_analysis.StructuralMechanicsAnalysis(model, ProjectParameters)
            self.test.Initialize()

    def modify_parameters(self, project_parameters):
        """This function can be used in derived classes to modify existing parameters
        before the execution of the test (e.g. switch to MPI)
        """
        pass

    def test_execution(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.test.RunSolutionLoop()

    def tearDown(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.test.Finalize()

class SimpleMeshMovingTest(StructuralMechanicsTestFactory):
    file_name = "mesh_moving_test/simple_mesh_moving_test"

class SDTwoDShearQuaPatchTest(StructuralMechanicsTestFactory):
    file_name = "patch_test/small_disp/patch_test_2D_shear_qua"

class SDTwoDShearTriPatchTest(StructuralMechanicsTestFactory):
    file_name = "patch_test/small_disp/patch_test_2D_shear_tri"

class SDTwoDTensionQuaPatchTest(StructuralMechanicsTestFactory):
    file_name = "patch_test/small_disp/patch_test_2D_tension_qua"

class SDTwoDTensionTriPatchTest(StructuralMechanicsTestFactory):
    file_name = "patch_test/small_disp/patch_test_2D_tension_tri"

class SDThreeDShearHexaPatchTest(StructuralMechanicsTestFactory):
    file_name = "patch_test/small_disp/patch_test_3D_shear_hexa"

class SDThreeDShearTetraPatchTest(StructuralMechanicsTestFactory):
    file_name = "patch_test/small_disp/patch_test_3D_shear_tetra"

class SDThreeDTensionHexaPatchTest(StructuralMechanicsTestFactory):
    file_name = "patch_test/small_disp/patch_test_3D_tension_hexa"

class SDThreeDTensionTetraPatchTest(StructuralMechanicsTestFactory):
    file_name = "patch_test/small_disp/patch_test_3D_tension_tetra"

class TLTwoDShearQuaPatchTest(StructuralMechanicsTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_2D_shear_qua"

class TLTwoDShearTriPatchTest(StructuralMechanicsTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_2D_shear_tri"

class TLTwoDTensionQuaPatchTest(StructuralMechanicsTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_2D_tension_qua"

class TLTwoDTensionTriPatchTest(StructuralMechanicsTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_2D_tension_tri"

class TLThreeDShearHexaPatchTest(StructuralMechanicsTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_3D_shear_hexa"

class TLThreeDShearTetraPatchTest(StructuralMechanicsTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_3D_shear_tetra"

class TLThreeDTensionHexaPatchTest(StructuralMechanicsTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_3D_tension_hexa"

class TLThreeDTensionTetraPatchTest(StructuralMechanicsTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_3D_tension_tetra"

class ULTwoDShearQuaPatchTest(StructuralMechanicsTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_2D_shear_qua"

class ULTwoDShearTriPatchTest(StructuralMechanicsTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_2D_shear_tri"

class ULTwoDTensionQuaPatchTest(StructuralMechanicsTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_2D_tension_qua"

class ULTwoDTensionTriPatchTest(StructuralMechanicsTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_2D_tension_tri"

class ULThreeDShearHexaPatchTest(StructuralMechanicsTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_3D_shear_hexa"

class ULThreeDShearTetraPatchTest(StructuralMechanicsTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_3D_shear_tetra"

class ULThreeDTensionHexaPatchTest(StructuralMechanicsTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_3D_tension_hexa"

class ULThreeDTensionTetraPatchTest(StructuralMechanicsTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_3D_tension_tetra"

class SprismMembranePatchTests(StructuralMechanicsTestFactory):
    file_name = "sprism_test/patch_membrane_test"

class SprismBendingPatchTests(StructuralMechanicsTestFactory):
    file_name = "sprism_test/patch_bending_test"

class EigenQ4Thick2x2PlateTests(StructuralMechanicsTestFactory):
    file_name = "eigen_test/Eigen_Q4_Thick_2x2_Plate_test"

class EigenTL3D8NCubeTests(StructuralMechanicsTestFactory):
    file_name = "eigen_test/Eigen_TL_3D8N_Cube_test"

class Eigen3D3NThinCircleTests(StructuralMechanicsTestFactory):
    file_name = "eigen_test/Eigen_3D3N_Thin_Circle_test"

class Fofi4PointTentnoCableTests(StructuralMechanicsTestFactory):
    file_name = "formfinding_test/Fofi_4Point_Tent_noCable_test"

class Fofi4PointTentCableTests(StructuralMechanicsTestFactory):
    file_name = "formfinding_test/Fofi_4Point_Tent_Cable_test"

class MembraneQ4PointLoadTests(StructuralMechanicsTestFactory):
    file_name = "membrane_test/Membrane_Q4_PointLoad_test"

class MembraneQ4TrussPointLoadTests(StructuralMechanicsTestFactory):
    file_name = "membrane_test/Membrane_Q4_Truss_PointLoad_test"

class Simple3D2NTrussTest(StructuralMechanicsTestFactory):
    file_name = "truss_test/nonlinear_3D2NTruss_test"

class Simple3D2NTrussLinearTest(StructuralMechanicsTestFactory):
    file_name = "truss_test/linear_3D2NTruss_test"

class Simple3D2NTrussDynamicTest(StructuralMechanicsTestFactory):
    file_name = "truss_test/dynamic_3D2NTruss_test"

class Simple3D2NTrussLinearCompressionPlasticTest(StructuralMechanicsTestFactory):
    file_name = "truss_test/linear_3D2NTruss_plastic_compression_test"

class Simple3D2NTrussLinearTensionPlasticTest(StructuralMechanicsTestFactory):
    file_name = "truss_test/linear_3D2NTruss_plastic_tension_test"

class Simple3D2NTrussNonLinearSnapthroughPlasticTest(StructuralMechanicsTestFactory):
    file_name = "truss_test/nonlinear_3D2NTruss_plastic_snapthrough_test"

class Simple3D2NTrussNonLinearTensionPlasticTest(StructuralMechanicsTestFactory):
    file_name = "truss_test/nonlinear_3D2NTruss_plastic_tension_test"

class Simple3D2NBeamCrTest(StructuralMechanicsTestFactory):
    file_name = "beam_test/nonlinear_3D2NBeamCr_test"

class Simple3D2NBeamCrLinearTest(StructuralMechanicsTestFactory):
    file_name = "beam_test/linear_3D2NBeamCr_test"

class Simple3D2NBeamCrDynamicTest(StructuralMechanicsTestFactory):
    file_name = "beam_test/dynamic_3D2NBeamCr_test"

class Simple2D2NBeamCrTest(StructuralMechanicsTestFactory):
    file_name = "beam_test/nonlinear_2D2NBeamCr_test"

class IsotropicDamageSimoJuPSTest(StructuralMechanicsTestFactory):
    file_name = "cl_test/IsotropicDamageSimoJu/PlaneStress_FourPointShear_test"

class SimpleSmallDeformationPlasticityMCTest(StructuralMechanicsTestFactory):
    file_name = "cl_test/SimpleSmallDeformationPlasticity/simple_small_deformation_plasticity_MC_test"

class SimpleSmallDeformationPlasticityVMTest(StructuralMechanicsTestFactory):
    file_name = "cl_test/SimpleSmallDeformationPlasticity/simple_small_deformation_plasticity_VM_test"

class SimpleSmallDeformationPlasticityDPTest(StructuralMechanicsTestFactory):
    file_name = "cl_test/SimpleSmallDeformationPlasticity/simple_small_deformation_plasticity_DP_test"

class SimpleSmallDeformationPlasticityTTest(StructuralMechanicsTestFactory):
    file_name = "cl_test/SimpleSmallDeformationPlasticity/simple_small_deformation_plasticity_T_test"

class SmallDeformationPlasticityTest(StructuralMechanicsTestFactory):
    file_name = "cl_test/SmallDeformationPlasticity/small_deformation_plasticity_test"

class ShellT3IsotropicLinearStaticStructScordelisLoRoofTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3_isotropic_linear_static_struct_scordelis_lo_roof"

class ShellT3AndQ4LinearStaticStructScordelisLoRoofTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_linear_static_struct_scordelis_lo_roof"

class ShellT3AndQ4LinearStaticStructPinchedCylinderTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_linear_static_struct_pinched_cylinder"

class ShellT3AndQ4LinearStaticStructPinchedHemisphereTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_linear_static_struct_pinched_hemisphere"

class ShellT3AndQ4LinearStaticStructClampedCylinderOrthotropicTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_linear_static_struct_clamped_cylinder_orthotropic"

class ShellT3AndQ4NonLinearStaticStructHingedCylRoofSnapthroughTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_nonlinear_static_struct_hinged_cyl_roof_snapthrough"

class ShellT3AndQ4NonLinearStaticStructHingedCylRoofSnapthroughOrthotropicTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_nonlinear_static_struct_hinged_cyl_roof_snapthrough_orthotropic"

class ShellT3AndQ4NonLinearDynamicStructOscillatingPlateTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_nonlinear_dynamic_struct_oscillating_plate"

class ShellT3AndQ4NonLinearDynamicStructOscillatingPlateLumpedTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_nonlinear_dynamic_struct_oscillating_plate_lumped"

class RigidFaceTestWithImposeRigidMovementProcess(StructuralMechanicsTestFactory):
    file_name = "rigid_test/rigid_test"

### OLD Tests Start, will be removed soon, Philipp Bucher, 31.01.2018 |---
class ShellQ4ThickBendingRollUpTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_Q4_Thick__BendingRollUp_test"
class ShellQ4ThickDrillingRollUpTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_Q4_Thick__DrillingRollUp_test"
class ShellQ4ThickOrthotropicLaminateLinearStaticTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_Q4_Thick_orthotropic_laminate_linear_static_test"
class ShellT3ThinBendingRollUpTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3_Thin__BendingRollUp_test"
class ShellT3ThinDrillingRollUpTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3_Thin__DrillingRollUp_test"
class ShellT3ThinOrthotropicLaminateLinearStaticTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3_Thin_orthotropic_laminate_linear_static_test"
class ShellT3IsotropicScordelisTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3_Isotropic_Scordelis_test"
class ShellT3ThickLinearStaticTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3_Thick_linear_static_test"
class ShellT3ThickNonLinearStaticTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3_Thick_nonlinear_static_test"
class ShellT3ThickLinearDynamicTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3_Thick_linear_dynamic_test"
class ShellT3ThickNonLinearDynamicTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3_Thick_nonlinear_dynamic_test"
class ShellT3ThickOrthotropicLaminateLinearStaticTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3_Thick_orthotropic_laminate_linear_static_test"
class ShellQ4ThinLinearStaticTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_Q4_Thin_linear_static_test"
class ShellQ4ThinNonLinearStaticTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_Q4_Thin_nonlinear_static_test"
class ShellQ4ThinLinearDynamicTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_Q4_Thin_linear_dynamic_test"
class ShellQ4ThinNonLinearDynamicTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_Q4_Thin_nonlinear_dynamic_test"
class ShellQ4ThinOrthotropicLaminateLinearStaticTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_Q4_Thin_orthotropic_laminate_linear_static_test"
### ---| OLD Tests End

class SprismPanTests(StructuralMechanicsTestFactory):
    file_name = "sprism_test/pan_test"

class PendulusTLTest(StructuralMechanicsTestFactory):
    file_name = "pendulus_test/pendulus_TL_test"

class PendulusULTest(StructuralMechanicsTestFactory):
    file_name = "pendulus_test/pendulus_UL_test"

class ShellT3AndQ4LinearStaticUnstructScordelisLoRoofTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_linear_static_unstruct_scordelis_lo_roof"

class ShellT3AndQ4LinearStaticUnstructUnstructPinchedCylinderTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_linear_static_unstruct_pinched_cylinder"

class ShellT3AndQ4LinearStaticUnstructPinchedHemisphereTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_linear_static_unstruct_pinched_hemisphere"

class ShellT3AndQ4LinearStaticUnstructClampedCylinderOrthotropicTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_linear_static_unstruct_clamped_cylinder_orthotropic"

class ShellT3AndQ4NonLinearStaticUnstructHingedCylRoofSnapthroughTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_nonlinear_static_unstruct_hinged_cyl_roof_snapthrough"

class ShellT3AndQ4NonLinearStaticUnstructHingedCylRoofSnapthroughOrthotropicTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_nonlinear_static_unstruct_hinged_cyl_roof_snapthrough_orthotropic"

class ShellT3AndQ4NonLinearDynamicUnstructOscillatingPlateTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_nonlinear_dynamic_unstruct_oscillating_plate"

class ShellT3AndQ4NonLinearDynamicUnstructOscillatingPlateLumpedTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_nonlinear_dynamic_unstruct_oscillating_plate_lumped"

class ShellT3AndQ4NonLinearDynamicStructPendulusTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_nonlinear_dynamic_struct_pendulus"

class ShellT3AndQ4NonLinearDynamicStructPendulusLumpedTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_nonlinear_dynamic_struct_pendulus_lumped"

class ShellT3AndQ4NonLinearDynamicUnstructPendulusTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_nonlinear_dynamic_unstruct_pendulus"

class ShellT3AndQ4NonLinearDynamicUnstructPendulusLumpedTests(StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_nonlinear_dynamic_unstruct_pendulus_lumped"

if __name__ == '__main__':
    KratosUnittest.main()
