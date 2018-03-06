import structural_mechanics_test_factory

class SDTwoDShearQuaPatchTest(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "patch_test/small_disp/patch_test_2D_shear_qua"

class SDTwoDShearTriPatchTest(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "patch_test/small_disp/patch_test_2D_shear_tri"

class SDTwoDTensionQuaPatchTest(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "patch_test/small_disp/patch_test_2D_tension_qua"

class SDTwoDTensionTriPatchTest(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "patch_test/small_disp/patch_test_2D_tension_tri"

class SDThreeDShearHexaPatchTest(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "patch_test/small_disp/patch_test_3D_shear_hexa"

class SDThreeDShearTetraPatchTest(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "patch_test/small_disp/patch_test_3D_shear_tetra"

class SDThreeDTensionHexaPatchTest(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "patch_test/small_disp/patch_test_3D_tension_hexa"

class SDThreeDTensionTetraPatchTest(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "patch_test/small_disp/patch_test_3D_tension_tetra"

class TLTwoDShearQuaPatchTest(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_2D_shear_qua"

class TLTwoDShearTriPatchTest(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_2D_shear_tri"

class TLTwoDTensionQuaPatchTest(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_2D_tension_qua"

class TLTwoDTensionTriPatchTest(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_2D_tension_tri"

class TLThreeDShearHexaPatchTest(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_3D_shear_hexa"

class TLThreeDShearTetraPatchTest(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_3D_shear_tetra"

class TLThreeDTensionHexaPatchTest(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_3D_tension_hexa"

class TLThreeDTensionTetraPatchTest(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "patch_test/total_lagrangian/patch_test_3D_tension_tetra"

class ULTwoDShearQuaPatchTest(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_2D_shear_qua"

class ULTwoDShearTriPatchTest(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_2D_shear_tri"

class ULTwoDTensionQuaPatchTest(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_2D_tension_qua"

class ULTwoDTensionTriPatchTest(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_2D_tension_tri"

class ULThreeDShearHexaPatchTest(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_3D_shear_hexa"

class ULThreeDShearTetraPatchTest(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_3D_shear_tetra"

class ULThreeDTensionHexaPatchTest(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_3D_tension_hexa"

class ULThreeDTensionTetraPatchTest(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "patch_test/updated_lagrangian/patch_test_3D_tension_tetra"

class SprismMembranePatchTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "sprism_test/patch_membrane_test"

class SprismBendingPatchTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "sprism_test/patch_bending_test"

class EigenQ4Thick2x2PlateTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "eigen_test/Eigen_Q4_Thick_2x2_Plate_test"

class EigenTL3D8NCubeTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "eigen_test/Eigen_TL_3D8N_Cube_test"

class Eigen3D3NThinCircleTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "eigen_test/Eigen_3D3N_Thin_Circle_test"

class Fofi4PointTentnoCableTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "formfinding_test/Fofi_4Point_Tent_noCable_test"

class Fofi4PointTentCableTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "formfinding_test/Fofi_4Point_Tent_Cable_test"

class MembraneQ4PointLoadTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "membrane_test/Membrane_Q4_PointLoad_test"

class MembraneQ4TrussPointLoadTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "membrane_test/Membrane_Q4_Truss_PointLoad_test"

class Simple3D2NTrussTest(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "truss_test/nonlinear_3D2NTruss_test"

class Simple3D2NTrussLinearTest(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "truss_test/linear_3D2NTruss_test"

class Simple3D2NTrussDynamicTest(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "truss_test/dynamic_3D2NTruss_test"

class Simple3D2NBeamCrTest(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "beam_test/nonlinear_3D2NBeamCr_test"

class Simple3D2NBeamCrLinearTest(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "beam_test/linear_3D2NBeamCr_test"

class Simple3D2NBeamCrDynamicTest(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "beam_test/dynamic_3D2NBeamCr_test"

class Simple2D2NBeamCrTest(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "beam_test/nonlinear_2D2NBeamCr_test"

class IsotropicDamageSimoJuPSTest(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "cl_test/IsotropicDamageSimoJu/PlaneStress_FourPointShear_test"

class ShellT3IsotropicLinearStaticStructScordelisLoRoofTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3_isotropic_linear_static_struct_scordelis_lo_roof"

class ShellT3AndQ4LinearStaticStructScordelisLoRoofTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_linear_static_struct_scordelis_lo_roof"

class ShellT3AndQ4LinearStaticStructPinchedCylinderTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_linear_static_struct_pinched_cylinder"

class ShellT3AndQ4LinearStaticStructPinchedHemisphereTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_linear_static_struct_pinched_hemisphere"

class ShellT3AndQ4LinearStaticStructClampedCylinderOrthotropicTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_linear_static_struct_clamped_cylinder_orthotropic"

class ShellT3AndQ4NonLinearStaticStructHingedCylRoofSnapthroughTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_nonlinear_static_struct_hinged_cyl_roof_snapthrough"

class ShellT3AndQ4NonLinearStaticStructHingedCylRoofSnapthroughOrthotropicTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_nonlinear_static_struct_hinged_cyl_roof_snapthrough_orthotropic"

class ShellT3AndQ4NonLinearDynamicStructOscillatingPlateTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_nonlinear_dynamic_struct_oscillating_plate"

class ShellT3AndQ4NonLinearDynamicStructOscillatingPlateLumpedTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3andQ4_nonlinear_dynamic_struct_oscillating_plate_lumped"

### OLD Tests Start, will be removed soon, Philipp Bucher, 31.01.2018 |---
class ShellQ4ThickBendingRollUpTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_Q4_Thick__BendingRollUp_test"
class ShellQ4ThickDrillingRollUpTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_Q4_Thick__DrillingRollUp_test"
class ShellQ4ThickOrthotropicLaminateLinearStaticTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_Q4_Thick_orthotropic_laminate_linear_static_test"
class ShellT3ThinBendingRollUpTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3_Thin__BendingRollUp_test"
class ShellT3ThinDrillingRollUpTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3_Thin__DrillingRollUp_test"
class ShellT3ThinOrthotropicLaminateLinearStaticTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3_Thin_orthotropic_laminate_linear_static_test"
class ShellT3IsotropicScordelisTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3_Isotropic_Scordelis_test"
class ShellT3ThickLinearStaticTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3_Thick_linear_static_test"
class ShellT3ThickNonLinearStaticTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3_Thick_nonlinear_static_test"
class ShellT3ThickLinearDynamicTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3_Thick_linear_dynamic_test"
class ShellT3ThickNonLinearDynamicTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3_Thick_nonlinear_dynamic_test"
class ShellT3ThickOrthotropicLaminateLinearStaticTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_T3_Thick_orthotropic_laminate_linear_static_test"
class ShellQ4ThinLinearStaticTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_Q4_Thin_linear_static_test"
class ShellQ4ThinNonLinearStaticTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_Q4_Thin_nonlinear_static_test"
class ShellQ4ThinLinearDynamicTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_Q4_Thin_linear_dynamic_test"
class ShellQ4ThinNonLinearDynamicTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_Q4_Thin_nonlinear_dynamic_test"
class ShellQ4ThinOrthotropicLaminateLinearStaticTests(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "shell_test/Shell_Q4_Thin_orthotropic_laminate_linear_static_test"
### ---| OLD Tests End
