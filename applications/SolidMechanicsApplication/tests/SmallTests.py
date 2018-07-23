# Definition of the classes for the SMALL TESTS

# Import TestFactory
import TestFactory as TF

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Check external dependencies
try:
  import KratosMultiphysics
  import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication
  missing_external_dependencies = False
  missing_application = ''
except ImportError as e:
    missing_external_dependencies = True
    # extract name of the missing application from the error message
    import re
    missing_application = re.search(r'''.*'KratosMultiphysics\.(.*)'.*''','{0}'.format(e)).group(1)


# Tests for elements:

# Small displacement elements SDE
class SD_Element2D4N_ShearTest(TF.TestFactory):
    file_name = "element_tests/small_displacement_elements/patch_test_2D4N_shear"
    file_parameters = "element_tests/shear_2D_inputs.json"

class SD_Element2D3N_ShearTest(TF.TestFactory):
    file_name = "element_tests/small_displacement_elements/patch_test_2D3N_shear"
    file_parameters = "element_tests/shear_2D_parameters_material.json"

class SD_Element2D4N_TensionTest(TF.TestFactory):
    file_name = "element_tests/small_displacement_elements/patch_test_2D4N_tension"
    file_parameters = "element_tests/tension_2D_parameters.json"

class SD_Element2D3N_TensionTest(TF.TestFactory):
    file_name = "element_tests/small_displacement_elements/patch_test_2D3N_tension"
    file_parameters = "element_tests/tension_2D_parameters.json"

class SD_Element3D8N_ShearTest(TF.TestFactory):
    file_name = "element_tests/small_displacement_elements/patch_test_3D8N_shear"
    file_parameters = "element_tests/shear_3D_parameters.json"

class SD_Element3D4N_ShearTest(TF.TestFactory):
    file_name = "element_tests/small_displacement_elements/patch_test_3D4N_shear"
    file_parameters = "element_tests/shear_3D_parameters.json"

class SD_Element3D8N_TensionTest(TF.TestFactory):
    file_name = "element_tests/small_displacement_elements/patch_test_3D8N_tension"
    file_parameters = "element_tests/tension_3D_parameters.json"

class SD_Element3D4N_TensionTest(TF.TestFactory):
    file_name = "element_tests/small_displacement_elements/patch_test_3D4N_tension"
    file_parameters = "element_tests/tension_3D_parameters.json"


# Total lagrangian elements TLE
class TL_Element2D4N_ShearTest(TF.TestFactory):
    file_name = "element_tests/total_lagrangian_elements/patch_test_2D4N_shear"
    file_parameters = "element_tests/shear_2D_parameters.json"

class TL_Element2D3N_ShearTest(TF.TestFactory):
    file_name = "element_tests/total_lagrangian_elements/patch_test_2D3N_shear"
    file_parameters = "element_tests/shear_2D_parameters.json"

class TL_Element2D4N_TensionTest(TF.TestFactory):
    file_name = "element_tests/total_lagrangian_elements/patch_test_2D4N_tension"
    file_parameters = "element_tests/tension_2D_parameters.json"

class TL_Element2D3N_TensionTest(TF.TestFactory):
    file_name = "element_tests/total_lagrangian_elements/patch_test_2D3N_tension"
    file_parameters = "element_tests/tension_2D_parameters.json"

class TL_Element3D8N_ShearTest(TF.TestFactory):
    file_name = "element_tests/total_lagrangian_elements/patch_test_3D8N_shear"
    file_parameters = "element_tests/shear_3D_parameters.json"

class TL_Element3D4N_ShearTest(TF.TestFactory):
    file_name = "element_tests/total_lagrangian_elements/patch_test_3D4N_shear"
    file_parameters = "element_tests/shear_3D_parameters.json"

class TL_Element3D8N_TensionTest(TF.TestFactory):
    file_name = "element_tests/total_lagrangian_elements/patch_test_3D8N_tension"
    file_parameters = "element_tests/tension_3D_parameters.json"

class TL_Element3D4N_TensionTest(TF.TestFactory):
    file_name = "element_tests/total_lagrangian_elements/patch_test_3D4N_tension"
    file_parameters = "element_tests/tension_3D_parameters.json"

# Updated lagrangian elements ULE
class UL_Element2D4N_ShearTest(TF.TestFactory):
    file_name = "element_tests/updated_lagrangian_elements/patch_test_2D4N_shear"
    file_parameters = "element_tests/shear_2D_parameters.json"

class UL_Element2D3N_ShearTest(TF.TestFactory):
    file_name = "element_tests/updated_lagrangian_elements/patch_test_2D3N_shear"
    file_parameters = "element_tests/shear_2D_parameters.json"

class UL_Element2D4N_TensionTest(TF.TestFactory):
    file_name = "element_tests/updated_lagrangian_elements/patch_test_2D4N_tension"
    file_parameters = "element_tests/tension_2D_parameters.json"

class UL_Element2D3N_TensionTest(TF.TestFactory):
    file_name = "element_tests/updated_lagrangian_elements/patch_test_2D3N_tension"
    file_parameters = "element_tests/tension_2D_parameters.json"

class UL_Element3D8N_ShearTest(TF.TestFactory):
    file_name = "element_tests/updated_lagrangian_elements/patch_test_3D8N_shear"
    file_parameters = "element_tests/shear_3D_parameters.json"

class UL_Element3D4N_ShearTest(TF.TestFactory):
    file_name = "element_tests/updated_lagrangian_elements/patch_test_3D4N_shear"
    file_parameters = "element_tests/shear_3D_parameters.json"

class UL_Element3D8N_TensionTest(TF.TestFactory):
    file_name = "element_tests/updated_lagrangian_elements/patch_test_3D8N_tension"
    file_parameters = "element_tests/tension_3D_parameters.json"

class UL_Element3D4N_TensionTest(TF.TestFactory):
    file_name = "element_tests/updated_lagrangian_elements/patch_test_3D4N_tension"
    file_parameters = "element_tests/tension_3D_parameters.json"

# Large displacement dynamic tests DYN
class Dynamic_Test_Bossak_TL_3D(TF.TestFactory):
    file_name = "dynamic_tests/solid_elements/dynamic_bossak_TL3D"
    file_parameters = None

class Dynamic_Test_Simo_TL_3D(TF.TestFactory):
    file_name = "dynamic_tests/solid_elements/dynamic_bossak_TL3D"
    file_parameters  = "dynamic_tests/solid_elements/dynamic_simo_TL3D_input.json"


# Large displacement beam elements BEM
class LD_Beam_BendingRollUpTest(TF.TestFactory):
    file_name = "element_tests/beam_elements/static_beam_bending"
    file_parameters = "element_tests/beam_elements/beam_bending_parameters.json"

class EMC_Beam_BendingRollUpTest(TF.TestFactory):
    file_name = "element_tests/beam_elements/emc_static_beam_bending"
    file_parameters = "element_tests/beam_elements/beam_bending_parameters.json"

class LD_Beam_DynamicRotation(TF.TestFactory):
    file_name = "dynamic_tests/beam_elements/dynamic_beam"
    file_parameters = None

class EMC_Beam_DynamicRotation(TF.TestFactory):
    file_name = "dynamic_tests/beam_elements/emc_dynamic_beam"
    file_parameters = None

# Large displacement shells SHE
class Thick_Shell3D4N_BendingRollUpTest(TF.TestFactory):
    file_name = "element_tests/shell_elements/thick_quadrilateral_shell_bending"
    file_parameters = None

class Thick_Shell3D4N_DrillingRollUpTest(TF.TestFactory):
    file_name = "element_tests/shell_elements/thick_quadrilateral_shell_drilling"
    file_parameters = None

class Thin_Shell3D3N_BendingRollUpTest(TF.TestFactory):
    file_name = "element_tests/shell_elements/thin_triangular_shell_bending"
    file_parameters = None

class Thin_Shell3D3N_DrillingRollUpTest(TF.TestFactory):
    file_name = "element_tests/shell_elements/thin_triangular_shell_drilling"
    file_parameters = None

# Eigen modal analysis tests
class EigenQ4Thick2x2PlateTests(TF.TestFactory):
    file_name = "eigen_tests/Eigen_Q4_Thick_2x2_Plate"
    file_parameters = None

class EigenTL3D8NCubeTests(TF.TestFactory):
    file_name = "eigen_tests/Eigen_TL_3D8N_Cube"
    file_parameters = None

class Eigen3D3NThinCircleTests(TF.TestFactory):
    file_name = "eigen_tests/Eigen_3D3N_Thin_Circle"
    file_parameters = None


def SetTestSuite(suites):
    small_suite = suites['small']

    small_suite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            #SDE
            SD_Element2D4N_ShearTest,
            SD_Element2D3N_ShearTest,
            SD_Element2D4N_TensionTest,
            SD_Element2D3N_TensionTest,
            SD_Element3D8N_ShearTest,
            SD_Element3D4N_ShearTest,
            SD_Element3D8N_TensionTest,
            SD_Element3D4N_TensionTest,
            #TLE
            TL_Element2D4N_ShearTest,
            TL_Element2D3N_ShearTest,
            TL_Element2D4N_TensionTest,
            TL_Element2D3N_TensionTest,
            TL_Element3D8N_ShearTest,
            TL_Element3D4N_ShearTest,
            TL_Element3D8N_TensionTest,
            TL_Element3D4N_TensionTest,
            #ULE
            UL_Element2D4N_ShearTest,
            UL_Element2D3N_ShearTest,
            UL_Element2D4N_TensionTest,
            UL_Element2D3N_TensionTest,
            UL_Element3D8N_ShearTest,
            UL_Element3D4N_ShearTest,
            UL_Element3D8N_TensionTest,
            UL_Element3D4N_TensionTest,
            #DYN
            Dynamic_Test_Bossak_TL_3D,
            Dynamic_Test_Simo_TL_3D,
            #BEM
            LD_Beam_BendingRollUpTest,
            EMC_Beam_BendingRollUpTest,
            LD_Beam_DynamicRotation,
            EMC_Beam_DynamicRotation,
            #SHE
            Thick_Shell3D4N_BendingRollUpTest,
            Thick_Shell3D4N_DrillingRollUpTest,
            Thin_Shell3D3N_BendingRollUpTest,
            Thin_Shell3D3N_DrillingRollUpTest

        ])
    )

    if (missing_external_dependencies == False):
        if( hasattr(KratosMultiphysics.ExternalSolversApplication,  "FEASTSolver") ):
            small_suite.addTests(
                KratosUnittest.TestLoader().loadTestsFromTestCases([
                    EigenQ4Thick2x2PlateTests,
                    EigenTL3D8NCubeTests,
                    Eigen3D3NThinCircleTests
                ])
            )
        else:
            print("FEASTSolver solver is not included in the compilation of the External Solvers Application")

    return small_suite
