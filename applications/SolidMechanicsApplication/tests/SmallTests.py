# Definition of the classes for the SMALL TESTS

# Import TestFactory
import TestFactory as TF

# Tests for elements:

# Small displacement elements SDE
class SD_Element2D4N_ShearTest(TF.TestFactory):
    file_name = "element_tests/small_displacement_elements/patch_test_2D4N_shear"
    file_parameters = "element_tests/shear_2D_parameters.json"
    
class SD_Element2D3N_ShearTest(TF.TestFactory):
    file_name = "element_tests/small_displacement_elements/patch_test_2D3N_shear"
    file_parameters = "element_tests/shear_2D_parameters.json"
    
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
    
# Large displacement shells SHE        
class Thick_Shell3D4N_BendingRollUpTest(TF.TestFactory):
    file_name = "element_tests/shell_elements/Shell_Q4_Thick__BendingRollUp"
    file_parameters = None

class Thick_Shell3D4N_DrillingRollUpTest(TF.TestFactory):
    file_name = "element_tests/shell_elements/Shell_Q4_Thick__DrillingRollUp"
    file_parameters = None

class Thin_Shell3D3N_BendingRollUpTest(TF.TestFactory):
    file_name = "element_tests/shell_elements/Shell_T3_Thin__BendingRollUp"
    file_parameters = None

class Thin_Shell3D3M_DrillingRollUpTest(TF.TestFactory):
    file_name = "element_tests/shell_elements/Shell_T3_Thin__DrillingRollUp"
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
