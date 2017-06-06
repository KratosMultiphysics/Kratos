# Definition of the classes for the SMALL TESTS

# Import TestFactory
import TestFactory as TF

# Tests for elements:

# Small displacement elements SDE
class SD_Element2D4N_ShearTest(TF.TestFactory):
    file_name = "element_tests/small_displacement_elements/patch_test_2D4N_shear"
    
class SD_Element2D3N_ShearTest(TF.TestFactory):
    file_name = "element_tests/small_displacement_elements/patch_test_2D3N_shear"
    
class SD_Element2D4N_TensionTest(TF.TestFactory):
    file_name = "element_tests/small_displacement_elements/patch_test_2D4N_tension"
    
class SD_Element2D3N_TensionTest(TF.TestFactory):
    file_name = "element_tests/small_displacement_elements/patch_test_2D3N_tension"
    
class SD_Element3D8N_ShearTest(TF.TestFactory):
    file_name = "element_tests/small_displacement_elements/patch_test_3D8N_shear"
    
class SD_Element3D4N_ShearTest(TF.TestFactory):
    file_name = "element_tests/small_displacement_elements/patch_test_3D4N_shear"
    
class SD_Element3D8N_TensionTest(TF.TestFactory):
    file_name = "element_tests/small_displacement_elements/patch_test_3D8N_tension"
    
class SD_Element3D4N_TensionTest(TF.TestFactory):
    file_name = "element_tests/small_displacement_elements/patch_test_3D4N_tension"
    

# Total lagrangian elements TLE
class TL_Element2D4N_ShearTest(TF.TestFactory):
    file_name = "element_tests/total_lagrangian_elements/patch_test_2D4N_shear"
    
class TL_Element2D3N_ShearTest(TF.TestFactory):
    file_name = "element_tests/total_lagrangian_elements/patch_test_2D3N_shear"
    
class TL_Element2D4N_TensionTest(TF.TestFactory):
    file_name = "element_tests/total_lagrangian_elements/patch_test_2D4N_tension"
    
class TL_Element2D3N_TensionTest(TF.TestFactory):
    file_name = "element_tests/total_lagrangian_elements/patch_test_2D3N_tension"
    
class TL_Element3D8N_ShearTest(TF.TestFactory):
    file_name = "element_tests/total_lagrangian_elements/patch_test_3D8N_shear"
    
class TL_Element3D4N_ShearTest(TF.TestFactory):
    file_name = "element_tests/total_lagrangian_elements/patch_test_3D4N_shear"
    
class TL_Element3D8N_TensionTest(TF.TestFactory):
    file_name = "element_tests/total_lagrangian_elements/patch_test_3D8N_tension"
    
class TL_Element3D4N_TensionTest(TF.TestFactory):
    file_name = "element_tests/total_lagrangian_elements/patch_test_3D4N_tension"


# Updated lagrangian elements ULE
class UL_Element2D4N_ShearTest(TF.TestFactory):
    file_name = "element_tests/updated_lagrangian_elements/patch_test_2D4N_shear"
    
class UL_Element2D3N_ShearTest(TF.TestFactory):
    file_name = "element_tests/updated_lagrangian_elements/patch_test_2D3N_shear"
    
class UL_Element2D4N_TensionTest(TF.TestFactory):
    file_name = "element_tests/updated_lagrangian_elements/patch_test_2D4N_tension"
    
class UL_Element2D3N_TensionTest(TF.TestFactory):
    file_name = "element_tests/updated_lagrangian_elements/patch_test_2D3N_tension"
    
class UL_Element3D8N_ShearTest(TF.TestFactory):
    file_name = "element_tests/updated_lagrangian_elements/patch_test_3D8N_shear"
    
class UL_Element3D4N_ShearTest(TF.TestFactory):
    file_name = "element_tests/updated_lagrangian_elements/patch_test_3D4N_shear"
    
class UL_Element3D8N_TensionTest(TF.TestFactory):
    file_name = "element_tests/updated_lagrangian_elements/patch_test_3D8N_tension"
    
class UL_Element3D4N_TensionTest(TF.TestFactory):
    file_name = "element_tests/updated_lagrangian_elements/patch_test_3D4N_tension"

    
# Large displacement shells SHE        
class Thick_Shell3D4N_BendingRollUpTest(TF.TestFactory):
    file_name = "element_tests/shell_elements/Shell_Q4_Thick__BendingRollUp_test"


class Thick_Shell3D4N_DrillingRollUpTest(TF.TestFactory):
    file_name = "element_tests/shell_elements/Shell_Q4_Thick__DrillingRollUp_test"


class Thin_Shell3D3N_BendingRollUpTest(TF.TestFactory):
    file_name = "element_tests/shell_elements/Shell_T3_Thin__BendingRollUp_test"


class Thin_Shell3D3M_DrillingRollUpTest(TF.TestFactory):
    file_name = "element_tests/shell_elements/Shell_T3_Thin__DrillingRollUp_test"


