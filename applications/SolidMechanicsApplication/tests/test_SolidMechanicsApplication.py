# import Kratos
import KratosMultiphysics 
import KratosMultiphysics.ExternalSolversApplication 
import KratosMultiphysics.SolidMechanicsApplication 

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

try:
  import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication
  missing_external_dependencies = False
  missing_application = ''
except ImportError as e:
    missing_external_dependencies = True
    # extract name of the missing application from the error message
    import re
    missing_application = re.search(r'''.*'KratosMultiphysics\.(.*)'.*''','{0}'.format(e)).group(1)

# Import the tests o test_classes to create the suits:

## SMALL TESTS

# Tests for elements:

# Small displacement elements SDE
from SmallTests import SD_Element2D4N_ShearTest     as SDE_2D4N_S_TEST
from SmallTests import SD_Element2D3N_ShearTest     as SDE_2D3N_S_TEST
from SmallTests import SD_Element2D4N_TensionTest   as SDE_2D4N_T_TEST
from SmallTests import SD_Element2D3N_TensionTest   as SDE_2D3N_T_TEST
from SmallTests import SD_Element3D8N_ShearTest     as SDE_3D8N_S_TEST
from SmallTests import SD_Element3D4N_ShearTest     as SDE_3D4N_S_TEST
from SmallTests import SD_Element3D8N_TensionTest   as SDE_3D8N_T_TEST
from SmallTests import SD_Element3D4N_TensionTest   as SDE_3D4N_T_TEST
# Total lagrangian elements TLE
from SmallTests import TL_Element2D4N_ShearTest     as TLE_2D4N_S_TEST
from SmallTests import TL_Element2D3N_ShearTest     as TLE_2D3N_S_TEST
from SmallTests import TL_Element2D4N_TensionTest   as TLE_2D4N_T_TEST
from SmallTests import TL_Element2D3N_TensionTest   as TLE_2D3N_T_TEST
from SmallTests import TL_Element3D8N_ShearTest     as TLE_3D8N_S_TEST
from SmallTests import TL_Element3D4N_ShearTest     as TLE_3D4N_S_TEST
from SmallTests import TL_Element3D8N_TensionTest   as TLE_3D8N_T_TEST
from SmallTests import TL_Element3D4N_TensionTest   as TLE_3D4N_T_TEST
# Updated lagrangian elements ULE
from SmallTests import UL_Element2D4N_ShearTest     as ULE_2D4N_S_TEST
from SmallTests import UL_Element2D3N_ShearTest     as ULE_2D3N_S_TEST
from SmallTests import UL_Element2D4N_TensionTest   as ULE_2D4N_T_TEST
from SmallTests import UL_Element2D3N_TensionTest   as ULE_2D3N_T_TEST
from SmallTests import UL_Element3D8N_ShearTest     as ULE_3D8N_S_TEST
from SmallTests import UL_Element3D4N_ShearTest     as ULE_3D4N_S_TEST
from SmallTests import UL_Element3D8N_TensionTest   as ULE_3D8N_T_TEST
from SmallTests import UL_Element3D4N_TensionTest   as ULE_3D4N_T_TEST
# Large displacement shells SHE
from SmallTests import Thick_Shell3D4N_BendingRollUpTest  as SHE_3D4N_B_TEST
from SmallTests import Thick_Shell3D4N_DrillingRollUpTest as SHE_3D4N_D_TEST
from SmallTests import Thin_Shell3D3N_BendingRollUpTest   as SHE_3D3N_B_TEST
from SmallTests import Thin_Shell3D3M_DrillingRollUpTest  as SHE_3D3N_D_TEST
# Eigenvalues tests
from SmallTests import EigenQ4Thick2x2PlateTests    as TEigenQ4Thick2x2PlateTests
from SmallTests import EigenTL3D8NCubeTests         as TEigenTL3D8NCubeTests
from SmallTests import Eigen3D3NThinCircleTests     as TEigen3D3NThinCircleTests


# Tests for constitutive models: (see ConstitutiveModelsApplication)

## NIGTHLY TESTS

# Constitutive small strains law CSS
from NightlyTests import SmallStrains_IsotropicDamage_SimoJu_Test  as CSS_IsoDamage_SimoJu_TEST


## VALIDATION TESTS
# from ValidationTests import ...




def AssambleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should populate the suites:
    "small", "nighlty" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''
    suites = KratosUnittest.KratosSuites

    ## SMALL TESTS
    #Create a test suit with the selected tests (Small tests):
    smallSuite = suites['small']

    # Patch test small displacement elements SDE
    smallSuite.addTest(SDE_2D4N_S_TEST('test_execution'))
    smallSuite.addTest(SDE_2D3N_S_TEST('test_execution'))
    smallSuite.addTest(SDE_2D4N_T_TEST('test_execution'))
    smallSuite.addTest(SDE_2D3N_T_TEST('test_execution'))
    smallSuite.addTest(SDE_3D8N_T_TEST('test_execution'))
    smallSuite.addTest(SDE_3D4N_T_TEST('test_execution'))
    smallSuite.addTest(SDE_3D8N_T_TEST('test_execution'))
    smallSuite.addTest(SDE_3D4N_T_TEST('test_execution'))   

    # Patch test total lagrangian elements TLE
    smallSuite.addTest(TLE_2D4N_S_TEST('test_execution'))
    smallSuite.addTest(TLE_2D3N_S_TEST('test_execution'))
    smallSuite.addTest(TLE_2D4N_T_TEST('test_execution'))
    smallSuite.addTest(TLE_2D3N_T_TEST('test_execution'))
    smallSuite.addTest(TLE_3D8N_T_TEST('test_execution'))
    smallSuite.addTest(TLE_3D4N_T_TEST('test_execution'))
    smallSuite.addTest(TLE_3D8N_T_TEST('test_execution'))
    smallSuite.addTest(TLE_3D4N_T_TEST('test_execution'))   

    # Patch test updated lagrangian elements ULE
    smallSuite.addTest(ULE_2D4N_S_TEST('test_execution'))
    smallSuite.addTest(ULE_2D3N_S_TEST('test_execution'))
    smallSuite.addTest(ULE_2D4N_T_TEST('test_execution'))
    smallSuite.addTest(ULE_2D3N_T_TEST('test_execution'))
    smallSuite.addTest(ULE_3D8N_T_TEST('test_execution'))
    smallSuite.addTest(ULE_3D4N_T_TEST('test_execution'))
    smallSuite.addTest(ULE_3D8N_T_TEST('test_execution'))
    smallSuite.addTest(ULE_3D4N_T_TEST('test_execution'))   

    # Tests for shell elements SHE
    smallSuite.addTest(SHE_3D4N_B_TEST('test_execution'))
    smallSuite.addTest(SHE_3D4N_D_TEST('test_execution'))
    smallSuite.addTest(SHE_3D3N_B_TEST('test_execution'))
    smallSuite.addTest(SHE_3D3N_D_TEST('test_execution'))

    #...
    if (missing_external_dependencies == False):
        if( hasattr(KratosMultiphysics.ExternalSolversApplication,  "FEASTSolver") ):
            # Eigenvalues tests
            smallSuite.addTest(TEigenQ4Thick2x2PlateTests('test_execution'))
            smallSuite.addTest(TEigen3D3NThinCircleTests('test_execution'))
            smallSuite.addTest(TEigenTL3D8NCubeTests('test_execution'))
        else:
            print("FEASTSolver solver is not included in the compilation of the External Solvers Application")
    

    ## NIGTHLY TESTS
    # Create a test suit with the selected tests for nightly run:   
    nightSuite = suites['nightly']

    # inlude Small Tests in nightly suite
    nightSuite.addTests(smallSuite)
    
    # Test constitutive model
    nightSuite.addTest(CSS_IsoDamage_SimoJu_TEST('test_execution'))

    #...


    
    ## VALIDATION TESTS   
    # Create a test suit with the selected tests for validation: (long tests do not run at night)
    validationSuite = suites['validation']
    
    #...

    

    ## ALL TESTS     
    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    
    allSuite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            #SDE
            SDE_2D4N_S_TEST,
            SDE_2D3N_S_TEST,
            SDE_2D4N_T_TEST,
            SDE_2D3N_T_TEST,
            SDE_3D8N_T_TEST,
            SDE_3D4N_T_TEST,
            SDE_3D8N_T_TEST,
            SDE_3D4N_T_TEST,
            #TLE
            TLE_2D4N_S_TEST,
            TLE_2D3N_S_TEST,
            TLE_2D4N_T_TEST,
            TLE_2D3N_T_TEST,
            TLE_3D8N_T_TEST,
            TLE_3D4N_T_TEST,
            TLE_3D8N_T_TEST,
            TLE_3D4N_T_TEST,
            #ULE
            ULE_2D4N_S_TEST,
            ULE_2D3N_S_TEST,
            ULE_2D4N_T_TEST,
            ULE_2D3N_T_TEST,
            ULE_3D8N_T_TEST,
            ULE_3D4N_T_TEST,
            ULE_3D8N_T_TEST,
            ULE_3D4N_T_TEST,
            #SHE
            SHE_3D4N_B_TEST,
            SHE_3D4N_D_TEST,
            SHE_3D3N_B_TEST,
            SHE_3D3N_D_TEST,
            #CSS
            CSS_IsoDamage_SimoJu_TEST
            ##...
        ])
    )

    if (missing_external_dependencies == False):
        if( hasattr(KratosMultiphysics.ExternalSolversApplication,  "FEASTSolver") ):
            allSuite.addTests(
                KratosUnittest.TestLoader().loadTestsFromTestCases([
                    TEigenQ4Thick2x2PlateTests,
                    TEigenTL3D8NCubeTests,
                    TEigen3D3NThinCircleTests
                ])
            )
        else:
            print("FEASTSolver solver is not included in the compilation of the External Solvers Application")
    
    
    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssambleTestSuites())
