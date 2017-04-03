# import Kratos
import KratosMultiphysics 
import KratosMultiphysics.ExternalSolversApplication 
import KratosMultiphysics.SolidMechanicsApplication 
import KratosMultiphysics.StructuralMechanicsApplication 

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
## SMALL TESTS
from SmallTests import SimpleMeshMovingTest as TSimpleMeshMovingTest
from SmallTests import DynamicBossakTests as TDynamicBossakTests
from SmallTests import DynamicNewmarkTests as TDynamicNewmarkTests
from SmallTests import SprismMembranePatchTests as TSprismMembranePatchTests
from SmallTests import SprismBendingPatchTests as TSprismBendingPatchTests
from SmallTests import ShellQ4ThickBendingRollUpTests as TShellQ4ThickBendingRollUpTests
from SmallTests import ShellQ4ThickDrillingRollUpTests as TShellQ4ThickDrillingRollUpTests
from SmallTests import ShellT3ThinBendingRollUpTests as TShellT3ThinBendingRollUpTests
from SmallTests import ShellT3ThinDrillingRollUpTests as TShellT3ThinDrillingRollUpTests
from SmallTests import EigenQ4Thick2x2PlateTests as TEigenQ4Thick2x2PlateTests
from SmallTests import EigenTL3D8NCubeTests as TEigenTL3D8NCubeTests
from SmallTests import Simple3D2NTrussTest as T3D2NTrussTest
from SmallTests import Simple3D2NBeamCrTest as T3D2NBeamCrTest

## NIGTHLY TESTS
from NightlyTests import ShellT3IsotropicScordelisTests as TShellT3IsotropicScordelisTests

## VALIDATION TESTS
from ValidationTests import SprismPanTests as TSprismPanTests

def AssambleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should pupulate the suites:
    "small", "nighlty" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''
    suites = KratosUnittest.KratosSuites

    # Create a test suit with the selected tests (Small tests):
    smallSuite = suites['small']
    smallSuite.addTest(TSimpleMeshMovingTest('test_execution'))
    smallSuite.addTest(TDynamicBossakTests('test_execution'))
    smallSuite.addTest(TDynamicNewmarkTests('test_execution'))
    smallSuite.addTest(TSprismMembranePatchTests('test_execution'))
    smallSuite.addTest(TSprismBendingPatchTests('test_execution'))
    smallSuite.addTest(TShellQ4ThickBendingRollUpTests('test_execution'))
    smallSuite.addTest(TShellQ4ThickDrillingRollUpTests('test_execution'))
    smallSuite.addTest(TShellT3ThinBendingRollUpTests('test_execution'))
    smallSuite.addTest(TShellT3ThinDrillingRollUpTests('test_execution'))
    smallSuite.addTest(TEigenQ4Thick2x2PlateTests('test_execution'))
    smallSuite.addTest(TEigenTL3D8NCubeTests('test_execution'))
    smallSuite.addTest(T3D2NTrussTest('test_execution'))  
    smallSuite.addTest(T3D2NBeamCrTest('test_execution'))    

    # Create a test suit with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)
    nightSuite.addTest(TShellT3IsotropicScordelisTests('test_execution'))
    
    # For very long tests that should not be in nighly and you can use to validate 
    validationSuite = suites['validation']
    validationSuite.addTest(TSprismPanTests('test_execution'))

    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            #TSimpleMeshMovingTest,
            #TDynamicBossakTests,
            #TDynamicNewmarkTests,
            #TSprismMembranePatchTests,
            #TSprismBendingPatchTests,
            #TShellQ4ThickBendingRollUpTests,
            #TShellQ4ThickDrillingRollUpTests,
            #TShellT3ThinBendingRollUpTests,
            #TShellT3ThinDrillingRollUpTests,
            #TShellT3IsotropicScordelisTests,
            T3D2NTrussTest,
            T3D2NBeamCrTest
            ######TSprismPanTests
        ])
    )
    
    if( hasattr(KratosMultiphysics.ExternalSolversApplication,  "FEASTSolver") ):
        allSuite.addTests(
            KratosUnittest.TestLoader().loadTestsFromTestCases([
                TEigenQ4Thick2x2PlateTests,
                TEigenTL3D8NCubeTests
            ])
        )
    else:
        print("FEASTSolver solver is not included in the compilation of the External Solvers Application")

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssambleTestSuites())
